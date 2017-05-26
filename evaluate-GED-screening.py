#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import argparse
import os
import json
import logging
import time
import math
import rdkit
import rdkit.Chem
from rdkit.ML.Scoring import Scoring
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

__author__ = 'Aneta Šťastná'
__license__ = 'X11'
__version__ = '1.0.1'


# region Representation and similarity computation
def _ged_similarity(query, active, ged_data):
    """ Gets the Graph Edit Distance from parsed json and normalizes it to <0,1>
    :return: number from <0,1>
    """
    
    # Get number of atoms and bonds of the molecules
    qatoms = query.GetNumAtoms() 
    qbonds = query.GetNumBonds()
    
    aatoms = active.GetNumAtoms()
    abonds = active.GetNumBonds()

    # Count the largest possible GED score
    bondCoef = ged_data["properties"]["edge delete/insert"]
    atomCoef = ged_data["properties"]["node delete/insert"]
    alpha = ged_data["properties"]["node-edge weight ratio"]
    maxged = alpha * atomCoef * (qatoms + aatoms) + (1 - alpha) * bondCoef * (qbonds + abonds)

    # Get the real GED score
    activeId = ged_data["results"]["matrix"]["targets"].index(active.GetProp("_Name"))
    ged = ged_data["results"]["matrix"][query.GetProp("_Name")][activeId]
    
    return (1 - ged/maxged)


def _load_molecules(sdf_path):
    """Load molecules from SDF file.

    :param sdf_path:
    :return: Molecules in dictionary under their names.
    """
    result = {}
    # Convert path to str os in come cases (Python 2.7, Win 7) the sdf_path
    # can be of type unicode and tha call would fail on invalid
    # argument type.
    for molecule in rdkit.Chem.SDMolSupplier(str(sdf_path)):
        if molecule is None:
            logging.error("Can't load molecule.")
            continue
        result[molecule.GetProp('_Name')] = molecule
    return result


def screening(input_path, input_directory, ged_results_file, output_path=None):
    """Perform a virtual screening.

    :param input_path: input .json file with basic screening params
    :param input_directory: directory with .sdf files
    :param ged_results_file: .json file with GED results and parameters
    :param output_path:
    :return:
    """
    with open(input_path) as input_stream:
        input_data = json.load(input_stream)
    # Load molecules.
    logging.info('Loading molecules ...')
    molecules = {}
    for file_item in input_data['files']:
        path = input_directory + file_item + '.sdf'
        if not os.path.exists(path):
            logging.error('Missing file: %s' % path)
            raise Exception('Missing file.')
        molecules.update(_load_molecules(path))
    # Parse ged results file
    with open(ged_results_file) as ged_stream:
        ged_data = json.load(ged_stream)
    # Create representation of active molecules.
    actives = []
    for active in input_data['data']['train']['ligands']:
        if active['name'] not in molecules:
            continue
        actives.append(molecules[active['name']])
    # Screening.
    logging.info('Screening ...')
    scores = []
    counter = 0
    counter_max = len(input_data['data']['test'])
    counter_step = math.floor(counter_max / 100.0) + 1
    time_begin = time.clock()
    for item in input_data['data']['test']:
        if item['name'] not in molecules:
            continue
        query = molecules[item['name']]
        # Counting similarity and searching for most similar active molecule
        similarity = 0
        similarMol = query
        for active in actives:
            currentSimilarity = _ged_similarity(query, active, ged_data)
            if (currentSimilarity > similarity):
                similarity = currentSimilarity
                similarMol = active
        scores.append({
            'name': item['name'],
            'similarity': similarity,
            'activity': item['activity'],
            'most-similar-active': similarMol.GetProp("_Name")
        })
        #if (item['activity'] == 1) create_picture(query, similar-active)
        if counter % counter_step == 0:
            logging.debug('%d/%d', counter, counter_max)            
        counter += 1
        logging.debug('counter: ' + str(counter))
    time_end = time.clock()
    # Evaluate screening.
    scores = sorted(scores,
                    key=lambda m: m['similarity'],
                    reverse=True)
    auc = Scoring.CalcAUC(scores, 'activity')
    ef = Scoring.CalcEnrichment(scores, 'activity', [0.005, 0.01, 0.02, 0.05])
    # Print results.
    print('AUC : ', auc)
    print('EF (0.5%, 1.0%, 2.0%, 5.0%) : ', ef)
    total_time = float(ged_data["properties"]["time"])/1000
    total_time += (time_end - time_begin)
    print('Execution time : %.2fs' % total_time)
    # Write result to a file.
    if not output_path is None and not output_path == '':
        if not os.path.exists(os.path.dirname(output_path)):
            os.makedirs(os.path.dirname(output_path))
        with open(output_path, 'w') as output_stream:
            json.dump({
                'properties': ged_data["properties"],
                'data': scores,
                'metadata': {
                    'auc': auc,
                    'ef': {
                        '0.005': ef[0],
                        '0.01': ef[1],
                        '0.02': ef[2],
                        '0.05': ef[3]
                    },
                    'fileName': os.path.basename(__file__),
                    'executionTime': total_time,
                    
                'definition': {
                        'selection': input_data['info']['selection'],
                        'molecules': input_data['info']['molecules'],
                        'index': input_data['info']['index'],
                        'dataset': input_data['info']['dataset'],
                        'method': input_data['info']['method'],
                        'config': 'config_file'
                    }
                }
            }, output_stream, indent=2)


def _read_configuration():
    """Read command line arguments and return configuration object.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='Perform a simple virtual screening.')
    parser.add_argument('-i', type=str, dest='input_path',
                        help='Path to input data (training and test in .json).',
                        required=True)
    parser.add_argument('-j', type=str, dest='input_directory',
                        help='Path to directory with input molecules (directory with .sdf files).',
                        required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='Path to the output file (output.json).',
                        required=True)
    parser.add_argument('-g', type=str, dest='ged_results',
                        help='Path to the file with counted graph edit distance (.json).',
                        required=True)

    return vars(parser.parse_args())


def _main():
    """Application entry point.

    """
    # Initialize logging.
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    #
    config = _read_configuration()
    screening(config['input_path'], config['input_directory'], config['ged_results'], config['output'])


if __name__ == '__main__':
    _main()
