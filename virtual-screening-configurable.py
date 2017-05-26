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
import mcsutils
from os.path import basename
from rdkit.ML.Scoring import Scoring
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

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

def _flush_results(output_path, scores):
    """ Flushes the current results to the output file.
    
    :param output_path: path to output file
    :param scores: current results
    """
    if not output_path is None and not output_path == '':
        if not os.path.exists(os.path.dirname(output_path)):
            os.makedirs(os.path.dirname(output_path))
        with open(output_path, 'w') as output_stream:
            json.dump({
                'data': scores
            }, output_stream, indent=2)
    return

def screening(input_dir, input_directory, config_file, output_path=None):
    """Perform a virtual screening.

    :param input_dir: path to input data (training and test in .json)
    :param input_directory: path to sdf files
    :param config_file: configuration file of mcs
    :param output_path: directory to save the results
    :return:
    """
    
    with open(input_dir) as input_stream:
        input_data = json.load(input_stream)
    # Load molecules.
    logging.info('Loading molecules ...')
    molecules = {}
    for file_item in input_data['files']:
        path = input_directory + file_item + '.sdf'
        if not os.path.exists(path):
            logging.error('Missing file: %s' % file_item)
            raise Exception('Missing file.')
        molecules.update(_load_molecules(path))
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
    inexact = 0
    counter_max = len(input_data['data']['test'])
    counter_step = math.floor(counter_max / 100.0) + 1
    params = mcsutils._parse_config(config_file)
    time_begin = time.clock()
    for item in input_data['data']['test']:
        if item['name'] not in molecules:
            continue
        query = molecules[item['name']]
        similarity = max([mcsutils._similarity(query, active, inexact, input_data['info'], params) for active in actives])
        scores.append({
            'name': item['name'],
            'similarity': similarity,
            'activity': item['activity']
        })
        if counter % counter_step == 0:
            logging.debug('%d/%d', counter, counter_max)            
            #_flush_results(output_path, scores)
        counter += 1
        #logging.debug('counter: ' + str(counter))
    time_end = time.clock()
    # Evaluate screening.
    scores = sorted(scores,
                    key=lambda m: m['similarity'],
                    reverse=True)
    auc = Scoring.CalcAUC(scores, 'activity')
    ef = Scoring.CalcEnrichment(scores, 'activity', [0.005, 0.01, 0.02, 0.05])
    # Print results.
    print('Input file: ', input_dir)
    print('Difficulty: ', input_directory)
    print('AUC : ', auc)
    print('EF (0.5%, 1.0%, 2.0%, 5.0%) : ', ef)
    print('Execution time : %.2fs' % (time_end - time_begin))
    # Write result to a file.
    if not output_path is None and not output_path == '':
        if not os.path.exists(os.path.dirname(output_path)):
            os.makedirs(os.path.dirname(output_path))
        with open(output_path, 'w') as output_stream:
            json.dump({
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
                    'executionTime': time_end - time_begin,
                    'inexactMolecules': inexact,
                    
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
    parser.add_argument('-i', type=str, dest='input_dir',
                        help='Path to input data (training and test in .json).',
                        required=True)
    parser.add_argument('-j', type=str, dest='input_directory',
                        help='Path to directory with input molecules (directory with .sdf files).',
                        required=True)
    parser.add_argument('-o', type=str, dest='output',
                        help='Path to the output file (output.json).',
                        required=True)
    parser.add_argument('-c', type=str, dest='config',
                        help='Path to the config file (.json).',
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
    screening(config['input_dir'], config['input_directory'], config['config'], config['output'])


if __name__ == '__main__':
    _main()
