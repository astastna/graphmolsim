#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import argparse
import os
import json
import logging
import rdkit
import rdkit.Chem
from rdkit.ML.Scoring import Scoring
from rdkit.Chem import rdFMCS
#from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit import Geometry
from rdkit import RDConfig
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib


__author__ = 'Aneta Šťastná'
__license__ = ''
__version__ = '1.0'

def _load_molecules(sdf_path, output_directory):
    """Load molecules from SDF file.

    :param sdf_path:
    :return: Molecules in dictionary under their names.
    """
    result = {}
    # Convert path to str os in come cases (Python 2.7, Win 7) the sdf_path
    # can be of type unicode and that call would fail on invalid
    # argument type.
    for molecule in rdkit.Chem.SDMolSupplier(str(sdf_path)):
        if molecule is None:
            logging.error("Can't load molecule.")
            continue
        result[molecule.GetProp('_Name')] = molecule
    return result

def _create_xml_header(xmlfile):
   xmlfile.write('<?xml version="1.0"?>\n')
   xmlfile.write('<GraphCollection>\n')
   xmlfile.write('<graphs>\n')

def _create_xml_footer(xmlfile):
    xmlfile.write('</graphs>\n')
    xmlfile.write('</GraphCollection>\n')

def generate_mol_sets(input_path, input_directory, output_directory=None):
    """Generates two xml files which can be handed to the GraphMatchingToolkit.
    source.cxl - file with all molecules from test
    target.cxl - file with all active molecules

    :param input_path:
    :param input_directory:
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
	logging.debug(path)
        if not os.path.exists(path):
            logging.error('Missing file: %s' % file_item)
            raise Exception('Missing file.')
        molecules.update(_load_molecules(path, output_directory))
    
    # Create representation of active molecules.
    actives = []
    for active in input_data['data']['train']['ligands']:
        if active['name'] not in molecules:
            continue
        actives.append(molecules[active['name']])

    # Preparing the screening by producing the .cxl files for GED algorithm.
    logging.info('Converting to cxl ...')
   
    # Creating the source cxl file with all the molecules for testing 
    input_file_name = os.path.basename(input_path).split(".")[0] 
    source_set = open(output_directory + "/source.cxl", "w")
    _create_xml_header(source_set)
    for item in input_data['data']['test']:
        if item['name'] not in molecules:
            continue
        query = molecules[item['name']]
        source_set.write('<print file="' + query.GetProp("_Name") + '.gxl"/>\n')
    _create_xml_footer(source_set)
    source_set.close()

    # Creating the target cxl file with all active molecules
    tartget_set = open(output_directory + "/target.cxl", "w")
    _create_xml_header(tartget_set)
    for active in actives:
        tartget_set.write('<print file="' + active.GetProp("_Name") + '.gxl"/>\n')
    _create_xml_footer(tartget_set)
    tartget_set.close()
    # GED algorithm then calculates GED between each pair (souce, target)

def _read_configuration():
    """Read command line arguments and return configuration object.

    :return:
    """
    parser = argparse.ArgumentParser(
        description='Generates molecule sets for running GraphEditDistance.')
    parser.add_argument('-i', type=str, dest='input_path',
                        help='Path to input data (training and test in .json).',
                        required=True)
    parser.add_argument('-j', type=str, dest='input_directory',
                        help='Path to directory with input molecules (directory with .sdf files).',
                        required=True)
    parser.add_argument('-o', type=str, dest='output_directory',
                        help='Path to the output directory.',
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
    generate_mol_sets(config['input_path'], config['input_directory'], config['output_directory'])

if __name__ == '__main__':
    _main()


