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

def _getFeatureFamily(mol):
    FEATURE_DEF_FILE = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    feat_factory = ChemicalFeatures.BuildFeatureFactory(FEATURE_DEF_FILE)
    hmol = rdkit.Chem.AddHs(mol)
    AllChem.EmbedMolecule(hmol,useRandomCoords=True)
    rc = rdkit.Chem.AllChem.EmbedMolecule(hmol)
    logging.debug("Getting features for mol "+ mol.GetProp("_Name"))
    if rc < 0 :
        rc = rdkit.Chem.AllChem.EmbedMolecule(hmol,useRandomCoords=True)
    if rc == 0 :
        try :
            if rdkit.Chem.AllChem.UFFOptimizeMolecule(hmol) != 0 :
                rdkit.Chem.AllChem.UFFOptimizeMolecule(hmol,maxIters=1000)
        except ValueError :
             logging.error("Problem with 3D version of molecule "+ hmol.GetProp("_Name"))
             pass
    feats = feat_factory.GetFeaturesForMol(hmol)
    atomFeatures = [["" for feature in range(len(feats))] for atom in range(hmol.GetNumAtoms())]
    for feature in feats:
        for atomId in feature.GetAtomIds():
            if feature.GetFamily() not in atomFeatures[atomId]:
                atomFeatures[atomId].append(feature.GetFamily());
    return atomFeatures
    

def _getElectronegativity(symbol):
    if (symbol == "C"):
        return 2.55
    elif (symbol == "N"):
        return 3.04
    elif (symbol == "Br"):
        return 2.96
    elif (symbol == "O"):
        return 3.44
    elif (symbol == "S"):
        return 2.58
    elif (symbol == "Cl"):
        return 3.16
    elif (symbol == "F"):
        return 3.98
    elif (symbol == "P"):
        return 2.19
    elif (symbol == "I"):
        return 2.66
    elif (symbol == "Si"):
        return 1.9
    elif (symbol == "As"):
        return 2.18
    elif (symbol == "Hg"):
        return 2
    elif (symbol == "Au"):
        return 2.54
    elif (symbol == "Fe"):
        return 1.83
    elif (symbol == "Se"):
        return 2.55
    elif (symbol == "Ga"):
        return 1.81
    else:
        logging.error("Undefined element type:" + symbol)

def _node_to_gxl(gxl, atom, featureList):
    gxl.write('<node id="_' + str(atom.GetIdx()) + '">\n')
    gxl.write('<attr name="valence"><int>' + str(atom.GetExplicitValence()) + '</int></attr>\n')
    gxl.write('<attr name="symbol"><string>' + atom.GetSymbol() + '</string></attr>\n')
    gxl.write('<attr name="electronegativity"><float>' + str(_getElectronegativity(atom.GetSymbol())) + '</float></attr>\n')
    featureString = ""
    for atomFeature in featureList[atom.GetIdx()]:
        if str(atomFeature) != "": featureString += str(atomFeature) + ","
    featureString = featureString[:-1]
    #print(featureString)
    gxl.write('<attr name="pharmacophores"><string>' + featureString + '</string></attr>\n')
    gxl.write('</node>\n')
    return

def _bond_type_to_num(bond_type):
    """ Converts: 
    SINGLE -> 1
    DOUBLE -> 2
    AROMATIC -> 1.5
    """
    if (str(bond_type) == "SINGLE"):
        return 1
    elif (str(bond_type) == "DOUBLE"):
        return 2
    elif (str(bond_type) == "AROMATIC"):
        return 1.5
    elif (str(bond_type) == "TRIPLE"):
        return 3
    else:
        print(bond_type)
        return 0



def _edge_to_gxl(gxl, edge):
    gxl.write('<edge from="_' + str(edge.GetBeginAtomIdx()) + '" to="_' + str(edge.GetEndAtomIdx()) + '">\n')
    gxl.write('<attr name="valence"><float>' + str(_bond_type_to_num(edge.GetBondType())) + '</float></attr>\n') #TODO check how the bond type is written - STRING -> convert to 1, 1.5, 2
    gxl.write('</edge>\n')
    return

def mol_to_gxl(molecule, output_directory):
    """ Transform the loaded molecule into the gxl file."""

    #Generate the gxl file.
    gxlpath = output_directory + "molecules/" + molecule.GetProp("_Name") + ".gxl"
    if (os.path.exists(gxlpath)):
    	return # the file has been already generated 
        #WARNING: this prevents the program to generate the file uselessly, but could cause problems when anybody changed the content of the file or created own with the same name
    gxl = open(gxlpath, 'w+')
    gxl.write('<?xml version="1.0"?>\n')
    gxl.write('<!DOCTYPE gxl SYSTEM "http://www.gupro.de/GXL/gxl-1.0.dtd">\n')
    gxl.write('<gxl>\n')
    gxl.write('<graph id="' + molecule.GetProp("_Name") + '.gxl">\n')
    fl = _getFeatureFamily(molecule)
    for a in molecule.GetAtoms():
	    _node_to_gxl(gxl, a, fl)        
    for e in molecule.GetBonds():
	    _edge_to_gxl(gxl, e)
    gxl.write('</graph>\n')
    gxl.write('</gxl>\n')
    gxl.close();

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
        mol_to_gxl(molecule, output_directory)
    return result

def _create_xml_header(xmlfile):
   xmlfile.write('<?xml version="1.0"?>\n')
   xmlfile.write('<GraphCollection>\n')
   xmlfile.write('<graphs>\n')

def _create_xml_footer(xmlfile):
    xmlfile.write('</graphs>\n')
    xmlfile.write('</GraphCollection>\n')

def generate_mol_sets(input_path, input_directory, output_directory=None):
    """ Generates cxl files for molecules in sdf files.
    The sdf files are specified in .json file in the input path.
    
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
        # During this loading the gxl files are created.
        molecules.update(_load_molecules(path, output_directory))

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


