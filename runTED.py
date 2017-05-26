#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
import sys
sys.path.append("zhang-shasha/zss/")
import argparse
import os
import json
import logging
import rdkit
import math
import time
import rdkit.Chem
from rdkit.ML.Scoring import Scoring
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit import Geometry
from rdkit import RDConfig
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib
import compare
import simple_tree

__author__ = 'Aneta Šťastná'
__license__ = ''
__version__ = '1.0'



def _getFeatureFamily(mol):
    FEATURE_DEF_FILE = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    feat_factory = ChemicalFeatures.BuildFeatureFactory(FEATURE_DEF_FILE)
    hmol = rdkit.Chem.AddHs(mol)
    AllChem.EmbedMolecule(hmol,useRandomCoords=True)
    rc = rdkit.Chem.AllChem.EmbedMolecule(hmol)
    #logging.debug("Getting features for mol "+ mol.GetProp("_Name"))
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

def absolute(a, b):
    return abs(a - b)

def any_in_common(a, b):
    for pharmacofore in a:
        if pharmacofore in b:
            return 1
    return 0

def _ted(query, active, prop):
    """ Gets the TED.

    :query: tree
    :active: tree
    :return: ted
    """
    # Setting the correct distance function and label values for empty node
    if (prop == 'symbol'):
        dist = compare.strdist
        empty = ""
    elif (prop == 'valence'):
        dist = absolute
        empty = 0
    elif (prop == 'pharmacophore'):
        dist = any_in_common 
        empty = [] 
    elif (prop == 'electronegativity'):
        dist = absolute
        empty = 0
    else:
        logging.error('Wrong prop type! Behaving like prop=symbol')
        logging.debug('prop='+ prop)
        dist = compare.strdist
        empty = ""

    ted = compare.distance(
        query, active, simple_tree.Node.get_children,
        insert_cost=lambda node: dist(empty, simple_tree.Node.get_label(node)),
        remove_cost=lambda node: dist(simple_tree.Node.get_label(node), empty),
        update_cost=lambda a, b: dist(simple_tree.Node.get_label(a), simple_tree.Node.get_label(b)) 
    )

    logging.debug("Tree edit distance: "+ str(ted))    
    return ted

def mol_to_tree(mol):
    """ Converts the given molecule to image and saves it into a file.
    :return: 
    """
    molName = mol.GetProp('_Name')

    # canonize the given molecule and create read-write molecule from the original one
    mol = rdkit.Chem.MolFromSmiles(rdkit.Chem.MolToSmiles(mol))
    tree = rdkit.Chem.RWMol(mol)
    # create a tree from the molecule
    atomIds = [ a.GetIdx() for a in tree.GetAtoms()]
    atomIds.sort()
    # atom ids are sorted, we want to start with the highest id
    for atomId in atomIds:
        atom = tree.GetAtomWithIdx(atomId)
        # skip atoms not on rings
        if not atom.IsInRing():
            continue
        # check whether neighbours with lower ids are also on ring
        neighbors = [ n.GetIdx() for n in atom.GetNeighbors()]
        neighbors.sort(reverse=True)
        neighborsId = [n for n in neighbors]
        for neigh in neighbors:
          if atom.GetIdx() > neigh and tree.GetAtomWithIdx(neigh).IsInRing():
                # and remove the bond if so
                tree = rdkit.Chem.RWMol(tree)
                tree.RemoveBond(atom.GetIdx(), neigh)
                tree = rdkit.Chem.AddHs(tree)
                # and update the molecule
                #rdkit.Chem.SanitizeMol(tree)

    tree.SetProp('_Name', molName)
    # preparation for picture creation
    AllChem.Compute2DCoords(tree)
    AllChem.Compute2DCoords(mol)

#    img=Draw.MolsToGridImage([mol, tree], molsPerRow=2, subImgSize=(300,300))
#    img.save('img_tree_test2/'+ molName +'-tree.png')
   
    return tree

def tree_to_graph(moltree, prop):
    """ This function converts the Python Mol object into the MolTree object
        which can be used in Zhang Sasha algorithm implementation.
        The function implicitly expects that the molecular graph is connected.

        :moltree: rdKit Mol object
        :return: SimpleTree object 
    """
    nodes = {}
    parent = {}
    # initialize the values for root
    root = moltree.GetAtomWithIdx(0)
    # creating the node according to given property
    if (prop == 'symbol'):
        root_node = simple_tree.Node(root.GetSymbol())
    elif (prop == 'valence'):
        root_node = simple_tree.Node(root.GetTotalValence())
    elif (prop == 'pharmacophore'):
        atomPharmacophores = _getFeatureFamily(moltree)
        root_node = simple_tree.Node(atomPharmacophores[0])
    elif (prop == 'electronegativity'):    
        root_node = simple_tree.Node(_getElectronegativity(root.GetSymbol()))
    else:
        logging.error('Wrong prop type! Behaving like prop=symbol')
        root_node = simple_tree.Node(root.GetSymbol())
    parent[0] = 0
    nodes[0] = root_node
    stack = [0]
    while (len(stack) > 0):
        # get the number of atom for processing and the corresponding Node object
        curr_id = stack.pop()
        curr_atom = moltree.GetAtomWithIdx(curr_id)
        curr_node = nodes[curr_id]
        
        # add all neighboring atoms as sons and add them to the stack for processing
        for n in curr_atom.GetNeighbors():
            neigh = n.GetIdx()
            # we mustn't add the parent to avoid neverending cycle
            if (neigh == parent[curr_id]):
                continue
            neigh_atom = moltree.GetAtomWithIdx(neigh)
            # creating the node according to given property
            if (prop == 'symbol'):
                neigh_node = simple_tree.Node(neigh_atom.GetSymbol())
            elif (prop == 'valence'):
                neigh_node = simple_tree.Node(neigh_atom.GetTotalValence())
            elif (prop == 'pharmacophore'):
                neigh_node = simple_tree.Node(atomPharmacophores[neigh])
            elif (prop == 'electronegativity'):
                neigh_node = simple_tree.Node(_getElectronegativity(neigh_atom.GetSymbol()))
            else:
                logging.error('Wrong prop type! Behaving like prop=symbol')
                logging.debug('prop='+ prop)
                neigh_node = simple_tree.Node(neigh_atom.GetSymbol())
            parent[neigh] = curr_id
            nodes[neigh] = neigh_node
            curr_node.addkid(neigh_node)
            stack.append(neigh)
    return nodes[0]            
    
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

def _load_molecules(sdf_path, sizes, bondSizes, prop):
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
        tree = mol_to_tree(molecule)
        graph = tree_to_graph(tree, prop)
        result[molecule.GetProp('_Name')] = graph
        sizes[molecule.GetProp('_Name')] = molecule.GetNumAtoms()
        bondSizes[molecule.GetProp('_Name')] = molecule.GetNumBonds()
    return result

def run_ted(input_path, input_directory, prop, output_path):
    """ Loads .sdf file, converts the molecules into trees with graph annotations, runs
    the TED, evaluates the results and saves them into a file.

    :param input_path:
    :param input_directory:
    :param output_path:
    :return:
    """
    with open(input_path) as input_stream:
        input_data = json.load(input_stream)
        
    # Load molecules and convert them to tree graphs.
    logging.info('Loading molecules ...')
    molecules = {}
    sizes = {}
    bondSizes = {}
    for file_item in input_data['files']:
        path = input_directory + file_item + '.sdf'
        logging.debug(path)
        if not os.path.exists(path):
            logging.error('Missing file: %s' % file_item)
            raise Exception('Missing file.')
        molecules.update(_load_molecules(path, sizes, bondSizes, prop))
    
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
        query_size = sizes[item['name']]
        query_bonds = bondSizes[item['name']]
        # Count pairwise similarity with all actives and choose the maximum.
        maxsim = 0
        for active in input_data['data']['train']['ligands']:
            if active['name'] not in molecules:
                continue
            active_graph = molecules[active['name']]
            active_size = sizes[active['name']]
            active_bonds = bondSizes[active['name']]
            ted = _ted(query, active_graph, prop)
            sim = 1.00 - ted / float(query_size + active_size + query_bonds + active_bonds)
            if (sim > maxsim):
                maxsim = sim
                minted = ted
        scores.append({
            'name': item['name'],
            'similarity': maxsim,
            'activity': item['activity'],
            'ted': minted 
        })
        if counter % counter_step == 0:
            logging.debug('%d/%d', counter, counter_max)
            _flush_results(output_path, scores)
        counter += 1
        logging.debug('counter: ' + str(counter))
    time_end = time.clock()
    logging.debug("Reached the end.")

    # Evaluate screening.
    scores = sorted(scores,
                    key=lambda m: m['similarity'],
                    reverse=True)
    auc = Scoring.CalcAUC(scores, 'activity')
    ef = Scoring.CalcEnrichment(scores, 'activity', [0.005, 0.01, 0.02, 0.05])

    # Print results.
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
        description='Runs and evaluates the TreeEditDistance.')
    parser.add_argument('-i', type=str, dest='input_path',
                        help='Path to input data (training and test in .json).',
                        required=True)
    parser.add_argument('-j', type=str, dest='input_directory',
                        help='Path to directory with input molecules (directory with .sdf files).',
                        required=True)
    parser.add_argument('-p', type=str, dest='prop',
                        help='The type of label used on atoms: symbol, valence, pharmacophore, electronegativity.',
                        required=True)
    parser.add_argument('-o', type=str, dest='output_path',
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
    if (config['prop'] not in ['valence', 'symbol', 'pharmacophore','electronegativity']):
        logging.error('Wrong prop argument type! Must be one of: symbol, valence, pharmacophore, electronegativity. Ending.')
        return 
    run_ted(config['input_path'], config['input_directory'], config['prop'], config['output_path'])

if __name__ == '__main__':
    _main()


