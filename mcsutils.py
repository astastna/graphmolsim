#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


author__ = 'Aneta Šťastná'
__license__ = ''
__version__ = '1.0'

def _similarity(left, right, inexact, info, params):
    """Computes molecule similarity through the mcs
    with the parameters given.

    :param left: molecule
    :param right: molecule
    :param path: path for output pictures
    :param info: information about current dataset
    :param params: parameters for calling
    :return: Value from range <0,1>
    """
    total = max(left.GetNumAtoms(), right.GetNumAtoms())
    mols = [left, right]

    maxCommonSubstructure = rdFMCS.FindMCS(mols, atomCompare=params[0], bondCompare=params[1], maximizeBonds=params[2], matchValences=params[3], ringMatchesRingOnly=params[4], completeRingsOnly=params[5], timeout=20)
    in_common = maxCommonSubstructure.numAtoms
    if maxCommonSubstructure.canceled:
        inexact = inexact+1

    return float(in_common) / float(total)

def _parse_config(config_path):
    """ Parse the mcs configuration from config file (.json)
    :param config_path

    :return list with objects for calling similarity
    """
    with open(config_path) as config_stream:
        config = json.load(config_stream)

    # Create the objects for calling of the MCS function according to the config file
    
    # Creating atomCompare object
    if config['atomCompare'] == 'CompareAny':
        atomCompareConf = rdFMCS.AtomCompare.CompareAny
    elif config['atomCompare'] == 'CompareElements':
        atomCompareConf = rdFMCS.AtomCompare.CompareElements
    elif config['atomCompare'] == 'CompareIsotopes':
        atomCompareConf = rdFMCS.AtomCompare.CompareIsotopes
    else:
        logging.error('Wrong config file option for atomCompare. Using default one.')
        atomCompareConf = rdFMCS.AtomCompare.CompareElements
    # Creating bondCompare object
    if config['bondCompare'] == 'CompareAny':
        bondCompareConf = rdFMCS.BondCompare.CompareAny
    elif config['bondCompare'] == 'CompareOrder':
        bondCompareConf = rdFMCS.BondCompare.CompareOrder
    elif config['bondCompare'] == 'CompareOrderExact':
        bondCompareConf = rdFMCS.BondCompare.CompareOrderExact
    else:
        logging.error('Wrong config file option for bondCompare. Using default one.')
        bondCompareConf = rdFMCS.BondCompare.CompareOrder
    maximizeBondsConf = config['maximizeBonds'] #boolean
    matchValencesConf = config['matchValences'] #boolean
    ringMatchesRingOnlyConf = config['ringMatchesRingOnly'] #boolean
    completeRingsOnlyConf = config['completeRingsOnly'] #boolean
    
    # Adding the objects to output list
    params = [atomCompareConf, bondCompareConf, maximizeBondsConf, matchValencesConf, ringMatchesRingOnlyConf, completeRingsOnlyConf]

    return params

