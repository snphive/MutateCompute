#!/usr/bin/env python3
"""
Enum for amino acid lists and dictionaries. Read from yaml file.
"""

import yaml
import os
from enum import Enum

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


# Name of enum may change from AminoAcids.AA if other types of data are stored here in future.
class AA(Enum):

    cwd = os.getcwd()
    YAML_FILE = '/Users/u0120577/PycharmProjects/MutateCompute/config/pathsAndDictionaries.yaml' if not \
        ('switchlab' in cwd or 'shazib' in cwd) else \
        '/switchlab/group/shazib/SnpEffect/config/pathsAndDictionaries.yaml'

    with open(YAML_FILE, 'r') as stream:

        try:

            paths_and_dictionaries = yaml.load(stream)
            DICT_AA_1TO3 = paths_and_dictionaries['ROOT']['dict_aa_1to3']
            DICT_AA_3TO1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']
            LIST_ALL_20_AA = paths_and_dictionaries['ROOT']['list_all_20_aa']

        except yaml.YAMLError as exc:
            print(exc)

