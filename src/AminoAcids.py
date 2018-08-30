import yaml
import os
from enum import Enum


# Name of enum may change from AminoAcids.AA if other types of data are stored here in future.
class AA(Enum):

    cwd = os.getcwd()
    YAML_FILE = '/Users/u0120577/PycharmProjects/MutateCompute/configuration/pathsAndDictionaries.yaml' if not \
        ('switchlab' in cwd or 'shazib' in cwd) else \
        '/switchlab/group/shazib/SnpEffect/configuration/pathsAndDictionaries.yaml'

    with open(YAML_FILE, 'r') as stream:

        try:

            paths_and_dictionaries = yaml.load(stream)
            DICT_AA_1TO3 = paths_and_dictionaries['ROOT']['dict_aa_1to3']
            DICT_AA_3TO1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']
            LIST_ALL_20_AA = paths_and_dictionaries['ROOT']['list_all_20_aa']

        except yaml.YAMLError as exc:
            print(exc)

