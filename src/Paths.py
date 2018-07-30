import yaml
from enum import Enum


# All paths are constants and stored in the config directory in 'pathsAndDictionaries' yaml file.
# They are read into the Paths enum from which all other classes source all path strings.
#
# All paths are absolute paths unless prefix with an "REL_" which indicates they are relative
# otherwise they are prefixed with "PATH_"
class Paths(Enum):

    PATH_YAML_FILE = '/Users/u0120577/PycharmProjects/MutateCompute/config/pathsAndDictionaries.yaml'

    with open(PATH_YAML_FILE, 'r') as stream:

        try:

            paths_and_dictionaries = yaml.load(stream)
            PATH_ZEUS_R_EXE = paths_and_dictionaries['ROOT']['path_zeus_r_exe']
            PATH_ZEUS_FOLDX_EXE = paths_and_dictionaries['ROOT']['path_zeus_fx_exe']
            PATH_ZEUS_AGADIR_EXE = paths_and_dictionaries['ROOT']['path_zeus_agad_exe']
            PATH_ZEUS_QSUB_EXE = paths_and_dictionaries['ROOT']['path_zeus_qsub_exe']
            PATH_LOCAL_R_EXE = paths_and_dictionaries['ROOT']['path_local_r_exe']
            PATH_LOCAL_FOLDX_EXE = paths_and_dictionaries['ROOT']['path_local_fx_exe']
            PATH_LOCAL_AGADIR_EXE = paths_and_dictionaries['ROOT']['path_local_agad_exe']

            PATH_ZEUS_SNPEFFECT = paths_and_dictionaries['ROOT']['path_zeus_snpeffect']
            PATH_LOCAL_MUTATECOMPUTE = paths_and_dictionaries['ROOT']['path_local_mutatecompute']

            REL_PATH_CONFIG = paths_and_dictionaries['ROOT']['path_rel_config']
            REL_PATH_CONFIG_AGAD = paths_and_dictionaries['ROOT']['path_rel_config_agadconfig']
            REL_PATH_CONFIG_FX = paths_and_dictionaries['ROOT']['path_rel_config_fxconfig']
            REL_PATH_CONFIG_GLOBAL_OPTIONS = paths_and_dictionaries['ROOT']['path_rel_config_globoptions']
            REL_PATH_CONFIG_FX_ACRUNSCRIPT = paths_and_dictionaries['ROOT']['path_rel_config_fxconf_acrunscript']
            REL_PATH_CONFIG_FX_BMRUNSCRIPT = paths_and_dictionaries['ROOT']['path_rel_config_fxconf_bmrunscript']
            REL_PATH_CONFIG_FX_REPRUNSCRIPT = paths_and_dictionaries['ROOT']['path_rel_config_fxconf_reprunscript']
            REL_PATH_SRC = paths_and_dictionaries['ROOT']['path_rel_src']
            REL_PATH_INPUT = paths_and_dictionaries['ROOT']['path_rel_input']
            REL_PATH_OUTPUT = paths_and_dictionaries['ROOT']['path_rel_output']
            REL_PATH_FX_AC = paths_and_dictionaries['ROOT']['path_rel_fx_ac']
            REL_PATH_FX_BM = paths_and_dictionaries['ROOT']['path_rel_fx_bm']
            REL_PATH_FX_REP = paths_and_dictionaries['ROOT']['path_rel_fx_rep']
            REL_PATH_MUTATE_FASTA = paths_and_dictionaries['ROOT']['path_rel_mutatefasta']

            DICT_AA_1TO3 = paths_and_dictionaries['ROOT']['dict_aa_1to3']
            DICT_AA_3TO1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']
            LIST_ALL_20_AA = paths_and_dictionaries['ROOT']['list_all_20_aa']

        except yaml.YAMLError as exc:
            print(exc)

