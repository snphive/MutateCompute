import yaml
from enum import Enum


# All paths are constants and stored in the config directory in 'pathsAndDictionaries' yaml file.
# They are read into the Paths enum from which all other classes source all path strings.
#
# All paths are absolute paths unless prefix with an "REL_" which indicates they are relative paths.
class Paths(Enum):

    YAML_FILE = '/Users/u0120577/PycharmProjects/MutateCompute/config/pathsAndDictionaries.yaml'

    with open(YAML_FILE.value, 'r') as stream:

        try:

            paths_and_dictionaries = yaml.load(stream)
            ZEUS_R_EXE = paths_and_dictionaries['ROOT']['path_zeus_r_exe']
            ZEUS_FOLDX_EXE = paths_and_dictionaries['ROOT']['path_zeus_fx_exe']
            ZEUS_AGADIR_EXE = paths_and_dictionaries['ROOT']['path_zeus_agad_exe']
            ZEUS_QSUB_EXE = paths_and_dictionaries['ROOT']['path_zeus_qsub_exe']
            LOCAL_R_EXE = paths_and_dictionaries['ROOT']['path_local_r_exe']
            LOCAL_FOLDX_EXE = paths_and_dictionaries['ROOT']['path_local_fx_exe']
            LOCAL_AGADIR_EXE = paths_and_dictionaries['ROOT']['path_local_agad_exe']

            ZEUS_SNPEFFECT = paths_and_dictionaries['ROOT']['path_zeus_snpeffect']
            LOCAL_MUTATECOMPUTE = paths_and_dictionaries['ROOT']['path_local_mutatecompute']

            REL_CONFIG = paths_and_dictionaries['ROOT']['path_rel_config']
            REL_CONFIG_AGAD = paths_and_dictionaries['ROOT']['path_rel_config_agadconfig']
            REL_CONFIG_JOBQ = paths_and_dictionaries['ROOT']['path_rel_config_clustjobq']
            REL_CONFIG_FX = paths_and_dictionaries['ROOT']['path_rel_config_fxconfig']
            REL_CONFIG_GLOBAL_OPTIONS = paths_and_dictionaries['ROOT']['path_rel_config_globoptions']
            REL_CONFIG_FX_ACRUNSCRIPT = paths_and_dictionaries['ROOT']['path_rel_config_fxconf_acrunscript']
            REL_CONFIG_FX_BMRUNSCRIPT = paths_and_dictionaries['ROOT']['path_rel_config_fxconf_bmrunscript']
            REL_CONFIG_FX_REPRUNSCRIPT = paths_and_dictionaries['ROOT']['path_rel_config_fxconf_reprunscript']
            REL_SRC = paths_and_dictionaries['ROOT']['path_rel_src']
            REL_INPUT = paths_and_dictionaries['ROOT']['path_rel_input']
            REL_OUTPUT = paths_and_dictionaries['ROOT']['path_rel_output']
            REL_FX_AC = paths_and_dictionaries['ROOT']['path_rel_fx_ac']
            REL_FX_BM = paths_and_dictionaries['ROOT']['path_rel_fx_bm']
            REL_FX_REP = paths_and_dictionaries['ROOT']['path_rel_fx_rep']
            REL_MUTATE_FASTA = paths_and_dictionaries['ROOT']['path_rel_mutatefasta']

            DICT_AA_1TO3 = paths_and_dictionaries['ROOT']['dict_aa_1to3']
            DICT_AA_3TO1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']
            LIST_ALL_20_AA = paths_and_dictionaries['ROOT']['list_all_20_aa']

            # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
            MC_CONFIG = LOCAL_MUTATECOMPUTE.value + REL_CONFIG.value
            MC_CONFIG_AGAD = LOCAL_MUTATECOMPUTE.value + REL_CONFIG_AGAD.value
            MC_CONFIG_JOBQ = LOCAL_MUTATECOMPUTE.value + REL_CONFIG_JOBQ.value
            MC_CONFIG_FX = LOCAL_MUTATECOMPUTE.value + REL_CONFIG_FX.value
            MC_CONFIG_GLOBAL_OPTIONS = LOCAL_MUTATECOMPUTE.value + REL_CONFIG_GLOBAL_OPTIONS.value
            MC_CONFIG_FX_ACRUNSCRIPT = LOCAL_MUTATECOMPUTE.value + REL_CONFIG_FX_ACRUNSCRIPT.value
            MC_CONFIG_FX_BMRUNSCRIPT = LOCAL_MUTATECOMPUTE.value + REL_CONFIG_FX_BMRUNSCRIPT.value
            MC_INPUT = LOCAL_MUTATECOMPUTE.value + REL_INPUT.value
            MC_OUTPUT = LOCAL_MUTATECOMPUTE.value + REL_OUTPUT.value

            SE_CONFIG = ZEUS_SNPEFFECT.value + REL_CONFIG.value
            SE_CONFIG_AGAD = ZEUS_SNPEFFECT.value + REL_CONFIG_AGAD.value
            SE_CONFIG_JOBQ = ZEUS_SNPEFFECT.value + REL_CONFIG_JOBQ.value
            SE_CONFIG_FX = ZEUS_SNPEFFECT.value + REL_CONFIG_FX.value
            SE_CONFIG_GLOBAL_OPTIONS = ZEUS_SNPEFFECT.value + REL_CONFIG_GLOBAL_OPTIONS.value
            SE_CONFIG_FX_ACRUNSCRIPT = ZEUS_SNPEFFECT.value + REL_CONFIG_FX_ACRUNSCRIPT.value
            SE_CONFIG_FX_BMRUNSCRIPT = ZEUS_SNPEFFECT.value + REL_CONFIG_FX_BMRUNSCRIPT.value
            SE_INPUT = ZEUS_SNPEFFECT.value + REL_INPUT.value
            SE_OUTPUT = ZEUS_SNPEFFECT.value + REL_OUTPUT.value

        except yaml.YAMLError as exc:
            print(exc)

