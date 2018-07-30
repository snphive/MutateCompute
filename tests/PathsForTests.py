from enum import Enum


# All paths are constants and stored in the config directory in 'pathsAndDictionaries' yaml file.
# They are read into the Paths enum from which all other classes source all path strings.
# However for tests, the absolute paths are directed towards the "tests" directory, so do not access the yaml config.
#
# All paths are absolute paths unless prefix with an "REL" which indicates they are relative.
# All paths are prefixed with "PATH"
class PFT(Enum):

    # ABSOLUTE PATHS FOR TESTS ROOT
    PATH_TESTS = "/Users/u0120577/PycharmProjects/MutateCompute/tests"

    # CONFIGURATION FILES FROM MAIN DIRECTORY (i.e. real data)
    PATH_CONFIG_FOR_READ_ONLY = "/Users/u0120577/PycharmProjects/MutateCompute/config"
    # INPUT FILES FROM MAIN DIRECTORY (i.e. real data)
    PATH_INPUT_FOR_READ_ONLY = "/Users/u0120577/PycharmProjects/MutateCompute/input_data"

    # The following might not even be used as entire config directory is being copied straight into /tests in one go.
    REL_PATH_CONFIG = "/config"
    REL_PATH_CONFIG_AGAD = "/config/agadir_config"
    REL_PATH_CONFIG_FX = "/config/foldx_config"
    REL_PATH_CONFIG_GLOBAL_OPTIONS = "/config/global_options"
    REL_PATH_CONFIG_FX_ACRUNSCRIPT = "/config/foldx_config/ac_runscript"
    REL_PATH_CONFIG_FX_BMRUNSCRIPT = "/config/foldx_config/bm_runscript"
    # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
    PATH_TESTS_CONFIG = PATH_TESTS.value + REL_PATH_CONFIG.value
    PATH_TESTS_CONFIG_AGAD = PATH_TESTS.value + REL_PATH_CONFIG_AGAD.value
    PATH_TESTS_CONFIG_FX = PATH_TESTS.value + REL_PATH_CONFIG_FX.value
    PATH_TESTS_CONFIG_GLOBAL_OPTIONS = PATH_TESTS.value + REL_PATH_CONFIG_GLOBAL_OPTIONS.value
    PATH_TESTS_CONFIG_FX_ACRUNSCRIPT = PATH_TESTS.value + REL_PATH_CONFIG_FX_ACRUNSCRIPT.value
    PATH_TESTS_CONFIG_FX_BMRUNSCRIPT = PATH_TESTS.value + REL_PATH_CONFIG_FX_BMRUNSCRIPT.value

    # REFERENCE FILES ONLY
    REL_PATH_REF_FILES = "/reference_files"
    # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
    PATH_TESTS_REF_FILES = PATH_TESTS.value + REL_PATH_REF_FILES.value

    # INPUT PATHS ONLY
    # e.g. /tests/input_data/<pdbname>/all_mutants/<fxmutantchainname>
    REL_PATH_INPUT = "/input_data"
    REL_PATH_ALL_MUTANTS = "/all_mutants"
    # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
    PATH_TESTS_INPUT = PATH_TESTS.value + REL_PATH_INPUT.value

    # OUTPUT PATHS ONLY
    # e.g.#1: tests/output_data/<pdbname>/foldx/build_model/<fxmutantchainname>
    # e.g.#2: tests/output_data/<pdbname>/mutate_fasta
    REL_PATH_OUTPUT = "/output_data"
    REL_PATH_FX_AC = "/foldx/analyse_complex"
    REL_PATH_FX_BM = "/foldx/build_model"
    REL_PATH_MUTATE_FASTA = "/mutate_fasta"
    # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
    PATH_TESTS_OUTPUT = PATH_TESTS.value + REL_PATH_OUTPUT.value

    DICT_AA_1TO3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
                    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}
    DICT_AA_3TO1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'H1S': 'H',
                    'H2S': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
                    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
    LIST_ALL_20_AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                      'Y']
