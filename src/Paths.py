import os
import yaml
from enum import Enum


# All paths are constants and stored in the configuration directory in 'pathsAndDictionaries' yaml file.
# They are read into the Paths enum from which all other classes source all path strings.
#
# All paths are absolute paths unless prefix with an "REL_" which indicates they are relative paths.
class Paths(Enum):

    YAML_FILE = '/Users/u0120577/PycharmProjects/MutateCompute/configuration/pathsAndDictionaries.yaml'

    with open(YAML_FILE, 'r') as stream:

        try:

            paths_and_dictionaries = yaml.load(stream)

            # ABSOLUTE PATHS
            ZEUS_R_EXE = paths_and_dictionaries['ROOT']['path_zeus_r_exe']
            ZEUS_FOLDX_EXE = paths_and_dictionaries['ROOT']['path_zeus_fx_exe']
            ZEUS_AGADIR_EXE = paths_and_dictionaries['ROOT']['path_zeus_agad_exe']
            ZEUS_QSUB_EXE = paths_and_dictionaries['ROOT']['path_zeus_qsub_exe']
            LOCAL_R_EXE = paths_and_dictionaries['ROOT']['path_local_r_exe']
            LOCAL_FOLDX_EXE = paths_and_dictionaries['ROOT']['path_local_fx_exe']
            LOCAL_AGADIR_EXE = paths_and_dictionaries['ROOT']['path_local_agad_exe']

            ZEUS_SNPEFFECT = paths_and_dictionaries['ROOT']['path_zeus_snpeffect']
            LOCAL_MUTATECOMPUTE = paths_and_dictionaries['ROOT']['path_local_mutatecompute']

            LOCAL_REPO_PDB_FASTA = paths_and_dictionaries['ROOT']['path_local_pdb_fasta_repository']

            #  DIRECTORY NAMES (USED TO CONSTRUCT ABSOLYTE PATHS BELOW)
            DIR_CONFIG = paths_and_dictionaries['ROOT']['dir_config']
            DIR_AGAD = paths_and_dictionaries['ROOT']['dir_agadconfig']
            DIR_JOBQ = paths_and_dictionaries['ROOT']['dir_clustjobq']
            DIR_GLOBAL_OPTIONS = paths_and_dictionaries['ROOT']['dir_globoptions']
            DIR_FXCONFIG = paths_and_dictionaries['ROOT']['dir_fxconfig']
            DIR_ACRUNSCRIPT = paths_and_dictionaries['ROOT']['dir_acrunscript']
            DIR_BMRUNSCRIPT = paths_and_dictionaries['ROOT']['dir_bmrunscript']
            DIR_REPRUNSCRIPT = paths_and_dictionaries['ROOT']['dir_reprunscript']
            DIR_SRC = paths_and_dictionaries['ROOT']['dir_src']
            DIR_INPUT = paths_and_dictionaries['ROOT']['dir_input']
            DIR_MUTANTS = paths_and_dictionaries['ROOT']['dir_mutants']
            DIR_OUTPUT = paths_and_dictionaries['ROOT']['dir_output']
            DIR_AC = paths_and_dictionaries['ROOT']['dir_ac']
            DIR_BM = paths_and_dictionaries['ROOT']['dir_bm']
            DIR_REP = paths_and_dictionaries['ROOT']['dir_rep']
            DIR_FASTAS = paths_and_dictionaries['ROOT']['dir_fastas']
            DIR_SEQUENCES = paths_and_dictionaries['ROOT']['dir_sequences']
            DIR_REPO_PDB_FASTA = paths_and_dictionaries['ROOT']['dir_pdb_fasta_repo']
            DIR_PDBS = paths_and_dictionaries['ROOT']['dir_pdbs']
            DIR_BLASTP = paths_and_dictionaries['ROOT']['dir_blastp']

            # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
            MC_CONFIG = os.path.join(LOCAL_MUTATECOMPUTE, DIR_CONFIG)
            MC_CONFIG_AGAD = os.path.join(MC_CONFIG, DIR_AGAD)
            MC_CONFIG_JOBQ = os.path.join(MC_CONFIG, DIR_JOBQ)
            MC_CONFIG_FX = os.path.join(MC_CONFIG, DIR_FXCONFIG)
            MC_CONFIG_GLOBAL_OPTIONS = os.path.join(MC_CONFIG, DIR_GLOBAL_OPTIONS)
            MC_CONFIG_FX_ACRUNSCRIPT = os.path.join(MC_CONFIG, DIR_ACRUNSCRIPT)
            MC_CONFIG_FX_BMRUNSCRIPT = os.path.join(MC_CONFIG, DIR_BMRUNSCRIPT)
            MC_INPUT = os.path.join(LOCAL_MUTATECOMPUTE, DIR_INPUT)
            MC_OUTPUT = os.path.join(LOCAL_MUTATECOMPUTE, DIR_OUTPUT)
            MC_OUTPUT_BLASTP = os.path.join(MC_OUTPUT, DIR_BLASTP)

            MC_REPO_PDBS = os.path.join(LOCAL_REPO_PDB_FASTA, DIR_PDBS)
            MC_REPO_FASTAS = os.path.join(LOCAL_REPO_PDB_FASTA, DIR_FASTAS)

            SE_CONFIG = os.path.join(ZEUS_SNPEFFECT, DIR_CONFIG)
            SE_CONFIG_AGAD = os.path.join(SE_CONFIG, DIR_AGAD)
            SE_CONFIG_JOBQ = os.path.join(SE_CONFIG, DIR_JOBQ)
            SE_CONFIG_FXCONFIG = os.path.join(SE_CONFIG, DIR_FXCONFIG)
            SE_CONFIG_GLOBAL_OPTIONS = os.path.join(SE_CONFIG, DIR_GLOBAL_OPTIONS)
            SE_CONFIG_FX_ACRUNSCRIPT = os.path.join(SE_CONFIG, DIR_ACRUNSCRIPT)
            SE_CONFIG_FX_BMRUNSCRIPT = os.path.join(SE_CONFIG, DIR_BMRUNSCRIPT)
            SE_INPUT = os.path.join(ZEUS_SNPEFFECT, DIR_INPUT)
            SE_OUTPUT = os.path.join(ZEUS_SNPEFFECT, DIR_OUTPUT)
            SE_REPO_PDB_FASTA = os.path.join(ZEUS_SNPEFFECT, DIR_REPO_PDB_FASTA)
            SE_REPO_PDBS = os.path.join(SE_REPO_PDB_FASTA, DIR_PDBS)
            SE_REPO_FASTAS = os.path.join(SE_REPO_PDB_FASTA, DIR_FASTAS)
            SE_OUTPUT_BLASTP = os.path.join(SE_OUTPUT, DIR_BLASTP)

        except yaml.YAMLError as exc:
            print(exc)

