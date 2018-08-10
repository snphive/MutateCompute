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
            DIR_CONFIG_AGAD = paths_and_dictionaries['ROOT']['dir_agadconfig']
            DIR_CONFIG_JOBQ = paths_and_dictionaries['ROOT']['dir_clustjobq']
            DIR_CONFIG_GLOBAL_OPTIONS = paths_and_dictionaries['ROOT']['dir_globoptions']
            DIR_CONFIG_FX = paths_and_dictionaries['ROOT']['dir_fxconfig']
            DIR_CONFIG_FX_ACRUNSCRIPT = paths_and_dictionaries['ROOT']['dir_acrunscript']
            DIR_CONFIG_FX_BMRUNSCRIPT = paths_and_dictionaries['ROOT']['dir_bmrunscript']
            DIR_CONFIG_FX_REPRUNSCRIPT = paths_and_dictionaries['ROOT']['dir_reprunscript']
            DIR_SRC = paths_and_dictionaries['ROOT']['dir_src']
            DIR_INPUT = paths_and_dictionaries['ROOT']['dir_input']
            DIR_MUTANTS = paths_and_dictionaries['ROOT']['dir_mutants']
            DIR_OUTPUT = paths_and_dictionaries['ROOT']['dir_output']
            DIR_FX_AC = paths_and_dictionaries['ROOT']['dir_ac']
            DIR_FX_BM = paths_and_dictionaries['ROOT']['dir_bm']
            DIR_FX_REP = paths_and_dictionaries['ROOT']['dir_rep']
            DIR_FASTAS = paths_and_dictionaries['ROOT']['dir_fastas']
            DIR_SEQUENCES = paths_and_dictionaries['ROOT']['dir_sequences']
            DIR_REPO_PDB_FASTA = paths_and_dictionaries['ROOT']['dir_pdb_fasta_repo']
            DIR_PDBS = paths_and_dictionaries['ROOT']['dir_pdbs']
            DIR_BLASTP = paths_and_dictionaries['ROOT']['dir_blastp']

            # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
            MC_CONFIG = os.path.join(LOCAL_MUTATECOMPUTE.value, DIR_CONFIG.value)
            MC_CONFIG_AGAD = os.path.join(LOCAL_MUTATECOMPUTE.va, DIR_CONFIG_AGAD.value)
            MC_CONFIG_JOBQ = os.path.join(LOCAL_MUTATECOMPUTE.value, DIR_CONFIG_JOBQ.value)
            MC_CONFIG_FX = os.path.join(LOCAL_MUTATECOMPUTE.value, DIR_CONFIG_FX.value)
            MC_CONFIG_GLOBAL_OPTIONS = os.path.join(LOCAL_MUTATECOMPUTE.value, DIR_CONFIG_GLOBAL_OPTIONS.value)
            MC_CONFIG_FX_ACRUNSCRIPT = os.path.join(LOCAL_MUTATECOMPUTE.value, DIR_CONFIG_FX_ACRUNSCRIPT.value)
            MC_CONFIG_FX_BMRUNSCRIPT = os.path.join(LOCAL_MUTATECOMPUTE.value, DIR_CONFIG_FX_BMRUNSCRIPT.value)
            MC_INPUT = os.path.join(LOCAL_MUTATECOMPUTE.value, DIR_INPUT.value)
            MC_OUTPUT = os.path.join(LOCAL_MUTATECOMPUTE.value, DIR_OUTPUT.value)
            MC_OUTPUT_BLASTP = os.path.join(MC_OUTPUT.value, DIR_BLASTP.value)
            MC_REPO_PDBS = os.path.join(LOCAL_REPO_PDB_FASTA.value, DIR_PDBS.value)
            MC_REPO_FASTAS = os.path.join(LOCAL_REPO_PDB_FASTA.value, DIR_FASTAS.value)

            SE_CONFIG = os.path.join(ZEUS_SNPEFFECT.value, DIR_CONFIG.value)
            SE_CONFIG_AGAD = os.path.join(ZEUS_SNPEFFECT.value, DIR_CONFIG_AGAD.value)
            SE_CONFIG_JOBQ = os.path.join(ZEUS_SNPEFFECT.value, DIR_CONFIG_JOBQ.value)
            SE_CONFIG_FX = os.path.join(ZEUS_SNPEFFECT.value, DIR_CONFIG_FX.value)
            SE_CONFIG_GLOBAL_OPTIONS = os.path.join(ZEUS_SNPEFFECT.value, DIR_CONFIG_GLOBAL_OPTIONS.value)
            SE_CONFIG_FX_ACRUNSCRIPT = os.path.join(ZEUS_SNPEFFECT.value, DIR_CONFIG_FX_ACRUNSCRIPT.value)
            SE_CONFIG_FX_BMRUNSCRIPT = os.path.join(ZEUS_SNPEFFECT.value, DIR_CONFIG_FX_BMRUNSCRIPT.value)
            SE_INPUT = os.path.join(ZEUS_SNPEFFECT.value, DIR_INPUT.value)
            SE_OUTPUT = os.path.join(ZEUS_SNPEFFECT.value, DIR_OUTPUT.value)
            SE_REPO_PDB_FASTA = os.path.join(ZEUS_SNPEFFECT.value, DIR_REPO_PDB_FASTA.value)
            SE_REPO_PDBS = os.path.join(SE_REPO_PDB_FASTA.value, DIR_PDBS.value)
            SE_REPO_FASTAS = os.path.join(SE_REPO_PDB_FASTA.value, DIR_FASTAS.value)
            SE_OUTPUT_BLASTP = os.path.join(SE_OUTPUT.value, DIR_BLASTP.value)

        except yaml.YAMLError as exc:
            print(exc)

