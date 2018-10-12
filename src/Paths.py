import os
import yaml
from enum import Enum
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)


# All paths are constants and stored in the configuration directory in 'pathsAndDictionaries' yaml file.
# They are read into the Paths enum from which all other classes source all path strings.
#
# All paths are absolute paths unless prefix with an "REL_" which indicates they are relative paths.
class Paths(Enum):

    cwd = os.getcwd()
    YAML_FILE = '/Users/u0120577/PycharmProjects/MutateCompute/configuration/pathsAndDictionaries.yaml' if not \
        ('switchlab' in cwd or 'shazib' in cwd) else \
        '/switchlab/group/shazib/SnpEffect/configuration/pathsAndDictionaries.yaml'

    with open(YAML_FILE, 'r') as stream:

        try:

            paths_and_dictionaries = yaml.load(stream)

            # ABSOLUTE PATHS - EXECUTABLES IN LOCAL
            LOCAL_R_EXE = paths_and_dictionaries['ROOT']['path_local_r_exe']
            LOCAL_FOLDX_EXE = paths_and_dictionaries['ROOT']['path_local_fx_exe']
            LOCAL_AGADIR_EXE = paths_and_dictionaries['ROOT']['path_local_agad_exe']
            # ABSOLUTE PATHS - EXECUTABLES IN ZEUS CLUSTER
            ZEUS_R_EXE = paths_and_dictionaries['ROOT']['path_zeus_r_exe']
            ZEUS_FOLDX_EXE = paths_and_dictionaries['ROOT']['path_zeus_fx_exe']
            ZEUS_AGADIR_EXE = paths_and_dictionaries['ROOT']['path_zeus_agad_exe']
            ZEUS_QSUB_EXE = paths_and_dictionaries['ROOT']['path_zeus_qsub_exe']

            # ABSOLUTE PATHS - ROOTS LOCAL AND CLUSTER
            LOCAL_MUTATECOMPUTE = paths_and_dictionaries['ROOT']['path_local_mutatecompute']
            LOCAL_IO_ROOT = paths_and_dictionaries['ROOT']['path_local_io_root']
            ZEUS_SNPEFFECT = paths_and_dictionaries['ROOT']['path_zeus_snpeffect']

            #  DIRECTORY NAMES (USED TO CONSTRUCT ABSOLUTE PATHS BELOW)
            DIR_CONFIG = paths_and_dictionaries['ROOT']['dir_config']
            DIR_MUTFASTA_JOBQ = paths_and_dictionaries['ROOT']['dir_mf_jobq']
            DIR_AGAD_CONFIG = paths_and_dictionaries['ROOT']['dir_agadconfig']
            DIR_AGAD_JOBQ = paths_and_dictionaries['ROOT']['dir_agadir_jobq']
            DIR_BM_JOBQ = paths_and_dictionaries['ROOT']['dir_bm_jobq']
            DIR_AC_JOBQ = paths_and_dictionaries['ROOT']['dir_ac_jobq']
            DIR_BLST_JOBQ = paths_and_dictionaries['ROOT']['dir_blst_jobq']
            DIR_GLOBAL_OPTIONS = paths_and_dictionaries['ROOT']['dir_globoptions']
            DIR_FXCONFIG = paths_and_dictionaries['ROOT']['dir_fxconfig']
            DIR_ACRUNSCRIPT = paths_and_dictionaries['ROOT']['dir_acrunscript']
            DIR_BMRUNSCRIPT = paths_and_dictionaries['ROOT']['dir_bmrunscript']
            DIR_REPRUNSCRIPT = paths_and_dictionaries['ROOT']['dir_reprunscript']
            DIR_SRC = paths_and_dictionaries['ROOT']['dir_src']
            DIR_INPUT = paths_and_dictionaries['ROOT']['dir_input']
            DIR_MUTANTS = paths_and_dictionaries['ROOT']['dir_mutants']
            DIR_MUTANTS_FASTAS = paths_and_dictionaries['ROOT']['dir_mutants_fastas']
            DIR_MUTANTS_MULTIFASTAS = paths_and_dictionaries['ROOT']['dir_mutants_multifastas']
            DIR_FOLDX = paths_and_dictionaries['ROOT']['dir_fx']
            DIR_AGADIR = paths_and_dictionaries['ROOT']['dir_agadir']
            DIR_OUTPUT = paths_and_dictionaries['ROOT']['dir_output']
            DIR_AC = paths_and_dictionaries['ROOT']['dir_ac']
            DIR_BM = paths_and_dictionaries['ROOT']['dir_bm']
            DIR_REP = paths_and_dictionaries['ROOT']['dir_rep']
            DIR_FASTAS = paths_and_dictionaries['ROOT']['dir_fastas']
            DIR_29611_FASTAS_1000 = paths_and_dictionaries['ROOT']['dir_29611_fastas_1000']
            DIR_SEQS_TXT_CSV = paths_and_dictionaries['ROOT']['dir_seqs_txt_csv']
            DIR_REPO_PDB_FASTA = paths_and_dictionaries['ROOT']['dir_pdb_fasta_repo']
            DIR_PDBS = paths_and_dictionaries['ROOT']['dir_pdbs']
            DIR_BLASTP = paths_and_dictionaries['ROOT']['dir_blastp']
            DIR_BASH = paths_and_dictionaries['ROOT']['dir_bash']
            #  CLUSTER JOB NAME PREFIXES
            PREFIX_MUTFSTA = paths_and_dictionaries['ROOT']['prefix_mutatefasta']
            PREFIX_BLSTP = paths_and_dictionaries['ROOT']['prefix_blastp']
            PREFIX_AGADIR = paths_and_dictionaries['ROOT']['prefix_agadir']
            PREFIX_FX_BM = paths_and_dictionaries['ROOT']['prefix_foldx_buildmodel']
            PREFIX_FX_AC = paths_and_dictionaries['ROOT']['prefix_foldx_analysecomplex']

            # ABSOLUTE PATHS BUILT FROM LOCAL ROOT AND DIRECTORY NAMES
            MC_CONFIG = os.path.join(LOCAL_MUTATECOMPUTE, DIR_CONFIG)
            MC_CONFIG_AGAD = os.path.join(MC_CONFIG, DIR_AGAD_CONFIG)
            MC_CONFIG_FX = os.path.join(MC_CONFIG, DIR_FXCONFIG)
            MC_CONFIG_GLOBAL_OPTIONS = os.path.join(MC_CONFIG, DIR_GLOBAL_OPTIONS)
            MC_CONFIG_ACRUNSCRIPT = os.path.join(MC_CONFIG, DIR_ACRUNSCRIPT)
            MC_CONFIG_BMRUNSCRIPT = os.path.join(MC_CONFIG, DIR_BMRUNSCRIPT)
            MC_INPUT = os.path.join(LOCAL_MUTATECOMPUTE, DIR_INPUT)
            MC_INPUT_FASTAS = os.path.join(MC_INPUT, DIR_FASTAS)
            MC_OUTPUT = os.path.join(LOCAL_MUTATECOMPUTE, DIR_OUTPUT)
            IO_OUTPUT = os.path.join(LOCAL_IO_ROOT, DIR_OUTPUT)
            IO_INPUT = os.path.join(LOCAL_IO_ROOT, DIR_INPUT)
            MC_BASH = os.path.join(LOCAL_MUTATECOMPUTE, DIR_BASH)
            MC_SRC = os.path.join(LOCAL_MUTATECOMPUTE, DIR_SRC)
            MC_OUTPUT_FASTAS = os.path.join(MC_OUTPUT, DIR_FASTAS)
            MC_OUTPUT_AGADIR = os.path.join(MC_OUTPUT, DIR_AGADIR)
            IO_OUTPUT_FASTAS = os.path.join(IO_OUTPUT, DIR_FASTAS)
            MC_OUTPUT_BLASTP = os.path.join(MC_OUTPUT, DIR_BLASTP)
            IO_OUTPUT_BLASTP = os.path.join(IO_OUTPUT, DIR_BLASTP)
            IO_OUTPUT_AGADIR = os.path.join(IO_OUTPUT, DIR_AGADIR)
            # ABSOLUTE PATHS TO FASTA & PDB REPOSITORIES BUILT FROM LOCAL ROOT AND DIRECTORY NAMES
            MC_REPO_PDB_FASTA = paths_and_dictionaries['ROOT']['path_local_pdb_fasta_repository']
            MC_REPO_PDBS = os.path.join(MC_REPO_PDB_FASTA, DIR_PDBS)
            MC_REPO_FASTAS = os.path.join(MC_REPO_PDB_FASTA, DIR_FASTAS)

            # ABSOLUTE PATHS BUILT FROM ZEUS ROOT AND DIRECTORY NAMES
            SE_CONFIG = os.path.join(ZEUS_SNPEFFECT, DIR_CONFIG)
            SE_CONFIG_MUTFASTA_JOBQ = os.path.join(SE_CONFIG, DIR_MUTFASTA_JOBQ)
            SE_CONFIG_AGAD = os.path.join(SE_CONFIG, DIR_AGAD_CONFIG)
            SE_CONFIG_AGAD_JOBQ = os.path.join(SE_CONFIG, DIR_AGAD_JOBQ)
            SE_CONFIG_FX = os.path.join(SE_CONFIG, DIR_FXCONFIG)
            SE_CONFIG_GLOBAL_OPTIONS = os.path.join(SE_CONFIG, DIR_GLOBAL_OPTIONS)
            SE_CONFIG_ACRUNSCRIPT = os.path.join(SE_CONFIG, DIR_ACRUNSCRIPT)
            SE_CONFIG_AC_JOBQ = os.path.join(SE_CONFIG, DIR_AC_JOBQ)
            SE_CONFIG_BMRUNSCRIPT = os.path.join(SE_CONFIG, DIR_BMRUNSCRIPT)
            SE_CONFIG_BM_JOBQ = os.path.join(SE_CONFIG, DIR_BM_JOBQ)
            SE_CONFIG_BLST_JOBQ = os.path.join(SE_CONFIG, DIR_BLST_JOBQ)
            SE_INPUT = os.path.join(ZEUS_SNPEFFECT, DIR_INPUT)
            SE_INPUT_FASTAS = os.path.join(SE_INPUT, DIR_FASTAS)
            SE_OUTPUT = os.path.join(ZEUS_SNPEFFECT, DIR_OUTPUT)
            SE_BASH = os.path.join(ZEUS_SNPEFFECT, DIR_BASH)
            SE_SRC = os.path.join(ZEUS_SNPEFFECT, DIR_SRC)
            SE_OUTPUT_FASTAS = os.path.join(SE_OUTPUT, DIR_FASTAS)
            SE_OUTPUT_BLASTP = os.path.join(SE_OUTPUT, DIR_BLASTP)
            SE_OUTPUT_AGADIR = os.path.join(SE_OUTPUT, DIR_AGADIR)
            # ABSOLUTE PATHS TO FASTA & PDB REPOSITORIES BUILT FROM ZEUS ROOT AND DIRECTORY NAMES
            SE_REPO_PDB_FASTA = os.path.join(ZEUS_SNPEFFECT, DIR_REPO_PDB_FASTA)
            SE_REPO_PDBS = os.path.join(SE_REPO_PDB_FASTA, DIR_PDBS)
            SE_REPO_FASTAS = os.path.join(SE_REPO_PDB_FASTA, DIR_FASTAS)

        except yaml.YAMLError as exc:
            print(exc)

    # Will set absolute paths to zeus cluster if wanted, otherwise left as default paths which are local.
    #
    # use_cluster   Boolean     True if specified operations are to run on cluster.
    @staticmethod
    def set_up_paths(use_cluster):
        Paths.CONFIG = Paths.SE_CONFIG.value if use_cluster else Paths.MC_CONFIG.value
        Paths.CONFIG_AGAD = Paths.SE_CONFIG_AGAD.value if use_cluster else Paths.MC_CONFIG_AGAD.value
        Paths.CONFIG_AGAD_JOBQ = Paths.SE_CONFIG_AGAD_JOBQ.value
        Paths.CONFIG_BM_JOBQ = Paths.SE_CONFIG_BM_JOBQ.value
        Paths.CONFIG_FX = Paths.SE_CONFIG_FX.value if use_cluster else Paths.MC_CONFIG_FX.value
        Paths.CONFIG_GLOBAL_OPTIONS = Paths.SE_CONFIG_GLOBAL_OPTIONS.value if use_cluster else Paths.MC_CONFIG_GLOBAL_OPTIONS.value
        Paths.CONFIG_ACRUNSCRIPT = Paths.SE_CONFIG_ACRUNSCRIPT.value if use_cluster else Paths.MC_CONFIG_ACRUNSCRIPT.value
        Paths.CONFIG_AC_JOBQ = Paths.SE_CONFIG_AC_JOBQ.value
        Paths.CONFIG_BLST_JOBQ = Paths.SE_CONFIG_BLST_JOBQ.value
        Paths.CONFIG_BMRUNSCRIPT = Paths.SE_CONFIG_BMRUNSCRIPT.value if use_cluster else Paths.MC_CONFIG_BMRUNSCRIPT.value
        # Paths.INPUT = Paths.SE_INPUT.value if use_cluster else Paths.MC_INPUT.value
        Paths.INPUT = Paths.SE_INPUT.value if use_cluster else Paths.IO_INPUT.value
        Paths.INPUT_FASTAS = Paths.SE_INPUT_FASTAS.value if use_cluster else Paths.MC_INPUT_FASTAS.value
        # Paths.OUTPUT = Paths.SE_OUTPUT.value if use_cluster else Paths.MC_OUTPUT.value
        Paths.OUTPUT = Paths.SE_OUTPUT.value if use_cluster else Paths.IO_OUTPUT.value
        # Paths.OUTPUT_FASTAS = Paths.SE_OUTPUT_FASTAS.value if use_cluster else Paths.MC_OUTPUT_FASTAS.value
        Paths.OUTPUT_FASTAS = Paths.SE_OUTPUT_FASTAS.value if use_cluster else Paths.IO_OUTPUT_FASTAS.value
        # Paths.OUTPUT_BLASTP = Paths.SE_OUTPUT_BLASTP.value if use_cluster else Paths.MC_OUTPUT_BLASTP.value
        Paths.OUTPUT_BLASTP = Paths.SE_OUTPUT_BLASTP.value if use_cluster else Paths.IO_OUTPUT_BLASTP.value
        Paths.OUTPUT_AGADIR = Paths.SE_OUTPUT_AGADIR.value if use_cluster else Paths.IO_OUTPUT_AGADIR.value
        Paths.REPO_PDB_FASTA = Paths.SE_REPO_PDB_FASTA.value if use_cluster else Paths.MC_REPO_PDB_FASTA.value
        Paths.REPO_PDBS = Paths.SE_REPO_PDBS.value if use_cluster else Paths.MC_REPO_PDBS.value
        Paths.REPO_FASTAS = Paths.SE_REPO_FASTAS.value if use_cluster else Paths.MC_REPO_FASTAS.value
        Paths.R_EXE = Paths.ZEUS_R_EXE.value if use_cluster else Paths.LOCAL_R_EXE.value
        Paths.FOLDX_EXE = Paths.ZEUS_FOLDX_EXE.value if use_cluster else Paths.LOCAL_FOLDX_EXE.value
        Paths.AGADIR_EXE = Paths.ZEUS_AGADIR_EXE.value if use_cluster else Paths.LOCAL_AGADIR_EXE.value
        Paths.BASH = Paths.SE_BASH.value if use_cluster else Paths.MC_BASH.value
        Paths.SRC = Paths.SE_SRC.value if use_cluster else Paths.MC_SRC.value

# pydevd.stoptrace()
