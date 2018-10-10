import os
import sys
import glob
import natsort
import time
from src.Main import Main
from src.Paths import Paths
from src.AminoAcids import AA
from src.Str import Str
from src.Agadir import Agadir
from src.Agadir import AgadCndtns
from src.Cluster import Cluster
from src.GeneralUtilityMethods import GUM
import multiprocessing as mp
import pydevd
pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

if len(sys.argv) < 2:
    use_cluster = False
else:
    use_cluster = True if sys.argv[1] == 'use_cluster=True' else False

Paths.set_up_paths(use_cluster)
# path_repo_pdbs = Paths.REPO_PDBS + '_10'

# NOTE: Ideally, if use_cluster is True, this script would set up an ssh tunnel to zeus and copy all src and config
# files to SnpEffect in zeus cluster before running the code therein - but this is not currently set up.
# You must manually copy over the scripts to the cluster, make an ssh tunnel and run the code on the cluster.

# for i in range(1):

# THE FOLLOWING ARE ONLY NECESSARY IF YOU WANT TO RUN THE WHOLE THING FROM VALUES IN THE GLOBAL OPTIONS TEXT FILE.
# path_repo_fastas = os.path.join(Paths.REPO_FASTAS, 'fastas_1000', str(startnum) + Str.DOTS3.value + str(endnum))
# globaloptions_lines = Main._read_global_options(Paths.CONFIG_GLOBAL_OPTIONS + '/global_options.txt')
# wanted_pdbfile_list = Main._build_filelist_for_analysis(globaloptions_lines, path_repo_pdbs)
# wanted_fastafile_list = Main._build_filelist_for_analysis(globaloptions_lines, path_repo_fastas)
# operations = Main._determine_which_operations_to_perform(globaloptions_lines)
# mutant_aa_list = Main._determine_residues_to_mutate_to(globaloptions_lines)
# path_dst = Paths.INPUT
# path_wanted_pdbfile_list = GUM.copy_files_from_repo_to_input_dirs(path_repo_pdbs, path_dst, wanted_pdbfile_list)
# path_wanted_fastafile_list = GUM.copy_files_from_repo_to_input_dirs(path_repo_fastas, path_dst, wanted_fastafile_list)

# path_input = Paths.INPUT
# path_output = Paths.OUTPUT

# startnum = 1
# endnum = 10
# dir_3dots = str(startnum) + Str.DOTS3.value + str(endnum)
# path_multifastas_3dots = os.path.join(Paths.INPUT, Paths.DIR_MUTANTS_MULTIFASTAS.value, dir_3dots)
# path_multifastafiles = path_multifastas_3dots + '/**/*' + Str.FSTAEXT.value
#
# operations = {'do_mutate_fasta': False, 'do_agadir': True, 'do_foldx_repair': False,'do_foldx_buildmodel': True,
#              'do_foldx_stability': False,'do_foldx_analysecomplex': False}
# use_multithread = False
# path_pdbfile_list = []
# path_fastafile_list = natsort.natsorted(glob.glob(path_multifastafiles, recursive=True))
# main = Main(operations, use_multithread, Paths.INPUT, Paths.OUTPUT, path_pdbfile_list, path_fastafile_list,
#             AA.LIST_ALL_20_AA.value)

# main = Main(operations, use_multithread, path_input, path_output, path_wanted_pdbfile_list,
# path_wanted_fastafile_list, mutant_aa_list)

#
# startnum = 1
# endnum = 10
# dir_3dots = str(startnum) + Str.DOTS3.value + str(endnum)
# path_pdbs_3dots = os.path.join(Paths.INPUT, Paths.DIR_PDBS.value, dir_3dots)
# path_pdbfiles = path_pdbs_3dots + '/**/*' + Str.PDBEXT.value

operations = {'do_mutate_fasta': False, 'do_agadir': False,
              'do_foldx_repair': False, 'do_foldx_buildmodel': True,
              'do_foldx_stability': False, 'do_foldx_analysecomplex': False}
use_multithread = False
# path_pdbfiles = natsort.natsorted(glob.glob(path_pdbfiles, recursive=True))
if use_cluster:
    path_pdbfiles = ['/switchlab/group/shazib/SnpEffect/tests/input_data/pdbs/RepairPDBtest_1.pdb']
else:
    path_pdbfiles = ['/Users/u0120577/PycharmProjects/MutateCompute/tests/input_data/pdbs/RepairPDBtest_1.pdb']
path_fastafiles = []
main = Main(operations, use_multithread, Paths.INPUT, Paths.OUTPUT, path_pdbfiles, path_fastafiles, ['A', 'C'])
            # AA.LIST_ALL_20_AA.value)

# if i == 50:
#     break
# time.sleep(10)
# startnum += 1000
# endnum += 1000

#
# def run_agadir_on_1000_fastas():
#     path_output_mutants = os.path.join(Paths.OUTPUT, Paths.DIR_MUTANTS_FASTAS.value,
# str(startnum) + Str.DOTS3.value + str(endnum))
#     path_to_fastas = path_output_mutants + '/**/*' + Str.FSTAEXT.value
#     path_fastafile_list = natsort.natsorted(glob.glob(path_to_fastas, recursive=True))
#     agadir = Agadir(AgadCndtns.INCELL_MAML.value)
#     for path_fastafile in path_fastafile_list:
#         time.sleep(1)
#         if use_cluster:
#             jobname = 'wr_' + path_fastafile.split('/')[-1]
#             path_to_script = os.path.join(Paths.SRC, 'run_write_1fastafile_per_fasta_from_multifastafile_zeus.py')
#             Cluster.write_job_q_bash(jobname, path_job_q_dir=Paths.CONFIG_JOBQ,
# python_script_with_paths=path_to_script + Str.SPCE.value + path_fastafile)
#             Cluster.run_job_q(path_job_q_dir=Paths.CONFIG_JOBQ)
#         else:
#             GUM.write_1_fastafile_per_fasta_from_multifastafile(path_dst=Paths.INPUT, path_fastafile=path_fastafile)

# In all cases, the start of operations begins with identifying which or how many pdb and/or fasta files are to be
# analysed. These must then be copied over from either the repository directory (REPO_PDB_FASTA) or the output_data
# directory, to the input_data directory.
# (In all cases, the necessary destination directory tree for the input files being copied over will be created if not
# already present. This involves a subdirectory bearing the name of the file, which for fastafiles is typically '1_A'
# and for pdbfiles is typically 'RepairPDB_1', (unchanged from RvdK's naming scheme)).
#
# E.g. '1_A.fasta' copied from ~/REPO_PDB_FASTA/FASTAs_100/ to ~/PycharmProjects/MutateCompute/input_data/1_A/
# 'RepairPDB_1' copied from ~/REPO_PDB_FASTA/PDBs_100/ to ~/PycharmProjects/MutateCompute/input_data/RepairPDB_1/
#
# The destination for output files is within the output_data directory tree.
# and FoldX BuildModel are creating mutants precisely for further analysis (by Agadir and FoldX)
#
#
# PROGRAM       INPUT FILES & DIRECTORIES                           OUTPUT FILES & DIRECTORIES
# -------------------------------------------------------------------------------------------------------------------
# MutateFasta   input_data/<fastafilename>/                         output_data/<fastafilename>/mutants/
# (completed calculations 05.09.18)
#
# Agadir        configuration/agadir_config/Options.txt
#               input_data/<fastafilename>/                         output_data/<fastafilename>/agadir/
#               input_data/<fastafilename>/mutants/                 output_data/<fastafilename>/agadir/<mutantfilename>/
# Tested on PSB cluster and locally on Macbook - both work. 13:35 10th Oct 2018 (Git revision number:
# 3f56bc782eecab941aeb29377541166f0c636c08)
#
# FoldX BM      configuration/foldx_config/commands_buildmodel.txt
#               configuration/foldx_config/options_buildmodel.txt
#               configuration/foldx_config/rotabase.txt
#               input_data/<pdbfilename>/                           output_data/<pdbfilename>/fx_bm/mutants/
#               input_data/<pdbfilename>/mutants/
#
# FoldX AC      configuration/foldx_config/rotabase.txt
#               input_data/<pdbfilename>/                           output_data/<pdbfilename>/fx_ac/
#               input_data/<pdbfilename>/mutants/                   output_data/<pdbfilename>/fx_ac/mutants/
#
#
# To start this script from cmd line, sh KickOff.sh

pydevd.stoptrace()
