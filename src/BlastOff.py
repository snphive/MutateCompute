import sys
import natsort
import glob
from src.Paths import Paths
from src.Scheduler import Scheduler
from src.IdentifyProtein import IdProt
from src.GeneralUtilityMethods import GUM
import pydevd
import threading
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

if len(sys.argv) < 2:
    use_cluster = False
else:
    use_cluster = True if sys.argv[1] == 'use_cluster=True' else False

Paths.set_up_paths(use_cluster=use_cluster)
path_input_fastas = IdProt._build_dir_tree_with_intermed_dir(path_root=Paths.INPUT,
                                                             intermed_dir=Paths.DIR_FASTAS.value, fastadir=None)
# path_repo = Paths.MC_REPO_FASTAS.value + '/fastas_1000_rest'
# wanted_file_list = GUM.copy_files_from_repo_to_input_dirs(path_repo_pdbs_or_fastas=path_repo,
#                                                           path_dst_dir=path_input_fastas, wanted_file_list=None)
# path_fastafiles = '/Users/u0120577/REPO_PDB_FASTA/fastas_10/5_C.fasta'

path_input_fastas_dir = Paths.MC_REPO_FASTAS.value + '/fastas_1000_rest'
# /Users/u0120577/REPO_PDB_FASTA/fastas/fastas_1000_rest
path_fastafiles = natsort.natsorted(glob.glob(path_input_fastas_dir + '/*.fasta'))

if not path_fastafiles:
    raise ValueError('No fasta files to process. Check paths are correct and check files are where you expect.')
Scheduler.start_blast(path_input_fastafiles=path_fastafiles, path_output=Paths.OUTPUT, write_idmaps_for_mysldb=True,
                      write_csv=True, write_xml=True, write_json=False)

# IdProt.map_seq_to_swsprt_acc_id_and_write_files(path_input_fastas_dir=path_input_fastas, path_output=Paths.OUTPUT,
#                                                 write_idmaps_for_mysldb=True, write_csv=True, write_xml=True,
#                                                 write_json=False)


# As with the MutateFasta, Agadir and FoldX programs, the start of operations begins with identifying which or how many
# fasta files are to be analysed. These must then be copied over from the repository directory (REPO_PDB_FASTA/FASTAs_),
# to the input_data directory.
#
# (And as with the other programs, the necessary destination directory tree for the input files being copied over will
# be created if not already present. This involves a subdirectory bearing the name of the file, which for fastafiles is
# typically '1_A' (unchanged from RvdK's naming scheme).)
#
# E.g. '1_A.fasta' copied from ~/REPO_PDB_FASTA/fastas_100/ to ~/PycharmProjects/MutateCompute/input_data/1_A/
#
#
# PROGRAM           INPUT FILES & DIRECTORIES                           OUTPUT FILES & DIRECTORIES
# -------------------------------------------------------------------------------------------------------------------
# IdentifyProtein   input_data/<fastafilename>/                         output_data/blastp/
#                                                                       output_data/blastp/<fastafilename>_idmaps/

# pydevd.stoptrace()
