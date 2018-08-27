from src.Main import Main
from src.Paths import Paths

path_repo_pdbs = Paths.REPO_PDBS + '_10'
path_repo_fastas = Paths.REPO_FASTAS + '_10'
main = Main(use_cluster=False, use_multithread=False, path_repo_pdbs=path_repo_pdbs, path_repo_fastas=path_repo_fastas)

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
#
# Agadir        configuration/agadir_config/Options.txt
#               input_data/<fastafilename>/                         output_data/<fastafilename>/agadir/
#               input_data/<fastafilename>/mutants/                 output_data/<fastafilename>/agadir/<mutantfilename>/
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
