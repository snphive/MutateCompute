from src.Main import Main
from src.Paths import Paths
from src.GeneralUtilityMethods import GUM
import os
import time

# use_cluster = True if sys.argv[1] == 'use_cluster=True' else False
use_cluster = False
use_multithread = False
Paths.set_up_paths(use_cluster=use_cluster)
path_repo_pdbs = Paths.REPO_PDBS + '_10'
startnum = 1
endnum = 1000

for i in range(79):

    path_repo_fastas = os.path.join(Paths.REPO_FASTAS, 'fastas_1000', str(startnum) + '...' + str(endnum))
    globaloptions_lines = Main._read_global_options(Paths.CONFIG_GLOBAL_OPTIONS + '/global_options.txt')
    wanted_pdbfile_list = Main._build_filelist_for_analysis(globaloptions_lines, path_repo_pdbs)
    wanted_fastafile_list = Main._build_filelist_for_analysis(globaloptions_lines, path_repo_fastas)
    operations = Main._determine_which_operations_to_perform(globaloptions_lines)
    mutant_aa_list = Main._determine_residues_to_mutate_to(globaloptions_lines)
    path_dst = Paths.INPUT
    path_wanted_pdbfile_list = GUM.copy_files_from_repo_to_input_dirs(path_repo_pdbs, path_dst, wanted_pdbfile_list)
    path_wanted_fastafile_list = GUM.copy_files_from_repo_to_input_dirs(path_repo_fastas, path_dst, wanted_fastafile_list)
    path_input = Paths.INPUT
    path_output = Paths.OUTPUT
    main = Main(use_cluster, operations, use_multithread, path_input, path_output, path_wanted_pdbfile_list,
                path_wanted_fastafile_list, mutant_aa_list)
    if i == 79:
        break
    time.sleep(120)
    startnum += 1000
    endnum += 1000


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
