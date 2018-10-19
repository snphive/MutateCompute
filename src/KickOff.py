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

"""
This module is for doing batch runs, locally or on cluster, typically running a folder of 1,000 input files through a
computation. 
This instead of reading in > 10,000 files at a time, as the latter seems to cause problems with iterating through such 
a large list in a for loop - (I experienced a delay on looping through a list of > 70,000, whereby the computations 
missed out about half of the inputs, generating only half of the outputs).   
"""

Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))
operations = {'do_mutate_fasta': False, 'do_agadir': True, 'do_foldx_repair': False, 'do_foldx_buildmodel': False,
              'do_foldx_stability': False, 'do_foldx_analysecomplex': False}
use_multithread = False
path_pdbfiles = []
startnum = 1
endnum = 1000

for i in range(29):
    dir_3dots = str(startnum) + Str.DOTS3.value + str(endnum)
    path_input_fastas_3dots_dir = os.path.join(Paths.INPUT_MUTS_MULTIFASTAS_29611_1000, dir_3dots)
    path_fastafiles = sorted(glob.glob(path_input_fastas_3dots_dir + '/**/*.fasta', recursive=True))
    if not path_fastafiles:
        raise ValueError('No fasta files to process. Check paths are correct and check files are where you expect.')
    main = Main(operations, use_multithread, Paths.INPUT, Paths.OUTPUT, path_pdbfiles, path_fastafiles,
                AA.LIST_ALL_20_AA.value)
    if i == 0:
        break
    startnum += 1000
    endnum += 1000


# PROGRAM       INPUT FILES & DIRECTORIES                           OUTPUT FILES & DIRECTORIES
# -------------------------------------------------------------------------------------------------------------------
# MutateFasta   input_data/<fastafilename>/                         output_data/<fastafilename>/mutants/
# Run on PSB cluster (and completed for pub_hg_fasta dataset in 29611_fastas_1000) and locally - both worked well
# ~6:30pm Thu 11th Oct 2018, and locally ~13:00 Fri 11th Oct 2018.
# latest commit#: 667a7ccf921dc3d700de3e7616cabffacd9c2de8

# Agadir runs on PSB cluster & locally (my Macbook). 13:35 10Oct18 commit#: 3f56bc782eecab941aeb29377541166f0c636c08
# On cluster it took 1507 seconds (25 mins) to complete ITERATING through 1...1000 and sending the jobs to the cluster,
# however it then took a further.. to complete all the computations.
# 1...1000 took up 14 Gigabytes..
# Agadir        configuration/agadir_config/Options.txt
#               input_data/mutants_multifastas/                     output_data/<fastafilename>/agadir/3dots_dir/
#               29611_fastas_1000/3dots_dir/<fastafilename>/        <fastafilename>/
#               mutants/fastafilename_mutants.fasta
#
#
# FoldX BM (version 3) runs on PSB cluster. 15:25 11Oct18 commit#: 6c2c266c918be564bb5f19edb113c4f43adb4538
# Does not run locally (have only tried FoldX version 4 locally)
# FoldX BM      configuration/foldx_config/commands_buildmodel.txt
#               configuration/foldx_config/options_buildmodel.txt
#               configuration/foldx_config/rotabase.txt
#               input_data/pdbs/<pdbfilename>/                      output_data/<pdbfilename>/fx_bm/mutants/
#               input_data/<pdbfilename>/mutants/
#
# FoldX AC      configuration/foldx_config/rotabase.txt
#               input_data/pdbs/<pdbfilename>/                      output_data/<pdbfilename>/fx_ac/
#               input_data/pdbs/<pdbfilename>/mutants/              output_data/<pdbfilename>/fx_ac/mutants/
#
#
# To start this script from cmd line, sh KickOff.sh

pydevd.stoptrace()

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
