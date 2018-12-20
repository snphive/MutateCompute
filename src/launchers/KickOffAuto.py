#!/usr/bin/env python3
"""
Script for running programs on files in different folders, locally or on cluster. Typically to run a folder of 1000
input files through a particular algorithm. This would be instead of reading in more than 10000 files at a time, as the
latter seems to cause problems with iterating through such a large list in a for loop - (I experienced a delay on
looping through a list of more than 70000, whereby computations missed out about half of the inputs, thereby generating
only about half of the outputs).

To run this script from:
- local terminal: Enter command sh KickOff.sh
- local Pycharm IDE: Select "KickOffManPy" in Run/Edit Configurations/Python: Click green 'play' button icon.
- remote (cluster) terminal: Enter command sh KickOffZeus.sh in /switchlab/group/shazib/SnpEffect/bash directory.
- remote (cluster) Pycharm IDE debugger: Select "KickOffRemoteDebug" in Run/Edit Configurations/Python Remote Debug. Click green
'bug' button icon. Make sure the script has "import pydevd" and "pydevd.settrace('localhost', port=51234, stdoutToServer=True,
stderrToServer=True)" at top of file and "pydevd.stoptrace()" at end. Log in to cluster with "ssh -p 7788
shazib@zeus.psb.ugent.be -R 51234:localhost:51234". Enter command sh KickOffZeus.sh in /switchlab/group/shazib/SnpEffect/bash
directory.
"""
import os
import sys
import glob

from src.Main import Main
from src.enums.Paths import Paths
from src.enums.AminoAcids import AA
from src.enums.Str import Str

# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

"""
Set up paths.
"""
Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))

"""
Set the value(s) of whichever operation(s) you want to run to True in the operations dict. 
operations is passed to Main via its constructor.
"""
operations = {'do_mutate_fasta': False, 'do_agadir': True, 'do_foldx_repair': False, 'do_foldx_buildmodel': False,
              'do_foldx_stability': False, 'do_foldx_analysecomplex': False}

"""
Multithreading is not up and running yet (Nov18). Set to False for now.
"""
use_multithread = False

"""
Get the pdb files you want to run FoldX on.
"""
path_pdbfiles = []

"""
Get the fasta files you want to run mutate_fasta or agadir on.
Select range of folder names (which have a numeric component) from which fasta files are to be read. 
For each folder and list of fasta files, run your selected program(s) via Main's constructor.
"""
startnum = 1
endnum = 1000

for i in range(29):
    dir_3dots = str(startnum) + Str.DOTS3.value + str(endnum)
    path_input_fastas_3dots_dir = os.path.join(Paths.INPUT_MUTS_MULTIFASTAS_29611_1000, dir_3dots)
    path_fastafiles = sorted(glob.glob(path_input_fastas_3dots_dir + '/**/*.fasta', recursive=True))
    if not path_fastafiles:
        raise ValueError('No fasta files to process. Check paths are correct and check files are where you expect.')
    """
    Kick off the program(s) via the constructor or Main class.
    """
    main = Main(operations, use_multithread, Paths.INPUT, Paths.OUTPUT, path_pdbfiles, path_fastafiles,
                AA.LIST_ALL_20_AA.value, write_to_csv_dumpfile_after_each_mutant=False)
    if i == 0:
        break
    startnum += 1000
    endnum += 1000

# pydevd.stoptrace()

import natsort
import time
import multiprocessing as mp
from src.enums.Conditions import Cond
from src.tools.GeneralUtilityMethods import GUM
from src.Agadir import Agadir
from src.Cluster import Cluster
"""
The 9 lines immediately below are only necessary if you want to run the programs via the global options file.
path_repo_fastas = os.path.join(Paths.REPO_FASTAS, 'fastas_1000', str(startnum) + Str.DOTS3.value + str(endnum))
globaloptions_lines = Main._read_global_options(Paths.CONFIG_GLOBAL_OPTIONS + '/global_options.txt')
wanted_pdbfile_list = Main._build_filelist_for_analysis(globaloptions_lines, path_repo_pdbs)
wanted_fastafile_list = Main._build_filelist_for_analysis(globaloptions_lines, path_repo_fastas)
operations = Main._determine_which_operations_to_perform(globaloptions_lines)
mutant_aa_list = Main._determine_residues_to_mutate_to(globaloptions_lines)
path_dst = Paths.INPUT
path_wanted_pdbfile_list = GUM.copy_files_from_repo_to_input_dirs(path_repo_pdbs, path_dst, wanted_pdbfile_list)
path_wanted_fastafile_list = GUM.copy_files_from_repo_to_input_dirs(path_repo_fastas, path_dst, wanted_fastafile_list)
"""

def run_agadir_on_1000_fastas():
    path_output_mutants = os.path.join(Paths.OUTPUT, Paths.DIR_MUTANTS_FASTAS.value, str(startnum) + Str.DOTS3.value + str(
        endnum))
    path_to_fastas = path_output_mutants + '/**/*' + Str.FSTAEXT.value
    path_fastafile_list = natsort.natsorted(glob.glob(path_to_fastas, recursive=True))
    agadir = Agadir(Cond.INCELL_MAML.value)
    for path_fastafile in path_fastafile_list:
        time.sleep(1)
        if GUM.using_cluster():
            jobname = Paths.PREFIX_WRITE.value + path_fastafile.split('/')[-1]
            path_to_script = os.path.join(Paths.SRC, 'write_1fastafile_per_fasta_from_multifastafile_zeus.py')
            Cluster.write_job_q_bash(jobname, path_job_q_dir=Paths.CONFIG_JOBQ,
                                     python_script_with_paths=path_to_script + Str.SPCE.value + path_fastafile)
            Cluster.run_job_q(path_job_q_dir=Paths.CONFIG_JOBQ)
        else:
            GUM.write_1_fastafile_per_fasta_from_multifastafile(path_dst=Paths.INPUT, path_fastafile=path_fastafile)
