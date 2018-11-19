#!/usr/bin/env python3
"""
Script for running any part of the codebase on relatively small numbers of files typically housed in a single folder.
Those to be run on multiple folders can be done through the KickOff.py module.

KickOffManual.py can be run locally or on cluster.
"""
import sys
import os
import glob
import natsort
from src.Main import Main
from src.Paths import Paths
from src.GeneralUtilityMethods import GUM
from src.AminoAcids import AA
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
Paths for the entire codebase are set accordingly. "use_cluster" is set to False by default.  
"""
Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))

if 'switchlab/group' in os.getcwd() and sys.argv[1].strip(' ') != 'use_cluster=True':
    raise ValueError('Program seems to be called from within /switchlab/group path. This is likely to be the psb '
                     'cluster path but the "use_cluster" parameter is not set to True. Please check you are starting '
                     'the program with bash/KickOffZeus.sh which includes the "use_cluster=True" parameter. \n'
                     'Alternatively, uncomment the Paths.set_up_paths(use_cluster=True) line below to manually set it' 
                     'here. (Edit True to False if you are running it locally and not on the cluster)')
# Paths.set_up_paths(use_cluster=True)

if 'switchlab/group' not in os.getcwd() and sys.argv[1].strip(' ') == 'use_cluster=True':
    raise ValueError('Program seems to be called from outside of /switchlab/group path. Hence it is likely you are not '
                     'running the program on the psb cluster. However the "use_cluster" parameter is set to True. '
                     'Please check you are either starting the program with bash/KickOff.sh which includes the '
                     '"use_cluster=False" parameter. \n'
                     'Alternatively, uncomment the Paths.set_up_paths(use_cluster=False) line below to manually set it' 
                     'here.')
# Paths.set_up_paths(use_cluster=False)

"""
Set the value(s) of whichever operation(s) you want to run to True in the operations dict. 
operations is passed to Main via its constructor.
"""
operations = {'do_mutate_fasta': False, 'do_agadir': False, 'do_foldx_repair': False, 'do_foldx_buildmodel': False,
              'do_foldx_stability': False, 'do_foldx_analysecomplex': False}

"""
Multithreading is not up and running yet (Nov18). Set to False for now.
"""
use_multithread = False

"""
Paths and other parameters for programs run outside of Main() and the algorithms. This is just to do a quick test-run
on some functionality.
"""
path_dir_files_to_pack = Paths.SE_OUTPUT_AGADIR + '/1...250/02064053-3f32-32e6-9660-aaaffc30db87'

"""
Get the fasta files you want to run mutate_fasta or agadir on.
"""
path_fastafiles = []
# path_input_fastas_dir = Paths.INPUT_MUTS_MULTIFASTAS_29611_1000 + '/1...250/'
# path_fastafiles = sorted(glob.glob(path_input_fastas_dir + '/**/*.fasta', recursive=True))
# path_fastafiles = sorted(glob.glob(path_input_fastas_dir + '/*.fasta'))
# if not path_fastafiles:
#     raise ValueError('No fasta files to process. Check paths are correct and check files are where you expect.')

"""
Get the pdb files you want to run FoldX on.
"""
path_pdbfiles = []
# path_input_pdbs_dir = '/switchlab/group/shazib/SnpEffect/output_data/analyse_complex'
# path_input_pdbs_dir = Paths.OUTPUT_AC
# path_pdbfiles = sorted(glob.glob(path_input_pdbs_dir + '/**/*.pdb', recursive=True))
# if not path_pdbfiles:
#     raise ValueError('No pdb files to process. Check paths are correct and check files are where you expect.')

"""
Kick off the program(s) via the constructor or Main class.
"""
main = Main(operations, use_multithread, Paths.INPUT, Paths.OUTPUT, path_pdbfiles, path_fastafiles,
            AA.LIST_ALL_20_AA.value)


# pydevd.stoptrace()

