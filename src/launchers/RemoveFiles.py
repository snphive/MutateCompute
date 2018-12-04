#!/usr/bin/env python3
"""
Script for removing specific (unnecessary) files from particular output directories.

NOTE: When running this script with remote debugger, there is a recursion error that occurs which ultimately interupts the
program prematurely.
"""
import os
import sys
import glob
from os import walk
from src.enums.Paths import Paths
from src.enums.Conditions import Cond
from src.tools.GeneralUtilityMethods import GUM
from src.FoldX import FoldX
from src.tools.OutputsParser import Parser
import pydevd
pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

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
Select Agadir-related files for deletion.
"""
# path_output_agadir_globresidfiles = Paths.OUTPUT_AGADIR + '/1...250/**/' + Str.GLOBRESIDUE_OUT.value
# path_globalresidueout_files = sorted(glob.glob(path_output_agadir_globresidfiles, recursive=True))
# if not path_globalresidueout_files:
#     raise ValueError('No PSX_globalresidue.out files to delete. Check paths are correct and check files are where
#       you expect.')
"""
Delete Agadir-related files.
"""
# for path_globalresidueout_file in path_globalresidueout_files:
#     # GUM.linux_remove_file(path_globalresidueout_file)
#     if not os.path.exists(path_globalresidueout_file):
#         print('Path: ' + path_globalresidueout_file + ' does not exist.')
#     else:
#         cmd = 'rm' + Str.SPCE.value + path_globalresidueout_file
#         try:
#             subprocess.call(cmd, shell=True)
#         except FileNotFoundError as fnf:
#             print('File not found: ' + path_globalresidueout_file)
#             print(fnf)
#         except OSError:
#             print(Str.PROBLNXCMD_MSG.value + cmd)

"""
Select pdb file(s).
"""
pdbnames = ['RepairPDB_14']

"""
Select specific mutants if you are only interested in these.
BE SURE TO SET THIS TO EMPTY LIST IF YOU DON'T WANT ANY OF THE SUBSEQUENT ACTIONS BELOW TO BE SPECIFIC TO THIS/THESE MUTANTS ONLY.
"""
# specific_fxmutants = ['RB186Q']
specific_fxmutants = []

"""
Decide which code blocks to run.
"""
delete_from_bm_outputs = False
delete_from_ac_outputs = True

"""
Define module-level function (before using it below).
"""
def remove_config_pdb_log_and_fxoutfiles(path_output_ac_or_bm):
    """
    Delete FoldX-related files - either from BuildModel or AnalyseComplex output folders.
    """
    path_output_ac_or_bm_pdb_fxmutantnames_dirs = []
    for pdbname in pdbnames:
        if specific_fxmutants:
            for specific_fxmutant in specific_fxmutants:
                path_output_ac_or_bm_pdb_fxmutantnames_dirs.append(os.path.join(path_output_ac_or_bm, pdbname, specific_fxmutant))
        else:
            path_output_ac_or_bm_pdb_fxmutantnames_dirs = glob.glob(os.path.join(path_output_ac_or_bm, pdbname, '*'))
        for path_output_ac_or_bm_pdb_fxmutantnames_dir in path_output_ac_or_bm_pdb_fxmutantnames_dirs:
            # fx.remove_config_files(path_output_ac_pdb_fxmutantnames_dir)
            # fx.remove_pdbfiles(path_output_ac_pdb_fxmutantnames_dir)
            # fx.remove_cluster_logfiles(path_output_ac_pdb_fxmutantnames_dir)
            # fx.remove_unnecessary_foldxfiles(path_output_ac_pdb_fxmutantnames_dir)
            fx.remove_all_sumry_except_1_0(path_output_ac_or_bm_pdb_fxmutantnames_dir)


fx = FoldX()
if delete_from_bm_outputs:
    remove_config_pdb_log_and_fxoutfiles(Paths.OUTPUT_BM)
if delete_from_ac_outputs:
    remove_config_pdb_log_and_fxoutfiles(Paths.OUTPUT_AC)

pydevd.stoptrace()
