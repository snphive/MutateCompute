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
from src.enums.Str import Str
from src.enums.Conditions import Cond
from src.Agadir import Agadir
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
1. Set up paths.
"""
Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))

"""
2. Select Agadir-related files for deletion:
"""
path_output_agadir_3dots_dir = os.path.join(Paths.OUTPUT_AGADIR, '1...250')

"""
3. Select pdb file(s).
"""
# pdbnames = ['RepairPDB_1', 'RepairPDB_3', 'RepairPDB_4', 'RepairPDB_5', 'RepairPDB_6', 'RepairPDB_7', 'RepairPDB_8',
#             'RepairPDB_9', 'RepairPDB_10']
pdbnames = ['RepairPDB_1']

"""
4. Select specific mutants if you are only interested in these.
BE SURE TO SET THIS TO EMPTY LIST IF YOU DON'T WANT ANY OF THE SUBSEQUENT ACTIONS BELOW TO BE SPECIFIC TO THIS/THESE MUTANTS ONLY.
"""
# specific_fxmutants = ['RB186Q']
specific_fxmutants = []

"""
5. This is the function encapsulting all calls to FoldX file removal methods. It is called in no.7 according to flags set in no.6
(Note you need to define this function (before using it below) because it is a module-level function rather than an object or 
class method).
"""
def remove_config_pdb_log_and_fxoutfiles(path_output_ac_or_bm):
    """
    Delete FoldX-related files - either from BuildModel or AnalyseComplex output folders.
    """
    path_output_ac_or_bm_pdb_fxmutant_dirs = []
    for pdbname in pdbnames:
        if specific_fxmutants:
            for specific_fxmutant in specific_fxmutants:
                path_output_ac_or_bm_pdb_fxmutant_dirs.append(os.path.join(path_output_ac_or_bm, pdbname, specific_fxmutant))
        else:
            path_output_ac_or_bm_pdb_fxmutant_dirs = glob.glob(os.path.join(path_output_ac_or_bm, pdbname, '*'))
        for path_output_ac_or_bm_pdb_fxmutant_dir in path_output_ac_or_bm_pdb_fxmutant_dirs:
            fx.remove_config_files(path_output_ac_or_bm_pdb_fxmutant_dir)
            fx.remove_pdbfiles(path_output_ac_or_bm_pdb_fxmutant_dir)
            fx.remove_cluster_logfiles(path_output_ac_or_bm_pdb_fxmutant_dir)
            fx.remove_unnecessary_foldxfiles(path_output_ac_or_bm_pdb_fxmutant_dir)
            fx.AnalyseComplex(Cond.INCELL_MAML_FX.value).remove_all_sumry_except_1_0(path_output_ac_or_bm_pdb_fxmutant_dir)


"""
6. Decide which code blocks to run.
"""
delete_from_agadir_outputs = False
delete_from_bm_outputs = True
delete_from_ac_outputs = False

"""
7. Call the function above according to the flags set in 
"""
fx = FoldX()
if delete_from_bm_outputs:
    remove_config_pdb_log_and_fxoutfiles(Paths.OUTPUT_BM)
if delete_from_ac_outputs:
    remove_config_pdb_log_and_fxoutfiles(Paths.OUTPUT_AC)

"""
8. Remove unwanted Agadir-related output files if flag set to true. 
"""
agadir = Agadir(Cond.INCELL_MAML_AG.value)
if delete_from_agadir_outputs:
    agadir.remove_unwanted_agad_output_files(path_output_agadir_3dots_dir)


pydevd.stoptrace()
