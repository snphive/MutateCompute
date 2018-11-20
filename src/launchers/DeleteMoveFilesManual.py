#!/usr/bin/env python3
"""
Script for removing specific (unnecessary) files from particular output directories.
"""
from src.enums.Paths import Paths
import os
import sys
from src.tools.OutputsParser import Parser
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
Select files for deletion.
"""
# path_output_agadir_globresidfiles = Paths.OUTPUT_AGADIR + '/1...250/**/' + Str.GLOBRESIDUE_OUT.value
# path_globalresidueout_files = sorted(glob.glob(path_output_agadir_globresidfiles, recursive=True))
# if not path_globalresidueout_files:
#     raise ValueError('No PSX_globalresidue.out files to delete. Check paths are correct and check files are where you '
#                      'expect.')

"""
Delete files.
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
Select directory for packing into tar.
"""

path_dir_txt_files_to_pack = os.path.join(Paths.OUTPUT_AGADIR, '02064053-3f32-32e6-9660-aaaffc30db87')

"""
Pack files in directory into tar.
"""
Parser().make_tarfile(path_dir_txt_files_to_pack)


# Running FoldX BuildModel and deleting config files immediately after running.
#
# from src.FoldX import FoldX
# from src.GeneralUtilityMethods import GUM
# from src.Conditions import Cond
# from os import walk
#
# fx_mutant_dirs = []
# pdbname = 'RepairPDB_2'
# for (dirpath, dirnames, filenames) in walk(os.path.join(Paths.SE_OUTPUT_BM.value, pdbname)):
#     fx_mutant_dirs.extend(dirnames)
#     break
#
# for fx_mutant_dir in fx_mutant_dirs:
#     bm = FoldX().BuildModel(Cond.INCELL_MAML_FX.value)
#     path_output_pdbname_mutant = os.path.join(Paths.SE_OUTPUT_BM.value, pdbname, fx_mutant_dir)
#     bm.remove_config_files(path_output_pdbname_mutant)
#     path_pdbfile = os.path.join(path_output_pdbname_mutant, pdbname + '.pdb')
#     if os.path.exists(path_pdbfile):
#         GUM.linux_remove_file(path_pdbfile)

# pydevd.stoptrace()
