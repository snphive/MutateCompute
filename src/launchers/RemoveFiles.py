#!/usr/bin/env python3
"""
Script for removing specific (unnecessary) files from particular output directories.

NOTE: When running this script with remote debugger, there is a recursion error that occurs which ultimately interupts the
program prematurely.
"""
import os
import sys
import glob
import warnings
from src.enums.Paths import Paths
from src.enums.Conditions import Cond
from src.tools.GeneralUtilityMethods import GUM
from src.enums.Str import Str
from src.Agadir import Agadir
from src.Cluster import Cluster
from src.FoldX import FoldX
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
1. Set up paths.
"""
Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))

"""
2. Select Agadir-related files for deletion:
"""
path_output_agadir_3dots_dir = ''
# path_output_agadir_3dots_dir = os.path.join(Paths.OUTPUT_AGADIR, '1...250')

"""
3. Select pdb file(s).
"""
# pdbnames = ['RepairPDB_1', 'RepairPDB_3', 'RepairPDB_4', 'RepairPDB_5', 'RepairPDB_6', 'RepairPDB_7', 'RepairPDB_8',
#             'RepairPDB_9', 'RepairPDB_10']
pdbnames = ['RepairPDB_1']
path_pdbfiles = []
for pdbname in pdbnames:
    path_pdbfiles.append(os.path.join(Paths.INPUT_PDBS, pdbname + Str.PDBEXT.value))
if not path_pdbfiles:
    warnings.warn_explicit(message="No pdb files to process. Check paths are correct and check files are where you expect.",
                           category=RuntimeWarning, filename="KickOff", lineno=68)
"""
4. Select specific mutants if you are only interested in these.
BE SURE TO SET THIS TO EMPTY LIST IF YOU DON'T WANT ANY OF THE SUBSEQUENT ACTIONS BELOW TO BE SPECIFIC TO THIS/THESE MUTANTS ONLY.
"""
specific_fxmutants = ['CA498L']
# specific_fxmutants = []

"""
5. Select which algorithm's unwanted files to remove: 
"""
delete_from_buildmodel_outputs = True
delete_from_analysecomplex_outputs = False
delete_from_agadir_outputs = False

"""
6. Remove unwanted FoldX-related input/output files:
"""

using_cluster = GUM.using_cluster()
if delete_from_buildmodel_outputs:
    fx = FoldX()
    path_output_bm_pdb_fxmutant_dirs = []
    for path_pdbfile in path_pdbfiles:
        pdbfile = os.path.basename(path_pdbfile)
        pdbname = pdbfile.split('.')[0]
        for specific_fxmutant in specific_fxmutants:
            if specific_fxmutant[0] == specific_fxmutant[len(specific_fxmutant) - 1]:
                print('Mutations to wild-type residues are currently not being performed, therefore the results directory '
                      'should not exist for this: ' + str(pdbname) + '_' + specific_fxmutant)
                continue
            path_output_bm_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname, specific_fxmutant))
        if not specific_fxmutants:
            path_output_bm_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_BM, pdbname, '*'))
        for path_output_bm_pdb_fxmutant_dir in path_output_bm_pdb_fxmutant_dirs:
            os.chdir(path_output_bm_pdb_fxmutant_dir)
            bm = fx.BuildModel(Cond.INCELL_MAML_FX.value)
            if bm.has_already_generated_avg_bm_fxoutfile(path_output_bm_pdb_fxmutant_dir):
                fxmutantname = os.path.basename(path_output_bm_pdb_fxmutant_dir)
                if fx.has_already_removed_config_logs_fxoutfile(path_output_bm_pdb_fxmutant_dir):
                    print('No configs, logs, or unnecessary fxoutfiles found in ' + str(pdbname) + '_' + fxmutantname +
                          '. Hence, nothing to delete.')
                    continue
                if using_cluster:
                    path_jobq_dir = GUM.os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fxmutantname)
                    Cluster.write_job_q_bash(jobname=Paths.PREFIX_FX_RM.value + fxmutantname, path_job_q_dir=path_jobq_dir,
                                             python_script_with_paths=os.path.join(Paths.SE_SRC_CLSTR_PYSCRPTS.value,
                                                    'remove_files_zeus.py' + Str.SPCE.value + path_output_bm_pdb_fxmutant_dir))
                    Cluster.run_job_q(path_job_q_dir=path_jobq_dir)
                else:
                    fx.rm_config_files(path_output_bm_pdb_fxmutant_dir)
                    fx.rm_cluster_logfiles(path_output_bm_pdb_fxmutant_dir, rm_non_empty_err_files=True)
                    fx.rm_unnecessary_fxoutfiles(path_output_bm_pdb_fxmutant_dir)

if delete_from_analysecomplex_outputs:
    fx = FoldX()
    path_output_bm_pdb_fxmutant_dirs = []
    path_output_ac_pdb_fxmutant_dirs = []
    for path_pdbfile in path_pdbfiles:
        pdbfile = os.path.basename(path_pdbfile)
        pdbname = pdbfile.split('.')[0]
        for specific_fxmutant in specific_fxmutants:
            if specific_fxmutant[0] == specific_fxmutant[len(specific_fxmutant) - 1]:
                print('Computations on for wild-type residues are currently not being performed, therefore the results '
                      'directory should not exist for this: ' + str(pdbname) + '_' + specific_fxmutant)
                continue
            path_output_bm_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname, specific_fxmutant))
            path_output_ac_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_AC, pdbname, specific_fxmutant))
        if not specific_fxmutants:
            path_output_ac_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_AC, pdbname, '*'))
        for path_output_ac_pdb_fxmutant_dir in path_output_ac_pdb_fxmutant_dirs:
            os.chdir(path_output_ac_pdb_fxmutant_dir)
            ac = fx.AnalyseComplex(Cond.INCELL_MAML_FX.value)
            if ac.has_already_generated_summary_ac_fxoutfile(path_output_ac_pdb_fxmutant_dir):
                fxmutantname = os.path.basename(path_output_ac_pdb_fxmutant_dir)
                if fx.has_already_removed_config_logs_fxoutfile(path_output_ac_pdb_fxmutant_dir):
                    print('No configs, logs, or unnecessary fxoutfiles found in ' + str(pdbname) + '_' + fxmutantname +
                          '. Hence, nothing to delete.')
                if using_cluster:
                    path_jobq_dir = GUM.os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fxmutantname)
                    Cluster.write_job_q_bash(jobname=Paths.PREFIX_FX_RM.value + fxmutantname, path_job_q_dir=path_jobq_dir,
                                             python_script_with_paths=os.path.join(Paths.SE_SRC_CLSTR_PYSCRPTS.value,
                                                    'remove_files_zeus.py' + Str.SPCE.value + path_output_ac_pdb_fxmutant_dir))
                    Cluster.run_job_q(path_job_q_dir=path_jobq_dir)
                else:
                    fx.rm_config_files(path_output_ac_pdb_fxmutant_dir)
                    fx.rm_pdbfiles(os.path.join(Paths.OUTPUT_BM, pdbname, fxmutantname))
                    fx.rm_pdbfiles(os.path.join(path_output_ac_pdb_fxmutant_dir))
                    fx.rm_cluster_logfiles(path_output_ac_pdb_fxmutant_dir, rm_non_empty_err_files=True)
                    fx.rm_unnecessary_fxoutfiles(path_output_ac_pdb_fxmutant_dir)
                    ac.rm_all_sumry_except_1_0(path_output_ac_pdb_fxmutant_dir)


"""
7. Remove unwanted Agadir-related output files. 
"""
# if delete_from_agadir_outputs:
#     agadir = Agadir(Cond.INCELL_MAML_AG.value)
#     agadir.remove_unwanted_agad_output_files(path_output_agadir_3dots_dir)


# pydevd.stoptrace()
