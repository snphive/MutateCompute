#!/usr/bin/env python3
"""
Script for running Agadir and/or FoldX on fasta/pdb files.

Includes functionality for reducing disk space usage by deleting input files after computation completed and deleting output
files that are not needed.

Will include functionality for read data from output files and write to summary files and/or MySQL database.

KickOff.py can be run locally or on cluster.
"""
import sys
import os
import glob
import warnings
from src.tools.GeneralUtilityMethods import GUM
from src.enums.Paths import Paths
from src.enums.AminoAcids import AA
from src.enums.Str import Str
from src.enums.Conditions import Cond
from src.Main import Main
from src.FoldX import FoldX
from src.Cluster import Cluster
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
1. Set up paths. ("use_cluster" is set to False by default.)  
"""
Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))
# Paths.set_up_paths('use_cluster=True')
# Paths.set_up_paths('use_cluster=False')

"""
2. Set the value(s) of whichever operation(s) you want to run to True in the operations dict, which passed to Main via its 
constructor.
"""
operations = {Str.OPER_RUN_MUT_FSTA.value: False, Str.OPER_RUN_AGDR.value: False, Str.OPER_RUN_FX_RPR.value: False,
              Str.OPER_RUN_FX_BM.value: True, Str.OPER_RUN_FX_STAB.value: False, Str.OPER_RUN_FX_AC.value: False}

"""
3. Choose multithreading.
(Multithreading is not up and running yet (Nov18). Set to False for now.)
"""
use_multithread = False

"""
4. Get the pdb files you want to run FoldX on.
"""
# path_input_pdbs_dir = '/switchlab/group/shazib/SnpEffect/output_data/analyse_complex'
# path_input_pdbs_dir = Paths.OUTPUT_AC
# path_pdbfiles = sorted(glob.glob(path_input_pdbs_dir + '/**/*.pdb', recursive=True))
# path_pdbfiles = []
# pdbnames = ['RepairPDB_1', 'RepairPDB_3', 'RepairPDB_4', 'RepairPDB_5', 'RepairPDB_6', 'RepairPDB_7', 'RepairPDB_8',
#             'RepairPDB_9', 'RepairPDB_10', 'RepairPDB_11']
pdbnames = ['RepairPDB_1']
path_pdbfiles = []
for pdbname in pdbnames:
    path_pdbfiles.append(os.path.join(Paths.INPUT_PDBS, pdbname + Str.PDBEXT.value))
if not path_pdbfiles:
    warnings.warn_explicit(message="No pdb files to process. (No problem if that is what you expected).",
                           category=RuntimeWarning, filename="KickOff", lineno=68)
"""
5. Select specific mutants if you are only interested in these.
MAKE SURE TO SET THIS LIST TO EMPTY if you don't want any of the subsequent actions below to be for these mutants only.
"""
# specific_fxmutants = ['CA498L']
specific_fxmutants = []
for path_pdbfile in path_pdbfiles:
    for specific_fxmutant in specific_fxmutants:
        if not GUM.is_valid_fxmutant_for_pdb(path_pdbfile, specific_fxmutant):
            raise ValueError('The specified mutant ' + specific_fxmutant + ' is not valid for this pdb: ' +
                             os.path.basename(path_pdbfile))

"""
6. Get the fasta files you want to run mutate_fasta or agadir on.
"""
path_fastafiles = []
# path_input_fastas_dir = Paths.INPUT_MUTS_MULTIFASTAS_29611_1000 + '/1...250/'
# path_fastafiles = sorted(glob.glob(path_input_fastas_dir + '/**/*.fasta', recursive=True))
# path_fastafiles = sorted(glob.glob(path_input_fastas_dir + '/*.fasta'))
if not path_fastafiles:
    warnings.warn_explicit(message="No fasta files to process. (No problem if that is what you expected).",
                           category=RuntimeWarning, filename="KickOff", lineno=84)
"""
7. Kick off the program(s) via the constructor or Main class. (Decide whether to write to csv after each mutant or after
all mutant computations are completed.)
"""
do_not_write_at_all = True
write_to_csv_dumpfile_after_each_mutant = False
write_to_csv_dumpfile_after_all_computations = False if write_to_csv_dumpfile_after_each_mutant else True

# main = Main(operations, use_multithread, Paths.INPUT, Paths.OUTPUT, path_pdbfiles, path_fastafiles, specific_fxmutants,
#             AA.LIST_ALL_20_AA.value, write_to_csv_dumpfile_after_each_mutant)

"""
8. After computation completed, WRITE to csv files.  
"""
if not do_not_write_at_all:
    if write_to_csv_dumpfile_after_all_computations:
        if operations[Str.OPER_RUN_FX_BM.value]:
            fx = FoldX()
            using_cluster = GUM.using_cluster()
            path_output_bm_pdb_fxmutant_dirs = []
            for path_pdbfile in path_pdbfiles:
                pdbfile = os.path.basename(path_pdbfile)
                pdbname = pdbfile.split('.')[0]
                for specific_fxmutant in specific_fxmutants:
                    if specific_fxmutant[0] == specific_fxmutant[len(specific_fxmutant) - 1]:
                        print('Mutations to wild-type residues are currently not being performed, therefore the results'
                              'directory should not exist for this: ' + str(pdbname) + '_' + specific_fxmutant)
                        continue
                    path_output_bm_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname, specific_fxmutant))
                if not specific_fxmutants:
                    path_output_bm_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_BM, pdbname, '*'))
                for path_output_bm_pdb_fxmutant_dir in path_output_bm_pdb_fxmutant_dirs:
                    fxmutantname = os.path.basename(path_output_bm_pdb_fxmutant_dir)
                    os.chdir(path_output_bm_pdb_fxmutant_dir)
                    bm = fx.BuildModel(Cond.INCELL_MAML_FX.value)
                    while not bm.has_already_generated_avg_bm_fxoutfile(path_output_bm_pdb_fxmutant_dir):
                        print(fx.Strs.AVG_BMDL_.value + pdbname + fx.Strs.FXOUTEXT.value +
                              ' has not been created yet for ' + fxmutantname +
                              '. Therefore cannot write to csv file just yet.')
                        if using_cluster:
                            jobname = Paths.PREFIX_FX_BM.value + fxmutantname
                            Cluster.wait_for_grid_engine_job_to_complete(jobname)
                    fx.write_to_csvfile_for_db_dump(path_output_bm_pdb_fxmutant_dir)
                    if using_cluster:
                        path_jobq_dir = GUM.os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fxmutantname)
                        Cluster.write_job_q_bash(jobname=Paths.PREFIX_WRITE.value + fxmutantname, path_job_q_dir=path_jobq_dir,
                                                 python_script_with_paths=os.path.join(Paths.SE_SRC_CLSTR_PYSCRPTS.value,
                                                 'write_csvdumpfile_from_fxout_zeus.py' + Str.SPCE.value +
                                                  path_output_bm_pdb_fxmutant_dir))
                        Cluster.run_job_q(path_job_q_dir=path_jobq_dir)
                    else:
                        fx.write_to_csvfile_for_db_dump(path_output_bm_pdb_fxmutant_dir)

    if write_to_csv_dumpfile_after_all_computations:
        if operations[Str.OPER_RUN_FX_AC.value]:
            fx = FoldX()
            using_cluster = GUM.using_cluster()
            path_output_ac_pdb_fxmutant_dirs = []
            for path_pdbfile in path_pdbfiles:
                pdbfile = os.path.basename(path_pdbfile)
                pdbname = pdbfile.split('.')[0]
                for specific_fxmutant in specific_fxmutants:
                    if specific_fxmutant[0] == specific_fxmutant[len(specific_fxmutant) - 1]:
                        print('Mutations to wild-type residues are currently not being performed, therefore the results'
                              'directory should not exist for this: ' + str(pdbname) + '_' + specific_fxmutant)
                        continue
                    path_output_ac_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_AC, pdbname, specific_fxmutant))
                if not specific_fxmutants:
                    path_output_ac_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_ac, pdbname, '*'))
                for path_output_ac_pdb_fxmutant_dir in path_output_ac_pdb_fxmutant_dirs:
                    fxmutantname = os.path.basename(path_output_ac_pdb_fxmutant_dir)
                    os.chdir(path_output_ac_pdb_fxmutant_dir)
                    ac = fx.AnalyseComplex(Cond.INCELL_MAML_FX.value)
                    while not ac.has_already_generated_summary_ac_fxoutfile(path_output_ac_pdb_fxmutant_dir):
                        print(fx.Strs.SMRY_AC_.value + pdbname + fx.Strs.FXOUTEXT.value + ' has not been created yet for ' +
                              fxmutantname + '. Therefore cannot write to csv file just yet.')
                        if using_cluster:
                            jobname = Paths.PREFIX_FX_AC.value + fxmutantname
                            Cluster.wait_for_grid_engine_job_to_complete(jobname)
                    fx.write_to_csvfile_for_db_dump(path_output_ac_pdb_fxmutant_dir)
                    if using_cluster:
                        path_jobq_dir = GUM.os_makedirs(Paths.CONFIG_AC_JOBQ, pdbname, fxmutantname)
                        Cluster.write_job_q_bash(jobname=Paths.PREFIX_WRITE.value + fxmutantname, path_job_q_dir=path_jobq_dir,
                                                 python_script_with_paths=os.path.join(Paths.SE_SRC_CLSTR_PYSCRPTS.value,
                                                 'write_csvdumpfile_from_fxout_zeus.py' + Str.SPCE.value +
                                                  path_output_ac_pdb_fxmutant_dir))
                        Cluster.run_job_q(path_job_q_dir=path_jobq_dir)
                    else:
                        fx.write_to_csvfile_for_db_dump(path_output_ac_pdb_fxmutant_dir)

"""
9. Decide whether to delete unwanted files at all:
"""
do_not_delete_at_all = False

"""
10.After computation completed, DELETE config files no longer needed.  
"""
if not do_not_delete_at_all:
    # Seems to be some bug(s) with removing log files using cluster jobs. Therefore forcing it not to send remove jobs
    # to cluster with flag use_cluster.
    use_cluster = False
    if operations[Str.OPER_RUN_FX_BM.value]:
        fx = FoldX()
        path_output_bm_pdb_fxmutant_dirs = []
        for path_pdbfile in path_pdbfiles:
            pdbfile = os.path.basename(path_pdbfile)
            pdbname = pdbfile.split('.')[0]
            for specific_fxmutant in specific_fxmutants:
                if specific_fxmutant[0] == specific_fxmutant[len(specific_fxmutant) - 1]:
                    print('Mutations to wild-type residues are currently not being performed, therefore the results'
                          'directory should not exist for this: ' + str(pdbname) + '_' + specific_fxmutant)
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
                    if use_cluster:
                        path_jobq_dir = GUM.os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fxmutantname)
                        Cluster.write_job_q_bash(jobname=Paths.PREFIX_FX_RM.value + fxmutantname, path_job_q_dir=path_jobq_dir,
                                                 python_script_with_paths=os.path.join(Paths.SE_SRC_CLSTR_PYSCRPTS.value,
                                                 'remove_files_zeus.py' + Str.SPCE.value + path_output_bm_pdb_fxmutant_dir))
                        Cluster.run_job_q(path_job_q_dir=path_jobq_dir)
                    else:
                        fx.rm_config_files(path_output_bm_pdb_fxmutant_dir)
                        fx.rm_cluster_logfiles(path_output_bm_pdb_fxmutant_dir, rm_non_empty_err_files=True)
                        fx.rm_unnecessary_fxoutfiles(path_output_bm_pdb_fxmutant_dir)

    if operations[Str.OPER_RUN_FX_AC.value]:
        fx = FoldX()
        path_output_bm_pdb_fxmutant_dirs = []
        path_output_ac_pdb_fxmutant_dirs = []
        for path_pdbfile in path_pdbfiles:
            pdbfile = os.path.basename(path_pdbfile)
            pdbname = pdbfile.split('.')[0]
            for specific_fxmutant in specific_fxmutants:
                if specific_fxmutant[0] == specific_fxmutant[len(specific_fxmutant) - 1]:
                    print('Computations on for wild-type residues are currently not being performed, therefore the '
                          'results directory should not exist for this: ' + str(pdbname) + '_' + specific_fxmutant)
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
                    if use_cluster:
                        path_jobq_dir = GUM.os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fxmutantname)
                        Cluster.write_job_q_bash(jobname=Paths.PREFIX_FX_RM.value + fxmutantname, path_job_q_dir=path_jobq_dir,
                                                 python_script_with_paths=os.path.join(Paths.SE_SRC_CLSTR_PYSCRPTS.value,
                                                 'remove_files_zeus.py' + Str.SPCE.value + path_output_ac_pdb_fxmutant_dir))
                        Cluster.run_job_q(path_job_q_dir=path_jobq_dir)
                    else:
                        fx.rm_config_files(path_output_ac_pdb_fxmutant_dir)
                        fx.rm_pdbfiles(os.path.join(Paths.OUTPUT_BM, pdbname, fxmutantname))
                        fx.rm_pdbfiles(os.path.join(path_output_ac_pdb_fxmutant_dir))
                        fx.rm_cluster_logfiles(path_output_ac_pdb_fxmutant_dir, rm_non_empty_err_files=False)
                        fx.rm_unnecessary_fxoutfiles(path_output_ac_pdb_fxmutant_dir)
                        ac.rm_all_sumry_except_1_0(path_output_ac_pdb_fxmutant_dir)


# pydevd.stoptrace()
