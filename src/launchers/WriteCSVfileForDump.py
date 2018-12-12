#!/usr/bin/env python3
"""
Script for reading multiple csv output files for writing to a single csv file, (in order to facilitate mysql dumps).

WriteCSVfileForDump.py can be run locally or on cluster.
"""
import sys
import os
import glob
import re
from src.tools.GeneralUtilityMethods import GUM
from src.enums.Paths import Paths
from src.enums.Str import Str
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
"""
2. Get csv files you want to read and write to single file.
"""
using_cluster = GUM.using_cluster()
path_output_ac_or_bm_dir = ''
path_output_csvfiles = glob.glob(Paths.OUTPUT_BM + '/**/*' + Str.CSVEXT.value, recursive=True)
"""
(some private functions defined here)
"""


def _extract_pdbname(filename: str):
    regex = re.compile(r"RepairPDB_[0-9]+_")
    res = regex.search(filename)
    pdbnme = filename[res.start() : res.end() - 1]
    return pdbnme


def _extract_fxmutantname(filename: str):
    regex = re.compile(r"RepairPDB_[0-9]+_")
    res = regex.search(filename)
    regex2 = re.compile(r"RepairPDB_[0-9]+_[A-Y][A-Y][0-9]+[A-Z]")
    res2 = regex2.search(filename)
    fxmutname = filename[res.end():res2.end()]
    return fxmutname


"""
3. Read csv files and write to a single csvfile for mysql dump. 
"""
for path_output_csvfile in path_output_csvfiles:
    filename = os.path.basename(path_output_csvfile)
    pdbname = _extract_pdbname(filename)
    fxmutantname = _extract_fxmutantname(filename)
    if using_cluster:
        path_jobq_dir = GUM.os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fxmutantname)
        Cluster.write_job_q_bash(jobname=Paths.PREFIX_FX_RM.value + fxmutantname, path_job_q_dir=path_jobq_dir,
                                 python_script_with_paths=os.path.join(Paths.SE_SRC_CLSTR_PYSCRPTS.value,
                                 'run_write_1csvfile_from_csvPerMutantfiles_zeus.py' + Str.SPCE.value + path_output_csvfile +
                                 Str.SPCE.value + path_output_ac_or_bm_dir + Str.SPCE.value + pdbname))
        Cluster.run_job_q(path_job_q_dir=path_jobq_dir)
    else:
        GUM.write_1_csvfile_from_csv_per_mutants(path_output_csvfile, path_output_ac_or_bm_dir, pdbname)
