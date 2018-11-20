#!/usr/bin/env python3
"""
Script for launching Blastp runs.
"""

import sys
import os
import glob
from src.tools.Paths import Paths
from src.Scheduler import Scheduler
from src.enums.Str import Str
import logging

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
Set up logger. 
"""
logfile = os.path.join(Paths.LOGS, 'log' + Str.LOGEXT.value)
logging.basicConfig(filename=logfile, level=logging.DEBUG)

"""
Get the fasta files you want to run mutate_fasta or agadir on.
"""
startnum = 1
endnum = 1000
path_input_files = os.path.join(Paths.INPUT_29611_FASTAS_1000, str(startnum) + Str.DOTS3.value + str(endnum))
path_fastafiles = sorted(glob.glob(path_input_files + '/**/*.fasta', recursive=True))
if not path_fastafiles:
    raise ValueError('No fasta files to process. Check paths are correct and check files are where you expect.')

"""
Kick off the program(s) via Scheduler's method start_blast().
"""
Scheduler.start_blast(path_input_fastafiles=path_fastafiles, path_output=Paths.OUTPUT, write_idmaps_for_mysldb=True,
                      write_csv=True, write_xml=True, write_json=False)

# IdProt.map_seq_to_swsprt_acc_id_and_write_files(path_input_fastas_dir=path_input_fastas, path_output=Paths.OUTPUT,
#                                                 write_idmaps_for_mysldb=True, write_csv=True, write_xml=True,
#                                                 write_json=False)


# pydevd.stoptrace()
