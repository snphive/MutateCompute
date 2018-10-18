import sys
import os
import natsort
import glob
from src.Paths import Paths
from src.Scheduler import Scheduler
from src.Str import Str
from src.IdentifyProtein import IdProt
from src.GeneralUtilityMethods import GUM
import threading
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))

startnum = 1
endnum = 1000
path_input_files = os.path.join(Paths.INPUT_29611_FASTAS_1000, str(startnum) + Str.DOTS3.value + str(endnum))
path_fastafiles = sorted(glob.glob(path_input_files + '/**/*.fasta', recursive=True))

if not path_fastafiles:
    raise ValueError('No fasta files to process. Check paths are correct and check files are where you expect.')
Scheduler.start_blast(path_input_fastafiles=path_fastafiles, path_output=Paths.OUTPUT, write_idmaps_for_mysldb=True,
                      write_csv=True, write_xml=True, write_json=False)

# IdProt.map_seq_to_swsprt_acc_id_and_write_files(path_input_fastas_dir=path_input_fastas, path_output=Paths.OUTPUT,
#                                                 write_idmaps_for_mysldb=True, write_csv=True, write_xml=True,
#                                                 write_json=False)


#
# PROGRAM           INPUT FILES & DIRECTORIES                           OUTPUT FILES & DIRECTORIES
# -------------------------------------------------------------------------------------------------------------------
# IdentifyProtein   input_data/<fastafilename>/                         output_data/blastp/
#                                                                       output_data/blastp/<fastafilename>_idmaps/
# Runs locally. Problem with running it as it is (via Biopython on the cluster).

# pydevd.stoptrace()
