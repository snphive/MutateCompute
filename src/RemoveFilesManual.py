from src.FoldX import FoldX
from src.Paths import Paths
import os
from os import walk
import pydevd
pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

fx_mutant_dirs = []
for (dirpath, dirnames, filenames) in walk(os.path.join(Paths.SE_OUTPUT_BM.value, 'RepairPDB_1')):
    fx_mutant_dirs.extend(dirnames)
    break

for fx_mutant_dir in fx_mutant_dirs:
    FoldX.BuildModel.remove_config_files(fx_mutant_dir)

pydevd.stoptrace()
