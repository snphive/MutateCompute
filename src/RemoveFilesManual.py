from src.FoldX import FoldX
from src.Paths import Paths
from src.GeneralUtilityMethods import GUM
from src.Conditions import Cond
import os
from os import walk
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

fx_mutant_dirs = []
pdbname = 'RepairPDB_2'
for (dirpath, dirnames, filenames) in walk(os.path.join(Paths.SE_OUTPUT_BM.value, pdbname)):
    fx_mutant_dirs.extend(dirnames)
    break

for fx_mutant_dir in fx_mutant_dirs:
    bm = FoldX().BuildModel(Cond.INCELL_MAML_FX.value)
    path_output_pdbname_mutant = os.path.join(Paths.SE_OUTPUT_BM.value, pdbname, fx_mutant_dir)
    bm.remove_config_files(path_output_pdbname_mutant)
    path_pdbfile = os.path.join(path_output_pdbname_mutant, pdbname + '.pdb')
    if os.path.exists(path_pdbfile):
        GUM.linux_remove_file(path_pdbfile)

# pydevd.stoptrace()
