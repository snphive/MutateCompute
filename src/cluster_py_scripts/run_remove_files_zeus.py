import os
import sys
from src.enums.Paths import Paths
from src.FoldX import FoldX
from src.enums.Conditions import Cond

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

path_output_ac_or_bm_pdb_fxmutant_dir = sys.argv[1]

"""
Delete files no longer needed from the output folders in order to reduce disk space usage (particularly important on the cluster).

:sys.argv[1] path_output_ac_pdb_fxmutant_dir: Absolute path to each pdb/mutant folder.
"""

ac_or_bm = path_output_ac_or_bm_pdb_fxmutant_dir.split('/')[-3]
pdbname = path_output_ac_or_bm_pdb_fxmutant_dir.split('/')[-2]
fxmutantname = path_output_ac_or_bm_pdb_fxmutant_dir.split('/')[-1]
print('run_remove_files_zeus.py ###################################')
fx = FoldX()
if ac_or_bm == Paths.DIR_AC.value:
    ac = fx.AnalyseComplex(Cond.INCELL_MAML_FX.value)
    print('remove pdb files from bm folders')
    fx.rm_pdbfiles(os.path.join(Paths.OUTPUT_BM, pdbname, fxmutantname))
    print('remove pdb files from ac folders')
    fx.rm_pdbfiles(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir))
    print('remove any unwanted extra summary files')
    ac.rm_all_sumry_except_1_0(path_output_ac_or_bm_pdb_fxmutant_dir)
elif ac_or_bm == Paths.DIR_BM.value:
    bm = fx.BuildModel(Cond.INCELL_MAML_FX.value)
print('remove config files')
fx.rm_config_files(path_output_ac_or_bm_pdb_fxmutant_dir)
print('remove logfiles from ac folders')
fx.rm_logfiles(path_output_ac_or_bm_pdb_fxmutant_dir)
print('remove unnecessary fxout files from ac folders')
fx.rm_unnecessary_fxoutfiles(path_output_ac_or_bm_pdb_fxmutant_dir)