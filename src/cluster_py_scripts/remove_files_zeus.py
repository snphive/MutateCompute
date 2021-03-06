#!/usr/bin/env python3
"""
Script for removing unnecessary FoldX input & output files from within cluster, specified in job.q.
"""
import os
import sys
from src.enums.Paths import Paths
from src.enums.Str import Str
from src.FoldX import FoldX
from src.tools.GeneralUtilityMethods import GUM
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
fx = FoldX()
if ac_or_bm == Paths.DIR_AC.value:
    print(ac_or_bm + 'it is in the ac scope')
    ac = fx.AnalyseComplex(Cond.INCELL_MAML_FX.value)
    fx.rm_pdbfiles(os.path.join(Paths.OUTPUT_BM, pdbname, fxmutantname))
    fx.rm_pdbfiles(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir))
    ac.rm_all_sumry_except_1_0(path_output_ac_or_bm_pdb_fxmutant_dir)
elif ac_or_bm == Paths.DIR_BM.value:
    print(ac_or_bm + 'it is in the bm scope')
    bm = fx.BuildModel(Cond.INCELL_MAML_FX.value)
    path_pdbfile = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, pdbname + Str.PDBEXT.value)
    if GUM.get_num_of_chains(path_pdbfile) == 1:
        # BuildModel has no use for these mutant pdb files and because this pdb only has 1 chain, AnalyseComplex has
        # no use for them either. Hence, they can already be removed.
        fx.rm_pdbfiles(path_output_ac_or_bm_pdb_fxmutant_dir)
else:
    print(ac_or_bm + ' something is not working !')
fx.rm_config_files(path_output_ac_or_bm_pdb_fxmutant_dir)
fx.rm_cluster_logfiles(path_output_ac_or_bm_pdb_fxmutant_dir, rm_non_empty_err_files=False)
fx.rm_unnecessary_fxoutfiles(path_output_ac_or_bm_pdb_fxmutant_dir)

