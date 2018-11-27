#!/usr/bin/env python3
"""
Script for launching the searching output files for correct number of output folders and output files.
E.g. Check that the complete repetoire of FoldX mutants were generated and that the desired output files were produced.

CheckComputationsCompleted.py can be run locally or on cluster (the latter via a bash script).
"""
import sys
import os
from src.FoldX import FoldX
from src.enums.Conditions import Cond
from src.enums.Paths import Paths
from src.enums.AminoAcids import AA
from src.enums.Str import Str
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
Paths for the entire codebase are set accordingly. "use_cluster" is set to False by default.  
"""
Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))

"""
Select parameters required.
"""
pdbfile = 'RepairPDB_14' + Str.PDBEXT.value
path_pdbfile = os.path.join(Paths.INPUT_PDBS, pdbfile)
amino_acids = AA.LIST_ALL_20_AA.value

"""
Check for all expected outputs.
"""

if not FoldX().BuildModel(Cond.INCELL_MAML_FX.value).confirm_all_dif_bm_fxoutfiles_computed(path_pdbfile, amino_acids):
    print('Warning: BuildModel has not completed all computations for this pdb: ' + os.path.basename(
        path_pdbfile))
else:
    print('BuildModel has completed all computations for this pdb: ' + os.path.basename(path_pdbfile))


# pydevd.stoptrace()
