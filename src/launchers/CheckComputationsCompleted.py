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
Note: ac=True or bm=True strings can be passed as 2nd or 3rd arguments from the bash script. 
"""
Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))
check_buildmodel_completed = sys.argv[2].strip(' ') == 'bm=True'
if not check_buildmodel_completed:
    check_buildmodel_completed = sys.argv[3].strip(' ') == 'bm=True'
check_analysecomplex_completed = sys.argv[2].strip(' ') == 'ac=True'
if not check_analysecomplex_completed:
    check_analysecomplex_completed = sys.argv[3].strip(' ') == 'ac=True'

"""
Select parameters required.
"""
pdbfile = 'RepairPDB_14' + Str.PDBEXT.value
path_pdbfile = os.path.join(Paths.INPUT_PDBS, pdbfile)
amino_acids = AA.LIST_ALL_20_AA.value

"""
Check for dif_BuildModel_RepairPDB_x.fxout file.
"""
if check_buildmodel_completed:
    num_of_missing_dif_fxoutfiles = FoldX().BuildModel(Cond.INCELL_MAML_FX.value).\
        find_num_of_missing_avg_bm_fxoutfiles(path_pdbfile, amino_acids)
    if num_of_missing_dif_fxoutfiles != 0:
        print('Warning: BuildModel has not completed computations for this pdb: ' + os.path.basename(path_pdbfile))
        print('Number of mutants missing dif files = ' + str(num_of_missing_dif_fxoutfiles))
    else:
        print('BuildModel has completed all computations for this pdb: ' + os.path.basename(path_pdbfile))

"""
Check for AnalyseComplex_RepairPDB_x_1_x.fxout file.
"""
if check_analysecomplex_completed:
    num_of_missing_interaction_energies = FoldX().AnalyseComplex(Cond.INCELL_MAML_FX.value).\
        find_num_of_missing_interaction_energies(path_pdbfile, amino_acids)
    if num_of_missing_interaction_energies != 0:
        print('Warning: AnalyseComplex has not completed computations for this pdb: ' + os.path.basename(path_pdbfile))
        print('Number of mutants missing interaction energies = ' + str(num_of_missing_interaction_energies))
    else:
        print('AnalyseComplex has completed interaction energy computations for this pdb: ' + os.path.basename(path_pdbfile))


# pydevd.stoptrace()
