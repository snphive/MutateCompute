#!/usr/bin/env python3
"""
Script for running Agadir and/or FoldX on fasta/pdb files.

Includes functionality for reducing disk space usage by deleting input files after computation completed and deleting output
files that are not needed.

Will include functionality for read data from output files and write to summary files and/or MySQL database.

KickOffManual.py can be run locally or on cluster.
"""
import sys
import os
import glob
from src.Main import Main
from src.enums.Paths import Paths
from src.enums.AminoAcids import AA
from src.enums.Str import Str
from src.enums.Conditions import Cond
from src.FoldX import FoldX
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

if 'switchlab/group' in os.getcwd() and sys.argv[1].strip(' ') != 'use_cluster=True':

    raise ValueError('Current working directory includes "/switchlab/group" in path. Hence you may be running the '
                     'program on the psb cluster, while the "use_cluster" parameter is False. '
                     'If launching from one of the bash scripts, check you are using bash/KickOffZeus.sh which '
                     'includes ("use_cluster=True") parameter, as opposed to the bash/KickOff.sh script. \n'
                     'Alternatively, uncomment the Paths.set_up_paths(use_cluster=True) line below to manually set it' 
                     'here.')
# Paths.set_up_paths(use_cluster=True)

if 'switchlab/group' not in os.getcwd() and sys.argv[1].strip(' ') == 'use_cluster=True':
    raise ValueError('Current working directory does not include "/switchlab/group" in path. Hence you might not be ' 
                     'running the program on the psb cluster, while the "use_cluster" parameter is True. '
                     'If launching from one of the bash scripts, check you are using bash/KickOff.sh which includes '
                     '("use_cluster=False") parameter, as opposed to the bash/KickOffZeus.sh script. \n'
                     'Alternatively, uncomment the Paths.set_up_paths(use_cluster=False) line below to manually set it' 
                     'here.')
# Paths.set_up_paths(use_cluster=False)

"""
Set the value(s) of whichever operation(s) you want to run to True in the operations dict. 
operations is passed to Main via its constructor.
"""
operations = {'do_mutate_fasta': False, 'do_agadir': False, 'do_foldx_repair': False, 'do_foldx_buildmodel': True,
              'do_foldx_stability': False, 'do_foldx_analysecomplex': False}

"""
Multithreading is not up and running yet (Nov18). Set to False for now.
"""
use_multithread = False

"""
Get the pdb files you want to run FoldX on.
"""
# path_pdbfiles = []
# path_input_pdbs_dir = '/switchlab/group/shazib/SnpEffect/output_data/analyse_complex'
# path_input_pdbs_dir = Paths.OUTPUT_AC
# path_pdbfiles = sorted(glob.glob(path_input_pdbs_dir + '/**/*.pdb', recursive=True))
path_pdbfiles = [os.path.join(Paths.INPUT_PDBS, 'RepairPDB_14' + Str.PDBEXT.value)]
if not path_pdbfiles:
    raise ValueError('No pdb files to process. Check paths are correct and check files are where you expect.')
"""
Select specific mutants if you are only interested in these.

BE SURE TO SET THIS TO EMPTY LIST IF YOU DON'T WANT ANY OF THE SUBSEQUENT ACTIONS BELOW TO BE SPECIFIC TO THIS/THESE MUTANTS ONLY.
"""
specific_fxmutants = ['AA101A']
# specific_fxmutants = []
"""
Get the fasta files you want to run mutate_fasta or agadir on.
"""
path_fastafiles = []
# path_input_fastas_dir = Paths.INPUT_MUTS_MULTIFASTAS_29611_1000 + '/1...250/'
# path_fastafiles = sorted(glob.glob(path_input_fastas_dir + '/**/*.fasta', recursive=True))
# path_fastafiles = sorted(glob.glob(path_input_fastas_dir + '/*.fasta'))
# if not path_fastafiles:
#     raise ValueError('No fasta files to process. Check paths are correct and check files are where you expect.')

"""
Kick off the program(s) via the constructor or Main class.
"""
main = Main(operations, use_multithread, Paths.INPUT, Paths.OUTPUT, path_pdbfiles, path_fastafiles, specific_fxmutants,
            AA.LIST_ALL_20_AA.value)

"""
After computation completed, DELETE config files no longer needed.  
"""
# if operations['do_foldx_buildmodel']:
#     fx = FoldX()
#     path_output_bm_pdb_fxmutant_dirs = []
#     for path_pdbfile in path_pdbfiles:
#         pdbname = os.path.basename(path_pdbfile).split('.')[0]
#         for specific_fxmutant in specific_fxmutants:
#             path_output_bm_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname, specific_fxmutant))
#         if not specific_fxmutants:
#             path_output_bm_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_BM, pdbname, '*'))
#         for path_output_bm_pdb_fxmutant_dir in path_output_bm_pdb_fxmutant_dirs:
#             if fx.BuildModel(Cond.INCELL_MAML_FX.value).has_already_generated_dif_bm_fxoutfile(path_output_bm_pdb_fxmutant_dir):
#                 fx.remove_config_files(path_output_bm_pdb_fxmutant_dir)

if operations['do_foldx_analysecomplex']:
    fx = FoldX()
    path_output_bm_pdb_fxmutant_dirs = []
    path_output_ac_pdb_fxmutant_dirs = []
    for path_pdbfile in path_pdbfiles:
        pdbname = os.path.basename(path_pdbfile).split('.')[0]
        for specific_fxmutant in specific_fxmutants:
            path_output_bm_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname, specific_fxmutant))
            path_output_ac_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_AC, pdbname, specific_fxmutant))
        if not specific_fxmutants:
            path_output_ac_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_AC, pdbname, '*'))
        for path_output_ac_pdb_fxmutant_dir in path_output_ac_pdb_fxmutant_dirs:
            if fx.AnalyseComplex(Cond.INCELL_MAML_FX.value).has_already_generated_summary_ac_fxoutfile(
                    path_output_ac_pdb_fxmutant_dir):
                fx.remove_config_files(path_output_ac_pdb_fxmutant_dir)
                fx.remove_pdbfiles(os.path.join(Paths.OUTPUT_BM, pdbname, os.path.basename(path_output_ac_pdb_fxmutant_dir)))
                fx.remove_pdbfiles(os.path.join(path_output_ac_pdb_fxmutant_dir))

"""
After computation completed, write result values to summary and/or MySQL database.  
"""



# """
# After deletion of unnecessary files, COMPRESS remaining files, (that will be copied over to storage directories)
# """
# if self._has_already_generated_summary_ac_fxoutfile(path_output_ac_pdb_fxmutant_dir):
#     fx.remove_all_pdbfiles(path_output_bm_pdb_fxmutant_dir)
#




# pydevd.stoptrace()

