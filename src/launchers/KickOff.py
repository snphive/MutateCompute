#!/usr/bin/env python3
"""
Script for running Agadir and/or FoldX on fasta/pdb files.

Includes functionality for reducing disk space usage by deleting input files after computation completed and deleting output
files that are not needed.

Will include functionality for read data from output files and write to summary files and/or MySQL database.

KickOff.py can be run locally or on cluster.
"""
import sys
import os
import glob
from src.tools.GeneralUtilityMethods import GUM
from src.Main import Main
from src.enums.Paths import Paths
from src.enums.AminoAcids import AA
from src.enums.Str import Str
from src.enums.Conditions import Cond
from src.tools.OutputsParser import Parser
from src.FoldX import FoldX
import pydevd
pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

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
operations = {Str.OPER_RUN_MUT_FSTA.value: False, Str.OPER_RUN_AGDR.value: False, Str.OPER_RUN_FX_RPR.value: False, Str.OPER_RUN_FX_BM.value: False,
              Str.OPER_RUN_FX_STAB.value: False, Str.OPER_RUN_FX_AC.value: False}

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
pdbname = 'RepairPDB_1'
path_pdbfiles = [os.path.join(Paths.INPUT_PDBS, pdbname + Str.PDBEXT.value)]
if not path_pdbfiles:
    raise ValueError('No pdb files to process. Check paths are correct and check files are where you expect.')
"""
Select specific mutants if you are only interested in these.

BE SURE TO SET THIS TO EMPTY LIST IF YOU DON'T WANT ANY OF THE SUBSEQUENT ACTIONS BELOW TO BE SPECIFIC TO THIS/THESE MUTANTS ONLY.
"""
# specific_fxmutants = ['AA101A']
specific_fxmutants = []
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
if operations[Str.OPER_RUN_FX_BM.value]:
    fx = FoldX()
    path_output_bm_pdb_fxmutant_dirs = []
    for path_pdbfile in path_pdbfiles:
        pdbname = os.path.basename(path_pdbfile).split('.')[0]
        for specific_fxmutant in specific_fxmutants:
            path_output_bm_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname, specific_fxmutant))
        if not specific_fxmutants:
            path_output_bm_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_BM, pdbname, '*'))
        for path_output_bm_pdb_fxmutant_dir in path_output_bm_pdb_fxmutant_dirs:
            bm = fx.BuildModel(Cond.INCELL_MAML_FX.value)
            if bm.has_already_generated_avg_bm_fxoutfile(path_output_bm_pdb_fxmutant_dir):
                fx.remove_config_files(path_output_bm_pdb_fxmutant_dir)

if operations[Str.OPER_RUN_FX_AC.value]:
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
            ac = fx.AnalyseComplex(Cond.INCELL_MAML_FX.value)
            if ac.has_already_generated_summary_ac_fxoutfile(
                    path_output_ac_pdb_fxmutant_dir):
                fx.remove_config_files(path_output_ac_pdb_fxmutant_dir)
                # fx.remove_pdbfiles(os.path.join(Paths.OUTPUT_BM, pdbname, os.path.basename(path_output_ac_pdb_fxmutant_dir)))
                # fx.remove_pdbfiles(os.path.join(path_output_ac_pdb_fxmutant_dir))

"""
Choose post-computation actions:
"""
write_bm_to_csv = False
write_bm_to_db = False
write_ac_to_csv = False
write_ac_to_db = False
pack_compress_bm_outputs = False
pack_compress_ac_outputs = True

"""
After computations completed, WRITE RESULTS to CSV files.  
"""
if write_bm_to_csv:
    fx = FoldX()
    path_output_bm_pdb_fxmutant_dirs = []
    for path_pdbfile in path_pdbfiles:
        pdbname = os.path.basename(path_pdbfile).split('.')[0]
        for specific_fxmutant in specific_fxmutants:
            path_output_bm_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname, specific_fxmutant))
        if not specific_fxmutants:
            path_output_bm_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_BM, pdbname, '*'))
        for path_output_bm_pdb_fxmutant_dir in path_output_bm_pdb_fxmutant_dirs:
            bm = fx.BuildModel(Cond.INCELL_MAML_FX.value)
            bm.write_bm_avg_fxout_to_csvfile_up_1dirlevel(path_output_bm_pdb_fxmutant_dir)
            GUM.linux_remove_dir(path_output_bm_pdb_fxmutant_dir)

# NEED TO ESTABLISH WHETHER YOU SHOULD READ FROM BOTH THE WT AND MUTANT SUMMARY FILE AND TAKE THE DIFFERENCE.
# NEED TO KNOW WHAT VALUES ARE SIGNIFICANT, E.G. ANYTHING < 0.01 KCAL/MOL SIGINIFICANT ??
if write_ac_to_csv:
    fx = FoldX()
    path_output_ac_pdb_fxmutant_dirs = []
    for path_pdbfile in path_pdbfiles:
        pdbname = os.path.basename(path_pdbfile).split('.')[0]
        for specific_fxmutant in specific_fxmutants:
            path_output_ac_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_AC, pdbname, specific_fxmutant))
        if not specific_fxmutants:
            path_output_ac_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_AC, pdbname, '*'))
        for path_output_ac_pdb_fxmutant_dir in path_output_ac_pdb_fxmutant_dirs:
            ac = fx.AnalyseComplex(Cond.INCELL_MAML_FX.value)
            ac.write_ac_sumry_fxout_to_csvfile_up_1dirlevel(path_output_ac_pdb_fxmutant_dir)
            GUM.linux_remove_dir(path_output_ac_pdb_fxmutant_dir)

"""
After computations completed, WRITE RESULTS to DATABASE.  
"""


"""
After computations completed, and results written, COMPRESS RESULTS to one tar per pdb (per algorithm) for transfer and storage.  
"""

if pack_compress_bm_outputs:
    for path_pdbfile in path_pdbfiles:
        pdbname = os.path.basename(path_pdbfile).split('.')[0]
        path_files_to_pack_dir = os.path.join(Paths.OUTPUT_BM, pdbname)
        Parser().make_tarfile(path_files_to_pack_dir)


if pack_compress_ac_outputs:
    for path_pdbfile in path_pdbfiles:
        pdbname = os.path.basename(path_pdbfile).split('.')[0]
        path_files_to_pack_dir = os.path.join(Paths.OUTPUT_AC, pdbname)
        Parser().make_tarfile(path_files_to_pack_dir)


pydevd.stoptrace()
