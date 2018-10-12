import sys
import glob
import natsort
from src.Main import Main
from src.Paths import Paths
from src.AminoAcids import AA
import pydevd
pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

"""
This module is for doing manual runs, locally or on cluster, either to quickly check something works or just to do a 
single input rather than a batch which is more for KickOff.py.
I am manually setting the use_cluster and path_inputs.. variables here but this is not essential.  
"""
use_cluster = True
Paths.set_up_paths(use_cluster)
operations = {'do_mutate_fasta': False, 'do_agadir': False, 'do_foldx_repair': False, 'do_foldx_buildmodel': True,
              'do_foldx_stability': False, 'do_foldx_analysecomplex': False}
use_multithread = False
# path_input_fastas_dir = '/Users/u0120577/REPO_PDB_FASTA/fastas/29611_fastas_1000/1...1000'
# path_fastafiles = sorted(glob.glob(path_input_fastas_dir + '/**/*.fasta', recursive=True))
path_fastafiles = []
# if not path_fastafiles:
#     raise ValueError('No fasta files to process. Check paths are correct and check files are where you expect.')
# path_pdbfiles = []
# path_input_pdbs_dir = '/Users/u0120577/REPO_PDB_FASTA/pdbs/32431_repairedpdbs_1000/1...1000'
path_input_pdbs_dir = '/switchlab/group/shazib/SnpEffect/input_data/RepairPDB_1'
path_pdbfiles = sorted(glob.glob(path_input_pdbs_dir + '/**/*.pdb', recursive=True))
if not path_pdbfiles:
    raise ValueError('No pdb files to process. Check paths are correct and check files are where you expect.')
main = Main(operations, use_multithread, Paths.INPUT, Paths.OUTPUT, path_pdbfiles, path_fastafiles,
            AA.LIST_ALL_20_AA.value)

pydevd.stoptrace()

