import sys
import glob
import natsort
from src.Main import Main
from src.Paths import Paths
import pydevd
pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

if len(sys.argv) < 2:
    use_cluster = False
else:
    use_cluster = True if sys.argv[1] == 'use_cluster=True' else False

Paths.set_up_paths(use_cluster)

operations = {'do_mutate_fasta': False, 'do_agadir': True,
              'do_foldx_repair': False, 'do_foldx_buildmodel': False,
              'do_foldx_stability': False, 'do_foldx_analysecomplex': False}
use_multithread = False


# path_input_fastas_dir = Paths.MC_REPO_FASTAS.value + '/fastas_1000_rest'
if use_cluster:
    path_input_fastas_dir = '/switchlab/group/shazib/SnpEffect/input_data/mutants_multifastas/1...2'
else:
    path_input_fastas_dir = '/Users/u0120577/PycharmProjects/input_data/mutants_multifastas/1...2'
path_fastafiles = natsort.natsorted(glob.glob(path_input_fastas_dir + '/**/*.fasta', recursive=True))
if not path_fastafiles:
    raise ValueError('No fasta files to process. Check paths are correct and check files are where you expect.')

path_pdbfiles = []

main = Main(operations, use_multithread, Paths.INPUT, Paths.OUTPUT, path_pdbfiles, path_fastafiles, ['A', 'C'])
            # AA.LIST_ALL_20_AA.value)

pydevd.stoptrace()
