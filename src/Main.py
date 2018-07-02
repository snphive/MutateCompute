from src.MutateFasta import MutateFasta
from src.Agadir import Agadir
from src.FoldX import FoldX
from src.Scheduler import Scheduler

import glob
import yaml
import pydevd
pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)


mutateFasta = None
agadir = None
foldx = None
scheduler = None

number_of_pdbs_to_process = 0
number_of_FASTA_files_to_process = 0
input_pdbs = []
input_fastas = []
do_mutate = False
do_agadir = False
do_foldx_repair = False
do_foldx_buildmodel = False
do_foldx_stability = False
do_foldx_analysecomplex = False


with open("/switchlab/group/shazib/SnpEffect/SourceFiles/Scripts/pathsAndDictionaries.yaml",
          'r') as stream:
    try:

        paths_and_dictionaries = yaml.load(stream)
        path_R_exe = paths_and_dictionaries['ROOT']['path_R_exe']
        path_FoldX_exe = paths_and_dictionaries['ROOT']['path_FoldX_exe']
        path_Agadir_exe = paths_and_dictionaries['ROOT']['path_Agadir_exe']
        path_QSub_exe = paths_and_dictionaries['ROOT']['path_QSub_exe']

        path_SnpEffect_dir = paths_and_dictionaries['ROOT']['path_SnpEffect_dir']
        path_SE_SourceFiles_Scripts_dir = paths_and_dictionaries['ROOT'][' path_SourceFiles_Scripts_dir']
        path_SE_Inputs_PDBs_dir = paths_and_dictionaries['ROOT']['path_SE_Inputs_PDBs_dir']
        path_SE_Inputs_FASTAs_dir = paths_and_dictionaries['ROOT']['path_SE_Inputs_FASTAs_dir']
        path_SE_Outputs_dir = paths_and_dictionaries['ROOT']['path_SE_Outputs_dir']
        path_SE_Outputs_Agadir_dir = paths_and_dictionaries['ROOT']['path_SE_Outputs_Agadir_dir']
        path_SE_Outputs_FoldX_dir = paths_and_dictionaries['ROOT']['path_SE_Outputs_FoldX_dir']

        dict_aa_1to3 = paths_and_dictionaries['ROOT']['dict_aa_1to3']
        dict_aa_3to1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']
        list_all_20_aa = paths_and_dictionaries['ROOT']['list_all_20_aa']

    except yaml.YAMLError as exc:
        print(exc)

mutate_compute_options_file = open('../options/MutateCompute_Options.txt', 'r').readlines()
for line in mutate_compute_options_file:
    if '#' in line:
        continue
    if 'PDBs:' in line:
        pdb_option = line.split(':')[-1].strip(';\n').strip()
        if pdb_option == 'All':
            pdb_paths = glob.glob(path_SE_Inputs_PDBs_dir + '*.pdb')
            for pdb_path in pdb_paths:
                input_pdbs.append(pdb_path.split('/')[-1])
        elif pdb_option == '':
            number_of_pdbs_to_analyse = 0
        elif not isinstance(pdb_option, str):
            number_of_pdbs_to_analyse = int(pdb_option)
        for x in range(1, number_of_pdbs_to_analyse + 1):
            input_pdbs.append(pdb_paths[x].split('/')[-1])
    print('PDBs to analyse:\t\t' + ",\t".join(input_pdbs))

    if 'FASTAs:' in line:
        fasta_option = line.split(':')[-1].strip(';\n').strip()
        if fasta_option == 'All':
            fasta_paths = glob.glob(path_SE_Inputs_FASTAs_dir + '*.fasta')
            for fasta_path in fasta_paths:
                input_fastas.append(fasta_path.split('/')[-1])
        elif fasta_option == '':
            number_of_fastas_to_analyse = 0
        elif not isinstance(fasta_option, str):
            number_of_fastas_to_analyse = int(fasta_option)
        for x in range(1, number_of_fastas_to_analyse + 1):
            input_fastas.append(fasta_paths[x].split('/')[-1])
    print('FASTA files to analyse:\t\t' + ",\t".join(input_fastas))

    if 'MUTATE_FASTA:' in line:
        do_mutate_fasta = line.split(':')[-1].strip(';\n').strip()
    if 'AGADIR:' in line:
        do_agadir = line.split(':')[-1].strip(';\n').strip()
    if 'FOLDX_REPAIR:' in line:
        do_foldx_repair = line.split(':')[-1].strip(';\n').strip()
    if 'FOLDX_BUILDMODEL:' in line:
        do_foldx_buildmodel = line.split(':')[-1].strip(';\n').strip()
    if 'FOLDX_STABILITY:' in line:
        do_foldx_stability = line.split(':')[-1].strip(';\n').strip()
    if 'FOLDX_ANALYSECOMPLEX:' in line:
        do_foldx_analysecomplex = line.split(':')[-1].strip(';\n').strip()

operations = {
    'do_mutate_fasta': do_mutate_fasta,
    'do_agadir': do_agadir,
    'do_foldx_repair': do_foldx_repair,
    'do_foldx_buildmodel': do_foldx_buildmodel,
    'do_foldx_stability': do_foldx_stability,
    'do_foldx_analysecomplex': do_foldx_analysecomplex,
}

if do_mutate_fasta or do_agadir or do_foldx_repair or do_foldx_buildmodel or do_foldx_stability or do_foldx_analysecomplex:
    scheduler = Scheduler(operations, input_pdbs, input_fastas, list_all_20_aa)
    scheduler.start()
else:
    print('All of the operations are set to FALSE in MutateCompute_Options.txt file - nothing to do')
