
from src.Mutate import Mutate
from src.Agadir import Agadir
from src.FoldX import FoldX
import glob
import yaml
import pydevd
pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)


mutate_instance = ''
agadir_instance = ''
foldx_instance = ''
number_of_pdbs_to_process = 0
number_of_FASTA_files_to_process = 0
pdb_list = []

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
        path_SE_PDBs_Inputs_dir = paths_and_dictionaries['ROOT']['path_SE_PDBs_Inputs_dir']
        path_SE_Outputs_dir = paths_and_dictionaries['ROOT']['path_SE_Outputs_dir']
        path_SE_Outputs_Agadir_dir = paths_and_dictionaries['ROOT']['path_SE_Outputs_Agadir_dir']
        path_SE_Outputs_FoldX_dir = paths_and_dictionaries['ROOT']['path_SE_Outputs_FoldX_dir']

        dict_aa_1to3 = paths_and_dictionaries['ROOT']['dict_aa_1to3']
        dict_aa_3to1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']
        list_all_20_aa = paths_and_dictionaries['ROOT']['list_all_20_aa']

        mutate_instance = Mutate()
        agadir_instance = Agadir()
        foldx_instance = FoldX()

    except yaml.YAMLError as exc:
        print(exc)

    mutate_compute_options_file = open('../options/MutateCompute_Options.txt', 'r').readlines()
    for line in mutate_compute_options_file:
        if '#' in line:
            continue
        if 'PDBs:' in line:
            if line.split(':')[-1].strip(';\n').strip() == "All":
                pdb_paths = glob.glob(path_SE_PDBs_Inputs_dir + '*.pdb')
                for pdb_path in pdb_paths:
                    pdb_temp = pdb_path.split('/')[-1]
                    pdb_list.append(pdb_temp)
            pdb_string = line.split(':')[-1].strip(';\n')
            if ',' in pdb_string:
                pdbs_temp = pdb_string.split(',')
                for pdb_temp in pdbs_temp:
                    pdb_list.append(pdb_temp)
            else:
                    pdb_list.append(pdb_string)
    print ('PDBs to analyse:\t\t' + ",\t".join(pdb_list))



