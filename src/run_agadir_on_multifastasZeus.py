import sys
import os
from src.GeneralUtilityMethods import GUM
from src.Agadir import Agadir
from src.Agadir import AgadCndtns
from src.Str import Str

path_fastafile = sys.argv[1]
path_output_root = sys.argv[2]

path_output_agadir_3dots = GUM.make_root_agadir_3dots_dirs(path_output_root, path_fastafile)
dir_3dots = GUM.get_3dots_dir(path_output_agadir_3dots)
agadir = Agadir(AgadCndtns.INCELL_MAML.value)
path_dst = GUM.make_root_agadir_3dots_filename_mutants_dirs(path_output_root, path_fastafile)
with open(path_fastafile) as f:
    fasta_str = ''
    is_first_line = True
    mutantfastafile = ''
    for line in f.readlines():
        if '>' in line:
            if not is_first_line:
                path_dst_mutant_filename = GUM._os_makedirs(path_dst, mutantfastafile.split('.')[0])
                path_dst_mutant_file = os.path.join(path_dst_mutant_filename, mutantfastafile)
                with open(path_dst_mutant_file, 'w') as temp_fastafile:
                    temp_fastafile.write(fasta_str)
                agadir = Agadir(AgadCndtns.INCELL_MAML.value)
                agadir.compute(path_dst_mutant_file)
                GUM.linux_remove_file(path_dst_mutant_file)
            fasta_str = line
            is_first_line = False
            mutantfastafile = line.split('>')[-1].split('\n')[0] + Str.FSTAEXT.value
        else:
            fasta_str += line
