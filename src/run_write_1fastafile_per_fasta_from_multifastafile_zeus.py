import sys
import os
from src.GeneralUtilityMethods import GUM
from src.Paths import Paths
from src.Str import Str

path_fastafile = sys.argv[1]
path_dst = GUM.make_path_agadir_3dots_filename_mutants_dirs(Paths.OUTPUT, path_fastafile)

with open(path_fastafile) as f:
    fasta_str = ''
    is_first_line = True
    mutantfastafile = ''
    for line in f.readlines():
        if '>' in line:
            if not is_first_line:
                path_dst_mutant_file = os.path.join(path_dst, mutantfastafile)
                with open(path_dst_mutant_file, 'w') as temp_fastafile:
                    temp_fastafile.write(fasta_str)
            fasta_str = line
            is_first_line = False
            mutantfastafile = line.split('>')[-1].split('\n')[0] + Str.FSTAEXT.value
        else:
            fasta_str += line
