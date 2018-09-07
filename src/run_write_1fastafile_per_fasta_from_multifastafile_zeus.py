import sys
import os
from src.GeneralUtilityMethods import GUM
from src.Paths import Paths

path_fastafile = sys.argv[1]

path_output_filename = GUM.make_root_input_3dots_dirs(Paths.LOCAL_IO_ROOT.value, path_fastafile)
with open(path_fastafile) as f:
    fasta_str = ''
    is_first_line = True
    for line in f.readlines():
        if '>' in line:
            if not is_first_line:
                fastafile = line.split('>')[-1].split('\n')[0] + '.fasta'
                with open(os.path.join(path_output_filename, fastafile), 'w') as f_to_write:
                    f_to_write.write(fasta_str)
            fasta_str = line
            is_first_line = False
        else:
            fasta_str += line
