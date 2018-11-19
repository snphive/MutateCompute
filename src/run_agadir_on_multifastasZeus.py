import sys
import os
from src.GeneralUtilityMethods import GUM
from src.Agadir import Agadir
from src.Agadir import AgadCndtns
from src.Str import Str

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

path_fastafile = sys.argv[1]
path_dst = sys.argv[2]

"""
To run agadir on file that has multiple fasta sequences. Iterate through the list, write the sequence to an
individual text file passing this to agadir.compute(). This newly-written individual text file is then deleted,
crucial to prevent too much memory taken up by each file when running large numbers of sequences.
:sys.argv[1] path_fastafile: Abs path to fasta file.
:sys.argv[2] path_dst: Abs path to output root dir.
"""

print('run_agadir_on_multifastasZeus.py ###################################')
# agadir = Agadir(AgadCndtns.INCELL_MAML.value)
path_dst = GUM.make_path_agadir_3dots_filename_mutants_dirs(path_dst, path_fastafile, add_filename_subdir=True)
with open(path_fastafile) as f:
    is_first_line = True
    fasta_str = ''
    mutantfastafilename = ''
    mutantfastafile = ''
    for line in f.readlines():
        if '>' in line:
            if not is_first_line:
                path_dst_mutant_filename = GUM._os_makedirs(path_dst, mutantfastafilename)
                path_dst_mutant_file = os.path.join(path_dst_mutant_filename, mutantfastafile)
                with open(path_dst_mutant_file, 'w') as g:
                    g.write(fasta_str)
                agadir = Agadir(AgadCndtns.INCELL_MAML.value)
                agadir.compute(path_dst_mutant_file)
                GUM.linux_remove_file(path_dst_mutant_file)
            fasta_str = line
            is_first_line = False
            mutantfastafilename = line.split('>')[-1].split('\n')[0]
            mutantfastafile = mutantfastafilename + Str.FSTAEXT.value
        else:
            fasta_str += line
