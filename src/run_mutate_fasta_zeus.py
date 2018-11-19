import sys
from src.AminoAcids import AA
from src.MutateFasta import MutateFasta

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

path_fastafile = sys.argv[1]
write_1_fasta_only = (sys.argv[2] == 'True')
write_fasta_per_mut = (sys.argv[3] == 'True')
path_output_3dots = sys.argv[4]

print('path_fastafile: ' + path_fastafile)
print('write_1_fastafile: ' + str(write_1_fasta_only))
print('write_fasta_per_mut: ' + str(write_fasta_per_mut))
print('path_output_3dots: ' + path_output_3dots)

"""
Mutates FASTA (typically wild-type) sequences at every position, to every residue specified, (typically the 19
other residues.)
:sys.argv[1] path_fastafile: Absolute path to fasta file.
:sys.argv[2] write_1_fasta_only: True writes 1 fasta file containing multiple fasta sequences (i.e. the mutants).
:sys.argv[3] write_fasta_per_mut: True writes a fasta file for each fasta sequences (i.e. 1 per mutant).
:sys.argv[4] path_output_3dots: Absolute path to a '3dot' subdirectory (i.e. 1...1000).
"""

print('run_mutate_fasta_zeus.py ###################################')
mutatefasta = MutateFasta(amino_acids=AA.LIST_ALL_20_AA.value)
mutatefasta.mutate_every_residue(path_fastafile, write_1_fasta_only, write_fasta_per_mut, path_output_3dots)
