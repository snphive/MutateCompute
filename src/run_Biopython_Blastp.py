from src.Biopython import Biopy
import sys


# fasta_str = fasta_str=sys.argv[1]


def run_blastp(fasta_str=sys.argv[1]):
    return Biopy.run_blastp(fasta_str)
