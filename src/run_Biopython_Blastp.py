from src.Biopython import Biopy
import sys

Biopy.run_blastp(fasta_str=sys.argv[1])
