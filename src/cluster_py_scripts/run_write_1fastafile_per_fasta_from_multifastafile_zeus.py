#!/usr/bin/env python3
"""
Script for writing 1 fastafile per mutant, read from a multifastafile, from within cluster, specified in job.q. Each file is
deleted immediately after use (to save disk space).
"""
import sys
import os
from src.tools.GeneralUtilityMethods import GUM
from src.enums.Paths import Paths
from src.enums.Str import Str

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

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
