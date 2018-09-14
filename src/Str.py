from src.Paths import Paths

from enum import Enum


class Str(Enum):
    TRUE_lc = 'true'
    FALSE_lc = 'false'
    NEWLN = '\n'
    SEMICO_NL = ';\n'
    SPCE = ' '
    FSLSH = '/'
    FSTAEXT = '.fasta'
    PDBEXT = '.pdb'
    DOTS3 = '...'
    PDBS = Paths.DIR_PDBS.value
    FASTAS = Paths.DIR_FASTAS.value
    HASH = '#'
