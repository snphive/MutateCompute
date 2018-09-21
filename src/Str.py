from src.Paths import Paths

from enum import Enum


class Str(Enum):
    TRUE_lc = 'true'
    FALSE_lc = 'false'
    NEWLN = '\n'
    TAB = "\t"
    SEMICO_NL = ';' + NEWLN
    SPCE = ' '
    FSLSH = '/'
    FSLSH_ASTRX = FSLSH + '*'
    FSTAEXT = '.fasta'
    PDBEXT = '.pdb'
    DOTS3 = '...'
    PDBS = Paths.DIR_PDBS.value
    FASTAS = Paths.DIR_FASTAS.value
    HASH = '#'
    JOBQ = 'job.q'
    runscrpt_txt = 'runscript.txt'
    TEMP = 'temp'
    ION_STRGTH = 'ion_strgth'
    PH = 'ph'
    TFE = 'tfe'
    STAB = 'stab'
    CONC = 'conc'
