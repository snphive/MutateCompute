from src.Paths import Paths

from enum import Enum


class Str(Enum):
    TRUE_lc = 'true'
    FALSE_lc = 'false'
    NEWLN = '\n'
    TAB = "\t"
    SEMICO_NL = ';' + NEWLN
    SPCE = ' '
    FSLSH = "/"
    FSLSH_ASTRX = FSLSH + '*'
    FSLSH_DBLASTRX_FSLSH_ASTRX = FSLSH_ASTRX + '*' + FSLSH_ASTRX
    FSTAEXT = '.fasta'
    PDBEXT = '.pdb'
    DOTS3 = '...'
    PDBS = Paths.DIR_PDBS.value
    FASTAS = Paths.DIR_FASTAS.value
    HASH = '#'
    JOBQ = 'job.q'
    runscrpt_txt = 'runscript.txt'
    COND_TEMP = 'temp'
    COND_ION_STRGTH = 'ion_strgth'
    COND_PH = 'ph'
    COND_TFE = 'tfe'
    COND_STAB = 'stab'
    COND_CONC = 'conc'
    PARTALLPATHEXISTS_MSG = 'Part/all dir tree already exists.'
    PROBLNXCMD_MSG = 'Problem with linux command: '
