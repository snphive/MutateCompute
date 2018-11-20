#!/usr/bin/env python3
"""
Enum for strings used by more than 1 other class.
"""

from src.enums.Paths import Paths
from enum import Enum

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


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
    XMLEXT = '.xml'
    JSONEXT = '.json'
    CSVEXT = '.csv'
    PDBEXT = '.pdb'
    LOGEXT = '.log'
    TARGZEXT = '.tar.gz'
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
    REPAIRPDB_  = 'RepairPDB_'
    GLOBRESIDUE_OUT = 'PSX_globalresidue.out'
