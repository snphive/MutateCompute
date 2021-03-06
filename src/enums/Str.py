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
    TXTEXT = '.txt'
    FSTAEXT = '.fasta'
    XMLEXT = '.xml'
    JSONEXT = '.json'
    CSVEXT = '.csv'
    PDBEXT = '.pdb'
    Q_EXT = '.q'
    LOGEXT = '.log'
    TARGZEXT = '.tar.gz'
    FXOUTEXT = '.fxout'
    CLSTR_ERR_LOGEXT = '.e'
    CLSTR_OUT_LOGEXT = '.o'
    DOTS3 = '...'
    PDBS = Paths.DIR_PDBS.value
    FASTAS = Paths.DIR_FASTAS.value
    HASH = '#'
    RNSCRPT = 'runscript'
    COND_TEMP = 'temp'
    COND_ION_STRGTH = 'ion_strgth'
    COND_PH = 'ph'
    COND_TFE = 'tfe'
    COND_STAB = 'stab'
    COND_CONC = 'conc'
    PARTALLPATHEXISTS_MSG = 'Part/all dir tree already exists.'
    PROBLNXCMD_MSG = 'Problem with linux command: '
    REPAIRPDB_ = 'RepairPDB_'
    CHMOD777 = 'chmod 777'
    DEFAULT_NUM_FX_RUNS_1 = 1
    OPER_RUN_MUT_FSTA = 'run_mutate_fasta'
    OPER_RUN_AGDR = 'run_agadir'
    OPER_RUN_FX_RPR = 'run_fx_repair'
    OPER_RUN_FX_BM = 'run_fx_buildmodel'
    OPER_RUN_FX_STAB = 'run_fx_stability'
    OPER_RUN_FX_AC = 'run_fx_analysecomplex'
