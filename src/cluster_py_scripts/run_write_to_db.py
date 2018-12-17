"""
Script for writing data to database.
"""
import os
import glob
from src.enums.Paths import Paths
from src.database.DAL import DAL
from src.enums.Str import Str


__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


write_ac_to_db = False
write_bm_to_db = False

path_output_ac_ddG_csvfiles = glob.glob(os.path.join(Paths.OUTPUT_AC, '/**/' + '*' + Str.CSVEXT.value), recursive=True)
path_output_bm_ddG_csvfiles = glob.glob(os.path.join(Paths.OUTPUT_BM, '/**/' + '*' + Str.CSVEXT.value), recursive=True)

"""
Write results to database.  
"""
if write_ac_to_db:
    for path_output_ac_ddG_csvfile in path_output_ac_ddG_csvfiles:
        DAL().write_csv_to_db(path_output_ac_ddG_csvfile)

if write_bm_to_db:
    for path_output_bm_ddG_csvfile in path_output_bm_ddG_csvfiles:
        DAL().write_csv_to_db(path_output_bm_ddG_csvfile)
