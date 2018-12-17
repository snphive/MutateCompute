#!/usr/bin/env python3
"""
Script for connecting to and reading/writing from/to a MySQL database.
"""
import os
import glob
from src.enums.Paths import Paths
from src.enums.Str import Str
from src.database.DBConnectDDL import DBCnxDDL
from src.database.DAL import DAL
__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


"""
1. Make connection
"""
db = DBCnxDDL()
cnx = db.connect_mysql_practice_db()

"""
2. To test out connection by inserting some columns and values..    
"""
# db.add_new_columns('ddG', 'FLOAT')
DAL().insert_values_by_column(["pdb", "mutation", "ddG"], ["RepairPDB_2", "AA101C", 1.23456789])

"""
3. Select which results to write to database.  
"""
write_ac_to_db = False
write_bm_to_db = False

"""
4. Write to database.  
"""
if write_ac_to_db:
    path_output_ac_ddG_csvfiles = glob.glob(os.path.join(Paths.OUTPUT_AC, '/**/' + '*' + Str.CSVEXT.value), recursive=True)
    for path_output_ac_ddG_csvfile in path_output_ac_ddG_csvfiles:
        DAL().write_csv_to_db(path_output_ac_ddG_csvfile)

if write_bm_to_db:
    path_output_bm_ddG_csvfiles = glob.glob(os.path.join(Paths.OUTPUT_BM, '/**/' + '*' + Str.CSVEXT.value), recursive=True)
    for path_output_bm_ddG_csvfile in path_output_bm_ddG_csvfiles:
        DAL().write_csv_to_db(path_output_bm_ddG_csvfile)

"""
5. Terminate connection
"""
db.disconnect_mysql_db(cnx)

