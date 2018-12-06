#!/usr/bin/env python3
"""
Script for connecting to and reading/writing from/to a MySQL database.
"""
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
db.connect_mysql_practice_db()

# db.add_new_columns('ddG', 'FLOAT')
DAL().insert_values_by_column(["pdb", "mutation", "ddG"], ["RepairPDB_2", "AA101C", 1.23456789])

"""
3. Terminate connection
"""
db.disconnect_mysql_db()

