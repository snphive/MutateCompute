#!/usr/bin/env python3
"""
Script for connecting to and reading/writing from/to a MySQL database.
"""
import mysql.connector
from src.enums.DBServer import Server
from src.database.DBConnection import DBCnx
# import mysqlx.connection

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

#
# connection = mysql.connector.connect(user='snpeffect_v5', passwd='R34WKKGR', host='127.0.0.1')
# #                                      auth_plugin='mysql_native_password')
# # , database='snpeffect_v5', port='3306')
# # connection = mysqlx.connection(user='root', password='K0yGrJ8(', host='127.0.0.1', database='testForScope',
# #                                      port='3306', auth_plugin='mysql_native_password')
# cursor = connection.cursor()
# cursor.execute("SHOW DATABASES")
# for x in cursor:
#     print(x)

# cnx = mysql.connector.connect(user='root', password='K0yGrJ8(', host='127.0.0.1', database='testForScope', port='3306',
#                               auth_plugin='mysql_native_password')
# cur = cnx.cursor()
# cur.execute("SHOW DATABASES")
# for x in cur:
#     print(x)
#
# cnx.close()

db = DBCnx()
db.connect_mysql_practice_db()
db.add_new_columns('ddG', 'FLOAT')
db.insert_values_by_column(["pdb", "mutation", "ddG"], ["Repair_PDB2", "AA101C", 1.23456789])
db.disconnect_mysql_practice_db()

