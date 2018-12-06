import mysql.connector
from src.database.DBConnectDDL import DBCnxDDL


"""
An approximation of Data Access Object pattern.

"""


class DAO(object):

    def __init(self):
        try:
            DBCnxDDL.connect_mysql_practice_db()
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))

    from enum import Enum

    class SQLDML(Enum):

        sql_INSRT = 'INSERT'
        sql_SLCT = 'SELECT'
        sql_FRM = 'FROM'
        sql_WHR = 'WHERE'
        sql_GRPBY = 'GROUP_BY'

    class SQLDDL(Enum):
        sql_CRT_TBL = 'CREATE TABLE'
        sql_ALTR = 'ALTER TABLE'
        sql_DRP = 'DROP'


