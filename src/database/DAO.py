import mysql.connector
from src.database.DBConnection import DBCnx


"""
Class for Data Access Layer. Includes connection logic.

"""


class DAO(object):

    def __init(self):
        try:
            DBCnx.connect_mysql_practice_db()
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))

    from enum import Enum

    class SQL_DML(Enum):

        sql_SLCT = 'SELECT'
        sql_FRM = 'FROM'
        sql_WHR = 'WHERE'
        sql_GRPBY = 'GROUP_BY'

    class SQL_DDL(Enum):
        sql_CRT_TBL = 'CREATE TABLE'

