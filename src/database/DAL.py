
"""
An approximation of Data Access Layer pattern.
Class for connecting to relational databases and performing SQL DML commands: e.g. INSERT, UPDATE, DELETE.


NOTE: the use of f-strings for string interpolation was only introduced in Python 3.6. The current (7Dec2018) latest version of
Python on the cluster is 3.5. Therefore, while this is my prefered syntax for insert_value_by_column(), I have provided an
alternative to prevent the syntax error it causes which others prevents the codebase from starting.
"""
import mysql.connector
from src.database.DBConnectDDL import DBCnxDDL
__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class DAL(object):

    def write_csv_to_db(self, path_output_bm_or_ac_pdb_sumry_csvfile: str):
        with open(path_output_bm_or_ac_pdb_sumry_csvfile) as f:
            csvfile_lines = f.readlines()
        column_names = csvfile_lines[0].split(',')
        column_values = csvfile_lines[1].split(',')
        db = DBCnxDDL.connect_mysql_practice_db()
        self.insert_values_by_column(db, column_names, column_values)
        DBCnxDDL().disconnect_mysql_db(db)

    def insert_values_by_column(self, mysql_practice_cnx, column_names: list, column_values: list):
        cursor = mysql_practice_cnx.cursor()
        column_names = ','.join(column_names)
        column_values = ','.join(['\"' + str(value) + '\"' for value in column_values])
        sql_cmd = 'INSERT INTO fx_data{columns} VALUES{values}'.format(column_names=column_names, column_values=column_values)
        # sql_cmd = f"""INSERT INTO fx_data({columns}) VALUES({values})"""
        try:
            cursor.execute(sql_cmd)
            mysql_practice_cnx.commit()
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))
            mysql_practice_cnx.rollback()
