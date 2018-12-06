
"""
An approximation of Data Access Layer pattern.
Class for connecting to relational databases and performing SQL DML commands: e.g. INSERT, UPDATE, DELETE.

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

    def insert_values_by_column(self, mysql_practice_cnx, columns: list, values: list):
        cursor = mysql_practice_cnx.cursor()
        columns = ','.join(columns)
        # values = ','.join([str(value) for value in values])
        values = ','.join(['\"' + str(value) + '\"' for value in values])
        sql_cmd = f"""INSERT INTO fx_data({columns}) VALUES({values})"""
        try:
            cursor.execute(sql_cmd)
            mysql_practice_cnx.commit()
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))
            mysql_practice_cnx.rollback()