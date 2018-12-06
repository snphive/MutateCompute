import mysql.connector
from src.enums.DBServer import Server
"""
Class for connecting to a relational database (currently only MySQL set up). 
Class also performs SQL DDL commands: e.g. CREATE, ALTER, DROP.
"""

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class DBCnxDDL(object):

    def connect_mysql_practice_db(self):
        """
        This method works. The table was created in
        :return:
        """
        try:
            mysql_practice_cnx = mysql.connector.connect(user=Server.PRACTICE_DB_MYSQL_USER.value,
                                                         password=Server.PRACTICE_DB_MYSQL_PW.value,
                                                         host=Server.PRACTICE_DB_MYSQL_NETADD.value,
                                                         database=Server.PRACTICE_DB_MYSQL_DBNAME.value,
                                                         port=Server.PRACTICE_DB_MYSQL_PORT.value)
            cursor = mysql_practice_cnx.cursor()
            cursor.execute("CREATE TABLE IF NOT EXISTS fx_data ( pdb VARCHAR(200) NOT NULL PRIMARY KEY, mutation VARCHAR(20))")
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))
            mysql_practice_cnx.rollback()
        return mysql_practice_cnx

    def add_new_columns(self, mysql_practice_cnx, column_name: str, data_type_defin: str):
        cursor = mysql_practice_cnx.cursor()
        try:
            cursor.execute("ALTER TABLE fx_data ADD COLUMN ddG FLOAT AFTER mutation")
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))
            mysql_practice_cnx.rollback()

    def disconnect_mysql_db(self, mysql_practice_cnx):
        print("Disconnecting mysql database" + str(mysql_practice_cnx))
        mysql_practice_cnx.close()

    """
    CONNECTING AND PERFORMING SQL ON SNP_EFFECT_V5
    """
    def connect_mysql_snpeffect_v5_db(self, mysql_snpeffect_v5_cnx):
        try:
            mysql_snpeffect_v5_cnx = mysql.connector.connect(user=Server.SNPEFFECT_V5_MYSQL_USER.value,
                                                             password=Server.SNPEFFECT_V5_MYSQL_PW.value,
                                                             host=Server.SNPEFFECT_V5_MYSQL_NETADD.value,
                                                             database=Server.SNPEFFECT_V5_MYSQL_DBNAME.value,
                                                             port=Server.SNPEFFECT_V5_MYSQL_PORT.value)
            cursor = mysql_snpeffect_v5_cnx.cursor()
            cursor.execute("CREATE TABLE testDBConx ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")

        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))
            mysql_snpeffect_v5_cnx.rollback()