import mysql.connector
from src.enums.DBServer import Server
"""
Class for connecting to database and performing certain sql actions.
"""


class DBCnx(object):

    def connect_mysql_practice_db(self):
        """
        This method works. The table was created in
        :return:
        """

        global mysql_practice
        try:
            mysql_practice = mysql.connector.connect(user=Server.PRACTICE_DB_MYSQL_USER.value,
                                                     password=Server.PRACTICE_DB_MYSQL_PW.value,
                                                     host=Server.PRACTICE_DB_MYSQL_NETADD.value,
                                                     database=Server.PRACTICE_DB_MYSQL_DBNAME.value,
                                                     port=Server.PRACTICE_DB_MYSQL_PORT.value)
            cursor = mysql_practice.cursor()
            cursor.execute("CREATE TABLE IF NOT EXISTS fx_data ( pdb VARCHAR(200) NOT NULL PRIMARY KEY, mutation VARCHAR(20) "
                           "NOT NULL )")
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))
            mysql_practice.rollback()
        # finally:
        #     mysql_practice.close()

    def add_new_columns(self, column_name: str, data_type_defin: str):
        global mysql_practice
        cursor = mysql_practice.cursor()
        try:
            cursor.execute("ALTER TABLE fx_data ADD COLUMN ddG FLOAT AFTER mutation")
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))
            mysql_practice.rollback()

    def insert_values_by_column(self, columns: list, values: list):
        global mysql_practice
        cursor = mysql_practice.cursor()
        sql_cmd = ' INSERT INTO fx_data(pdb, mutation, ddG) VALUES(\"RepairPDB_5\", \"AA101A\", \"0.00023\") '
        try:
            # columns = ','.join(columns)
            # values = ','.join([str(value) for value in values])
            # values = [('"' + value + '"') for value in values])
            # sql_cmd = f'INSERT INTO testForScope.testDBCnx ({columns}) VALUES ({values});'
            # cursor.execute(sql_cmd)

            # abc = 'RepairPDB_3'
            # sql = f'INSERT INTO testForScope.testDBCnx (pdb) VALUES ({abc})'
            # sql_insery = """ sql """
            # cursor.execute(sql_insery)
            # pdbname =
            # cursor.execute("SHOW DATABASES")
            # for x in cur:
            #     print(x)
            # cursor.execute("SELECT * FROM ")
            cursor.execute(sql_cmd)
            mysql_practice.commit()

        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))
            mysql_practice.rollback()

    def disconnect_mysql_practice_db(self):
        global mysql_practice
        print("Disconnecting mysql practice db")
        mysql_practice.close()

    """
    CONNECTING AND PERFORMING SQL ON SNP_EFFECT_V5
    """
    def connect_mysql_snpeffect_v5_db(self):
        global mysql_snpeffect5
        try:
            mysql_snpeffect5 = mysql.connector.connect(user=Server.SNPEFFECT_V5_MYSQL_USER.value,
                                                       password=Server.SNPEFFECT_V5_MYSQL_PW.value,
                                                       host=Server.SNPEFFECT_V5_MYSQL_NETADD.value,
                                                       database=Server.SNPEFFECT_V5_MYSQL_DBNAME.value,
                                                       port=Server.SNPEFFECT_V5_MYSQL_PORT.value)
            cursor = mysql_snpeffect5.cursor()
            cursor.execute("CREATE TABLE testDBConx ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")

        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))
            mysql_snpeffect5.rollback()
        # finally:
        #     mysql_snpeffect5.close()

    def disconnect_mysql_snpeffect_v5_db(self):
        global mysql_snpeffect5
        print("Disconnecting mysql snpeffect_v5")
        mysql_snpeffect5.close()
