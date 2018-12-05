import mysql.connector
from src.enums.DBServer import Server


class DBCnx(object):

    def connect_mysql_practice_db(self):
        global mysql_practice
        try:
            # cnx = mysql.connector.connect(user='root', password='K0yGrJ8(', host='127.0.0.1', database='mydb', port='3306')
            mysql_practice = mysql.connector.connect(user=Server.PRACTICE_DB_MYSQL_USER.value,
                                                     password=Server.PRACTICE_DB_MYSQL_PW.value,
                                                     host=Server.PRACTICE_DB_MYSQL_NETADD.value,
                                                     database=Server.PRACTICE_DB_MYSQL_DBNAME.value,
                                                     port=Server.PRACTICE_DB_MYSQL_PORT.value)
            cur = mysql_practice.cursor()
            cur.execute("CREATE TABLE IF NOT EXISTS testDBCnx ( pdb VARCHAR(200) NOT NULL PRIMARY KEY, mutation VARCHAR(20) NOT "
                        "NULL )")
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))

    def add_new_columns(self, column_name: str, data_type_defin: str):
        global mysql_practice
        cur = mysql_practice.cursor()
        cur.execute("ALTER TABLE testDBCnx ADD COLUMN ddG FLOAT AFTER title")

    def insert_values(self):
        global mysql_practice
        cur = mysql_practice.cursor()
        cur.execute("INSERT INTO testDBCnx (pdb, mutation, ddG) VALUES ('RepairPDB_1', 'AA101A', 0.002)")
    #     may need to include schema name, e.g. testForScope.testDBCnx

    def disconnect_mysql_practice_db(self):
        global mysql_practice
        mysql_practice.close()

    def connect_mysql_snpeffect_v5_db(self):
        global mysql_snpeffect5
        try:
            mysql_snpeffect5 = mysql.connector.connect(user=Server.SNPEFFECT_V5_MYSQL_USER.value,
                                                       password=Server.SNPEFFECT_V5_MYSQL_PW.value,
                                                       host=Server.SNPEFFECT_V5_MYSQL_NETADD.value,
                                                       database=Server.SNPEFFECT_V5_MYSQL_DBNAME.value,
                                                       port=Server.SNPEFFECT_V5_MYSQL_PORT.value)
            cur = mysql_snpeffect5.cursor()
            cur.execute("CREATE TABLE testDBConx ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")

        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))

    def disconnect_mysql_snpeffect_v5_db(self):
        global mysql_snpeffect5
        mysql_snpeffect5.close()
