import mysql.connector
from src.enums.UNPW import Server


class DBCnx(object):

    cnx_mysql_practice = ''
    cnx_mysql_snpeffect5 = ''

    def connect_mysql_practice_db(self):
        global cnx

        try:
            # cnx = mysql.connector.connect(user='root', password='K0yGrJ8(', host='127.0.0.1', database='mydb', port='3306')
            cnx = mysql.connector.connect(user=Server.PRACTICE_DB_MYSQL_USER.value,
                                          password=Server.PRACTICE_DB_MYSQL_PW.value,
                                          host=Server.PRACTICE_DB__MYSQL_NETADD.value,
                                          database=Server.PRACTICE_DB_MYSQL_DBNAME.value,
                                          port=Server.PRACTICE_DB__MYSQL_PORT.value)
            cur = cnx.cursor()
            cur.execute("CREATE TABLE testDBCon2 ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")
            # cnx.close()
        except mysql.connector.Error as err:
            print("Something went wrong: {}".format(err))

    def disconnect_mysql_practice_db(self):
        global cnx
        cnx.close()


# """
# Read ddG from csv file and write to database.
# :param path_output_bm_pdb_fxmutant_ddG_csvfile:
# """
# connection = mysql.connector.connect(user=Server.SNPEFFECT_V5_MYSQL_UN.value,
#                                      password=Server.SNPEFFECT_V5_MYSQL_PW.value,
#                                      host=Server.SNPEFFECT_V5_MYSQL_NETADD.value)
# , database='snpeffect_v5', port='3306')
# cursor = connection.cursor()
# cursor.execute("SHOW DATABASES")
# for x in cursor:
#     print(x)
# cursor.execute("CREATE DATABASE SnpEffect_v5.0")
# cursor.execute(
#     "CREATE TABLE testingDBConnection ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")