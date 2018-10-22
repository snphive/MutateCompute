import site
import sys
import multiprocessing as mp
import time

print('user site-packages is at: ' + site.getusersitepackages())
print('global site-packages is at: ' + ''.join(site.getsitepackages()))
print('python version is at: ' + str(sys.version_info[0]))

print('This computer has ' + str(mp.cpu_count()) + ' cpus')

start_time = time.perf_counter()
elapsed = time.perf_counter() - start_time
print('Time taken to complete iterating through fasta_list (after methods have been called): ' + str(elapsed))


# The following 4 lines of commented-out code successfully connect to my mysql database and create a table.
# import mysql.connector
# cnx = mysql.connector.connect(user='root', password='K0yGrJ8(', host='127.0.0.1', database='mydb', port='3306')
# cur = cnx.cursor()
# cur.execute("CREATE TABLE testingDBConnection ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")

# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)
# pydevd.stoptrace()


# Might be better to log in manually, as this can hang if mulitplexing is used.
# @staticmethod
# def log_in_cluster_ssh(*modules):
# HOST = "zeus.psb.ugent.be"
# # Ports are handled in ~/.ssh/config since we use OpenSSH
# COMMAND = "-p 7788 shazib"
# Ports are handled in ~/.ssh/config since we use OpenSSH
# ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND], shell=False, stdout=subprocess.PIPE,
# stderr=subprocess.PIPE)
# result = ssh.stdout.readlines()
# if result == []:
#     error = ssh.stderr.readlines()
#     print("ERROR: %s" % error, file=sys.stderr)
# else:
#     print(result)

# or use fabric??
# import os
# from fabric.api import *
#
# HOST="www.example.org"
# COMMAND="uname -a"
#
# env.user = os.getenv('SSH_USER', 'vagrant')
# env.password = os.getenv('SSH_PASSWORD', 'vagrant')
#
# @hosts(HOST)
# def do_something():
#     run(COMMAND)