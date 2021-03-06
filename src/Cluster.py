#!/usr/bin/env python3
"""
Class for writing job.q bash and to run a job using the job.q bash script.

(Includes a nested Enum of cluster-specific strings located at the end of the class.)

Currently (2018) the cluster is made up of:
odin1: 64 cores, split over 4 cpus
odin2: 64 cores, split over 4 cpus
hodor1: 44 cores, split over 2 cpus
hodor2: 44 cores, split over 2 cpus
BUT hodor1 & 2 are intel cpus which support Hyperthreading, which means that you have 2 slots per core
odin1 & 2 are amd cpus and they don't have this feature. (based on Frederik Delaere email 09Oct2018)
"""

import os
from src.enums.Str import Str
import subprocess
from src.enums.Paths import Paths
import time
import signal

# The job.q is a script that includes all necessary information for the grid engine, in terms of which computations
# to perform as well as how to run these computations

# MEMORY ###############################################################################################################
# If you don't specify anything, your job will be able to use 1 GB of memory. If you need more, you need to specify it:
# for e.g. if pass max_memory_GB='4' to build_job_q_bash(), it will set the maximum memory to 4 GB.

# OPTIONS DEFINITELY AVAILABLE ON CLUSTER ##############################################################################
# https://psbweb05.psb.ugent.be/clusterwiki/doku.php?id=general:gridengine and those used by Rob van der Kant
#
# -N <job name>             Name of job.
# -V                        Specifies to pass the environment variables to the job.
# -q <queue>                Sets the queue. Without this, it will use the standard queue.
# -o <output_logfile>       Name of output log file
#
# -pe serial <n_slots>      Specifies parallel environment and number of slots you want to use on a single node.
# -pe serial <min>-<max>    With this syntax, the cluster will search for a host with a min 4 and max 8 slots free.
# -l mem_free               This value should be the TOTAL amount of memory that you EXPECT your job will need,
#                           mem_free is required when selecting multiple slots with the -pe option. This allows it to
#                           select a node that currently has memory available.
#                           BUT you must also specify memory limit (h_vmem) as well.
# -l h_vmem=<size>          Memory limited for your job (e.g. 3G or 3500M). NOTE: you must divide it by the number of
#                           slots because OGE multiplies this parameter by n_slots when you ask for multiple slots.
#                           This is written into the method below.
#                           E.g. If specifying 2 slots, and you want to set your memory limit to 2 GB, then you must
#                           pass in the h_vmem value as 1, because OGE will multiply 2 (slots) by 1 (GB) -> 2 GB.
#
# -l hostname=<nodename>    Name of cluster node to use.
# -cwd                      Run in current working directory
#
# -XX:ParallelGCThreads=1   Use this if you want to re-enable Java's parallel garbage collector, but it was disabled
#                           because it takes up a lot of memory.

# OPTIONS POSSIBLY AVAILABLE ON CLUSTER  ###############################################################################
# bioinformatics.mdc-berlin.de/intro2UnixandSGE/sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html
#
# -w e                      Verifies options and abort if there is an error
# -m ea                     Will send email when job ends or aborts.
# -M <emailaddress>         Email address to send email to.
# -l h_rt=<hh:mm:ss>        Maximum run time (hours, minutes and seconds)
# -l s_rt=<hh:mm:ss>        Soft run time limit - NB: both s_rt and h_rt must be set.
# -wd <dir>                 Set working directory for this job as <dir>

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class Cluster(object):


    @staticmethod
    def write_job_q_bash(jobname: str, path_job_q_dir: str, startnum='', endnum='', using_runscript=False, path_runscript_dir='',
                         python_script_with_paths='', queue='', n_slots='', total_memory_GB='', memory_limit_GB='',
                         cluster_node='', max_runtime_limit_hard='00:00:30', max_runtime_limit_soft='00:00:01',
                         email_address='shahin.zibaee@kuleuven.vib.be'):
        """
        Note that only 1st 2 arguments are required, the rest have default named arguments that can be overwritten.
        If you  want to run parallel jobs using shared-memory to the specified number of slots, you must change replace
        'serial' with 'smp'.

        :param jobname: N specifies job name, e.g. concatenation of mutant name + computation-specific prefix.
        :param path_job_q_dir: Name of destination dir for this job.q file. Root fixed to /config/cluster_jobq.
        :param startnum: Starting number for an array job. If either startnum or endnum are empty, no array job.
        :param endnum: End number for an array job. If either startnum or endnum are empty, no array job.
        :param using_runscript: True/False using runscript (hence running FoldX).
        :param path_runscript_dir: Which script to run followed by paths to relevant executables, such as to Qsub.
        :param python_script_with_paths:
        :param queue: q specifies which oge queue you want to use, e.g. 'all.q' for all queues.
        :param n_slots: Serial is number of slots you want your job to use. There are 8 slots per cluster node.
        :param total_memory_GB: mem_free is the total amount of memory (here as GB) you expect your job will need.
        :param memory_limit_GB: h_vmem is the max memory (here as GB) you want to allow your job to use.
        :param cluster_node: Hostname specifies a specific node on the cluster you want to use e.g. hodor1.vib.
        :param max_runtime_limit_hard: Maximum runtime that any single job is allowed to run. After this time the OGE aborts the
        job using the SIGKILL signal.
        :param max_runtime_limit_soft: Maximum runtime that any single job will run before the OGE warns the job by sending the
        job the SIGUSR1 signal. (This behaviour is effective only if the job catches and handles that warning signal.)
        :param email_address: To which emails will be sent if send_email_job_end is True.
        :return:
        """
        try:
            os.makedirs(path_job_q_dir)
        except FileExistsError:
            print(Str.PARTALLPATHEXISTS_MSG.value)
        job_q = []
        job_q.append('#!/bin/bash' + Str.NEWLN.value)
        job_q.append('#$ -N' + Str.SPCE.value + jobname + '' + Str.NEWLN.value)
        if startnum != '' and endnum != '':
            job_q.append('#$ -t' + Str.SPCE.value + str(startnum) + '-' + str(endnum) + Str.NEWLN.value)
        job_q.append('#$ -V' + Str.NEWLN.value)
        if queue != '':
            job_q.append('#$ -q' + Str.SPCE.value + queue + Str.NEWLN.value)

        multicore_memory_command = ''
        if n_slots != '' and total_memory_GB == '':
            raise ValueError('n_slots was specified but no total memory was given.')
        elif n_slots != '' and total_memory_GB != '':
            multicore_memory_command = '#$ -pe serial' + Str.SPCE.value + n_slots + Str.SPCE.value
            multicore_memory_command += '-l mem_free=' + total_memory_GB + 'G'

        if n_slots != '' and memory_limit_GB != '':
            memory_limit_GB = int(memory_limit_GB) / int(n_slots)

            if multicore_memory_command != '':
                multicore_memory_command += ',h_vmem=' + str(memory_limit_GB) + 'G' + Str.NEWLN.value
                job_q.append(multicore_memory_command)
            else:
                job_q.append('#$ -l h_vmem=' + str(memory_limit_GB) + 'G' + Str.NEWLN.value)

        # job_q.append('#$ -l h_rt=' + max_runtime_limit_hard + Str.NEWLN.value)
        # job_q.append('#$ -l s_rt=' + max_runtime_limit_soft + Str.NEWLN.value)

        if cluster_node != '':
            job_q.append('#$ -l hostname=' + cluster_node + Str.NEWLN.value)
        if 'agadir' in python_script_with_paths.lower():
            # job_q.append('#$ -wd /switchlab/group/shazib/SnpEffect/cluster_logfiles' + Str.NEWLN.value)
            job_q.append('#$ -wd' + Str.SPCE.value + Paths.SE_CLSTR_LOGFILES.value + Str.NEWLN.value)

        # else:
        #     job_q.append('#$ -wd' + Str.SPCE.value + path_dst_dir + Str.NEWLN.value)
        else:
            job_q.append('#$ -cwd' + Str.NEWLN.value)
        job_q.append('source ~/.bash_profile' + Str.NEWLN.value)

        # if email_address != '':
        #     job_q.append('#$ -m ea' + Str.NEWLN.value)
        #     job_q.append('#S -M ' + email_address + Str.NEWLN.value)
        if using_runscript:
            job_q.append(Paths.ZEUS_FOLDX_EXE.value + Str.SPCE.value + Cluster.CLSTR.RNFL.value + Str.SPCE.value +
                         os.path.join(path_runscript_dir, Cluster.CLSTR.RNSCPT_TXT.value) + Str.NEWLN.value)

        if python_script_with_paths != '':
            job_q.append('python3' + Str.SPCE.value + python_script_with_paths + Str.NEWLN.value)
            # not sure if path to qsub is also needed, as it is present in Rob's scripts (OptProt has path/to/python
            # script then space then the qsub path string, which seemed odd as the path/to/qsub is specified in the
            # prefix to qsub command in subprocess.call())

        with open(os.path.join(path_job_q_dir, Cluster.CLSTR.JOBQFILE.value), 'w') as job_q_file:
            job_q_str = ''.join(job_q)
            job_q_file.write(job_q_str)

        return job_q_str

    @staticmethod
    def signal_handler(signum, frame, jobname: str):
        print('got SIGUSR1 for job: ' + jobname + '. This is just a warning that the job has exceeded the soft runtime limit '
                                                  'of 1 second. It will not terminate, unless it exceeds the hard runtime limit '
                                                  'of 30 seconds'.format(signum))

    @staticmethod
    def run_job_q(path_job_q_dir: str, run_in_background=False):
        """
        :param path_job_q_dir: Absolute path of the directory holding the job.q file.
        :param run_in_background: True to run detached from the terminal.
        """
        signal.signal(signal.SIGUSR1, Cluster.signal_handler)
        nohup = ''
        ampersand = ''
        if run_in_background:
            nohup = 'nohup' + Str.SPCE.value
            ampersand = Str.SPCE.value + '&'
        # try:
        cmd = nohup + Paths.ZEUS_QSUB_EXE.value + Cluster.CLSTR.QSUB.value + Str.SPCE.value + os.path.join(path_job_q_dir,
                                                                                    Cluster.CLSTR.JOBQFILE.value) + ampersand
        # except:
        #     signal.signal(signal.SIGUSR1, Cluster.signal_term_handler)

        subprocess.call(cmd, shell=True)

    @staticmethod
    def wait_for_grid_engine_job_to_complete(grid_engine_job_prefix_or_full_name: str):
        """
        :param grid_engine_job_prefix_or_full_name: Typically two characters and underscore to be concatenated to a unique job
        identifier. The prefix is likely an abbreviation of the type of computation being performed, e.g. 'BM_pdb123' for FoldX
        BuildModel of 'pdb123', RM_ for the removing unnecessary input and output files after FoldX job run.
        """
        print('Cluster.wait_for_grid_ending..() is called.......')
        check_qstat = subprocess.Popen(Cluster.CLSTR.QSTAT.value, stdout=subprocess.PIPE)
        output_qstat = check_qstat.stdout.read()
        # output_qstat = output_qstat.decode('utf-8')
        while grid_engine_job_prefix_or_full_name.encode('utf-8') in output_qstat:
            print('Waiting ' + str(Cluster.CLSTR.JOB_WAIT_SECS.value) + 'secs for ' + grid_engine_job_prefix_or_full_name +
                  ' to finish.')
            time.sleep(Cluster.CLSTR.JOB_WAIT_SECS.value)
            check_qstat = subprocess.Popen(Cluster.CLSTR.QSTAT.value, stdout=subprocess.PIPE)
            output_qstat = check_qstat.stdout.read()

    from enum import Enum

    class CLSTR(Enum):
        ALLQ = 'all' + Str.Q_EXT.value
        NODE_HDR1 = 'hodor1.vib'
        NODE_HDR2 = 'hodor2.vib'
        NODE_ODN1 = 'odin1.vib'
        NODE_ODN2 = 'odin2.vib'
        JOBQFILE = 'job' + Str.Q_EXT.value
        QSTAT = 'qstat'
        QSUB = 'qsub'
        RNFL = '-runfile'
        RNSCPT_TXT = Str.RNSCRPT.value + Str.TXTEXT.value
        JOB_WAIT_SECS = 5

