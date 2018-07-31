import yaml
import os
import traceback
from src.Paths import Paths
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


class Cluster(object):

    with open("/Users/u0120577/PycharmProjects/MutateCompute/config/pathsAndDictionaries.yaml", 'r') as stream:

        try:
            paths_and_dictionaries = yaml.load(stream)
            path_zeus_FoldX_exe = paths_and_dictionaries['ROOT']['path_zeus_FoldX_exe']
            path_local_FoldX_exe = paths_and_dictionaries['ROOT']['path_local_FoldX_exe']

        except yaml.YAMLError as exc:
            print(exc)

    # Named arguments are used here to allow default values to be set.
    # It means that the caller of this method must also pass named parameters to the method.
    #
    # job_name                 N specifies job name, e.g. concatenation of mutant name + computation-specific prefix
    # using_runscript          True/False using runscript (hence running FoldX)
    # python_script_with_path  Which script to run with python command to run it
    # queue                    q specifies which oge queue you want to use, e.g. 'all.q' for all queues
    # n_slots                  serial is number of slots you want your job to use. There are 8 slots per cluster node.
    # total_memory_GB          mem_free is the total amount of memory (here as GB) you expect your job will need
    # memory_limit_GB          h_vmem is the max memory (here as GB) you want to allow your job to use
    # cluster_node             hostname specifies a specific node on the cluster you want to use e.g. hodor1.vib
    @staticmethod
    def write_job_q_bash(job_name, path_job_q_dest, using_runscript=True, python_script_with_paths='', queue='',
                         n_slots='', total_memory_GB='', memory_limit_GB='', cluster_node=''):
        if not os.path.exists(path_job_q_dest):
            traceback.print_exc()
            raise ValueError('Destination path for job.q is missing. It must be created first')
            return

        job_q = []
        job_q.append('#!/bin/bash\n'+'#$ -N '+job_name+'\n'+'#$ -V\n')

        if queue != '':
            job_q.append('#$ -q ' + queue + '\n')

        multicore_memory_command = ''
        if n_slots != '' and total_memory_GB == '':
            raise ValueError('n_slots was specified but no total memory was given.')
        elif n_slots != '' and total_memory_GB != '':
            multicore_memory_command = '#$ -pe serial ' + n_slots + ' '
            multicore_memory_command += '-l mem_free=' + total_memory_GB + 'G'

        if n_slots != '' and memory_limit_GB != '':
            memory_limit_GB = int(memory_limit_GB) / int(n_slots)

            if multicore_memory_command != '':
                multicore_memory_command += ',h_vmem=' + str(memory_limit_GB) + 'G\n'
                job_q.append(multicore_memory_command)
            else:
                job_q.append('#$ -l h_vmem=' + str(memory_limit_GB) + 'G\n')

        if cluster_node != '':
            job_q.append('#$ -l hostname=' + cluster_node + '\n')

        job_q.append('#$ -cwd\n' + 'source ~/.bash_profile\n')

        if using_runscript:
            job_q.append(Paths.PATH_ZEUS_FOLDX_EXE + ' -runfile runscript.txt\n')

        if python_script_with_paths != '':
            job_q.append('python ' + python_script_with_paths + '\n')

        with open(path_job_q_dest + '/job.q', 'w') as job_q_file:
            job_q_str = ''.join(job_q)
            job_q_file.write(job_q_str)

        return job_q_str
