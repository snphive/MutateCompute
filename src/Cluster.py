import src.Main

# The job.q is a script that includes all necessary information for the grid engine, in terms of which computations
# to perform as well as how to run these computations

# MEMORY ###############################################################################################################
# If you don't specify anything, your job will be able to use 1 GB of memory. If you need more, you need to specify it:
# for e.g. if pass max_memory_GB='4' to build_job_q_bash(), it will set the maximum memory to 4 GB.

# OPTIONS DEFINITELY AVAILABLE ON CLUSTER ##############################################################################
# https://psbweb05.psb.ugent.be/clusterwiki/doku.php?id=general:gridengine and those used by Rob van der Kant
# -N <job name>             Name of job.
# -V                        Specifies to pass the environment variables to the job.
# -q <queue>                Sets the queue. Without this, it will use the standard queue.
# -o <output_logfile>       Name of output log file
# -pe serial <n_slots>      Specifies parallel environment and number of slots you want to use on a single node.
# -l mem_free               Needed together with -pe option for it to select a node that currently has memory available
#                           This value should be the TOTAL amount of memory that you EXPECT you job will need,
#                           BUT you must also specify h_vmem memory limit as memory limit/n_slots since OGE multiplies
#                           it by n_slots when you ask for multiple slots.
#                           (Memory / processor slot. So if specifying 2 slots, total memory'll be 2X hvmem_value)
# -l h_vmem=<size>          Max memory required (e.g. 3G or 3500M)
# -l hostname=<nodename>    Name of cluster node to use.
# -cwd                      Run in current working directory

# -XX:ParallelGCThreads=1   Use this if you want to re-enable Java's parallel garbage collector, but it was disabled
#                           because it takes up a lot of memory.

# OPTIONS POSSIBLY AVAILABLE ON CLUSTER  ###############################################################################
# bioinformatics.mdc-berlin.de/intro2UnixandSGE/sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html

# -w e                          Verifies options and abort if there is an error
# -m ea                     Will send email when job ends or aborts.
# -M <emailaddress>         Email address to send email to.
# -l h_rt=<hh:mm:ss>        Maximum run time (hours, minutes and seconds)
# -l s_rt=<hh:mm:ss>        Soft run time limit - NB: both s_rt and h_rt must be set.
# -wd <dir>                 Set working directory for this job as <dir>


class Cluster(object):

    # Named arguments are used here to allow default values to be set.
    # It means that the caller of this method must also pass named parameters to the method.
    # using_runscript          True/False using runscript (hence running FoldX)
    # python_script_with_path  Which script to run with python command to run it.
    # job_name                 Typically concatenation of mutant name and computation-specific prefix.
    # queue                    Which oge queue e.g. 'all.q' for all queues.
    # n_slots                  There are 8 slots (cpu cores) per cluster node.
    # max_memory_GB            Maximum memory in gigabytes to use e.g. 3 for 3 gigabytes.
    # cluster_node             Specific node on cluster to use e.g. hostname=hodor1.vib
    @staticmethod
    def build_job_q_bash(using_runscript, python_script_with_path, job_name, queue='', n_slots='',
                         total_memory_GB='', max_memory_GB='', cluster_node=''):
        job_q_file = open('./job.q', 'w')
        job_q_file.write('#!/bin/bash\n')
        job_q_file.write('#$ -N ' + job_name + '\n')
        job_q_file.write('#$ -V\n')

        if queue != '':
            job_q_file.write('#$ -q ' + queue + '\n')

        if n_slots != '' and total_memory_GB == '':
            raise ValueError('n_slots was specified but no total memory was given.')

        else:
            job_q_file.write('#$ -pe serial ' + n_slots + '\n')
            job_q_file.write('#$ -l mem_free=' + total_memory_GB + 'G\n')

        if max_memory_GB != '' and n_slots != '':
            max_memory_GB = int(max_memory_GB) / int(n_slots)
            job_q_file.write('#$ -l h_vmem=' + max_memory_GB + 'G\n')

        if cluster_node != '':
            job_q_file.write('#$ -l hostname=' + cluster_node + '\n')

        job_q_file.write('#$ -cwd\n')
        job_q_file.write('source ~/.bash_profile\n')

        if using_runscript:
            job_q_file.write(src.Main.path_zeus_FoldX_exe + ' -runfile runscript.txt\n')

        if python_script_with_path != '':
            job_q_file.write('python ' + python_script_with_path + '\n')

        job_q_file.close()

