#!/usr/bin/env python3
"""
Class for running program(s) on each file from a list of files. This is where multithreading code is intended to be
(though not up and running).
"""

import os
import time
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing as mp
# import concurrent.futures.ThreadPoolExecutor
from threading import Thread
from src.MutateFasta import MutateFasta
from src.Agadir import Agadir
from src.FoldX import FoldX
from src.enums.Conditions import Cond
from src.IdentifyProtein import IdProt
from src.tools.GeneralUtilityMethods import GUM
from src.Cluster import Cluster
from src.enums.Paths import Paths
from src.enums.Str import Str

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class Scheduler(object):

    @staticmethod
    def start(operations: dict, use_multithread: bool, path_input: str, path_output: str, path_pdbfiles: list,
              path_fastafiles: list, specific_fxmutants: list, amino_acids: list, write_1_fasta_only: bool,
              write_fasta_per_mut: bool, write_to_csv_dumpfile_after_each_mutant: bool):
        """
        Iterate through a list of fasta files or pdb files and perform Agadir, or Foldx computations as specified by 'operations'.
        :param operations: Each operation paired with True/False flag to indicate whether or not to perform the operation.
        :param use_multithread: True to employ parallel processing.
        :param path_input: Absolute path to input_data root dir.
        :param path_output: Absolute path to output_data root dir.
        :param path_pdbfiles: Absolute path to pdb input files.
        :param path_fastafiles: Absolute path to fasta input files.
        :param specific_fxmutants: Given when specific mutants only should be calculated.
        :param amino_acids: Amino acids that mutation operations should use to mutate to.
        :param write_1_fasta_only: True to write any fasta output data to 1 fasta file, each separated by \n.
        :param write_fasta_per_mut: True to write any fasta output data as 1 fasta file per mutant. (Uses a lot of disk space).
        :param write_to_csv_dumpfile_after_each_mutant: True to write ddG values from fxout files to one csv file (for database
        dump).
        """
        if path_fastafiles:
            if operations[Scheduler.Strs.OPER_RUN_MUT_FSTA.value]:
                path_output_fastas_3dots = GUM.make_path_fastas_3dots_dirs(path_output, path_fastafiles[0])
                mutate_fasta = MutateFasta(amino_acids)
                for path_fastafile in path_fastafiles:
                    sleep_secs = 0 if len(path_fastafiles) < 200 else len(path_fastafiles) / 5000
                    time.sleep(sleep_secs)
                    if use_multithread:
                        # Scheduler._launch_thread(target=mutate_fasta.mutate_every_residue,
                        #                          args=[path_fastafile, write_1_fasta_only, write_fasta_per_mut,
                        #                                path_output_3dots])
                        Scheduler._launch_process(target=mutate_fasta.mutate_every_residue,
                                                  args=[path_fastafile, write_1_fasta_only, write_fasta_per_mut,
                                                        path_output_fastas_3dots])
                    elif not GUM.using_cluster():
                        mutate_fasta.mutate_every_residue(path_fastafile, write_1_fasta_only, write_fasta_per_mut,
                                                              path_output_fastas_3dots)
                    if GUM.using_cluster():
                        jobname = Paths.PREFIX_MUTFSTA.value + path_fastafile.split('/')[-1]
                        write_1_fasta_only = True
                        write_fasta_per_mut = False
                        Cluster.write_job_q_bash(jobname=jobname, path_job_q_dir=Paths.SE_CONFIG_MUTFASTA_JOBQ.value,
                                                 python_script_with_paths=os.path.join(Paths.SE_SRC.value,
                                                 'run_mutate_fasta_zeus.py') +
                                                 Str.SPCE.value + path_fastafile +
                                                 Str.SPCE.value + str(write_1_fasta_only) +
                                                 Str.SPCE.value + str(write_fasta_per_mut) +
                                                 Str.SPCE.value + path_output_fastas_3dots, queue='',
                                                 n_slots='', total_memory_GB='', memory_limit_GB='3', cluster_node='')
                        Cluster.run_job_q(path_job_q_dir=Paths.SE_CONFIG_MUTFASTA_JOBQ.value)

            if operations[Scheduler.Strs.OPER_RUN_AGDR.value]:
                agadir = Agadir(Cond.INCELL_MAML_FX.value)
                for path_fastafile in path_fastafiles:
                    sleep_secs = 0 if len(path_fastafiles) < 200 else len(path_fastafiles) / 1000
                    time.sleep(sleep_secs)
                    if GUM.using_cluster():
                        print('Calling scheduler.do_agadir using_cluster condition')
                        jobname = Paths.PREFIX_AGADIR.value + path_fastafile.split('/')[-1]
                        Cluster.write_job_q_bash(jobname=jobname, path_job_q_dir=Paths.SE_CONFIG_AGAD_JOBQ.value,
                                                 python_script_with_paths=os.path.join(Paths.SE_SRC.value,
                                                 'run_agadir_on_multifastas_zeus.py' +
                                                 Str.SPCE.value + path_fastafile +
                                                 Str.SPCE.value + Paths.SE_OUTPUT.value))
                        Cluster.run_job_q(path_job_q_dir=Paths.SE_CONFIG_AGAD_JOBQ.value)

                    path_dst = GUM.make_path_agadir_3dots_filename_mutants_dirs(path_output, path_fastafile,
                                                                                add_filename_subdir=True)
                    if use_multithread:
                        # Scheduler._launch_thread(target=agadir.run_agadir_on_multifastas,
                        #                          args=[path_fastafile, path_output])
                        Scheduler._launch_process(target=agadir.run_agadir_on_multifastas,
                                                  args=[path_fastafile, path_dst])
                    elif not GUM.using_cluster() and not use_multithread:
                        agadir.run_agadir_on_multifastas(path_fastafile, path_dst)
        if path_pdbfiles:
            for path_pdbfile in path_pdbfiles:
                if operations[Scheduler.Strs.OPER_RUN_FX_BM.value]:
                    buildmodel = FoldX().BuildModel(Cond.INCELL_MAML_FX.value)
                    if use_multithread:
                        Scheduler._launch_thread(target=buildmodel.mutate_protein_structure,
                                                 args=[path_pdbfile, amino_acids, specific_fxmutants])
                    else:
                        buildmodel.mutate_protein_structure(path_pdbfile, amino_acids, specific_fxmutants,
                                                write_to_csv_dumpfile_after_each_mutant=write_to_csv_dumpfile_after_each_mutant)
                if operations[Scheduler.Strs.OPER_RUN_FX_AC.value]:
                    analysecomplex = FoldX().AnalyseComplex(Cond.INCELL_MAML_FX.value)
                    if use_multithread:
                        Scheduler._launch_thread(target=analysecomplex.calculate_complex_energies, args=path_pdbfile)
                    else:
                        analysecomplex.calculate_complex_energies(path_pdbfile, specific_fxmutants,
                                                write_to_csv_dumpfile_after_each_mutant=write_to_csv_dumpfile_after_each_mutant)
                if operations[Scheduler.Strs.OPER_RUN_FX_RPR.value]:
                    repair = FoldX().Repair(Cond.INCELL_MAML_FX.value)
                    if use_multithread:
                        Scheduler._launch_thread(target=repair.do_repair, args=path_pdbfile)
                    else:
                        repair.do_repair(path_pdbfile)

        # if path_pdbfiles:
        #     for path_pdbfile in path_pdbfiles:
        #         if operations['do_foldx_buildmodel']:
        #             buildmodel = FoldX().BuildModel(Cond.INCELL_MAML_FX.value)
        #             if buildmodel.find_num_of_missing_avg_bm_fxoutfiles(path_pdbfile, amino_acids) != 0:
        #                 print('Warning: BuildModel has not completed all computations for this pdb: ' + os.path.basename(
        #                     path_pdbfile))
        #             else:
        #                 print('BuildModel has completed all computations for this pdb: ' + os.path.basename(path_pdbfile))

                # if operations['do_foldx_analysecomplex']:
                #     analysecomplex = FoldX().AnalyseComplex(Cond.INCELL_MAML_FX.value)
                #     if analysecomplex.find_num_of_missing_interaction_energies(path_pdbfile, amino_acids) != 0:
                #         print('Warning: AnalyseComplex has not completed all computations for this pdb: ' + os.path.basename(
                #             path_pdbfile))
                #     else:
                #         print('AnalyseComplex has completed all computations for this pdb: ' + os.path.basename(path_pdbfile))


    @staticmethod
    def start_blast(path_input_fastafiles: list, path_output: str, write_idmaps_for_mysldb=True, write_csv=True, write_xml=True,
                    write_json=False, use_multithread=False):
        """

        :param path_input_fastafiles:
        :param path_output:
        :param write_idmaps_for_mysldb:
        :param write_csv:
        :param write_xml:
        :param write_json:
        :param use_multithread:
        :return:
        """
        if use_multithread:
            print('')
            # Scheduler._launch_thread(target=IdProt.map_seq_to_swsprt_acc_id_and_write_files,
            #                          args=[path_input_fastas_dir, path_output, write_idmaps_for_mysldb, write_csv,
            #                                write_xml, write_json])

            # Scheduler._launch_process_threadpool(target=IdProt., args=)

        else:
            IdProt.map_seq_to_swsprt_acc_id_and_write_files(path_input_fastafiles=path_input_fastafiles,
                                                            path_output=path_output,
                                                            write_idmaps_for_mysqldb=write_idmaps_for_mysldb,
                                                            write_csv=write_csv, write_xml=write_xml,
                                                            write_json=write_json)

    @staticmethod
    def _launch_thread(target: callable(object), args: list, in_background=False):
        """

        :param target:
        :param args:
        :param in_background:
        :return:
        """
        t = Thread(target=target, args=args)
        t.daemon = in_background
        t.start()
        t.join()

    @staticmethod
    def _launch_process(target: callable(object), args: list):
        """

        :param target:
        :param args:
        :return:
        """
        p = mp.Process(target=target, args=args)
        p.start()
        p.join()

    # no pool size given to ThreadPool(x), defaults to number of CPUs (which is 4 on my macbook).
    @staticmethod
    def _launch_process_threadpool(methodname: str, arg_list: list):
        """

        :param methodname:
        :param arg_list:
        :return:
        """
        pool = ThreadPool(4)
        results = pool.map(func=methodname, iterable=arg_list)
        pool.close()
        pool.join()
        return results
    # Scheduler._launch_process_threadpool(mutate_fasta.mutate_every_residue_to_every_aa_write_1_file,
    #                                          path_fastafiles)

    # Note: try pool.imap and pool.map_async in place of pool.map

    from enum import Enum

    class Strs(Enum):
        OPER_RUN_MUT_FSTA = Str.OPER_RUN_MUT_FSTA.value
        OPER_RUN_AGDR = Str.OPER_RUN_AGDR.value
        OPER_RUN_FX_RPR = Str.OPER_RUN_FX_RPR.value
        OPER_RUN_FX_BM = Str.OPER_RUN_FX_BM.value
        OPER_RUN_FX_AC = Str.OPER_RUN_FX_AC.value
        OPER_RUN_FX_STAB = Str.OPER_RUN_FX_STAB.value
