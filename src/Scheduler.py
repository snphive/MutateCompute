import os
import time
from multiprocessing import Process
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing as mp
# import concurrent.futures.ThreadPoolExecutor
from threading import Thread
from src.MutateFasta import MutateFasta
from src.Agadir import Agadir
from src.Agadir import AgadCndtns
from src.FoldX import FoldX
from src.Conditions import Cond
from src.IdentifyProtein import IdProt
from src.GeneralUtilityMethods import GUM
from src.Cluster import Cluster
from src.Paths import Paths
from src.Str import Str


class Scheduler(object):

    # Not sure yet if I should instantiate the objects once at beginning rather than a new instance for every and every
    # fasta or pdb.
    # (NOTE: argument list is identical as those of the calling method _start_scheduler())
    #
    # operations            Dictionary  Operations each paired with a flag, True to run the operation.
    # use_multithread       Boolean     True to employ parallel processing.
    # path_input            String      Abs path to input_data root dir.
    # path_output           String      Abs path to output_data root dir.
    # path_pdbfile_list     List        Pdbfiles to run (incl. .pdb extension)
    # path_fastafile_list   List        Fastafiles to run (incl. .fasta extension)
    # mutant_aa_list        List        All amino acids that any mutations operations will use to mutate residues to.
    # write_1_fasta_only    Boolean     True to write any fasta output data to 1 fasta file, each separated by \n.
    # write_fasta_per_mut   Boolean     True to write any fasta output data as 1 fasta file per mutant.
    @staticmethod
    def start(operations, use_multithread, path_input, path_output, path_pdbfile_list, path_fastafile_list,
              mutant_aa_list, write_1_fasta_only, write_fasta_per_mut):
        print('STARTING SCHEDULER')
        start_time = time.perf_counter()
        if path_fastafile_list:
            if operations['do_mutate_fasta']:
                path_output_fastas_3dots = GUM.make_root_fastas_3dots_dirs(path_output, path_fastafile_list[0])
                mutate_fasta = MutateFasta(mutant_aa_list)
                for path_fastafile in path_fastafile_list:
                    sleep_secs = 0 if len(path_fastafile_list) < 200 else len(path_fastafile_list) / 1000
                    time.sleep(sleep_secs)
                    if use_multithread:
                        # Scheduler._launch_thread(target=mutate_fasta.mutate_every_residue,
                        #                          args=[path_fastafile, write_1_fasta_only, write_fasta_per_mut,
                        #                                path_output_3dots])
                        Scheduler._launch_process(target=mutate_fasta.mutate_every_residue,
                                                  args=[path_fastafile, write_1_fasta_only, write_fasta_per_mut,
                                                        path_output_fastas_3dots])
                    elif not use_multithread:
                        mutate_fasta.mutate_every_residue(path_fastafile, write_1_fasta_only, write_fasta_per_mut,
                                                              path_output_fastas_3dots)
                    if GUM.using_cluster():
                        jobname = path_fastafile.split('/')[-1]
                        # Cluster.write_job_q_bash(job_name=jobname, path_job_q_dir=Paths.SE_CONFIG_JOBQ.value,
                        #                          python_script_with_paths=Paths.ZEUS_SNPEFFECT +
                        #                           'src/run_MutCompZeus.py ' + path_fastafile_list, queue='',
                        #                           n_slots='', total_memory_GB='', memory_limit_GB='3',
                        #                           cluster_node='')
                        # Cluster.run_job_q(path_job_q_dir=Paths.SE_CONFIG_JOBQ.value)

            if operations['do_agadir']:
                agadir = Agadir(AgadCndtns.INCELL_MAML.value)
                for path_fastafile in path_fastafile_list:
                    sleep_secs = 0 if len(path_fastafile_list) < 200 else len(path_fastafile_list) / 1000
                    time.sleep(sleep_secs)
                    if GUM.using_cluster():
                        print('Calling scheduler.do_agadir using_cluster condition')
                        jobname = Paths.PREFIX_AGADIR.value + path_fastafile.split('/')[-1]
                        Cluster.write_job_q_bash(job_name=jobname, path_job_q_dir=Paths.SE_CONFIG_AGAD_JOBQ.value,
                                                 python_script_with_paths=os.path.join(Paths.SRC,
                'run_agadir_on_multifastasZeus.py' + Str.SPCE.value + path_fastafile + Str.SPCE.value + Paths.OUTPUT))
                        Cluster.run_job_q(path_job_q_dir=Paths.SE_CONFIG_AGAD_JOBQ.value)
                    if use_multithread:
                        # Scheduler._launch_thread(target=agadir.run_agadir_on_multifastas,
                        #                          args=[path_fastafile, path_output])
                        Scheduler._launch_process(target=agadir.run_agadir_on_multifastas,
                                                  args=[path_fastafile, path_output])
                    elif not GUM.using_cluster() and not use_multithread:
                        agadir.run_agadir_on_multifastas(path_fastafile, path_output)

        elapsed = time.perf_counter() - start_time
        print('Time taken to complete iterating through fasta_list (after methods have been called): ' + str(elapsed))
        if path_pdbfile_list:
            for path_pdbfile in path_pdbfile_list:
                if operations['do_foldx_repair']:
                    repair = FoldX().Repair(Cond.INCELL_MAML_FX.value)
                    if use_multithread:
                        Scheduler._launch_thread(target=repair.do_repair, args=path_pdbfile)
                    else:
                        repair.do_repair(path_pdbfile)
                if operations['do_foldx_buildmodel']:
                    buildmodel = FoldX().BuildModel(Cond.INCELL_MAML_FX.value)
                    if use_multithread:
                        Scheduler._launch_thread(target=buildmodel.mutate_protein_structure,
                                                 args=[path_pdbfile, mutant_aa_list])
                    else:
                        buildmodel.mutate_protein_structure(path_pdbfile, mutant_aa_list)
                if operations['do_foldx_stability']:
                    stability = FoldX().Stability(Cond.INCELL_MAML_FX.value)
                    if use_multithread:
                        Scheduler._launch_thread(target=stability.calculate_stability, args=path_pdbfile)
                    else:
                        stability.calculate_stability(path_pdbfile)
                if operations['do_foldx_analysecomplex']:
                    analysecomplex = FoldX().AnalyseComplex(Cond.INCELL_MAML_FX.value)
                    if use_multithread:
                        Scheduler._launch_thread(target=analysecomplex.calculate_complex_energies, args=path_pdbfile)
                    else:
                        analysecomplex.calculate_complex_energies(path_pdbfile)

    @staticmethod
    def start_blast(path_input_fastafiles, path_output, write_idmaps_for_mysldb=True, write_csv=True, write_xml=True,
                    write_json=False, use_multithread=False):
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
    def _launch_thread(target, args, in_background=False):
        t = Thread(target=target, args=args)
        t.daemon = in_background
        t.start()
        t.join()

    @staticmethod
    def _launch_process(target, args):
        p = mp.Process(target=target, args=args)
        p.start()
        p.join()

    # no pool size given to ThreadPool(x), defaults to number of CPUs (which is 4 on my macbook).
    @staticmethod
    def _launch_process_threadpool(methodname, arg_list):
        pool = ThreadPool(4)
        results = pool.map(func=methodname, iterable=arg_list)
        pool.close()
        pool.join()
        return results
    # Scheduler._launch_process_threadpool(mutate_fasta.mutate_every_residue_to_every_aa_write_1_file,
    #                                          path_fastafile_list)

    # Note: try pool.imap and pool.map_async in place of pool.map

