from src.MutateFasta import MutateFasta
from src.Agadir import Agadir
from src.FoldX import FoldX
from threading import Thread
import time
import multiprocessing as mp


class Scheduler(object):

    # Not sure yet if I should instantiate the objects once at beginning rather than a new instance for every and every
    # fasta or pdb
    #
    # operations            Dictionary
    # path_input            String
    # pdb_list              List
    # fastafile_list        List
    # mutant_aa_list        List
    # use_multithread       Boolean
    # write_1_fasta_only    Boolean
    # write_fasta_per_mut   Boolean
    # path_output           String
    @staticmethod
    def start(operations, path_input, pdb_list, fastafile_list, mutant_aa_list, use_multithread, write_1_fasta_only,
              write_fasta_per_mut, path_output):
        print('This computer has ' + str(mp.cpu_count()) + ' cpus')
        start_time = time.perf_counter()
        mutate_fasta = MutateFasta()
        for fastafile in fastafile_list:
            if operations['do_mutate_fasta']:
                if use_multithread:
                    Scheduler._launch_thread(target=mutate_fasta.mutate_every_residue,
                                             args=[path_input, fastafile_list, fastafile, mutant_aa_list,
                                                   write_1_fasta_only, write_fasta_per_mut, path_output])
                    # Scheduler._launch_process(target=mutate_fasta.mutate_every_residue_in_fasta_list,
                    #                          args=[path_input, fasta_list, mutant_aa_list])
                else:
                    mutate_fasta.mutate_every_residue(path_input, fastafile_list, fastafile, mutant_aa_list,
                                                      write_1_fasta_only, write_fasta_per_mut, path_output)
            if operations['do_agadir']:
                agadir = Agadir()
                if use_multithread:
                    Scheduler._launch_thread(target=agadir.compute, args=fastafile)
                else:
                    agadir.compute(fastafile)
        elapsed = time.perf_counter() - start_time
        print('Time taken to complete iterating through fasta_list (after methods have been called): ' + str(elapsed))
        for pdb in pdb_list:
            if operations['do_foldx_repair']:
                repair = FoldX().Repair()
                if use_multithread:
                    Scheduler._launch_thread(target=repair.do_repair, args=pdb)
                else:
                    repair.do_repair(pdb)
            if operations['do_foldx_buildmodel']:
                write_wt_fastafiles = True
                buildmodel = FoldX().BuildModel()
                if use_multithread:
                    Scheduler._launch_thread(target=buildmodel.mutate_residues_of_pdb,
                                             args=[pdb, mutant_aa_list, write_wt_fastafiles])
                else:
                    buildmodel.mutate_residues_of_pdb(pdb, mutant_aa_list, write_wt_fastafiles)
            if operations['do_foldx_stability']:
                stability = FoldX().Stability()
                if use_multithread:
                    Scheduler._launch_thread(target=stability.calculate_stability, args=pdb)
                else:
                    stability.calculate_stability(pdb)
            if operations['do_foldx_analysecomplex']:
                analysecomplex = FoldX().AnalyseComplex()
                if use_multithread:
                    Scheduler._launch_thread(target=analysecomplex.calculate_complex_energies, args=pdb)
                else:
                    analysecomplex.calculate_complex_energies(pdb)

    @staticmethod
    def _launch_thread(target, args):
        t = Thread(target=target, args=args)
        t.start()
        t.join()

    @staticmethod
    def _launch_process(target, args):
        p = mp.Process(target=target, args=args)
        p.start()
        p.join()
