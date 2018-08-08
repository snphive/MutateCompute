from src.MutateFasta import MutateFasta
from src.Agadir import Agadir
from src.FoldX import FoldX
from threading import Thread


class Scheduler(object):

    @staticmethod
    def start(operations, path_input, pdb_list, fasta_list, mutant_aa_list):
        process_started = False
        for fasta in fasta_list:
            if operations['do_mutate_fasta']:
                mutate_fasta = MutateFasta()
                thread = Thread(target=mutate_fasta.mutate_every_residue_in_fasta_list,
                                args=(path_input, fasta_list, mutant_aa_list))
                thread.start()
                process_started = True
                thread.join()
            if operations['do_agadir']:
                agadir = Agadir()
                thread = Thread(target=agadir.compute, args=fasta)
                thread.start()
                process_started = True
                thread.join()

        for pdb in pdb_list:
            if operations['do_foldx_repair']:
                repair = FoldX().Repair()
                thread = Thread(target=repair.do_repair, args=pdb)
                thread.start()
                process_started = True
                thread.join()

            if operations['do_foldx_buildmodel']:
                write_wt_fastafiles = True
                buildmodel = FoldX().BuildModel()
                thread = Thread(target=buildmodel.mutate_residues_of_pdb,
                                args=(pdb, mutant_aa_list, write_wt_fastafiles))
                thread.start()
                process_started = True
                thread.join()

            if operations['do_foldx_stability']:
                stability = FoldX().Stability()
                thread = Thread(target=stability.calculate_stability, args=pdb)
                thread.start()
                process_started = True
                thread.join()

            if operations['do_foldx_analysecomplex']:
                analysecomplex = FoldX().AnalyseComplex()
                thread = Thread(target=analysecomplex.calculate_complex_energies, args=pdb)
                thread.start()
                process_started = True
                thread.join()

        return process_started
