from src.MutateFasta import MutateFasta
from src.Agadir import Agadir
from src.FoldX import FoldX
from threading import Thread


class Scheduler(object):

    @staticmethod
    def start(operations, path_input, input_pdb_list, input_fasta_list, mutant_aa_list):
        for input_fasta in input_fasta_list:
            if operations['do_mutate_fasta'] == 'True':
                mutate_fasta = MutateFasta()
                thread = Thread(target=mutate_fasta.mutate_every_residue_fastas,
                                args=(path_input, input_fasta_list, mutant_aa_list))
                thread.start()
                thread.join()
            if operations['do_agadir']:
                agadir = Agadir()
                thread = Thread(target=agadir.compute, args=input_fasta)
                thread.start()
                thread.join()

        for input_pdb in input_pdb_list:
            if operations['do_foldx_repair']:
                repair = FoldX().Repair()
                repair.do_repair(input_pdb)

            if operations['do_foldx_buildmodel']:
                write_wt_fasta_files = True
                FoldX.BuildModel().mutate_residues_of_pdb(input_pdb, mutant_aa_list, write_wt_fasta_files)

            # if operations['do_foldx_stability']:
            #     stability = FoldX().Stability()
            #     stability.do_stuff(input_pdb)
            # if operations['do_foldx_analysecomplex']:
            #     analyseComplex = FoldX().AnalyseComplex()
            #     analyseComplex.do_stuff(input_pdb)

