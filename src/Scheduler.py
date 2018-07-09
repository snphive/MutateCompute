from src.MutateFasta import MutateFasta
from src.Agadir import Agadir
from src.FoldX import FoldX
from threading import Thread
from time import sleep


class Scheduler(object):

    def __init__(self, absolute_path_inputs, absolute_path_outputs):
        self.abs_path_inputs = absolute_path_inputs
        self.abs_path_outputs = absolute_path_outputs

    def start(self, operations, input_pdbs, input_fastas, list_of_all_20_aa):

        for input_fasta in input_fastas:
            if operations['do_mutate_fasta']:
                mutate_fasta = MutateFasta()
                mutate_all_residues = True
                thread = Thread(target=mutate_fasta.create_mutant_fastas, args=(input_fasta, mutate_all_residues,
                                                                                list_of_all_20_aa))
                thread.start()
                thread.join()
            if operations['do_agadir']:
                agadir = Agadir()
                thread = Thread(target=agadir.compute, args=input_fasta)
                thread.start()
                thread.join()

        for input_pdb in input_pdbs:
            if operations['do_foldx_repair']:
                FoldX().Repair(input_pdb)

            if operations['do_foldx_buildmodel']:
                # Code will be adapted so user can mutate to specific residues to, (rather than all or none)
                mutate_all_residues = True
                if mutate_all_residues:
                    mutate_to_residues = list_of_all_20_aa
                write_wt_fasta_files = True
                FoldX().BuildModel().mutate_residues_of_pdb(self.abs_path_inputs, self.abs_path_outputs, input_pdb,
                                                            mutate_to_residues, write_wt_fasta_files)
            if operations['do_foldx_stability']:
                FoldX().Stability(input_pdb)
            if operations['do_foldx_analysecomplex']:
                FoldX().AnalyseComplex(input_pdb)

