import subprocess
import os
from src.Str import Str
from src.Paths import Paths
from src.GeneralUtilityMethods import GUM
from src.Cluster import Cluster


class FoldX(object):

    def __init__(self):
        print('FoldX constructor')

    class Repair(object):

        def __init__(self):
            print('Repair constructor')

        def do_repair(self, input_pdb):
            print('do repair method')

    class BuildModel(object):

        def __init__(self):
            print('buildmodel constructor')

        # Mutate specified amino acids in this pdb to all listed amino acids using FoldX BuildModel,
        # FoldX uses a runscript file, which must be written here.
        #
        # path_pdbfile                   String   Input pdb to be mutated.
        # mutant_aa_list        List     List of amino acids that you want to mutate your pdb to.
        # write_wt_fastafiles   Boolean  True/False is you want to write the wild-type sequence of the input pdb out.
        def mutate_all_pdb(self, path_input, path_pdbfile, mutant_aa_list, path_output):
            pdbname = path_pdbfile.split('/')[-1].split('.')[0]
            path_input_pdbs_pdbname = GUM._os_makedirs(path_input, Paths.DIR_PDBS.value, pdbname)
            path_runscript = GUM._os_makedirs(path_input_pdbs_pdbname, Paths.DIR_BMRUNSCRIPT)
            pdbname_chain_fasta_dict = GUM.extract_pdbname_chain_fasta_from_pdbs(path_input, path_pdbfile,
                                                                                 path_output, )
            fx_mutant_name_list = self._make_fx_mutant_name_list(mutant_aa_list, pdbname_chain_fasta_dict)

            action = '<BuildModel>#,individual_list.txt'
            GUM.write_runscript_for_pdbs(path_runscript, path_pdbfile, action)

            for fx_mutant_name in fx_mutant_name_list:
                path_jobq_indivlist_dest = GUM._os_makedirs(path_input_pdbs_pdbname, fx_mutant_name)
                self._write_individual_list_for_mutant(fx_mutant_name, path_jobq_indivlist_dest)
                job_name = Paths.PREFIX_FX_BM.value + fx_mutant_name
                Cluster.write_job_q_bash(job_name, path_job_q_dir=path_jobq_indivlist_dest,
                                         path_runscript_dir=path_runscript)
                if os.path.exists(os.path.join(path_runscript, 'runscript.txt')):
                    Cluster.run_job_q(path_job_q_dir=path_jobq_indivlist_dest)
                else:
                    raise ValueError('No runscript file was found')

        def _make_fx_mutant_name_list(self, mutate_to_aa_list, pdbname_chain_fasta_dict):
            fx_mutant_name_list = []
            for pdbname_chain, fasta_sequence in pdbname_chain_fasta_dict.items():
                chain = pdbname_chain.split('_')[-1]
                for index, wt_aa in enumerate(fasta_sequence):
                    position = index + 1
                    for mutant_aa in mutate_to_aa_list:
                        fx_mutant_name_list.append(wt_aa + chain + str(position) + mutant_aa)
            return fx_mutant_name_list

        # Not tested yet.
        def _write_individual_list_for_mutant(self, fx_mutant_name, path_indivlist_dest):
            path_indivlist_dest_txtfile = os.path.join(path_indivlist_dest, 'individual_list.txt')
            with open(path_indivlist_dest_txtfile, 'w') as indiv_list_file_opened:
                indiv_list_file_opened.write(fx_mutant_name + ';\n')

    class Stability(object):

        def __init__(self):
            print('helloworld constructor')

        def calculate_stability(self, pdb):
            print('not implemented yet')

    class AnalyseComplex(object):

        def __init__(self):
            print('helloworld constructor')

        def calculate_complex_energies(self, pdbname):
            self._prepare_for_FoldX_AnalyseComplex(pdbname)
            print('not yet implemented')

        def _prepare_for_FoldX_AnalyseComplex(self, repair_pdbname):
            _0_1_2_pdbs = ['0.pdb,', '1.pdb,', '2.pdb,']
            repair_pdbname_1_ = repair_pdbname + '_1_'
            wt_repair_pdbname_1_ = 'WT_' + repair_pdbname_1_

            path_to_runscript = './'
            pdbs_to_analyse = repair_pdbname_1_ + _0_1_2_pdbs[0] + \
                              repair_pdbname_1_ + _0_1_2_pdbs[1] + \
                              repair_pdbname_1_ + _0_1_2_pdbs[2] + \
                              wt_repair_pdbname_1_ + _0_1_2_pdbs[0] + \
                              wt_repair_pdbname_1_ + _0_1_2_pdbs[1] + \
                              wt_repair_pdbname_1_ + _0_1_2_pdbs[2]
            action = '<AnalyseComplex>#'
            GUM.write_runscript_for_pdbs(path_to_runscript, pdbs_to_analyse, action)
