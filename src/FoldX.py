import subprocess
import os
from src.GeneralUtilityMethods import GUM
from src.Cluster import Cluster
from src.Paths import Paths


class FoldX(object):

    def __init__(self):
        print('FoldX constructor')

    class Repair(object):

        def __init__(self):
            print('Repair constructor')

        def do_repair(self, input_pdb):
            print('do repair method')

    class BuildModel(object):

        FXBM_jobname_prefix = 'FXBM_'

        def __init__(self):
            print('buildmodel constructor')

        # Mutate specified amino acids in this pdb to all listed amino acids using FoldX BuildModel,
        # FoldX uses a runscript file, which must be written here.
        #
        # pdb                   String   Input pdb to be mutated.
        # mutate_to_aa_list     List     List of amino acids that you want to mutate your pdb to.
        # write_wt_fastafiles   Boolean  True/False is you want to write the wild-type sequence of the input pdb out.
        def mutate_residues_of_pdb(self, pdb, mutate_to_aa_list, write_wt_fastafiles):
            pdbname = pdb.split('.')[0]
            path_input_PDBs_pdbname = GUM.create_dir_tree(Paths.MC_INPUT, 'PDBs', pdbname)
            path_runscript_dest = GUM.create_dir_tree(path_input_PDBs_pdbname, 'FX_BuildModel')
            pdbname_chain_fasta_dict = GUM.extract_pdbname_chain_fasta_from_pdbs(pdb, Paths.MC_INPUT,
                                                                                 write_wt_fastafiles, Paths.MC_OUTPUT)
            fx_mutant_name_list = self._make_fx_mutant_name_list(mutate_to_aa_list, pdbname_chain_fasta_dict)

            try:
                os.makedirs(path_runscript_dest)
            except FileExistsError:
                print('Part of all of path already exists. This is absolutely fine.')

            action = '<BuildModel>#,individual_list.txt'
            GUM.write_runscript_for_pdbs(path_runscript_dest, 'RepairPDB_' + pdb, action)

            for fx_mutant_name in fx_mutant_name_list:
                path_jobq_indivlist_dest = GUM.create_dir_tree(path_input_PDBs_pdbname, fx_mutant_name)

                try:
                    os.makedirs(path_jobq_indivlist_dest)
                except FileExistsError:
                    print('Part of all of path already exists. This is absolutely fine.')

                try:
                    os.chdir(path_jobq_indivlist_dest)
                except ValueError:
                    print('directory tree was not created')

                self._write_individual_list_for_mutant(fx_mutant_name, path_jobq_indivlist_dest)
                # path_foldx = self.path_zeus_foldx_exe if use_cluster else self.path_local_foldx_exe
                job_name = self.FXBM_jobname_prefix + fx_mutant_name
                Cluster.write_job_q_bash(job_name, path_jobq_indivlist_dest)

                if os.path.exists(os.path.join(path_runscript_dest, 'runscript.txt')):
                    subprocess.call('qsub job.q', shell=True)
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


