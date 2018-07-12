import subprocess
import os
from src.GeneralUtilityMethods import GUM


class FoldX(object):

    def __init__(self):
        print('FoldX constructor')

    class Repair(object):

        def __init__(self):
            print('Repair constructor')

        def do_repair(self):
            print('do repair method')

    class BuildModel(object):

        def __init__(self, path_zeus_foldx_exe, path_local_foldx_exe):
            print('buildmodel constructor')
            self.foldx_buildmodel_job_prefix = 'FXBM_'
            self.path_zeus_foldx_exe = path_zeus_foldx_exe
            self.path_local_foldx_exe = path_local_foldx_exe

        # Mutate specified amino acids in this pdb to all listed amino acids using FoldX BuildModel,
        # FoldX uses a runscript file, which must be written here.
        # abs_path_inputs          Absolute path of Inputs folder.
        # abs_path_outputs         Absolute path of Outputs folder.
        # input_pdb                Input pdb to be mutated.
        # mutate_to_aa_list        List of amino acids that you want to mutate your pdb to.
        # write_wt_fasta_files     True/False is you want to write the wild-type sequence of the input pdb out.
        def mutate_residues_of_pdb(self, abs_path_inputs, abs_path_outputs, input_pdb, mutate_to_aa_list, write_wt_fasta_files):
            pdb_name_chain_fasta_dict = GUM.extract_pdb_name_fasta_chains_from_pdb(input_pdb, abs_path_inputs, abs_path_outputs, write_wt_fasta_files)
            mutant_name_fx_list = self._build_list_of_mutant_names_for_foldx(mutate_to_aa_list, pdb_name_chain_fasta_dict)
            self._make_runscript_for_FoldX_BuildModel(input_pdb)
            abs_path_outputs_FX_BM_pdb = self._create_foldx_buildmodel_pdb_dir_tree(abs_path_outputs, input_pdb)
            for mutant_name_fx in mutant_name_fx_list:
                abs_path_outputs_FX_BM_pdb_mutant = GUM.create_directory_tree(abs_path_outputs_FX_BM_pdb, mutant_name_fx)
                try:
                    os.chdir(abs_path_outputs_FX_BM_pdb_mutant)
                except ValueError:
                    print('directory tree was not created')
                self._write_individual_list_for_mutant(mutant_name_fx, abs_path_outputs_FX_BM_pdb_mutant)
                no_python_script = ''
                use_cluster = True
                path_foldx = self.path_zeus_foldx_exe if use_cluster else self.path_local_foldx_exe
                if os.path.isfile('./runscript.txt'):
                    self._execute_job(use_cluster, mutant_name_fx, no_python_script, path_foldx)
                else:
                    raise ValueError('No runscript file was found')

        def _create_foldx_buildmodel_pdb_dir_tree(self, abs_path_outputs, pdb):
            return GUM.create_directory_tree(abs_path_outputs, 'FoldX', 'BuildModel', pdb.split('.')[0])

        def _execute_job(self, use_cluster, mutant_name_fx, python_script_with_path, path_foldx):
            job_name = self.foldx_buildmodel_job_prefix + mutant_name_fx
            using_runscript = False
            if use_cluster:
                no_specific_queue = ''
                no_max_memory = ''
                no_specific_cluster_node = ''
                using_runscript = True

            GUM.build_job_q_bash(use_cluster, job_name, no_specific_queue, no_max_memory, no_specific_cluster_node,
                                     using_runscript, path_foldx, python_script_with_path)
            subprocess.call('qsub job.q', shell=True)

        # Note FoldX expects mutant names to have format wtaa_chain_position_mutantaa.
        # Variable names for values with this format are denoted with '_fx_'
        def _build_list_of_mutant_names_for_foldx(self, mutate_to_aa_list, pdb_name_chain_fasta_dict):
            mutant_name_fx_list = []
            for pdb_name_chain, fasta_sequence in pdb_name_chain_fasta_dict.items():
                chain = pdb_name_chain.split('_')[-1]
                for index, wt_aa in enumerate(fasta_sequence):
                    position = index + 1
                    for mutant_aa in mutate_to_aa_list:
                        mutant_name_fx_list.append(wt_aa + chain + str(position) + mutant_aa)
            return mutant_name_fx_list

        # Note FoldX expects the mutant_name to have the form: WTaa_chain_position_mutantaa (no gaps or underscores)
        def _make_runscript_for_FoldX_BuildModel(self, name_of_repair_pdb, abs_path_runscript):
            abs_path_runscript = './'
            repair_pdbs = 'RepairPDB_' + name_of_repair_pdb + '.pdb'
            show_sequence_detail = False
            action = '<BuildModel>#,individual_list.txt'
            print_networks = False
            calculate_stability = False
            GUM.write_runscript_for_pdbs(abs_path_runscript, repair_pdbs, show_sequence_detail, action, print_networks,
                                                          calculate_stability)

        def _write_individual_list_for_mutant(self, mutant_name, absolute_path_output_dir):
            individual_list_for_this_mutant_only = open(absolute_path_output_dir + 'individual_list.txt', 'w')
            individual_list_for_this_mutant_only.write(mutant_name + ';\n')
            individual_list_for_this_mutant_only.close()

    class Stability(object):

        def __init__(self):
            print('helloworld constructor')

    class AnalyseComplex(object):

        def __init__(self):
            print('helloworld constructor')

        def _prepare_for_FoldX_AnalyseComplex(self, repair_pdb_name):
            _0_1_2_pdbs = ['0.pdb,', '1.pdb,', '2.pdb,']
            repair_pdb_name_1_ = repair_pdb_name + '_1_'
            wt_repair_pdb_name_1_ = 'WT_' + repair_pdb_name_1_

            path_to_runscript = './'
            pdbs_to_analyse = repair_pdb_name_1_ + _0_1_2_pdbs[0] + \
                              repair_pdb_name_1_ + _0_1_2_pdbs[1] + \
                              repair_pdb_name_1_ + _0_1_2_pdbs[2] + \
                              wt_repair_pdb_name_1_ + _0_1_2_pdbs[0] + \
                              wt_repair_pdb_name_1_ + _0_1_2_pdbs[1] + \
                              wt_repair_pdb_name_1_ + _0_1_2_pdbs[2]
            show_sequence_detail = False
            action = '<AnalyseComplex>#'
            print_networks = False
            calculate_stability = False
            GUM.write_runscript_for_pdbs(path_to_runscript, pdbs_to_analyse,
                                         show_sequence_detail, action, print_networks,
                                         calculate_stability)


