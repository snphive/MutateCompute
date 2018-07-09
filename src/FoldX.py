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

        def __init__(self):
            print('buildmodel constructor')

        def mutate_residues_of_pdb(self, abs_path_inputs, abs_path_outputs, input_pdb, mutate_to_aa_list, write_wt_fasta_files):
            pdb_name_chain_fasta_dict = GUM.extract_pdb_name_fasta_chains_from_pdb(input_pdb, abs_path_inputs, abs_path_outputs, write_wt_fasta_files)
            mutant_names_list = self._build_list_of_mutant_names_for_foldx(mutate_to_aa_list, pdb_name_chain_fasta_dict)
            for mutant_name in mutant_names_list:
                self._make_runscript_FoldX_BuildModel(mutant_name, input_pdb)
                self._write_indivdual_list_for_mutant(mutant_name)

        def _build_list_of_mutant_names_for_foldx(self, mutate_to_aa_list, pdb_name_chain_fasta_dict):
            mutant_name_list = []

            for pdb_name_chain, fasta_sequence in pdb_name_chain_fasta_dict.items():
                chain = pdb_name_chain.split('_')[-1]

                for index, wt_aa in enumerate(fasta_sequence):
                    position = index + 1

                    for mutant_aa in mutate_to_aa_list:
                        mutant_name_list.append(wt_aa + chain + str(position) + mutant_aa)

            return mutant_name_list

        # Note that FoldX expects the mutant_name to have the form: WTaa_chain_position_mutantaa (no gaps or underscores)
        def _make_runscript_FoldX_BuildModel(self, mutant_name, name_of_repair_pdb):
            path_to_runscript = './'
            repair_pdbs = 'RepairPDB_' + name_of_repair_pdb + '.pdb'
            show_sequence_detail = False
            action = '<BuildModel>#,individual_list.txt'
            print_networks = False
            calculate_stability = False
            GUM.write_runscript_for_pdbs(path_to_runscript, repair_pdbs, show_sequence_detail, action, print_networks,
                                                          calculate_stability)

        def _write_individual_list_for_mutant(self, mutant_name):
            individual_list_for_this_mutant_only = open('individual_list.txt', 'w')
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
