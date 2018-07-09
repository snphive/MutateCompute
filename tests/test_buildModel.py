from unittest import TestCase
from src.FoldX import FoldX
from src.GeneralUtilityMethods import GUM

import mock


class TestBuildModel(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.absolute_path_tests = '/Users/u0120577/PycharmProjects/MutateCompute/tests'
        cls.abs_path_tests_inputs = cls.absolute_path_tests + '/Inputs'
        cls.abs_path_tests_inputs_fasta = cls.abs_path_tests_inputs + '/Fasta'
        cls.abs_path_tests_inputs_pdbs = cls.abs_path_tests_inputs + '/PDBs'
        cls.abs_path_tests_outputs = cls.absolute_path_tests + '/Outputs'
        cls.abs_path_tests_outputs_agadir = cls.abs_path_tests_outputs + '/Agadir'
        cls.abs_path_tests_outputs_fasta = cls.abs_path_tests_outputs + '/Fasta'
        cls.abs_path_tests_outputs_foldx = cls.abs_path_tests_outputs + '/FoldX'
        cls.abs_path_tests_outputs_foldx_buildmodel = cls.abs_path_tests_outputs_foldx + '/BuildModel'

    @classmethod
    def tearDownClass(cls):
        cls.absolute_path_tests = ''
        cls.abs_path_tests_inputs = ''
        cls.abs_path_tests_inputs_fasta = ''
        cls.abs_path_tests_inputs_pdbs = ''
        cls.abs_path_tests_outputs = ''
        cls.abs_path_tests_outputs_agadir = ''
        cls.abs_path_tests_outputs_fasta = ''
        cls.abs_path_tests_outputs_foldx = ''
        cls.abs_path_tests_outputs_foldx_buildmodel = ''

    def test_mutate_residues(self):
        input_pdb = 'RepairPDB_Test1Residue.pdb'
        write_wt_fasta_files = False
        mutant_aa_list = ['A', 'C', 'D']
        FoldX().BuildModel().mutate_residues_of_pdb(self.abs_path_inputs, self.abs_path_outputs, input_pdb, mutant_aa_list, write_wt_fasta_files)
        # with mock.patch('buildModel._build_list_of_mutant_names_for_foldx', return_value='list_of_mutant_names') as abc_urandom_function:
            # assert abc_urandom(5) == 'abcpumpkins'

    # Note FoldX expects mutant names to have format wtaa_chain_position_mutantaa
    def test_build_list_of_mutant_names_for_foldx(self):
        input_pdb = 'RepairPDB_Test1Residue.pdb'
        write_wt_fasta_files = False
        pdb_name_chain_fasta_dict = GUM.extract_pdb_name_fasta_chains_from_pdb(input_pdb, self.abs_path_tests_inputs,
                                                self.abs_path_tests_outputs_foldx_buildmodel, write_wt_fasta_files)
        mutant_aa_list = ['A', 'C', 'D']
        mutant_name_fx_list = FoldX().BuildModel()._build_list_of_mutant_names_for_foldx(mutant_aa_list, pdb_name_chain_fasta_dict)
        expected_mutant_name_list = ['RA1A', 'RA1C', 'RA1D', 'RB1A', 'RB1C', 'RB1D']
        self.assertEqual(mutant_name_fx_list, expected_mutant_name_list)
