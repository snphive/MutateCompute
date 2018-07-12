from unittest import TestCase
from src.FoldX import FoldX
from src.GeneralUtilityMethods import GUM
import mock
from mock import patch


class TestBuildModel(TestCase):

    def setUp(self):
        print('bla')

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

    # @patch('FoldX.BuildModel.')
    def test_execute_job(self):
        # mock_GUM =
        input_pdb = 'RepairPDB_Test1Residue.pdb'
        write_wt_fasta_files = False
        mutant_aa_list = ['A', 'C', 'D']
        FoldX().BuildModel().mutate_residues_of_pdb(self.abs_path_tests_inputs, self.abs_path_tests_outputs, input_pdb,
                                                    mutant_aa_list, write_wt_fasta_files)

    # Note FoldX expects mutant names to have format wtaa_chain_position_mutantaa.
    # Variable names for values with this format are denoted with '_fx_'
    def test_build_list_of_mutant_names_for_foldx(self):
        pdb_name_chain_fasta_dict = {'Test1Residue_A': 'R', 'Test1Residue_B': 'R'}
        mutant_aa_list = ['A', 'C', 'D']
        mutant_name_fx_list = FoldX().BuildModel()._build_list_of_mutant_names_for_foldx(mutant_aa_list, pdb_name_chain_fasta_dict)
        expected_mutant_name_fx_list = ['RA1A', 'RA1C', 'RA1D', 'RB1A', 'RB1C', 'RB1D']
        self.assertEqual(mutant_name_fx_list, expected_mutant_name_fx_list)
