from unittest import TestCase
from src.FoldX import FoldX
from src.Cluster import Cluster
from src.GeneralUtilityMethods import GUM
from unittest.mock import patch


class TestBuildModel(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.absolute_path_tests = '/Users/u0120577/PycharmProjects/MutateCompute/tests'
        cls.abs_path_tests_inputs = cls.absolute_path_tests + '/Inputs'
        cls.abs_path_tests_inputs_foldx = cls.abs_path_tests_inputs + '/FoldX'
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
        cls.abs_path_tests_inputs_foldx = ''
        cls.abs_path_tests_inputs_fasta = ''
        cls.abs_path_tests_inputs_pdbs = ''
        cls.abs_path_tests_outputs = ''
        cls.abs_path_tests_outputs_agadir = ''
        cls.abs_path_tests_outputs_fasta = ''
        cls.abs_path_tests_outputs_foldx = ''
        cls.abs_path_tests_outputs_foldx_buildmodel = ''

    def setUp(self):
        self.path_zeus_foldx_exe = '/switchlab/group/tools/FoldX_2015/FoldX'
        self.path_local_foldx_exe = '/Users/u0120577/SNPEFFECT/executables/FoldX'
        # foldx = FoldX(path_zeus_foldx_exe, path_local_foldx_exe)
        foldx = FoldX()
        self.buildModel = foldx.BuildModel(self.path_zeus_foldx_exe, self.path_local_foldx_exe)

    def tearDown(self):
        self.buildModel = None

    @patch.object(FoldX.BuildModel, '_write_job_q_bash_to_run_on_cluster_using_runscript')
    @patch.object(FoldX.BuildModel, '_write_individual_list_for_mutant')
    @patch('os.chdir')
    @patch.object(FoldX.BuildModel, '_create_foldx_buildmodel_pdb_dir_tree')
    @patch.object(FoldX.BuildModel, '_write_runscript_for_FoldX_BuildModel')
    @patch.object(FoldX.BuildModel, '_make_list_of_mutant_names_for_foldx')
    @patch.object(GUM, 'extract_pdb_name_fasta_chains_from_pdb')
    def test_mutate_residues_of_pdb(self, mock_extract_pdb_name_fasta_chains_from_pdb,
                                    mock_make_list_of_mutant_names_for_foldx,
                                    mock_write_runscript_for_FoldX_BuildModel,
                                    mock_create_foldx_buildmodel_pdb_dir_tree,
                                    mock_os_chdir,
                                    mock_write_individual_list_for_mutant,
                                    mock__write_job_q_bash_to_run_on_cluster_using_runscript):
        # arrange
        mock_extract_pdb_name_fasta_chains_from_pdb.return_value = {'RepairPDB_Test1Residue_A': 'R',
                                                                    'RepairPDB_Test1Residue_B': 'R'}
        mock_make_list_of_mutant_names_for_foldx.return_value = ['RA1A', 'RA1C', 'RA1D', 'RB1A', 'RB1C', 'RB1D']
        input_pdb = 'RepairPDB_Test1Residue.pdb'
        input_pdb_name = input_pdb.split('_')[-1].split('.')[0]
        mock_create_foldx_buildmodel_pdb_dir_tree.return_value = self.abs_path_tests_outputs_foldx_buildmodel + '/' + \
                                                                 input_pdb_name
        mock_os_chdir.return_value = None
        write_wt_fasta_files = False
        mutant_aa_list = ['A', 'C', 'D']
        # action
        self.buildModel.mutate_residues_of_pdb(self.abs_path_tests_inputs, self.abs_path_tests_outputs, input_pdb,
                                               mutant_aa_list, write_wt_fasta_files)
        # assert
        expected_call_count_1 = 1
        expected_call_count_6 = 6
        self.assertEqual(expected_call_count_1, mock_extract_pdb_name_fasta_chains_from_pdb.call_count)
        self.assertEqual(expected_call_count_1, mock_make_list_of_mutant_names_for_foldx.call_count)
        self.assertEqual(expected_call_count_1, mock_write_runscript_for_FoldX_BuildModel.call_count)
        self.assertEqual(expected_call_count_1, mock_create_foldx_buildmodel_pdb_dir_tree.call_count)
        self.assertEqual(expected_call_count_6, mock_write_individual_list_for_mutant.call_count)
        self.assertEqual(expected_call_count_6, mock__write_job_q_bash_to_run_on_cluster_using_runscript.call_count)

    # Note FoldX expects mutant names to have format wtaa_chain_position_mutantaa.
    # Variable names for values with this format are denoted with '_fx_'
    def test_make_list_of_mutant_names_for_foldx(self):
        pdb_name_chain_fasta_dict = {'Test1Residue_A': 'R', 'Test1Residue_B': 'R'}
        mutant_aa_list = ['A', 'C', 'D']
        mutant_name_fx_list = self.buildModel._make_list_of_mutant_names_for_foldx(mutant_aa_list, pdb_name_chain_fasta_dict)
        expected_mutant_name_fx_list = ['RA1A', 'RA1C', 'RA1D', 'RB1A', 'RB1C', 'RB1D']
        self.assertEqual(mutant_name_fx_list, expected_mutant_name_fx_list)

    def test_write_job_q_bash_to_run_on_cluster_using_runscript_no_queue_slots_memlimit_or_clusternode(self):
        # arrange
        foldx_buildmodel_job_prefix = 'FXBM_'
        mutant_name_fx = 'RA1A'
        job_name = foldx_buildmodel_job_prefix + mutant_name_fx
        python_script_with_path = '/somepath/to/python/script.py'
        abs_path_job_q_file = self.abs_path_tests_inputs_foldx
        expected_job_q = '#!/bin/bash\n' + '#$ -N ' + job_name + '\n' + '#$ -V\n' + '#$ -cwd\n' + \
                              'source ~/.bash_profile\n' + self.path_zeus_foldx_exe + ' -runfile runscript.txt\n' + \
                              'python ' + python_script_with_path + '\n'
        single_space = ' '
        not_expected_job_q = '#!/bin/bash' + single_space + '\n' + '#$ -N ' + job_name + '\n' + '#$ -V\n' + \
                             '#$ -cwd\n' + 'source ~/.bash_profile\n' + self.path_zeus_foldx_exe + \
                             ' -runfile runscript.txt\n' + 'python ' + python_script_with_path + '\n'
        missing_new_line = ''
        not_expected_job_q_2 = '#!/bin/bash\n' + '#$ -N ' + job_name + missing_new_line + '#$ -V\n' + \
                             '#$ -cwd\n' + 'source ~/.bash_profile\n' + self.path_zeus_foldx_exe + \
                             ' -runfile runscript.txt\n' + 'python ' + python_script_with_path + '\n'
        # action
        actual_job_q = self.buildModel._write_job_q_bash_to_run_on_cluster_using_runscript(mutant_name_fx,
                                                                        python_script_with_path, abs_path_job_q_file)
        # assert
        self.assertEqual(actual_job_q, expected_job_q)
        self.assertNotEqual(actual_job_q, not_expected_job_q)
        self.assertNotEqual(actual_job_q, not_expected_job_q_2)

    # def test_write_job_q_filenames_check(self):
        # job_q_file = open('./job.q', 'w') - how to assert filenames are same?
        # job_q_file.write('#!/bin/bash\n') - not sure what this  does

    # unused file comparison which relies to some extent on the directory structure being in place which I didn't know how to mock
    # so I made the write_job.. method return the string that it writes to job.q file.
    #     with open(abs_path_job_q_file + '/job.q') as f:
    #         actual_job_q = f.read()
    #         self.assertEqual(actual_job_q, expected_job_q)
    #         self.assertNotEqual(actual_job_q, not_expected_job_q)
    #         self.assertNotEqual(actual_job_q, not_expected_job_q_2)
