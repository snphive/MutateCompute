from unittest import TestCase
from src.FoldX import FoldX
from src.GeneralUtilityMethods import GUM
from src.Cluster import Cluster
from unittest.mock import patch
from tests.testHelperMethods import THM


# Note: A test pdb is used here which is taken from the RepairPDBs folder, but includes only the first 5 residues of
# to facilitate manual creation of test variables, such as the list of all mutants.
# Note: Wherever a path is indicated, it is the absolute path. Relative paths are always prefixed with rel_
class TestBuildModel(TestCase):

    @classmethod
    def setUpClass(cls):
        # paths common to both input & output
        cls.path_tests = '/Users/u0120577/PycharmProjects/MutateCompute/tests'
        cls.rel_path_PDBs = '/PDBs'
        cls.rel_path_FoldX = '/FoldX'
        # input paths only
        cls.rel_path_Inputs = '/Inputs'
        cls.path_tests_Inputs = cls.path_tests + cls.rel_path_Inputs
        cls.path_tests_Inputs_PDBs = cls.path_tests_Inputs + cls.rel_path_PDBs
        cls.rel_path_Cluster = '/Cluster'
        cls.rel_path_BuildModel = '/BuildModel'
        cls.rel_path_Fasta = '/Fasta'
        cls.rel_path_Options_Agadir = '/Options/Agadir'
        cls.rel_path_Options_FoldX = '/Options/FoldX'
        cls.rel_path_Options_Cluster = '/Options/Cluster'
        # output paths only
        cls.rel_path_Outputs = '/Outputs'
        cls.path_tests_Outputs = cls.path_tests + cls.rel_path_Outputs
        cls.rel_path_BuildModel = '/BuildModel'
        cls.rel_path_AnalyseComplex = '/AnalyseComplex'

    def setUp(self):
        # THM.remove_tests_Inputs_Outputs_folders()
        self.path_zeus_foldx_exe = '/switchlab/group/tools/FoldX_2015/FoldX'
        self.path_local_foldx_exe = '/Users/u0120577/SNPEFFECT/executables/FoldX'
        foldx = FoldX()
        self.buildModel = foldx.BuildModel(self.path_zeus_foldx_exe, self.path_local_foldx_exe)

    def tearDown(self):
        self.buildModel = None
        # THM.remove_tests_Inputs_Outputs_folders()

    @patch('subprocess.call')
    @patch.object(Cluster, 'write_job_q_bash')
    @patch.object(FoldX.BuildModel, '_write_individual_list_for_mutant')
    @patch('os.chdir')
    @patch.object(GUM, 'create_dir_tree')
    @patch.object(GUM, 'write_runscript_for_pdbs')
    @patch.object(FoldX.BuildModel, '_make_fx_mutant_name_list')
    @patch.object(GUM, 'extract_pdbname_chain_fasta_from_pdb')
    def test_mutate_residues_of_pdb(self,
                                    mock_extract_pdbname_chain_fasta_from_pdb,
                                    mock_make_fx_mutant_name_list,
                                    mock_write_runscript_for_pdbs,
                                    mock_create_dir_tree,
                                    mock_os_chdir,
                                    mock_write_individual_list_for_mutant,
                                    mock_write_job_q_bash,
                                    mock_subprocess_call):
        # arrange
        mock_extract_pdbname_chain_fasta_from_pdb.return_value = {'RepairPDB_1_A': 'RVYLT', 'RepairPDB_1_B': 'RVYLT'}
        # Note FoldX expects mutant names to have format wtaa_chain_position_mutantaa.
        mock_make_fx_mutant_name_list.return_value = ['RA1A', 'RA1C', 'VA2A', 'VA2C', 'YA3A', 'YA3C', 'LA4A',
                                                          'LA4C', 'TA5A', 'TA5C', 'RB1A', 'RB1C', 'VB2A', 'VB2C',
                                                          'YB3A', 'YB3C', 'LB4A', 'LB4C', 'TB5A', 'TB5C']
        pdb = 'RepairPDB_1.pdb'

        mock_os_chdir.return_value = None
        write_wt_fasta_files = False
        mutant_aa_list = ['A', 'C']
        # action
        self.buildModel.mutate_residues_of_pdb(self.path_tests_Inputs, self.path_tests_Outputs, pdb, mutant_aa_list,
                                               write_wt_fasta_files)
        # assert
        expected_call_count_1 = 1
        expected_call_count_20 = len(mock_make_fx_mutant_name_list.return_value)
        expected_call_count_22 = 22
        self.assertEqual(expected_call_count_1, mock_extract_pdbname_chain_fasta_from_pdb.call_count)
        self.assertEqual(expected_call_count_1, mock_make_fx_mutant_name_list.call_count)
        self.assertEqual(expected_call_count_1, mock_write_runscript_for_pdbs.call_count)
        self.assertEqual(expected_call_count_22, mock_create_dir_tree.call_count)
        self.assertEqual(expected_call_count_20, mock_write_individual_list_for_mutant.call_count)
        self.assertEqual(expected_call_count_20, mock_write_job_q_bash.call_count)
        self.assertEqual(expected_call_count_20, mock_subprocess_call.call_count)

    # Note FoldX expects mutant names to have format wtaa_chain_position_mutantaa.
    # Variable names for values with this format are denoted with '_fx_'
    def test_make_fx_mutant_name_list(self):
        pdbname_chain_fasta_dict = {'Test1Residue_A': 'R', 'Test1Residue_B': 'R'}
        mutant_aa_list = ['A', 'C', 'D']
        fx_mutant_name_list = self.buildModel._make_fx_mutant_name_list(mutant_aa_list, pdbname_chain_fasta_dict)
        expected_fx_mutant_name_list = ['RA1A', 'RA1C', 'RA1D', 'RB1A', 'RB1C', 'RB1D']
        self.assertEqual(fx_mutant_name_list, expected_fx_mutant_name_list)



    # def test_write_job_q_filenames_check(self):
        # job_q_file = open('./job.q', 'w') - how to assert filenames are same?
        # job_q_file.write('#!/bin/bash\n') - not sure what this  does

    # unused file comparison which relies to some extent on the directory structure being in place which I didn't know how to mock
    # so I made the write_job.. method return the string that it writes to job.q file.
    #     with open(path_job_q_file + '/job.q') as f:
    #         actual_job_q = f.read()
    #         self.assertEqual(actual_job_q, expected_job_q)
    #         self.assertNotEqual(actual_job_q, not_expected_job_q)
    #         self.assertNotEqual(actual_job_q, not_expected_job_q_2)

    # pdbname = pdb.split('.')[0]
    # path_tests_Inputs_PDBs_pdbname = self.path_tests_Inputs_PDBs + '/' + pdbname
    # path_tests_Inputs_PDBs_pdbname_FX_BM = path_tests_Inputs_PDBs_pdbname + '/FX_BuildModel'
    # path_tests_Inputs_PDBs_pdbname_FXmutantname = path_tests_Inputs_PDBs_pdbname + '/RA1A'
    # mock_os_path_join.side_effect = [path_tests_Inputs_PDBs_pdbname, path_tests_Inputs_PDBs_pdbname_FX_BM,
    #                                  path_tests_Inputs_PDBs_pdbname_FXmutantname]
