from unittest import TestCase
from src.GeneralUtilityMethods import GUM
from tests.testHelperMethods import THM
import os
from unittest.mock import patch


class TestGUM(TestCase):

    @classmethod
    def setUpClass(cls):
        # paths common to both input & output
        cls.path_tests = '/Users/u0120577/PycharmProjects/MutateCompute/tests'
        cls.path_rel_PDBs = '/PDBs'
        cls.path_rel_FoldX = '/FoldX'
        # input paths only
        cls.path_rel_Inputs = '/Inputs'
        cls.path_tests_Inputs = cls.path_tests + cls.path_rel_Inputs
        cls.path_tests_Inputs_PDBs = cls.path_tests_Inputs + cls.path_rel_PDBs
        cls.path_rel_Cluster = '/Cluster'
        cls.path_rel_BuildModel = '/BuildModel'
        cls.path_rel_Fasta = '/Fasta'
        cls.path_rel_Options_Agadir = '/Options/Agadir'
        cls.path_rel_Options_FoldX = '/Options/FoldX'
        cls.path_rel_Options_Cluster = '/Options/Cluster'
        # output paths only
        cls.path_rel_Outputs = '/Outputs'
        cls.path_tests_Outputs = cls.path_tests + cls.path_rel_Outputs
        cls.path_rel_BuildModel = '/BuildModel'
        cls.path_rel_AnalyseComplex = '/AnalyseComplex'

    # @classmethod
    # def tearDownClass(cls):
    #     cls.path_local_tests, cls.path_rel_Inputs, cls.path_rel_Cluster,cls.path_rel_BuildModel, cls.path_rel_Fasta,\
    #     cls.path_rel_Options_Agadir, cls.path_rel_Options_FoldX, cls.path_rel_Options_Cluster, cls.path_rel_Outputs, \
    #     cls.path_rel_Agadir, cls.path_rel_PDBs, cls.path_rel_FoldX, cls.path_rel_FoldX_BuildModel, \
    #     cls.path_rel_FoldX_AnalyseComplex = ''

    # def setUp(self):
    #     # THM.remove_tests_Inputs_Outputs_folders()
    #
    # def tearDown(self):
    #     # THM.remove_tests_Inputs_Outputs_folders()

    # def test_wait_for_grid_engine_job_to_complete(self):
    #     self.fail()

    # write_runscript_for_pdbs() takes 6 arguments. The last 3 (namely show_sequence_detail, print_networks,
    # calculate_stability) are keyword (named) arguments. All have default values assigned in the method argument so
    # that when the method is called without supplying a (keyword-named) value, the default value is applied.
    # (The default values for all 3 are False).
    def test_write_runscript_for_pdbs(self):
        # arrange
        pdb = 'RepairPDB_1.pdb'
        path_runscript = self.path_tests + self.path_rel_Inputs + self.path_rel_PDBs + '/' + pdb.split('.')[0] + \
                         '/FX_BuildModel'
        if not os.path.exists(path_runscript):
            os.makedirs(path_runscript)
        action = '<BuildModel>#,individual_list.txt'
        expected_runscript = '<TITLE>FOLDX_runscript;\n' + '<JOBSTART>#;\n' + '<PDBS>' + pdb + ';\n' + \
                             '<BATCH>#;\n' + '<COMMANDS>FOLDX_commandfile;\n' + action + ';\n' + '<END>#;\n' + \
                             '<OPTIONS>FOLDX_optionfile;\n' + '<Temperature>298;\n' + '<IonStrength>0.05;\n' + \
                             '<ph>7;\n' + '<moveNeighbours>true;\n' + '<VdWDesign>2;\n' + '<numberOfRuns>3;\n' + \
                             '<OutPDB>#;\n' + '<END>#;\n' + '<JOBEND>#;\n' + '<ENDFILE>#;\n'
        single_space = ' '
        not_expected_runscript = single_space + expected_runscript
        # action
        actual_runscript = GUM.write_runscript_for_pdbs(path_runscript, pdb, action)
        # assert
        self.assertEqual(actual_runscript, expected_runscript)
        self.assertNotEqual(actual_runscript, not_expected_runscript)

    # @patch.object(GUM, '_make_subfoldername')
    def test_copy_and_move_files(self):
        # arrange
        path_src_dir = '/Users/u0120577/ROB_HOMOLOGY_HUMAN/RepairPDBstest'
        path_dst_dir = self.path_tests_Inputs_PDBs
        starting_num = 1
        # action
        GUM.copy_and_move_pdb_files(path_src_dir, path_dst_dir, starting_num, total_num_to_copy=9)
        # assert
        # expected_call_count_1 = 1
        # self.assertEqual(expected_call_count_1, mock_make_subfoldername.call_count)
