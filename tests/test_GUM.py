from unittest import TestCase
from src.GeneralUtilityMethods import GUM
from tests.HelperMethods import HM
import os
from unittest.mock import patch
from tests.TestPaths import TPaths
import subprocess


class TestGUM(TestCase):

    # Currently the tests are copying over all configuration and input data from the main directory into the tests before
    # running the tests (i.e. here in the setUpClass method).
    # The data in those main folders will be programmatically generated but is currently manually transferred.
    @classmethod
    def setUpClass(cls):
        if not os.path.exists(TPaths.MC_TESTS_CONFIG.value):
            GUM.linux_copy(path_src=TPaths.CONFIG_FOR_READ_ONLY.value, path_dst=TPaths.MC_TESTS_CONFIG.value,
                           do_recursively=True)

        if not os.path.exists(TPaths.MC_TESTS_INPUT.value):
            GUM.linux_copy(path_src=TPaths.INPUT_FOR_READ_ONLY.value, path_dst=TPaths.MC_TESTS_INPUT.value,
                           do_recursively=True)

    # Write_runscript_for_pdbs() takes 6 arguments. The last 3 (namely show_sequence_detail, print_networks,
    # calculate_stability) are keyword (named) arguments. All have default values assigned in the method argument so
    # that when the method is called without supplying a (keyword-named) value, the default value is applied.
    # (The default values for all 3 are False).
    #
    # 30.07.18 Have redesigned the directory structure such that runscripts will go in configuration/foldx/and maybe
    # another
    # level such as analyse_complex or build_model or stability etc.
    def test_write_runscript_for_pdbs(self):
        # arrange
        pdb = 'RepairPDB_1.pdb'
        path_runscript = TPaths.MC_TESTS_CONFIG_FX.value
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

    # CHECK WHAT HAPPENS IF THERE ARE FEWER FILES THAN SPECIFIED TO MOVE BY TOTAL_NUM_TO_COPY
    # @patch.object(GUM, '_make_subfoldername')
    def test_copy_and_move_pdb_files(self):
        # arrange
        path_src_dir = TPaths.MC_TESTS_INPUT.value
        path_dst_dir = TPaths.MC_TESTS_OUTPUT.value
        starting_num = 1
        # action
        GUM.copy_and_move_pdb_files(path_src_dir, path_dst_dir, starting_num, total_num_to_copy=9)
        # assert
        # expected_call_count_1 = 1
        # self.assertEqual(expected_call_count_1, mock_make_subfoldername.call_count)

    def test_copy_input_files_from_repo_to_input(self):
        # arrange
        path_src_dir = '/Users/u0120577/PDB_repo10'
        path_dst_dir = TPaths.MC_TESTS_INPUT.value
        files = os.listdir(path_src_dir)
        subprocess.call('cp -r /Users/u0120577/PDB_repo10 ' + TPaths.MC_TESTS_REFFILES.value, shell=True)
        # action
        GUM.copy_files_from_repo_to_input_filedir(path_src_dir, path_dst_dir)
        # assert
        for file in files:
            self.assertTrue(os.path.exists(TPaths.MC_TESTS_INPUT.value + '/' + file.split('.')[0] + '/' + file))
