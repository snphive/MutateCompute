from unittest import TestCase
from src.Agadir import Agadir
from tests.TestPathsAndListsSeqs import TPLS
from src.GeneralUtilityMethods import GUM
import os
import glob


class TestAgadir(TestCase):

    def setUp(self):
        self.agadir = Agadir()

    def tearDown(self):
        self.agadir = None

    def test_copy_files_from_output_fastafilemutants_dir_to_input(self):
        # arrange
        # path_repo_fastafiles = TPLS.REPO_PDB_FASTA.value + '/**/*.fasta'
        path_repo_fastafilelist = glob.glob(TPLS.REPO_PDB_FASTA.value + '/fastas_10/*.fasta')

        path_dst_arrange = os.path.join(TPLS.MC_TESTS_OUTPUT.value, TPLS.DIR_FASTAS.value)
        GUM.linux_copy_files_into_own_subdir(path_src_filelist=path_repo_fastafilelist, path_dst=path_dst_arrange)
        # path_dst_arrange = os.path.join(path_dst_arrange, TPLS.DIR_MUTANTS.value)
        expected_filelist = ['1_A.fasta', '1_B.fasta', '2_A.fasta', '3_A.fasta', '3_B.fasta',
                             '4_A.fasta', '5_A.fasta', '5_B.fasta', '5_C.fasta', '5_D.fasta']
        path_dst_input = os.path.join(TPLS.MC_TESTS_INPUT.value, TPLS.DIR_FASTAS.value)
        # act
        filelist_rtn = self.agadir.copy_files_from_output_to_input_subdirs(path_output=path_dst_arrange,
                                                                           path_input=path_dst_input)
        actual_filelist = os.listdir(path_dst_input)
        # assert
        self.assertListEqual(expected_filelist, actual_filelist)
        self.assertListEqual(expected_filelist, filelist_rtn)

    def test_compute(self):
        # arrange
        # act
        # assert
        self.assertTrue()
        self.assertEqual()

