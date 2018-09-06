from unittest import TestCase
from src.Agadir import Agadir
from src.Agadir import AgadCndtns
from tests.TestPathsAndListsSeqs import TPLS
from src.GeneralUtilityMethods import GUM
from src.Paths import Paths
import os
import glob


class TestAgadir(TestCase):

    def setUp(self):
        self.agadir = Agadir(AgadCndtns.INCELL_MAML.value)

    def tearDown(self):
        self.agadir = None

    def test_copy_files_output_to_input_fastasubdirs(self):
        # arrange
        path_repo_fastafilelist = glob.glob(TPLS.REPO_PDB_FASTA.value + '/fastas_10/*.fasta')
        path_output_fastas = os.path.join(TPLS.MC_TESTS_OUTPUT.value, TPLS.DIR_FASTAS.value)
        GUM.linux_copy_files(path_repo_fastafilelist, path_output_fastas, into_own_subdirs=True)
        expected_filelist = ['1_A.fasta', '1_B.fasta', '2_A.fasta', '3_A.fasta', '3_B.fasta', '4_A.fasta', '5_A.fasta',
                             '5_B.fasta', '5_C.fasta', '5_D.fasta']
        path_input_fastas = os.path.join(TPLS.MC_TESTS_INPUT.value, TPLS.DIR_FASTAS.value)
        expected_path_filelist = [path_input_fastas + '/' + file for file in expected_filelist]
        GUM.linux_remove_files_in_dir(path_input_fastas)
        # act
        self.agadir.copy_files_output_to_input_fastasubdirs(path_output_fastas, path_input_fastas, filename_list=None,
                                                            into_own_subdirs=True)
        actual_path_filelist = glob.glob(path_input_fastas + '/*.fasta')
        # assert
        self.maxDiff = None
        self.assertListEqual(expected_path_filelist, actual_path_filelist)

    def test_compute(self):
        # arrange

        # act
        self.agadir.compute(path_input_fastas=TPLS.MC_TESTS_INPUT_FASTAS.value, path_output=TPLS.MC_TESTS_OUTPUT.value)
        # assert

    def test_write_1_fastafile_per_fasta_from_multifastafile(self):
        # arrange
        Paths.set_up_paths(use_cluster=False)
        path_fastafile = '/Users/u0120577/PycharmProjects/output_data/mutants_fastas/1...1000/1_A/mutants/' \
                         '1_A_mutants.fasta'
        # act
        self.agadir.write_1_fastafile_per_fasta_from_multifastafile(path_fastafile)
