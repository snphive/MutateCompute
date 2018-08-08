from unittest import TestCase
from src.MutateFasta import MutateFasta
from tests.TestPathsAndLists import TPL


class TestMutateFasta(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.mutateFasta = MutateFasta()

    @classmethod
    def tearDownClass(cls):
        cls.mutateFasta = None

# MAKES SENSE TO WRITE TESTS FOR THE PRIVATE METHODS FIRST.
    # def test_mutate_every_residue_in_fasta_list(self):
    #     # arrange
    #     path_input = TPL.MC_TESTS_INPUT
    #     with open('/Users/u0120577/FASTA_repo100') as fastafiles
    #         fasta_list = fastafiles.readlines()
    #     mutant_aa_list = TPL.LIST_ALL_20_AA
    #
    #     # action
    #     MutateFasta().mutate_every_residue_in_fasta_list(path_input, fasta_list, mutant_aa_list)
    #     # assert

    def test__build_complete_paths_for_fastafiles(self):
        # arrange
        path_input = TPL.MC_TESTS_INPUT.value
        fastafile_list = ['1_A.fasta', '1_B.fasta', '2_A.fasta']
        expected_path_list = ['/Users/u0120577/PycharmProjects/MutateCompute/tests/input_data/fastas/1_A/1_A.fasta',
                              '/Users/u0120577/PycharmProjects/MutateCompute/tests/input_data/fastas/1_B/1_B.fasta',
                              '/Users/u0120577/PycharmProjects/MutateCompute/tests/input_data/fastas/2_A/2_A.fasta']
        mutateFasta = MutateFasta()
        # action
        path_list = mutateFasta._build_complete_paths_for_fastafiles(path_input, fastafile_list)
        # assert
        self.assertListEqual(expected_path_list, path_list)

    # def test__build_list_of_title_sequence_dict_from_fastafiles(self):
    #     self.fail()
    #
    # def test__mutate_sequences_in_list_of_dict(self):
    #     self.fail()
    #
    # def test__mutate_fasta(self):
    #     self.fail()
