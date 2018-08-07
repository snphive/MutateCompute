from unittest import TestCase
from src.MutateFasta import MutateFasta
from tests.TestPathsAndLists import TPL

class TestMutateFasta(TestCase):


# MAKES SENSE TO WRITE TESTS FOR THE PRIVATE METHODS FIRST.
    # def test_mutate_every_residue_in_fasta_list(self):
    #     # arrange
    #     path_input = TPL.MC_TESTS_INPUT
    #     with open('/Users/u0120577/FASTA_repo100') as fasta_files
    #         fasta_list = fasta_files.readlines()
    #     mutant_aa_list = TPL.LIST_ALL_20_AA
    #
    #     # action
    #     MutateFasta().mutate_every_residue_in_fasta_list(path_input, fasta_list, mutant_aa_list)
    #     # assert

    def test__build_path_for_fasta_file(self):
        # arrange
        path_input = TPL.MC_TESTS_INPUT
        fasta_list = ['1_A.fasta', '1_B.fasta', '2_A.fasta', '3_A.fasta']

        # action
        path_fasta_file_list = MutateFasta._build_path_for_fasta_file(path_input, fasta_list)
        # assert
        # assertEqual()

    def test__build_list_of_title_sequence_dict_from_fasta_files(self):
        self.fail()

    def test__mutate_sequences_in_list_of_dict(self):
        self.fail()

    def test__mutate_fasta(self):
        self.fail()
