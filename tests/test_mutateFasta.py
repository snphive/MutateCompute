from unittest import TestCase
from src.MutateFasta import MutateFasta
from tests.TestPathsAndListsSeqs import TPLS
from unittest.mock import patch
from unittest.mock import MagicMock
from src.GeneralUtilityMethods import GUM


class TestMutateFasta(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.mutateFasta = MutateFasta()

    @classmethod
    def tearDownClass(cls):
        cls.mutateFasta = None

    def test__add_mutantTitle_mutatedSeq_to_dict(self):
        # arrange
        titleSeqDict = {'WT': 'SCI'}
        wt_title = 'WT'
        wt_seq = TPLS.FASTA_TRIPEP_SEQ_WILDTYPE.value
        mutant_aa_list = TPLS.LIST_ALL_20_AA.value
        expected_mutantTitle_SeqDict = TPLS.FASTA_TRIPEP_TITLE_TITLE_SEQ_DICT_ALL_20_INCL_WT.value
        # action
        mutantTitle_SeqDict = self.mutateFasta._add_mutantTitle_mutatedSeq_to_dict(titleSeqDict, wt_title, wt_seq,
                                                                                   mutant_aa_list)
        # assert
        # self.maxDiff = None
        self.assertEqual(len(expected_mutantTitle_SeqDict.keys()), len(mutantTitle_SeqDict.keys()))
        self.assertDictEqual(expected_mutantTitle_SeqDict, mutantTitle_SeqDict)

    # I CAN'T WORK OUT HOW TO MAKE MOCK'S SIDE-EFFECT THING WORK!! HOW TO MAKE IT MOCK A SPECIFIC METHOD TO RETURN
    # ONE VALUE ON THE FIRST TIME IT'S CALLED AND ANOTHER VALUE THE 2ND TIME IT'S CALLED.
    def test__populate_titleTitleSeqDictDict_with_mutants(self):
        # arrange
        titleTitleSeqDictDict = {'WT1': {'WT1': 'SC'}, 'WT2': {'WT2': 'GE'}}
        mutant_aa_list = ['A', 'C']
        expected_titleTitleSeqDictDict_w_muts = {'WT1': {'WT1': 'SC', 'S1A': 'AC', 'S1C': 'CC', 'C2A': 'SA'},
                                                 'WT2': {'WT2': 'GE', 'G1A': 'AE', 'G1C': 'CE', 'E2A': 'GA',
                                                         'E2C': 'GC'}}
        # action
        titleTitleSeqDictDict_w_muts = self.mutateFasta._populate_titleTitleSeqDictDict_with_mutants(
            titleTitleSeqDictDict, mutant_aa_list)
        # assert
        self.assertDictEqual(expected_titleTitleSeqDictDict_w_muts, titleTitleSeqDictDict_w_muts)

    # TODO
    def test_mutate_every_residue_in_fasta_list(self):
        # arrange
        path_input = TPLS.MC_TESTS_INPUT.value
        fastafile_list = ['WT']
        mutant_aa_list = ['A']
        # expected_mutants = # TODO
        # action
        mutants = self.mutateFasta.mutate_every_residue_in_fasta_list(path_input, fastafile_list, mutant_aa_list)
        # assert
        # self.assertDictEqual(expected_mutants, mutants)

    # TODO
    def test_write_mutants_to_file(self):
        # arrange
        mutants = {}
        # action
        self.mutateFasta._write_mutants_to_file(mutants)
        # assert
        # compare with reference files