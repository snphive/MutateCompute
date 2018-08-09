import os
import glob
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
        titleSeqDict = {'WT_SCI': 'SCI'}
        wt_title = 'WT_SCI'
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
        titleTitleSeqDictDict = {'WT1_SC': {'WT1_SC': 'SC'}, 'WT2_GE': {'WT2_GE': 'GE'}}
        mutant_aa_list = ['A', 'C']
        expected_titleTitleSeqDictDict_w_muts = {'WT1_SC': {'WT1_SC': 'SC', 'S1A': 'AC', 'S1C': 'CC', 'C2A': 'SA'},
                                                 'WT2_GE': {'WT2_GE': 'GE', 'G1A': 'AE', 'G1C': 'CE', 'E2A': 'GA',
                                                         'E2C': 'GC'}}
        # action
        titleTitleSeqDictDict_w_muts = self.mutateFasta._populate_titleTitleSeqDictDict_with_mutants(
            titleTitleSeqDictDict, mutant_aa_list)
        # assert
        self.assertDictEqual(expected_titleTitleSeqDictDict_w_muts, titleTitleSeqDictDict_w_muts)

    def test_mutate_every_residue_in_fasta_list(self):
        # arrange
        path_input = TPLS.MC_TESTS_INPUT.value
        fastafile_list = ['WT_SCI.fasta']
        mutant_aa_list = ['A']
        expected_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        # action
        mutants = self.mutateFasta.mutate_every_residue_in_fasta_list(path_input, fastafile_list, mutant_aa_list)
        # assert
        self.assertDictEqual(expected_mutants, mutants)

    # Write one fastafile for all mutants. All other flags are set to false.
    def test_write_mutants_to_file_1(self):
        # arrange
        title_titleSeq_w_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        path_input = TPLS.MC_TESTS_INPUT.value
        make_fastafile_all_mutants = True
        make_fastafile_per_mutant = False
        path_output = ''
        expected_all_mutants_fastafile = '>S1A\nACI\n>C2A\nSAI\n>I3A\nSCA\n'
        # action
        self.mutateFasta._write_mutants_to_file(title_titleSeq_w_mutants, path_input, make_fastafile_all_mutants,
                               make_fastafile_per_mutant, path_output, make_csv_file=False, make_txt_file=False)
        # assert
        self.assertTrue(os.path.exists(TPLS.MC_TESTS_INPUT.value + '/fastas/WT_SCI/mutants'),
            msg='~/PycharmProjects/MutateCompute/tests/input_data/fastas/WT_SCI/mutants folder was not found')
        self.assertTrue(os.path.exists(TPLS.MC_TESTS_INPUT.value + '/fastas/WT_SCI/mutants/WT_SCI_mutants.fasta'),
            msg='~/PycharmProjects/MutateCompute/tests/input_data/fastas/WT_SCI/mutants/WT_SCI.fasta file not found')
        with open(TPLS.MC_TESTS_INPUT.value + '/fastas/WT_SCI/mutants/WT_SCI_mutants.fasta') as all_mutants_fastafile:
            all_mutants_fastafile = all_mutants_fastafile.read()
        self.assertEqual(expected_all_mutants_fastafile, all_mutants_fastafile)

    # Write one fastafile per mutant. All other flags are false.
    def test_write_mutants_to_file_2(self):
        # arrange
        title_titleSeq_w_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        path_input = TPLS.MC_TESTS_INPUT.value
        make_fastafile_all_mutants = False
        make_fastafile_per_mutant = True
        path_output = ''
        expected_mutant_fastafile_list = ['S1A.fasta', 'C2A.fasta', 'I3A.fasta']
        # action
        self.mutateFasta._write_mutants_to_file(title_titleSeq_w_mutants, path_input, make_fastafile_all_mutants,
                               make_fastafile_per_mutant, path_output, make_csv_file=False, make_txt_file=False)
        mutant_fastafile_list = glob.glob(TPLS.MC_TESTS_INPUT.value + '/fastas/WT_SCI/mutants/*.fasta')
        mutant_fastafile_list = [mutant_fastafile.split('/')[-1] for mutant_fastafile in mutant_fastafile_list]

        # assert
        self.assertTrue(os.path.exists(TPLS.MC_TESTS_INPUT.value + '/fastas/WT_SCI/mutants'),
            msg='~/PycharmProjects/MutateCompute/tests/input_data/fastas/WT_SCI/mutants folder not found')
        for exp_mutant_fastafile in expected_mutant_fastafile_list:
            self.assertTrue(exp_mutant_fastafile in mutant_fastafile_list)

    # Write one txt file (for wild-type and all mutants). All other flags are false.
    def test_write_mutants_to_file_3(self):
        # arrange
        title_titleSeq_w_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        path_input = TPLS.MC_TESTS_INPUT.value
        make_fastafile_all_mutants = False
        make_fastafile_per_mutant = False
        make_csv_file = False
        make_txt_file = True
        path_output = TPLS.MC_TESTS_OUTPUT.value
        expected_txt_file_str = 'WT_SCI:SCI\nS1A:ACI\nC2A:SAI\nI3A:SCA\n'
        # action
        self.mutateFasta._write_mutants_to_file(title_titleSeq_w_mutants, path_input, make_fastafile_all_mutants,
                    make_fastafile_per_mutant, path_output, make_csv_file=make_csv_file, make_txt_file=make_txt_file)
        # assert
        self.assertTrue(os.path.exists(TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences'),
            msg=TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences folder not found')
        self.assertTrue(os.path.exists(TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences/WT_SCI_mutants.txt'),
            msg=TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences/WT_SCI_mutants.txt file not found')
        with open(TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences/WT_SCI_mutants.txt') as txt_opened:
            txt_file_str = txt_opened.read()
        self.assertEqual(expected_txt_file_str, txt_file_str)

    # Write one csv file (for wild-type and all mutants). All other flags are false.
    def test_write_mutants_to_file_4(self):
        # arrange
        title_titleSeq_w_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        path_input = TPLS.MC_TESTS_INPUT.value
        make_fastafile_all_mutants = False
        make_fastafile_per_mutant = False
        make_csv_file = True
        make_txt_file = False
        path_output = TPLS.MC_TESTS_OUTPUT.value
        expected_csv_file_str = 'WT_SCI:SCI,S1A:ACI,C2A:SAI,I3A:SCA,'
        # action
        self.mutateFasta._write_mutants_to_file(title_titleSeq_w_mutants, path_input, make_fastafile_all_mutants,
                    make_fastafile_per_mutant, path_output, make_csv_file=make_csv_file, make_txt_file=make_txt_file)
        # assert
        self.assertTrue(os.path.exists(TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences'),
            msg=TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences folder not found')
        self.assertTrue(os.path.exists(TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences/WT_SCI_mutants.csv'),
            msg=TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences/WT_SCI_mutants.csv file not found')
        with open(TPLS.MC_TESTS_OUTPUT.value + '/WT_SCI/sequences/WT_SCI_mutants.csv') as csv_opened:
            csv_file_str = csv_opened.read()
        self.assertEqual(expected_csv_file_str, csv_file_str)
