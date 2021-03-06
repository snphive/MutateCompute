import os
from unittest import TestCase
from src.MutateFasta import MutateFasta
from tests.TestPathsAndListsSeqs import TPLS
from src.tools.GeneralUtilityMethods import GUM


class TestMutateFasta(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.mutateFasta = MutateFasta()

    @classmethod
    def tearDownClass(cls):
        cls.mutateFasta = None

    def test_convert_titleSeqDict_to_titleTitleSeqDictDict(self):
        # arrange
        title_sequence_dict = {'1_A': TPLS.FASTA_SEQ_1_A.value}
        expected_titleTitleSeqDictDict = {'1_A': {'1_A': TPLS.FASTA_SEQ_1_A.value}}
        # act
        titleTitleSeqDictDict = MutateFasta.convert_titleSeqDict_to_titleTitleSeqDictDict(title_sequence_dict)
        # assert
        self.assertDictEqual(expected_titleTitleSeqDictDict, titleTitleSeqDictDict)

    def test_make_titleSeqDict_from_fastafile(self):
        # arrange
        path_fastafile_list = '/Users/u0120577/PycharmProjects/MutateCompute/tests/input_data/fastas/1_A/1_A.fasta'
        expected_title_seq_dict = {'1_A': TPLS.FASTA_SEQ_1_A.value}
        # act
        title_seq_dict = self.mutateFasta.make_titleSeqDict_from_fastafile(path_fastafile_list)
        # assert
        self.assertDictEqual(expected_title_seq_dict, title_seq_dict)

    def test__build_complete_paths_for_fastafiles(self):
        # arrange
        path_input = TPLS.MC_TESTS_INPUT.value
        fastafile_list = ['1_A.fasta', '1_B.fasta', '2_A.fasta']
        expected_path_list = ['/Users/u0120577/PycharmProjects/MutateCompute/tests/input_data/fastas/1_A/1_A.fasta',
                              '/Users/u0120577/PycharmProjects/MutateCompute/tests/input_data/fastas/1_B/1_B.fasta',
                              '/Users/u0120577/PycharmProjects/MutateCompute/tests/input_data/fastas/2_A/2_A.fasta']
        # act
        path_list = self.mutateFasta.build_complete_paths_for_fastafiles(path_input, fastafile_list)
        # assert
        self.assertListEqual(expected_path_list, path_list)

    def test__add_mutantTitle_mutatedSeq_to_dict(self):
        # arrange
        titleSeqDict = {'WT_SCI': 'SCI'}
        wt_title = 'WT_SCI'
        wt_seq = TPLS.FASTA_TRIPEP_SEQ_WILDTYPE.value
        amino_acids = TPLS.LIST_ALL_20_AA.value
        expected_mutantTitle_SeqDict = TPLS.FASTA_TRIPEP_TITLE_TITLE_SEQ_DICT_ALL_20_INCL_WT.value
        # act
        mutantTitle_SeqDict = self.mutateFasta._add_mutantTitle_mutatedSeq_to_dict(titleSeqDict, wt_title, wt_seq,
                                                                                   amino_acids)
        # assert
        # self.maxDiff = None
        self.assertEqual(len(expected_mutantTitle_SeqDict.keys()), len(mutantTitle_SeqDict.keys()))
        self.assertDictEqual(expected_mutantTitle_SeqDict, mutantTitle_SeqDict)

    # I CAN'T WORK OUT HOW TO MAKE MOCK'S SIDE-EFFECT THING WORK!! HOW TO MAKE IT MOCK A SPECIFIC METHOD TO RETURN
    # ONE VALUE ON THE FIRST TIME IT'S CALLED AND ANOTHER VALUE THE 2ND TIME IT'S CALLED.
    def test__populate_title_titleSeq_with_mutants(self):
        # arrange
        titleTitleSeqDictDict = {'WT1_SC': {'WT1_SC': 'SC'}, 'WT2_GE': {'WT2_GE': 'GE'}}
        mutant_aa_list = ['A', 'C']
        expected_titleTitleSeqDictDict_w_muts = {'WT1_SC': {'WT1_SC': 'SC', 'S1A': 'AC', 'S1C': 'CC', 'C2A': 'SA'},
                                                 'WT2_GE': {'WT2_GE': 'GE', 'G1A': 'AE', 'G1C': 'CE', 'E2A': 'GA',
                                                         'E2C': 'GC'}}
        # act
        titleTitleSeqDictDict_w_muts = self.mutateFasta._populate_title_titleSeq_with_mutants(
            titleTitleSeqDictDict, mutant_aa_list)
        # assert
        self.assertDictEqual(expected_titleTitleSeqDictDict_w_muts, titleTitleSeqDictDict_w_muts)

    def test_mutate_every_residue(self):
        # arrange
        path_input = TPLS.MC_TESTS_INPUT.value
        fastafile_list = ['WT_SCI.fasta']
        fastafile = ['WT_SCI.fasta']
        mutant_aa_list = ['A']
        write_1_fasta_only = True
        write_fasta_per_mut = False
        path_output = TPLS.MC_TESTS_OUTPUT.value
        expected_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        # act
        mutants = self.mutateFasta.mutate_every_residue(fastafile, mutant_aa_list,
                                                        write_fasta_per_mut=write_fasta_per_mut,
                                                        path_output_3dots=path_output)
        # assert
        self.assertDictEqual(expected_mutants, mutants)

    # Write one fastafile for all mutants. All other flags are set to false.
    def test_write_mutants_1(self):
        # arrange
        path_tests_output_fastas_wtsci_mutants = os.path.join(TPLS.MC_TESTS_OUTPUT_FASTAS.value, 'WT_SCI',
                                                             TPLS.DIR_MUTANTS.value)
        GUM.linux_remove_all_files_in_dir(path_tests_output_fastas_wtsci_mutants)
        title_titleSeq_w_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        expected_all_mutants_single_fastafile = '>WT_SCI\nSCI\n>S1A\nACI\n>C2A\nSAI\n>I3A\nSCA\n'
        # act
        self.mutateFasta._write_mutants(title_titleSeq_w_mutants, write_1_fasta_only=True, write_fasta_per_mut=False,
                                        path_output_3dots=TPLS.MC_TESTS_OUTPUT.value, write_csv=False, write_txt=False)
        # assert
        self.assertTrue(os.path.exists(path_tests_output_fastas_wtsci_mutants),
            msg=path_tests_output_fastas_wtsci_mutants + ' folder was not found.')

        path_tests_output_fastas_sci_mutants_wtscifilemuts = os.path.join(path_tests_output_fastas_wtsci_mutants,
                                                                            'WT_SCI_mutants.fasta')
        self.assertTrue(os.path.exists(path_tests_output_fastas_sci_mutants_wtscifilemuts),
            msg=path_tests_output_fastas_sci_mutants_wtscifilemuts + ' file not found.')

        with open(path_tests_output_fastas_sci_mutants_wtscifilemuts) as all_mutants_fastafile:
            all_mutants_fastafile = all_mutants_fastafile.read()
        self.assertEqual(expected_all_mutants_single_fastafile, all_mutants_fastafile)

    # Write one fastafile per mutant. All other flags are false.
    def test_write_mutants_2(self):
        # arrange
        path_tests_output_fastas_wtsci = os.path.join(TPLS.MC_TESTS_OUTPUT_FASTAS.value, 'WT_SCI')
        GUM.linux_remove_all_files_in_dir(path_tests_output_fastas_wtsci)
        title_titleSeq_w_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        expected_mutant_fastafile_list = ['WT_SCI.fasta', 'S1A.fasta', 'C2A.fasta', 'I3A.fasta']
        # act
        self.mutateFasta._write_mutants(title_titleSeq_w_mutants, write_1_fasta_only=False, write_fasta_per_mut=True,
                                        path_output_3dots=TPLS.MC_TESTS_OUTPUT.value, write_csv=False, write_txt=False)
        path_tests_output_fastas_wtsci_mutants = os.path.join(TPLS.MC_TESTS_OUTPUT.value, TPLS.DIR_FASTAS.value,
                                                             'WT_SCI', TPLS.DIR_MUTANTS.value)
        # mutant_fastafile_list = glob.glob(path_tests_output_fastas_wtsci + '/**/*.fasta', recursive=True)
        # mutant_fastafile_list = [mutant_fastafile.split('/')[-1] for mutant_fastafile in mutant_fastafile_list]
        mutant_fastafile_list = GUM.make_filelist_in_current_dir_and_subdirs_recursively(path_tests_output_fastas_wtsci,
                                                                                     '.fasta')
        # assert
        self.assertTrue(path_tests_output_fastas_wtsci_mutants, msg=path_tests_output_fastas_wtsci_mutants +
                                                                    ' folder not found.')
        for exp_mutant_fastafile in expected_mutant_fastafile_list:
            self.assertTrue(exp_mutant_fastafile in mutant_fastafile_list)

    # Write one txt file (for wild-type and all mutants). All other flags are false.
    def test_write_mutants_3(self):
        # arrange
        path_tests_output_seqs_wtsci = os.path.join(TPLS.MC_TESTS_OUTPUT.value, TPLS.DIR_SEQS_TXT_CSV.value, 'WT_SCI')
        GUM.linux_remove_all_files_in_dir(path_tests_output_seqs_wtsci)
        title_titleSeq_w_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        expected_txt_file_str = 'WT_SCI:SCI\nS1A:ACI\nC2A:SAI\nI3A:SCA\n'
        # act
        self.mutateFasta._write_mutants(title_titleSeq_w_mutants, write_1_fasta_only=False, write_fasta_per_mut=False,
                                        path_output_3dots=TPLS.MC_TESTS_OUTPUT.value, write_csv=False, write_txt=True)
        # assert
        self.assertTrue(path_tests_output_seqs_wtsci, msg=path_tests_output_seqs_wtsci + ' folder not found.')
        path_tests_output_seqs_wtsci_mutstxt = os.path.join(path_tests_output_seqs_wtsci, 'WT_SCI_mutants.txt')
        self.assertTrue(os.path.exists(path_tests_output_seqs_wtsci_mutstxt),
            msg=path_tests_output_seqs_wtsci_mutstxt + ' file not found.')
        with open(path_tests_output_seqs_wtsci_mutstxt) as txt_opened:
            txt_file_str = txt_opened.read()
        self.assertEqual(expected_txt_file_str, txt_file_str)

    # Write one csv file (for wild-type and all mutants). All other flags are false.
    def test_write_mutants_4(self):
        # arrange
        path_tests_output_seqs_wtsci = os.path.join(TPLS.MC_TESTS_OUTPUT.value, TPLS.DIR_SEQS_TXT_CSV.value, 'WT_SCI')
        GUM.linux_remove_all_files_in_dir(path_tests_output_seqs_wtsci)
        title_titleSeq_w_mutants = {'WT_SCI': {'WT_SCI': 'SCI', 'S1A': 'ACI', 'C2A': 'SAI', 'I3A': 'SCA'}}
        expected_csv_file_str = 'WT_SCI:SCI,S1A:ACI,C2A:SAI,I3A:SCA,'
        # act
        self.mutateFasta._write_mutants(title_titleSeq_w_mutants, write_1_fasta_only=False, write_fasta_per_mut=False,
                                        path_output_3dots=TPLS.MC_TESTS_OUTPUT.value, write_csv=True, write_txt=False)
        # assert
        self.assertTrue(path_tests_output_seqs_wtsci, msg=path_tests_output_seqs_wtsci + ' folder not found')
        path_tests_output_seqs_wtscimutscsv = os.path.join(path_tests_output_seqs_wtsci, 'WT_SCI_mutants.csv')
        self.assertTrue(os.path.exists(path_tests_output_seqs_wtscimutscsv),
                        msg=path_tests_output_seqs_wtscimutscsv + ' file not found.')
        with open(path_tests_output_seqs_wtscimutscsv) as csv_opened:
            csv_file_str = csv_opened.read()
        self.assertEqual(expected_csv_file_str, csv_file_str)
