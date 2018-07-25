from unittest import TestCase
from Bio.Blast import NCBIWWW
from src.Biopython import Biopy
from unittest.mock import patch
from src.GeneralUtilityMethods import GUM
import os
import filecmp


class TestBiopython(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.maxDiff = None
        cls.PATH_FASTA_INPUT_FILES = '/Users/u0120577/ROB_HOMOLOGY_HUMAN/Fasta'
        cls.SWSPRT_DB = 'swissprot'
        cls.SWSPRT_PROTS_NUM = 20341
        cls.ZERO_GAP = 0
        cls.QSTRT_1 = 1
        cls.PATH_OUTPUT_ROOT = "/Users/u0120577/SNPEFFECT/Outputs"
        cls.PATH_TESTS_OUTPUTS = '/Users/u0120577/PycharmProjects/MutateCompute/tests/Outputs'

        cls.FASTA_1_A_FILE = '1_A.fasta'
        cls.PATH_FASTA_1_A = cls.PATH_FASTA_INPUT_FILES + '/' + cls.FASTA_1_A_FILE
        cls.FASTA_1_A_SEQ = 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGN' \
                             'FVRVIQTFNRTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRS' \
                             'LTKRNAVRTDQHNSKWLSEPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARL' \
                             'VCSVTDEDGPETHFDELEDVFLLETDNPRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRP' \
                             'GTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSIYPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVL' \
                             'PTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVSQVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRR' \
                             'SRRQDVRHGNPLTQCR'
        cls.FASTA_1_A_NAME = cls.FASTA_1_A_FILE.split('.')[0]
        cls.FASTA_1_A_XML = cls.FASTA_1_A_NAME + '.xml'
        cls.FASTA_1_A_STR = ">" + cls.FASTA_1_A_NAME + "\n" + cls.FASTA_1_A_SEQ

        # Create any test files needed for unit tests and write to tests/Files folder
        cls.PATH_TESTS_FILES = '/Users/u0120577/PycharmProjects/MutateCompute/tests/Files'
        cls.BLAST_RUN_DIR = 'blast_run'
        path_tests_Files_blastrun_1_A_xml = cls.PATH_TESTS_FILES + '/' + cls.BLAST_RUN_DIR + '/' + cls.FASTA_1_A_XML
        if not os.path.exists(path_tests_Files_blastrun_1_A_xml):
            result_handle_1_A = NCBIWWW.qblast(program=Biopy.BlastParam.BLST_P.value,
                                               database=Biopy.BlastParam.SWSPRT.value,
                                               sequence=cls.FASTA_1_A_STR,
                                               entrez_query=Biopy.BlastParam.HOMSAP_ORG.value,
                                               alignments=Biopy.BlastParam.MAX_ALIGN_20.value,
                                               hitlist_size=Biopy.BlastParam.MAX_HIT_20.value)
            TestBiopython._write_blast_run_to_tests_dir(cls.PATH_TESTS_FILES, result_handle_1_A, cls.FASTA_1_A_XML)

    @classmethod
    def tearDownClass(cls):
        cls.PATH_FASTA_INPUT_FILES = None
        cls.SWSPRT_DB = None
        cls.SWSPRT_PROTS_NUM = None
        cls.ZERO_GAP = None
        cls.QSTRT_1 = None
        cls.PATH_OUTPUT_ROOT = None
        cls.FASTA_1_A_FILE = None
        cls.PATH_FASTA_1_A = None
        cls.FASTA_1_A_SEQ = None
        cls.FASTA_1_A_NAME = None
        cls.FASTA_1_A_STR = None
        cls.blast_res_1_A = None

    # This test is not ideal. It relies on runnning the actual qblast and writing the result to an xml
    # file, if only once. It then runs the qblast method in src.Biopython.Biopy but also has to write the
    # result to an xml in the test method. Hence it does not exclusively test the _run_blastp method. It
    # can fail if writing the result or reading the pre-written result fails. (I cannot find a more direct
    # mechanism for file comparison when the file is still in the io file format.)
    # NOTE: using filecmp.cmp() is useless because each blastp run has a unique query id which is written in the output
    def test_run_blastp_1_A(self):
        # arrange
        filecmp.clear_cache()
        path_actual = self.PATH_TESTS_OUTPUTS + '/' + self.BLAST_RUN_DIR + '/' + self.FASTA_1_A_XML
        path_expected = self.PATH_TESTS_FILES + '/' + self.BLAST_RUN_DIR + '/' + self.FASTA_1_A_XML
        # action
        result_handle_1_A = Biopy._run_blastp(self.FASTA_1_A_STR)
        TestBiopython._write_blast_run_to_tests_dir(self.PATH_TESTS_OUTPUTS, result_handle_1_A, self.FASTA_1_A_XML)
        # assert
        with open(path_actual) as actual, open(path_expected) as expected:
            for actual_line, exp_line in zip(actual, expected):
                if 'Iteration_query-ID' in exp_line or 'BlastOutput_query-ID' in exp_line:
                    self.assertNotEqual(exp_line, actual_line)
                else:
                    self.assertEqual(exp_line, actual_line)




    # path_tests        tests/Outputs where actual test results are stored
    #                   tests/Files where reference results are stored
    # result_handle     the direct return object/handle of Biopython's qblast method
    @staticmethod
    def _write_blast_run_to_tests_dir(path_tests, result_handle, file):
        path_blast_run = GUM.create_dir_tree(path_tests, TestBiopython.BLAST_RUN_DIR)
        path_blastrun_xml = path_blast_run + '/' + file
        with open(path_blastrun_xml, 'w') as resx:
            resx.write(result_handle.read())
            result_handle.close()

    @patch.object(Biopy, "_write_qblast_xml_result")
    @patch.object(NCBIWWW, "qblast")
    @patch.object(Biopy, "_read_input_fasta_seq")
    def test_find_identical_blastp_hit_swissprot_for_fasta_input_1_A(self, mock__read_input_fasta_seq, mock_qblast,
                                                                 mock__write_qblast_xml_result):
        # arrange
        len_1_A = len(self.FASTA_1_A_SEQ)
        query_len = align_len = idents = q_end = len_1_A
        acc = 'Q99985'
        hit_len = 751
        hit_start_pos = 28
        hit_end_pos = hit_start_pos + query_len - 1
        hit_def = 'RecName: Full=Semaphorin-3C; AltName: Full=Semaphorin-E; Short=Sema E; Flags: Precursor'
        expected_align_dict_1_A = {'accession_num': acc, 'length': hit_len, 'hit_def': hit_def,
                                       'hsp_dict': {'align_length': align_len, 'gaps': self.ZERO_GAP,
                                                    'identities': idents, 'query_end': q_end,
                                                    'query_start': self.QSTRT_1, 'sbjct_end': hit_end_pos,
                                                    'sbjct_start': hit_start_pos}}
        expected_qblast_dict_1_A = {'database': self.SWSPRT_DB, 'database_seqs_num': self.SWSPRT_PROTS_NUM,
                                        'query_length': query_len, 'query_seq_id': self.FASTA_1_A_NAME,
                                        'alignment_dict': expected_align_dict_1_A}
        mock__read_input_fasta_seq.return_value = self.FASTA_1_A_STR
        with open(self.PATH_OUTPUT_ROOT + "/test_blastp_1_A.xml") as test_1_A_xml:
            mock_qblast.return_value = test_1_A_xml
        mock__write_qblast_xml_result.return_value = self.PATH_OUTPUT_ROOT + "/blastp_swissprot/1_A.xml"
        # action
        result = Biopy.find_identical_blastp_hit_swissprot_for_fasta_input(self.PATH_FASTA_1_A, self.PATH_OUTPUT_ROOT)
        # assert
        self.assertEqual(expected_qblast_dict_1_A, result)

    # @patch.object(Biopy, "_write_qblast_xml_result")
    # @patch.object(NCBIWWW, "qblast")
    # @patch.object(Biopy, "_read_input_fasta_seq")
    # def test_find_identical_blastp_hit_swissprot_for_fasta_input_1_B(self, mock__read_input_fasta_seq, mock_qblast,
    #                                                                  mock__write_qblast_xml_result):
    #     # arrange
    #     print('todo')
    #     # action
    #     # assert
    #
    # @patch.object(Biopy, "_write_qblast_xml_result")
    # @patch.object(NCBIWWW, "qblast")
    # @patch.object(Biopy, "_read_input_fasta_seq")
    # def test_find_identical_blastp_hit_swissprot_for_fasta_input_2_A(self, mock__read_input_fasta_seq, mock_qblast,
    #                                                                  mock__write_qblast_xml_result):
    #     # arrange
    #     print('todo')
    #     # action
    #     # assert
    #
    # @patch.object(Biopy, "_write_qblast_xml_result")
    # @patch.object(NCBIWWW, "qblast")
    # @patch.object(Biopy, "_read_input_fasta_seq")
    # def test_find_identical_blastp_hit_swissprot_for_fasta_input_3_A(self, mock__read_input_fasta_seq, mock_qblast,
    #                                                                  mock__write_qblast_xml_result):
    #     # arrange
    #     print('todo')
    #     # action
    #     # assert
    #
    # @patch.object(Biopy, "_write_qblast_xml_result")
    # @patch.object(NCBIWWW, "qblast")
    # @patch.object(Biopy, "_read_input_fasta_seq")
    # def test_find_identical_blastp_hit_swissprot_for_fasta_input_3_B(self, mock__read_input_fasta_seq, mock_qblast,
    #                                                                  mock__write_qblast_xml_result):
    #     # arrange
    #     print('todo')
    #     # action
    #     # assert
