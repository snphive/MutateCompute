import os
import io
from unittest import TestCase
from unittest import mock
from unittest.mock import patch
from src.Biopython import Biopy
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from tests.HelperMethods import HM
from src.IdentifyProtein import IdProt
from tests.TestPathsAndListsSeqs import TPLS


class TestIdProt(TestCase):

    # constants related to blastp parameters - for all mutants
    SWSPRT_DB = 'swissprot'
    SWSPRT_PROTS_NUM = 20341
    ZERO_GAP = 0
    QSTRT_1 = 1
    # constants related to output - for all mutants
    DIR_BLASTP = 'blastp'
    REL_BLASTP = '/' + DIR_BLASTP
    PATH_TESTS_REFFILES_BLASTP = TPLS.MC_TESTS_REFFILES.value + REL_BLASTP

    @classmethod
    def setUpClass(cls):

        #     mutant 1_A     ##################################################################################

        # constants related to output - mutant 1_A
        cls.XML_1_A_BLASTP_OUTPUT_FILE = '1_A.xml'

        # constants related to input - mutant 1_A
        cls.FASTAFILE_1_A = '1_A.fasta'
        cls.DIR_1_A = '1_A'
        cls.PATH_FASTA_1_A = os.path.join(TPLS.MC_TESTS_INPUT_FASTAS.value, cls.DIR_1_A, cls.FASTAFILE_1_A)
        cls.FASTA_SEQ_1_A = TPLS.FASTA_SEQ_1_A.value
        cls.NAME_1_A = cls.FASTAFILE_1_A.split('.')[0]
        cls.XML_FILE_1_A = cls.NAME_1_A + '.xml'
        cls.FASTA_STR_1_A = ">" + cls.NAME_1_A + "\n" + cls.FASTA_SEQ_1_A

        #     mutant 1_B     ##################################################################################

        cls.FASTAFILE_1_B = '1_B.fasta'
        cls.DIR_1_B = '1_B'
        cls.PATH_FASTA_1_B = os.path.join(TPLS.MC_TESTS_INPUT_FASTAS.value, cls.DIR_1_B, cls.FASTAFILE_1_B)
        cls.FASTA_SEQ_1_B = TPLS.FASTA_SEQ_1_B.value

        cls.NAME_1_B = cls.FASTAFILE_1_B.split('.')[0]
        cls.XML_FILE_1_B = cls.NAME_1_B + '.xml'
        cls.FASTA_STR_1_B = ">" + cls.NAME_1_B + "\n" + cls.FASTA_SEQ_1_B

        #     mutant 2_A     ##################################################################################

        cls.FASTAFILE_2_A = '2_A.fasta'
        cls.DIR_2_A = '2_A'
        cls.PATH_FASTA_2_A = os.path.join(TPLS.MC_TESTS_INPUT_FASTAS.value, cls.DIR_2_A, cls.FASTAFILE_2_A)
        cls.FASTA_SEQ_2_A = TPLS.FASTA_SEQ_2_A.value
        cls.NAME_2_A = cls.FASTAFILE_2_A.split('.')[0]
        cls.XML_FILE_2_A = cls.NAME_2_A + '.xml'
        cls.FASTA_STR_2_A = ">" + cls.NAME_2_A + "\n" + cls.FASTA_SEQ_2_A

        #     mutant 3_A     ##################################################################################

        cls.FASTAFILE_3_A = '3_A.fasta'
        cls.DIR_3_A = '3_A'
        cls.PATH_FASTA_3_A = os.path.join(TPLS.MC_TESTS_INPUT_FASTAS.value, cls.DIR_3_A, cls.FASTAFILE_3_A)
        cls.FASTA_SEQ_3_A = TPLS.FASTA_SEQ_3_A.value
        cls.NAME_3_A = cls.FASTAFILE_3_A.split('.')[0]
        cls.XML_FILE_3_A = cls.NAME_3_A + '.xml'
        cls.FASTA_STR_3_A = ">" + cls.NAME_3_A + "\n" + cls.FASTA_SEQ_3_A

        #     mutant 3_B     ##################################################################################

        cls.FASTAFILE_3_B = '3_B.fasta'
        cls.DIR_3_B = '3_B'
        cls.PATH_FASTA_3_B = os.path.join(TPLS.MC_TESTS_INPUT_FASTAS.value, cls.DIR_3_B, cls.FASTAFILE_3_B)
        cls.FASTA_SEQ_3_B = TPLS.FASTA_SEQ_3_B.value
        cls.NAME_3_B = cls.FASTAFILE_3_B.split('.')[0]
        cls.XML_FILE_3_B = cls.NAME_3_B + '.xml'
        cls.FASTA_STR_3_B = ">" + cls.NAME_3_B + "\n" + cls.FASTA_SEQ_3_B

        #     mutant 10_B     ##################################################################################

        cls.FASTAFILE_10_B = '10_B.fasta'
        cls.DIR_10_B = '10_B'
        cls.PATH_FASTA_10_B = os.path.join(TPLS.MC_TESTS_INPUT_FASTAS.value, cls.DIR_10_B, cls.FASTAFILE_10_B)
        cls.FASTA_SEQ_10_B = TPLS.FASTA_SEQ_10_B.value
        cls.NAME_10_B = cls.FASTAFILE_10_B.split('.')[0]
        cls.XML_FILE_10_B = cls.NAME_10_B + '.xml'
        cls.FASTA_STR_10_B = ">" + cls.NAME_10_B + "\n" + cls.FASTA_SEQ_10_B
        ###############################################################################################################

        sequence_dict = {cls.XML_FILE_1_A: cls.FASTA_STR_1_A, cls.XML_FILE_1_B: cls.FASTA_STR_1_B,
                         cls.XML_FILE_2_A: cls.FASTA_STR_2_A, cls.XML_FILE_3_A: cls.FASTA_STR_3_A,
                         cls.XML_FILE_3_B: cls.FASTA_STR_3_B, cls.XML_FILE_10_B: cls.FASTA_STR_10_B}
        TestIdProt._build_reference_blastp_output_xml_files(sequence_dict)

        cls.PATH_TESTS_REFFILES_BLASTP_1_A_XML = os.path.join(TestIdProt.PATH_TESTS_REFFILES_BLASTP, cls.XML_FILE_1_A)

    # This is a helper method, called from setUpClass()
    @staticmethod
    def _build_reference_blastp_output_xml_files(sequence_dict):
        for blast_out_xml, fasta_input in sequence_dict.items():
            if not os.path.exists(os.path.join(TestIdProt.PATH_TESTS_REFFILES_BLASTP, blast_out_xml)):
                result_handle = NCBIWWW.qblast(program=Biopy.BlastParam.BLST_P.value,
                                               database=Biopy.BlastParam.SWSPRT.value,
                                               sequence=fasta_input,
                                               entrez_query=Biopy.BlastParam.HOMSAP_ORG.value,
                                               alignments=Biopy.BlastParam.MAX_ALIGN_20.value,
                                               hitlist_size=Biopy.BlastParam.MAX_HIT_20.value)
                HM.write_blastp_to_tests_dir(TPLS.MC_TESTS_REFFILES.value, TestIdProt.DIR_BLASTP,
                                             result_handle, blast_out_xml)

    def test_start(self):
        # arrange
        path_fastafile = os.path.join(TPLS.MC_TESTS_INPUT_FASTAS.value, '1_A', '1_A.fasta')
        write_blastp_json = False
        build_idmap_csv = False
        path_output = TPLS.MC_TESTS_OUTPUT.value
        # action
        IdProt.start(path_fastafile, write_blastp_json, build_idmap_csv, path_output)
        # assert

    # I wanted to mock open() read(). Not quite managed it yet, but even so it prevents the method throwing an
    # exception if it doesn't find a file at the location it expects. As a result the method can continue and it
    # uses the subsequent qblast mock method anyway so it doesn't matter about mock_open as long as it doesn't
    # stop the method execution. It is disabled atm, so tester must be sure there is a file at the expected location
    # with the expected name.
    # mock_open(read_data="ABC")
    @patch.object(IdProt, "_write_blastp_xml_result")
    @patch.object(NCBIWWW, "qblast")
    # @patch('builtins.open', create=True) mock_open,
    def test_start_1_A(self,  mock_qblast, mock__write_blastp_xml_result):
        # arrange
        len_1_A = len(self.FASTA_SEQ_1_A)
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
                                        'query_length': query_len, 'query_seq_id': self.NAME_1_A,
                                        'alignment_dict': expected_align_dict_1_A}
        with open(self.PATH_TESTS_REFFILES_BLASTP_1_A_XML) as test_1_A_xml:
            mock_qblast.return_value = test_1_A_xml
        mock__write_blastp_xml_result.return_value = self.PATH_TESTS_REFFILES_BLASTP_1_A_XML
        path_output = os.path.join(TPLS.MC_TESTS_OUTPUT_BLASTP.value, self.DIR_1_A)
        # action
        result_dict = IdProt.start(path_fastafile=self.PATH_FASTA_1_A, write_blastp_json=False,
                                   build_idmap_csv=False, path_output=path_output)
        # assert
        self.assertEqual(expected_qblast_dict_1_A, result_dict)

    def test__write_blastp_xml_result(self):
        # arrange
        path_output = TPLS.MC_TESTS_OUTPUT.value
        fastafile_name = self.DIR_1_A
        with open(os.path.join(TPLS.MC_TESTS_REFFILES.value, TPLS.DIR_BLASTP.value, self.XML_FILE_1_A)) as xml_1_A_open:
            blastp_result = xml_1_A_open
            expected_path_output_blastp_xmlfile = os.path.join(TPLS.MC_TESTS_OUTPUT_BLASTP.value, self.XML_FILE_1_A)
            # action
            path_output_blastp_xmlfile = IdProt._write_blastp_xml_result(path_output, fastafile_name, blastp_result)
        # assert
        self.assertTrue(os.path.exists(expected_path_output_blastp_xmlfile))
        self.assertEqual(expected_path_output_blastp_xmlfile, path_output_blastp_xmlfile)

    def test__write_dict_to_json_file(self):
        # arrange
        path_output = TPLS.MC_TESTS_OUTPUT.value
        filename = ''
        blastp_result_dict = {}
        expected_path_res_file = ''
        # action
        path_result_file = IdProt._write_dict_to_json_file(path_output, filename, blastp_result_dict)
        # assert
        self.assertEqual()
        self.assertTrue(os.path.exists(path_result_file))

    def test__build_idmap_for_mysqldb(self):
        # arrange
        blastp_result_dict = {}
        filename = ''
        # action
        IdProt._build_idmap_for_mysqldb(blastp_result_dict, filename)
        # assert
        self.assertEqual()

    def test__write_idmap_csv(self):
        # arrange
        path_output = TPLS.MC_TESTS_OUTPUT.value
        filename = ''
        id_map = {}
        # action
        path_csv_file = IdProt._write_idmap_csv(path_output, filename, id_map)
        # assert
        self.assertEqual()
        self.assertTrue(os.path.exists(path_csv_file))