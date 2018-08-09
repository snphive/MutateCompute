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
        cls.DIR_PDB_1_A = '1_A'
        cls.PATH_FASTA_1_A = TPLS.MC_TESTS_INPUT.value + '/' + cls.DIR_PDB_1_A + '/' + cls.FASTAFILE_1_A
        cls.FASTA_SEQ_1_A = TPLS.FASTA_SEQ_1_A.value
        cls.NAME_1_A = cls.FASTAFILE_1_A.split('.')[0]
        cls.XML_FILE_1_A = cls.NAME_1_A + '.xml'
        cls.FASTA_STR_1_A = ">" + cls.NAME_1_A + "\n" + cls.FASTA_SEQ_1_A

        #     mutant 1_B     ##################################################################################

        cls.FASTAFILE_1_B = '1_B.fasta'
        cls.DIR_PDB_1_B = '1_B'
        cls.PATH_FASTA_1_B = TPLS.MC_TESTS_INPUT.value + '/' + cls.DIR_PDB_1_B + '/' + cls.FASTAFILE_1_B
        cls.FASTA_SEQ_1_B = TPLS.FASTA_SEQ_1_B.value

        cls.NAME_1_B = cls.FASTAFILE_1_B.split('.')[0]
        cls.XML_FILE_1_B = cls.NAME_1_B + '.xml'
        cls.FASTA_STR_1_B = ">" + cls.NAME_1_B + "\n" + cls.FASTA_SEQ_1_B

        #     mutant 2_A     ##################################################################################

        cls.FASTAFILE_2_A = '2_A.fasta'
        cls.DIR_PDB_2_A = '2_A'
        cls.PATH_FASTA_2_A = TPLS.MC_TESTS_INPUT.value + '/' + cls.DIR_PDB_2_A + '/' + cls.FASTAFILE_2_A
        cls.FASTA_SEQ_2_A = TPLS.FASTA_SEQ_2_A.value
        cls.NAME_2_A = cls.FASTAFILE_2_A.split('.')[0]
        cls.XML_FILE_2_A = cls.NAME_2_A + '.xml'
        cls.FASTA_STR_2_A = ">" + cls.NAME_2_A + "\n" + cls.FASTA_SEQ_2_A

        #     mutant 3_A     ##################################################################################

        cls.FASTAFILE_3_A = '3_A.fasta'
        cls.DIR_PDB_3_A = '3_A'
        cls.PATH_FASTA_3_A = TPLS.MC_TESTS_INPUT.value + '/' + cls.DIR_PDB_3_A + '/' + cls.FASTAFILE_3_A
        cls.FASTA_SEQ_3_A = TPLS.FASTA_SEQ_3_A.value
        cls.NAME_3_A = cls.FASTAFILE_3_A.split('.')[0]
        cls.XML_FILE_3_A = cls.NAME_3_A + '.xml'
        cls.FASTA_STR_3_A = ">" + cls.NAME_3_A + "\n" + cls.FASTA_SEQ_3_A

        #     mutant 3_B     ##################################################################################

        cls.FASTAFILE_3_B = '3_B.fasta'
        cls.DIR_PDB_3_B = '3_B'
        cls.PATH_FASTA_3_B = TPLS.MC_TESTS_INPUT.value + '/' + cls.DIR_PDB_3_B + '/' + cls.FASTAFILE_3_B
        cls.FASTA_SEQ_3_B = TPLS.FASTA_SEQ_3_B.value
        cls.NAME_3_B = cls.FASTAFILE_3_B.split('.')[0]
        cls.XML_FILE_3_B = cls.NAME_3_B + '.xml'
        cls.FASTA_STR_3_B = ">" + cls.NAME_3_B + "\n" + cls.FASTA_SEQ_3_B

        #     mutant 10_B     ##################################################################################

        cls.FASTAFILE_10_B = '10_B.fasta'
        cls.DIR_PDB_10_B = '10_B'
        cls.PATH_FASTA_10_B = TPLS.MC_TESTS_INPUT.value + '/' + cls.DIR_PDB_10_B + '/' + cls.FASTAFILE_10_B
        cls.FASTA_SEQ_10_B = TPLS.FASTA_SEQ_10_B.value
        cls.NAME_10_B = cls.FASTAFILE_10_B.split('.')[0]
        cls.XML_FILE_10_B = cls.NAME_10_B + '.xml'
        cls.FASTA_STR_10_B = ">" + cls.NAME_10_B + "\n" + cls.FASTA_SEQ_10_B
        ###############################################################################################################

        sequence_dict = {cls.XML_FILE_1_A: cls.FASTA_STR_1_A, cls.XML_FILE_1_B: cls.FASTA_STR_1_B,
                         cls.XML_FILE_2_A: cls.FASTA_STR_2_A, cls.XML_FILE_3_A: cls.FASTA_STR_3_A,
                         cls.XML_FILE_3_B: cls.FASTA_STR_3_B, cls.XML_FILE_10_B: cls.FASTA_STR_10_B}
        TestIdProt._build_reference_blastp_output_xml_files(sequence_dict)

        cls.PATH_TESTS_REFFILES_BLASTP_1_A_XML = TestIdProt.PATH_TESTS_REFFILES_BLASTP + '/' + cls.XML_FILE_1_A


    @staticmethod
    def _build_reference_blastp_output_xml_files(sequence_dict):
        for blast_out_xml, fasta_input in sequence_dict.items():
            if not os.path.exists(TestIdProt.PATH_TESTS_REFFILES_BLASTP + '/' + blast_out_xml):
                result_handle = NCBIWWW.qblast(program=Biopy.BlastParam.BLST_P.value,
                                               database=Biopy.BlastParam.SWSPRT.value,
                                               sequence=fasta_input,
                                               entrez_query=Biopy.BlastParam.HOMSAP_ORG.value,
                                               alignments=Biopy.BlastParam.MAX_ALIGN_20.value,
                                               hitlist_size=Biopy.BlastParam.MAX_HIT_20.value)
                HM.write_blastp_to_tests_dir(TPLS.MC_TESTS_REFFILES.value, TestIdProt.DIR_BLASTP,
                                             result_handle, blast_out_xml)

    @patch.object(IdProt, "_write_blastp_xml_result")
    @patch.object(NCBIWWW, "qblast")
    @patch('builtins.open', create=True)
    def test_start_1_A(self, mock_open, mock_qblast, mock__write_blastp_xml_result):
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

        # I was attempting to mock the open() and read() functions to simply return the fasta_str but this doesn't seem
        # to quite achieve it. Nonetheless it's useful to mock this as the method won't crash if the fasta file
        # it is looking for doesn't exist at the specified path. Furthermore the fasta string value is used in the
        # start() method to pass to biopython's qblast method which is mocked here anyway so it doesn't matter what
        # value is passed to it.
        mock.mock_open(read_data="ABC")
        with open(self.PATH_TESTS_REFFILES_BLASTP_1_A_XML) as test_1_A_xml:
            mock_qblast.return_value = test_1_A_xml
        mock__write_blastp_xml_result.return_value = self.PATH_TESTS_REFFILES_BLASTP_1_A_XML
        path_output = TPLS.MC_TESTS_OUTPUT.value + '/' + self.DIR_PDB_1_A + TestIdProt.DIR_BLASTP
        # action
        result_dict = IdProt.start(path_fastafile=self.PATH_FASTA_1_A, write_blastp_json=False, build_idmap_csv=False,
                                   path_output=path_output)
        # assert
        self.assertEqual(expected_qblast_dict_1_A, result_dict)


    def test__write_blastp_xml_result(self):
        # arrange
        path_output = TPLS.MC_TESTS_OUTPUT.value
        filename = ''
        blastp_result = {}
        # action
        IdProt._write_blastp_xml_result(path_output, filename, blastp_result)
        # assert

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

    def test__build_idmap(self):
        # arrange
        blastp_result_dict = {}
        filename = ''
        # action
        IdProt._build_idmap(blastp_result_dict, filename)
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