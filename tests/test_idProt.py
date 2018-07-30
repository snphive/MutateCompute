from unittest import TestCase
from unittest.mock import patch
from src.Biopython import Biopy
import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from tests.HelperMethods import HM
from src.IdentifyProtein import IdProt
from tests.PathsForTests import PFT


class TestIdProt(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.PATH_FASTA_INPUT_FILES = '/Users/u0120577/PycharmProjects/MutateCompute/tests/Inputs/Fasta'
        cls.PATH_TESTS_OUTPUTS = '/Users/u0120577/PycharmProjects/MutateCompute/tests/Outputs'
        cls.DIR_BLAST_SP_HS_20_20 = 'blastp_sp_hs_20_20'

        cls.SWSPRT_DB = 'swissprot'
        cls.SWSPRT_PROTS_NUM = 20341
        cls.ZERO_GAP = 0
        cls.QSTRT_1 = 1

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

        cls.FASTA_1_B_FILE = '1_B.fasta'
        cls.PATH_FASTA_1_B = cls.PATH_FASTA_INPUT_FILES + '/' + cls.FASTA_1_B_FILE
        cls.FASTA_1_B_SEQ = 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGN' \
                            'FVRVIQTFNRTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRS' \
                            'LTKRNAVRTDQHNSKWLSEPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARL' \
                            'VCSVTDEDGPETHFDELEDVFLLETDNPRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRP' \
                            'GTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSIYPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVL' \
                            'PTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVSQVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRR' \
                            'SRRQDVRHGNPLTQCR'
        cls.FASTA_1_B_NAME = cls.FASTA_1_B_FILE.split('.')[0]
        cls.FASTA_1_B_XML = cls.FASTA_1_B_NAME + '.xml'
        cls.FASTA_1_B_STR = ">" + cls.FASTA_1_B_NAME + "\n" + cls.FASTA_1_B_SEQ

        cls.FASTA_2_A_FILE = '2_A.fasta'
        cls.PATH_FASTA_2_A = cls.PATH_FASTA_INPUT_FILES + '/' + cls.FASTA_2_A_FILE
        cls.FASTA_2_A_SEQ = 'EIVQYGVKNNTTFLECAPKSPQASIKWLLQKDKDRRKEVKLNERIIATSQGLLIRSVQGSDQGLYHCIATENSFKQTIAKINFKVLD'
        cls.FASTA_2_A_NAME = cls.FASTA_2_A_FILE.split('.')[0]
        cls.FASTA_2_A_XML = cls.FASTA_2_A_NAME + '.xml'
        cls.FASTA_2_A_STR = ">" + cls.FASTA_2_A_NAME + "\n" + cls.FASTA_2_A_SEQ

        cls.FASTA_3_A_FILE = '3_A.fasta'
        cls.PATH_FASTA_3_A = cls.PATH_FASTA_INPUT_FILES + '/' + cls.FASTA_3_A_FILE
        cls.FASTA_3_A_SEQ = 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGN' \
                            'FVRVIQTFNRTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRS' \
                            'LTKRNAVRTDQHNSKWLSEPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARL' \
                            'VCSVTDEDGPETHFDELEDVFLLETDNPRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRP' \
                            'GTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSIYPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVL' \
                            'PTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVSQVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRR' \
                            'SRRQDVRHGNPLTQCR'
        cls.FASTA_3_A_NAME = cls.FASTA_3_A_FILE.split('.')[0]
        cls.FASTA_3_A_XML = cls.FASTA_3_A_NAME + '.xml'
        cls.FASTA_3_A_STR = ">" + cls.FASTA_3_A_NAME + "\n" + cls.FASTA_3_A_SEQ

        cls.FASTA_3_B_FILE = '3_B.fasta'
        cls.PATH_FASTA_3_B = cls.PATH_FASTA_INPUT_FILES + '/' + cls.FASTA_3_B_FILE
        cls.FASTA_3_B_SEQ = 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGN' \
                            'FVRVIQTFNRTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRS' \
                            'LTKRNAVRTDQHNSKWLSEPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARL' \
                            'VCSVTDEDGPETHFDELEDVFLLETDNPRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRP' \
                            'GTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSIYPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVL' \
                            'PTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVSQVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRR' \
                            'SRRQDVRHGNPLTQCR'
        cls.FASTA_3_B_NAME = cls.FASTA_3_B_FILE.split('.')[0]
        cls.FASTA_3_B_XML = cls.FASTA_3_B_NAME + '.xml'
        cls.FASTA_3_B_STR = ">" + cls.FASTA_3_B_NAME + "\n" + cls.FASTA_3_B_SEQ

        cls.FASTA_10_B_FILE = '10_B.fasta'
        cls.PATH_FASTA_10_B = cls.PATH_FASTA_INPUT_FILES + '/' + cls.FASTA_10_B_FILE
        cls.FASTA_10_B_SEQ = 'SCIQFTRHASDVLLNLNRLRSRDILTDVVIVVSREQFRAHKTVLMACSGLFYSIFTDQLKCNLSVINLDPEINPEGFCILLDFMYT' \
                             'SRLNLREGNIMAVMATAMYLQMEHVVDTCRKFIKAS'
        cls.FASTA_10_B_NAME = cls.FASTA_10_B_FILE.split('.')[0]
        cls.FASTA_10_B_XML = cls.FASTA_10_B_NAME + '.xml'
        cls.FASTA_10_B_STR = ">" + cls.FASTA_10_B_NAME + "\n" + cls.FASTA_10_B_SEQ

        cls.PATH_TESTS_FILES_BLASTRUN = PFT.PATH_TESTS_REF_FILES + '/' + cls.DIR_BLAST_SP_HS_20_20
        cls.PATH_TESTS_FILES_BLASTRUN_1_A_XML = cls.PATH_TESTS_FILES_BLASTRUN + '/' + cls.FASTA_1_A_XML
        sequence_dict = {cls.FASTA_1_A_XML: cls.FASTA_1_A_STR, cls.FASTA_1_B_XML: cls.FASTA_1_B_STR,
                         cls.FASTA_2_A_XML: cls.FASTA_2_A_STR, cls.FASTA_3_A_XML: cls.FASTA_3_A_STR,
                         cls.FASTA_3_B_XML: cls.FASTA_3_B_STR, cls.FASTA_10_B_XML: cls.FASTA_10_B_STR}
        TestIdProt._build_reference_blastp_output_xml_files(sequence_dict)

    @classmethod
    def tearDownClass(cls):
        cls.PATH_FASTA_INPUT_FILES = None
        cls.SWSPRT_DB = None
        cls.SWSPRT_PROTS_NUM = None
        cls.ZERO_GAP = None
        cls.QSTRT_1 = None
        cls.PATH_TESTS_OUTPUTS = None
        cls.FASTA_1_A_FILE = None
        cls.PATH_FASTA_1_A = None
        cls.FASTA_1_A_SEQ = None
        cls.FASTA_1_A_NAME = None
        cls.FASTA_1_A_XML = None
        cls.FASTA_1_A_STR = None
        cls.PATH_TESTS_FILES = None
        cls.DIR_BLAST_SP_HS_20_20 = None

    @staticmethod
    def _build_reference_blastp_output_xml_files(sequence_dict):
        for blast_out_xml, fasta_input in sequence_dict.items():
            if not os.path.exists(TestIdProt.PATH_TESTS_FILES_BLASTRUN + '/' + blast_out_xml):
                result_handle = NCBIWWW.qblast(program=Biopy.BlastParam.BLST_P.value,
                                               database=Biopy.BlastParam.SWSPRT.value,
                                               sequence=fasta_input,
                                               entrez_query=Biopy.BlastParam.HOMSAP_ORG.value,
                                               alignments=Biopy.BlastParam.MAX_ALIGN_20.value,
                                               hitlist_size=Biopy.BlastParam.MAX_HIT_20.value)
                THM.write_blast_run_to_tests_dir(TestIdProt.PATH_TESTS_FILES, TestIdProt.DIR_BLAST_SP_HS_20_20,
                                                 result_handle, blast_out_xml)



    @patch.object(Biopy, "_write_qblast_xml_result")
    @patch.object(NCBIWWW, "qblast")
    @patch.object(Biopy, "_read_input_fasta_seq")
    def test_start_1_A(self, mock__read_input_fasta_seq, mock_qblast,
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
        with open(self.PATH_TESTS_FILES_BLASTRUN_1_A_XML) as test_1_A_xml:
            mock_qblast.return_value = test_1_A_xml
        mock__write_qblast_xml_result.return_value = self.PATH_TESTS_FILES_BLASTRUN_1_A_XML
        # action
        result_dict = IdProt.start(self.PATH_FASTA_1_A, self.PATH_TESTS_OUTPUTS)
        # assert
        self.assertEqual(expected_qblast_dict_1_A, result_dict)


    def test__read_fasta_file(self):
        self.fail()

    def test__write_blastp_xml_result(self):
        self.fail()

    def test__write_dict_to_json_file(self):
        self.fail()

    def test__build_idmap(self):
        self.fail()

    def test__write_idmap_csv(self):
        self.fail()
