from unittest import TestCase
from Bio.Blast import NCBIWWW
from src.Biopython import Biopy
from unittest.mock import patch


class TestBiopython(TestCase):

    def setUp(self):
        self.path_fasta_input_test_seq1 = '/Users/u0120577/ROB_HOMOLOGY_HUMAN/Fasta/1_A.fasta'
        self.path_fasta_input_test_seq2 = '/Users/u0120577/ROB_HOMOLOGY_HUMAN/Fasta/1_B.fasta'
        self.path_fasta_input_test_seq3 = '/Users/u0120577/ROB_HOMOLOGY_HUMAN/Fasta/2_A.fasta'
        self.path_fasta_input_test_seq4 = '/Users/u0120577/ROB_HOMOLOGY_HUMAN/Fasta/3_A.fasta'
        self.path_fasta_input_test_seq5 = '/Users/u0120577/ROB_HOMOLOGY_HUMAN/Fasta/3_B.fasta'

        self.path_output_root = '/Users/u0120577/SNPEFFECT/Outputs'

    def tearDown(self):
        self.fasta_repaired_1_A = None

    @patch.object(Biopy, '_write_qblast_xml_result')
    @patch.object(NCBIWWW, 'qblast')
    @patch.object(Biopy, '_read_input_fasta_seq')
    def test_find_identical_blastp_hit_swissprot_for_fasta_input(self, mock__read_input_fasta_seq, mock_qblast,
                                                                 mock__write_qblast_xml_result):
        # arrange
        mock__read_input_fasta_seq.return_value = '>1_A\nRVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNIS' \
                                                  'QEALSVFWPASTIKVEECKMAGKDPTHGCGNFVRVIQTFNRTHLYVCGSGAFSPVCTYLNRG' \
                                                  'RRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRSLTKRNA' \
                                                  'VRTDQHNSKWLSEPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPND' \
                                                  'TGGLRSLVNKWTTFLKARLVCSVTDEDGPETHFDELEDVFLLETDNPRTTLVYGIFTTSSSV' \
                                                  'FKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRPGTCPGGAFTPNMRTTKEF' \
                                                  'PDDVVTFIRNHPLMYNSIYPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGT' \
                                                  'VQKVVVLPTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVSQVSLHRCHIY' \
                                                  'GTACADCCLARDPYCAWDGHSCSRFYPTGKRRSRRQDVRHGNPLTQCR'

        with open('/Users/u0120577/SNPEFFECT/Outputs/test_blastp_1_A.xml') as test_1_A_xml:
            mock_qblast.return_value = test_1_A_xml
        mock__write_qblast_xml_result.return_value = '/Users/u0120577/SNPEFFECT/Outputs/blastp_swissprot/1_A.xml'

        path_fasta_input_query_seq = self.path_fasta_input_test_seq1
        path_output_root = self.path_output_root

        # action
        result = Biopy.find_identical_blastp_hit_swissprot_for_fasta_input(path_fasta_input_query_seq,
                                                                          path_output_root)
        print('nothing: ' + result)
        # assert
        # self.assertEqual()

