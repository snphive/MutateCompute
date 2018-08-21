import os
from unittest import TestCase
from src.Biopython import Biopy
from unittest.mock import patch
from tests.HelperMethods import HM
from tests.TestPathsAndListsSeqs import TPLS
from Bio.Blast import NCBIXML


class TestBiopython(TestCase):

    @classmethod
    def setUpClass(cls):
        # constants related to input
        cls.FASTAFILE_1_A = '1_A.fasta'
        cls.DIR_PDB_1_A = '1_A'
        cls.PATH_FASTA_1_A = os.path.join(TPLS.MC_TESTS_INPUT.value, cls.DIR_PDB_1_A, cls.FASTAFILE_1_A)
        cls.FASTA_SEQ_1_A = 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGN' \
                            'FVRVIQTFNRTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRS' \
                            'LTKRNAVRTDQHNSKWLSEPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARL' \
                            'VCSVTDEDGPETHFDELEDVFLLETDNPRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRP' \
                            'GTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSIYPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVL' \
                            'PTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVSQVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRR' \
                            'SRRQDVRHGNPLTQCR'
        cls.NAME_1_A = cls.FASTAFILE_1_A.split('.')[0]
        cls.XML_FILE_1_A = cls.NAME_1_A + '.xml'
        cls.FASTA_STR_1_A = ">" + cls.NAME_1_A + "\n" + cls.FASTA_SEQ_1_A

        # constants related to output
        cls.DIR_BLASTP = 'blastp'
        cls.REL_BLASTP = '/' + cls.DIR_BLASTP
        cls.XML_1_A_BLASTP_OUTPUT_FILE = '1_A.xml'

        # constants related to blastp parameters
        cls.SWSPRT_DB = 'swissprot'
        cls.SWSPRT_PROTS_NUM = 20341
        cls.ZERO_GAP = 0
        cls.QSTRT_1 = 1

    # This test is not ideal. It relies on running the actual qblast and writing the result to an xml
    # file, if only once. It then runs the qblast method in src.Biopython.Biopy but also has to write the
    # result to an xml in the test method. Hence it does not exclusively test the _run_blastp method. It
    # can fail if writing the result or reading the pre-written result fails. (I cannot find a more direct
    # mechanism for file comparison when the file is still in the io file format.)
    # NOTE: using filecmp.cmp() is useless because each blastp run has a unique query id which is written in the output
    def test_run_blastp_1_A(self):
        # arrange
        Hsp_eval = "Hsp_evalue"
        path_actual = os.path.join(TPLS.MC_TESTS_OUTPUT.value + self.REL_BLASTP, self.XML_1_A_BLASTP_OUTPUT_FILE)
        path_expected = os.path.join(TPLS.MC_TESTS_REF_FILES.value, self.REL_BLASTP, self.XML_1_A_BLASTP_OUTPUT_FILE)
        # act
        result_handle_1_A = Biopy.run_blastp(fasta_str=self.FASTA_STR_1_A)
        blastp_output_xml = self.XML_1_A_BLASTP_OUTPUT_FILE
        HM.write_blastp_to_tests_dir(path_tests=TPLS.MC_TESTS_OUTPUT.value, blastp_dir=self.DIR_BLASTP,
                                     result_handle=result_handle_1_A, blastp_output_xmlfile=blastp_output_xml)
        # assert
        with open(path_actual) as actual, open(path_expected) as expected:
            for actual_line, exp_line in zip(actual, expected):
                if 'Iteration_query-ID' in exp_line or 'BlastOutput_query-ID' in exp_line:
                    self.assertNotEqual(exp_line, actual_line)
                elif 'Statistics' in exp_line:
                    break
                elif Hsp_eval in exp_line:
                    self.assertTrue(TestBiopython._has_insignificant_diff_Hsp_evalues(exp_line, actual_line, Hsp_eval))
                else:
                    self.assertEqual(exp_line, actual_line)

    def test_parse_and_filter_blastp_xml_to_dict(self):
        self.fail()

    # NOTE: The difference cut-off of 0.0001 was arbitrarily chosen, may change this at some point.
    @staticmethod
    def _has_insignificant_diff_Hsp_evalues(exp_evalue_line, actual_evalue_line, str_tag):
        exp_evalue = TestBiopython.__extract_value_from_tagged_line(line=exp_evalue_line, tag=str_tag)
        actual_evalue = TestBiopython.__extract_value_from_tagged_line(line=actual_evalue_line, tag=str_tag)
        return abs(exp_evalue - actual_evalue) < 0.0001

    @staticmethod
    def __extract_value_from_tagged_line(line, tag):
        start_tag = "<" + tag + ">"
        end_tag = "</" + tag + ">"
        value_start_index = line.find(start_tag) + len(start_tag)
        if value_start_index == -1:
            raise ValueError(start_tag + " not found in: " + line)
        value_end_index = line.find(end_tag)
        if value_end_index == -1:
            raise ValueError(end_tag + " not found in: " + line)
        return float(line[value_start_index:value_end_index])

    # def test__make_list_of_dicts_of_hsps_with_0gaps_and_queryLen_equal_alignLen(self):
    #     # arrange
    #     # might need to mock these objects, this is not working
    #     class Hsp(object):
    #         def __init(self):
    #             self.align_length = self.gaps = self.identities = self.query_end = self.query_start = self.sbjct_end
    # = self.sbjct_start = 0
    #             self.expect = 0.0
    #
    #     class Alignment(object):
    #         def __init(self):
    #             self.accession = ''
    #             self.hit_def = ''
    #             self.length = 0
    #             self.hsps = []
    #
    #
    #     hsp1 = Hsp()
    #     hsp1.align_length = hsp1.identities = hsp1.query_end = hsp1.sbjct_end = 100
    #     hsp1.query_start = hsp1.sbjct_start = 1
    #     alignment1 = Alignment()
    #     alignment1.hsps.append(Hsp)
    #     alignment1.accession = 'X1234'
    #     alignment1.hit_def = 'hit_X1234'
    #     alignment1.length = 100
    #     alignments = []
    #     alignments.append(alignment1)
    #
    #     hsp2 = Hsp()
    #     hsp2.align_length = hsp2.identities = hsp2.query_end = hsp2.sbjct_end = 101
    #     hsp2.query_start = hsp2.sbjct_start = 2
    #     alignment2 = Alignment()
    #     alignment2.hsps.append(Hsp)
    #     alignment2.accession = 'A4321'
    #     alignment2.hit_def = 'hit_A4321'
    #     alignment2.length = 100
    #     alignments.append(alignment2)
    #
    #     query_length = 100
    #     # act
    #     Biopy._make_list_of_dicts_of_hsps_with_0gaps_and_queryLen_equal_alignLen(query_length, alignments)
    #     # assert
