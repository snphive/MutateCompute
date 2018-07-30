from unittest import TestCase
from src.Biopython import Biopy
from unittest.mock import patch
from tests.HelperMethods import HM


class TestBiopython(TestCase):

    # This test is not ideal. It relies on running the actual qblast and writing the result to an xml
    # file, if only once. It then runs the qblast method in src.Biopython.Biopy but also has to write the
    # result to an xml in the test method. Hence it does not exclusively test the _run_blastp method. It
    # can fail if writing the result or reading the pre-written result fails. (I cannot find a more direct
    # mechanism for file comparison when the file is still in the io file format.)
    # NOTE: using filecmp.cmp() is useless because each blastp run has a unique query id which is written in the output
    def test_run_blastp_1_A(self):
        # arrange
        Hsp_eval = "Hsp_evalue"
        path_actual = self.PATH_TESTS_OUTPUTS + '/' + self.DIR_BLAST_SP_HS_20_20 + '/' + self.FASTA_1_A_XML
        path_expected = self.PATH_TESTS_FILES + '/' + self.DIR_BLAST_SP_HS_20_20 + '/' + self.FASTA_1_A_XML
        # action
        result_handle_1_A = Biopy._run_blastp(self.FASTA_1_A_STR)
        blast_output_xml = self.FASTA_1_A_XML
        HM.write_blast_run_to_tests_dir(self.PATH_TESTS_OUTPUTS, TestBiopython.DIR_BLAST_SP_HS_20_20,
                                         result_handle_1_A, blast_output_xml)
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

    # NOTE: The difference cut-off of 0.0001 was arbitrarily chosen, may change this at some point.
    @staticmethod
    def _has_insignificant_diff_Hsp_evalues(exp_evalue_line, actual_evalue_line, str_tag):
        exp_evalue = TestBiopython.__extract_value_from_tagged_line(exp_evalue_line, str_tag)
        actual_evalue = TestBiopython.__extract_value_from_tagged_line(actual_evalue_line, str_tag)
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

