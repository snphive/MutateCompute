from unittest import TestCase
from src.GeneralUtilityMethods import GUM
from tests.HelperMethods import HM
import os
from unittest.mock import patch
from tests.TestPathsAndLists import TPL
import subprocess


class TestGUM(TestCase):

    # Copies over all configuration & input data from the main directory into the /tests/ before running tests.
    @classmethod
    def setUpClass(cls):
        if not os.path.exists(TPL.MC_TESTS_CONFIG.value):
            GUM.linux_copy(path_src=TPL.CONFIG_FOR_READ_ONLY.value, path_dst=TPL.MC_TESTS.value,
                           do_recursively=True)

        if not os.path.exists(TPL.MC_TESTS_INPUT.value):
            GUM.linux_copy(path_src=TPL.INPUT_FOR_READ_ONLY.value, path_dst=TPL.MC_TESTS.value,
                           do_recursively=True)

    # Write_runscript_for_pdbs() takes 6 arguments. The last 3 (namely show_sequence_detail, print_networks,
    # calculate_stability) are keyword (named) arguments. All have default values assigned in the method argument so
    # that when the method is called without supplying a (keyword-named) value, the default value is applied.
    # (The default values for all 3 are False).
    #
    # 30.07.18 Have redesigned the directory structure such that runscripts will go in configuration/foldx/and maybe
    # another level such as analyse_complex or build_model or stability etc.
    def test_write_runscript_for_pdbs(self):
        # arrange
        pdb = 'RepairPDB_1.pdb'
        path_runscript = TPL.MC_TESTS_CONFIG_FX.value
        if not os.path.exists(path_runscript):
            os.makedirs(path_runscript)
        action = '<BuildModel>#,individual_list.txt'
        expected_runscript = '<TITLE>FOLDX_runscript;\n' + '<JOBSTART>#;\n' + '<PDBS>' + pdb + ';\n' + \
                             '<BATCH>#;\n' + '<COMMANDS>FOLDX_commandfile;\n' + action + ';\n' + '<END>#;\n' + \
                             '<OPTIONS>FOLDX_optionfile;\n' + '<Temperature>298;\n' + '<IonStrength>0.05;\n' + \
                             '<ph>7;\n' + '<moveNeighbours>true;\n' + '<VdWDesign>2;\n' + '<numberOfRuns>3;\n' + \
                             '<OutPDB>#;\n' + '<END>#;\n' + '<JOBEND>#;\n' + '<ENDFILE>#;\n'
        single_space = ' '
        not_expected_runscript = single_space + expected_runscript
        # action
        actual_runscript = GUM.write_runscript_for_pdbs(path_runscript, pdb, action)
        # assert
        self.assertEqual(actual_runscript, expected_runscript)
        self.assertNotEqual(actual_runscript, not_expected_runscript)

    # CHECK WHAT HAPPENS IF THERE ARE FEWER FILES THAN SPECIFIED TO MOVE BY TOTAL_NUM_TO_COPY
    # @patch.object(GUM, '_make_subfoldername')
    # def test_copy_and_move_pdb_files(self):
    #     # arrange
    #     path_src_dir = TPL.MC_TESTS_INPUT.value
    #     path_dst_dir = TPL.MC_TESTS_OUTPUT.value
    #     starting_num = 1
    #     # action
    #     GUM.copy_and_move_pdb_files(path_src_dir, path_dst_dir, starting_num, total_num_to_copy=9)
    #     # assert
    #     # expected_call_count_1 = 1
    #     # self.assertEqual(expected_call_count_1, mock_make_subfoldername.call_count)
    def test_copy_input_files_from_repo_to_input(self):
        # arrange
        path_src_repo_dir = TPL.REPO_PDB_FASTA.value
        path_dst_dir = TPL.MC_TESTS_INPUT.value
        src_file_list = ['RepairPDB_1.pdb', 'RepairPDB_2.pdb', 'RepairPDB_3.pdb', 'RepairPDB_4.pdb', 'RepairPDB_5.pdb']
        wanted_file_list = ['RepairPDB_1.pdb', 'RepairPDB_2.pdb']
        # action
        copied_wanted_file_list = GUM.copy_files_from_repo_to_input_filedir(path_src_repo_dir, path_dst_dir,
                                                                            src_file_list, wanted_file_list,
                                                                            copy_all_files=False)
        path_copied_file_list = [path_dst_dir + '/' + x for x in copied_wanted_file_list]

        # assert
        self.assertEqual(wanted_file_list, copied_wanted_file_list)
        for path_copied_file in path_copied_file_list:
            self.assertTrue(path_copied_file)

    # Same as test_copy_input_files_from_repo_to_input() above but for fasta files.
    def test_copy_input_files_from_repo_to_input_2(self):
        # arrange
        path_src_repo_dir = TPL.REPO_PDB_FASTA.value
        path_dst_dir = TPL.MC_TESTS_INPUT.value
        src_file_list = ['1_A.fasta', '1_B.fasta', '2_A.fasta', '3_A.fasta', '3_B.fasta']
        wanted_file_list = ['1_A.fasta', '3_A.fasta']
        # action
        copied_wanted_file_list = GUM.copy_files_from_repo_to_input_filedir(path_src_repo_dir, path_dst_dir,
                                                                            src_file_list, wanted_file_list,
                                                                            copy_all_files=False)
        path_copied_file_list = [path_dst_dir + '/' + x for x in copied_wanted_file_list]
        # assert
        self.assertEqual(wanted_file_list, copied_wanted_file_list)
        for path_copied_file in path_copied_file_list:
            self.assertTrue(path_copied_file)

    def test_extract_pdbname_chain_fasta_from_pdb(self):
        # arrange
        pdbfiles = ['RepairPDB_1.pdb', 'RepairPDB_2.pdb']
        path_input = TPL.MC_TESTS_INPUT.value
        write_fastafile = False
        path_to_write_fastafile = ''
        expected_pdbname_chain_fastaseq_dict = {
        'RepairPDB_1_A': 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGNFVR'
                         'VIQTFNRTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRSLTKRNA'
                         'VRTDQHNSKWLSEPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARLVCSVTDEDG'
                         'PETHFDELEDVFLLETDNPRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRPGTCPGGAFTPNM'
                         'RTTKEFPDDVVTFIRNHPLMYNSIYPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVLPTNNSVSGELILEEL'
                         'EVFKNHAPITTMKISSKKQQLYVSSNEGVSQVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRRSRRQDVRHGNPLTQCR',
        'RepairPDB_1_B': 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGNFVR'
                         'VIQTFNRTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRSLTKRNA'
                         'VRTDQHNSKWLSEPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARLVCSVTDEDG'
                         'PETHFDELEDVFLLETDNPRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRPGTCPGGAFTPNM'
                         'RTTKEFPDDVVTFIRNHPLMYNSIYPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVLPTNNSVSGELILEEL'
                         'EVFKNHAPITTMKISSKKQQLYVSSNEGVSQVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRRSRRQDVRHGNPLTQCR',
        'RepairPDB_2_A': 'EIVQYGVKNNTTFLECAPKSPQASIKWLLQKDKDRRKEVKLNERIIATSQGLLIRSVQGSDQGLYHCIATENSFKQTIAKINFKVLD'}
        # action
        pdbname_chain_fastaseq_dict = GUM.extract_pdbname_chain_fasta_from_pdbs(pdbfiles, path_input, write_fastafile,
                                                                            path_to_write_fastafile)
        # assert
        self.assertDictEqual(expected_pdbname_chain_fastaseq_dict, pdbname_chain_fastaseq_dict)

    def test_extract_all_chains_from_pdb(self):
        # arrange
        pdbfile = 'RepairPDB_1.pdb'
        path_pdbfile = TPL.MC_TESTS_INPUT.value + '/' + pdbfile.split('.')[0]
        expected_protein_chains = ['A', 'B']
        # action
        protein_chains = GUM.extract_all_chains_from_pdb(pdbfile, path_pdbfile)
        # assert
        self.assertEqual(expected_protein_chains, protein_chains)

    def test_remove_prefix_and_suffix(self):
        # arrange
        name = 'RepairPDB_1_A.pdb'
        prefix = 'RepairPDB_'
        suffix = '.pdb'
        expected_word = '1_A'
        # action
        word = GUM._remove_prefix_and_suffix(prefix, suffix, name)
        # assert
        self.assertEqual(expected_word, word)

    def test_get_sequenceonly_from_fastafile(self):
        # arrange
        path_fastafile = TPL.MC_TESTS_INPUT.value + '/fastas/1_A/1_A.fasta'
        expected_seq = \
        'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGNFVRVIQTFNRTHLYVCGSGAF' \
        'SPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRSLTKRNAVRTDQHNSKWLSEPMFVDAHVIPDGTDPNDAKVYFF' \
        'FKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARLVCSVTDEDGPETHFDELEDVFLLETDNPRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNG' \
        'PFAHKEGPNHQLISYQGRIPYPRPGTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSIYPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKV' \
        'VVLPTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVSQVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRRSRRQDVRHGNPLTQCR'
        # action
        sequence = GUM.get_sequenceonly_from_fastafile(path_fastafile)
        # assert
        self.assertEqual(expected_seq, sequence)
