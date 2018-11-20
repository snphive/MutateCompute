from unittest import TestCase
from src.tools.GeneralUtilityMethods import GUM
import os
import glob
from src.enums.Paths import Paths
from unittest.mock import patch
from tests.TestPathsAndListsSeqs import TPLS
from src.enums.Str import Str


class TestGUM(TestCase):

    # Copies over all configuration & input data from the main directory into the /tests/ before running tests.
    @classmethod
    def setUpClass(cls):
        if not os.path.exists(TPLS.MC_TESTS_CONFIG.value):
            GUM.linux_copy_all_files_in_dir(path_src_dir=TPLS.CONFIG_FOR_READ_ONLY.value,
                                            path_dst_dir=TPLS.MC_TESTS.value, recursively=True)

        if not os.path.exists(TPLS.MC_TESTS_INPUT.value):
            GUM.linux_copy_all_files_in_dir(path_src_dir=TPLS.INPUT_FOR_READ_ONLY.value,
                                            path_dst_dir=TPLS.MC_TESTS.value, recursively=True)

    # CHECK WHAT HAPPENS IF THERE ARE FEWER FILES THAN SPECIFIED TO MOVE BY TOTAL_NUM_TO_COPY
    # @patch.object(GUM, '_make_subfoldername')
    # def test_copy_and_move_pdb_files(self):
    #     # arrange
    #     path_src_dir = TPL.MC_TESTS_INPUT.value
    #     path_dst_dir = TPL.MC_TESTS_OUTPUT.value
    #     starting_num = 1
    #     # act
    #     GUM.copy_and_move_pdb_files(path_src_dir, path_dst_dir, starting_num, total_num_to_copy=9)
    #     # assert
    #     # expected_call_count_1 = 1
    #     # self.assertEqual(expected_call_count_1, mock_make_subfoldername.call_count)
    @patch.object(GUM, 'get_pdb_or_fastafile_list_from_subdirs')
    def test_copy_files_from_repo_to_input_dst_dir(self, mock_get_pdb_or_fastafile_list_from_subdirs):
        # arrange
        path_repo = TPLS.REPO_PDB_FASTA.value
        path_dst = TPLS.MC_TESTS_INPUT.value
        wanted_file_list = ['RepairPDB_1.pdb', 'RepairPDB_2.pdb']
        mock_get_pdb_or_fastafile_list_from_subdirs.return_value = ['RepairPDB_1.pdb', 'RepairPDB_2.pdb',
                                                                'RepairPDB_3.pdb', 'RepairPDB_4.pdb', 'RepairPDB_5.pdb']
        copy_all_files_in_dir =  False
        # act
        copied_wanted_file_list = GUM.copy_files_from_repo_to_input_dirs(path_repo, path_dst, wanted_file_list)
        path_copied_file_list = [path_dst + '/' + x for x in copied_wanted_file_list]

        # assert
        self.assertEqual(wanted_file_list, copied_wanted_file_list)
        for path_copied_file in path_copied_file_list:
            self.assertTrue(path_copied_file)

    # Same as test_copy_input_files_from_repo_to_input() above but for fasta files.
    @patch.object(GUM, 'get_pdb_or_fastafile_list_from_subdirs')
    def test_copy_input_files_from_repo_to_input_2(self, mock_get_pdb_or_fastafile_list_from_subdirs):
        # arrange
        path_repo = TPLS.REPO_PDB_FASTA.value
        path_dst = TPLS.MC_TESTS_INPUT.value
        wanted_file_list = ['1_A.fasta', '3_A.fasta']
        mock_get_pdb_or_fastafile_list_from_subdirs.return_value = ['1_A.fasta', '1_B.fasta', '2_A.fasta', '3_A.fasta',
                                                                    '3_B.fasta']
        copy_all_files_in_dir = False
        # act
        copied_wanted_file_list = GUM.copy_files_from_repo_to_input_dirs(path_repo, path_dst, wanted_file_list)
        path_copied_file_list = [path_dst + '/' + x for x in copied_wanted_file_list]
        # assert
        self.assertEqual(wanted_file_list, copied_wanted_file_list)
        for path_copied_file in path_copied_file_list:
            self.assertTrue(path_copied_file)

    def test_extract_pdbname_chain_fasta_from_pdbs(self):
        # arrange
        path_pdbfiles = [os.path.join(TPLS.MC_TESTS_INPUT.value, 'RepairPDB_1', 'RepairPDB_1.pdb'),
                         os.path.join(TPLS.MC_TESTS_INPUT.value, 'RepairPDB_2', 'RepairPDB_2.pdb')]
        path_to_write_fastafile = ''
        # Need at some point to make a decision about whether to include 'RepairPDB_' in the names.
        # expected_pdbname_chain_fastaseq_dict = {
        # 'RepairPDB_1_A': TPLS.FASTA_SEQ_1_A.value,
        # 'RepairPDB_1_B': TPLS.FASTA_SEQ_1_B.value,
        # 'RepairPDB_2_A': TPLS.FASTA_SEQ_2_A.value}
        expected_pdbname_chain_fastaseq_dict = {
        '1_A_1': TPLS.FASTA_SEQ_1_A.value,
        '1_B_540': TPLS.FASTA_SEQ_1_B.value,
        '2_A_1': TPLS.FASTA_SEQ_2_A.value}
        # act
        pdbname_chain_fastaseq_dict = GUM.extract_pdbname_chain_startpos_fasta_from_pdbs(path_pdbfiles=path_pdbfiles)
        # assert
        # self.maxDiff = None
        self.assertDictEqual(expected_pdbname_chain_fastaseq_dict, pdbname_chain_fastaseq_dict)

    def test_extract_all_chains_from_pdb(self):
        # arrange
        pdbfile = 'RepairPDB_1.pdb'
        path_pdbfile = os.path.join(TPLS.MC_TESTS_INPUT.value, pdbfile.split('.')[0])
        expected_protein_chains = ['A', 'B']
        # act
        protein_chains = GUM.extract_all_chains_from_pdb(pdbfile, path_pdbfile)
        # assert
        self.maxDiff = None
        self.assertEqual(expected_protein_chains, protein_chains)

    def test_remove_prefix_and_suffix(self):
        # arrange
        input_str = 'RepairPDB_1_A.pdb'
        prefix = 'RepairPDB_'
        suffix = '.pdb'
        expected_trimmed = '1_A'
        # act
        trimmed = GUM._remove_prefix_and_suffix(input_str, prefix, suffix)
        # assert
        self.assertEqual(expected_trimmed, trimmed)

    def test_get_sequenceonly_from_fastafile(self):
        # arrange
        path_fastafile = TPLS.MC_TESTS_INPUT.value + '/fastas/1_A/1_A.fasta'
        expected_seq_1_A = TPLS.FASTA_SEQ_1_A.value
        # act
        sequence = GUM.get_sequenceOnly_from_fastafile(path_fastafile)
        # assert
        self.assertEqual(expected_seq_1_A, sequence)

    # NEED TO ADD AN ASSERT (BUT THE METHOD DOES WORK)
    def test_move_files_into_own_subdirs(self):
        # arrange
        path_repo_fastafilelist = glob.glob(TPLS.REPO_PDB_FASTA.value + '/fastas_10/*.fasta')
        GUM.linux_copy_files(path_repo_fastafilelist, path_dst=TPLS.MC_TESTS_INPUT_FASTAS.value, into_own_subdirs=False)
        # act
        GUM._move_files_into_own_subdirs(path_dir=TPLS.MC_TESTS_INPUT_FASTAS.value)
        # assert
        print('nothing')

    def test_copy_files_to_3dot_dir_1(self):
        # act
        Paths.set_up_paths(use_cluster=False)
        # GUM.copy_files_to_3dot_dir(path_src_dir='/Users/u0120577/ROB_HOMOLOGY_HUMAN/Fasta',
        GUM.copy_files_to_3dot_dir(path_src_dir='/Users/u0120577/ROB_HOMOLOGY_HUMAN/pub_hg_fasta',
                                   path_dst_dir=Paths.REPO_FASTAS + '/29611_fastas_1000',
                                   file_extension=Str.FSTAEXT.value, starting_num=1, num_to_copy_per_subdir=1000)
        # assert
        self.assertTrue(os.path.exists(os.path.join(Paths.REPO_FASTAS, '29611_fastas_1000')),
                        'No folder at  ~/REPO_PDB_FASTA/fastas/29611_fastas_1000 found')
        self.assertTrue(os.path.exists(os.path.join(Paths.REPO_FASTAS, '29611_fastas_1000', '1...1000')),
                        'No folder at  ~/REPO_PDB_FASTA/fastas/29611_fastas_1000/1...1000 found')
        self.assertTrue(os.path.exists(os.path.join(Paths.REPO_FASTAS, '29611_fastas_1000', '2001...3000')),
                        'No folder at  ~/REPO_PDB_FASTA/fastas/29611_fastas_1000/2001...3000 found')
        self.assertTrue(os.path.exists(os.path.join(Paths.REPO_FASTAS, '29611_fastas_1000', '1...1000',
                                                    '0001c4ad-0f86-3c74-9010-11fe302e765c.fasta')),
                        'No file at  ~/REPO_PDB_FASTA/fastas/29611_fastas_1000/1...1000/'
                        '0001c4ad-0f86-3c74-9010-11fe302e765c.fasta found')
        self.assertTrue(os.path.exists(os.path.join(Paths.REPO_FASTAS, '29611_fastas_1000', '2001...3000',
                                                    '048b6fff-4495-3872-bccf-3d6d77cad903.fasta')),
                        'No file at  ~/REPO_PDB_FASTA/fastas/29611_fastas_1000/2001...3000/'
                        '048b6fff-4495-3872-bccf-3d6d77cad903.fasta found')

    def test_copy_files_to_3dot_dir_2(self):
        # act
        Paths.set_up_paths(use_cluster=False)
        GUM.copy_files_to_3dot_dir(path_src_dir='/Users/u0120577/ROB_HOMOLOGY_HUMAN/RepairPDBs',
                                   path_dst_dir=Paths.REPO_PDBS, file_extension='.pdb', starting_num=1,
                                   num_to_copy_per_subdir=1000)
        # assert
        self.assertTrue(os.path.exists(os.path.join(Paths.REPO_PDBS, 'pdbs_1000', '2001...3000', 'RepairPDB_2031.pdb')),
                        'No file at  ~/REPO_PDB_FASTA/pdbs/pdbs_1000/2001...3000/RepairPDB_2031.pdb found')

    def test_write_1_fastafile_per_fasta_from_multifastafile(self):
        # arrange
        path_dst = TPLS.MC_TESTS_INPUT.value
        Paths.set_up_paths(use_cluster=False)
        path_fastafile = os.path.join(TPLS.MC_TESTS_INPUT.value, 'mutants_fastas_multifastafiles', '1...10', '1_A',
                                      'mutants', '1_A_mutants.fasta')
        # act
        GUM.write_1_fastafile_per_fasta_from_multifastafile(path_dst, path_fastafile)
