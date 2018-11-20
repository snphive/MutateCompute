import os
from unittest import TestCase
from unittest.mock import patch
from src.FoldX import FoldX
from src.enums.Conditions import Cond
from src.tools.GeneralUtilityMethods import GUM
from src.Cluster import Cluster
from src.enums.Paths import Paths
from tests.TestPathsAndListsSeqs import TPLS


# Note: A test pdb is used here which is taken from the RepairPDBs folder, but includes only the first 5 residues of
# to facilitate manual creation of test variables, such as the list of all mutants.
# Note: Wherever a path is indicated, it is the absolute path. Relative paths are always prefixed with rel_
class TestBuildModel(TestCase):

    def setUp(self):
        self.foldx = FoldX()
        self.buildModel = self.foldx.BuildModel(Cond.INCELL_MAML_FX.value)

    def tearDown(self):
        self.buildModel = None

    @patch.object(FoldX.BuildModel, '_make_fx_mutant_name_list')
    @patch.object(GUM, 'extract_pdbname_chain_startpos_fasta_from_pdbs')
    def test_mutate_protein_structure_Not_Cluster(self, mock_extract_pdbname_chain_startpos_fasta_from_pdbs,
                                                  mock__make_fx_mutant_name_list):
        # arrange
        Paths.set_up_paths(use_cluster=False)
        mock_extract_pdbname_chain_startpos_fasta_from_pdbs.return_value = {'RepairPDB_1A': 'RVYLT',
                                                                            'RepairPDB_1B': 'RVYLT'}
        # Note FoldX expects mutant names to have format wtaa_chain_position_mutant aa.
        mock__make_fx_mutant_name_list.return_value = ['RA1A', 'RA1C', 'VA2A', 'VA2C', 'YA3A', 'YA3C', 'LA4A', 'LA4C',
                                                       'TA5A', 'TA5C', 'RB540A', 'RB540C', 'VB541A', 'VB541C',
                                                       'YB542A', 'YB542C', 'LB543A', 'LB543C', 'TB544A', 'TB544C']
        amino_acids = ['A', 'C']
        path_pdb = '/Users/u0120577/PycharmProjects/MutateCompute/tests/input_data/pdbs/RepairPDB_1_first5aa.pdb'
        # act
        self.buildModel.mutate_protein_structure(path_pdbfile=path_pdb, amino_acids=amino_acids)
        # assert
        expected_call_count_1 = 1
        self.assertEqual(expected_call_count_1, mock_extract_pdbname_chain_startpos_fasta_from_pdbs.call_count)
        self.assertEqual(expected_call_count_1, mock__make_fx_mutant_name_list.call_count)

    @patch('subprocess.call')
    @patch.object(Cluster, 'write_job_q_bash')
    @patch.object(FoldX.BuildModel, '_write_individual_list_for_mutant')
    @patch('os.chdir')
    @patch.object(GUM, '_os_makedirs')
    @patch.object(FoldX, 'write_runscript_file')
    @patch.object(FoldX.BuildModel, '_make_fx_mutant_name_list')
    @patch.object(GUM, 'extract_pdbname_chain_startpos_fasta_from_pdbs')
    def test_mutate_protein_structure_Cluster(self, mock_extract_pdbname_chain_startpos_fasta_from_pdbs,
                                              mock__make_fx_mutant_name_list, mock_write_runscript_file,
                                              mock__os_makedirs, mock_os_chdir, mock__write_individual_list_for_mutant,
                                              mock_write_job_q_bash, mock_subprocess_call):
        # arrange
        Paths.set_up_paths(use_cluster=True)
        mock_extract_pdbname_chain_startpos_fasta_from_pdbs.return_value = {'RepPDB_1A_test': 'RVYLT', 'RepPDB_1B_test': 'RVYLT'}
        # Note FoldX expects mutant names to have format wtaa_chain_position_mutant aa.
        mock__make_fx_mutant_name_list.return_value = ['RA1A', 'RA1C', 'VA2A', 'VA2C', 'YA3A', 'YA3C', 'LA4A', 'LA4C',
                                                       'TA5A', 'TA5C', 'RB540A', 'RB540C', 'VB541A', 'VB541C',
                                                       'YB542A', 'YB542C', 'LB543A', 'LB543C', 'TB544A', 'TB544C']
        mock_os_chdir.return_value = None
        write_wt_fastafiles = False
        amino_acids = ['A', 'C']
        path_pdb = '/Users/u0120577/PycharmProjects/input_data/pdbs/1...10/RepairPDB_1.pdb'
        # act
        self.buildModel.mutate_protein_structure(path_pdbfile=path_pdb, amino_acids=amino_acids)
        # assert
        expected_call_count_1 = 1
        expected_call_count_20 = len(mock__make_fx_mutant_name_list.return_value)
        expected_call_count_22 = 22
        self.assertEqual(expected_call_count_1, mock_extract_pdbname_chain_startpos_fasta_from_pdbs.call_count)
        self.assertEqual(expected_call_count_1, mock__make_fx_mutant_name_list.call_count)
        self.assertEqual(expected_call_count_1, mock_write_runscript_file.call_count)
        self.assertEqual(expected_call_count_22, mock__os_makedirs.call_count)
        self.assertEqual(expected_call_count_20, mock__write_individual_list_for_mutant.call_count)
        self.assertEqual(expected_call_count_20, mock_write_job_q_bash.call_count)
        self.assertEqual(expected_call_count_20, mock_subprocess_call.call_count)

    # Note FoldX expects mutant names to have format wtaa_chain_position_mutantaa.
    # Variable names for values with this format are denoted with '_fx_'
    def test_make_fx_mutant_name_list(self):
        pdbname_chain_startpos_fasta_dict = {'Test1Residue_A': 'R', 'Test1Residue_B': 'R'}
        amino_acids = ['A', 'C', 'D']
        fx_mutant_name_list = self.buildModel._make_fx_mutant_name_list(amino_acids, pdbname_chain_startpos_fasta_dict)
        expected_fx_mutant_name_list = ['RA1A', 'RA1C', 'RA1D', 'RB1A', 'RB1C', 'RB1D']
        self.assertEqual(fx_mutant_name_list, expected_fx_mutant_name_list)

    # def test_write_job_q_filenames_check(self):
        # job_q_file = open('./job.q', 'w') - how to assert filenames are same?
        # job_q_file.write('#!/bin/bash\n') - not sure what this  does

    # file comparison scripts, currently unused. It relies to some extent on the directory structure being in place
    # which I didn't know how to mock so I made the write_job.. method return the string that it writes to job.q file.
    # with open(os.path.join(path_job_q_file, 'job.q')) as f:
    #   actual_job_q = f.read()
    #   self.assertEqual(actual_job_q, expected_job_q)
    #   self.assertNotEqual(actual_job_q, not_expected_job_q)
    #   self.assertNotEqual(actual_job_q, not_expected_job_q_2)

    # pdbname = pdb.split('.')[0]
    # path_tests_Inputs_PDBs_pdbname = os.path.join(PFT.PATH_TESTS_INPTS_PDBS, pdbname)
    # path_tests_Inputs_PDBs_pdbname_FX_BM = os.path.join(path_tests_Inputs_PDBs_pdbname, FX_BuildModel')
    # path_tests_Inputs_PDBs_pdbname_FXmutant = os.path.join(path_tests_Inputs_PDBs_pdbname, 'RA1A')
    # mock_os_path_join.side_effect = [path_tests_Inputs_PDBs_pdbname, path_tests_Inputs_PDBs_pdbname_FX_BM,
    #                                  path_tests_Inputs_PDBs_pdbname_FXmutant]

    # Write_runscript_for_pdbs() takes 6 arguments. The last 3 (namely show_sequence_detail, print_networks,
    # calculate_stability) are keyword (named) arguments. All have default values assigned in the method argument so
    # that when the method is called without supplying a (keyword-named) value, the default value is applied.
    # (The default values for all 3 are False).
    #
    # 30.07.18 Have redesigned the directory structure such that runscripts will go in /configuration /foldx_config
    # /ac_runscript or /bm_runscript (and others).
    def test_write_runscript_file(self):
        # arrange
        pdb = 'RepairPDB_1.pdb'
        path_runscript = TPLS.MC_TESTS_CONFIG_FXCONFIG_BMRUNSCRIPT.value
        try:
            os.makedirs(path_runscript)
        except FileExistsError:
            print('Part/all dir tree already exists. No probs.')
        action = '<BuildModel>#,individual_list.txt'
        expected_runscript = '<TITLE>FOLDX_runscript;\n' + '<JOBSTART>#;\n' + '<PDBS>' + pdb + ';\n' + \
                             '<BATCH>#;\n' + '<COMMANDS>FOLDX_commandfile;\n' + action + ';\n' + '<END>#;\n' + \
                             '<OPTIONS>FOLDX_optionfile;\n' + '<Temperature>298;\n' + '<IonStrength>0.05;\n' + \
                             '<ph>7;\n' + '<moveNeighbours>true;\n' + '<VdWDesign>2;\n' + '<numberOfRuns>3;\n' + \
                             '<OutPDB>#;\n' + '<END>#;\n' + '<JOBEND>#;\n' + '<ENDFILE>#;\n'
        single_space = ' '
        not_expected_runscript = single_space + expected_runscript
        # act
        actual_runscript = self.foldx.write_runscript_file(
            path_runscript=TPLS.MC_TESTS_CONFIG_FXCONFIG_BMRUNSCRIPT.value, pdb=pdb,
            conditions=Cond.INCELL_MAML_FX.value, action=action)
        # assert
        self.assertEqual(actual_runscript, expected_runscript)
        self.assertNotEqual(actual_runscript, not_expected_runscript)

    def test_write_individual_list_for_mutant(self):
        # arrange
        fx_mutant_name = 'RA1A'
        path_file_dest_dir = TPLS.MC_TESTS_CONFIG_FXCONFIG_BMRUNSCRIPT.value
        # act
        self.buildmodel._write_individual_list_for_mutant(fx_mutant_name, path_file_dest_dir)
        # assert

    def test_make_fx_mutant_name_list(self):
        # arrange
        amino_acids = TPLS.LIST_ALL_20_AA.value
        pdbname_chain_startpos_fasta = 'PDB1_X_1: TESTSEQ'
        # act
        self.buildModel._make_fx_mutant_name_list(amino_acids, pdbname_chain_startpos_fasta)
        # assert

    def test_make_output_dir_and_copy_fxconfig_files_in(self):
         # arrange
         path_output = TPLS.MC_TESTS_OUTPUT.value
         pdbname = 'test_pdb_dir'
         fx_mutant_name = 'wtaa_chain_pos_mut_name'
         # act
         self.buildModel._make_output_dir_and_copy_fxconfig_files_in(path_output, pdbname, fx_mutant_name)
         # assert
