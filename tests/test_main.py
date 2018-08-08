from unittest import TestCase
from unittest.mock import patch
from src.Main import Main
from tests.TestPathsAndListsSeqs import TPLS
from src.GeneralUtilityMethods import GUM
from src.Scheduler import Scheduler


class TestMain(TestCase):

    @classmethod
    def setUpClass(cls):
        GUM.linux_copy(path_src=TPLS.CONFIG_FOR_READ_ONLY.value, path_dst=TPLS.MC_TESTS.value, do_recursively=True)

    def test__read_global_options(self):
        # arrange
        path_globoptions_file = TPLS.MC_TESTS_CONFIG_GLOBAL_OPTIONS.value + '/MutateCompute_Options.txt'
        expected_global_options = ['######## NUMBER OF PDBs TO ANALYSE ##############################################'
                                   '#################################\n', '#\n', '## All means all PDBs in the '
                                    'specified source directory\n', '#\n', 'PDBs: 4;\n', '#\n', '######## NUMBER OF '
                                    'FASTA FILES TO ANALYSE #########################################################'
                                    '###############\n', '#\n', '## All means all FASTA files in specified source '
                                    'directory\n', '#\n', 'FASTAs: 10;\n', '#\n', '######## OPERATIONS YOU WANT TO '
                                    'RUN ############################################################################'
                                    '##\n', '#\n', '## MUTATE_FASTA = generate FASTA files of every possible point '
                                    'mutant for a given FASTA file.\n', '## AGADIR = run Agadirwrapper. The Options.txt'
                                    ' file specifies which of Agadirwrapper algorithms to perform.\n', '## FOLDX_REPAIR'
                                    ' = repair pdb with FoldX repair algorithm\n', '## FOLDX_BUILDMODEL = mutate '
                                    'sequence, remodel and determine DDG of mutated pdb with FoldX BuildModel '
                                    'algorithm\n', '## FOLDX_STABILITY = determine DDG of pdb with FoldX Stability '
                                    'algorithm\n', '## FOLDX_ANALYSECOMPLEX = determine DDG of pdb with FoldX '
                                    'Stability algorithm\n', '#\n', 'MUTATE_FASTA: TRUE\n', 'AGADIR: FALSE\n',
                                   'FOLDX_REPAIR: FALSE\n', 'FOLDX_BUILDMODEL: FALSE\n', 'FOLDX_STABILITY: FALSE\n',
                                   'FOLDX_ANALYSECOMPLEX: FALSE\n', '#\n', '######## AMINO ACIDS YOU WANT TO MUTATE '
                                   'TO ########################################################################\n',
                                   '#\n', '## All means all 20 amino acids. Otherwise each amino acid in FASTA format,'
                                    ' e.g. "ACDEFGHIKLMNOP;"\n', '#\n', 'RESIDUES: All;']
        # action
        global_options = Main._read_global_options(path_globoptions_file)
        # assert
        self.assertEqual(expected_global_options, global_options)

    def test__build_filelist_for_analysis(self):
        # arrange
        globaloptions_lines = ['#\n', 'PDBs: 4;\n', '#']
        PDB_or_FASTA = 'PDBs'
        path_input = TPLS.MC_TESTS_INPUT.value
        expected_pdb_list = ['RepairPDB_1.pdb', 'RepairPDB_2.pdb', 'RepairPDB_3.pdb', 'RepairPDB_4.pdb']
        # action
        pdb_list = Main._build_filelist_for_analysis(globaloptions_lines, PDB_or_FASTA, path_input)
        # assert
        self.assertNotEqual('', pdb_list)
        self.assertListEqual(expected_pdb_list, pdb_list)

        # arrange
        globaloptions_lines = ['#\n', 'FASTAs: 4;\n', '#']
        PDB_or_FASTA = 'FASTAs'
        expected_fasta_list = ['1_A.fasta', '1_B.fasta', '2_A.fasta', '3_A.fasta']
        # action
        fasta_list = Main._build_filelist_for_analysis(globaloptions_lines, PDB_or_FASTA, path_input)
        # assert
        self.assertNotEqual('', fasta_list)
        self.assertListEqual(expected_fasta_list, fasta_list)

    def test__determine_which_operations_to_perform(self):
        # arrange
        globaloptions_lines = ['## FOLDX_ANALYSECOMPLEX = determine DDG of pdb with FoldX Stability algorithm\n',
                               '#\n', 'MUTATE_FASTA: TRUE\n', 'AGADIR: FALSE\n', 'FOLDX_REPAIR: FALSE\n',
                               'FOLDX_BUILDMODEL: FALSE\n', 'FOLDX_STABILITY: FALSE\n', 'FOLDX_ANALYSECOMPLEX: FALSE\n',
                               '#\n']
        expected_operations = {'do_mutate_fasta': True, 'do_agadir': False, 'do_foldx_repair': False,
                               'do_foldx_buildmodel': False, 'do_foldx_stability': False,
                               'do_foldx_analysecomplex': False}
        # action
        operations = Main._determine_which_operations_to_perform(globaloptions_lines)
        # assert
        self.assertNotEqual({}, operations)
        self.assertDictEqual(expected_operations, operations)

    def test__determine_residues_to_mutate_to(self):
        # arrange
        globaloptions_lines = ['#\n', 'RESIDUES: All;\n']
        expected_residues = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                             'W', 'Y']
        # action
        residues = Main._determine_residues_to_mutate_to(globaloptions_lines)
        # assert
        self.assertListEqual(expected_residues, residues)

        # arrange
        globaloptions_lines = ['#\n', 'RESIDUES: ACDEFGHI;\n']
        expected_residues = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
        # action
        residues = Main._determine_residues_to_mutate_to(globaloptions_lines)
        # assert
        self.assertListEqual(expected_residues, residues)

        # arrange
        globaloptions_lines = ['#\n', 'RESIDUES: AXZ;\n']
        expected_residues = ['A']
        # action
        residues = Main._determine_residues_to_mutate_to(globaloptions_lines)
        # assert
        self.assertListEqual(expected_residues, residues)



    @patch.object(Scheduler, 'start')
    def test__start_scheduler(self, mock_start):
        # arrange
        operations = {'do_mutate_fasta': True}
        # operations = {} would need to add something to here to check for raising a ValueError (mock sideeffect?)
        path_input = TPLS.MC_TESTS_INPUT
        pdb_list = ['RepairPDB_1.pdb']
        fasta_list = ['1_A.fasta']
        list_of_mutant_aa = ['A']
        mock_start.return_value = True
        expected = True
        # action
        actual = Main._start_scheduler(operations, path_input, pdb_list, fasta_list, list_of_mutant_aa)
        # assert
        self.assertEqual(expected, actual)
