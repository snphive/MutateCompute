from unittest import TestCase
from src.FoldX import FoldX
import threading as Thread
from src.MutateFasta import MutateFasta
from src.Scheduler import Scheduler
from tests.TestPathsAndListsSeqs import TPLS

class TestScheduler(TestCase):


    def test_start_multithreading(self):
        Thread.activeCount()


    # @patch.object(Scheduler, 'start')
    def test__start_scheduler(self, mock_start):
        mock_start.return_value = True
        expected = True
        # act
        actual = Scheduler.start(operations={'do_mutate_fasta': True}, path_input=TPLS.MC_TESTS_INPUT.value,
                                 path_output=TPLS.MC_TESTS_OUTPUT, pdb_list=['RepairPDB_1.pdb'],
                                 fastafile_list=['1_A.fasta'], amino_acids=['A'], use_multithread=False,
                                 write_1_fasta_only=True, write_fasta_per_mut=True)
        # assert
        self.assertEqual(expected, actual)