from unittest import TestCase
from src.FoldX import FoldX
import threading as Thread
from src.MutateFasta import MutateFasta


class TestScheduler(TestCase):


    def test_start_multithreading(self):
        Thread.activeCount()
