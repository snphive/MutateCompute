import os
from unittest import TestCase
from src.Cluster import Cluster
from tests.testHelperMethods import THM
from src.GeneralUtilityMethods import GUM


class TestCluster(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.path_zeus_FoldX_exe = '/switchlab/group/tools/FoldX_2015/FoldX'
        cls.path_local_FoldX_exe = '/Users/u0120577/SNPEFFECT/executables/FoldX'
        # paths common to both input & output
        cls.path_tests = '/Users/u0120577/PycharmProjects/MutateCompute/tests'
        cls.rel_path_PDBs = '/PDBs'
        cls.rel_path_FoldX = '/FoldX'
        # input paths only
        cls.rel_path_Inputs = '/Inputs'
        cls.path_tests_Inputs = cls.path_tests + cls.rel_path_Inputs
        cls.path_tests_Inputs_PDBs = cls.path_tests_Inputs + cls.rel_path_PDBs
        cls.rel_path_Cluster = '/Cluster'
        cls.rel_path_BuildModel = '/BuildModel'
        cls.rel_path_Fasta = '/Fasta'
        cls.rel_path_Options_Agadir = '/Options/Agadir'
        cls.rel_path_Options_FoldX = '/Options/FoldX'
        cls.rel_path_Options_Cluster = '/Options/Cluster'
        # output paths only
        cls.rel_path_Outputs = '/Outputs'
        cls.path_tests_Outputs = cls.path_tests + cls.rel_path_Outputs
        cls.rel_path_BuildModel = '/BuildModel'
        cls.rel_path_AnalyseComplex = '/AnalyseComplex'

    def setUp(self):
        self.cluster = Cluster()

    def tearDown(self):
        self.cluster = None
        THM.remove_tests_Inputs_Outputs_folders()

    # The method has named keyword argument set to '' for following: python_script_with_paths, queue, n_slots,
    # total_memory_GB, memory_limit_GB, cluster_node. None of these values is given here so the default will be applied
    # such that the resulting bash string does not include them.
    def test_write_job_q_bash(self):
        # arrange
        pdbname = 'RepairPDB_1_A'
        fxbm_jobname_prefix = 'FXBM_'
        fx_mutant_name = 'RA1A'
        jobname = fxbm_jobname_prefix + fx_mutant_name
        path_job_q_file = GUM.create_dir_tree(self.path_tests_Inputs_PDBs, pdbname, fx_mutant_name)
        if not os.path.exists(path_job_q_file):
            os.makedirs(path_job_q_file)
        expected_job_q = '#!/bin/bash\n' + '#$ -N ' + jobname + '\n' + '#$ -V\n' + '#$ -cwd\n' + \
                              'source ~/.bash_profile\n' + self.path_zeus_FoldX_exe + ' -runfile runscript.txt\n' + \
                         '\n'
        single_space = ' '
        not_expected_job_q = '#!/bin/bash' + single_space + '\n' + '#$ -N ' + jobname + '\n' + '#$ -V\n' + \
                             '#$ -cwd\n' + 'source ~/.bash_profile\n' + self.path_zeus_FoldX_exe + \
                             ' -runfile runscript.txt\n' + '\n'
        missing_new_line = ''
        not_expected_job_q_2 = '#!/bin/bash\n' + '#$ -N ' + jobname + missing_new_line + '#$ -V\n' + \
                             '#$ -cwd\n' + 'source ~/.bash_profile\n' + self.path_zeus_FoldX_exe + \
                             ' -runfile runscript.txt\n' + '\n'
        # action
        actual_job_q = self.cluster.write_job_q_bash(jobname, path_job_q_file)
        # assert
        self.assertEqual(expected_job_q, actual_job_q)
        self.assertNotEqual(not_expected_job_q, actual_job_q)
        self.assertNotEqual(not_expected_job_q_2, actual_job_q)
        self._test_job_q_bash_file_created(path_job_q_file)

    def _test_job_q_bash_file_created(self, path_job_q):
        self.assertTrue(os.path.exists(path_job_q), 'path to job.q path does not exist: ' + path_job_q)
        if os.path.exists(path_job_q):
            self.assertTrue(os.path.isfile('job.q'))

