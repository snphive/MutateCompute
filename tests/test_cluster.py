import os
from unittest import TestCase
from src.Cluster import Cluster
from src.tools.GeneralUtilityMethods import GUM
from tests.TestPathsAndListsSeqs import TPLS


class TestCluster(TestCase):

    # Currently the tests are copying over all configuration and input data from the main directory into the tests
    # before running the tests (i.e. here in the setUpClass method).
    # The data in those main folders will be programmatically generated but is currently manually transferred.
    @classmethod
    def setUpClass(cls):
        if not os.path.exists(TPLS.MC_TESTS_CONFIG.value):
            GUM.linux_copy_all_files_in_dir(path_src_dir=TPLS.CONFIG_FOR_READ_ONLY.value,
                                            path_dst_dir=TPLS.MC_TESTS_CONFIG.value)

        if not os.path.exists(TPLS.MC_TESTS_INPUT.value):
            GUM.linux_copy_all_files_in_dir(path_src_dir=TPLS.INPUT_FOR_READ_ONLY.value,
                                            path_dst_dir=TPLS.MC_TESTS_INPUT.value)

    # @classmethod
    # def tearDownClass(cls):
    #     # HM.remove_config_folders()

    def setUp(self):
        self.cluster = Cluster()

    def tearDown(self):
        self.cluster = None

    # The method has named keyword argument set to '' for following: python_script_with_paths, queue, n_slots,
    # total_memory_GB, memory_limit_GB, cluster_node. None of these values is given here so the default will be applied
    # such that the resulting bash string does not include them.
    def test_write_job_q_bash(self):
        # arrange
        pdbname = 'RepairPDB_1_A'
        fxbm_jobname_prefix = 'FXBM_'
        fx_mutant_name = 'RA1A'
        jobname = fxbm_jobname_prefix + fx_mutant_name
        expected_job_q = '#!/bin/bash\n' + '#$ -N ' + jobname + '\n' + '#$ -V\n' + '#$ -cwd\n' + \
                              'source ~/.bash_profile\n' + TPLS.ZEUS_FOLDX_EXE.value + ' -runfile runscript.txt\n'
        # not_expected_job_q is just the same as expected_job_q but has an extra single single
        single_space = ' '
        not_expected_job_q = '#!/bin/bash' + single_space + '\n' + '#$ -N ' + jobname + '\n' + '#$ -V\n' + \
                             '#$ -cwd\n' + 'source ~/.bash_profile\n' + TPLS.ZEUS_FOLDX_EXE.value + \
                             ' -runfile runscript.txt\n'
        # not_expected_job_q is just the same as expected_job_q but is missing a \n
        missing_new_line = ''
        not_expected_job_q_2 = '#!/bin/bash\n' + '#$ -N ' + jobname + missing_new_line + '#$ -V\n' + \
                             '#$ -cwd\n' + 'source ~/.bash_profile\n' + TPLS.ZEUS_FOLDX_EXE.value + \
                             ' -runfile runscript.txt\n'
        # act
        actual_job_q = self.cluster.write_job_q_bash(jobname, TPLS.MC_TESTS_CONFIG_JOBQ.value, using_runscript=True)
        # assert
        self.assertEqual(expected_job_q, actual_job_q)
        self.assertNotEqual(not_expected_job_q, actual_job_q)
        self.assertNotEqual(not_expected_job_q_2, actual_job_q)
        self._test_job_q_bash_file_created(os.path.join(TPLS.MC_TESTS_CONFIG_JOBQ.value, 'job.q'))

    def _test_job_q_bash_file_created(self, path_job_q_file):
        path_job_q_dir = path_job_q_file.strip(path_job_q_file.split('/')[-1])
        self.assertTrue(os.path.exists(path_job_q_dir), 'path to job.q directory (which should contain the job.q file) '
                                                        'does not exist: ' + path_job_q_dir)
        if os.path.exists(path_job_q_dir):
            self.assertTrue(os.path.exists(path_job_q_file), '(Absolute path to) job.q file does not exist: '
                            + path_job_q_file)
            self.assertTrue(os.path.isfile(path_job_q_file))

