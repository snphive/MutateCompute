import os
import shutil
from src.GeneralUtilityMethods import GUM
from tests.TestPaths import TPaths


# A collection of "helper" methods to be used by unit tests.
class HM(object):

    @staticmethod
    def remove_tests_input_output_dirs():
        if os.path.exists(TPaths.MC_TESTS.value):
            if os.path.exists(TPaths.MC_TESTS_INPUT.value):
                HM.__delete_subdirectory_tree_of_tests(TPaths.MC_TESTS_INPUT.value)
            if os.path.exists(TPaths.MC_TESTS_OUTPUT.value):
                HM.__delete_subdirectory_tree_of_tests(TPaths.MC_TESTS_OUTPUT.value)

    @staticmethod
    def remove_config_folders():
        if os.path.exists(TPaths.MC_TESTS.value):
            if os.path.exists(TPaths.MC_TESTS_CONFIG.value):
                HM.__delete_subdirectory_tree_of_tests(TPaths.REL_CONFIG.value)

    @staticmethod
    def __delete_subdirectory_tree_of_tests(path_to_delete):
        os.chdir(TPaths.MC_TESTS.value)
        if not os.getcwd() == TPaths.MC_TESTS.value:
            raise ValueError('Current working directory is not MutateCompute/tests. '
                             '\nNot proceeding with deletion of ' + path_to_delete)
        else:
            try:
                shutil.rmtree('/' + path_to_delete)
            except OSError as e:
                print("Error removing: %s - %s." % (e.filename, e.strerror))

    # path_tests        tests/Outputs where actual test results are stored
    #                   tests/Files where reference results are stored
    # result_handle     the direct return object/handle of Biopython's qblast method
    @staticmethod
    def write_blast_run_to_tests_dir(path_tests, blast_dir, result_handle, blast_output_xml_file):
        path_blast_run = GUM.create_dir_tree(path_tests, blast_dir)
        path_blastrun_xml = path_blast_run + '/' + blast_output_xml_file
        with open(path_blastrun_xml, 'w') as resx:
            resx.write(result_handle.read())
            result_handle.close()
