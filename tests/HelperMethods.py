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

    # path_tests        /tests/output_data where actual test results are stored or
    #                   /tests/reference_files where reference results are stored.
    # blastp_dir        Directory name of blastp data. Typically "blastp"
    # result_handle     the direct return object/handle of Biopython's qblast method
    @staticmethod
    def write_blastp_to_tests_dir(path_tests, blastp_dir, result_handle, blastp_output_xml_file):
        path_blastp = GUM.create_dir_tree(path_tests, blastp_dir)
        path_blastp_xml = path_blastp + '/' + blastp_output_xml_file
        with open(path_blastp_xml, 'w') as resx:
            resx.write(result_handle.read())
            result_handle.close()
