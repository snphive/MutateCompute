import os
import shutil


class THM(object):

    # paths common to both input & output
    path_tests = '/Users/u0120577/PycharmProjects/MutateCompute/tests'
    rel_path_PDBs = '/PDBs'
    rel_path_FoldX = '/FoldX'
    # input paths only
    rel_path_Inputs = '/Inputs'
    path_tests_Inputs = path_tests + rel_path_Inputs
    path_tests_Inputs_PDBs = path_tests_Inputs + rel_path_PDBs
    rel_path_Cluster = '/Cluster'
    rel_path_BuildModel = '/BuildModel'
    rel_path_Fasta = '/Fasta'
    rel_path_Options_Agadir = '/Options/Agadir'
    rel_path_Options_FoldX = '/Options/FoldX'
    rel_path_Options_Cluster = '/Options/Cluster'
    # output paths only
    rel_path_Outputs = '/Outputs'
    path_tests_Outputs = path_tests + rel_path_Outputs
    rel_path_BuildModel = '/BuildModel'
    rel_path_AnalyseComplex = '/AnalyseComplex'

    @staticmethod
    def remove_tests_Inputs_Outputs_folders():
        if os.path.exists(THM.path_tests):
            if os.path.exists(THM.path_tests_Inputs):
                THM.__delete_directory_tree_of_tests_InputsOutputs(THM.rel_path_Inputs.strip('/'))
            if os.path.exists(THM.path_tests_Outputs):
                THM.__delete_directory_tree_of_tests_InputsOutputs(THM.rel_path_Outputs.strip('/'))

    @staticmethod
    def __delete_directory_tree_of_tests_InputsOutputs(rel_path_InputsOutputs):
        os.chdir(THM.path_tests)
        print(os.getcwd())
        if not os.getcwd() == '/Users/u0120577/PycharmProjects/MutateCompute/tests':
            raise ValueError('Current working directory is not MutateCompute/tests. Not proceeding with deletion of '
                             'inputs & outputs subfolders')
        else:
            try:
                shutil.rmtree(rel_path_InputsOutputs)
            except OSError as e:
                print("Error removing: %s - %s." % (e.filename, e.strerror))

