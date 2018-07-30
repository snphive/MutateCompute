from behave import *
import yaml
import src.Main
from src.MutateFasta import MutateFasta

with open("/switchlab/group/shazib/SnpEffect/SourceFiles/Scripts/test_pathsAndDictionaries.yaml",
          'r') as stream:
    try:

        paths_and_dictionaries = yaml.load(stream)
        path_R_exe = paths_and_dictionaries['ROOT']['path_R_exe']
        path_FoldX_exe = paths_and_dictionaries['ROOT']['path_FoldX_exe']
        path_Agadir_exe = paths_and_dictionaries['ROOT']['path_Agadir_exe']
        path_QSub_exe = paths_and_dictionaries['ROOT']['path_QSub_exe']

        path_SnpEffect = paths_and_dictionaries['ROOT']['path_SnpEffect']
        path_SE_SourceFiles_Scripts = paths_and_dictionaries['ROOT'][' path_SourceFiles_Scripts']
        path_SE_Inputs_PDBs = paths_and_dictionaries['ROOT']['path_SE_Inputs_PDBs']
        path_SE_Inputs_FASTAs = paths_and_dictionaries['ROOT']['path_SE_Inputs_FASTAs']
        path_SE_Outputs = paths_and_dictionaries['ROOT']['path_SE_Outputs']
        path_SE_Outputs_Agadir = paths_and_dictionaries['ROOT']['path_SE_Outputs_Agadir']
        path_SE_Outputs_FoldX = paths_and_dictionaries['ROOT']['path_SE_Outputs_FoldX']

        dict_aa_1to3 = paths_and_dictionaries['ROOT']['dict_aa_1to3']
        dict_aa_3to1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']
        list_all_20_aa = paths_and_dictionaries['ROOT']['list_all_20_aa']

    except yaml.YAMLError as exc:
        print(exc)


@given('{fasta} input files')
def step_impl(context, fasta):
    mutate_all_residues = True
    context.mutateFasta = MutateFasta(fasta, mutate_all_residues, list_all_20_aa)
    pass


@when("MutateCompute_Options file specifies to perform the Mutate operation")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@then("result should be FASTA files for every possible amino acid substitution")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@given("FASTA input files")
def step_impl(context):
    context.main = src.Main
    context.mutateFasta = MutateFasta()
    pass


@when("MutateCompute_Options file specifies to perform the Agadir operation")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@then("result should be Agadir PSX output files")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass