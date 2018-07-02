from behave import *
import yaml
import src.Main
from src.MutateFasta import MutateFasta

with open("/switchlab/group/shazib/SnpEffect/SourceFiles/Scripts/pathsAndDictionaries.yaml",
          'r') as stream:
    try:

        paths_and_dictionaries = yaml.load(stream)
        path_R_exe = paths_and_dictionaries['ROOT']['path_R_exe']
        path_FoldX_exe = paths_and_dictionaries['ROOT']['path_FoldX_exe']
        path_Agadir_exe = paths_and_dictionaries['ROOT']['path_Agadir_exe']
        path_QSub_exe = paths_and_dictionaries['ROOT']['path_QSub_exe']

        path_SnpEffect_dir = paths_and_dictionaries['ROOT']['path_SnpEffect_dir']
        path_SE_SourceFiles_Scripts_dir = paths_and_dictionaries['ROOT'][' path_SourceFiles_Scripts_dir']
        path_SE_Inputs_PDBs_dir = paths_and_dictionaries['ROOT']['path_SE_Inputs_PDBs_dir']
        path_SE_Inputs_FASTAs_dir = paths_and_dictionaries['ROOT']['path_SE_Inputs_FASTAs_dir']
        path_SE_Outputs_dir = paths_and_dictionaries['ROOT']['path_SE_Outputs_dir']
        path_SE_Outputs_Agadir_dir = paths_and_dictionaries['ROOT']['path_SE_Outputs_Agadir_dir']
        path_SE_Outputs_FoldX_dir = paths_and_dictionaries['ROOT']['path_SE_Outputs_FoldX_dir']

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