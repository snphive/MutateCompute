from behave import *

use_step_matcher("re")


@given("pdb input files")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@when("MutateCompute_Options file specifies to perform the FoldX BuildModel operation")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@then("result should be FoldX DDG output files for all specified amino acid substitutions")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass





@given("pdb input files")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@when("MutateCompute_Options file specifies to perform the FoldX Stability operation")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@then("result should be FoldX DDG output files")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass




@given("pdb input files of separate proteins and/or protein chains that form a complex")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@when("MutateCompute_Options file specifies to perform the FoldX AnalyseComplex operation")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@step("different proteins of a protein complex are included")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass


@then("result should be FoldX DDG output files for all complexes for which there is a structure")
def step_impl(context):
    """
    :type context: behave.runner.Context
    """
    pass

