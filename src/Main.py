#!/usr/bin/env python3
"""
Class for initiating all programs except blast.

Calling the constuctor can be done in the IDE through the short script called KickOffAuto.py.
Alternatively it can be run from the command line via KickOff.sh, which does the following 3 things:
1. Sources .profile script that sets the python environment variable to this MutateCompute project folder.
2. Sets path to the virtual environment folder which includes all the downloaded libraries used here.
3. Executes KickOffAuto.py via python3 command to indicate that version 3 and not 2 should be used here.
"""

import glob
import natsort
from src.Scheduler import Scheduler
from src.enums.Str import Str
from src.enums.Paths import Paths
from src.enums.AminoAcids import AA

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class Main(object):

    def __init__(self, operations: dict, use_multithread: bool, path_input: str, path_output: str, path_pdbfiles: list,
                 path_fastafiles: list, specific_fxmutants: list, amino_acids: list):
        """
        All algorithmic analyses start from here. IdentifyProtein.py which runs Blast analyses is the only program that
        currently runs independently.
        :param operations: Operations flagged True to run, False to not run.
        :param use_multithread: True if program(s) should be spawned on separate threads/processes.
        :param path_input: Abs path to input root dir
        :param path_output: Abs path to output root dir
        :param path_pdbfiles: Abs path to repo of pdbfiles.
        :param path_fastafiles: Abs path to repo of fastafiles.
        :param amino_acids: All amino acids that mutation operations will mutate residues to.
        """
        if operations == {} or Main._all_ops_are_false(operations):
            raise ValueError("All options in 'operations' were either set to FALSE or some typographical error. "
                             "Check /config/global_options/global_options.txt was written correctly")
        elif operations[Str.OPER_RUN_MUT_FSTA.value] or operations[Str.OPER_RUN_AGDR.value] or \
                operations[Str.OPER_RUN_FX_RPR.value] or operations[Str.OPER_RUN_FX_BM.value] or \
                operations[Str.OPER_RUN_FX_STAB.value] or operations[Str.OPER_RUN_FX_AC.value]:
            Scheduler.start(operations, use_multithread, path_input, path_output, path_pdbfiles, path_fastafiles,
                            specific_fxmutants, amino_acids, write_1_fasta_only=True, write_fasta_per_mut=False)

    @staticmethod
    def _all_ops_are_false(operations: dict):
        """
        :param operations:
        :return: True indicates all oprations are set to False in global_options.txt.
        """
        if not operations[Str.OPER_RUN_MUT_FSTA.value] and not operations[Str.OPER_RUN_AGDR.value] \
                and not operations[Str.OPER_RUN_FX_RPR.value] and not operations[Str.OPER_RUN_FX_BM.value] \
                and not operations[Str.OPER_RUN_FX_STAB.value] and not operations[Str.OPER_RUN_FX_AC.value]:
            return True
        else:
            return False

    @staticmethod
    def _read_global_options(path_globaloptions_file: str):
        """
        Reads the global_options.txt in /config/global_options breaking the text up according to newlines.
        :param path_globaloptions_file: Abs path to global_options file (incl. extension).
        :return: Text of the file as a list of strings for each line ending with \n.
        """
        with open(path_globaloptions_file, 'r') as globaloptions:
            globaloptions_lines = globaloptions.readlines()
        return globaloptions_lines

    # Note the risk of inefficiency where the program spends a lot of time sorting the source directory to build the
    # list of available target files, but only a relatively small number of files will be used E.g. source dir sorts
    # about 30,000 pdbs for user to get access to first 5 only!
    @staticmethod
    def _build_filelist_for_analysis(globaloptions_lines: list, path_repo_pdbs_or_fastas: str):
        """
        Takes the /config/global_options/global_options.txt text as a list of \n separated lines.
        Reads the option given at either "PDBs" and "FASTAs".
        The option is expected to be in one of four forms: "", "all", "<filename(s)>", or "<number>".
        '' indictes no files - returns None.
        'all' indicates all files in specified input source directory.
        '<filename(s)>' indicates specific files by name (if found in the specified input source directory).
        '<number>' indicates the maximum number of files from the naturally-sorted list that is based on what files are
        present in the specified input_data source directory.
        E.g. line = "PDBs: 4;\n", returns ['RepairPDB_1.pdb', 'RepairPDB_2.pdb', 'RepairPDB_3.pdb', 'RepairPDB_4.pdb'].
        :param globaloptions_lines: Alphanumeric text lines of the global options ending with "\n".
        :param path_repo_pdbs_or_fastas: Path to the repository files (including pdbfiles and fastafiles).
        :return: list of pdbfiles or fastafiles according to global_options.txt. The list can only include files that
        are actually in the path_repo (typically referred to here as a repository dir)
        """
        file_list = []
        file_extension = Str.PDBEXT.value
        pdbs_or_fastas_option = Main.Strs.OPT_PDBs.value
        if Paths.DIR_FASTAS.value in path_repo_pdbs_or_fastas:
            file_extension = Str.FSTAEXT.value
            pdbs_or_fastas_option = Main.Strs.OPT_FASTAs.value
        path_repo_files = path_repo_pdbs_or_fastas + Str.FSLSH_ASTRX.value + file_extension
        for line in globaloptions_lines:
            if '#' in line:
                continue
            if pdbs_or_fastas_option in line:
                # Main.__validate_option_text(line)
                pdb_or_fasta_option = Main.__get_text_after_colon_before_semi(line)
                if pdb_or_fasta_option == '' or pdb_or_fasta_option == '0':
                    file_list = None
                    break
                elif pdb_or_fasta_option.isdigit() and not pdb_or_fasta_option.isalpha():
                    # check that the number is reasonable, i.e. less than 1000 ?
                    num_of_files_to_analyse = int(pdb_or_fasta_option)
                    path_file_natsorted_list = natsort.natsorted(glob.glob(path_repo_files))
                    if not path_file_natsorted_list:
                        raise ValueError('No files were found in the input_data directory tree. Files must be manually '
                                         'placed there')
                    for i, path_file_natsorted in enumerate(path_file_natsorted_list):
                        if i == num_of_files_to_analyse:
                            break
                        file_list.append(path_file_natsorted.split('/')[-1])
                elif pdb_or_fasta_option.isalpha():
                    if pdb_or_fasta_option.lower() == Main.Strs.OPT_ALL.value:
                        path_file_list = glob.glob(path_repo_files)
                        for path_file in path_file_list:
                            file_list.append(path_file.split('/')[-1])
                    # else:
                        # if there is a comma-separated list of pdb or fasta names
                else:
                    raise ValueError('The option for pdbs or fastas is not All, a number, a list of sequence names or '
                                     'blank. Something is wrong with the PDBs and/or FASTAs global_options option')
                break
        return file_list

    @staticmethod
    def _determine_which_operations_to_perform(globaloptions_lines: list):
        """
        Determines which operations to perform, according to the global options file - True or False indicates to
        perform this operation or not, respectively. Format expected is "True" or "False" comes after ":" and before
        ";\n"
        :param globaloptions_lines: Alphanumeric text lines of the global options ending with "\n".
        :return: dict with operation str as key, True/False boolean as value.
        """
        operations = {}
        for line in globaloptions_lines:
            if "#" in line:
                continue
            if Main.Strs.OPT_MUT_FSTA.value in line:
                operations[Str.OPER_RUN_MUT_FSTA.value] = Main.__is_true(line)
            if Main.Strs.OPT_AGDR.value in line:
                operations[Str.OPER_RUN_AGDR.value] = Main.__is_true(line)
            if Main.Strs.OPT_FX_RPR.value in line:
                operations[Str.OPER_RUN_FX_RPR.value] = Main.__is_true(line)
            if Main.Strs.OPT_FX_BM.value in line:
                operations[Str.OPER_RUN_FX_BM.value] = Main.__is_true(line)
            if Str.OPER_RUN_FX_STAB.value in line:
                operations[Str.OPER_RUN_FX_STAB.value] = Main.__is_true(line)
            if Str.OPER_RUN_FX_AC.value in line:
                operations[Str.OPER_RUN_FX_AC.value] = Main.__is_true(line)
        return operations

    @staticmethod
    def _determine_residues_to_mutate_to(globaloptions_lines: list):
        """
        Determines which residues to mutate to according to the global_options file. Given either a list of residues
        in FASTA) format, or the word "All" to indicate all 20 residues.
        :param globaloptions_lines: Alphanumeric text lines of the global options ending with "\n".
        :return: list of amino acids that will be used to mutate protein sequence to (at every position).
        """
        mutant_aa_list = []
        for line in globaloptions_lines:
            if '#' in line:
                continue
            if Main.Strs.OPT_RESDS.value in line:
                aa_option = Main.__get_text_after_colon_before_semi(line)
                if aa_option.lower() == Main.Strs.OPT_ALL.value:
                    mutant_aa_list = AA.LIST_ALL_20_AA.value
                else:
                    for aa in aa_option:
                        if aa in AA.LIST_ALL_20_AA.value:
                            mutant_aa_list.append(aa)
        return mutant_aa_list

    @staticmethod
    def __is_true(line: str):
        """
        :param line: Alphanumeric text of global options ending with "\n". Expecting "option name:<option value>;\n"
        :return: True boolean if text string is "True", False if text is "False".
        """
        return Main.__get_text_after_colon_before_semi(line).lower() == 'true'

    @staticmethod
    def __get_text_after_colon_before_semi(line: str):
        """
        :param line: Alphanumeric text of global options ending with "\n". Expecting "option name:<option value>;\n"
        :return: all text string after a colon and before a semi-colon, stripping out any white space.
        """
        return line.split(':')[-1].strip(';\n').strip()

    from enum import Enum

    class Strs(Enum):
        OPER_RUN_MUT_FSTA = Str.OPER_RUN_MUT_FSTA.value
        OPER_RUN_AGDR = Str.OPER_RUN_AGDR.value
        OPER_RUN_FX_RPR = Str.OPER_RUN_FX_RPR.value
        OPER_RUN_FX_BM = Str.OPER_RUN_FX_BM.value
        OPER_RUN_FX_AC = Str.OPER_RUN_FX_AC.value
        OPER_RUN_FX_STAB = Str.OPER_RUN_FX_STAB.value
        OPT_MUT_FSTA = 'MUTATE_FASTA:'
        OPT_AGDR = 'AGADIR:'
        OPT_FX_RPR = 'FOLDX_REPAIR:'
        OPT_FX_BM = 'FOLDX_BUILDMODEL:'
        OPT_FX_STAB = 'FOLDX_STABILITY:'
        OPT_FX_AC = 'FOLDX_ANALYSECOMPLEX:'
        OPT_RESDS = 'RESIDUES'
        OPT_ALL = 'all'
        OPT_PDBs = 'PDBs'
        OPT_FASTAs = 'FASTAs'
