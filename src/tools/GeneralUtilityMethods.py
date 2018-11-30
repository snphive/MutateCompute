#!/usr/bin/env python3
"""
Class GeneralUtilityMethods for "general purpose" methods used by more than 1 other class. Includes for linux commands.

Some variable naming conventions that I employ here are (the following are all always strings):
'fasta': FASTA sequence only, hence without any name given on the line above, prefixed with '>'.
'fasta_str': FASTA name and sequence. Hence '>synuclein\nMDVFMKGLS...'
'fastafile': name of the fasta file, including the file extension which is typically '.fasta'
'fastafilename': name of the fasta file, not including the file extension. This might be used to make a file specific
 subdir (i.e. bearing the same name.)
"""
import os
import glob
import shutil
import subprocess
import time
import warnings
from src.enums.Str import Str
from src.enums.AminoAcids import AA
from src.enums.Paths import Paths
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class GUM(object):

    # def __init__(self):

    @staticmethod
    def extract_pdbname_chain_startpos_fasta_from_pdbs(path_pdbfiles, path_fastafiles_to_be_written=''):
        """
        Assumes standard pdb format with 'ATOM' as the first string at start of each line of atomic coordinates
        and with the chain at the 22nd character (index position 21) and the residue number within index 22 to 26.
        NOTE: If these very specific aspects are not exactly matching, the method will fail, i.e. it is not very robust.
        Return value should look like:
        {'RepairPDB_1_A_1': 'RVYLTFDELRETK..etc', 'RepairPDB_2_A': 'RVYLTFDELRETK..etc', etc}

        Extract the pdbname, chain, position of starting residue and fasta sequence from a pdb file. Can also write
        fasta file to specified dst.
        :param path_pdbfiles: str OR List of str, Absolute path to pdbfile(s) (including '.pdb') .
        :param path_fastafiles_to_be_written: Absolute path of fastafiles to be written. Empty string indicates no
        fasta files will be written.
        :return: dict, 'pdbname_chain_startpos' is key, fasta sequence is value.
        """
        # pdbname_chain_fastaseq_dict: dict = {} version 3 only
        pdbname_chain_startpos_fastaseq_dict = {}
        if isinstance(path_pdbfiles, str):
            path_pdbfiles = [path_pdbfiles]
        for path_pdbfile in path_pdbfiles:
            with open(path_pdbfile) as f:
                pdbfile_lines = f.readlines()
            pdbfile_lines_with_atom = []
            protein_chains = []
            protein_chains_startpos = {}
            for pdbfile_line in pdbfile_lines:
                if 'ATOM' == pdbfile_line[0:4]:
                    protein_chain = pdbfile_line[21]
                    pdbfile_lines_with_atom.append(pdbfile_line)
                    if protein_chain not in protein_chains:
                        protein_chains.append(protein_chain)
                        protein_chains_startpos[protein_chain] = pdbfile_line[22:26].strip(' ')
            for protein_chain in protein_chains:
                fasta = []
                residue_num = '0'
                for pdbfile_line in pdbfile_lines_with_atom:
                    protein_chain_on_this_line = pdbfile_line[21]
                    residue_num_on_this_line = pdbfile_line[22:26].strip(' ')
                    if protein_chain_on_this_line == protein_chain and residue_num_on_this_line != residue_num:
                        residue_num = residue_num_on_this_line
                        amino_acid_on_this_line = pdbfile_line[17:20]
                        if amino_acid_on_this_line in AA.DICT_AA_3TO1.value.keys():
                            fasta.append(AA.DICT_AA_3TO1.value[amino_acid_on_this_line])
                        else:
                            print('This 3-letter word is not recognised as 1 of the 20 amino acids! '
                                  'Cannot extract FASTA from ' + path_pdbfile + ' !')

                fasta = "".join(fasta)
                pdbname = GUM._remove_prefix_and_suffix(path_pdbfile.split('/')[-1].split('.')[0], 'RepairPDB_', '_1_0')
                pdbname_chain_startpos = pdbname + '_' + protein_chain + '_' + protein_chains_startpos[protein_chain]
                print(pdbname_chain_startpos + ' : ' + fasta)
                pdbname_chain_startpos_fastaseq_dict[pdbname_chain_startpos] = fasta
                if path_fastafiles_to_be_written != '':
                    GUM.write_fastafile_to_name_chain_dir(pdbname_chain_startpos_fastaseq_dict,
                                                          path_fastafiles_to_be_written)
        return pdbname_chain_startpos_fastaseq_dict

    @staticmethod
    def write_fastafile_to_name_chain_dir(pdbname_chain_fasta: dict, path_to_write_fastafile_root: str):
        """
        Write the fasta file (with >title) to a new subdirectory in the specified directory. The subdirectory has the
        same name as the fasta file (which includes the protein chain).
        :param pdbname_chain_fasta: The pdbname_chain is the key, the amino acid sequence is the value.
        :param path_to_write_fastafile_root: Absolute path of root for fasta sequence file (without >title).
        :return:
        """
        for pdbname_chain, fasta_sequence in pdbname_chain_fasta.items():
            path_output_pdb_pdbchain = GUM.create_dir_tree(path_to_write_fastafile_root, pdbname_chain.split('_')[0],
                                                           pdbname_chain)
            with open(path_output_pdb_pdbchain + pdbname_chain + Str.FSTAEXT.value, 'w') as fastafile:
                fastafile.write('>' + pdbname_chain + '\n')
                fastafile.write(fasta_sequence)

    @staticmethod
    def extract_all_chains_from_pdb(path_pdbfile: str):
        """
        For the pdb passed here, all chains are read from the pdb file.
        :param pdbfile: String name of pdb file (incl. '.pdb' extension).
        :param path_pdbfile: Absolute path to directory where the target pdbfile is.
        :return: list of str, protein chains
        """
        with open(path_pdbfile) as f:
            pdbfile_lines = f.readlines()
        protein_chains = []
        for pdbfile_line in pdbfile_lines:
            if 'ATOM' == pdbfile_line[0:4]:
                protein_chain = pdbfile_line[21]
                if protein_chain not in protein_chains:
                    protein_chains.append(protein_chain)
        return protein_chains

    @staticmethod
    def get_num_of_chains(path_pdbfile: str):
        """

        :param path_pdbfile:
        :return: Number of protein chains (aka molecules) within the specified pdb file.
        """
        return len(GUM.extract_all_chains_from_pdb(path_pdbfile))


    @staticmethod
    def _remove_prefix_and_suffix(input_text: str, prefix: str, suffix: str):
        """
        :param input_text: Any string that you want to trim.
        :param prefix: The prefix you want to remove from the input string.
        :param suffix: The suffix you want to remove from the input string.
        :return: str, input with specified prefix/suffix removed.
        """
        trimmed = input_text
        if input_text.startswith(prefix):
            trimmed = input_text.replace(prefix, '')
        if trimmed.endswith(suffix):
            trimmed = trimmed.replace(suffix, '')
        return trimmed

    @staticmethod
    def get_sequenceOnly_from_fastafile(path_fastafile: str):
        """
        Originally used in solubis.py but not used at the moment in this project. Keeping as may be useful.
        Reads a fastafile for the sequence part only (not including the >title).
        :param path_fastafile: Absolute path to the fasta file, including fasta file itself (incl. .fasta ext).
        :return: str, FASTA sequence
        """
        with open(path_fastafile, 'r') as fastafile_opened:
            fastafile_lines = fastafile_opened.readlines()
            fasta = fastafile_lines[1] if len(fastafile_lines) == 2 else fastafile_lines[0]
        return fasta

    ###################################################################################################################
    # Builds a directory tree with any number of child nodes. Each new node is only ever one level down (siblings).
    # Adding 3 child nodes produces root-leaf1; root-leaf2; root-leaf3
    #
    #     --------- root ----------
    #     |           |           |
    #   leaf1       leaf2       leaf3
    #
    #####################################################################################################
    @staticmethod
    def create_dir_tree_one_level(path_root: str, *args):
        """
        See above for description.
        :param path_root: The path to the root (ideally absolute), "/path/of/root"
        :param args: list or list of str. Name or list of names of directory to add to root only.
        :return: list of str, complete paths of newly made directory tree.
        """
        complete_paths = []
        if not os.path.exists(path_root):
            os.makedirs(path_root)
        for leaf in args:
            path_root_leaf = os.path.join(path_root, leaf)
            complete_paths.append(path_root_leaf)
            if not os.path.exists(path_root_leaf):
                os.mkdir(path_root_leaf)
        if not complete_paths:
            warnings.warn_explicit(message='No new directory trees were created. Returning an empty list.',
                                   category=RuntimeWarning, filename='GeneralUtilityMethods.py', lineno=288)
        return complete_paths

    ###################################################################################################################
    # Builds a directory tree with any number of child nodes. If more than 1 node is added, only the first is a direct
    # child of the root. Subsequent child nodes are added to the last added node, so each is a child of the previous.
    # Adding 3 leaves produces root-leaf1-leaf2-leaf3:
    #
    #    root
    #      |
    #    leaf1
    #      |
    #    leaf2
    #      |
    #    leaf3
    #
    ###################################################################################################################
    #
    # NOTE: os.makedirs(os.path.join(path, *paths)) already does this, but I found it problematic for mocking, possibly
    # due to some bug when combining unittest.mock.patch with the pydevd debugger.
    # It does not create duplicates and does not raise any errors or exceptions.
    #
    # NOTE: os.makedirs(path) does the same thing but throws an exception if the path already exists.
    @staticmethod
    def create_dir_tree(path_root: str, *leaves):
        """
        See above for description.
        :param path_root: The path to the root (ideally absolute), "/path/of/root"
        :param leaves: list or list of str. Name or list of names of directory to add to root, then root-leaf.
        :return: str, complete path of newly made directory tree.
        """
        if not os.path.exists(path_root):
            os.makedirs(path_root)
        for leaf in leaves:
            path_root_leaf = os.path.join(path_root, leaf)
            if not os.path.exists(path_root_leaf):
                os.mkdir(path_root_leaf)
            path_root = path_root_leaf
        return path_root

    # From 07Aug, started phasing in use os.makedirs() instead of my own create_dir_tree() method.
    @staticmethod
    def os_makedirs(path_root: str, *new_dirs):
        """
        :param path_root:
        :param new_dirs:
        :return: str, newly-made directory tree.
        """
        path_root_newdirs = os.path.join(path_root, *new_dirs)
        try:
            os.makedirs(path_root_newdirs)
        except FileExistsError:
            print(Str.PARTALLPATHEXISTS_MSG.value)
        return path_root_newdirs

    @staticmethod
    def copy_files_to_3dot_dir(path_src_dir: str, path_dst_dir: str, file_extension: str, starting_num: int,
                               num_to_copy_per_subdir=1000):
        """
        Copy a chosen number of files (1000 by default) of specified file type from a single, named src dir to a named
        dst dir, but inside newly-made subfolders named '1...1000' for example to store the first 1000 files.
        :param path_src_dir: Absolute path of source dir for files to copy from.
        :param path_dst_dir: Absolute path of destination dir for files to copy to.
        :param file_extension: File extension of files to copy.
        :param starting_num: Number of first file to copy from sorted src dir to dst subdir.
        :param num_to_copy_per_subdir: Number of files to copy. Default of 1000.
        """
        # path_sorted_file_list = natsort.natsorted(glob.glob(path_src_dir + '/*' + file_extension))
        path_sorted_file_list = sorted(glob.glob(path_src_dir + '/*' + file_extension))
        subdir = file_extension.strip('.') + 's_'
        path_dst_dir = os.path.join(path_dst_dir, subdir + str(num_to_copy_per_subdir))
        if not path_sorted_file_list:
            raise FileExistsError('No files found at this path')
        for i, path_file in enumerate(path_sorted_file_list):
            if i % num_to_copy_per_subdir == 0:
                path_dst_3dot_dir = GUM._make_path_3dot_dir(path_dst_dir, starting_num + i, num_to_copy_per_subdir + i)
            cmd = 'cp -n' + Str.SPCE.value + path_file + Str.SPCE.value + path_dst_3dot_dir
            try:
                subprocess.call(cmd, shell=True)
            except OSError:
                print(Str.PROBLNXCMD_MSG.value + cmd)

    @staticmethod
    def _make_path_3dot_dir(path_dst_dir: str, start_num: int, end_num: int):
        """
        Builds a subfolder to house the specified number of pdbs. E.g. if total num to copy is 100 and starting_num is 1:
        Folder will be name "1...100". But if starting_num is 20025, folder will get name "20025...20125"
        :param path_dst_dir: Absolute path to destination directory.
        :param start_num: Starting number of 3dot dir, e.g. 1 in the dir name 1...100
        :param end_num: Ending number of 3dot dir, e.g. 100 in the dir name 1...100
        :return: str, Absolute path of newly made 3dot directory.
        """
        return GUM.os_makedirs(path_dst_dir, str(start_num) + Str.DOTS3.value + str(end_num))

    @staticmethod
    def make_path_fastas_3dots_dirs(path_dst_root: str, path_fastafile: str):
        """
        Copies the 3dot subdir from the source dir (path_fastafile) and builds a dst dir with the same 3dots subdir.
        This method is currently only used for building output dir for MutateFasta.
        Expects path_fastafile to be: ~/PycharmProjects/MutateCompute/input_data/fastas/1...1000/1_A/1_A.fasta
        or ~/PycharmProjects/MutateCompute/input_data/29611_fastas_1000/1...1000/00bdkjhi09234kjn3349_4234.fasta, for
        example.
        As such it builds path of child subdirs of /fastas upto the 1...1000 3dots folder name, not including any
        subdirs from there down, such as /<filename>/filename.fasta or just /filename.fasta.
        :param path_dst_root: Absolute path of root to new directory being made (destination dir).
        :param path_fastafile: Absolute path of file that contain the subdir names to be copied to dst dir.
        :return: str, Absolute path of newly-made directory tree.
        """
        dirs_in_path_fastafile = path_fastafile.split('/')
        path_fastas_3dots_dirs = []
        copy_from_here = False
        for dir_in_path_fastafile in dirs_in_path_fastafile[:-1]:
            if dir_in_path_fastafile == '':
                continue
            if copy_from_here or dir_in_path_fastafile == Paths.DIR_FASTAS.value:
                copy_from_here = True
                if dir_in_path_fastafile == Paths.DIR_FASTAS.value:
                    dir_in_path_fastafile = 'mutants_' + dir_in_path_fastafile
                path_fastas_3dots_dirs.append(dir_in_path_fastafile)
                if Str.DOTS3.value in dir_in_path_fastafile:
                    break
        return GUM.os_makedirs(path_dst_root, '/'.join(path_fastas_3dots_dirs))

    @staticmethod
    def make_path_agadir_3dots_dirs(path_root: str, path_fastafile: str):
        """
        Copies the 3dot subdir from the source dir and builds a dst dir with the same 3dot subdir.
        This method is currently only used for building output dir for Agadir.
        Takes path_fastafile: /any_dir/any_dir/any_dir/<3dots dir>/<filename>[-3]/<mutants>[-2]/<file.fasta>[-1]
        It builds: path_root/agadir/<3dots dir>/<any other subdirs between 3dots_dir and [-3] although none are expected.
        """
        path_fastafile_split_to_list = path_fastafile.split('/')
        path_agadir_3dots_dirs = []
        copy_from_here = False
        for path_dir in path_fastafile_split_to_list[:-3]:
            if path_dir == '':
                continue
            if Str.DOTS3.value in path_dir:
                path_agadir_3dots_dirs.append(Paths.DIR_AGADIR.value)
                path_agadir_3dots_dirs.append(path_dir)
                copy_from_here = True
                continue
            if copy_from_here:
                path_agadir_3dots_dirs.append(path_dir)
        return GUM.os_makedirs(path_root, '/'.join(path_agadir_3dots_dirs))

    @staticmethod
    def make_path_agadir_3dots_filename_mutants_dirs(path_root: str, path_fastafile: str, add_filename_subdir: bool):
        """
        Copies the 3dot subdir from the source dir and builds a dst dir with the same 3dot subdir.
        This method is currently only used for building output dir for Agadir.
        Takes path_fastafile: /any_dir/any_dir/any_dir/<3dots dir>/<filename>[-3]/<mutants>[-2]/<file.fasta>[-1]
        It builds: path_root/agadir/<3dots dir>/<filename> into which all the agadir files would go
        (includes <any other subdirs between 3dots_dir and [-3] although none are expected.)
        :param path_root:
        :param path_fastafile:
        :param add_filename_subdir: True if the dir path should include a dedicated dir name for each file, which is
        essential for Agadir because the dir names of mutants will be the same between different proteins, and the
        output files of Agadir do not state the name of the protein in the filename. (This was a new parameter not
        needed before because the dir structure for input_data used to include filenamedir/mutants/fastafiles
        :return:
        """
        path_fastafile_split_to_list = path_fastafile.split('/')
        path_agadir_3dots_dirs = []
        copy_from_here = False
        for path_dir in path_fastafile_split_to_list[:-1]:
            if path_dir == '':
                continue
            if Str.DOTS3.value in path_dir:
                path_agadir_3dots_dirs.append(Paths.DIR_AGADIR.value)
                path_agadir_3dots_dirs.append(path_dir)
                copy_from_here = True
                continue
            if copy_from_here:
                path_agadir_3dots_dirs.append(path_dir)
        if add_filename_subdir:
            fastafilename = path_fastafile.split('/')[-1].split('.')[0]
            suffix = '_mutants'
            if fastafilename.endswith(suffix):
                n = len(suffix)
                fastafilename = fastafilename[:-n]
            path_agadir_3dots_dirs.append(fastafilename)
        return GUM.os_makedirs(path_root, '/'.join(path_agadir_3dots_dirs))

    @staticmethod
    def make_path_3dots_dirs(path_root: str, path_fastafile: str):
        """
        Copies the 3dot subdir from the source dir and builds a dst dir with the same 3dot subdir.
        This method is currently only used for building input dir for individual fastafiles - hence not currently used.
        Expects path_fastafile to be: ~/PycharmProjects/output_data/mutants_fastas/1...1000/1_A/
        mutants/1_A_mutants.fasta
        As such it builds path of child subdirs of output_data/mutant_fastas upto but not including
        /1_A/mutants/1_A_mutants.fasta
        :param path_root:
        :param path_fastafile:
        :return:
        """
        path_fastafile_split_to_list = path_fastafile.split('/')
        path_input_3dots_dirs = []
        copy_from_here = False
        for path_dir in path_fastafile_split_to_list[:-1]:
            if path_dir == '':
                continue
            if path_dir == Paths.DIR_MUTANTS.value:
                break
            if copy_from_here:
                path_input_3dots_dirs.append(path_dir)
            elif path_dir == Paths.DIR_OUTPUT.value:
                copy_from_here = True
        return GUM.os_makedirs(path_root, '/'.join(path_input_3dots_dirs))

    # THESE LINUX COPY METHODS WILL BE COMBINED INTO ONE WITH FLAGS FOR CREATING OWN SUBDIRS BASED ON FILE NAMES AND
    # AND FOR DOING THE COPY RECURSIVELY. A FILELIST WILL ALSO BE INCLUDED BUT CAN BE EMPTY. IT MIGHT ALSO BE WORTH
    # ENFORCING THAT THERE IS A PATH TO THE LOCATION OF THE FILES WITHOUT THE ACTUAL ABS PATH INCLUDING THE FILES
    # THEMSELVES, I.E. THE GLOB.GLOB COULD BE DONE WITHIN THE LINUX COPY METHOD RATHER THAN OUTSIDE BUT THIS COULD
    # ALSO INDICATE THAT HAVING TWO SEPARATE LINUX COPY METHODS MIGHT STILL BE WORTH HAVING RATHER THAN JUST ONE.

    @staticmethod
    def linux_copy_all_files_in_dir(path_src_dir: str, path_dst_dir: str, recursively=False, files_only=False):
        """
        :param path_src_dir: Absolute path of source directory to copy files from.
        :param path_dst_dir: Absolute path of destination directory to copy files to.
        :param recursively: True to copy all files in current directory and all files in all subdirectories.
        :param files_only: True to copy all files only, not with folders (not sure this is functional).
        """
        recurse_option = '-r' if recursively else ''
        path_src_dir = path_src_dir + Str.FSLSH_ASTRX.value if files_only else path_src_dir
        cmd = 'cp' + Str.SPCE.value + recurse_option + Str.SPCE.value + path_src_dir + Str.SPCE.value + path_dst_dir
        try:
            subprocess.call(cmd, shell=True)
        except OSError:
            print(Str.PROBLNXCMD_MSG.value + cmd)

    # NOTE: current cp command specifies not to overwrite existing files.
    @staticmethod
    def linux_copy_files(path_src_files, path_dst: str, into_own_subdirs: bool):
        """
        :param path_src_files: str OR List of str, Absolute path(s) to src file(s) to copy.
        :param path_dst: Absolute path(s) to dst dir for files to be copied to.
        :param into_own_subdirs: True to make individual dirs for each file (bearing same name as file).
        """
        if isinstance(path_src_files, str):
            path_src_files = [path_src_files]
        for path_src_file in path_src_files:
            if into_own_subdirs:
                filename = path_src_file.split('/')[-1].split('.')[0]
                path_dst_filenamedir = GUM.os_makedirs(path_dst, filename)
                cmd = 'cp -n' + Str.SPCE.value + path_src_file + Str.SPCE.value + path_dst_filenamedir
            else:
                cmd = 'cp -n' + Str.SPCE.value + path_src_file + Str.SPCE.value + path_dst
            try:
                subprocess.call(cmd, shell=True)
            except OSError:
                print(Str.PROBLNXCMD_MSG.value + cmd)

    @staticmethod
    def linux_copy_specified_files(path_src_files, path_dst_dir: str):
        """
        :param path_src_files: str or list of str, Absolute path(s) to src file(s) to copy.
        :param path_dst_dir: Absolute path(s) to dst dir for files to be copied to.
        """
        if isinstance(path_src_files, str):
            path_src_files = [path_src_files]
        for path_src_file in path_src_files:
            cmd = 'cp -f' + Str.SPCE.value + path_src_file + Str.SPCE.value + path_dst_dir
            try:
                subprocess.call(cmd, shell=True)
            except OSError:
                print(Str.PROBLNXCMD_MSG.value + cmd)

    @staticmethod
    def _move_files_into_own_subdirs(path_dir: str):
        """
        Finds files in the immediate specified directory, creates a subdir with same name as the file and moves the
        file into this new subdirectory.
        :param path_dir: Absolute path of files that will each be moved into their own subdir with same name.
        """
        path_files_in_dir = glob.glob(path_dir + '/*.*')
        for path_file_in_dir in path_files_in_dir:
            if os.path.isfile(path_file_in_dir):
                filename = path_file_in_dir.split('/')[-1].split('.')[0]
                path_dir_filename = GUM.os_makedirs(path_dir, filename)
                cmd = 'mv' + Str.SPCE.value + path_file_in_dir + Str.SPCE.value + path_dir_filename
                try:
                    subprocess.call(cmd, shell=True)
                except OSError:
                    print(Str.PROBLNXCMD_MSG.value + cmd)

    @staticmethod
    def linux_remove_all_files_in_dir(path_dir: str):
        """
        :param path_dir:
        """
        cmd = 'rm' + Str.SPCE.value + path_dir + "/*"
        try:
            subprocess.call(cmd, shell=True)
        except OSError:
            print(Str.PROBLNXCMD_MSG.value + cmd)

    @staticmethod
    def linux_remove_file(path_file_to_remove: str):
        """
        Removes the specified file.
        :param path_file_to_remove: Absolute path to file to be deleted. This can also be an expression including '*'
        """
        cmd = 'rm -f' + Str.SPCE.value + path_file_to_remove
        try:
            subprocess.call(cmd, shell=True)
        except FileNotFoundError as fnf:
            print('File not found: ' + path_file_to_remove)
            print(fnf)
        except OSError:
            print(Str.PROBLNXCMD_MSG.value + cmd)

    @staticmethod
    def copy_files_from_repo_to_input_dirs(path_repo_pdbs_or_fastas: str, path_dst_dir: str, wanted_files: list):
        """
        Copy files from a source repository subdirectory to destination directory.
        Src dir typically               ~/REPO_PDB_FASTA/fastas/fastas_1000/1...1000/filename.fasta
        Dst dir typically ~/PycharmProjects/MutateCompute/input_data/fastas/1...1000/<filename>/filename.fasta

        1. Removes any files from wanted list that are not found in the source directory.
        2. Creates new <filename> subdir for each file: e.g. /input_data/<filename> subdirs in dst dir for each file
        3. Copies each file to its new corresponding subdir.
        :param path_repo_pdbs_or_fastas: Absolute path of repository subdir from which to copy pdb or fasta files. It 
        is a subdir of REPO_PDB_FASTA which itself contains only directories of files (not files).
        :param path_dst_dir: Path of dest dir to which the specified pdbfiles or fastafiles are copied, via (creating)
        individual subdirs for each, bearing same name. Typically /input_data/fastas or just /input_data. If the latter,
        then /fastas will be created as subdir of /input_data.
        :param wanted_files: A subset of files that you want, specified by name (incl. file ext).
        :return: wanted_files containing only those files that were found in, and copied from, the path_repo src dir.
        """
        if not wanted_files:
            return
        path_repo_pdbs_or_fastas_split_to_list = path_repo_pdbs_or_fastas.split('/')
        path_dst_subdirs = []
        for path_repo_pdbs_or_fastas_dir in path_repo_pdbs_or_fastas_split_to_list:
            if '' == path_repo_pdbs_or_fastas_dir:
                continue
            if Str.DOTS3.value in path_repo_pdbs_or_fastas_dir:
                path_dst_subdirs.append(path_repo_pdbs_or_fastas_dir)
            if Paths.DIR_FASTAS.value == path_repo_pdbs_or_fastas_dir or \
                    Paths.DIR_PDBS.value == path_repo_pdbs_or_fastas_dir:
                path_dst_subdirs.append(path_repo_pdbs_or_fastas_dir)
        path_dst_dir = GUM.os_makedirs(path_dst_dir, *path_dst_subdirs)
        path_input_fastafiles = []
        path_available_files = GUM.get_pdb_or_fastafiles_from_subdir(path_repo_pdbs_or_fastas)
        for wanted_file in wanted_files:
            path_wantedfile = os.path.join(path_repo_pdbs_or_fastas, wanted_file)
            if path_wantedfile not in path_available_files:
                wanted_files.remove(wanted_file)
            else:
                path_file_to_copy = os.path.join(path_repo_pdbs_or_fastas, wanted_file)
                into_own_subdirs = True
                GUM.linux_copy_files(path_file_to_copy, path_dst_dir, into_own_subdirs)
                path_input_fastafiles.append(GUM._build_path_file(path_dst_dir, wanted_file, into_own_subdirs))
        return path_input_fastafiles

    @staticmethod
    def write_1_fastafile_per_fasta_from_multifastafile(path_dst: str, path_fastafile: str):
        """
        Currently only used for Agadir.
        Reads the multi-fastafiles that are in output_data (to which MutateFasta wrote them), and writes the fastafiles
        as individual fastafiles in the (because that's how I assume agadirwrapper expects them, though I could be
        wrong). The dst dir (path_dst_filename) is ~/PycharmProjects/input_data/mutants_fastas/1...1000/1_A if local;
         /switchlab/group/shazib/SnpEffect/input_data/mutants_fastas/1...1000/1_A if on zeus cluster.
        :param path_dst:
        :param path_fastafile:
        """
        fasta_str = ''
        fastafile = ''
        is_first_line = True
        path_dst_filename = GUM.make_path_agadir_3dots_dirs(path_dst, path_fastafile)
        with open(path_fastafile) as f:
            for line in f.readlines():
                if '>' in line:
                    if not is_first_line:
                        with open(os.path.join(path_dst_filename, fastafile), 'w') as j:
                            j.write(fasta_str)
                    fasta_str = line
                    is_first_line = False
                    fastafile = line.split('>')[-1].split('\n')[0] + Str.FSTAEXT.value
                else:
                    fasta_str += line


    @staticmethod
    def _build_path_file(path_root: str, filename: str, into_own_subdirs: bool):
        """
        :param path_root:
        :param filename:
        :param into_own_subdirs:
        :return:
        """
        filename = filename.split('.')[0] if into_own_subdirs else ''
        return GUM.os_makedirs(path_root, filename, filename)

    @staticmethod
    def get_subdirname_starting_with(path: str, starting_with: str):
        """
        :param path:
        :param starting_with:
        :return:
        """
        dirnames = []
        dirname = ''
        for dir_tuple in os.walk(path):
            if dir_tuple[0].split('/')[-1].startswith(starting_with):
                dirnames.append(dir_tuple[0])
                dirname = ''.join(dirnames[0])
                break
        if dirname == '':
            raise ValueError('No directory with a name starting with ' + starting_with + ' was found in ' + path)
        return dirname

    @staticmethod
    def get_pdb_or_fastafiles_from_subdir(path_repo_pdbs_or_fastas: str):
        """
        The repo directory structure is expected to be /REPO_PDB_FASTA/pdbs_<number> and /REPO_PDB_FASTA/fastas_<number>.
        Note: It assumes only 1 PDB and 1 FASTA subdirectory. If there is more than 1, they won't be seen.
        :param path_repo_pdbs_or_fastas: Absolute path of the source repository directory.
        :return: Full list of pdbfiles or fastafiles from one of two subdirectories that should be in REPO_PDB_FASTA.
        """
        file_ext = Str.PDBEXT.value if (Paths.DIR_PDBS.value in path_repo_pdbs_or_fastas) else Str.FSTAEXT.value
        return glob.glob(path_repo_pdbs_or_fastas + '/*' + file_ext)

    @staticmethod
    def get_3dots_dir(path_with_3dots: str):
        """
        :param path_with_3dots:
        :return:
        """
        dir_3dots = ''
        for this_dir in path_with_3dots.split('/'):
            if Str.DOTS3.value in this_dir:
                dir_3dots = this_dir
                break
        return dir_3dots

    # Note the use of list comprehension.
    @staticmethod
    def make_filelist_in_current_dir_and_subdirs_recursively(path_current_dir: str, dot_file_ext: str):
        """
        :param path_current_dir:
        :param dot_file_ext:
        :return:
        """
        my_list = glob.glob(path_current_dir + '/**/*' + dot_file_ext, recursive=True)
        return [my_file.split('/')[-1] for my_file in my_list]

    @staticmethod
    def using_cluster():
        """
        :return: True if current working directory is the cluster, based on presence of string "switchlab".
        """
        return 'switchlab' in os.getcwd()

    ###################################################################################################################

    @staticmethod
    def remove_inputdata_dir_tree():
        """
        Permanently removes input_data and all contents!
        """
        path_to_delete = Paths.INPUT.value
        GUM.give_userWarning_deleting_dir(directoryname=path_to_delete, lineno=332)
        if os.path.exists(path_to_delete):
            GUM.__delete_subdirectory_tree_of_inputdata()

    @staticmethod
    def remove_output_dirs():
        """
        Permanently removes output_data and all contents!
        """
        path_to_delete = Paths.OUTPUT.value
        GUM.give_userWarning_deleting_dir(directoryname=path_to_delete, lineno=340)
        if os.path.exists(path_to_delete):
            GUM.__delete_subdirectory_tree_of_outputdata()

    @staticmethod
    def remove_config_folders():
        """
        Permanently removes config and all contents!
        """
        path_to_delete = Paths.CONFIG.value
        GUM.give_userWarning_deleting_dir(directoryname=path_to_delete, lineno=348)
        if os.path.exists(path_to_delete):
            GUM.__delete_subdirectory_tree_of_outputdata()

    @staticmethod
    def give_userWarning_deleting_dir(directoryname: str, lineno: int):
        """
        :param directoryname:
        :param lineno:
        :return:
        """
        warnings.warn_explicit(message='You are about to delete ' + directoryname + 'dir tree. You have 10 secs to '
                                                                                    'abort!', category=UserWarning,
                               filename='GeneralUtilityMethods.py', lineno=lineno)
        time.sleep(10)

    # PRIVATE METHOD
    @staticmethod
    def __delete_subdirectory_tree_of_config():
        """
        Permanently deletes /config and all its contents.
        """
        try:
            shutil.rmtree(Paths.CONFIG.value)
        except OSError as e:
            print("Error removing: %s - %s." % (e.filename, e.strerror))

    # PRIVATE METHOD
    @staticmethod
    def __delete_subdirectory_tree_of_inputdata():
        """
        Permanently deletes /input_data and all its contents.
        """
        try:
            shutil.rmtree(Paths.INPUT.value)
        except OSError as e:
            print("Error removing: %s - %s." % (e.filename, e.strerror))

    # PRIVATE METHOD
    @staticmethod
    def __delete_subdirectory_tree_of_outputdata():
        """
        Permanently deletes /output_data and all its contents.
        """
        try:
            shutil.rmtree(Paths.OUTPUT.value)
        except OSError as e:
            print("Error removing: %s - %s." % (e.filename, e.strerror))

    # @staticmethod
    # def perform_local_multithreading(process_method_name, file_to_process):
    #     # init objects
    #     cores = mp.cpu_count()
    #     pool = mp.Pool(cores)
    #     jobs = []
    #     # create jobs
    #     with open(file_to_process) as f:
    #         for line in f:
    #             jobs.append(pool.apply_async(process_method_name(line)))
    #     # wait for all jobs to finish
    #     for job in jobs:
    #         job.get()
    #     # clean up
    #     pool.close()
