import os
import glob
import sys
import shutil
import subprocess
import time
import warnings
from src.AminoAcids import AA
from src.Paths import Paths
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)


class GUM(object):

    space = ' '
    fslash = '/'

    # The runscript.txt is an input file for FoldX indicating which pdbs to analyse and which programs to run on them.
    # path_runscript        String      Absolute path for runscript.txt file being written.
    # pdbs                  String      pdb(s) (incl. .pdb extension) inputs for FoldX.
    # show_sequence_detail  Boolean     True will provide extra information in output.
    # action                String      FoldX computation to be performed (e.g. BuildModel, AnalyseComplex, etc).
    # print_networks        Boolean     True to output network data on .. ?
    # calculate_stability   Boolean     True include a stability calculation.
    #
    # 30.07.18 Redesigned the directory structure such that runscripts will go in configuration/foldx/ & maybe another
    # level such as analyse_complex or build_model or stability etc.
    @staticmethod
    def write_runscript_for_pdbs(path_runscript, pdbs, action, show_sequence_detail=False, print_networks=False,
                                 calculate_stability=False):
        runscript = []
        runscript.append('<TITLE>FOLDX_runscript;\n')
        runscript.append('<JOBSTART>#;\n')
        runscript.append('<PDBS>' + pdbs + ';\n')
        runscript.append('<BATCH>#;\n')
        runscript.append('<COMMANDS>FOLDX_commandfile;\n')
        if show_sequence_detail:
            runscript.append('<SequenceDetail>#;\n')
        runscript.append(action + ';\n')
        if print_networks:
            runscript.append('<PrintNetworks>#;\n')
        if calculate_stability:
            runscript.append('<Stability>#;\n')
        runscript.append('<END>#;\n')
        runscript.append('<OPTIONS>FOLDX_optionfile;\n')
        runscript.append('<Temperature>298;\n')
        runscript.append('<IonStrength>0.05;\n')
        runscript.append('<ph>7;\n')
        runscript.append('<moveNeighbours>true;\n')
        runscript.append('<VdWDesign>2;\n')
        runscript.append('<numberOfRuns>3;\n')
        runscript.append('<OutPDB>#;\n')
        runscript.append('<END>#;\n')
        runscript.append('<JOBEND>#;\n')
        runscript.append('<ENDFILE>#;\n')
        with open(path_runscript + '/runscript.txt', 'w') as runscript_file:
            runscript_str = ''.join(runscript)
            runscript_file.write(runscript_str)
        return runscript_str

    # Extracts and writes a FASTA file for each chain described in the pdb.
    # Assumes standard pdb format with 'ATOM' as the first string at start of each line of atomic coordinates
    # and with the chain at the 22nd character (index position 21) and the residue number within index 22 to 26.
    # if these very specific aspects are not exactly matching, the method will fail, i.e. it is not very robust.
    # E.g. return value might look like:
    # {'RepairPDB_1_A': 'RVYLTFDELRETK..etc', 'RepairPDB_2_A': 'RVYLTFDELRETK..etc', etc}
    #
    # pdbfiles                  String OR List of Strings   The pdb file (including ".pdb") or files (no path).
    # path_input                String                      Path to where the pdb subdirectories can be found.
    # write_fastafile           Boolean                     True to write fasta string incl. >title to new fasta file.
    # path_fastafile            String                      Absolute path of root for fasta file (without >title) that
    #                                                       will be written if write_fastafile is True.
    #
    # Returns a dictionary of the pdbname_chain and the protein sequence in FASTA format.
    @staticmethod
    def extract_pdbname_chain_fasta_from_pdbs(pdbfiles, path_input, write_fastafile, path_fastafile):
        # pdbname_chain_fastaseq_dict: dict = {} version 3 only
        pdbname_chain_fastaseq_dict = {}
        if isinstance(pdbfiles, str):
            pdbfiles = [pdbfiles]
        for pdbfile in pdbfiles:
            with open(os.path.join(path_input, pdbfile.split('.')[0], pdbfile)) as pdbfile_opened:
                pdbfile_lines = pdbfile_opened.readlines()
            pdbfile_lines_with_atom = []
            protein_chains = []
            for pdbfile_line in pdbfile_lines:
                if 'ATOM' == pdbfile_line[0:4]:
                    protein_chain = pdbfile_line[21]
                    pdbfile_lines_with_atom.append(pdbfile_line)
                    if protein_chain not in protein_chains:
                        protein_chains.append(protein_chain)
            for protein_chain in protein_chains:
                fasta_list = []
                residue_num = '0'
                for pdbfile_line in pdbfile_lines_with_atom:
                    protein_chain_on_this_line = pdbfile_line[21]
                    residue_num_on_this_line = pdbfile_line[22:26].strip(' ')
                    if protein_chain_on_this_line == protein_chain and residue_num_on_this_line != residue_num:
                        residue_num = residue_num_on_this_line
                        amino_acid_on_this_line = pdbfile_line[17:20]
                        if amino_acid_on_this_line in AA.DICT_AA_3TO1.value.keys():
                            fasta_list.append(AA.DICT_AA_3TO1.value[amino_acid_on_this_line])
                        else:
                            print('This 3-letter word is not recognised as 1 of the 20 amino acids! '
                                  'Cannot extract FASTA from ' + pdbfile + ' !')

                fasta_sequence = "".join(fasta_list)
                pdbname = GUM._remove_prefix_and_suffix(pdbfile.split('/')[-1].split('.')[0], 'RepairPDB_', '_1_0')
                pdbname_chain = pdbname + '_' + protein_chain
                print(pdbname_chain + ' : ' + fasta_sequence)
                pdbname_chain_fastaseq_dict[pdbname_chain] = fasta_sequence
                if write_fastafile:
                    GUM.write_fastafile_to_name_chain_dir(pdbname_chain_fastaseq_dict, path_fastafile)
        return pdbname_chain_fastaseq_dict

    # Write the fasta file (with >title) to a new subdirectory in the specified directory. The subdirectory has the same
    # name as the fasta file (which includes the protein chain).
    #
    # pdbname_chain_fasta_dict  Dictionary  The pdbname_chain is the key, the amino acid sequence is the value.
    # path_fasta_root           String      Absolute path of root for fasta sequence file (without >title).
    @staticmethod
    def write_fastafile_to_name_chain_dir(pdbname_chain_fasta_dict, path_to_write_fastafile_root):
        for pdbname_chain, fasta_sequence in pdbname_chain_fasta_dict.items():
            path_output_pdb_pdbchain = GUM.create_dir_tree(path_to_write_fastafile_root, pdbname_chain.split('_')[0],
                                                           pdbname_chain)
            with open(path_output_pdb_pdbchain + pdbname_chain + '.fasta', 'w') as fastafile:
                fastafile.write('>' + pdbname_chain + '\n')
                fastafile.write(fasta_sequence)

    # For the pdb passed here, all chains are read from the pdb file.
    #
    # pdbfile           String      String name of pdb file (incl. '.pdb' extension).
    # path_pdbfile      String      Absolute path to directory where the target pdbfile is.
    #
    # Returns list of protein chains in the pdbfile (as alphabetic single-letter characters).
    @staticmethod
    def extract_all_chains_from_pdb(pdbfile, path_pdbfile):
        with open(os.path.join(path_pdbfile, pdbfile)) as pdbfile_opened:
            pdbfile_lines = pdbfile_opened.readlines()
        protein_chains = []
        for pdbfile_line in pdbfile_lines:
            if 'ATOM' == pdbfile_line[0:4]:
                protein_chain = pdbfile_line[21]
                if protein_chain not in protein_chains:
                    protein_chains.append(protein_chain)
        return protein_chains

    # input_string      String      Any string that you want to trim.
    # prefix            String      The prefix you want to remove from the input string.
    # suffix            String      The suffix you want to remove from the input string.
    #
    # Returns the input string with the specified prefix and suffix trimmed off.
    @staticmethod
    def _remove_prefix_and_suffix(input_string, prefix, suffix):
        trimmed = input_string
        if input_string.startswith(prefix):
            trimmed = input_string.replace(prefix, '')
        if trimmed.endswith(suffix):
            trimmed = trimmed.replace(suffix, '')
        return trimmed

    # Originally used in solubis.py but not used at the moment in this project. Keeping as may be useful.
    # Reads a fasta file for the sequence part only (not including the >title).
    #s
    # path_fastafile     String      Absolute path to the fasta file, including fasta file itself (incl. .fasta ext).
    #
    # Returns the (amino acid) sequence of the fasta file.
    @staticmethod
    def get_sequenceOnly_from_fastafile(path_fastafile):
        with open(path_fastafile, 'r') as fastafile_opened:
            fastafile_lines = fastafile_opened.readlines()
            fasta_seq = fastafile_lines[1] if len(fastafile_lines) == 2 else fastafile_lines[0]
        return fasta_seq

    ###################################################################################################################
    # Builds a directory tree with any number of child nodes. Each new node is only ever one level down (siblings).
    # Adding 3 child nodes produces root-leaf1; root-leaf2; root-leaf3
    #
    #     --------- root ----------
    #     |           |           |
    #   leaf1       leaf2       leaf3
    #
    #####################################################################################################
    #
    # path_root     String                      The path to the root (ideally absolute), "/path/of/root"
    # *args         String or list of Strings   Name or list of names of directory to add to root only.
    @staticmethod
    def create_dir_tree_one_level(path_root, *args):
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
    # path_root     String                      The path to the root (ideally absolute), "/path/of/root"
    # *leaves       String or list of Strings   Name or list of names of directory to add to root, then root-leaf.
    #
    # NOTE: os.makedirs(path) does the same thing but throws an exception if the path already exists.
    @staticmethod
    def create_dir_tree(path_root, *leaves):
        if not os.path.exists(path_root):
            os.makedirs(path_root)
        for leaf in leaves:
            path_root_leaf = os.path.join(path_root, leaf)
            if not os.path.exists(path_root_leaf):
                os.mkdir(path_root_leaf)
            path_root = path_root_leaf
        return path_root

    # From 07Aug, started phasing in use os.makedirs() instead of my own create_dir_tree method (after discovering it!).
    @staticmethod
    def _os_makedirs(path_root, *new_dirs):
        path_root_newdirs = os.path.join(path_root, *new_dirs)
        try:
            os.makedirs(path_root_newdirs)
        except FileExistsError:
            print('Part or all of path already exists. No problemo.')
        return path_root_newdirs

    # path_src_dir          String      The path of source for pdb files to copy from.
    # path_dst_dir          String      The path of destination for pdb files to copy to.
    # starting_num          int         The starting number part of name for destination directory.
    # total_num_to_copy     int         The total number of pdbs to copy&move (also used to destination directory name).
    #                                   It is given as a keyword argument so that a default of 100 can be set when no
    #                                   argument is given.
    #
    # NB: I DON'T THINK I'VE COMPLETED THIS METHOD BECAUSE IT DOESN'T LOOK LIKE IT USES THE STARTING NUMBER TO DECIDE
    # WHICH 100 PDBS TO COPY AND MOVE.***** *** **
    @staticmethod
    def copy_and_move_pdb_files(path_src_dir, path_dst_dir, starting_num, total_num_to_copy=100):
        if os.path.exists(path_src_dir):
            file_list = os.listdir(path_src_dir)
            sorted_file_list = sorted(file_list)
            threedots_subdirname = GUM._make_3dots_subdirname(starting_num, total_num_to_copy)
            path_src_dir_subfoldername = path_src_dir + GUM.fslash + threedots_subdirname
            path_dst_dir_subfoldername = path_dst_dir + GUM.fslash + threedots_subdirname
            target_file_list = []
            for n in range(total_num_to_copy):
                num = n + 1
                target_file_list.append('RepairPDB_' + str(num) + '.pdb')
            if not os.path.exists(path_src_dir_subfoldername):
                os.mkdir(path_src_dir_subfoldername)
            for file in sorted_file_list:
                if file.endswith('.pdb') and file in target_file_list:
                    shutil.move(path_src_dir + GUM.fslash + file, path_src_dir_subfoldername + GUM.fslash + file)
            GUM.linux_copy_all_files_in_dir(path_src_dir_subfoldername, path_dst_dir_subfoldername, recursively=True)

    # Builds a subfolder to house the specified number of pdbs. E.g. if total num to copy is 100 and starting_num is 1:
    # Folder will be name "1...100". But if starting_num is 20025, folder will get name "20025...20125"
    @staticmethod
    def _make_3dots_subdirname(starting_num, total_num_to_copy):
        return str(starting_num) + '...' + str(total_num_to_copy)

    # THESE LINUX COPY METHODS WILL BE COMBINED INTO ONE WITH FLAGS FOR CREATING OWN SUBDIRS BASED ON FILE NAMES AND
    # AND FOR DOING THE COPY RECURSIVELY. A FILELIST WILL ALSO BE INCLUDED BUT CAN BE EMPTY. IT MIGHT ALSO BE WORTH
    # ENFORCING THAT THERE IS A PATH TO THE LOCATION OF THE FILES WITHOUT THE ACTUAL ABS PATH INCLUDING THE FILES
    # THEMSELVES, I.E. THE GLOB.GLOB COULD BE DONE WITHIN THE LINUX COPY METHOD RATHER THAN OUTSIDE BUT THIS COULD
    # ALSO INDICATE THAT HAVING TWO SEPARATE LINUX COPY METHODS MIGHT STILL BE WORTH HAVING RATHER THAN JUST ONE.

    # Copies all files inside specified src dir. Does not include the folder itself.
    #
    # path_src_dir      String      Abs path of source directory to copy files from.
    # path_dst_dir      String      Abs path of destination directory to copy files to.
    # recursively       Boolean     True to copy all files in current directory and all files in all subdirectories.
    #
    # NOTE: current cp command specifies not to overwrite existing files
    @staticmethod
    def linux_copy_all_files_in_dir(path_src_dir, path_dst_dir, recursively=False):
        recurse_cmd = (GUM.space + '-r') if recursively else ''
        cmd = 'cp -n' + recurse_cmd + GUM.space + path_src_dir + '/.' + GUM.space + path_dst_dir
        try:
            subprocess.call(cmd, shell=True)
        except OSError:
            print('Problem with linux cp command.')

    # NOTE: current cp command specifies not to overwrite existing files.
    @staticmethod
    def linux_copy_files(path_src_filelist, path_dst, into_own_subdirs):
        if isinstance(path_src_filelist, str):
            path_src_filelist = [path_src_filelist]
        for path_src_file in path_src_filelist:
            if into_own_subdirs:
                filename = path_src_file.split('/')[-1].split('.')[0]
                path_dst_filenamedir = GUM._os_makedirs(path_dst, filename)
                cmd = 'cp -n' + GUM.space + path_src_file + GUM.space + path_dst_filenamedir
            else:
                cmd = 'cp -n' + GUM.space + path_src_file + GUM.space + path_dst
            try:
                subprocess.call(cmd, shell=True)
            except OSError:
                print('Problem with linux cp command.')

    # Finds files in the immediate specified directory, creates a subdir with same name as the file and moves the file
    # into this new subdirectory.
    #
    # path_dir    String      Abs path of files that will each be moved into their own subdir with same name.
    @staticmethod
    def _move_files_into_own_subdirs(path_dir):
        path_files_in_dir = glob.glob(path_dir + '/*.*')
        for path_file_in_dir in path_files_in_dir:
            if os.path.isfile(path_file_in_dir):
                filename = path_file_in_dir.split('/')[-1].split('.')[0]
                path_dir_filename = GUM._os_makedirs(path_dir, filename)
                cmd = 'mv' + GUM.space + path_file_in_dir + GUM.space + path_dir_filename
                try:
                    subprocess.call(cmd, shell=True)
                except OSError:
                    print('Problem with linux cp command.')

    @staticmethod
    def linux_remove_files_in_dir(path_dir):
        cmd = 'rm' + GUM.space + path_dir + "/*"
        try:
            subprocess.call(cmd, shell=True)
        except OSError:
            print('Problem with linux cp command.')

    # Might be better to log in manually, as this can hang if mulitplexing is used.
    # @staticmethod
    # def log_in_cluster_ssh(*modules):
        # HOST = "zeus.psb.ugent.be"
        # # Ports are handled in ~/.ssh/config since we use OpenSSH
        # COMMAND = "-p 7788 shazib"
        # Ports are handled in ~/.ssh/config since we use OpenSSH
        # ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND], shell=False, stdout=subprocess.PIPE,
        # stderr=subprocess.PIPE)
        # result = ssh.stdout.readlines()
        # if result == []:
        #     error = ssh.stderr.readlines()
        #     print("ERROR: %s" % error, file=sys.stderr)
        # else:
        #     print(result)

    # or use fabric??
    # import os
    # from fabric.api import *
    #
    # HOST="www.example.org"
    # COMMAND="uname -a"
    #
    # env.user = os.getenv('SSH_USER', 'vagrant')
    # env.password = os.getenv('SSH_PASSWORD', 'vagrant')
    #
    # @hosts(HOST)
    # def do_something():
    #     run(COMMAND)

    # Copy files from a source repository subdirectory to destination directory.
    # (Src dir is typically dirs like ~/REPO_PDB_FASTA/pdbs_10 and ~/REPO_PDB_FASTA/fastas_100.
    # Dst dirs are typically MutateCompute/input_data/<pdbname> directories and
    # MutateCompute/input_data/fastas/<fastaname> directories).
    #
    # 1. Remove any files from wanted list that are not found in the source directory.
    # 2. Create /input_data/<pdbname> or input_data/fastas/<fastaname> subdirectories in dst dir for each pdb or fasta
    #    to be copied in to.
    # 3. Copy the pdbfiles or fastafiles to their corresponding subdirectories.
    #
    # E.g. file1.pdb will be copied to /input_data/file1/file.pdb
    # E.g. file1_A.fasta will be copied to /input_data/fastas/file1_A/file1_A.fasta
    #
    # path_repo_subdir              String      Abs path of repository subdir from which to copy pdb or fasta files.
    #                                           It is a subdir of REPO_PDB_FASTA which itself contains only directories
    #                                           of files (not files).
    # path_dst_dir                  String      Path of dest dir to which the specified pdbfiles or fastafiles are
    #                                           copied, via (creating) individual subdirs for each, bearing same name.
    #                                           Typically /input_data/fastas or just /input_data. If the latter, then
    #                                           /fastas will be created as subdir of /input_data.
    # wanted_file_list              List        A subset of files that you want, specified by name (incl. file ext).
    #                                           If this is None or '', it means you want all files in the specified dir.
    #
    # Returns wanted_file_list containing only those files that were found in, and copied from, the path_repo src dir.
    @staticmethod
    def copy_files_from_repo_to_input_dirs(path_repo_pdbs_or_fastas, path_dst_dir, wanted_file_list):
        # determine if
        pdbs_or_fastas = Paths.DIR_PDBS
        if Paths.DIR_FASTAS.value in path_repo_pdbs_or_fastas:
            pdbs_or_fastas = 'FASTAs'
        # add /fastas if path_dst_dir doesn't already have this at the end of this path string.
        if pdbs_or_fastas == Paths.DIR_FASTAS.value:
            if path_dst_dir.split('/')[-1] != Paths.DIR_FASTAS.value:
                path_dst_dir = GUM._os_makedirs(path_dst_dir, 'fastas')
        available_file_list = GUM.get_pdb_or_fastafile_list_from_subdir(path_repo_pdbs_or_fastas)
        if not wanted_file_list or wanted_file_list == '':
            GUM.linux_copy_all_files_in_dir(path_src_dir=path_repo_pdbs_or_fastas, path_dst_dir=path_dst_dir,
                                            recursively=True)
        else:
            for wanted_file in wanted_file_list:
                path_wantedfile = os.path.join(path_repo_pdbs_or_fastas, wanted_file)
                if path_wantedfile not in available_file_list:
                    wanted_file_list.remove(wanted_file)
                else:
                    wanted_file_dst_dirname = wanted_file.split('.')[0]
                    path_dst_dir = GUM._os_makedirs(path_dst_dir, wanted_file_dst_dirname)
                    path_file_to_copy = os.path.join(path_repo_pdbs_or_fastas, wanted_file)
                    GUM.linux_copy_files(path_file_to_copy, path_dst_dir)
        return wanted_file_list

    @staticmethod
    def get_subdirname_starting_with(path, starting_with):
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

    # The repo directory structure is expected to be /REPO_PDB_FASTA/pdbs_<number> and /REPO_PDB_FASTA/fastas_<number>.
    # Note: It assumes only 1 PDB and 1 FASTA subdirectory. If there is more than 1, they won't be seen.
    #
    # path_repo_pdbs_or_fastas  String          Absolute path of the source repository directory.
    #
    # Returns full list of pdbfiles or fastafiles from one of two subdirectories that should be in REPO_PDB_FASTA.
    @staticmethod
    def get_pdb_or_fastafile_list_from_subdir(path_repo_pdbs_or_fastas):
        file_ext = '.pdb' if (Paths.DIR_FASTAS.value in path_repo_pdbs_or_fastas) else '.fasta'
        filelist = glob.glob(path_repo_pdbs_or_fastas + '/*' + file_ext)
        return filelist

    # Note use of list comprehension.
    @staticmethod
    def make_filelist_in_current_dir_and_subdirs_recursively(path_current_dir, dot_file_ext):
        my_list = glob.glob(path_current_dir + '/**/*' + dot_file_ext, recursive=True)
        return [my_file.split('/')[-1] for my_file in my_list]

    @staticmethod
    def using_cluster():
        return 'switchlab' in os.getcwd()

    ######################################################################################################################
    # Permanently removes input_data and all contents!
    @staticmethod
    def remove_inputdata_dir_tree():
        path_to_delete = Paths.INPUT.value
        GUM.do_userWarning_deleting_dir(dir=path_to_delete, lineno=332)
        if os.path.exists(path_to_delete):
            GUM.__delete_subdirectory_tree_of_inputdata()

    # Permanently removes output_data and all contents!
    @staticmethod
    def remove_output_dirs():
        path_to_delete = Paths.OUTPUT.value
        GUM.do_userWarning_deleting_dir(dir=path_to_delete, lineno=340)
        if os.path.exists(path_to_delete):
            GUM.__delete_subdirectory_tree_of_outputdata()

    # Permanently removes config and all contents!
    @staticmethod
    def remove_config_folders():
        path_to_delete = Paths.CONFIG.value
        GUM.do_userWarning_deleting_dir(dir=path_to_delete, lineno=348)
        if os.path.exists(path_to_delete):
            GUM.__delete_subdirectory_tree_of_outputdata()

    @staticmethod
    def do_userWarning_deleting_dir(dir, lineno):
        warnings.warn_explicit(message='You are about to delete ' + dir + 'dir tree. You have 10 secs to abort!',
                               category=UserWarning, filename='GeneralUtilityMethods.py', lineno=lineno)
        time.sleep(10)

    # PRIVATE METHOD
    # Permanently deletes /config and all its contents.
    @staticmethod
    def __delete_subdirectory_tree_of_config():
        try:
            shutil.rmtree(Paths.CONFIG.value)
        except OSError as e:
            print("Error removing: %s - %s." % (e.filename, e.strerror))

    # PRIVATE METHOD
    # Permanently deletes /input_data and all its contents.
    @staticmethod
    def __delete_subdirectory_tree_of_inputdata():
        try:
            shutil.rmtree(Paths.INPUT.value)
        except OSError as e:
            print("Error removing: %s - %s." % (e.filename, e.strerror))

    # PRIVATE METHOD
    # Permanently deletes /output_data and all its contents.
    @staticmethod
    def __delete_subdirectory_tree_of_outputdata():
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
