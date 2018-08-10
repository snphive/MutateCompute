import os
import glob
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

    @staticmethod
    def wait_for_grid_engine_job_to_complete(grid_engine_job_prefix, message_to_print):
        check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
        output_qstat = check_qstat.stdout.read()
        while grid_engine_job_prefix in output_qstat:
            print('Waiting for ' + message_to_print)
            time.sleep(10)
            check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
            output_qstat = check_qstat.stdout.read()

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
        pdbname_chain_fastaseq_dict: dict = {}
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

    # Combines the absolute path given to every fasta file in the given list. Builds ayn missing subdirectories such as
    # /fastas/<fastaname>
    # This separation of path to input directory and fasta files was done to provide more flexibility (it makes it
    # more testable too).
    #
    # path_input        String       Absolute path of directory holding list of fasta files.
    # fastafile_list    List         List of strings which are the target fasta files (including extensions).
    #
    # Returns list of fasta files with absolute path to them, e.g. ~/.../input_data/fastas/1_A/1_A.fasta
    @staticmethod
    def build_complete_paths_for_fastafiles(path_input, fastafile_list):
        return [path_input + '/fastas/' + fastafile.split('.')[0] + '/' + fastafile for fastafile in fastafile_list]

    # Splits a fastafile's text contents into a dictionary. The key is fastafile's title, value is (amino acid)
    # sequence.
    # E.g.
    # ">1_A\n'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVG...' etc"
    #
    # is saved as
    #
    # {
    #   '1_A': <seq of 1_A>,
    #   '2_A': <seq of 2_A>,
    #    etc..
    # }
    #
    #
    # The fastafile is assumed to be either two lines: 1st line ">title", 2nd line sequence, but the method also allows
    # there to be either no title (i.e. sequence only) or an empty title (">" with no title after it). In these both
    # scenarios, the title is created based on the name of the file itself.
    #
    # path_fastafiles   String OR List of Strings   Absolute path to the fasta file, including the file (incl. .fasta)
    #                                               e.g. ~/../input_data/fastas/1_A/1_A.fasta
    #
    # Returns title:sequence dictionary.
    @staticmethod
    def make_titleSeqDict_from_fastafile(path_fastafiles):
        title_sequence_dict = {}
        if isinstance(path_fastafiles, str):
            path_fastafiles     = [path_fastafiles]
        for path_fastafile in path_fastafiles:
            with open(path_fastafile, 'r') as path_fastafile_opened:
                fastafile_lines = path_fastafile_opened.readlines()
                if len(fastafile_lines) > 2 or len(fastafile_lines) < 1:
                    raise ValueError('There is either an extra unexpected carriage return or no sequence data at all')
                elif len(fastafile_lines) == 2:
                    title = GUM._remove_prefix_and_suffix(fastafile_lines[0],  ">", "\n")
                if title == "" or len(fastafile_lines) == 1:
                    title = os.path.splitext(path_fastafiles)[0]
                fasta_seq = fastafile_lines[1] if len(fastafile_lines) == 2 else fastafile_lines[0]
                title_sequence_dict[title] = fasta_seq
        return title_sequence_dict

    # Converts a title:sequence dictionary to a title: {title:seq dictionary} dictionary. One use of this is to add
    # 19 more titles per residue for different titles which represent point mutants.
    # E.g.
    #
    # {
    #   '1_A': <seq of 1_A>,
    #   '2_A': <seq of 2_A>,
    #    etc..
    # }
    #
    # is converted to
    #
    # {
    #   '1_A':
    #           {'1_A'      : <seq of 1_A>,
    #            '1_A_R1A'  : <seq of 1_A_R1A>, ... etc },
    #   '2_A':
    #           {'2_A'      : <seq of 2_A>,
    #            '2_A_P1A'  : <seq of 2_A_P1A>, ... etc },
    #    etc..
    # }
    #
    # title_sequence_dict   Dictionary   Fasta title (i.e. sequence id) is key. Sequence in FASTA format is value.
    #
    # Returns dictionary of title: title:sequence dictionaries.
    @staticmethod
    def convert_titleSeqDict_to_titleTitleSeqDictDict(title_sequence_dict):
        title_titleSeqDict_dict = {}
        for title in title_sequence_dict.keys():
            title_titleSeqDict_dict[title] = title_sequence_dict
        return title_titleSeqDict_dict


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
    # *args         String or list of Strings   Name or list of names of directory to add to root, then root-leaf.
    #
    # NOTE: os.makedirs(path) does the same thing but throws an exception if the path already exists.
    @staticmethod
    def create_dir_tree(path_root, *args):
        if not os.path.exists(path_root):
            os.makedirs(path_root)
        for leaf in args:
            path_root_leaf = os.path.join(path_root, leaf)
            if not os.path.exists(path_root_leaf):
                os.mkdir(path_root_leaf)
            path_root = path_root_leaf
        return path_root_leaf

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
            subfoldername = GUM._make_subfoldername(starting_num, total_num_to_copy)
            path_src_dir_subfoldername = path_src_dir + GUM.fslash + subfoldername
            path_dst_dir_subfoldername = path_dst_dir + GUM.fslash + subfoldername
            target_file_list = []

            for n in range(total_num_to_copy):
                num = n + 1
                target_file_list.append('RepairPDB_' + str(num) + '.pdb')

            if not os.path.exists(path_src_dir_subfoldername):
                os.mkdir(path_src_dir_subfoldername)

            for file in sorted_file_list:

                if file.endswith('.pdb') and file in target_file_list:
                    shutil.move(path_src_dir + GUM.fslash + file, path_src_dir_subfoldername + GUM.fslash + file)

            GUM.linux_copy(path_src_dir_subfoldername, path_dst_dir_subfoldername, do_recursively=True)

    # Builds a subfolder to house the specified number of pdbs. E.g. if total num to copy is 100 and starting_num is 1:
    # Folder will be name "1...100". But if starting_num is 20025, folder will get name "20025...20125"
    @staticmethod
    def _make_subfoldername(starting_num, total_num_to_copy):
        return str(starting_num) + '...' + str(total_num_to_copy)

    @staticmethod
    def linux_copy(path_src, path_dst, do_recursively):
        recurse_cmd = (GUM.space + '-r') if do_recursively else ''
        cmd = 'cp' + recurse_cmd + GUM.space + path_src + GUM.space + path_dst
        subprocess.call(cmd, shell=True)

    # Copy files from a source repository directory to destination directories.
    # (Src dir is typically ~/REPO_PDB_FASTA/PDBs_10 and ~/REPO_PDB_FASTA/FASTAs_100.
    # Dst dirs are typically MutateCompute/input_data/<pdbname> directories and
    # MutateCompute/input_data/fasta/<fastaname> directories).
    # 1. Identify which subfolder (/PDBs_ or /FASTAs_) is src subdirectory in REPO_PDB_FASTA src directory.
    #    Get exact name of the subfolder and set this in the path of src directory (i.e. /PDBs_100).
    # 2. Remove any files from wanted list that are not found in the source directory.
    # 3. Create /input_data/<pdbname> or input_data/fastas/<fastaname> subdirectories in dst dir for each pdb or fasta
    #    to be copied in to.
    # 4. Copy the pdb or fasta files in to their corresponding subdirectories.
    #
    # E.g. file1.pdb will be copied to /input_data/file1/file.pdb
    # E.g. file1_A.fasta will be copied to /input_data/fastas/file1_A/file1_A.fasta
    #
    # path_src_repo_dir     String      Path of repository directory from which to copy pdb files (pdb, fasta, other).
    # path_dst_dir          String      Path of destination directory to which the specified pdbfiles are copied, via
    #                                   (creating) individual subdirectories for each file, bearing the same name.
    # src_file_list
    # wanted_file_list      List        List of strings of pdb files or a list of fasta files (incl. extension) to copy.
    # copy_all_files        Boolean     False by default, otherwise will recursively copy all files from
    #                                   path_src_pdb_dir to path_dst_dir_pdbname
    #                                   (If this is set to true, the wanted_pdb_list is ignored.)
    #
    # Returns wanted_file_list containing only those files that were found in, and copied from, the source directory.
    @staticmethod
    def copy_files_from_repo_to_input_filedir(path_src_repo_dir, path_dst_dir, src_file_list, wanted_file_list,
                                              copy_all_files=False):
        PDBs_or_FASTAs = 'PDBs' if wanted_file_list[0].endswith('.pdb') else 'FASTAs'
        path_src_repo_dir = GUM.get_subdirname_starting_with(path_src_repo_dir, PDBs_or_FASTAs)
        path_root_dst_dir = path_dst_dir if PDBs_or_FASTAs == 'PDBs' else path_dst_dir + '/fastas'
        for wanted_file in wanted_file_list:
            if wanted_file not in src_file_list:
                wanted_file_list.remove(wanted_file)
            else:
                path_dst_dir = GUM.create_dir_tree(path_root_dst_dir, wanted_file.split('.')[0])
                path_file_to_copy = os.path.join(path_src_repo_dir, wanted_file)
                GUM.linux_copy(path_file_to_copy, path_dst_dir, do_recursively=copy_all_files)
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

    # The repo directory structure is expected to be /REPO_PDB_FASTA/PDBs_<number> and /REPO_PDB_FASTA/FASTAs_<number>.
    # Note: It assumes only 1 PDB and 1 FASTA subdirectory. If there is more than 1, the rest will not be seen.
    #
    # path_src_repo_dir     String          Absolute path of the source repository directory.
    # PDBs_or_FASTAs        String          Either 'PDBs' or 'FASTAs'
    #
    # Returns full list of pdb files or fasta files from one of two subdirectories that should be in REPO_PDB_FASTA.
    @staticmethod
    def get_filelist_from_subdirs(path_src_repo_dir, PDBs_or_FASTAs):
        file_ext = '.pdb' if (PDBs_or_FASTAs == 'PDBs') else '.fasta'
        filelist = glob.glob(path_src_repo_dir + '/' + PDBs_or_FASTAs + '*/*' + file_ext)
        return filelist

######################################################################################################################
    # Permanently removes input_data and all contents!
    @staticmethod
    def remove_inputdata_dir_tree():
        path_to_delete = Paths.MC_INPUT.value
        GUM.do_userWarning_deleting_dir(dir=path_to_delete, lineno=332)
        if os.path.exists():
            GUM.__delete_subdirectory_tree_of_inputdata()

    # Permanently removes output_data and all contents!
    @staticmethod
    def remove_output_dirs():
        path_to_delete = Paths.MC_OUTPUT.value
        GUM.do_userWarning_deleting_dir(dir=path_to_delete, lineno=340)
        if os.path.exists(path_to_delete):
            GUM.__delete_subdirectory_tree_of_outputdata()

    # Permanently removes config and all contents!
    @staticmethod
    def remove_config_folders():
        path_to_delete = Paths.MC_CONFIG.value
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
            shutil.rmtree(Paths.MC_CONFIG.value)
        except OSError as e:
            print("Error removing: %s - %s." % (e.filename, e.strerror))

    # PRIVATE METHOD
    # Permanently deletes /input_data and all its contents.
    @staticmethod
    def __delete_subdirectory_tree_of_inputdata():
        try:
            shutil.rmtree(Paths.MC_INPUT.value)
        except OSError as e:
            print("Error removing: %s - %s." % (e.filename, e.strerror))

    # PRIVATE METHOD
    # Permanently deletes /output_data and all its contents.
    @staticmethod
    def __delete_subdirectory_tree_of_outputdata():
        try:
            shutil.rmtree(Paths.MC_OUTPUT.value)
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
