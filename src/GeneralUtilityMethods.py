import subprocess
import time
import os
import shutil
from src.AminoAcids import AA
import warnings
from src.Paths import Paths
import multiprocessing as mp
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
    #
    # pdbs              String or List of Strings   The pdb file (including ".pdb") or files
    # path_input        String                      Path to where the pdbs can be found.
    # write_fasta       Boolean                     True to write out as a side effect of the method.
    # path_fasta_root   String                      Absolute path of root for fasta sequence file (without >title).
    @staticmethod
    def extract_pdbname_chain_fasta_from_pdb(pdbs, path_input_pdbs, write_fasta_seq, path_fasta_root):
        pdbname_chain_fasta_dict: dict = {}

        if isinstance(pdbs, str):
            pdbs = [pdbs]

        for pdb in pdbs:
            pdb_file = open(path_input_pdbs + pdb).readlines()
            atom_lines = []
            protein_chains = []

            for line in pdb_file:

                if 'ATOM' == line[0:4]:
                    protein_chain = line[21]
                    atom_lines.append(line)

                    if protein_chain not in protein_chains:
                        protein_chains.append(protein_chain)

            for protein_chain in protein_chains:
                fasta_list = []
                resnum = '0'

                for line in atom_lines:

                    if line[21] == protein_chain and resnum != line[22:26].strip(' '):
                        resnum = line[22:26].strip(' ')
                        amino_acid = line[17:20]

                        if amino_acid in AA.DICT_AA_3TO1.value.keys():
                            fasta_list.append(AA.DICT_AA_3TO1.value[amino_acid])
                        else:
                            print('This 3-letter word is not recognised as 1 of the 20 amino acids! '
                                  'Cannot extract FASTA from ' + pdb + ' !')

                fasta_sequence = "".join(fasta_list)
                pdbname = GUM._remove_prefix_and_suffix('RepairPDB_', '_1_0', pdb.split('/')[-1].split('.')[0])
                pdbname_chain = pdbname + '_' + protein_chain
                print(pdbname_chain + ' : ' + fasta_sequence)
                pdbname_chain_fasta_dict[pdbname_chain] = fasta_sequence
                if write_fasta_seq:
                    GUM.write_fasta(pdbname_chain_fasta_dict, path_fasta_root)
        return pdbname_chain_fasta_dict

    # pdbname_chain_fasta_dict  Dictionary  The pdbname_chain is the key, the amino acid sequence is the value.
    # path_fasta_root           String      Absolute path of root for fasta sequence file (without >title).
    @staticmethod
    def write_fasta(pdbname_chain_fasta_dict, path_fasta_root):
        for pdbname_chain, fasta_sequence in pdbname_chain_fasta_dict.items():
            path_output_pdb_pdbchain = GUM.create_dir_tree(path_fasta_root, pdbname_chain.split('_')[0], pdbname_chain)
            fasta_file = open(path_output_pdb_pdbchain + pdbname_chain + '.fasta', 'w')
            fasta_file.write('>' + pdbname_chain + '\n')
            fasta_file.write(fasta_sequence)
            fasta_file.close()

    # For the pdb passed here, all chains are read from the pdb file and returned as a list of chain character(s).
    #
    # pdb               pdb (io file?)      The pdb file (with ",pdb" extension).
    # rel_path_to_pdb   String              The relative path to the pdb.
    @staticmethod
    def extract_all_chains_from_pdb(pdb, rel_path_to_pdb):
        pdb_file = open(rel_path_to_pdb + pdb).readlines()
        protein_chains = []

        for line in pdb_file:

            if 'ATOM' == line[0:4]:
                protein_chain = line[21]

                if protein_chain not in protein_chains:
                    protein_chains.append(protein_chain)

        return protein_chains

    # prefix    String      The prefix to remove
    # suffix    String      The suffix to remove
    # pdbname   String      The pdb name (without ".pdb" extension)
    @staticmethod
    def _remove_prefix_and_suffix(prefix, suffix, pdbname):

        if prefix in pdbname:
            str.replace(pdbname, prefix, '')

        if suffix in pdbname:
            str.replace(pdbname, suffix, '')

        return pdbname

    # Originally used in solubis.py but removed as it was not necessary.
    # Keeping this code for now should it be needed in near future
    # Reads a FASTA file and returns the sequence only (not including the >title).
    #
    # path_to_fasta     String      The path to the fasta file.
    @staticmethod
    def get_sequence_only_from_fasta_file(path_to_fasta):
        with open(path_to_fasta, 'r').readlines() as fasta_lines:
            if len(fasta_lines) > 2 or len(fasta_lines) < 1:
                raise ValueError('There is either an extra unexpected carriage return or no sequence data at all')
            fasta_seq = fasta_lines[1] if len(fasta_lines) == 2 else fasta_lines[0]
        return fasta_seq

    @staticmethod
    def get_title_sequence_dict_from_fasta_file(path_to_fasta):
        title_sequence_dict = {}
        with open(path_to_fasta, 'r').readlines() as fasta_lines:
            if len(fasta_lines) > 2 or len(fasta_lines) < 1:
                raise ValueError('There is either an extra unexpected carriage return or no sequence data at all')
            title = fasta_lines[0] if len(fasta_lines) == 2 else os.path.splitext(path_to_fasta)[0]
            fasta_seq = fasta_lines[1] if len(fasta_lines) == 2 else fasta_lines[0]
            title_sequence_dict[title] = fasta_seq
        return title_sequence_dict

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
            path_root_leaf = path_root + '/' + leaf
            complete_paths.append(path_root_leaf)

            if not os.path.exists(path_root_leaf):
                os.mkdir(path_root_leaf)

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
    @staticmethod
    def create_dir_tree(path_root, *args):
        path_root_leaf = path_root

        if not os.path.exists(path_root):
            os.makedirs(path_root)

        for leaf in args:
            path_root_leaf = path_root + '/' + leaf

            if not os.path.exists(path_root_leaf):
                os.mkdir(path_root_leaf)
            path_root = path_root_leaf

        if path_root_leaf == path_root:
            warnings.warn_explicit(message='No new directory tree was created. Returning the unchanged root.',
                                   category=RuntimeWarning, filename='GeneralUtilityMethods.py', lineno=260)
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
    #
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
    # Folder will be name "1...100". If starting_num is 20025, it will name it "20025...20125"
    @staticmethod
    def _make_subfoldername(starting_num, total_num_to_copy):
        return str(starting_num) + '...' + str(total_num_to_copy)

    # Reads FASTA input file in and returns the text including the >name, the /n, and the amino acid sequence.
    #
    # path_fasta_file       String      FASTA file (with .fasta extension) and absolute path to the file.
    @staticmethod
    def read_fasta_file(path_fasta_file):
        with open(path_fasta_file) as fasta_io:
            fasta_str = fasta_io.read()
        return fasta_str

    @staticmethod
    def linux_copy(path_src, path_dst, do_recursively):
        recurse_cmd = (GUM.space + '-r') if do_recursively else ''
        cmd = 'cp' + recurse_cmd + GUM.space + path_src + GUM.space + path_dst
        subprocess.call(cmd, shell=True)

    # Moves files from one directory to another list of directories that bear the same name.
    # Typically this method is used for getting a list of pdbs from a directory and copying them over to the
    # directories in the input_data folder which should already have a folder for each of the listed directories
    # which are passed as
    @staticmethod
    def copy_input_files_from_repo_to_input(path_src_dir, path_dst_dir_list):
        file_list = os.listdir(path_src_dir)
        for path_dst_dir in path_dst_dir_list:

            GUM.linux_copy()

            GUM.create_dir_tree_one_level(path_root_input, input_pdb_list, input_fasta_list)








    # BEWARE OF THE REMOVE METHODS AS THEY WILL PERMANENTLY DELETE THE SPECIFIED INPUT, OUTPUT OR CONFIG FOLDER!
    @staticmethod
    def remove_input_dirs():
        time
        if os.path.exists(Paths.MC_INPUT.value):
            GUM.__delete_subdirectory_tree_of_tests(Paths.MC_INPUT.value)

    # BEWARE OF THE REMOVE METHODS AS THEY WILL PERMANENTLY DELETE THE SPECIFIED INPUT, OUTPUT OR CONFIG FOLDER!
    @staticmethod
    def remove_output_dirs():
        if os.path.exists(Paths.MC_OUTPUT.value):
            GUM.__delete_subdirectory_tree_of_tests(Paths.MC_OUTPUT.value)

    # BEWARE OF THE REMOVE METHODS AS THEY WILL PERMANENTLY DELETE THE SPECIFIED INPUT, OUTPUT OR CONFIG FOLDER!
    @staticmethod
    def remove_config_folders():
        if os.path.exists(Paths.MC_CONFIG.value):
            GUM.__delete_subdirectory_tree_of_tests(Paths.REL_CONFIG.value)

    # NEVER EVER CALL THIS PRIVATE METHOD FROM ANYWHERE. IT IS DESIGNED TO ONLY BE CALLED BY THE REMOVE METHODS ABOVE.
    @staticmethod
    def __delete_subdirectory_tree_of_tests(path_to_delete):
        os.chdir(Paths.MC_TESTS.value)
        if not os.getcwd() == Paths.MC_TESTS.value:
            raise ValueError('Current working directory is not MutateCompute/tests. '
                             '\nNot proceeding with deletion of ' + path_to_delete)
        else:
            try:
                shutil.rmtree('/' + path_to_delete)
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
