import subprocess
import time
import yaml
import os
import shutil
import multiprocessing as mp
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)


class GUM(object):

    space = ' '
    fslash = '/'
    dict_aa_3to1 = {}
    path_zeus_FoldX_exe: ''
    path_local_FoldX_exe: ''

    with open("/Users/u0120577/PycharmProjects/MutateCompute/config/pathsAndDictionaries.yaml", 'r') as stream:

        try:

            paths_and_dictionaries = yaml.load(stream)
            path_zeus_Foldx_exe = paths_and_dictionaries['ROOT']['path_zeus_FoldX_exe']
            path_local_Foldx_exe = paths_and_dictionaries['ROOT']['path_local_FoldX_exe']

            path_zeus_SnpEffect = paths_and_dictionaries['ROOT']['path_zeus_SnpEffect']
            path_local_MutateCompute = paths_and_dictionaries['ROOT']['path_local_MutateCompute']

            path_rel_src = paths_and_dictionaries['ROOT']['path_rel_src']
            path_rel_Inputs_PDBs = paths_and_dictionaries['ROOT']['path_rel_Inputs_PDBs']
            path_rel_Inputs_Fasta = paths_and_dictionaries['ROOT']['path_rel_Inputs_Fasta']
            path_rel_Inputs_Options = paths_and_dictionaries['ROOT']['path_rel_Inputs_Options']
            path_rel_Outputs_PDBs = paths_and_dictionaries['ROOT']['path_rel_Outputs_PDBs']

            path_rel = paths_and_dictionaries['ROOT']['path_rel']
            path_rel_Cluster = paths_and_dictionaries['ROOT']['path_rel_Cluster']
            path_rel_Agadir = paths_and_dictionaries['ROOT']['path_rel_Agadir']
            path_rel_FoldX = paths_and_dictionaries['ROOT']['path_rel_FoldX']
            path_rel_FoldX_BuildModel = paths_and_dictionaries['ROOT']['path_rel_FoldX_BuildModel']
            path_rel_FoldX_Repair = paths_and_dictionaries['ROOT']['path_rel_FoldX_Repair']
            path_rel_FoldX_AnalyseComplex = paths_and_dictionaries['ROOT']['path_rel_FoldX_AnalyseComplex']

            dict_aa_3to1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']

        except yaml.YAMLError as exc:

            print(exc)

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
    # abs_path_runscript    String      Absolute path for runscript.txt file being written.
    # pdbs                  String      pdb(s) (incl. .pdb extension) inputs for FoldX.
    # show_sequence_detail  Boolean     True will provide extra information in output.
    # action                String      FoldX computation to be performed (e.g. BuildModel, AnalyseComplex, etc).
    # print_networks        Boolean     True to output network data on .. ?
    # calculate_stability   Boolean     True include a stability calculation.
    @staticmethod
    def write_runscript_for_pdbs(path_runscript, pdb, action, show_sequence_detail=False, print_networks=False,
                                 calculate_stability=False):
        runscript = []
        runscript.append('<TITLE>FOLDX_runscript;\n')
        runscript.append('<JOBSTART>#;\n')

        runscript.append('<PDBS>' + pdb + ';\n')
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
    # if these very specific aspects are not exactly matching, the method will fail.
    @staticmethod
    def extract_pdbname_chain_fasta_from_pdb(pdbs, path_Inputs, path_Outputs, write_wt_fasta_files):
        pdbname_chain_fasta_dict = {}

        if isinstance(pdbs, str):
            pdbs = [pdbs]

        for pdb in pdbs:
            pdb_file = open(path_Inputs + '/PDBs/' + pdb).readlines()
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

                        if amino_acid in GUM.dict_aa_3to1.keys():  # else throw some kind of exception or error message?
                            fasta_list.append(GUM.dict_aa_3to1[amino_acid])
                        else:
                            print('This 3-letter word is not recognised as 1 of the 20 amino acids! '
                                  'Cannot extract FASTA from ' + pdb + ' !')

                fasta_sequence = "".join(fasta_list)
                pdbname = GUM._remove_prefix_and_suffix('RepairPDB_', '_1_0', pdb.split('/')[-1].split('.')[0])
                pdbname_chain = pdbname + '_' + protein_chain
                print(pdbname_chain + ' : ' + fasta_sequence)
                pdbname_chain_fasta_dict[pdbname_chain] = fasta_sequence
                if write_wt_fasta_files:
                    GUM.write_fasta(pdbname_chain_fasta_dict, path_Outputs)
        return pdbname_chain_fasta_dict

    @staticmethod
    def write_fasta(pdbname_chain_fasta_dict, path_Outputs):
        for pdbname_chain, fasta_sequence in pdbname_chain_fasta_dict.iteritems():
            fasta_file = open(path_Outputs + 'Fasta/' + pdbname_chain + '.fasta', 'w')
            fasta_file.write('>' + pdbname_chain + '\n')
            fasta_file.write(fasta_sequence)
            fasta_file.close()

    # Extracts All chains that are included in a pdb
    # returns them in a list
    @staticmethod
    def extract_all_chains_from_pdb(pdb, relative_path_to_pdb):
        pdb_file = open(relative_path_to_pdb + pdb).readlines()
        protein_chains = []

        for line in pdb_file:

            if 'ATOM' == line[0:4]:
                protein_chain = line[21]

                if protein_chain not in protein_chains:
                    protein_chains.append(protein_chain)

        return protein_chains

    @staticmethod
    def _remove_prefix_and_suffix(prefix, suffix, pdbname):

        if prefix in pdbname:
            str.replace(pdbname, prefix, '')

        if suffix in pdbname:
            str.replace(pdbname, suffix, '')

        return pdbname

    # Originally used in solubis.py but removed as it was not necessary.
    # Keeping this code for now should it be needed in near future
    @staticmethod
    def get_fasta_sequence(path_to_fasta):

        fasta_lines = open(path_to_fasta, 'r').readlines()
        fasta_list = []

        if len(fasta_lines) == 2:

            for amino_acid in fasta_lines[1]:
                fasta_list.append(amino_acid)

        else:

            for fasta_line in fasta_lines[1:]:

                for amino_acid in fasta_line:
                    fasta_list.append(amino_acid)

        fasta = ''.join(fasta_list)
        return fasta

    # Build a directory tree composed of absolute path of the root and any number of child nodes.
    # 3 leaves => root-leaf; root- leaf; root-leaf
    #
    #     --------- root ----------
    #     |           |           |
    #   leaf1       leaf2       leaf3
    #
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

    # Builds a directory tree with no branches. root-leaf-leaf-leaf
    #
    #    root
    #      |
    #    leaf1
    #      |
    #    leaf2
    #      |
    #    leaf3
    #
    # NOTE: os.makedirs(os.path.join(path, *paths)) already does this, but I found it annoying to mock, possibly
    # some bug when combining mocking with the debugger.
    # it does not create duplicates and does not raise any errors or exceptions.
    @staticmethod
    def create_dir_tree(path_root, *args):

        if not os.path.exists(path_root):
            os.makedirs(path_root)

        for leaf in args:
            path_root_leaf = path_root + '/' + leaf

            if not os.path.exists(path_root_leaf):
                os.mkdir(path_root_leaf)

        return path_root_leaf

    @staticmethod
    def build_input_output_directory_trees(build_local_dir_tree, build_zeus_dir_tree):
        rel_path_list = [GUM.path_rel_src, GUM.path_rel_Inputs_PDBs, GUM.path_rel_Inputs_Fasta,
                         GUM.path_rel_Inputs_Options, GUM.path_rel_Outputs_PDBs]

        if not build_local_dir_tree and not build_zeus_dir_tree:
            raise ValueError('Both options false - neither local or cluster trees will be created')

        else:

            if build_local_dir_tree:
                GUM.create_dir_tree_one_level(GUM.path_local_MutateCompute, rel_path_list)

            if build_zeus_dir_tree:
                GUM.create_dir_tree_one_level(GUM.path_zeus_SnpEffect, rel_path_list)


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

            copy_linux_cmd = 'cp' + GUM.space + '-r' + GUM.space + path_src_dir_subfoldername + GUM.space + \
                             path_dst_dir_subfoldername
            subprocess.call(copy_linux_cmd, shell=True)

    @staticmethod
    def _make_subfoldername(starting_num, total_num_to_copy):
        return str(starting_num) + '...' + str(total_num_to_copy)

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

