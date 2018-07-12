import subprocess
import time
import yaml
import os
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)


class GUM(object):

    dict_aa_3to1 = {}
    path_foldx_exe = ''

    with open("/Users/u0120577/PycharmProjects/MutateCompute/config/pathsAndDictionaries.yaml", 'r') as stream:

        try:

            paths_and_dictionaries = yaml.load(stream)
            dict_aa_3to1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']
            path_foldx_exe = paths_and_dictionaries['ROOT']['path_FoldX_exe']

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

    # The runscript.txt is an input file for FoldX that informs the which pdbs to analyse and which programs to run
    @staticmethod
    def write_runscript_for_pdbs(path_runscript, pdbs, show_sequence_detail, action, print_networks,
                                 calculate_stability):
        runscript_file = open(path_runscript + 'runscript.txt', 'w')
        runscript_file.write('<TITLE>FOLDX_runscript;\n')
        runscript_file.write('<JOBSTART>#;\n')
        runscript_file.write('<PDBS>' + pdbs + ';\n')
        runscript_file.write('<BATCH>#;\n')
        runscript_file.write('<COMMANDS>FOLDX_commandfile;\n')
        if show_sequence_detail:
            runscript_file.write('<SequenceDetail>#;\n')
        runscript_file.write(action + ';\n')
        if print_networks:
            runscript_file.write('<PrintNetworks>#;\n')
        if calculate_stability:
            runscript_file.write('<Stability>#;\n')
        runscript_file.write('<END>#;\n')
        runscript_file.write('<OPTIONS>FOLDX_optionfile;\n')
        runscript_file.write('<Temperature>298;\n')
        runscript_file.write('<IonStrength>0.05;\n')
        runscript_file.write('<ph>7;\n')
        runscript_file.write('<moveNeighbours>true;\n')
        runscript_file.write('<VdWDesign>2;\n')
        runscript_file.write('<numberOfRuns>3;\n')
        runscript_file.write('<OutPDB>#;\n')
        runscript_file.write('<END>#;\n')
        runscript_file.write('<JOBEND>#;\n')
        runscript_file.write('<ENDFILE>#;\n')
        runscript_file.close()

    # Extracts and writes a FASTA file for each chain described in the pdb.
    # Assumes standard pdb format with 'ATOM' as the first string at start of each line of atomic coordinates
    # and with the chain at the 22nd character (index position 21) and the residue number within index 22 to 26.
    # if these very specific aspects are not exactly matching, the method will fail.
    @staticmethod
    def extract_pdb_name_fasta_chains_from_pdb(pdbs, abs_path_inputs, abs_path_outputs, write_wt_fasta_files):
        pdb_name_chain_fasta_dict = {}
        if isinstance(pdbs, str):
            pdbs = [pdbs]
        for pdb in pdbs:
            pdb_file = open(abs_path_inputs + '/PDBs/' + pdb).readlines()
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
                pdb_name = GUM._remove_prefix_suffix('RepairPDB_', '_1_0', pdb.split('/')[-1].split('.')[0])
                pdb_name_chain = pdb_name + '_' + protein_chain
                print(pdb_name_chain + ' : ' + fasta_sequence)
                pdb_name_chain_fasta_dict[pdb_name_chain] = fasta_sequence
                if write_wt_fasta_files:
                    GUM.write_fasta(pdb_name_chain_fasta_dict, abs_path_outputs)
        return pdb_name_chain_fasta_dict

    @staticmethod
    def write_fasta(pdb_name_chain_fasta_dict, abs_path_outputs):
        for pdb_name_chain, fasta_sequence in pdb_name_chain_fasta_dict.iteritems():
            fasta_file = open(abs_path_outputs + 'Fasta/' + pdb_name_chain + '.fasta', 'w')
            fasta_file.write('>' + pdb_name_chain + '\n')
            fasta_file.write(fasta_sequence)
            fasta_file.close()

    # Extracts All chains that are included in a pdb
    # returns them in a list
    @staticmethod
    def extract_all_chains_from_pdb(pdb, relative_path_to_pdb):
        cwd = os.getcwd()  # for debugging
        pdb_file = open(relative_path_to_pdb + pdb).readlines()
        protein_chains = []
        for line in pdb_file:
            if 'ATOM' == line[0:4]:
                protein_chain = line[21]
                if protein_chain not in protein_chains:
                    protein_chains.append(protein_chain)
        return protein_chains

    @staticmethod
    def _remove_prefix_suffix(prefix, suffix, pdb_name):
        if prefix in pdb_name:
            str.replace(pdb_name, prefix, '')
        if suffix in pdb_name:
            str.replace(pdb_name, suffix, '')
        return pdb_name

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
        fasta = "".join(fasta_list)
        return fasta

    # Build a directory tree composed of absolute path of the root and any number of child nodes.
    @staticmethod
    def create_directory_tree(absolute_path_root, *args):
        if not os.path.exists(absolute_path_root):
            os.mkdir(absolute_path_root)
        os.chdir(absolute_path_root)
        for child_path in args:
            if not os.path.exists(child_path):
                os.mkdir(child_path)
            os.chdir(child_path)
        return os.getcwd()
