import glob
from src.Scheduler import Scheduler
from src.GeneralUtilityMethods import GUM
from src.Paths import Paths
from src.AminoAcids import AA
import natsort

# The following 4 lines of code successfully connect to
# import mysql.connector
# cnx = mysql.connector.connect(user='root', password='K0yGrJ8(', host='127.0.0.1', database='mydb', port='3306')
# cur = cnx.cursor()
# cur.execute("CREATE TABLE testingDBConnection ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")
#
# Currently, the entire set of programs is kicked off by calling the constructor of Main()
# This might be done through a short script like the run_agadir.py script I wrote for the refactored Solubis code.
# Hence run_Main.py could then be called by a cmd line or bash script "python run_Main.py".
class Main(object):

    # At the moment this is set up to run either in local or on cluster, not both. This is temporary.
    def __init__(self, cluster):
        self.pdbs = 'PDBs'
        self.fastas = 'FASTAs'
        globaloptions_lines = Main._read_global_options(Paths.MC_CONFIG_GLOBAL_OPTIONS.value +
                                                        '/MutateCompute_Options.txt')
        wanted_pdbfile_list = Main._build_filelist_for_analysis(globaloptions_lines, self.pdbs, Paths.MC_INPUT.value)
        wanted_fastfile_list = Main._build_filelist_for_analysis(globaloptions_lines, self.fastas, Paths.MC_INPUT.value)
        operations = Main._determine_which_operations_to_perform(globaloptions_lines)
        list_of_mutant_aa = Main._determine_residues_to_mutate_to()
        path_dst_dir = Paths.SE_INPUT.value if cluster else Paths.MC_INPUT.value
        path_src_repo_dir = Paths.SE_REPO_PDB_FASTA.value if cluster else Paths.LOCAL_REPO_PDB_FASTA.value
        src_pdbfile_list = GUM.get_filelist_from_subdirs(path_src_repo_dir, self.pdbs)
        src_fastafile_list = GUM.get_filelist_from_subdirs(path_src_repo_dir, self.fastas)
        wanted_pdbfile_list = GUM.copy_files_from_repo_to_input_filedir(path_src_repo_dir, path_dst_dir, src_pdbfile_list,
                                                                        wanted_pdbfile_list)
        wanted_fastafile_list = GUM.copy_files_from_repo_to_input_filedir(path_src_repo_dir, path_dst_dir,
                                                                          src_fastafile_list, wanted_fastfile_list)
        Main._start_scheduler(operations, Paths.MC_INPUT, wanted_pdbfile_list, wanted_fastafile_list, list_of_mutant_aa)

    # Reads the MutateCompute_Options.txt in /configuration/global_options and .
    # returns the text in a list of strings each ending with a /n (this is the functionality of .readlines()
    # path_globaloptions_file   String      Name of the options file (with extension) and absolute path to the file.
    @staticmethod
    def _read_global_options(path_globaloptions_file):
        with open(path_globaloptions_file, 'r') as globaloptions:
            globaloptions_lines = globaloptions.readlines()
        return globaloptions_lines

    # Takes the /configuration/global_options/MutateCompute_Options.txt text as a list of \n separated
    # lines.
    # Reads the option given at either "PDBs" and "FASTAs".
    # The option is expected to be in one of four forms: "", "all", "<filename(s)>", or "<number>".
    # '' indictes no files - returns None.
    # 'all' indicates all files in specified input source directory.
    # '<filename(s)>' indicates specific files by name (if found in the specified input source directory).
    # '<number>' indicates the maximum number of files from the naturallyy sorted specified input source directory.
    # Note that the sorting of the source directory before building the list of target files could be inefficient if
    # the number of files in the source directory is much larger than the number of files the user has asked for via
    # the options file. I.e. if the source directory contains 30000 pdbs and the user wants the first 5, it will be
    # sorting 30000 pdbs just to get the first 5 of them.
    # E.g. line = "PDBs: 4;\n", returns ['RepairPDB_1.pdb', 'RepairPDB_2.pdb', 'RepairPDB_3.pdb', 'RepairPDB_4.pdb'].
    #
    # globaloptions_lines   List        List of strings, each element is a line of the options text ending with '/n'.
    # PDB_or_FASTA          String      Either "PDBs" or "FASTAs" to indicate which files to retrieve.
    # path_input            String      Path to the input_data files (including PDB and FASTA files).
    # returns a list of pdb or fasta files according to the global_options PDBs or FASTAs. (Files include extensions).
    @staticmethod
    def _build_filelist_for_analysis(globaloptions_lines, PDB_or_FASTA, path_input):
        file_list = []
        for line in globaloptions_lines:
            if '#' in line:
                continue
            if PDB_or_FASTA in line:
                # Main._validate_option_text(line)
                pdb_or_fasta_option = line.split(':')[-1].strip(';\n').strip()
                file_extension = '.pdb' if PDB_or_FASTA == 'PDBs' else '.fasta'
                path_files = path_input + '/**/*' + file_extension
                if pdb_or_fasta_option.lower() == 'all':
                    path_file_list = glob.glob(path_files)
                    for path_file in path_file_list:
                        file_list.append(path_file.split('/')[-1])
                elif pdb_or_fasta_option == '':
                    file_list = None
                else:
                    num_of_files_to_analyse = int(pdb_or_fasta_option)
                    path_file_natsorted_list = natsort.natsorted(glob.glob(path_files))
                    for i, path_file_natsorted in enumerate(path_file_natsorted_list):
                        if i == num_of_files_to_analyse:
                            break
                        file_list.append(path_file_natsorted.split('/')[-1])
                break
        return file_list

    # Validates the "PDBs:" or "FASTAs:" has expected values: "", "all", "filename", "<number within range 1 to 1000>"
    # If the input for PDBs: or FASTAs happens to be "<name>.fasta" or "<name>.pdb", respectively, it is still accepted
    # as the name of the file is the important thing and the extension is ignore anyway.
    # It is also accepted if there was no ";" at the end of the value.
    # Strangely enough, strip(';\n') removes any number of consecutive semi-colons from either side of "\n". However,
    # the validation will fail is there are semi-colons or other non-dot punctuation marks within the value.
    # @staticmethod
    # def _validate_option_text(line):
    #     re.match(['a-zA-Z0-9_'], line)
    #
    # Determines which operations to performed as listed in the global options file, by setting associating it to
    # True or False to indicate whether this operation is to be performed or not, respectively.
    #
    # globaloptions_lines   List        List of strings, each element is a line of the global options ending with /n.
    # returns a dictionary of operation text as key and a boolean as the value.
    @staticmethod
    def _determine_which_operations_to_perform(globaloptions_lines):
        operations = {}
        for line in globaloptions_lines:
            if "#" in line:
                continue
            if "MUTATE_FASTA:" in line:
                operations['do_mutate_fasta'] = line.split(':')[-1].strip(';\n').strip() == 'TRUE'
            if "AGADIR:" in line:
                operations['do_agadir'] = line.split(':')[-1].strip(';\n').strip() == 'TRUE'
            if "FOLDX_REPAIR:" in line:
                operations['do_foldx_repair'] = line.split(':')[-1].strip(';\n').strip() == 'TRUE'
            if "FOLDX_BUILDMODEL:" in line:
                operations['do_foldx_buildmodel'] = line.split(':')[-1].strip(';\n').strip() == 'TRUE'
            if "FOLDX_STABILITY:" in line:
                operations['do_foldx_stability'] = line.split(':')[-1].strip(';\n').strip() == 'TRUE'
            if "FOLDX_ANALYSECOMPLEX:" in line:
                operations['do_foldx_analysecomplex'] = line.split(':')[-1].strip(';\n').strip() == 'TRUE'
        return operations

    # Determines which residues to performed as listed in the global options file, by setting associating it to
    # True or False to indicate whether this operation is to be performed or not, respectively.
    #
    # globaloptions_lines   List        List of strings, each element is a line of the global options ending with /n.
    # returns a list of amino acids that will be used to mutate any protein sequences to (at every position).
    @staticmethod
    def _determine_residues_to_mutate_to(globaloptions_lines):
        mutant_aa_list = []
        for line in globaloptions_lines:
            if '#' in line:
                continue
            if 'RESIDUES' in line:
                aa_option = line.split(':')[-1].strip(';\n').strip()
                if aa_option.lower() == 'all':
                    mutant_aa_list = AA.LIST_ALL_20_AA.value
                else:
                    for aa in aa_option:
                        if aa in AA.LIST_ALL_20_AA.value:
                            mutant_aa_list.append(aa)
        return mutant_aa_list

    @staticmethod
    def _start_scheduler(operations, path_input, pdb_list, fasta_list, list_of_mutant_aa):
        process_started = False
        if operations == {}:
            raise ValueError('All operations options were either set to FALSE or did not read TRUE/FALSE at all. '
                             'Check /configuration/global_options/MutateCompute_Options.txt was written correctly')
        elif operations['do_mutate_fasta'] or operations['do_agadir'] or operations['do_foldx_repair'] \
                or operations['do_foldx_buildmodel'] or operations['do_foldx_stability'] \
                or operations['do_foldx_analysecomplex']:
            process_started = Scheduler.start(operations, path_input, pdb_list, fasta_list, list_of_mutant_aa)
        return process_started

# cnx is the mysql connector (see top of script)
# cnx.close()
