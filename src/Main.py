import glob
from src.Scheduler import Scheduler
from src.GeneralUtilityMethods import GUM
from src.Paths import Paths
from src.AminoAcids import AA
import natsort
# The following 4 lines of commented-out code successfully connect to my mysql database and create a table.
# import mysql.connector
# cnx = mysql.connector.connect(user='root', password='K0yGrJ8(', host='127.0.0.1', database='mydb', port='3306')
# cur = cnx.cursor()
# cur.execute("CREATE TABLE testingDBConnection ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")


# Everything, other than the Blast program, runs from here, via the constructor of Main().
#
# Calling the constuctor can be done in the IDE through the short script called KickOff.py.
# Alternatively it can be run from the command line via KickOff.sh, which does 3 essential things:
# 1. Sources .profile script that sets the python environment variable to this MutateCompute project folder.
# 2. Sets path to the virtual environment folder which includes all the downloaded libraries used here.
# 3. Executes KickOff.py via python3 command to indicate that version 3 and not 2 should be used here.
class Main(object):

    # All algorithmic analyses start from here.
    # IdentifyProtein.py which runs Blast analyses is the only program that currently runs independently.
    #
    # use_cluster   Boolean     True if program(s) should run on the cluster (whereupon paths are set to zeus paths).
    #                           Currently set to false by default.
    #
    # At the moment, this is set up to run either in local or on cluster, not both.
    def __init__(self, use_cluster, use_multithreading):
        Paths.set_up_paths(use_cluster)
        self.pdbs = 'PDBs'
        self.fastas = 'FASTAs'
        globaloptions_lines = Main._read_global_options(Paths.CONFIG_GLOBAL_OPTIONS + '/global_options.txt')
        wanted_pdbfile_list = Main._build_filelist_for_analysis(globaloptions_lines, self.pdbs, Paths.REPO_PDB_FASTA)
        wanted_fastafile_list = Main._build_filelist_for_analysis(globaloptions_lines, self.fastas, Paths.REPO_PDB_FASTA)
        operations = Main._determine_which_operations_to_perform(globaloptions_lines)
        list_of_mutant_aa = Main._determine_residues_to_mutate_to(globaloptions_lines)
        path_dst_dir = Paths.INPUT
        path_src_repo_dir = Paths.REPO_PDB_FASTA
        available_pdbfile_list = GUM.get_filelist_from_subdirs(path_src_repo_dir, self.pdbs)
        available_fastafile_list = GUM.get_filelist_from_subdirs(path_src_repo_dir, self.fastas)
        wanted_pdbfile_list = GUM.copy_files_from_repo_to_input_dst_dir(path_src_repo_dir, path_dst_dir,
                                                                        available_pdbfile_list, wanted_pdbfile_list)
        wanted_fastafile_list = GUM.copy_files_from_repo_to_input_dst_dir(path_src_repo_dir, path_dst_dir,
                                                                        available_fastafile_list, wanted_fastafile_list)
        Main._start_scheduler(operations, Paths.INPUT, wanted_pdbfile_list, wanted_fastafile_list, list_of_mutant_aa,
                              use_multithreading)

    # Reads the global_options.txt in /configuration/global_options breaking the text up according to newlines.
    #
    # path_globaloptions_file   String      Abs path to global_options file (incl. extension).
    #
    # Returns the text of the file as a list of strings for each line ending with \n.
    @staticmethod
    def _read_global_options(path_globaloptions_file):
        with open(path_globaloptions_file, 'r') as globaloptions:
            globaloptions_lines = globaloptions.readlines()
        return globaloptions_lines

    # Takes the /configuration/global_options/global_options.txt text as a list of \n separated
    # lines.
    #
    # Reads the option given at either "PDBs" and "FASTAs".
    # The option is expected to be in one of four forms: "", "all", "<filename(s)>", or "<number>".
    # '' indictes no files - returns None.
    # 'all' indicates all files in specified input source directory.
    # '<filename(s)>' indicates specific files by name (if found in the specified input source directory).
    # '<number>' indicates the maximum number of files from the naturally-sorted list that is based on what files are
    # present in the specified input_data source directory.
    # E.g. line = "PDBs: 4;\n", returns ['RepairPDB_1.pdb', 'RepairPDB_2.pdb', 'RepairPDB_3.pdb', 'RepairPDB_4.pdb'].
    #
    # Note the risk of inefficiency where the program spends a lot of time sorting the source directory to build the
    # list of available target files, but only a relatively small number of files will be used E.g. source dir sorts
    # about 30,000 pdbs for user to get access to first 5 only!
    #
    # globaloptions_lines   List        Alphanumeric text lines of the global options ending with "\n".
    # PDB_or_FASTA          String      Either "PDBs" or "FASTAs" to indicate which files to retrieve.
    # path_input            String      Path to the input_data files (including PDB and FASTA files).
    #
    # Returns a list of pdbfiles or fastafiles according to that specified in the global_options/Options file at
    # values given for 'PDBs' and 'FASTAs'. (pdbiles and fastafiles are strings including .pdb and .fasta extensions).
    @staticmethod
    def _build_filelist_for_analysis(globaloptions_lines, PDB_or_FASTA, path_input):
        file_list = []
        file_extension = '.pdb' if PDB_or_FASTA == 'PDBs' else '.fasta'
        path_files = path_input + '/**/*' + file_extension
        for line in globaloptions_lines:
            if '#' in line:
                continue
            if PDB_or_FASTA in line:
                # Main.__validate_option_text(line)
                pdb_or_fasta_option = Main.__get_text_after_colon_before_semi(line)
                if pdb_or_fasta_option == '':
                    file_list = None
                    break
                elif pdb_or_fasta_option.isdigit() and not pdb_or_fasta_option.isalpha():
                    # check that the number is reasonable, i.e. less than 1000 ?
                    num_of_files_to_analyse = int(pdb_or_fasta_option)
                    path_file_natsorted_list = natsort.natsorted(glob.glob(path_files))
                    if not path_file_natsorted_list:
                        raise ValueError('No files were found in the input_data directory tree. Files must be manually '
                                         'placed there')
                    for i, path_file_natsorted in enumerate(path_file_natsorted_list):
                        if i == num_of_files_to_analyse:
                            break
                        file_list.append(path_file_natsorted.split('/')[-1])
                elif pdb_or_fasta_option.isalpha():
                    if pdb_or_fasta_option.lower() == 'all':
                        path_file_list = glob.glob(path_files)
                        for path_file in path_file_list:
                            file_list.append(path_file.split('/')[-1])
                    # else:
                        # if there is a comma-separated list of pdb or fasta names
                else:
                    raise ValueError('The option for pdbs or fastas is not All, a number, a list of sequence names or '
                                     'blank. Something is wrong with the PDBs and/or FASTAs global_options option')
                break
        return file_list

    # Determines which operations to perform, according to the global options file - True or False indicates to perform
    # this operation or not, respectively.
    # Format expected is "True" or "False" comes after ":" and before ";\n"
    #
    # globaloptions_lines   List    Alphanumeric text lines of the global options ending with "\n".
    #
    # Returns dictionary. Key is operation text, Value is boolean.
    @staticmethod
    def _determine_which_operations_to_perform(globaloptions_lines):
        operations = {}
        for line in globaloptions_lines:
            if "#" in line:
                continue
            if "MUTATE_FASTA:" in line:
                operations['do_mutate_fasta'] = Main.__is_true(line)
            if "AGADIR:" in line:
                operations['do_agadir'] = Main.__is_true(line)
            if "FOLDX_REPAIR:" in line:
                operations['do_foldx_repair'] = Main.__is_true(line)
            if "FOLDX_BUILDMODEL:" in line:
                operations['do_foldx_buildmodel'] = Main.__is_true(line)
            if "FOLDX_STABILITY:" in line:
                operations['do_foldx_stability'] = Main.__is_true(line)
            if "FOLDX_ANALYSECOMPLEX:" in line:
                operations['do_foldx_analysecomplex'] = Main.__is_true(line)
        return operations

    # Determines which residues to mutate to according to the global_options file.
    # Given either a list of residues in FASTA) format, or the word "All" to indicate all 20 residues.
    #
    # globaloptions_lines   List    Alphanumeric text lines of the global options ending with "\n".
    #
    # Returns list of amino acids that will be used to mutate protein sequence to (at every position).
    @staticmethod
    def _determine_residues_to_mutate_to(globaloptions_lines):
        mutant_aa_list = []
        for line in globaloptions_lines:
            if '#' in line:
                continue
            if 'RESIDUES' in line:
                aa_option = Main.__get_text_after_colon_before_semi(line)
                if aa_option.lower() == 'all':
                    mutant_aa_list = AA.LIST_ALL_20_AA.value
                else:
                    for aa in aa_option:
                        if aa in AA.LIST_ALL_20_AA.value:
                            mutant_aa_list.append(aa)
        return mutant_aa_list

    @staticmethod
    def _start_scheduler(operations, path_input, pdb_list, fasta_list, list_of_mutant_aa, use_multithreading):
        process_started = False
        if operations == {}:
            raise ValueError("All options in 'operations' were either set to FALSE or there was a typo of some form. "
                             "Check that /configuration/global_options/global_options.txt was written correctly")
        elif operations['do_mutate_fasta'] or operations['do_agadir'] or operations['do_foldx_repair'] \
                or operations['do_foldx_buildmodel'] or operations['do_foldx_stability'] \
                or operations['do_foldx_analysecomplex']:
            process_started = Scheduler.start(operations, path_input, pdb_list, fasta_list, list_of_mutant_aa,
                                              use_multithreading)
        return process_started

    # line   String    Alphanumeric text line of the global options ending with "\n".
    #                  Must have format "option name:<option value>;\n"
    #
    # Returns Boolean True if text is "True", False if text is "False".
    @staticmethod
    def __is_true(line):
        return Main.__get_text_after_colon_before_semi(line) == 'TRUE'

    # line   String    Alphanumeric text line of the global options ending with "\n".
    #                  Must have format "option name:<option value>;\n"
    #
    # Returns all text string after a colon and before a semi-colon, stripping out any white space.
    @staticmethod
    def __get_text_after_colon_before_semi(line):
        return line.split(':')[-1].strip(';\n').strip()

    # Validates the "PDBs:" or "FASTAs:" has expected values: "", "all", "filename", "<number within range 1 to 1000>"
    # If the input for PDBs: or FASTAs happens to be "<name>.fasta" or "<name>.pdb", respectively, it is still accepted
    # as the name of the file is the important thing and the extension is ignore anyway.
    # It is also accepted if there was no ";" at the end of the value.
    # Strangely enough, strip(';\n') removes any number of consecutive semi-colons from either side of "\n". However,
    # the validation will fail is there are semi-colons or other non-dot punctuation marks within the value.
    # @staticmethod
    # def __validate_option_text(line):
    #     re.match(['a-zA-Z0-9_'], line)


# cnx is the mysql connector (see top of script)
# cnx.close()
