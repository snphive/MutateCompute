import glob
from src.Scheduler import Scheduler
from src.GeneralUtilityMethods import GUM
from src.Paths import Paths
from src.AminoAcids import AA


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
        globaloptions = Main._read_global_options()
        input_pdb_list = Main._determine_which_pdbs_to_analyse(globaloptions)
        input_fasta_list = Main._determine_which_fasta_to_analyse(globaloptions)
        operations = Main._determine_which_operations_to_perform(globaloptions)
        list_of_mutant_aa = Main._determine_residues_to_mutate_to()
        path_dst_pdb_dir = Paths.SE_INPUT if cluster else Paths.MC_INPUT
        path_src_pdb_dir = Paths.SE_PDB_REPO if cluster else Paths.LOCAL_PDB_REPO
        GUM.copy_input_files_from_repo_to_input(path_src_pdb_dir, path_dst_pdb_dir, input_pdb_list, input_fasta_list)
        Main._start_scheduler(operations, Paths.MC_INPUT, input_pdb_list, input_fasta_list, list_of_mutant_aa)

    @staticmethod
    def _read_global_options():
        return open(Paths.MC_CONFIG_GLOBAL_OPTIONS + '/MutateCompute_Options.txt', 'r').readlines()

    @staticmethod
    def _determine_which_pdbs_to_analyse(global_options):
        pdb_path_list = []
        number_of_pdbs_to_analyse = 0
        input_pdb_list = []
        for line in global_options:
            if '#' in line:
                continue
            if 'PDBs:' in line:
                pdb_option = line.split(':')[-1].strip(';\n').strip().lower()
                if pdb_option == 'all'  or pdb_option == '':
                    pdb_path_list = glob.glob(Paths.MC_INPUT + '*.pdb')
                    for pdb_path in pdb_path_list:
                        input_pdb_list.append(pdb_path.split('/')[-1])
                elif pdb_option == '':
                    number_of_pdbs_to_analyse = 0
                elif not isinstance(pdb_option, str):
                    number_of_pdbs_to_analyse = int(pdb_option)
                for x in range(1, number_of_pdbs_to_analyse + 1):
                    input_pdb_list.append(pdb_path_list[x].split('/')[-1])
        return input_pdb_list

    @staticmethod
    def _determine_which_fasta_to_analyse(global_options):
        fasta_path_list = []
        number_of_fastas_to_analyse = 0
        input_fasta_list = []
        for line in global_options:
            if "#" in line:
                continue
            if 'FASTAs:' in line:
                fasta_option = line.split(':')[-1].strip(';\n').strip().lower()
                if fasta_option == 'all' or fasta_option == '':
                    fasta_path_list = glob.glob(Paths.MC_INPUT + '*.fasta')
                    for fasta_path in fasta_path_list:
                        input_fasta_list.append(fasta_path.split('/')[-1])
                elif fasta_option == '':
                    number_of_fastas_to_analyse = 0
                elif not isinstance(fasta_option, str):
                    number_of_fastas_to_analyse = int(fasta_option)
                for x in range(1, number_of_fastas_to_analyse + 1):
                    input_fasta_list.append(fasta_path_list[x].split('/')[-1])
        return input_fasta_list

    @staticmethod
    def _determine_which_operations_to_perform(global_options):
        operations = {}
        for line in global_options:
            if "#" in line:
                continue
            if "MUTATE_FASTA:" in line:
                operations['do_mutate_fasta'] = line.split(':')[-1].strip(';\n').strip() == 'True'
            if "AGADIR:" in line:
                operations['do_agadir'] = line.split(':')[-1].strip(';\n').strip() == 'True'
            if "FOLDX_REPAIR:" in line:
                operations['do_foldx_repair'] = line.split(':')[-1].strip(';\n').strip() == 'True'
            if "FOLDX_BUILDMODEL:" in line:
                operations['do_foldx_buildmodel'] = line.split(':')[-1].strip(';\n').strip() == 'True'
            if "FOLDX_STABILITY:" in line:
                operations['do_foldx_stability'] = line.split(':')[-1].strip(';\n').strip() == 'True'
            if "FOLDX_ANALYSECOMPLEX:" in line:
                operations['do_foldx_analysecomplex'] = line.split(':')[-1].strip(';\n').strip() == 'True'
        return operations

    @staticmethod
    def _determine_residues_to_mutate_to(global_options):
        list_of_mutant_aa = []
        for line in global_options:
            if '#' in line:
                continue
            if 'RESIDUES' in line:
                aa_option = line.split(':')[-1].strip(';\n').strip()
                if aa_option == 'all' or aa_option == '':
                    list_of_mutant_aa = AA.LIST_ALL_20_AA
                else:
                    for aa in aa_option:
                        list_of_mutant_aa.append(aa)
        return list_of_mutant_aa

    @staticmethod
    def _start_scheduler(operations, path_input, input_pdb_list, input_fasta_list, list_of_mutant_aa):
        if operations['do_mutate_fasta'] or operations['do_agadir'] or operations['do_foldx_repair'] \
                or operations['do_foldx_buildmodel'] or operations['do_foldx_stability'] \
                or operations['do_foldx_analysecomplex']:
            Scheduler.start(operations, path_input, input_pdb_list, input_fasta_list, list_of_mutant_aa)
        else:
            print("All of the operations are set to FALSE in MutateCompute_Options.txt file - hence nothing to do")



# cnx is the mysql connector (see top of script)
# cnx.close()
