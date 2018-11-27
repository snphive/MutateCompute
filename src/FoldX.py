#!/usr/bin/env python3
"""
Class for launching FoldX algorithms (BuildModel, AnalyseComplex, Repair). Includes code for writing FoldX's runscript file
and processing output files.

(Includes a nested Enum of FoldX-specific strings located at the end of the class.)

Note: FoldX output files are written to the working directory which must also contain certain config files and pdbs. Thus these
output & input files are located together - in output_data/build_model & output_data/analyse_complex dirs.

Note: "fxmutant" names refers specifically to the FoldX naming convention for mutants: wt_aa|chain|position|mutant_aa.
"""
import subprocess
import os
import glob
from src.enums.Str import Str
from src.enums.Paths import Paths
from src.enums.UNPW import Server
from src.tools.GeneralUtilityMethods import GUM
from src.Cluster import Cluster
import mysql.connector
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class FoldX(object):

    def write_runscript_file(self, path_runscript: str, pdbs: str, conditions: dict, action: str, num_of_runs=3,
                             show_sequence_detail=False, print_networks=False, calculate_stability=False):
        """
        The runscript.txt is an input file for FoldX indicating which pdbs to analyse & which programs to run on them.
        :param path_runscript: Absolute path for runscript.txt file being written.
        :param pdbs: pdb(s) (including .pdb extension) inputs for FoldX.
        :param conditions: temperature (Kelvin), pH and ionic strength (Molar) of solute.
        :param action: FoldX computation to be performed (e.g. BuildModel, AnalyseComplex, etc).
        :param num_of_runs: Number of times FoldX runs the pdb through algorithm. Lower to 2 or 1 for speed/disk space.
        :param show_sequence_detail: True will provide extra information in output. 
        :param print_networks: True to output network data on .. ?
        :param calculate_stability: True include a stability calculation.
        :return: The string representation of the runscript file that is written.
        """
        runscript = []
        runscript.append('<TITLE>' + self.Strs.FXRUNSCRPT.value + Str.SEMICO_NL.value)
        runscript.append('<JOBSTART>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<PDBS>' + pdbs + Str.SEMICO_NL.value)
        runscript.append('<BATCH>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<COMMANDS>' + self.Strs.FXCMDFILE.value + Str.SEMICO_NL.value)
        if show_sequence_detail:
            runscript.append('<SequenceDetail>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append(action + Str.SEMICO_NL.value)
        if print_networks:
            runscript.append('<PrintNetworks>' + Str.HASH.value + Str.SEMICO_NL.value)
        if calculate_stability:
            runscript.append('<Stability>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<END>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<OPTIONS>' + self.Strs.FXOPTFILE.value + Str.SEMICO_NL.value)
        runscript.append('<Temperature>' + str(conditions[Str.COND_TEMP.value]) + Str.SEMICO_NL.value)
        runscript.append('<IonStrength>' + str(conditions[Str.COND_ION_STRGTH.value]) + Str.SEMICO_NL.value)
        runscript.append('<ph>' + str(conditions[Str.COND_PH.value]) + Str.SEMICO_NL.value)
        runscript.append('<moveNeighbours>' + Str.TRUE_lc.value + Str.SEMICO_NL.value)
        runscript.append('<VdWDesign>2' + Str.SEMICO_NL.value)
        runscript.append('<numberOfRuns>' + str(num_of_runs) + Str.SEMICO_NL.value)
        runscript.append('<OutPDB>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<END>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<JOBEND>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<ENDFILE>' + Str.HASH.value + Str.SEMICO_NL.value)
        path_runscript_file = os.path.join(GUM.os_makedirs(path_runscript), self.Strs.runscrpt_txt.value)
        with open(path_runscript_file, 'w') as f:
            runscript_str = ''.join(runscript)
            f.write(runscript_str)
        return runscript_str

    def write_fxmutantnames(self, amino_acids: list, pdbname_chain_startpos_fasta: dict):
        """
        :param amino_acids: Amino acids that (every residue in) the protein will be mutated to.
        :param pdbname_chain_startpos_fasta: pdbname_underscore_chain and fasta sequence as key-value pair.
        :return: List of fx_style mutant names: concatenation of wtaa, chain, position and mutantaa.
        """
        fxmutantnames = []
        for pdbname_chain_startpos, fasta_sequence in pdbname_chain_startpos_fasta.items():
            chain = pdbname_chain_startpos.split('_')[-2]
            startpos = pdbname_chain_startpos.split('_')[-1]
            for index, wt_aa in enumerate(fasta_sequence):
                position = index + int(startpos)
                for mutant_aa in amino_acids:
                    fxmutantnames.append(wt_aa + chain + str(position) + mutant_aa)
        return fxmutantnames

    def remove_repaired_pdbfiles(self, path_output_bm_pdb_fxmutant_dir):
        pdbname = path_output_bm_pdb_fxmutant_dir.split('/')[-2]
        path_pdbfiles = glob.glob(path_output_bm_pdb_fxmutant_dir + pdbname + FoldX.Strs.UNDRSCR1_.value + '*' + Str.PDBEXT.value)
        for path_pdbfile in path_pdbfiles:
            GUM.linux_remove_file(path_pdbfile)

    def _get_num_of_repaired_pdbfiles(self, path_output_pdb_fxmutant_dir: str):
        """
        Finds pdb files in the specified mutant folder, with pdbname followed by "_1_". Hence does not include those for
        the wild-type pdb files.
        These pdb files (named with the _1_ suffix) are generated by FoldX's BuildModel. They are subsequently copied over to
        AnalyseComplex output folders as they are required as inputs to FoldX's AnalyseComplex computation.
        :param path_output_pdb_fxmutant_dir:
        :return: The number of repaired pdb files.
        """
        path_output_pdb_fxmutant_dir = path_output_pdb_fxmutant_dir.split('/')
        pdbname = path_output_pdb_fxmutant_dir[-2]
        path_output_pdb_fxmutant_dir = ''.join(path_output_pdb_fxmutant_dir)
        filename = pdbname + FoldX.Strs.UNDRSCR1_.value + '*' + Str.PDBEXT.value
        repaired_pdbs = glob.glob(os.path.join(path_output_pdb_fxmutant_dir, filename))
        return len(repaired_pdbs)

    class BuildModel(object):

        def __init__(self, cond: dict):
            """
            :param cond: temperature, pH, ionic strength of solute, as key-value pairs.
            """
            self.conditions = cond

        def mutate_protein_structure(self, path_pdbfile: str, amino_acids: list):
            """
            Mutate all amino acids in this pdb to all listed amino acids using FoldX BuildModel.
            Note: individual_list needs to go in the same dir as runscript.txt.
            Note: "fxmutant" names refers specifically to the FoldX naming convention for mutants: wt_aa|chain|position|mutant_aa.
            :param path_pdbfile: Absolute path to pdbfile.
            :param amino_acids: Amino acids that proteins will be mutated to.
            """
            pdbfile = path_pdbfile.split('/')[-1]
            pdbname = pdbfile.split('.')[0]
            pdbname_chain_startpos_fasta_dict = GUM.extract_pdbname_chain_startpos_fasta_from_pdbs(path_pdbfile)
            fx = FoldX()
            fxmutantnames = fx.write_fxmutantnames(amino_acids, pdbname_chain_startpos_fasta_dict)
            action = fx.Strs.BM.value + Str.HASH.value + ',' + fx.Strs.indiv_list_txt.value
            path_runscript_dir = os.path.join(Paths.CONFIG_BMRUNSCRIPT, pdbname)
            fx.write_runscript_file(path_runscript_dir, pdbfile, self.conditions, action)
            path_runscript_file = os.path.join(path_runscript_dir, fx.Strs.runscrpt_txt.value)
            has_one_chain_only = GUM.get_num_of_chains(path_pdbfile) == 1

            for fxmutantname in fxmutantnames:
                path_output_bm_pdb_fxmutant_dir = GUM.os_makedirs(Paths.OUTPUT, Paths.DIR_BM.value, pdbname, fxmutantname)
                if self._has_already_generated_dif_bm_fxoutfile(path_output_bm_pdb_fxmutant_dir):
                    print('BuildModel output data (Dif_BuildModel_RepairPDB_X.fxout) for this pdb (mutant) already exists in '
                          'output folder.')
                    continue
                GUM.linux_copy_all_files_in_dir(path_src_dir=Paths.CONFIG_FX, path_dst_dir=path_output_bm_pdb_fxmutant_dir,
                                                files_only=True)
                self._write_individual_list_for_mutant(path_dst_dir=path_output_bm_pdb_fxmutant_dir)
                os.chdir(path_output_bm_pdb_fxmutant_dir)
                path_files_to_copy = [path_runscript_file, path_pdbfile]
                GUM.linux_copy_specified_files(path_src_files=path_files_to_copy, path_dst_dir=path_output_bm_pdb_fxmutant_dir)
                if GUM.using_cluster():
                    path_jobq = GUM.os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fxmutantname)
                    Cluster.write_job_q_bash(jobname=Paths.PREFIX_FX_BM.value + fxmutantname, path_job_q_dir=path_jobq,
                                             using_runscript=True, path_runscript_dir=path_runscript_dir)
                    # path_jobq_file = os.path.join(path_jobq, Str.JOBQ.value)
                    # GUM.linux_copy_specified_files(path_jobq_file, path_dst_dir=path_output_pdbname_mutant)
                    if os.path.exists(path_runscript_file):
                        Cluster.run_job_q(path_job_q_dir=path_jobq)
                    else:
                        raise ValueError(fx.Strs.NO_RUNSCRPT_FILE_MSG.value)
                else:
                    cmd = Str.CHMOD777.value + Str.SPCE.value + Paths.LOCAL_FOLDX_EXE.value
                    try:
                        subprocess.call(cmd, shell=True)
                    except OSError:
                        print(Str.PROBLNXCMD_MSG.value + cmd)
                    cmd = Paths.FOLDX_EXE + Str.SPCE.value + fx.Strs.DSH_RUNFILE.value + Str.SPCE.value + path_runscript_file
                    if os.path.exists(path_runscript_file):
                        subprocess.call(cmd, shell=True)
                    else:
                        raise ValueError(fx.Strs.NO_RUNSCRPT_FILE_MSG.value)
                if has_one_chain_only:
                    # BuildModel has no use for these mutant pdb files and because this pdb only has 1 chain, AnalyseComplex has
                    # no use for them either. Hence, they can already be deleted, saving disk space (crucial for cluster).
                    fx.remove_repaired_pdbfiles(path_output_bm_pdb_fxmutant_dir)

            # for fxmutantname in fxmutantnames:
            #     path_output_bm_pdb_fxmutant_dir = os.path.join(Paths.OUTPUT_BM, pdbname, fxmutantname)
            #     ddG_average = self._write_ddG_csv_file(path_output_bm_pdb_fxmutant_dir, pdbname, fxmutantname)
            #     self.write_ddG_to_DB(ddG_average)

        def _has_already_generated_dif_bm_fxoutfile(self, path_output_bm_pdb_fxmutant_dir: str):
            """
            :param path_output_bm_pdb_fxmutant_dir:
            :return: True if there are atleast as many AnalyseComplex_RepairPDB_x_1_...fxout files as expected by the number of
            repaired pdb files that would be used for the computation.
            """
            path_output_bm_pdb_fxmutant_dir = path_output_bm_pdb_fxmutant_dir.split('/')
            fxmutant_dir = path_output_bm_pdb_fxmutant_dir[-1]
            pdbname = path_output_bm_pdb_fxmutant_dir[-2]
            fx = FoldX()
            path_dif_buildmodel_fxoutfile = os.path.join(Paths.OUTPUT_BM, pdbname, fxmutant_dir, fx.Strs.DIF_BLDMDL_.value +
                                                         pdbname + fx.Strs.FXOUTEXT.value)
            return os.path.exists(path_dif_buildmodel_fxoutfile)

        def write_ddG_to_DB(self, ddG_average: float):
            """
            Write the average ddG to database - TODO
            :param ddG_average:
            :return:
            """
            connection = mysql.connector.connect(user=Server.SNPEFFECT_V5_MYSQL_UN.value,
                                                 password=Server.SNPEFFECT_V5_MYSQL_PW.value,
                                                 host=Server.SNPEFFECT_V5_MYSQL_NETADD.value)
                                                #, database='snpeffect_v5', port='3306')
            cursor = connection.cursor()
            cursor.execute("SHOW DATABASES")
            for x in cursor:
                print(x)
            # cursor.execute("CREATE DATABASE SnpEffect_v5.0")
            # cursor.execute(
            #     "CREATE TABLE testingDBConnection ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")

        def _write_ddG_csv_file(self, path_output_pdb_fxmutant_dir: str, pdbname: str, fxmutantname: str):
            """
            :param path_output_pdb_fxmutant_dir:
            :param pdbname:
            :param fxmutantname:
            :return:
            """
            ddG_average = 0.0
            ddG_1 = ''
            ddG_2 = ''
            ddG_3 = ''
            fx = FoldX()
            path_output_dif_file = os.path.join(path_output_pdb_fxmutant_dir, fx.Strs.DIF_FXOUTFILE.value)
            if not os.path.exists(path_output_dif_file):
                print(path_output_dif_file + ' has not been written yet')
            else:
                self.remove_config_files(path_output_pdb_fxmutant_dir)
                GUM.linux_remove_file(os.path.join(path_output_pdb_fxmutant_dir, pdbname + Str.PDBEXT.value))
                with open(path_output_dif_file) as f:
                    dif_fxoutfile_lines = f.readlines()
                for line in dif_fxoutfile_lines:
                    if ddG_average != '':
                        break
                    if not pdbname + '_' in line:
                        continue
                    else:
                        ddG = line.split(Str.TAB.value)[1]
                        if ddG_1 == '':
                            ddG_1 = ddG
                        elif ddG_2 == '':
                            ddG_2 = ddG
                        elif ddG_3 == '':
                            ddG_3 = ddG
                ddG_average = (float(ddG_1) + float(ddG_2) + float(ddG_3)) / 3
                with open(os.path.join(path_output_pdb_fxmutant_dir, fx.Strs.DDG_CSV.value), 'w') as f:
                    f.write(pdbname + ',' + fxmutantname + ',' + str(ddG_average))
            return ddG_average

        def remove_config_files(self, path_output_pdb_fxmutant_dir: str):
            """
            Remove the copies of config files that had to be copied to the FoldX output dir, but once used serve no
            purpose and only take up disk space.
            :param path_output_pdb_fxmutant_dir: Absoluate path to output dir for each mutant housing the individual
            copies of the config files to be deleted.
            """
            fx = FoldX()
            path_runscript_txt = os.path.join(path_output_pdb_fxmutant_dir, fx.Strs.runscrpt_txt.value)
            path_rotabase_txt = os.path.join(path_output_pdb_fxmutant_dir, fx.Strs.ROTABASE_TXT.value)
            path_cmds_bm_txt = os.path.join(path_output_pdb_fxmutant_dir, fx.Strs.CMNDS_BM_TXT.value)
            path_cmds_stab_txt = os.path.join(path_output_pdb_fxmutant_dir, fx.Strs.CMNDS_STAB_TXT.value)
            path_options_bm_txt = os.path.join(path_output_pdb_fxmutant_dir, fx.Strs.OPTNS_BM_TXT.value)
            path_options_stab_txt = os.path.join(path_output_pdb_fxmutant_dir, fx.Strs.OPTNS_STAB_TXT.value)
            path_indiv_list_txt = os.path.join(path_output_pdb_fxmutant_dir, fx.Strs.indiv_list_txt.value)
            path_files_to_remove = [path_runscript_txt, path_rotabase_txt, path_cmds_bm_txt, path_cmds_stab_txt,
                                    path_options_bm_txt, path_options_stab_txt, path_indiv_list_txt]
            for path_file_to_remove in path_files_to_remove:
                if os.path.exists(path_file_to_remove):
                    GUM.linux_remove_file(path_file_to_remove)

        # Not unit tested yet.
        def _write_individual_list_for_mutant(self, path_dst_dir: str):
            """
            Writes the 'individual_list.txt' file that is used by fx/runscript to identify the protein to process.
            :param path_dst_dir: Abs path to the output dir of foldx computation for each mutant
            """
            fxmutantname = path_dst_dir.split('/')[-1]
            with open(os.path.join(path_dst_dir, FoldX().Strs.indiv_list_txt.value), 'w') as f:
                f.write(fxmutantname + ';\n')

    class AnalyseComplex(object):

        def __init__(self, cond: dict):
            """
            :param cond: temperature, pH, ionic strength of solute, as key-value pairs.
            """
            self.conditions = cond

        def calculate_complex_energies(self, path_pdbfile: str):
            """
            Calculates interaction energies for a complex of protein chains for a specified pdb and list of amino acids.
            The repaired pdbs are read from each output_data/buildmodel/<pdbname>/<fxmutantname> folder.
            :param path_pdbfile:
            """
            pdbfile = path_pdbfile.split('/')[-1]
            pdbname = pdbfile.split('.')[0]
            wt_pdbname = FoldX.Strs.WT_.value + pdbname
            fx = FoldX()
            path_output_bm_pdb_fxmutant_dirs = sorted(glob.glob(os.path.join(Paths.OUTPUT_BM, pdbname, '*')))
            action = fx.Strs.AC.value + Str.HASH.value
            # action += '0'  # According to FoldX manual (version3), write 0 for all possible interation in a complex.
            # Write AB for example if only interested in iteraction energies of chain A with chain B.
            has_one_chain_only = GUM.get_num_of_chains(path_pdbfile) == 1
            for path_output_bm_pdb_fxmutant_dir in path_output_bm_pdb_fxmutant_dirs:
                if has_one_chain_only:
                    print('Only one protein chain. No interactions to be quantified for this pdbfile.')
                    break
                fxmutant_dir = os.path.basename(path_output_bm_pdb_fxmutant_dir)
                path_output_ac_pdb_fxmutant_dir = GUM.os_makedirs(Paths.OUTPUT, Paths.DIR_AC.value, pdbname, fxmutant_dir)
                if self._has_already_generated_all_ac_fxoutfile(path_output_ac_pdb_fxmutant_dir):
                    print('AnalyseComplex_RepairPDB_X.fxout for this pdb (mutant) already exists in output folder.')
                    break
                pdbs_to_analyse = [pdbname + FoldX.Strs._1_012_SUFFIX_PDBS.value[0], ',', wt_pdbname +
                                   FoldX.Strs._1_012_SUFFIX_PDBS.value[0]]
                num_of_repaired_pdbs = fx._get_num_of_repaired_pdbfiles(path_output_bm_pdb_fxmutant_dir)
                for i in range(1, num_of_repaired_pdbs):
                    pdbs_to_analyse.append(',')
                    pdbs_to_analyse.append(pdbname)
                    pdbs_to_analyse.append(FoldX.Strs._1_012_SUFFIX_PDBS.value[i])
                    pdbs_to_analyse.append(',')
                    pdbs_to_analyse.append(wt_pdbname)
                    pdbs_to_analyse.append(FoldX.Strs._1_012_SUFFIX_PDBS.value[i])
                pdbs_to_analyse = ''.join(pdbs_to_analyse)
                fx.write_runscript_file(path_output_ac_pdb_fxmutant_dir, pdbs_to_analyse, self.conditions, action)
                path_runscript_file = os.path.join(path_output_ac_pdb_fxmutant_dir, fx.Strs.runscrpt_txt.value)
                GUM.linux_copy_all_files_in_dir(Paths.CONFIG_FX, path_output_ac_pdb_fxmutant_dir, files_only=True)
                os.chdir(path_output_ac_pdb_fxmutant_dir)
                files_to_copy = pdbs_to_analyse.split(',')
                path_files_to_copy = []
                for file_to_copy in files_to_copy:
                    path_files_to_copy.append(os.path.join(path_output_bm_pdb_fxmutant_dir, file_to_copy))
                GUM.linux_copy_specified_files(path_files_to_copy, path_dst_dir=path_output_ac_pdb_fxmutant_dir)

                if GUM.using_cluster():
                    Cluster.write_job_q_bash(jobname=Paths.PREFIX_FX_AC.value + fxmutant_dir,
                                             path_job_q_dir=path_output_ac_pdb_fxmutant_dir, using_runscript=True,
                                             path_runscript_dir=path_output_ac_pdb_fxmutant_dir)
                    if os.path.exists(path_runscript_file):
                        Cluster.run_job_q(path_job_q_dir=path_output_ac_pdb_fxmutant_dir)
                    else:
                        raise ValueError(fx.Strs.NO_RUNSCRPT_FILE_MSG.value)
                else:
                    cmd = Str.CHMOD777.value + Str.SPCE.value + Paths.LOCAL_FOLDX_EXE.value
                    try:
                        subprocess.call(cmd, shell=True)
                    except OSError:
                        print(Str.PROBLNXCMD_MSG.value + cmd)
                    cmd = Paths.FOLDX_EXE + Str.SPCE.value + fx.Strs.DSH_RUNFILE.value + Str.SPCE.value + path_runscript_file
                    if os.path.exists(path_runscript_file):
                        subprocess.call(cmd, shell=True)
                    else:
                        raise ValueError(fx.Strs.NO_RUNSCRPT_FILE_MSG.value)
                if self._has_already_generated_all_ac_fxoutfile(path_output_ac_pdb_fxmutant_dir):
                    fx.remove_repaired_pdbfiles(path_output_bm_pdb_fxmutant_dir)

        def _has_already_generated_all_ac_fxoutfile(self, path_output_ac_pdb_fxmutant_dir: str):
            """
            :param path_output_ac_pdb_fxmutant_dir:
            :return: True if there are atleast as many AnalyseComplex_RepairPDB_x_1_...fxout files as expected by the number of
            repaired pdb files that would be used for the computation.
            """
            path_output_ac_pdb_fxmutant_dir = path_output_ac_pdb_fxmutant_dir.split('/')
            fxmutant_dir = path_output_ac_pdb_fxmutant_dir[-1]
            pdbname = path_output_ac_pdb_fxmutant_dir[-2]
            path_output_ac_pdb_fxmutant_dir = ''.join(path_output_ac_pdb_fxmutant_dir)
            fx = FoldX()
            num_of_repaired_pdbfiles = fx._get_num_of_repaired_pdbfiles(path_output_ac_pdb_fxmutant_dir)
            num_of_analysecomplex_fxoutfiles = len(glob.glob(os.path.join(Paths.OUTPUT_AC, pdbname, fxmutant_dir,
                                                                          fx.Strs.ANLYSCMPLX_.value + pdbname +
                                                                          fx.Strs.UNDRSCR1_.value + '*' +
                                                                          fx.Strs.FXOUTEXT.value)))
            return num_of_repaired_pdbfiles <= num_of_analysecomplex_fxoutfiles

    class Repair(object):

        def __init__(self, cond: dict):
            """
            :param cond: temperature, pH, ionic strength of solute, as key-value pairs.
            """
            self.conditions = cond

        def do_repair(self, input_pdb):
            print('TODO')

    from enum import Enum

    class Strs(Enum):

        AC = '<AnalyseComplex>'
        BM = '<BuildModel>'
        indiv_list_txt = 'individual_list.txt'
        FXRUNSCRPT = 'FOLDX_runscript'
        FXCMDFILE = 'FOLDX_commandfile'
        FXOPTFILE = 'FOLDX_optionfile'
        runscrpt_txt = Str.runscrpt_txt.value
        NO_RUNSCRPT_FILE_MSG = 'No runscript file was found'
        DSH_RUNFILE = '-runfile'
        DIF_FXOUTFILE = 'Dif_BuildModel_RepairPDBtest_1.fxout'
        DDG_CSV = 'ddG' + Str.CSVEXT.value
        FXOUTEXT = Str.FXOUTEXT.value
        ROTABASE_TXT = 'rotabase.txt'
        CMNDS_BM_TXT = 'commands_buildmodel.txt'
        CMNDS_STAB_TXT = 'commands_stability.txt'
        OPTNS_BM_TXT = 'options_buildmodel.txt'
        OPTNS_STAB_TXT = 'options_stability.txt'
        _1_012_SUFFIX_PDBS = ['_1_0' + Str.PDBEXT.value, '_1_1' + Str.PDBEXT.value, '_1_2' + Str.PDBEXT.value]
        UNDRSCR1_ = '_1_'
        WT_ = 'WT_'
        ANLYSCMPLX_ = 'AnalyseComplex_'
        DIF_BLDMDL_ = 'Dif_BuildModel_'


# pydevd.stoptrace()
