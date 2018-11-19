#!/usr/bin/env python3
"""
Class for launching FoldX algorithms (BuildModel, AnalyseComplex, Stability). Includes code for writing FoldX'' runscript
file and processing output files.
Includes a nested Enum of FoldX-specific strings located at the end of the class.
"""
import subprocess
import os
from src.Str import Str
from src.Paths import Paths
from src.GeneralUtilityMethods import GUM
from src.Cluster import Cluster
import mysql.connector
from src.Conditions import Cond

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class FoldX(object):

    def __init__(self):
        self.BuildModel(Cond.INCELL_MAML_FX.value).write_ddG_to_DB()

    # 30.07.18 Redesigned the directory structure such that runscripts will go in configuration/foldx/ & maybe another
    # level such as analyse_complex or build_model or stability etc.
    def write_runscript_file(self, path_runscript: str, pdbs: str, conditions: dict, action: str,
                             show_sequence_detail=False, print_networks=False, calculate_stability=False):
        """
        The runscript.txt is an input file for FoldX indicating which pdbs to analyse & which programs to run on them.
        :param path_runscript: Absolute path for runscript.txt file being written.
        :param pdbs: pdb(s) (including .pdb extension) inputs for FoldX.
        :param conditions: temperature (Kelvin), pH and ionic strength (Molar) of solute.
        :param action: FoldX computation to be performed (e.g. BuildModel, AnalyseComplex, etc).
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
        runscript.append('<numberOfRuns>3' + Str.SEMICO_NL.value)
        runscript.append('<OutPDB>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<END>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<JOBEND>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<ENDFILE>' + Str.HASH.value + Str.SEMICO_NL.value)
        path_runscript_file = os.path.join(GUM._os_makedirs(path_runscript), self.Strs.runscrpt_txt.value)
        with open(path_runscript_file, 'w') as f:
            runscript_str = ''.join(runscript)
            f.write(runscript_str)
        return runscript_str

    class Repair(object):

        def __init__(self, cond: dict):
            """
            :param cond: temperature, pH, ionic strength of solute, as key-value pairs.
            """
            self.conditions = cond

        def do_repair(self, input_pdb):
            print('do repair method')

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
            :param path_pdbfile: Absolute path to pdbfile.
            :param amino_acids: Amino acids that proteins will be mutated to.
            """
            pdbfile = path_pdbfile.split('/')[-1]
            pdbname = pdbfile.split('.')[0]
            pdbname_chain_startpos_fasta_dict = GUM.extract_pdbname_chain_startpos_fasta_from_pdbs(path_pdbfile)
            fx_mutant_name_list = FoldX._make_fx_mutant_name_list(amino_acids, pdbname_chain_startpos_fasta_dict)
            action = FoldX().Strs.BM.value + Str.HASH.value + ',' + FoldX().Strs.indiv_list_txt.value
            path_runscript_dir = os.path.join(Paths.CONFIG_BMRUNSCRIPT, pdbname)
            FoldX().write_runscript_file(path_runscript_dir, pdbfile, self.conditions, action)
            path_runscript_file = os.path.join(path_runscript_dir, FoldX().Strs.runscrpt_txt.value)

            for fx_mutant_name in fx_mutant_name_list:
                path_output_pdbname_mutant = FoldX._make_output_dir_and_copy_fxconfig_files_in(Paths.OUTPUT, pdbname,
                                                                                              fx_mutant_name)
                self._write_individual_list_for_mutant(path_output_pdbname_mutant)
                os.chdir(path_output_pdbname_mutant)
                path_files_to_copy = [path_runscript_file, path_pdbfile]
                GUM.linux_copy_specified_files(path_files_to_copy, path_dst_dir=path_output_pdbname_mutant)
                if GUM.using_cluster():
                    path_jobq = GUM._os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fx_mutant_name)
                    Cluster.write_job_q_bash(jobname=Paths.PREFIX_FX_BM.value + fx_mutant_name,
                                             path_job_q_dir=path_jobq, path_dst_dir=path_output_pdbname_mutant,
                                             path_runscript_dir=path_runscript_dir, using_runscript=True)
                    # path_jobq_file = os.path.join(path_jobq, Str.JOBQ.value)
                    # GUM.linux_copy_specified_files(path_jobq_file, path_dst_dir=path_output_pdbname_mutant)
                    if os.path.exists(path_runscript_file):
                        Cluster.run_job_q(path_job_q_dir=path_jobq)
                    else:
                        raise ValueError(FoldX().Strs.NO_RUNSCRPT_FILE_MSG.value)
                else:
                    cmd = 'chmod 777' + Str.SPCE.value + Paths.LOCAL_FOLDX_EXE.value
                    try:
                        subprocess.call(cmd, shell=True)
                    except OSError:
                        print(Str.PROBLNXCMD_MSG.value + cmd)
                    cmd = Paths.FOLDX_EXE + Str.SPCE.value + FoldX().Strs.DSH_RUNFILE.value + Str.SPCE.value + \
                          path_runscript_file
                    if os.path.exists(path_runscript_file):
                        subprocess.call(cmd, shell=True)
                    else:
                        raise ValueError(FoldX().Strs.NO_RUNSCRPT_FILE_MSG.value)

            for fx_mutant_name in fx_mutant_name_list:
                path_output_pdbname_mutant = os.path.join(Paths.OUTPUT, pdbname, fx_mutant_name)
                self._write_ddG_csv_file(path_output_pdbname_mutant, pdbname, fx_mutant_name)

        def write_ddG_to_DB(self):
            connection = mysql.connector.connect(user='snpeffect_v5',
                                                 password='R34WKKGR',
                                                 host='127.0.0.1')
                                                #, database='snpeffect_v5', port='3306')
            cursor = connection.cursor()
            cursor.execute("SHOW DATABASES")
            for x in cursor:
                print(x)
            # cursor.execute("CREATE DATABASE SnpEffect_v5.0")
            # cursor.execute(
            #     "CREATE TABLE testingDBConnection ( id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT, title TEXT NOT NULL )")

        def _write_ddG_csv_file(self, path_output_pdbname_mutant: str, pdbname: str, fx_mutant_name: str):
            """
            :param path_output_pdbname_mutant:
            :param pdbname:
            :param fx_mutant_name:
            :return:
            """
            ddG_average = ''
            ddG_1 = ''
            ddG_2 = ''
            ddG_3 = ''
            path_output_dif_file = os.path.join(path_output_pdbname_mutant, FoldX().Strs.DIF_FXOUTFILE.value)
            if not os.path.exists(path_output_dif_file):
                print(path_output_dif_file + ' has not be written yet')
            else:
                self.remove_config_files(path_output_pdbname_mutant)
                GUM.linux_remove_file(os.path.join(path_output_pdbname_mutant, pdbname + '.pdb'))
                with open(path_output_dif_file) as f:
                    dif_fxoutfile_lines = f.readlines()
                for line in dif_fxoutfile_lines:
                    if ddG_average != '':
                        break
                    if not pdbname + '_' in line:
                        continue
                    else:
                        ddG = float(line.split(Str.TAB.value)[1])
                        if ddG_1 == '':
                            ddG_1 = ddG
                        elif ddG_2 == '':
                            ddG_2 = ddG
                        elif ddG_3 == '':
                            ddG_3 = ddG
                ddG_average = (ddG_1 + ddG_2 + ddG_3) / 3
                with open(os.path.join(path_output_pdbname_mutant, FoldX().Strs.DDG_CSV.value), 'w') as f:
                    f.write(pdbname + ',' + fx_mutant_name + ',' + str(ddG_average))
            return ddG_average

        def remove_config_files(self, path_output_pdbname_mutant: str):
            """
            Remove the copies of config files that had to be copied to the FoldX output dir, but once used serve no
            purpose and only take up disk space.
            :param path_output_pdbname_mutant: Absoluate path to output dir for each mutant housing the individual
            copies of the config files to be deleted.
            """
            path_runscript_txt = os.path.join(path_output_pdbname_mutant, FoldX().Strs.runscrpt_txt.value)
            path_rotabase_txt = os.path.join(path_output_pdbname_mutant, FoldX().Strs.ROTABASE_TXT.value)
            path_cmds_bm_txt = os.path.join(path_output_pdbname_mutant, FoldX().Strs.CMNDS_BM_TXT.value)
            path_cmds_stab_txt = os.path.join(path_output_pdbname_mutant, FoldX().Strs.CMNDS_STAB_TXT.value)
            path_options_bm_txt = os.path.join(path_output_pdbname_mutant, FoldX().Strs.OPTNS_BM_TXT.value)
            path_options_stab_txt = os.path.join(path_output_pdbname_mutant, FoldX().Strs.OPTNS_STAB_TXT.value)
            path_indiv_list_txt = os.path.join(path_output_pdbname_mutant, FoldX().Strs.indiv_list_txt.value)
            path_files_to_remove = [path_runscript_txt, path_rotabase_txt, path_cmds_bm_txt, path_cmds_stab_txt,
                                    path_options_bm_txt, path_options_stab_txt, path_indiv_list_txt]
            for path_file_to_remove in path_files_to_remove:
                if os.path.exists(path_file_to_remove):
                    GUM.linux_remove_file(path_file_to_remove)



        # Not tested yet.
        def _write_individual_list_for_mutant(self, path_dest_mutant_dir: str):
            """
            Writes the 'individual_list.txt' file that is used by fx/runscript to identify the protein to process.
            :param path_dest_mutant_dir: Abs path to the output dir of foldx computation for each mutant 
            """
            fx_mutant_name = path_dest_mutant_dir.split('/')[-1]
            with open(os.path.join(path_dest_mutant_dir, FoldX().Strs.indiv_list_txt.value), 'w') as f:
                f.write(fx_mutant_name + ';\n')

    class Stability(object):

        def __init__(self, cond: dict):
            """
            :param cond: temperature, pH, ionic strength of solute, as key-value pairs.
            """
            self.conditions = cond

        def calculate_stability(self, pdb):
            print('TODO')

    class AnalyseComplex(object):

        def __init__(self, cond: dict):
            """
            :param cond: temperature, pH, ionic strength of solute, as key-value pairs.
            """
            self.conditions = cond

        def calculate_complex_energies(self, path_pdbfile: str, amino_acids: list):
            pdbfile = path_pdbfile.split('/')[-1]
            pdbname = pdbfile.split('.')[0]
            pdbname_chain_startpos_fasta_dict = GUM.extract_pdbname_chain_startpos_fasta_from_pdbs(path_pdbfile)
            fx_mutant_name_list = FoldX._make_fx_mutant_name_list(amino_acids, pdbname_chain_startpos_fasta_dict)

            path_runscript_dir = os.path.join(Paths.CONFIG_ACRUNSCRIPT, pdbname)
            wt_pdbname = FoldX.Strs.WT_.value + pdbname
            pdbs_to_analyse = pdbname + FoldX.Strs._1_0_1_2_SUFFIX_PDBS[0] + ',' + \
                              pdbname + FoldX.Strs._1_0_1_2_SUFFIX_PDBS[1] + ',' + \
                              pdbname + FoldX.Strs._1_0_1_2_SUFFIX_PDBS[2] + ',' + \
                              wt_pdbname + FoldX.Strs._1_0_1_2_SUFFIX_PDBS[0] + ',' + \
                              wt_pdbname + FoldX.Strs._1_0_1_2_SUFFIX_PDBS[1] + ',' + \
                              wt_pdbname + FoldX.Strs._1_0_1_2_SUFFIX_PDBS[2]
            action = FoldX().Strs.AC.value + Str.HASH.value
            FoldX().write_runscript_file(path_runscript_dir, pdbs_to_analyse, self.conditions, action)
            path_runscript_file = os.path.join(path_runscript_dir, FoldX().Strs.runscrpt_txt.value)

            for fx_mutant_name in fx_mutant_name_list:
                path_output_pdbname_mutant = FoldX._make_output_dir_and_copy_fxconfig_files_in(Paths.OUTPUT, pdbname,
                                                                                              fx_mutant_name)
                os.chdir(path_output_pdbname_mutant)
                path_files_to_copy = [path_runscript_file, path_pdbfile]
                GUM.linux_copy_specified_files(path_files_to_copy, path_dst_dir=path_output_pdbname_mutant)
                if GUM.using_cluster():
                    path_jobq = GUM._os_makedirs(Paths.CONFIG_AC_JOBQ, pdbname, fx_mutant_name)
                    Cluster.write_job_q_bash(jobname=Paths.PREFIX_FX_BM.value + fx_mutant_name,
                                             path_job_q_dir=path_jobq, path_dst_dir=path_output_pdbname_mutant,
                                             path_runscript_dir=path_runscript_dir, using_runscript=True)
                    # path_jobq_file = os.path.join(path_jobq, Str.JOBQ.value)
                    # GUM.linux_copy_specified_files(path_jobq_file, path_dst_dir=path_output_pdbname_mutant)
                    if os.path.exists(path_runscript_file):
                        Cluster.run_job_q(path_job_q_dir=path_jobq)
                    else:
                        raise ValueError(FoldX().Strs.NO_RUNSCRPT_FILE_MSG.value)
                else:
                    cmd = 'chmod 777' + Str.SPCE.value + Paths.LOCAL_FOLDX_EXE.value
                    try:
                        subprocess.call(cmd, shell=True)
                    except OSError:
                        print(Str.PROBLNXCMD_MSG.value + cmd)
                    cmd = Paths.FOLDX_EXE + Str.SPCE.value + FoldX().Strs.DSH_RUNFILE.value + Str.SPCE.value + \
                          path_runscript_file
                    if os.path.exists(path_runscript_file):
                        subprocess.call(cmd, shell=True)
                    else:
                        raise ValueError(FoldX().Strs.NO_RUNSCRPT_FILE_MSG.value)

    def _make_output_dir_and_copy_fxconfig_files_in(self, path_output: str, pdbname: str, fx_mutant_name: str):
        """
        :param path_output:
        :param pdbname:
        :param fx_mutant_name:
        :return:
        """
        path_output_pdbname_mutant = GUM._os_makedirs(path_output, Paths.DIR_BM.value, pdbname, fx_mutant_name)
        GUM.linux_copy_all_files_in_dir(Paths.CONFIG_FX, path_output_pdbname_mutant, files_only=True)
        return path_output_pdbname_mutant

    def _make_fx_mutant_name_list(self, amino_acids: list, pdbname_chain_startpos_fasta: dict):
        """
        :param amino_acids: Amino acids that (every residue in) the protein will be mutated to.
        :param pdbname_chain_startpos_fasta: pdbname_underscore_chain and fasta sequence as key-value pair.
        :return: List of fx_style mutant names: concatenation of wtaa, chain, position and mutantaa.
        """
        fx_mutant_name_list = []
        for pdbname_chain_startpos, fasta_sequence in pdbname_chain_startpos_fasta.items():
            chain = pdbname_chain_startpos.split('_')[-2]
            startpos = pdbname_chain_startpos.split('_')[-1]
            for index, wt_aa in enumerate(fasta_sequence):
                position = index + int(startpos)
                for mutant_aa in amino_acids:
                    fx_mutant_name_list.append(wt_aa + chain + str(position) + mutant_aa)
        return fx_mutant_name_list

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
        DDG_CSV = 'ddG.csv'
        ROTABASE_TXT = 'rotabase.txt'
        CMNDS_BM_TXT = 'commands_buildmodel.txt'
        CMNDS_STAB_TXT = 'commands_stability.txt'
        OPTNS_BM_TXT = 'options_buildmodel.txt'
        OPTNS_STAB_TXT = 'options_stability.txt'
        _1_0_1_2_SUFFIX_PDBS = ['_1_0' + Str.PDBEXT.value, '_1_1' + Str.PDBEXT.value, '_1_2' + Str.PDBEXT.value]
        WT_ = 'WT_'



