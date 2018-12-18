#!/usr/bin/env python3
"""
Class for launching FoldX algorithms (BuildModel, AnalyseComplex, Repair). Includes code for writing FoldX's runscript file
and processing output files.

(Includes a nested Enum of FoldX-specific strings located at the end of the class.)

Note: FoldX output files are written to the working directory which must also contain certain config files and pdbs (i.e. the
input files). As such, this means that there is only one mutation folder for both the input and the output files.
Currently, this is only the case with FoldX computations, in all other parts of the MutateCompute project, input files and
output files are read and written from separate paths (input_data/ and output_data/).

Note: "fxmutantname" refers specifically to the FoldX naming convention for mutants: wt_aa|chain|position|mutant_aa.
Hence, the FoldX format for an Alanine to Valine substitution, at position 100, on protein chain B is "AB100V".
"""
import subprocess
import os
import glob
from src.enums.Str import Str
from src.enums.Paths import Paths
from src.tools.GeneralUtilityMethods import GUM
from src.Cluster import Cluster
import pydevd
pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class FoldX(object):

    def write_runscript_file(self, path_runscript: str, pdbs: str, conditions: dict, action: str,
                             num_of_runs=Str.DEFAULT_NUM_FX_RUNS_1.value, show_sequence_detail=False, print_networks=False,
                             calculate_stability=False):
        """
        The runscript.txt is an input file for FoldX indicating which pdbs to analyse & which programs to run on them.
        :param path_runscript: Absolute path for runscript.txt file being written.
        :param pdbs: pdb(s) (including .pdb extension) inputs for FoldX.
        :param conditions: temperature (Kelvin), pH and ionic strength (Molar) of solute.
        :param action: FoldX computation to be performed (e.g. BuildModel, AnalyseComplex, etc).
        :param num_of_runs: Number of times FoldX runs the pdb through algorithm. Default in FoldXv3 manual is 1. (Note: It is
        not possible to access FoldX.Strs enum in method argument so resorted to storing DEFAULT_NUM_FX_RUNS in Str enum.)
        :param show_sequence_detail: True will provide extra information in output. 
        :param print_networks: True to output network data on .. ?
        :param calculate_stability: True to include a stability calculation.
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
        runscript.append('<VdWDesign>' + str(self.Strs.DEFAULT_VDW_DSGN.value) + Str.SEMICO_NL.value)
        runscript.append('<numberOfRuns>' + str(num_of_runs) + Str.SEMICO_NL.value)
        runscript.append('<OutPDB>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<END>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<JOBEND>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<ENDFILE>' + Str.HASH.value + Str.SEMICO_NL.value)
        path_runscript_file = os.path.join(GUM.os_makedirs(path_runscript), self.Strs.RNSCRPT_TXT.value)
        with open(path_runscript_file, 'w') as f:
            runscript_str = ''.join(runscript)
            f.write(runscript_str)
        return runscript_str

    def make_fxmutantnames(self, amino_acids: list, pdbname_chain_startpos_fasta: dict):
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

    def rm_pdbfiles(self, path_output_ac_or_bm_pdb_fxmutant_dir: str):
        """
        Remove all pdb files from specified directory directory.
        :param path_output_ac_or_bm_pdb_fxmutant_dir: Absolute path for pdb files to be deleted.
        """
        path_pdbfiles = glob.glob(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, '*' + Str.PDBEXT.value))
        for path_pdbfile in path_pdbfiles:
            GUM.linux_remove_file(path_pdbfile)

    def rm_logfiles(self, path_output_ac_or_bm_pdb_fxmutant_dir: str):
        """
        Remove all e.1234567 and o.1234567 cluster log files from specified output directory.
        :param path_output_ac_or_bm_pdb_fxmutant_dir: Absolute path for log files to be deleted.
        """
        fxmutantname = os.path.basename(path_output_ac_or_bm_pdb_fxmutant_dir)
        path_output_ac_or_bm_pdb_fxmutant_ologfiles = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, '*_' + fxmutantname +
                                                                   Str.CLSTR_OUT_LOGEXT.value + '*')
        GUM.linux_remove_file(path_output_ac_or_bm_pdb_fxmutant_ologfiles)

        path_output_ac_or_bm_pdb_fxmutant_elogfiles = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, '*_' + fxmutantname +
                                                                   Str.CLSTR_ERR_LOGEXT.value + '*')
        path_output_ac_or_bm_pdb_fxmutant_elogfiles = glob.glob(path_output_ac_or_bm_pdb_fxmutant_elogfiles)

        for path_output_ac_or_bm_pdb_fxmutant_elogfile in path_output_ac_or_bm_pdb_fxmutant_elogfiles:
            if os.stat(path_output_ac_or_bm_pdb_fxmutant_elogfile).st_size != 0:
                with open(path_output_ac_or_bm_pdb_fxmutant_elogfile) as f:
                    path_output_ac_or_bm_pdb_fxmutant_elogfile = path_output_ac_or_bm_pdb_fxmutant_elogfiles.split('/')
                    pdbname = path_output_ac_or_bm_pdb_fxmutant_elogfile[-3]
                    fxmutantname = path_output_ac_or_bm_pdb_fxmutant_elogfile[-2]
                    print('This error logfile is not 0 bytes in size (' + pdbname + ':' + fxmutantname + ')')
                    print('The error logfile reads: ' + f.read())
            else:
                GUM.linux_remove_file(path_output_ac_or_bm_pdb_fxmutant_elogfile)

    def rm_unnecessary_fxoutfiles(self, path_output_ac_or_bm_pdb_fxmutant_dir: str):
        """
        Remove all unnecessary FoldX output files.
        For BuildModel, this includes Raw_BuildModel_pdbname.fxout, Dif_BuildModel_pdbname.fxout, BuildModel_pdbname.fxout and
        PdbList_BuildModel_pdbname.fxout.
        For AnalyseComplex, this includes AnalyseComplex_pdbname.fxout, Indiv_energies_AnalyseComplex_pdbname.fxout,
        Interaction_AnalyseComplex_pdbname.fxout and Interface_Residues_AnalyseComplex_pdbname.fxout.
        :param path_output_ac_or_bm_pdb_fxmutant_dir: Absolute path for unnecessary FoldX files to be deleted.
        """
        path_output_ac_or_bm_pdb_fxmutant_dir = path_output_ac_or_bm_pdb_fxmutant_dir.split('/')
        ac_or_bm = path_output_ac_or_bm_pdb_fxmutant_dir[-3]
        pdbname = path_output_ac_or_bm_pdb_fxmutant_dir[-2]
        path_output_ac_or_bm_pdb_fxmutant_dir = '/'.join(path_output_ac_or_bm_pdb_fxmutant_dir)
        fx = FoldX()
        prefixes_of_files_to_rm = fx.Strs.BM_FILES_TO_DELETE.value if ac_or_bm == Paths.DIR_BM.value else fx.Strs.AC_FILES_TO_DELETE.value
        for prefix_of_file_to_rm in prefixes_of_files_to_rm:
            GUM.linux_remove_file(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, prefix_of_file_to_rm + '*' + pdbname + '*' +
                                               fx.Strs.FXOUTEXT.value))

    def rm_config_files(self, path_output_ac_or_bm_pdb_fxmutant_dir: str):
        """
        Remove all config files from specified output directory. This includes runscript.txt, rotabase.txt,
        commands_stability.txt, commands_buildmodel.txt, options_stability.txt, options_buildmodel.txt, individual_list.txt,
        and job.q files.
        :param path_output_ac_or_bm_pdb_fxmutant_dir: Absolute path for config files to be deleted.
        """
        fx = FoldX()
        path_runscript_txt = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, fx.Strs.RNSCRPT_TXT.value)
        path_rotabase_txt = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, fx.Strs.ROTABASE_TXT.value)
        path_cmds_bm_txt = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, fx.Strs.CMNDS_BM_TXT.value)
        path_cmds_stab_txt = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, fx.Strs.CMNDS_STAB_TXT.value)
        path_options_bm_txt = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, fx.Strs.OPTNS_BM_TXT.value)
        path_options_stab_txt = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, fx.Strs.OPTNS_STAB_TXT.value)
        path_indiv_list_txt = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, fx.Strs.INDIV_LST_TXT.value)
        path_job_q = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, fx.Strs.JOB_Q.value)
        path_files_to_remove = [path_runscript_txt, path_rotabase_txt, path_cmds_bm_txt, path_cmds_stab_txt,
                                path_options_bm_txt, path_options_stab_txt, path_indiv_list_txt, path_job_q]
        for path_file_to_remove in path_files_to_remove:
            if os.path.exists(path_file_to_remove):
                GUM.linux_remove_file(path_file_to_remove)

    def _get_num_of_repaired_pdbfiles(self, path_output_pdb_fxmutant_dir: str):
        """
        Finds pdb files in the specified mutant folder, with pdbname followed by "_1_". NOTE: This does not include those for
        the wild-type pdb files, so it should return half the total number of repaired pdbs.
        These pdb files (named with the _1_ suffix) are generated by FoldX's BuildModel. They are subsequently copied over to
        AnalyseComplex output folders where they serve as input files to the AnalyseComplex computation.
        :param path_output_pdb_fxmutant_dir: Absolute path to a pdb mutant folder.
        :return: The number of repaired pdb files.
        """
        path_output_pdb_fxmutant_dir = path_output_pdb_fxmutant_dir.split('/')
        pdbname = path_output_pdb_fxmutant_dir[-2]
        path_output_pdb_fxmutant_dir = '/'.join(path_output_pdb_fxmutant_dir)
        filename = pdbname + FoldX.Strs.UNDRSCR1.value + '*' + Str.PDBEXT.value
        repaired_pdbs = glob.glob(os.path.join(path_output_pdb_fxmutant_dir, filename))
        return len(repaired_pdbs)

    def has_already_removed_config_logs_fxoutfile(self, path_output_ac_or_bm_pdb_fxmutant_dir: str):
        """
        :param path_output_ac_or_bm_pdb_fxmutant_dir: Absolute path to the folder of interest.
        :return: True if unnecessary files are absent from specified location.
        """
        has_no_files_to_rm = True
        fx = FoldX()
        has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                                                    fx.Strs.RNSCRPT_TXT.value))
        has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                                                    fx.Strs.ROTABASE_TXT.value))
        has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                                                    fx.Strs.CMNDS_BM_TXT.value))
        has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                                                    fx.Strs.CMNDS_STAB_TXT.value))
        has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                                                    fx.Strs.OPTNS_BM_TXT.value))
        has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                                                    fx.Strs.OPTNS_STAB_TXT.value))
        has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                                                    fx.Strs.INDIV_LST_TXT.value))
        has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                                                    fx.Strs.JOB_Q.value))

        path_output_ac_or_bm_pdb_fxmutant_dir = path_output_ac_or_bm_pdb_fxmutant_dir.split('/')
        ac_or_bm = path_output_ac_or_bm_pdb_fxmutant_dir[-3]
        pdbname = path_output_ac_or_bm_pdb_fxmutant_dir[-2]
        fxmutantname = path_output_ac_or_bm_pdb_fxmutant_dir[-1]
        path_output_ac_or_bm_pdb_fxmutant_dir = '/'.join(path_output_ac_or_bm_pdb_fxmutant_dir)
        prefixes_of_files_to_rm = fx.Strs.BM_FILES_TO_DELETE.value if ac_or_bm == Paths.DIR_BM.value else \
            fx.Strs.AC_FILES_TO_DELETE.value
        for prefix_of_file_to_rm in prefixes_of_files_to_rm:
            path_fxout_files = glob.glob(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, prefix_of_file_to_rm + '*' +
                                                      pdbname + '*' + fx.Strs.FXOUTEXT.value))
            for path_fxout_file in path_fxout_files:
                has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(path_fxout_file)
        path_output_ac_or_bm_pdb_fxmutant_ologfiles = glob.glob(os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                                '*_' + fxmutantname + Str.CLSTR_OUT_LOGEXT.value + '*'))
        for path_ologfile in path_output_ac_or_bm_pdb_fxmutant_ologfiles:
            has_no_files_to_rm = has_no_files_to_rm and not os.path.exists(path_ologfile)
        return has_no_files_to_rm

    def write_to_csvfile_for_db_dump(self, path_output_ac_or_bm_pdb_fxmutant_dir: str):
        """
        Write ddG data from specified fxout output mutant folder to single csv file (intended for database dump).
        :param path_output_ac_or_bm_pdb_fxmutant_dir: Absolute path to pdb/fxmutant folder from where the ddG data is read.
        :return: Absolute path to csv dumpfile which is the file that ddG data is written to.
        """
        ddG = ''
        csvfile_header = ''
        fxout_prefix = ''
        path_output_ac_or_bm = ''
        path_output_ac_or_bm_ddG_csv_dumpfile = ''
        is_writing_ddG = False
        is_writing_interaction_ddG = False
        ddG_lines_row_index = 0
        ddG_value_cell_index = 0
        if Paths.DIR_BM.value in path_output_ac_or_bm_pdb_fxmutant_dir:
            is_writing_ddG = True
            ddG = '_ddG'
            csvfile_header = 'pdb,fxmutant,ddG'
            fxout_prefix = self.Strs.AVG_BMDL_.value
            path_output_ac_or_bm = Paths.OUTPUT_BM
            ddG_lines_row_index = self.Strs.AVG_BM_FILE_DDG_LINE_INDEX.value
            ddG_value_cell_index = self.Strs.AVG_BM_FILE_DDG_CELL_INDEX.value
        elif Paths.DIR_AC.value in path_output_ac_or_bm_pdb_fxmutant_dir:
            is_writing_interaction_ddG = True
            ddG = '_interaction_ddG'
            fxout_prefix = self.Strs.SMRY_AC_.value
            csvfile_header = 'pdb,fxmutant,interaction ddG'
            path_output_ac_or_bm = Paths.OUTPUT_AC
            ddG_lines_row_index = self.Strs.SMRY_AC_FILE_INTER_ENERGY_LINE_INDEX.value
            ddG_value_cell_index = self.Strs.SMRY_AC_FILE_INTER_ENERGY_CELL_INDEX.value
        else:
            print(
                "The output directory path does not contain build_model or analyse_complex. If you are trying to read/write "
                "Agadir-related files, you'll need to complete the logic for it here (in FoldX class)")
        pdbname = path_output_ac_or_bm_pdb_fxmutant_dir.split('/')[-2]
        fxmutantname = path_output_ac_or_bm_pdb_fxmutant_dir.split('/')[-1]
        path_output_ac_or_bm_pdb_fxmutant_ddG_fxoutfile = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir, fxout_prefix +
                                                                       pdbname + self.Strs.UNDRSCR1.value +
                                                                       self.Strs.FXOUTEXT.value)
        if not os.path.exists(path_output_ac_or_bm_pdb_fxmutant_ddG_fxoutfile):
            print('Warning: ' + path_output_ac_or_bm_pdb_fxmutant_ddG_fxoutfile + ' not found. Cannot find ddG value to write '
                  'to csvfile.')
        else:
            avg_ddG = 0
            path_output_ac_or_bm_ddG_csv_dumpfile = os.path.join(path_output_ac_or_bm, pdbname + '_' + fxmutantname + ddG +
                                                                 Str.CSVEXT.value)
            if not os.path.exists(path_output_ac_or_bm_ddG_csv_dumpfile):
                with open(path_output_ac_or_bm_ddG_csv_dumpfile, 'w') as ddG_csv:
                    ddG_csv.write(csvfile_header + Str.NEWLN.value)
            if is_writing_interaction_ddG:
                path_output_ac_or_bm_pdb_fxmutant_ddG_fxoutfile_wt = os.path.join(path_output_ac_or_bm_pdb_fxmutant_dir,
                                                        fxout_prefix + self.Strs.WT_.value + pdbname + self.Strs.UNDRSCR1.value
                                                                                  + self.Strs.FXOUTEXT.value)
                with open(path_output_ac_or_bm_pdb_fxmutant_ddG_fxoutfile, 'r') as ddG_f_mut:
                    ddG_mut_lines = ddG_f_mut.readlines()
                with open(path_output_ac_or_bm_pdb_fxmutant_ddG_fxoutfile_wt, 'r') as ddG_f_wt:
                    ddG_wt_lines = ddG_f_wt.readlines()
                    avg_ddG_mut_values = ddG_mut_lines[ddG_lines_row_index]
                    avg_ddG_mut_values = avg_ddG_mut_values.split()
                    avg_ddG_mut = float(avg_ddG_mut_values[ddG_value_cell_index])
                    avg_ddG_wt_values = ddG_wt_lines[ddG_lines_row_index]
                    avg_ddG_wt_values = avg_ddG_wt_values.split()
                    avg_ddG_wt = float(avg_ddG_wt_values[ddG_value_cell_index])
                    avg_ddG = str(avg_ddG_wt - avg_ddG_mut)
            elif is_writing_ddG:
                with open(path_output_ac_or_bm_pdb_fxmutant_ddG_fxoutfile, 'r') as ddG_f:
                    ddG_lines = ddG_f.readlines()
                    avg_ddG_values = ddG_lines[ddG_lines_row_index]
                    avg_ddG_values = avg_ddG_values.split()
                    avg_ddG = avg_ddG_values[ddG_value_cell_index]
            with open(path_output_ac_or_bm_ddG_csv_dumpfile, 'a+') as csv_f:
                csv_f.write(pdbname + ',' + fxmutantname + ',' + avg_ddG + Str.NEWLN.value)
        return path_output_ac_or_bm_ddG_csv_dumpfile

    class BuildModel(object):

        def __init__(self, cond: dict):
            """
            :param cond: temperature, pH, ionic strength of solute, as key-value pairs.
            """
            self.conditions = cond

        def mutate_protein_structure(self, path_pdbfile: str, amino_acids: list, specific_fxmutants=None, write_to_csvfile=True):
            """
            Mutate all amino acids in this pdb to all listed amino acids using FoldX BuildModel.
            Note: individual_list needs to go in the same dir as runscript.txt.
            Note: "fxmutant" names refers specifically to the FoldX naming convention for mutants: wt_aa|chain|position|mutant_aa.
            :param path_pdbfile: Absolute path to pdbfile.
            :param amino_acids: Amino acids that proteins will be mutated to.
            :param specific_fxmutants: Specific mutants to perform BuildModel computation on. Remains empty if all mutations
            should be calculated (specified by amino_acids list).
            :param write_to_csvfile: True to write ddG to 1 csv dump file for all mutants in this pdb.
            """
            pdbfile = path_pdbfile.split('/')[-1]
            pdbname = pdbfile.split('.')[0]
            pdbname_chain_startpos_fasta_dict = GUM.extract_pdbname_chain_startpos_fasta_from_pdbs(path_pdbfile)
            fx = FoldX()
            fxmutantnames = specific_fxmutants if specific_fxmutants else fx.make_fxmutantnames(amino_acids,
                                                                                                pdbname_chain_startpos_fasta_dict)
            for specific_fxmutant in specific_fxmutants:
                if specific_fxmutant not in fxmutantnames:
                    raise ValueError('The specified mutant ' + specific_fxmutant + ' that is in your specific_fxmutants list is '
                                     ' not valid for this pdb!')
            action = fx.Strs.BM_ACTN.value + Str.HASH.value + ',' + fx.Strs.INDIV_LST_TXT.value
            path_runscript_dir = os.path.join(Paths.CONFIG_BMRUNSCRIPT, pdbname)
            fx.write_runscript_file(path_runscript_dir, pdbfile, self.conditions, action)
            path_runscript_file = os.path.join(path_runscript_dir, fx.Strs.RNSCRPT_TXT.value)
            using_cluster = GUM.using_cluster()
            for fxmutantname in fxmutantnames:
                if fxmutantname[0] == fxmutantname[len(fxmutantname) - 1]:
                    continue
                path_output_bm_pdb_fxmutant_dir = GUM.os_makedirs(Paths.OUTPUT, Paths.DIR_BM.value, pdbname, fxmutantname)
                if not self.has_already_generated_avg_bm_fxoutfile(path_output_bm_pdb_fxmutant_dir):
                    GUM.linux_copy_all_files_in_dir(path_src_dir=Paths.CONFIG_FX, path_dst_dir=path_output_bm_pdb_fxmutant_dir,
                                                    files_only=True)
                    self._write_individual_list_for_mutant(path_dst_dir=path_output_bm_pdb_fxmutant_dir)
                    os.chdir(path_output_bm_pdb_fxmutant_dir)
                    path_files_to_copy = [path_runscript_file, path_pdbfile]
                    GUM.linux_copy_specified_files(path_src_files=path_files_to_copy, path_dst_dir=path_output_bm_pdb_fxmutant_dir)
                    if using_cluster:
                        path_jobq = GUM.os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fxmutantname)
                        jobname = Paths.PREFIX_FX_BM.value + fxmutantname
                        Cluster.write_job_q_bash(jobname=jobname, path_job_q_dir=path_jobq, using_runscript=True,
                                                 path_runscript_dir=path_runscript_dir)
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
                # The following section can be run in a separate thread, subject to checking a list of completed jobs,
                # For example, a producer/consumer model.
                else:
                    print(fx.Strs.AVG_BMDL_.value + pdbname + fx.Strs.FXOUTEXT.value + ' already exists for: ' + fxmutantname)
                if write_to_csvfile:
                    print('Writing values to csv dump file for: ' + pdbname + '_' + fxmutantname)
                    while not self.has_already_generated_avg_bm_fxoutfile(path_output_bm_pdb_fxmutant_dir):
                        print(fx.Strs.AVG_BMDL_.value + pdbname + fx.Strs.FXOUTEXT.value + ' has not been created yet for ' +
                              fxmutantname + '. Therefore cannot write to csv file just yet.')
                        if using_cluster:
                            jobname = Paths.PREFIX_FX_BM.value + fxmutantname
                            Cluster.wait_for_grid_engine_job_to_complete(jobname)
                    fx.write_to_csvfile_for_db_dump(path_output_bm_pdb_fxmutant_dir)

        def has_already_generated_avg_bm_fxoutfile(self, path_output_bm_pdb_fxmutant_dir: str):
            """
            :param path_output_bm_pdb_fxmutant_dir:
            :return: True if there are atleast as many AnalyseComplex_RepairPDB_x_1_...fxout files as expected by the number of
            repaired pdb files that would be used for the computation.
            """
            path_output_bm_pdb_fxmutant_dir = path_output_bm_pdb_fxmutant_dir.split('/')
            fxmutant_dir = path_output_bm_pdb_fxmutant_dir[-1]
            pdbname = path_output_bm_pdb_fxmutant_dir[-2]
            fx = FoldX()
            path_avg_buildmodel_fxoutfile = os.path.join(Paths.OUTPUT_BM, pdbname, fxmutant_dir, fx.Strs.AVG_BMDL_.value +
                                                         pdbname + fx.Strs.FXOUTEXT.value)
            return os.path.exists(path_avg_buildmodel_fxoutfile)

        def find_num_of_missing_avg_bm_fxoutfiles(self, path_pdbfile: str, amino_acids: list):
            """
            Checks that BuildModel ddGs (in avg_bm_repairedpdb_fxout) have been calculated for all expected mutants of the
            specified pdb, according to the amino acid list.
            :param path_pdbfile: Absolute path to pdb file.
            :param amino_acids: List of amino acids that the pdb is mutated to at every position, thereby indicating which
            repaired mutant pdbs should have generated (previously by BuildModel).
            :return: True if all interactions energies for the specified pdb and amino acid mutations have been calculated.
            """
            num_of_missing_mutant_files = 0
            pdbfile = path_pdbfile.split('/')[-1]
            pdbname = pdbfile.split('.')[0]
            pdbname_chain_startpos_fasta_dict = GUM.extract_pdbname_chain_startpos_fasta_from_pdbs(path_pdbfile)
            fx = FoldX()
            expected_fxmutantnames = fx.make_fxmutantnames(amino_acids, pdbname_chain_startpos_fasta_dict)
            path_output_bm_pdb_fxmutant_dirs = sorted(glob.glob(os.path.join(Paths.OUTPUT_BM, pdbname, '*')))
            num_of_fxmutantnames_dirs = len(path_output_bm_pdb_fxmutant_dirs)
            num_of_expected_fxmutantnames = len(expected_fxmutantnames)
            if num_of_expected_fxmutantnames != num_of_fxmutantnames_dirs:
                print('You do not have the expected number of output folders for this pdb file.\nExpected number of fxmutant '
                      'folders for ' + pdbname + ' = ' + str(num_of_expected_fxmutantnames) + '\nActual number of fxmutant '
                      'folders = ' + str(num_of_fxmutantnames_dirs))
                for expected_fxmutantname in expected_fxmutantnames:
                    path_output_bm_pdb_expected_fxmutantname = os.path.join(Paths.OUTPUT_BM, pdbname, expected_fxmutantname)
                    if path_output_bm_pdb_expected_fxmutantname not in path_output_bm_pdb_fxmutant_dirs:
                        num_of_missing_mutant_files += 1
                        print('Specifically, you are missing ' + expected_fxmutantname + ' for pdb: ' + pdbname)
            else:
                print('Correct number of output_data/buildmodel/' + pdbname + '/<fxmutantnames> folders')
                for path_output_bm_pdb_fxmutant_dir in path_output_bm_pdb_fxmutant_dirs:
                    if not self.has_already_generated_avg_bm_fxoutfile(path_output_bm_pdb_fxmutant_dir):
                        num_of_missing_mutant_files += 1
                        fxmutantname = os.path.basename(path_output_bm_pdb_fxmutant_dir)
                        print('However, ' + fx.Strs.AVG_BMDL_.value + pdbname + fx.Strs.FXOUTEXT.value + ' file is missing from '
                              + fxmutantname)
            return num_of_missing_mutant_files

        def write_bm_avg_fxout_to_1csvfile_up_2dirlevels(self, path_output_bm_pdb_fxmutant_dir):
            """
            Reads the Average_BuildModel_..fxout file and writes a csv file of the data lines only (excluding top 8 lines (
            FoldX version and Consortium info, pdb name, output type..). The csv filename incorporates the pdbname and
            mutation, so there is no longer the need for the extra pdb folder and the many fxmutant subfolders, hence it is
            written two levels up, to /output_data/analyse_complex/. The fxmtuant directory is then deleted, saving about 20KB
            per mutant.
            :param path_output_bm_pdb_fxmutant_dir: Absolute path for mutant directory holding Average_BuildModel_..fxout file.
            :return: New csv output file (with its absolute path).
            """
            path_output_bm_pdb_fxmutant_dir = path_output_bm_pdb_fxmutant_dir.split('/')
            pdbname = path_output_bm_pdb_fxmutant_dir[-2]
            fxmutantname = path_output_bm_pdb_fxmutant_dir[-1]
            path_output_bm_pdb_fxmutant_dir = '/'.join(path_output_bm_pdb_fxmutant_dir)
            fx = FoldX()
            path_output_bm_pdb_fxmutant_avg_fxoutfile = os.path.join(path_output_bm_pdb_fxmutant_dir, fx.Strs.AVG_BMDL_.value +
                                                                     pdbname + fx.Strs.FXOUTEXT.value)
            path_output_bm_avg_csvfile = os.path.join(Paths.OUTPUT_BM, pdbname + '_' + fxmutantname + '_' +
                                                      fx.Strs.AVG_BMDL_CSV.value)
            if not os.path.exists(path_output_bm_avg_csvfile):
                with open(path_output_bm_avg_csvfile, 'w') as csv_f:
                    csv_f.write('pdb,fxmutant,ddG\n')
            if not os.path.exists(path_output_bm_pdb_fxmutant_avg_fxoutfile):
                if not os.path.exists(path_output_bm_avg_csvfile):
                    print('Warning: ' + path_output_bm_pdb_fxmutant_avg_fxoutfile + ' not found. But ' +
                          path_output_bm_avg_csvfile + ' also not found. It seems the Average_BuildModel_ value for this '
                                                           'mutant has not been calculated yet. Otherwise the data files have '
                                                           'been transferred to a new location or lost.')
                    return ''
                else:
                    print(path_output_bm_avg_csvfile + ' has already been written.')
                    return path_output_bm_avg_csvfile
            else:
                with open(path_output_bm_pdb_fxmutant_avg_fxoutfile) as avg_f:
                    avg_file_lines = avg_f.readlines()
                with open(path_output_bm_avg_csvfile, 'w') as csv_f:
                    # csv_f.write('pdb,fxmutant,ddG\n')
                    csv_f.write(pdbname + ',' + fxmutantname + ',')
                    avg_ddG_values = avg_file_lines[fx.Strs.AVG_BM_FILE_DDG_LINE_INDEX.value]
                    avg_ddG_values = avg_ddG_values.split()
                    csv_f.write(avg_ddG_values[fx.Strs.AVG_BM_FILE_DDG_CELL_INDEX.value])
            return path_output_bm_avg_csvfile

        def _write_individual_list_for_mutant(self, path_dst_dir: str):
            """
            Writes the 'individual_list.txt' file that is used by fx/runscript to identify the protein to process.
            :param path_dst_dir: Abs path to the output dir of foldx computation for each mutant
            """
            fxmutantname = path_dst_dir.split('/')[-1]
            with open(os.path.join(path_dst_dir, FoldX().Strs.INDIV_LST_TXT.value), 'w') as f:
                f.write(fxmutantname + ';\n')

    class AnalyseComplex(object):

        def __init__(self, cond: dict):
            """
            :param cond: temperature, pH, ionic strength of solute, as key-value pairs.
            """
            self.conditions = cond

        def calculate_complex_energies(self, path_pdbfile: str, specific_fxmutants=None, write_to_csvfile=True):
            """
            Calculates interaction energies for a complex of protein chains for a specified pdb and list of amino acids.
            The repaired pdbs are read from each output_data/buildmodel/<pdbname>/<fxmutantname> folder.
            :param path_pdbfile:
            :param specific_fxmutants: Given when only specific mutants should be computed (i.e. when some error/interruption has
            resulted in  a few mutants not being computed. There is no time wasted checking thousands of folders for output
            files when you only want to run a handful of mutants).
            :param write_to_csvfile: True to write interaction ddG to 1 csv dump file for all mutants in this pdb.
            """
            jobname = ''
            pdbfile = path_pdbfile.split('/')[-1]
            pdbname = pdbfile.split('.')[0]
            wt_pdbname = FoldX.Strs.WT_.value + pdbname
            fx = FoldX()
            using_cluster = GUM.using_cluster()
            path_output_bm_pdb_fxmutant_dirs = []
            if specific_fxmutants is not None:
                for specific_fxmutant in specific_fxmutants:
                    if specific_fxmutant[0] == specific_fxmutant[len(specific_fxmutant) - 1]:
                        continue
                    path_output_bm_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname, specific_fxmutant))
            else:
                path_output_bm_pdb_fxmutant_dirs = sorted(glob.glob(os.path.join(Paths.OUTPUT_BM, pdbname, '*')))
                if not path_output_bm_pdb_fxmutant_dirs:
                    print('No BuildModel output folders found for this pdb. AnalyseComplex requires repaired pdbs files for each '
                          'mutant. Checking the AnalyseComplex mutant output folders for the necessary repaired pdbs .....')
                path_output_ac_pdb_fxmutant_dirs = sorted(glob.glob(os.path.join(Paths.OUTPUT_AC, pdbname, '*')))
                if not path_output_ac_pdb_fxmutant_dirs:
                    print('No AnalyseComplex output folders found for this pdb either. Try running BuildModel first.')
                else:
                    for path_output_ac_pdb_fxmutant in path_output_ac_pdb_fxmutant_dirs:
                        num_of_repaired_pdbs = fx._get_num_of_repaired_pdbfiles(path_output_ac_pdb_fxmutant)
                        if num_of_repaired_pdbs == 0:
                            print('Necessary repaired pdbs (with format _1_0.fxout) are not present in the AnalyseComplex output '
                                  'folder either for: ' + path_output_ac_pdb_fxmutant)

            action = fx.Strs.AC_ACTN.value + Str.HASH.value
            # action += '0'  # According to FoldX manual (version3), write 0 for all possible interation in a complex.
            # Write AB for example if only interested in iteraction energies of chain A with chain B.
            has_one_chain_only = GUM.get_num_of_chains(path_pdbfile) == 1
            for path_output_bm_pdb_fxmutant_dir in path_output_bm_pdb_fxmutant_dirs:
                if has_one_chain_only:
                    print('Only one protein chain. No interactions to be quantified for this pdbfile.')
                    return
                path_output_bm_pdb_fxmutant_repaired = glob.glob(os.path.join(path_output_bm_pdb_fxmutant_dir,
                                                                              pdbname + fx.Strs.UNDRSCR1.value + '*' +
                                                                              Str.PDBEXT.value))
                path_output_bm_pdb_fxmutant_wt_repaired = glob.glob(os.path.join(path_output_bm_pdb_fxmutant_dir, wt_pdbname +
                                                                                 fx.Strs.UNDRSCR1.value + '*' + Str.PDBEXT.value))
                if not os.path.exists(path_output_bm_pdb_fxmutant_repaired[0]) or not os.path.exists(str(
                        path_output_bm_pdb_fxmutant_wt_repaired[0])):
                    print('Repaired pdb files for this mutant cannot be found at the expected BuildModel output location: ' +
                    path_output_bm_pdb_fxmutant_dir + '. AnalyseComplex cannot be run for this pdb/mutation. Check BuildModel '
                                                      'has been run for this pdb/mutation.')
                    return
                fxmutantname = os.path.basename(path_output_bm_pdb_fxmutant_dir)
                path_output_ac_pdb_fxmutant = GUM.os_makedirs(Paths.OUTPUT, Paths.DIR_AC.value, pdbname, fxmutantname)
                pdbs_to_analyse = []
                if not self.has_already_generated_summary_ac_fxoutfile(path_output_ac_pdb_fxmutant):
                    num_of_repaired_pdbs = fx._get_num_of_repaired_pdbfiles(path_output_bm_pdb_fxmutant_dir)
                    if num_of_repaired_pdbs == 1:
                        pdbs_to_analyse = [pdbname + fx.Strs.UNDRSCR1.value + Str.PDBEXT.value, ',', wt_pdbname +
                                           fx.Strs.UNDRSCR1.value + Str.PDBEXT.value]
                    elif num_of_repaired_pdbs > 1:
                        for i in range(1, num_of_repaired_pdbs):
                            pdbs_to_analyse.append(',')
                            pdbs_to_analyse.append(pdbname)
                            pdbs_to_analyse.append(FoldX.Strs.UNDRSCR1_012_SUFFIX_PDBS.value[i])
                            pdbs_to_analyse.append(',')
                            pdbs_to_analyse.append(wt_pdbname)
                            pdbs_to_analyse.append(FoldX.Strs.UNDRSCR1_012_SUFFIX_PDBS.value[i])
                    pdbs_to_analyse = ''.join(pdbs_to_analyse)
                    fx.write_runscript_file(path_output_ac_pdb_fxmutant, pdbs_to_analyse, self.conditions, action)
                    path_runscript_file = os.path.join(path_output_ac_pdb_fxmutant, fx.Strs.RNSCRPT_TXT.value)
                    GUM.linux_copy_all_files_in_dir(Paths.CONFIG_FX, path_output_ac_pdb_fxmutant, files_only=True)
                    os.chdir(path_output_ac_pdb_fxmutant)
                    files_to_copy = pdbs_to_analyse.split(',')
                    path_files_to_copy = []
                    for file_to_copy in files_to_copy:
                        path_files_to_copy.append(os.path.join(path_output_bm_pdb_fxmutant_dir, file_to_copy))
                    GUM.linux_copy_specified_files(path_files_to_copy, path_dst_dir=path_output_ac_pdb_fxmutant)
                    if using_cluster:
                        jobname = Paths.PREFIX_FX_AC.value + fxmutantname
                        Cluster.write_job_q_bash(jobname=jobname, path_job_q_dir=path_output_ac_pdb_fxmutant,
                                                 using_runscript=True, path_runscript_dir=path_output_ac_pdb_fxmutant)
                        if os.path.exists(path_runscript_file):
                            Cluster.run_job_q(path_job_q_dir=path_output_ac_pdb_fxmutant)
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
                else:
                    print(fx.Strs.SMRY_AC_.value + pdbname + fx.Strs.FXOUTEXT.value + ' already exists for: ' + fxmutantname)
                # The following section can be run in a separate thread, subject to checking a list of completed jobs,
                # For example, a producer/consumer model.
                if write_to_csvfile:
                    i = 0
                    while not self.has_already_generated_summary_ac_fxoutfile(path_output_ac_pdb_fxmutant):
                        print(fx.Strs.SMRY_AC_.value + pdbname + fx.Strs.FXOUTEXT.value + ' has not been created yet for ' +
                              fxmutantname + '. Therefore cannot write to csv file just yet... (waiting 10 seconds to try again)')
                        i += 6
                        if i > 10:
                            raise ValueError('AnalyseComplex has failed to produce the required ' + fx.Strs.SMRY_AC_.value +
                                             pdbname + fx.Strs.FXOUTEXT.value + ' file for: ' + fxmutantname + ' after 60 '
                                             'seconds. Something is wrong, aborting with a ValueError. Perhaps check the '
                                             'cluster logfiles??')
                    print('Writing values to csv dump file for: ' + pdbname + '_' + fxmutantname)
                    if using_cluster:
                        jobname = Paths.PREFIX_FX_AC.value + fxmutantname
                        Cluster.wait_for_grid_engine_job_to_complete(jobname)
                    fx.write_to_csvfile_for_db_dump(path_output_ac_pdb_fxmutant)

        def has_already_generated_summary_ac_fxoutfile(self, path_output_ac_pdb_fxmutant_dir: str):
            """
            :param path_output_ac_pdb_fxmutant_dir: Absolute path to the output folder for a specific mutant.
            :return: True if Summary_AnalyseComplex_RepairPDB_x_1_...fxout file is found in the mutant directory indicating
            that this mutant has already been run through the Foldx AnalyseComplex algorithm successfully.
            """
            path_output_ac_pdb_fxmutant_dir = path_output_ac_pdb_fxmutant_dir.split('/')
            pdbname = path_output_ac_pdb_fxmutant_dir[-2]
            path_output_ac_pdb_fxmutant_dir = '/'.join(path_output_ac_pdb_fxmutant_dir)
            fx = FoldX()
            path_smry_ac_pdbname_1_fxout = glob.glob(os.path.join(path_output_ac_pdb_fxmutant_dir, fx.Strs.SMRY_AC_.value +
                                                                  pdbname + fx.Strs.UNDRSCR1.value + '*' +
                                                                  fx.Strs.FXOUTEXT.value))
            path_smry_ac_wtpdbname_1_fxout = glob.glob(os.path.join(path_output_ac_pdb_fxmutant_dir, fx.Strs.SMRY_AC_.value +
                                                                    fx.Strs.WT_.value + pdbname + fx.Strs.UNDRSCR1.value + '*'
                                                                    + fx.Strs.FXOUTEXT.value))
            return len(path_smry_ac_pdbname_1_fxout) > 0 and len(path_smry_ac_pdbname_1_fxout) == len(
                path_smry_ac_wtpdbname_1_fxout)

        def find_num_of_missing_interaction_energies(self, path_pdbfile: str, amino_acids: list):
            """
            Checks that AnalyseComplex interaction energies have been calculated for all expected mutants of the specified pdb,
            according to the amino acid list.
            :param path_pdbfile: Absolute path to pdb file.
            :param amino_acids: List of amino acids that the pdb is mutated to at every position, thereby indicating which
            repaired mutant pdbs should have generated (previously by BuildModel).
            :return: True if all interactions energies for the specified pdb and amino acid mutations have been calculated.
            """
            num_of_missing_mutant_files = 0
            pdbfile = path_pdbfile.split('/')[-1]
            pdbname = pdbfile.split('.')[0]
            pdbname_chain_startpos_fasta_dict = GUM.extract_pdbname_chain_startpos_fasta_from_pdbs(path_pdbfile)
            fx = FoldX()
            expected_fxmutantnames = []
            num_of_expected_fxmutantnames = 0
            if GUM.get_num_of_chains(path_pdbfile) != 1:
                expected_fxmutantnames = fx.make_fxmutantnames(amino_acids, pdbname_chain_startpos_fasta_dict)
                num_of_expected_fxmutantnames = len(expected_fxmutantnames)
            path_output_ac_pdb_fxmutant_dirs = sorted(glob.glob(os.path.join(Paths.OUTPUT_AC, pdbname, '*')))
            num_of_fxmutantnames_dirs = len(path_output_ac_pdb_fxmutant_dirs)
            if num_of_expected_fxmutantnames != num_of_fxmutantnames_dirs:
                print('You do not have the expected number of output folders for this pdb file.\nExpected number of fxmutant '
                      'folders for ' + pdbname + ' = ' + str(num_of_expected_fxmutantnames) + '\nActual number of fxmutant '
                      'folders = ' + str(num_of_fxmutantnames_dirs))
                for expected_fxmutantname in expected_fxmutantnames:
                    path_output_ac_pdb_expected_fxmutantname = os.path.join(Paths.OUTPUT_AC, pdbname, expected_fxmutantname)
                    if path_output_ac_pdb_expected_fxmutantname not in path_output_ac_pdb_fxmutant_dirs:
                        num_of_missing_mutant_files += 1
                        print('Specifically, you are missing ' + expected_fxmutantname + ' for pdb: ' + pdbname)
            else:
                print('You have the correct number of output_data/analysecomplex/' + pdbname + '/<fxmutantnames> folders')
                for path_output_ac_pdb_fxmutant_dir in path_output_ac_pdb_fxmutant_dirs:
                    if not self.has_already_generated_summary_ac_fxoutfile(path_output_ac_pdb_fxmutant_dir):
                        num_of_missing_mutant_files += 1
                        fxmutantname = os.path.basename(path_output_ac_pdb_fxmutant_dir)
                        print('However, ' + fx.Strs.SMRY_AC_.value + pdbname + fx.Strs.FXOUTEXT.value + ' file is missing from ' +
                              fxmutantname)
            return num_of_missing_mutant_files

        def rm_all_sumry_except_1_0(self, path_output_ac_pdb_fxmutant_dir: str):
            """
            In cases where FoldX has been run more than once, producing multiple outputs of all fxout files including the
            Summary_AnalyseComplex_ fxout files. All are deleted except for one (_1_0.fxout), one for the mutant and one for the WT.
            :param path_output_ac_pdb_fxmutant_dir:
            """
            fx = FoldX()
            path_sumry_files = glob.glob(os.path.join(path_output_ac_pdb_fxmutant_dir, fx.Strs.SMRY_AC_.value + '*'))
            for path_sumry_file in path_sumry_files:
                if fx.Strs.UNDRSCR1_0_FXOUT.value in path_sumry_file:
                    continue
                else:
                    GUM.linux_remove_file(path_sumry_file)

        def write_ac_sumry_fxout_to_1csvfile_up_2dirlevels(self, path_output_ac_pdb_fxmutant_dir: str):
            """
            Reads the Summary_AnalyseComplex_..fxout file and writes a csv file of the data lines only (excluding top 8
            lines (FoldX version and Consortium info, pdb name, output type..). The csv filename incorporates the pdbname and
            mutation, so there is no longer the need for the extra pdb folder and the many fxmutant subfolders, hence it is
            written two levels up, to /output_data/analyse_complex/. The fxmtuant directory is then deleted, saving about 20KB
            per mutant.
            :param path_output_ac_pdb_fxmutant_dir: absolute path to mutant directory holding the
            Summary_AnalyseComplex_..fxout file.
            :return: New csv output file (with its absolute path)
            """
            path_output_ac_pdb_fxmutant_dir = path_output_ac_pdb_fxmutant_dir.split('/')
            pdbname = path_output_ac_pdb_fxmutant_dir[-2]
            fxmutantname = path_output_ac_pdb_fxmutant_dir[-1]
            path_output_ac_pdb_fxmutant_dir = '/'.join(path_output_ac_pdb_fxmutant_dir)
            fx = FoldX()
            path_output_ac_pdb_fxmutant_sumry_fxoutfile = os.path.join(path_output_ac_pdb_fxmutant_dir, fx.Strs.SMRY_AC_.value
                                                                       + pdbname + fx.Strs.UNDRSCR1_0_FXOUT.value)
            path_output_ac_pdb_fxmutant_sumry_fxoutfile_wt = os.path.join(path_output_ac_pdb_fxmutant_dir,
                                                                          fx.Strs.SMRY_AC_.value + fx.Strs.WT_.value + pdbname
                                                                          + fx.Strs.UNDRSCR1_0_FXOUT.value)
            path_output_ac_sumry_csvfile = os.path.join(Paths.OUTPUT_AC, pdbname + '_' + fxmutantname + '_' +
                                                        fx.Strs.SMRY_AC_CSV.value)
            if not os.path.exists(path_output_ac_sumry_csvfile):
                with open(path_output_ac_sumry_csvfile, 'w') as csv_f:
                    csv_f.write('pdb,fxmutant,interaction ddG\n')
            if not os.path.exists(path_output_ac_pdb_fxmutant_sumry_fxoutfile):
                print(path_output_ac_pdb_fxmutant_sumry_fxoutfile + ' has not been written yet')
            else:
                with open(path_output_ac_pdb_fxmutant_sumry_fxoutfile) as sumry_f:
                    sumry_file_lines = sumry_f.readlines()
                with open(path_output_ac_pdb_fxmutant_sumry_fxoutfile_wt) as sumry_wt_f:
                    sumry_wt_file_lines = sumry_wt_f.readlines()
                with open(path_output_ac_sumry_csvfile, 'w') as csv_f:
                    # csv_f.write('pdb,fxmutant,interaction energy\n')
                    csv_f.write(pdbname + ',' + fxmutantname + ',')
                    sumry_values = sumry_file_lines[fx.Strs.SMRY_AC_FILE_INTER_ENERGY_LINE_INDEX.value]
                    sumry_wt_values = sumry_wt_file_lines[fx.Strs.SMRY_AC_FILE_INTER_ENERGY_LINE_INDEX.value]
                    sumry_values = sumry_values.split()
                    sumry_wt_values = sumry_wt_values.split()
                    int_energy = float(sumry_values[fx.Strs.SMRY_AC_FILE_INTER_ENERGY_CELL_INDEX.value])
                    int_energy_wt = float(sumry_wt_values[fx.Strs.SMRY_AC_FILE_INTER_ENERGY_CELL_INDEX.value])
                    csv_f.write(str(int_energy_wt - int_energy))
            return path_output_ac_sumry_csvfile

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

        FXOUTEXT = Str.FXOUTEXT.value
        AC_ACTN = '<AnalyseComplex>'
        BM_ACTN = '<BuildModel>'
        INDIV_LST_TXT = 'individual_list' + Str.TXTEXT.value
        FXRUNSCRPT = 'FOLDX_runscript'
        FXCMDFILE = 'FOLDX_commandfile'
        FXOPTFILE = 'FOLDX_optionfile'
        RNSCRPT_TXT = Str.RNSCRPT.value + Str.TXTEXT.value
        DEFAULT_VDW_DSGN = 2
        JOB_Q = Cluster.CLSTR.JOBQFILE.value
        NO_RUNSCRPT_FILE_MSG = 'No runscript file was found'
        DSH_RUNFILE = '-runfile'
        DIF_FXOUTFILE = 'Dif_BuildModel_RepairPDBtest_1' + FXOUTEXT
        # this is temporary, just for trying something out
        DDG_CSV = 'ddG' + Str.CSVEXT.value
        ROTABASE_TXT = 'rotabase' + Str.TXTEXT.value
        CMNDS_BM_TXT = 'commands_buildmodel' + Str.TXTEXT.value
        CMNDS_STAB_TXT = 'commands_stability' + Str.TXTEXT.value
        OPTNS_BM_TXT = 'options_buildmodel' + Str.TXTEXT.value
        OPTNS_STAB_TXT = 'options_stability' + Str.TXTEXT.value
        UNDRSCR1_012_SUFFIX_PDBS = ['_1_0' + Str.PDBEXT.value, '_1_1' + Str.PDBEXT.value, '_1_2' + Str.PDBEXT.value]
        UNDRSCR1 = '_1'
        UNDRSCR1_ = UNDRSCR1 + '_'
        UNDRSCR1_0 = UNDRSCR1_ + '0'
        UNDRSCR1_0_PDB = UNDRSCR1_0 + Str.PDBEXT.value
        UNDRSCR1_0_FXOUT = UNDRSCR1_0 + FXOUTEXT
        WT_ = 'WT_'
        BMDL_ = 'BuildModel_'
        AVG_BMDL_ = 'Average_' + BMDL_
        RAW_BMDL_ = 'Raw_' + BMDL_
        DIF_BMDL_ = 'Dif_' + BMDL_
        PDBLST_BMDL = 'PdbList_' + BMDL_
        ACPLX_ = 'AnalyseComplex_'
        SMRY_AC_ = 'Summary_' + ACPLX_
        INDIV_ENRG_ACPLX_ = 'Indiv_energies_' + ACPLX_
        INTACTN_ACPLX_ = 'Interaction_' + ACPLX_
        INTRFC_RESDS_ACPLX_ = 'Interface_Residues_' + ACPLX_
        BM_FILES_TO_DELETE = [BMDL_, RAW_BMDL_, DIF_BMDL_, PDBLST_BMDL]
        AC_FILES_TO_DELETE = [ACPLX_, INDIV_ENRG_ACPLX_, INTACTN_ACPLX_, INTRFC_RESDS_ACPLX_]
        SMRY_AC_CSV = 'summary_ac' + Str.CSVEXT.value
        AVG_BMDL_CSV = 'avg_bm' + Str.CSVEXT.value
        Pdb_FXOUTFILE_HEADER = 'Pdb'
        TOTL_ddG = 'Total energy (ddG)'
        AVG_BM_FILE_DDG_LINE_INDEX = 9
        AVG_BM_FILE_DDG_CELL_INDEX = 2
        SMRY_AC_FILE_INTER_ENERGY_CELL_INDEX = 5
        SMRY_AC_FILE_INTER_ENERGY_LINE_INDEX = 9


# pydevd.stoptrace()
