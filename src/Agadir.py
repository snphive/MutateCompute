#!/usr/bin/env python3
"""
Class for launching Agadir algorithms. Includes code for writing Agadir's options file and processing output files.
"""

import os
import subprocess
from src.enums.Str import Str
from src.tools.GeneralUtilityMethods import GUM
from src.enums.Paths import Paths
from src.enums.Conditions import Cond
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class Agadir(object):

    def __init__(self, conditions):
        self.temp = str(conditions[Str.COND_TEMP.value])
        self.ph = str(conditions[Str.COND_PH.value])
        self.ion_strgth = str(conditions[Str.COND_ION_STRGTH.value])
        self.tfe = str(conditions[Str.COND_TFE.value])
        self.stab = str(conditions[Str.COND_STAB.value])
        self.conc = str(conditions[Str.COND_CONC.value])

    def run_agadir_on_multifastas(self, path_fastafile: str, path_dst: str):
        """
        To run agadir on file that has multiple fasta sequences. Iterate through the list, write the sequence to an
        individual text file passing this to agadir.compute(). This newly-written individual text file is then deleted,
        crucial to prevent too much memory taken up by each file when running large numbers of sequences.
        :param path_fastafile: Abs path to fasta file.
        :param path_dst: Abs path to output root dir.
        """
        with open(path_fastafile) as f:
            is_first_line = True
            fasta_str = ''
            mutantfastafilename = ''
            mutantfastafile = ''
            for line in f.readlines():
                if '>' in line:
                    if not is_first_line:
                        path_dst_mutant_filename = GUM.os_makedirs(path_dst, mutantfastafilename)
                        path_dst_mutant_file = os.path.join(path_dst_mutant_filename, mutantfastafile)
                        with open(path_dst_mutant_file, 'w') as g:
                            g.write(fasta_str)
                        agadir = Agadir(Cond.INCELL_MAML.value)
                        agadir.compute(path_dst_mutant_file)
                        GUM.linux_remove_file(path_dst_mutant_file)
                    fasta_str = line
                    is_first_line = False
                    mutantfastafilename = line.split('>')[-1].split('\n')[0]
                    mutantfastafile = mutantfastafilename + Str.FSTAEXT.value
                else:
                    fasta_str += line
        self.remove_unwanted_globalresidue_output_files(path_dst_mutant_filename)

    def remove_unwanted_globalresidue_output_files(self, path_dst_mutant_filename: str):
        """
        Deemed unnecessary (in JS meeting 24.10.18), the PSX_globalresidue.output file is to be deleted to save space.
        :param path_dst_mutant_filename: Absolute path of the file to be deleted.
        :return:
        """
        path_out_mutantfilename_globalresidue_out = os.path.join(path_dst_mutant_filename, Str.GLOBRESIDUE_OUT.value)
        GUM.linux_remove_file(path_out_mutantfilename_globalresidue_out)

    def compute(self, path_dst_fastafile: str):
        """
        Runs Agadirwrapper according to the options specified in /configuration/agadir_config/Options.txt file.
        Expecting path_fastafile to potentially contain more than 1 sequence.
        :param path_dst_fastafile: Abs path to fasta file.
        """
        print('Agadir.compute is called on: ' + path_dst_fastafile)
        path_agadir_dst_dir = '/'.join(path_dst_fastafile.split('/')[:-1])
        os.chdir(path_agadir_dst_dir)
        if GUM.using_cluster():
            # Cluster.wait_for_grid_engine_job_to_complete(Paths.PREFIX_AGADIR.value)
            path_agadir_exe = Paths.ZEUS_AGADIR_EXE.value + 'agadirwrapper'
            path_config_agad = Paths.SE_CONFIG_AGAD.value
        else:
            print('Local agadir exe should be: /Users/u0120577/SNPEFFECT/executables/agadir_10042012')
            cmd = 'chmod 777 /Users/u0120577/SNPEFFECT/executables/agadir_10042012/agadir_wrapper'
            try:
                subprocess.call(cmd, shell=True)
            except OSError:
                print(Str.PROBLNXCMD_MSG.value + cmd)
            path_agadir_exe = Paths.LOCAL_AGADIR_EXE.value + 'agadir_wrapper'
            path_config_agad = Paths.MC_CONFIG_AGAD.value
        # cmd = Paths.AGADIR_EXE + agadir_exe + Str.SPCE.value + path_dst_fastafile + Str.SPCE.value + \
        #       Paths.CONFIG_AGAD + '/Options.txt'
        cmd = path_agadir_exe + Str.SPCE.value + path_dst_fastafile + Str.SPCE.value + path_config_agad + '/Options.txt'
        # cmd = '/switchlab/group/tools/agadir_10042012/' + agadir_exe + Str.SPCE.value + path_dst_fastafile +
        # Str.SPCE.value + '/switchlab/group/shazib/SnpEffect/configuration/agadir_config/Options.txt'
        subprocess.call(cmd, shell=True)

    def write_agadir_options_file(self, ntrm='#', ctrm='#', global_tot=Str.TRUE_lc.value, tango_wndw=Str.TRUE_lc.value,
                                  waltz_wndw=Str.TRUE_lc.value, limbo_wndw=Str.FALSE_lc.value,
                                  agadir_wndw=Str.FALSE_lc.value, csblnca_wndw=Str.FALSE_lc.value,
                                  complex_wndw=Str.FALSE_lc.value, repeat_wndw=Str.FALSE_lc.value,
                                  pat_tango_wndw=Str.FALSE_lc.value, tango_resid=Str.FALSE_lc.value,
                                  waltz_resid=Str.FALSE_lc.value, limbo_resid=Str.FALSE_lc.value,
                                  complex_resid=Str.FALSE_lc.value, agadir_resid=Str.FALSE_lc.value,
                                  csblnca_resid=Str.FALSE_lc.value, repeat_resid=Str.FALSE_lc.value,
                                  wndws_file_per_seq=Str.FALSE_lc.value, resid_file_per_seq=Str.FALSE_lc.value):
        """
        # There are 3 types of calculations and xx algorithms which this can be applied to.
        # 1) Global total; 2) Window calculation; 3) Residue calculation.
        Global means the sum of all the scores for the whole protein.
        Residue means the score for each residue in the whole protein.
        Window means .....
        :param ntrm:
        :param ctrm:
        :param global_tot:
        :param tango_wndw:
        :param waltz_wndw:
        :param limbo_wndw:
        :param agadir_wndw:
        :param csblnca_wndw:
        :param complex_wndw:
        :param repeat_wndw:
        :param pat_tango_wndw:
        :param tango_resid:
        :param waltz_resid:
        :param limbo_resid:
        :param complex_resid:
        :param agadir_resid:
        :param csblnca_resid:
        :param repeat_resid:
        :param wndws_file_per_seq:
        :param resid_file_per_seq:
        :return: Agadir options text file.
        """
        options = []
        options.append('<TITLE>AGADIR_optionfile' + Str.SEMICO_NL.value)
        options.append('<Temperature>' + self.temp + '.' + Str.SEMICO_NL.value)
        options.append('<pH>' + self.ph + Str.SEMICO_NL.value)
        options.append('<IonStrength>' + self.ion_strgth + Str.SEMICO_NL.value)
        options.append('<TfeConc>' + self.tfe + '.' + Str.SEMICO_NL.value)
        options.append('<Stability>' + self.stab + '.' + Str.SEMICO_NL.value)
        options.append('<Concentration>' + self.conc + '.' + Str.SEMICO_NL.value)
        options.append('<Nterm>' + ntrm + Str.SEMICO_NL.value)
        options.append('<Cterm>' + ctrm + Str.SEMICO_NL.value)
        options.append('<global_total>' + global_tot + Str.SEMICO_NL.value)
        options.append('<tango_window>' + tango_wndw + Str.SEMICO_NL.value)
        options.append('<waltz_window>' + waltz_wndw + Str.SEMICO_NL.value)
        options.append('<limbo_window>' + limbo_wndw + Str.SEMICO_NL.value)
        options.append('<agadir_window>' + agadir_wndw + Str.SEMICO_NL.value)
        options.append('<casablanca_window>' + csblnca_wndw + Str.SEMICO_NL.value)
        options.append('<complex_window>' + complex_wndw + Str.SEMICO_NL.value)
        options.append('<repeat_window>' + repeat_wndw + Str.SEMICO_NL.value)
        options.append('<patentTango_window>' + pat_tango_wndw + Str.SEMICO_NL.value)
        options.append('<tango_residue>' + tango_resid + Str.SEMICO_NL.value)
        options.append('<waltz_residue>' + waltz_resid + Str.SEMICO_NL.value)
        options.append('<limbo_residue>' + limbo_resid + Str.SEMICO_NL.value)
        options.append('<complex_residue>' + complex_resid + Str.SEMICO_NL.value)
        options.append('<agadir_residue>' + agadir_resid + Str.SEMICO_NL.value)
        options.append('<casablanca_residue>' + csblnca_resid + Str.SEMICO_NL.value)
        options.append('<repeat_residue>' + repeat_resid + Str.SEMICO_NL.value)
        options.append('<windows_file_per_sequence>' + wndws_file_per_seq + Str.SEMICO_NL.value)
        options.append('<residue_file_per_sequence>' + resid_file_per_seq + Str.SEMICO_NL.value)
        options.append('<END>' + Str.NEWLN.value)
        with open(Paths.CONFIG_AGAD + '/Options.txt', 'w') as f:
            options_str = ''.join(options)
            f.write(options_str)
        return options_str

    def copy_files_output_to_input_fastasubdirs(self, path_output_fastas: str, path_input_fastas: str,
                                                filename_list=None, into_own_subdirs=True):
        """
        Copies fastafiles from subdirectories of output_data (e.g. generated by MutateFasta) to subdirs of input_data.
        Typically, from ../output_data/fastas/<fastafilename>/mutants/. to ../input_data/fastas/<fastafilename>/mutants
        :param path_output_fastas: Source directory of files to be copied.
        :param path_input_fastas: Destination directory of files to be copied.
        :param filename_list: Names of files to copy (if not all output_fastas dir), None by default.
        :param into_own_subdirs: True to copy files into new subdir bearing same name.
        """
        if filename_list is not None:
            path_fastas_filelist = [os.path.join(path_output_fastas, filename, filename + Str.FSTAEXT.value) for
                                    filename in filename_list]
            GUM.linux_copy_files(path_fastas_filelist, path_input_fastas, into_own_subdirs)
        else:
            path_input_fastas_list = path_input_fastas.split('/')
            path_input_list = path_input_fastas_list[:-1]
            path_input = '/'.join(path_input_list)
            GUM.linux_copy_all_files_in_dir(path_output_fastas, path_input, recursively=True)
            if into_own_subdirs:
                GUM._move_files_into_own_subdirs(path_dir=path_input_fastas)


# pydevd.stoptrace()
