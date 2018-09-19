import subprocess
import os
from src.Str import Str
from src.Paths import Paths
from src.GeneralUtilityMethods import GUM
from src.Cluster import Cluster


class FoldX(object):

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
        runscript.append('<Temperature>' + str(conditions['temp']) + Str.SEMICO_NL.value)
        runscript.append('<IonStrength>' + str(conditions['ion_strgth']) + Str.SEMICO_NL.value)
        runscript.append('<ph>' + str(conditions['ph']) + Str.SEMICO_NL.value)
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

        # Mutate specified amino acids in this pdb to all listed amino acids using FoldX BuildModel,
        # FoldX uses a runscript file, which must be written here.
        def mutate_protein_structure(self, path_pdbfile: str, amino_acids: list):
            """
            Mutate all amino acids in this pdb to all listed amino acids using FoldX BuildModel.
            Note: individual_list needs to go in the same dir as runscript.txt.
            :param path_pdbfile: Absolute path to pdbfile.
            :param amino_acids: Amino acids that proteins will be mutated to.
            """
            pdbfile = path_pdbfile.split('/')[-1]
            pdbname = pdbfile.split('.')[0]
            pdbname_chain_fasta_dict = GUM.extract_pdbname_chain_fasta_from_pdbs(path_pdbfile)
            fx_mutant_name_list = self._make_fx_mutant_name_list(amino_acids, pdbname_chain_fasta_dict)
            action = FoldX().Strs.BM.value + Str.HASH.value + ',' + FoldX().Strs.indiv_list_txt.value
            FoldX().write_runscript_file(Paths.CONFIG_BMRUNSCRIPT, pdbfile, self.conditions, action)
            path_runscript_file = os.path.join(Paths.CONFIG_BMRUNSCRIPT, FoldX().Strs.runscrpt_txt.value)

            for fx_mutant_name in fx_mutant_name_list:
                path_output_pdbname_mutant = self._make_output_dir_and_copy_fxconfig_files_in(Paths.OUTPUT, pdbname,
                                                                                              fx_mutant_name)
                self._write_individual_list_for_mutant(path_output_pdbname_mutant)
                os.chdir(path_output_pdbname_mutant)
                path_files_to_copy = [path_runscript_file, path_pdbfile]
                GUM.linux_copy_specified_files(path_files_to_copy, path_output_pdbname_mutant)
                if GUM.using_cluster():
                    path_jobq = GUM._os_makedirs(Paths.CONFIG_BM_JOBQ, pdbname, fx_mutant_name)
                    Cluster.write_job_q_bash(job_name=Paths.PREFIX_FX_BM.value + fx_mutant_name,
                                             path_job_q_dir=path_jobq, path_runscript_dir=Paths.CONFIG_BMRUNSCRIPT,
                                             using_runscript=True)
                    if os.path.exists(os.path.join(Paths.CONFIG_BMRUNSCRIPT, FoldX().Strs.runscrpt_txt.value)):
                        Cluster.run_job_q(path_job_q_dir=path_jobq)
                    else:
                        raise ValueError('No runscript file was found')
                else:
                    cmd = 'chmod 777 /Users/u0120577/SNPEFFECT/executables/FoldX/foldx'
                    try:
                        subprocess.call(cmd, shell=True)
                    except OSError:
                        print('Problem with this linux cmd: ' + cmd)
                    cmd = Paths.FOLDX_EXE + Str.SPCE.value + FoldX().Strs.DSH_RUNFILE.value + Str.SPCE.value + \
                          path_runscript_file
                    if os.path.exists(path_runscript_file):
                        subprocess.call(cmd, shell=True)
                    else:
                        raise ValueError('No runscript file was found')

        def _make_output_dir_and_copy_fxconfig_files_in(self, path_output, pdbname, fx_mutant_name):
            path_output_pdbname_mutant = GUM._os_makedirs(path_output, Paths.DIR_BM.value, pdbname, fx_mutant_name)
            GUM.linux_copy_all_files_in_dir(Paths.CONFIG_FX, path_output_pdbname_mutant, overwrite=False)
            return path_output_pdbname_mutant

        def _make_fx_mutant_name_list(self, amino_acids: list, pdbname_chain_fasta: dict):
            """
            :param amino_acids: Amino acids that (every residue in) the protein will be mutated to.
            :param pdbname_chain_fasta: pdbname_underscore_chain and fasta sequence as key-value pair.
            :return: List of fx_style mutant names: concatenation of wtaa, chain, position and mutantaa.
            """
            fx_mutant_name_list = []
            for pdbname_chain, fasta_sequence in pdbname_chain_fasta.items():
                chain = pdbname_chain.split('_')[-1]
                for index, wt_aa in enumerate(fasta_sequence):
                    position = index + 1
                    for mutant_aa in amino_acids:
                        fx_mutant_name_list.append(wt_aa + chain + str(position) + mutant_aa)
            return fx_mutant_name_list

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

        def calculate_complex_energies(self, pdbname):
            self._prepare_for_FoldX_AnalyseComplex(pdbname)
            print('TODO')

        def _prepare_for_FoldX_AnalyseComplex(self, repair_pdbname):
            _0_1_2_pdbs = ['0.pdb,', '1.pdb,', '2.pdb,']
            repair_pdbname_1_ = repair_pdbname + '_1_'
            wt_repair_pdbname_1_ = 'WT_' + repair_pdbname_1_

            path_to_runscript = './'
            pdbs_to_analyse = repair_pdbname_1_ + _0_1_2_pdbs[0] + \
                              repair_pdbname_1_ + _0_1_2_pdbs[1] + \
                              repair_pdbname_1_ + _0_1_2_pdbs[2] + \
                              wt_repair_pdbname_1_ + _0_1_2_pdbs[0] + \
                              wt_repair_pdbname_1_ + _0_1_2_pdbs[1] + \
                              wt_repair_pdbname_1_ + _0_1_2_pdbs[2]
            action = FoldX().Strs.AC.value + Str.HASH.value
            FoldX().write_runscript_file(path_to_runscript, pdbs_to_analyse, self.conditions, action)

    from enum import Enum

    class Strs(Enum):

        AC = '<AnalyseComplex>'
        BM = '<BuildModel>'
        indiv_list_txt = 'individual_list.txt'
        FXRUNSCRPT = 'FOLDX_runscript'
        FXCMDFILE = 'FOLDX_commandfile'
        FXOPTFILE = 'FOLDX_optionfile'
        runscrpt_txt = 'runscript.txt'
        DSH_RUNFILE = '-runfile'




