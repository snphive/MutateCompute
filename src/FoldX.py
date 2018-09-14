import subprocess
import os
from src.Str import Str
from src.Paths import Paths
from src.GeneralUtilityMethods import GUM
from src.Cluster import Cluster


class FoldX(object):

    def __init__(self, conditions: object):
        """
        :param conditions: different associated conditions in dictionaries
        """
        self.temp = str(conditions['temp'])
        self.ph = str(conditions['ph'])
        self.ion_strgth = str(conditions['ion_strgth'])

    # 30.07.18 Redesigned the directory structure such that runscripts will go in configuration/foldx/ & maybe another
    # level such as analyse_complex or build_model or stability etc.
    def write_runscript_textfile(self, path_runscript: str, pdbs: str, action: str, show_sequence_detail=False,
                                 print_networks=False, calculate_stability=False):
        """
        The runscript.txt is an input file for FoldX indicating which pdbs to analyse & which programs to run on them.
        :param path_runscript: Absolute path for runscript.txt file being written.
        :param pdbs: pdb(s) (including .pdb extension) inputs for FoldX.
        :param action: FoldX computation to be performed (e.g. BuildModel, AnalyseComplex, etc).
        :param show_sequence_detail: True will provide extra information in output. 
        :param print_networks: True to output network data on .. ?
        :param calculate_stability: True include a stability calculation.
        :return: 
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
        runscript.append('<Temperature>' + self.temp + Str.SEMICO_NL.value)
        runscript.append('<IonStrength>' + self.ion_strgth + Str.SEMICO_NL.value)
        runscript.append('<ph>' + self.ph + Str.SEMICO_NL.value)
        runscript.append('<moveNeighbours>' + Str.TRUE_lc.value + Str.SEMICO_NL.value)
        runscript.append('<VdWDesign>2' + Str.SEMICO_NL.value)
        runscript.append('<numberOfRuns>3' + Str.SEMICO_NL.value)
        runscript.append('<OutPDB>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<END>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<JOBEND>' + Str.HASH.value + Str.SEMICO_NL.value)
        runscript.append('<ENDFILE>' + Str.HASH.value + Str.SEMICO_NL.value)
        with open(os.path.join(path_runscript, self.Strs.runscrpt_txt.value), 'w') as runscript_file:
            runscript_str = ''.join(runscript)
            runscript_file.write(runscript_str)
        return runscript_str

    class Repair(object):

        def __init__(self):
            print('Repair constructor')

        def do_repair(self, input_pdb):
            print('do repair method')

    class BuildModel(object):

        def __init__(self):
            self.fx = FoldX()

        # Mutate specified amino acids in this pdb to all listed amino acids using FoldX BuildModel,
        # FoldX uses a runscript file, which must be written here.
        def mutate_all_pdb(self, path_input: str, path_pdbfile: str, amino_acids: list, path_output: str):
            """
            Mutate all amino acids in this pdb to all listed amino acids using FoldX BuildModel.
            :param path_input: Absolute path of input root directory.
            :param path_pdbfile: Absolute path to pdbfile.
            :param amino_acids: Amino acids that proteins will be mutated to.
            :param path_output: Absolute path of output root directory.
            """
            pdbname = path_pdbfile.split('/')[-1].split('.')[0]
            path_input_pdbs_pdbname = GUM._os_makedirs(path_input, Paths.DIR_PDBS.value, pdbname)
            path_runscript = GUM._os_makedirs(path_input_pdbs_pdbname, Paths.DIR_BMRUNSCRIPT.value)
            pdbname_chain_fasta_dict = GUM.extract_pdbname_chain_fasta_from_pdbs(path_input, path_pdbfile,
                                                                                 path_output, )
            fx_mutant_name_list = self._make_fx_mutant_name_list(amino_acids, pdbname_chain_fasta_dict)

            action = self.fx.Strs.BM.value + Str.HASH.value + ',' + self.fx.Strs.indiv_list_txt.value
            FoldX().write_runscript_textfile(path_runscript, path_pdbfile, action)

            for fx_mutant_name in fx_mutant_name_list:
                path_jobq_indivlist_dest = GUM._os_makedirs(path_input_pdbs_pdbname, fx_mutant_name)
                self._write_individual_list_for_mutant(fx_mutant_name, path_jobq_indivlist_dest)
                job_name = Paths.PREFIX_FX_BM.value + fx_mutant_name
                if GUM.using_cluster():
                    Cluster.write_job_q_bash(job_name, path_job_q_dir=path_jobq_indivlist_dest,
                                             path_runscript_dir=path_runscript, using_runscript=True)
                    if os.path.exists(os.path.join(path_runscript, self.fx.Strs.runscrpt_txt.value)):
                        Cluster.run_job_q(path_job_q_dir=path_jobq_indivlist_dest)
                    else:
                        raise ValueError('No runscript file was found')
                else:
                    cmd = Paths.FOLDX_EXE + Str.SPCE.value + self.fx.Strs.DSH_RUNFILE.value + Str.SPCE.value + \
                          os.path.join(Paths.CONFIG_FX_BMRUNSCRIPT, self.fx.Strs.runscrpt_txt.value)
                    subprocess.call(cmd, shell=True)

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
        def _write_individual_list_for_mutant(self, fx_mutant_name: str, path_file_dest_dir: str):
            """
            Writes the 'individual_list.txt' file that is used by fx/runscript to identify the protein to process.
            :param fx_mutant_name: fx_style mutant name: concatenation of wtaa, chain, position and mutantaa.
            :param path_file_dest_dir:
            :return:
            """
            with open(os.path.join(path_file_dest_dir, self.fx.Strs.indiv_list_txt.value), 'w') as f:
                f.write(fx_mutant_name + ';\n')

    class Stability(object):

        def __init__(self):
            print('helloworld constructor')

        def calculate_stability(self, pdb):
            print('not implemented yet')

    class AnalyseComplex(object):

        def __init__(self):
            print('helloworld constructor')

        def calculate_complex_energies(self, pdbname):
            self._prepare_for_FoldX_AnalyseComplex(pdbname)
            print('not yet implemented')

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
            FoldX().write_runscript_textfile(path_to_runscript, pdbs_to_analyse, action)

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
