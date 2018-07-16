from unittest import TestCase
from src.FoldX import FoldX
import yaml


class TestScheduler(TestCase):

    global path_R_exe
    global path_FoldX_exe
    global path_Agadir_exe
    global path_QSub_exe
    global path_SnpEffect
    global path_SE_SourceFiles_Scripts
    global path_SE_Inputs_PDBs
    global path_SE_Inputs_FASTAs
    global path_SE_Outputs
    global path_SE_Outputs_Agadir
    global path_SE_Outputs_FoldX
    global dict_aa_1to3
    global dict_aa_3to1
    global list_all_20_aa

    with open("~/PycharmProjects/MutateCompute/config/pathsAndDictionaries.yaml", 'r') as stream:
        try:

            paths_and_dictionaries = yaml.load(stream)
            path_R_exe = paths_and_dictionaries['ROOT']['path_R_exe']
            path_FoldX_exe = paths_and_dictionaries['ROOT']['path_FoldX_exe']
            path_Agadir_exe = paths_and_dictionaries['ROOT']['path_Agadir_exe']
            path_QSub_exe = paths_and_dictionaries['ROOT']['path_QSub_exe']

            path_SnpEffect = paths_and_dictionaries['ROOT']['path_SnpEffect']
            path_SE_SourceFiles_Scripts = paths_and_dictionaries['ROOT'][' path_SourceFiles_Scripts']
            path_SE_Inputs_PDBs = paths_and_dictionaries['ROOT']['path_SE_Inputs_PDBs']
            path_SE_Inputs_FASTAs = paths_and_dictionaries['ROOT']['path_SE_Inputs_FASTAs']
            path_SE_Outputs = paths_and_dictionaries['ROOT']['path_SE_Outputs']
            path_SE_Outputs_Agadir = paths_and_dictionaries['ROOT']['path_SE_Outputs_Agadir']
            path_SE_Outputs_FoldX = paths_and_dictionaries['ROOT']['path_SE_Outputs_FoldX']

            dict_aa_1to3 = paths_and_dictionaries['ROOT']['dict_aa_1to3']
            dict_aa_3to1 = paths_and_dictionaries['ROOT']['dict_aa_3to1']
            list_all_20_aa = paths_and_dictionaries['ROOT']['list_all_20_aa']

        except yaml.YAMLError as exc:
            print(exc)

    @classmethod
    def setUpClass(cls):
        cls.absolute_path_tests = '~/PycharmProjects/MutateCompute/tests'
        cls.abs_path_tests_inputs = cls.absolute_path_tests + '/Inputs'
        cls.abs_path_tests_inputs_fasta = cls.abs_path_tests_inputs + '/Fasta'
        cls.abs_path_tests_inputs_pdbs = cls.abs_path_tests_inputs + '/PDBs'
        cls.abs_path_tests_outputs = cls.absolute_path_tests + '/Outputs'
        cls.abs_path_tests_outputs_agadir = cls.abs_path_tests_outputs + '/Agadir'
        cls.abs_path_tests_outputs_fasta = cls.abs_path_tests_outputs + '/Fasta'
        cls.abs_path_tests_outputs_foldx = cls.abs_path_tests_outputs + '/FoldX'
        cls.abs_path_tests_outputs_foldx_buildmodel = cls.abs_path_tests_outputs_foldx + '/BuildModel'

    def test_yaml_config_parameters_uploaded(self):
        msg = 'path not as expected'
        self.assertEqual(self, '/software/shared/apps/general/R/3.1.2/bin/Rscript', path_R_exe, 'R exe' + msg)
        self.assertEqual(self, '/switchlab/group/tools/FoldX_2015/FoldX', path_FoldX_exe, 'FoldX exe' + msg)
        self.assertEqual(self, '/switchlab/group/tools/agadir_10042012/agadirwrapper', path_Agadir_exe, 'Agadirwrapper exe' + msg)
        self.assertEqual(self, '/opt/sge/bin/lx-amd64/', path_QSub_exe, 'QSub exe' + msg)
        self.assertEqual(self, '/switchlab/group/shazib/SnpEffect', path_SnpEffect, 'SnpEffect directory' + msg)
        self.assertEqual(self, '/switchlab/group/shazib/SnpEffect/SourceFiles/Scripts', path_SE_SourceFiles_Scripts, 'SourceFiles/Scripts directory' + msg)
        self.assertEqual(self, '/switchlab/group/shazib/SnpEffect/Inputs/PDBs', path_SE_Inputs_PDBs, 'Inputs/PDBs directory' + msg)
        self.assertEqual(self, '/switchlab/group/shazib/SnpEffect/Inputs/FASTAs', path_SE_Inputs_FASTAs, 'Inputs/Fastas directory' + msg)
        self.assertEqual(self, '/switchlab/group/shazib/SnpEffect/Outputs', path_SE_Outputs, '/Outputs directory' + msg)
        self.assertEqual(self, '/switchlab/group/shazib/SnpEffect/Outputs/Agadir', path_SE_Outputs_Agadir, '/Outputs/Agadir directory' + msg)
        self.assertEqual(self, '/switchlab/group/shazib/SnpEffect/Outputs/FoldX', path_SE_Outputs_FoldX, '/Outputs/FoldX directory' + msg)
        self.assertEqual(self, {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE',
                                'K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER',
                                'T':'THR','V':'VAL','W':'TRP','Y':'TYR'}, dict_aa_1to3, 'dict_aa_1to3 not as expected')
        self.assertEqual(self, {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','H1S':'H',
                                'H2S':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
                                'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}, dict_aa_3to1,
                         'dict_aa_3to1 not as expected')
        self.assertEqual(self, ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'],
                         list_all_20_aa, 'list_all_20_aa not as expected')

    def test_start(self):

        self.fail()
