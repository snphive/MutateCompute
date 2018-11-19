#!/usr/bin/env python3
"""
Class for generating mutated protein sequences (FASTA format) at every position in the protein to any other
amino acid(s).
"""
import os
from src.GeneralUtilityMethods import GUM
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from src.Paths import Paths
from src.AminoAcids import AA
from src.Str import Str

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class MutateFasta(object):

    def __init__(self, amino_acids: list):
        """
        Constructor
        :param amino_acids: all amino acids that sequence(s) will be mutated to.
        """
        self.__AAs = amino_acids

    def mutate_every_residue(self, path_fastafile: str, write_1_fasta_only: bool, write_fasta_per_mut: bool,
                             path_output_3dots: str, write_csv=False, write_txt=False):
        """
        Mutates FASTA (typically wild-type) sequences at every position, to every residue specified, (typically the 19
        other residues.)
        :param path_fastafile: Absolute path to fasta file.
        :param write_1_fasta_only: True writes 1 fasta file containing multiple fasta sequences (i.e. the mutants).
        :param write_fasta_per_mut: True writes a fasta file for each fasta sequences (i.e. 1 per mutant).
        :param path_output_3dots: Absolute path to a '3dot' subdirectory (i.e. 1...1000).
        :param write_csv: True writes a csv file containing all mutants per wildtype fasta, (False by default).
        :param write_txt: True writes a txt file containing all mutants per wildtype fasta, (False by default).
        :return: Dictionary with the wild-type protein's name as key (e.g. 1_A), and another dictionary as the value.
        The second dictionary has each mutant for this protein, with its name (e.g. R1A) as the key, the mutant's
        sequence as the value.
        """
        titleSeq = self.make_titleSeqDict_from_fastafile(path_fastafile)
        title_titleSeq = self.convert_titleSeqDict_to_titleTitleSeqDictDict(titleSeq)
        title_titleSeq_w_mutants = self._populate_title_titleSeq_with_mutants(title_titleSeq)
        self._write_mutants(title_titleSeq_w_mutants, write_1_fasta_only, write_fasta_per_mut, path_output_3dots,
                            write_csv, write_txt)
        return title_titleSeq_w_mutants

    def mutate_every_residue_to_every_aa_write_1_file(self, path_fastafile: str):
        """

        :param path_fastafile:
        """
        path_output_3dots = GUM.make_path_fastas_3dots_dirs(Paths.OUTPUT, path_fastafile)
        self.mutate_every_residue(path_fastafile=path_fastafile, write_1_fasta_only=True, write_fasta_per_mut=False,
                                  path_output_3dots=path_output_3dots)

    def _populate_title_titleSeq_with_mutants(self, title_titleSeq: dict):
        """
    
        :param title_titleSeq: Fasta title is the key. The value is mutant-title: mutated-sequence dictionary. This
        already has the dict of wt-title: (wt-title:wt-sequence) dict.
        :return: wt-title: mutant-title:mutant-sequence
        """
        title_titleSeq_w_mutants = title_titleSeq
        for wt_title in list(title_titleSeq.keys()):
            wt_seq = title_titleSeq[wt_title][wt_title]
            title_titleSeq_w_mutants[wt_title] = self._add_mutantTitle_mutatedSeq_to_dict(
                title_titleSeq[wt_title], wt_title, wt_seq)
        return title_titleSeq_w_mutants

    def _add_mutantTitle_mutatedSeq_to_dict(self, titleSeq: dict, wt_title: str, wt_seq: str):
        """
        A sequence is mutated at every position to every residue that is included in the list (mutant_aa_list).
        E.g. input as {'230498_A': 'MDVFM'} output as {'230498_A': 'MDVFM','M1A': 'ADVFM','M1C': 'CDVFM', etc.. }
        :param titleSeq: Already should have wt-title: wt-title:wt-sequence dictionary.
        :param wt_title: Name of fasta sequence (same as the fastafile name).
        :param wt_seq: Amino acid sequence in the fastafile (in FASTA format).
        :return: Dictionary of sequence name as key, FASTA sequence itself as value, starting with wild-type and
        followed by every mutant according to a specified list of residues.
        """
        titleSeqDict_w_mutants = {wt_title: titleSeq[wt_title]}
        mutable_seq = MutableSeq(wt_seq, IUPAC.protein)
        for i, wt_aa in enumerate(mutable_seq):
            for mutant_aa in self.__AAs:
                mutant_title = wt_aa + str(i + 1) + mutant_aa
                wt_aa_at_i = mutable_seq[i]
                if not mutant_aa == wt_aa_at_i:
                    mutable_seq[i] = mutant_aa
                    titleSeqDict_w_mutants[mutant_title] = str(mutable_seq)
                    mutable_seq[i] = wt_aa_at_i
        return titleSeqDict_w_mutants

    def _write_mutants(self, title_titleSeq_w_mutants: dict, write_1_fasta_only: bool, write_fasta_per_mut: bool,
                       path_output_3dots: str, write_csv=False, write_txt=False):
        """
        Writes the mutants out to fastafiles and/or csv files and/or txt files. These can be written in one file
        containing all mutants or one file per mutant.
        The fastafiles are written to /output_data/<fastafilename>/mutants/.
        The csvfiles and txtfiles are written to /output_data/<fastafilename/sequences/.
        The reason for fastafile mutants being written to an input folder is that these sequences are generated as
        direct inputs for the mutation operation.
        :param title_titleSeq_w_mutants: Title of wild-type associated to every mutant title:sequence.
        :param write_1_fasta_only: True to write one fastafile containing all mutants, separated by \n.
        :param write_fasta_per_mut: True to write one fastafile per mutant.
        :param path_output_3dots: Absolute path of output_data dir (where fasta, txt, csv written), includes to
        subdirs: /fastas/xxxx...yyyy/ e.g. 1001...2000
        :param write_csv: True to write 1 csv file for wt & mutants. False by default.
        :param write_txt: True to write 1 txt file for wt & mutants. False by default.
        """
        path_1_fastafile = None
        path_1_fastafile_open = None
        path_fastafilepermut = None
        path_seqscsv = None
        path_seqscsv_open = None
        path_seqstxt = None
        path_seqstxt_open = None

        for wt_title, title_seq in title_titleSeq_w_mutants.items():
            for mut_title, mut_seq in title_seq.items():
                if write_1_fasta_only and path_1_fastafile is None:
                    path_1_fastafile = GUM._os_makedirs(path_output_3dots, wt_title, Paths.DIR_MUTANTS.value)
                    path_1_fastafile = os.path.join(path_1_fastafile, wt_title + '_mutants.fasta')
                    path_1_fastafile_open = open(path_1_fastafile, 'w')
                    path_1_fastafile_open.write('>' + wt_title + '\n' + title_seq[wt_title] + '\n')
                if write_fasta_per_mut and path_fastafilepermut is None:
                    path_fastafilepermut = GUM._os_makedirs(path_output_3dots, Paths.DIR_FASTAS.value, wt_title)
                    path_fastafilepermut = os.path.join(path_fastafilepermut, wt_title + Str.FSTAEXT.value)
                    path_fastafilepermut_open = open(path_fastafilepermut, 'w')
                    path_fastafilepermut_open.write('>' + wt_title + '\n' + title_seq[wt_title] + '\n')
                    path_fastafilepermut_open.close()
                if write_csv and path_seqscsv is None:
                    path_seqscsv = GUM._os_makedirs(path_output_3dots, Paths.DIR_SEQS_TXT_CSV.value, wt_title)
                    path_seqscsv = os.path.join(path_seqscsv, wt_title + '_mutants.csv')
                    path_seqscsv_open = open(path_seqscsv, 'w')
                    path_seqscsv_open.write(wt_title + ':' + title_seq[wt_title] + ',')
                if write_txt and path_seqstxt is None:
                    path_seqstxt = GUM._os_makedirs(path_output_3dots, Paths.DIR_SEQS_TXT_CSV.value, wt_title)
                    path_seqstxt = os.path.join(path_seqstxt, wt_title + '_mutants.txt')
                    path_seqstxt_open = open(path_seqstxt, 'w')
                    path_seqstxt_open.write(wt_title + ':' + title_seq[wt_title] + '\n')
                elif mut_title is not wt_title:
                    if write_1_fasta_only and path_1_fastafile_open is not None:
                        path_1_fastafile_open.write('>' + mut_title + '\n' + mut_seq + '\n')
                    if write_fasta_per_mut and path_fastafilepermut is not None:
                        path_fastafilepermut = os.path.join(path_output_3dots, Paths.DIR_FASTAS.value, wt_title,
                                                            Paths.DIR_MUTANTS.value, mut_title + Str.FSTAEXT.value)
                        path_fastafilepermut_open = open(path_fastafilepermut, 'w')
                        path_fastafilepermut_open.write('>' + mut_title + '\n' + mut_seq + '\n')
                        path_fastafilepermut_open.close()
                    if write_csv and path_seqscsv_open is not None:
                        path_seqscsv_open.write(mut_title + ':' + mut_seq + ',')
                    if write_txt and path_seqstxt_open is not None:
                        path_seqstxt_open.write(mut_title + ':' + mut_seq + '\n')

        if path_1_fastafile_open is not None:
            path_1_fastafile_open.close()
        if path_seqstxt_open is not None:
            path_seqstxt_open.close()
        if path_seqscsv_open is not None:
            path_seqscsv_open.close()

    def make_titleSeqDict_from_fastafile(self, path_input_fastafiles):
        """
        Splits a fastafile's text contents into a dictionary. The key is fastafile's title, value is (amino acid)
        sequence. E.g. ">1_A\n'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVG...' etc" is saved as
            {
               '1_A': <seq of 1_A>,
               '2_A': <seq of 2_A>,
               etc..
            }
        The fasta file is assumed to be either two lines: 1st line ">title", 2nd line sequence, but the method also
        allows there to be either no title (i.e. sequence only) or an empty title (">" with no title after it). In these
        both scenarios, the title is created based on the name of the file itself.
        :param path_input_fastafiles: Str or List of Str, Absolute path to the fastafile(s) (incl. .fasta ext)
        e.g. ~/../input_data/fastas/1_A/1_A.fasta
        :return: title:sequence dictionary.
        """
        title_sequence_dict = {}
        if isinstance(path_input_fastafiles, str):
            path_input_fastafiles = [path_input_fastafiles]
        for path_fastafile in path_input_fastafiles:
            with open(path_fastafile, 'r') as path_fastafile_opened:
                fastafile_lines = path_fastafile_opened.readlines()
                if len(fastafile_lines) > 2 or len(fastafile_lines) < 1:
                    raise ValueError('There is either an extra unexpected carriage return or no sequence data at all')
                elif len(fastafile_lines) == 2:
                    title = GUM._remove_prefix_and_suffix(fastafile_lines[0],  ">", "\n")
                if title == "" or len(fastafile_lines) == 1:
                    title = os.path.splitext(path_input_fastafiles)[0]
                fasta_seq = fastafile_lines[1] if len(fastafile_lines) == 2 else fastafile_lines[0]
                title_sequence_dict[title] = fasta_seq
        return title_sequence_dict

    def convert_titleSeqDict_to_titleTitleSeqDictDict(self, title_sequence_dict: dict):
        """
        Converts a title:sequence dictionary to a title: {title:seq dictionary} dictionary. One use of this is to add
        19 more titles per residue for different titles - representing point mutants.
        E.g.
            {
              '1_A': <seq of 1_A>,
              '2_A': <seq of 2_A>,
               etc..
            }
            is converted to
            {
              '1_A':
                      {'1_A'      : <seq of 1_A>,
                       '1_A_R1A'  : <seq of 1_A_R1A>, ... etc },
              '2_A':
                      {'2_A'      : <seq of 2_A>,
                       '2_A_P1A'  : <seq of 2_A_P1A>, ... etc },
               etc..
            }
        :param title_sequence_dict: Fasta title (i.e. sequence id) is key. Sequence in FASTA format is value.
        :return: Dictionary of title: title:sequence dictionaries.
        """
        title_titleSeqDict_dict = {}
        for title in title_sequence_dict.keys():
            title_titleSeqDict_dict[title] = title_sequence_dict
        return title_titleSeqDict_dict

    def build_complete_paths_for_input_fastas(self, path_input, fastafiles):
        """
        Builds list of absolute path(s) for fastafile or fastafiles in list, by adding an extra directory between the
        path_input and the file itself, named exactly the same as the file.
        E.g. input_data/fastas/<fastaname>/fastafile.

        :param path_input: Absolute path of directory holding list of fasta files.
        :param fastafiles: Target fasta files (including extensions).
        :return: List of fasta files with absolute path to them, e.g. ~/.../input_data/fastas/1_A/1_A.fasta
        """
        if isinstance(fastafiles, str):
            fastafiles = [fastafiles]
        return [path_input + '/' + Paths.DIR_FASTAS.value + '/' + fastafile.split('.')[0] + '/' + fastafile for
                fastafile in fastafiles]
