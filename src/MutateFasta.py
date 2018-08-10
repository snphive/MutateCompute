import os
from src.GeneralUtilityMethods import GUM
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC


class MutateFasta(object):

    # def __init__(self, path_input, fastafile_list, mutant_aa_list):
    #     mutants = self.mutate_every_residue_in_fasta_list(path_input, fastafile_list, mutant_aa_list)
    #     self._write_mutants_to_file(fastamutants, path_input, make_fastafile_for_all_mutants,
    #                            make_fastafile_per_mutant,  path_output, make_csv_file=True, make_txt_file=True)

    def __init__(self):
        print('')

    # Mutates FASTA (typically wild-type) sequences at every position, to every residue specified (typically all other
    # 19 residues.)
    #
    # path_input        String
    # fastafile_list    List
    # mutant_aa_list    List
    #
    # Returns ..
    def mutate_every_residue_in_fasta_list(self, path_input, fastafile_list, mutant_aa_list):
        path_fastafile_list = GUM.build_complete_paths_for_fastafiles(path_input, fastafile_list)
        titleSeqDict = GUM.make_titleSeqDict_from_fastafile(path_fastafile_list)
        titleTitleSeqDictDict = GUM.convert_titleSeqDict_to_titleTitleSeqDictDict(titleSeqDict)
        titleTitleSeqDictDict_w_mutants = self._populate_titleTitleSeqDictDict_with_mutants(
            titleTitleSeqDictDict, mutant_aa_list)
        return titleTitleSeqDictDict_w_mutants

    # titleTitleSeqDictDict   Dictionary        Fasta title is key. Value is mutant-title: mutated-sequence dictionary.
    #                                           This already has the dict of wt-title: (wt-title:wt-sequence) dict.
    # mutant_aa_list                List        All residues that the sequence(s) should be mutated to.
    #
    # Returns wt-title: mutant-title:mutant-sequence
    def _populate_titleTitleSeqDictDict_with_mutants(self, titleTitleSeqDictDict, mutant_aa_list):
        titleTitleSeqDictDicts_w_mutants = titleTitleSeqDictDict
        for wt_title in list(titleTitleSeqDictDict.keys()):
            wt_seq = titleTitleSeqDictDict[wt_title][wt_title]
            titleTitleSeqDictDicts_w_mutants[wt_title] = self._add_mutantTitle_mutatedSeq_to_dict(
                titleTitleSeqDictDict[wt_title], wt_title, wt_seq, mutant_aa_list)
        return titleTitleSeqDictDicts_w_mutants

    # A sequence is mutated at every position to every residue that is included in the list (mutant_aa_list).
    #
    # E.g.
    # input as {'230498_A': 'MDVFM'}
    # output as {'230498_A': 'MDVFM','M1A': 'ADVFM','M1C': 'CDVFM', etc.. }
    #
    # titleSeqDict   Dictionary      Already should have wt-title: wt-title:wt-sequence dictionary.
    # wt_seq                  String
    # mutant_aa_list                List
    #
    def _add_mutantTitle_mutatedSeq_to_dict(self, titleSeqDict, wt_title,  wt_seq, mutant_aa_list):
        titleSeqDict_w_mutants = {wt_title: titleSeqDict[wt_title]}
        mutable_seq = MutableSeq(wt_seq, IUPAC.protein)
        for i, wt_aa in enumerate(mutable_seq):
            for mutant_aa in mutant_aa_list:
                mutant_title = wt_aa + str(i + 1) + mutant_aa
                wt_aa_at_i = mutable_seq[i]
                if not mutant_aa == wt_aa_at_i:
                    mutable_seq[i] = mutant_aa
                    titleSeqDict_w_mutants[mutant_title] = str(mutable_seq)
                    mutable_seq[i] = wt_aa_at_i
        return titleSeqDict_w_mutants

    # Writes the mutants out either to fastafiles. There is the option to write each mutant out to its own fastafile,
    # and/or all mutants to one fastafile, and/or all mutants to one csv file and/or all mutants to a txt file.
    #
    # title_titleSeq_w_mutants      Dictionary  Title of wild-type associated to every mutant title:sequence
    #                                           starting with the wild-type title:sequence pair.
    # path_input                    String      Absolute path of input_data dir (where fastafiles will be written).
    # make_fastafile_all_mutants    Boolean     True to write one fastafile containing all mutants, separated by \n.
    # make_fastafile_per_mutant     Boolean     True to write one fastafile per mutant.
    # path_output                   String      Absolute path of input_data dir (where txt & csv seqs will be written).
    # make_csv_file                 Boolean     True to write one csv file for all mutants including wild-type.
    # make_txt_file                 Boolean     True to write one txt file for all mutants including wild-type.
    def _write_mutants_to_file(self, title_titleSeq_w_mutants, path_input, make_fastafile_all_mutants,
                               make_fastafile_per_mutant, path_output, make_csv_file=False, make_txt_file=False):
        path_all_mutants_fastafile = None
        path_seqs_w_mutants = None
        for wt_title, title_seq in title_titleSeq_w_mutants.items():
            path_input_fastas_title_mutants = os.path.join(path_input, 'fastas', wt_title, 'mutants')
            try:
                os.makedirs(path_input_fastas_title_mutants)
            except FileExistsError:
                print('Part of all of path already exists. This is absolutely fine.')
            if make_fastafile_all_mutants:
                path_all_mutants_fastafile = open(os.path.join(path_input_fastas_title_mutants,
                                                               wt_title + '_mutants.fasta'), 'w')
            if make_csv_file or make_txt_file:
                path_output_title_seqs = os.path.join(path_output, wt_title, 'sequences')
                try:
                    os.makedirs(path_output_title_seqs)
                except FileExistsError:
                    print('Part of all of path already exists. This is absolutely fine.')
            if make_csv_file:
                path_seqs_w_mutants = open(os.path.join(path_output_title_seqs, wt_title + '_mutants.csv'), 'w')
                path_seqs_w_mutants.write(wt_title + ':' + title_seq[wt_title] + ',')
            if make_txt_file:
                path_seqs_w_mutants = open(os.path.join(path_output_title_seqs, wt_title + '_mutants.txt'), 'w')
                path_seqs_w_mutants.write(wt_title + ':' + title_seq[wt_title] + '\n')
            for mut_title, mut_seq in title_seq.items():
                if mut_title is not wt_title:
                    if path_all_mutants_fastafile is not None:
                        path_all_mutants_fastafile.write('>' + mut_title + '\n' + mut_seq + '\n')
                    if path_seqs_w_mutants is not None:
                        if make_csv_file:
                            path_seqs_w_mutants.write(mut_title + ':' + mut_seq + ',')
                        if make_txt_file:
                            path_seqs_w_mutants.write(mut_title + ':' + mut_seq + '\n')
                    if make_fastafile_per_mutant:
                        path_mutant_fastafile = open(os.path.join(path_input_fastas_title_mutants,
                                                                  mut_title + '.fasta'), 'w')
                        path_mutant_fastafile.write('>' + mut_title + '\n' + mut_seq)
                        path_mutant_fastafile.close()
            if path_all_mutants_fastafile is not None:
                path_all_mutants_fastafile.close()
            if path_seqs_w_mutants is not None:
                path_seqs_w_mutants.close()
