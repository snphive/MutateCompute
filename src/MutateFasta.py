from src.GeneralUtilityMethods import GUM
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC


class MutateFasta(object):

    def __init__(self):
        mutants = self.mutate_every_residue_in_fasta_list()
        self._write_mutants_to_file(mutants)

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
        titleSeqDict_w_mutants = {}
        titleSeqDict_w_mutants[wt_title] = titleSeqDict[wt_title]
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

    # TODO
    # Writes the mutants data structure out to a csv file and txt file.
    def _write_mutants_to_file(self, mutants):
        print('')
