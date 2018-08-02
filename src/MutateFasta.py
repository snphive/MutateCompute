from src.GeneralUtilityMethods import GUM
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC


class MutateFasta(object):

    def __init__(self):
        print('hello world MutateFasta constructor')

    # Mutates FASTA (typically wild-type) sequences at every position, to every residue specified (typically all other
    # 19 residues.)
    def mutate_every_residue_fastas(self, path_input, fasta_list, mutant_aa_list):
        path_fasta_file_list = self._build_path_to_fasta_add_extension(path_input, fasta_list)
        fasta_title_sequence_dict_list = self._read_fasta_files_for_sequences(path_fasta_file_list)
        mutant_fasta_list = self._mutate_fasta_list(fasta_title_sequence_dict_list, mutant_aa_list)
        return mutant_fasta_list

    def _build_path_to_fasta_add_extension(self, path_input, fasta_list):
        path_fasta_file_list = []
        for fasta in fasta_list:
            path_fasta_file_list.append(path_input + '/' + fasta + '/' + fasta + '.fasta')
        return path_fasta_file_list

    # returns a list of fasta title:sequence dictionaries
    def _read_fasta_files_for_sequences(self, path_fasta_file_list):
        fasta_title_sequence_dict_list = []
        for path_fasta_file in path_fasta_file_list:
            fasta_title_sequence_dict_list.append(GUM.get_title_sequence_dict_from_fasta_file(path_fasta_file))
        return fasta_title_sequence_dict_list

    def _mutate_fasta_list(self, fasta_title_sequence_dict_list, mutant_aa_list):
        result_dict = {}
        for fasta_title_sequence_dict in fasta_title_sequence_dict_list:
            for title, sequence in fasta_title_sequence_dict.items():
                mutant_fasta_list = self._mutate_fasta(sequence, mutant_aa_list)
                result_dict[title] = mutant_fasta_list
                self._write_list_to_file(result_dict)

    def _mutate_fasta(self, sequence, mutant_aa_list):
        mutated_sequences_list = []
        # my_seq = Seq(sequence, IUPAC.protein)
        # mutable_seq = my_seq.tomutable()
        mutable_seq = MutableSeq(sequence, IUPAC.protein)
        for i, wt_aa in enumerate(mutable_seq):
            for mutant_aa in mutant_aa_list:
                wt_aa_at_i = mutable_seq[i]
                if not mutant_aa == wt_aa_at_i:
                    mutable_seq[i] = mutant_aa
                    mutated_sequences_list.append(mutable_seq)
                    mutable_seq[i] = wt_aa_at_i
        return mutated_sequences_list
