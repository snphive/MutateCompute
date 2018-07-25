# Biopython packages NCBIWWW and NCBIXML are used here to perform a blast of query fasta sequences and to parse the
# resulting output xml file, respectively.
#######################################################################################################################
# NCBIXML parses the xml to an object called blast_record which has 49 attributes:
#
# blast_record
#   |
#   alignments
#   |
#   database
#   |
#   database_sequences
#   |
#   query_length
#   |
#   query_letters
#######################################################################################################################
# Nested in alignments are the following attributes:
#
# <Alignment> (unnamed element in zero-index-based list of alignments)
#   |
#   accession (string, identifier, e.g. Q99985)
#   |
#   hit_def (name of protein hit; alternative names)
#   |
#   hit_id (different identifiers for protein hit)
#   |
#   hsps (list of high scoring pairs)
#   |
#   length (full length of protein hit, not the query protein)
#   |
#   title (concatenation of hit_id and hit_def)
#######################################################################################################################
# Nested in each element of hsps list are 17 attributes, including the following:
#
# <HSP> (unnamed element in zero-index-based list of hsps)
#   |
#   align_length (int, length of this alignment.)
#   |
#   expect (float, 0-to-1 score that indicates how unlikely-to-likely this query sequence is to produce a random hit)
#   |
#   gaps (int, gaps in the alignment)
#   |
#   identities (int, number of identical residues)
#   |
#   positives (int, either identical or similar residues)
#   |
#   query_end (position of end of this alignment on query sequence)
#   |
#   query_start (position of start of this alignment on query sequence)
#   |
#   sbjct_end (position of end of this alignment on hit sequence)
#   |
#   sbjct_start (position of start of this alignment on hit sequence)
#
#######################################################################################################################
import os
import time
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from src.GeneralUtilityMethods import GUM


class Biopy(object):

    OUTPUT_LEAF = "blastp_swissprot"
    TWO_HUNDRED_KB = 200000

    @staticmethod
    def _read_input_fasta_seq(path_fasta_input_query_seq):
        with open(path_fasta_input_query_seq) as fasta_io:
            fasta_str = fasta_io.read()
        return fasta_str

    @staticmethod
    def _write_qblast_xml_result(qblast_result, path_output_root, output_leaf, output_filename):
        path_output = GUM.create_dir_tree(path_output_root, output_leaf)
        path_result_file = path_output + "/" + output_filename + ".xml"
        with open(path_result_file, 'w') as out_handle:
            out_handle.write(qblast_result.read())
            qblast_result.close()
        return path_result_file

    # NOTE: Time taken for this remote qblast for sequence 1_A.fasta was about 20 seconds.
    # NOTE: The query sequence id is also the filename and is used here for the name of the output xml.
    @staticmethod
    def find_identical_blastp_hit_swissprot_for_fasta_input(path_fasta_input_query_seq, path_output_root):
        fasta_input_seq = Biopy._read_input_fasta_seq(path_fasta_input_query_seq)
        qblast_result = Biopy._run_blastp(fasta_input_seq)
        filename = path_fasta_input_query_seq.split("/")[-1].split('.')[0]
        path_qblast_result = Biopy._write_qblast_xml_result(qblast_result, path_output_root, Biopy.OUTPUT_LEAF,
                                                            filename)
        minimum_filesize_for_parsing = Biopy.TWO_HUNDRED_KB
        qblast_xml_size = os.stat(path_qblast_result).st_size
        if qblast_xml_size > minimum_filesize_for_parsing:
            raise ValueError("Qblast xml output size is over 200 KB in size. Something may be wrong.")

        with open(path_qblast_result) as f:
            blast_record = NCBIXML.read(f)
            query_length = blast_record.query_length
            query_letters = blast_record.query_letters
            database_used = blast_record.database
            database_seqs_num = blast_record.database_sequences
            alignments = blast_record.alignments

        qblast_result_dict = {'query_seq_id': filename,
                              'query_length': query_length,
                              'database': database_used,
                              'database_seqs_num': database_seqs_num,
                              'alignment_dict': Biopy._make_alignmentDict_zeroGaps_queryLen_equal_alignLen(query_length,
                                                                                                           alignments)}
        Biopy.__print_discrepancies_in_query_sequence_length(query_length, query_letters, path_fasta_input_query_seq)
        return qblast_result_dict

    # Using alignments restricted to 5 and hitlist restricted to 5, the time qblast takes to run sequence 1_A.fasta
    # remotely is only reduced from 22 to 21 seconds. However, the output xml is reduced in size from 77 KB to 15 KB.
    # The remaining data blast_record are not used so nothing is lost by reducing the scope of the analysis. We are
    # only looking for 100 % identity anyway.
    # NOTE: occasionally the qblast took over 4 minutes!
    @staticmethod
    def _run_blastp(fasta_input_str):
        return NCBIWWW.qblast(program=Biopy.BlastParam.BLST_P.value,
                              database=Biopy.BlastParam.SWSPRT.value,
                              sequence=fasta_input_str,
                              entrez_query=Biopy.BlastParam.HOMSAP_ORG.value,
                              alignments=Biopy.BlastParam.MAX_ALIGN_20.value,
                              hitlist_size=Biopy.BlastParam.MAX_HIT_20.value)

    @staticmethod
    def _make_alignmentDict_zeroGaps_queryLen_equal_alignLen(query_length, alignments):
        alignment_dict = {'accession_num': 0, 'length': 0, 'hit_def': '', 'hsp_dict': {}}
        for alignment in alignments:
            for hsp in alignment.hsps:
                if hsp.expect == 0.0 and hsp.gaps == 0 and query_length == hsp.align_length:
                    alignment_dict['accession_num'] = alignment.accession
                    alignment_dict['length'] = alignment.length
                    alignment_dict['hit_def'] = alignment.hit_def
                    alignment_dict['hsp_dict']['align_length'] = hsp.align_length
                    alignment_dict['hsp_dict']['gaps'] = hsp.gaps
                    alignment_dict['hsp_dict']['identities'] = hsp.identities
                    alignment_dict['hsp_dict']['query_end'] = hsp.query_end
                    alignment_dict['hsp_dict']['query_start'] = hsp.query_start
                    alignment_dict['hsp_dict']['sbjct_end'] = hsp.sbjct_end
                    alignment_dict['hsp_dict']['sbjct_start'] = hsp.sbjct_start
        return alignment_dict

    @staticmethod
    def __print_discrepancies_in_query_sequence_length(query_length, query_letters, path_fasta_input_query_seq):
        if not query_length == query_letters:
            print("query_length of blast_record is not same as query letters of blast_record")

        if not query_length == len(SeqIO.read(path_fasta_input_query_seq, 'fasta').seq):
            print("Length of fasta input sequence not same as query length of blast_record")

# Note this enum class is created to hold constant values.
    from enum import Enum

    class BlastParam(Enum):
        BLST_P = "blastp"
        SWSPRT = "swissprot"
        HOMSAP_ORG = "Homo Sapiens[Organism]"
        MAX_ALIGN_5 = 5
        MAX_HIT_5 = 5
        MAX_ALIGN_20 = 20
        MAX_HIT_20 = 20
