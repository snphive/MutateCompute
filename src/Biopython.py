# Biopython packages NCBIWWW and NCBIXML are used here to perform a blast of query fasta sequences and to parse the
# resulting output xml file, respectively.
#######################################################################################################################
# NCBIXML parses the xml to an object called blast_record which has 49 attributes:
#
# blast_record: alignments; database; database_sequences; query_length; query_letters
#######################################################################################################################
# Nested in alignments are the following attributes:
#
# <Alignment> (unnamed element in zero-index-based list of alignments)
#   accession (string, identifier, e.g. Q99985)
#   hit_def (name of protein hit; alternative names)
#   hit_id (different identifiers for protein hit)
#   hsps (list of high scoring pairs)
#   length (full length of protein hit, not the query protein)
#   title (concatenation of hit_id and hit_def)
#######################################################################################################################
# Nested in each element of hsps list are 17 attributes, including the following:
#
# <HSP> (unnamed element in zero-index-based list of hsps)
#   align_length (int, length of this alignment.)
#   expect (float, 0-to-1 score that indicates how unlikely-to-likely this query sequence is to produce a random hit)
#   gaps (int, gaps in the alignment)
#   identities (int, number of identical residues)
#   positives (int, either identical or similar residues)
#   query_end (position of end of this alignment on query sequence)
#   query_start (position of start of this alignment on query sequence)
#   sbjct_end (position of end of this alignment on hit sequence)
#   sbjct_start (position of start of this alignment on hit sequence)
#
#######################################################################################################################
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO


class Biopy(object):

    # Using alignments restricted to 5 and hitlist restricted to 5, the time qblast takes to run sequence 1_A.fasta
    # remotely is only reduced from 22 to 21 seconds. However, the output xml is reduced in size from 77 KB to 15 KB.
    # The remaining data blast_record are not used so nothing is lost by reducing the scope of the analysis. We are
    # only looking for 100 % identity anyway. NOTE: occasionally the qblast took over 4 minutes!
    @staticmethod
    def run_blastp(fasta_input_str):
        return NCBIWWW.qblast(program=Biopy.BlastParam.BLST_P.value,
                              database=Biopy.BlastParam.SWSPRT.value,
                              sequence=fasta_input_str,
                              entrez_query=Biopy.BlastParam.HOMSAP_ORG.value,
                              alignments=Biopy.BlastParam.MAX_ALIGN_20.value,
                              hitlist_size=Biopy.BlastParam.MAX_HIT_20.value)

    # NOTE: Time taken for this remote qblast for sequence 1_A.fasta was about 20 seconds.
    # NOTE: The query sequence id is also the filename and is used here for the name of the output xml.
    @staticmethod
    def parse_filter_blastp_xml_to_dict(path_qblast_result, filename, path_fasta_file):
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
                              'alignment_dict': Biopy._make_dict_of_hsps_with_0Gaps_queryLen_is_alignLen(query_length,
                                                                                                        alignments)}
        len_sbjct_prot = qblast_result_dict['alignment_dict']['hsp_dict']['sbjct_end'] - \
                         qblast_result_dict['alignment_dict']['hsp_dict']['sbjct_start']

        Biopy.__print_discrepancies(query_length, query_letters, path_fasta_file, len_sbjct_prot)
        return qblast_result_dict

    # Builds a dictionary of relevant info for each high-scoring pair (hsp) alignment where there are zero gaps and
    # the alignment length is exactly the same as the query length, hence it is a 100% identity match. However, the hit
    # fasta may still be longer than the query fasta. This info is stored in sbjct_start and sbjct_end.
    # query_length      int     qblast output value, refers to length of the input fasta sequence.
    # alignments        String ? can find out what Type this is with Python)*******
    @staticmethod
    def _make_dict_of_hsps_with_0Gaps_queryLen_is_alignLen(query_length, alignments):
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

    # query_length      int     qblast output value, refers to length of the input fasta sequence.
    # query_letters     int     qblast letters value, refers to number of input fasta sequence characters.
    # path_fasta_file   String  FASTA file (with extension) and absolute path.
    # len_hit           int     length of protein that blast hit
    @staticmethod
    def __print_discrepancies_in_query_sequence_length(query_length, query_letters, path_fasta_file, len_hit):
        if not query_length == query_letters:
            print("Warning: query_length of blast_record is not same as query letters of blast_record.")
        if not query_length == len(SeqIO.read(path_fasta_file, 'fasta').seq):
            print("Warning: length of fasta input sequence not same as query length of blast_record.")
        if not query_length < len_hit:
            print("NB: length of fasta input (" + str(query_length) + ") is less than length of hit(" + str(len_hit))

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
