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
import warnings


class Biopy(object):

    # Using alignments restricted to a fixed number (currently 20) and hitlist restricted to a fixed number (currently
    # 20), it takes about 21 seconds to run this qblast remotely.
    # However, using smaller number of alignments than the default (of 50), the output which is parsed to an xml file
    # is reduced in size from ~77 KB to ~53 KB. Any reduction could be important when you consider that there are
    # 78,586 fasta files to run (that's likely over 1.5 GB in total).
    # We are only interested in 100% identity hits, so unless a sequence has over 20 isoforms which all include that
    # sequence, the remaining hits on the blast_record are not used anyway.
    #
    # NOTE: occasionally the qblast took about 4 minutes.
    #
    # fasta_Str     String     The text of a fasta file including >title/n and sequence (standard FASTA format).
    #
    # Returns what biopython's NCBIWWW.qblast returns. (Data type is _io.TextIOWrapper.)
    @staticmethod
    def run_blastp(fasta_str):
        return NCBIWWW.qblast(program=Biopy.BlastParam.BLST_P.value,
                              database=Biopy.BlastParam.SWSPRT.value,
                              sequence=fasta_str,
                              entrez_query=Biopy.BlastParam.HOMSAP_ORG.value,
                              alignments=Biopy.BlastParam.MAX_ALIGN_20.value,
                              hitlist_size=Biopy.BlastParam.MAX_HIT_20.value)

    # Parses qblast result and filters (assigns to a data structure) only those fields that are of interest.
    #
    # NCBIXML.read() takes the _io.TextIOTextWrapper (which is created when xml file is opened).
    # NCBIXML.read() returns a Bio.Blast.Record.Blast object that has assigned various values in the qblast record as
    # fields that can be directly accessed.
    # NOTE: Time taken for this remote qblast for sequence 1_A.fasta was about 20 seconds.
    # NOTE: The query sequence id is also the filename and is used here for the name of the output xml.
    #
    # path_qblast_result    String      Abs path to the output xml file.
    # fastafile_name        String      Name of fastafile (without .fasta extension).
    # path_fastafile        String      Absolute path to fastafile (used downstream to spot any seq length anomalies).
    #
    # Returns a data structure of the elements of the qblast result that are considered here to be pertinent.
    @staticmethod
    def parse_filter_blastp_xml_to_dict(path_qblast, fastafile_name, path_fastafile):
        with open(path_qblast) as f:
            blast_record = NCBIXML.read(f)
            query_length = blast_record.query_length
            query_letters = blast_record.query_letters
            database_used = blast_record.database
            database_seqs_num = blast_record.database_sequences
            alignments = blast_record.alignments

        qblast_dict = {'query_seq_id': fastafile_name,
                        'query_length': query_length,
                        'database': database_used,
                        'database_seqs_num': database_seqs_num,
                        'identical_aligns_list': Biopy._make_list_of_dicts_of_hsps_with_0gaps_and_queryLen_equal_alignLen(
                        query_length, alignments)}
        Biopy._warn_if_discrepancies_in_query_seq_length(qblast_dict, query_length, query_letters, path_fastafile)
        return qblast_dict

    # Builds a list_of_dictionaries of relevant info for each high-scoring pair (hsp) alignment where there are zero
    # gaps and the alignment length is exactly the same as the query length, hence it is a 100% identity match.
    #
    # However, the hit fasta may still be longer than the query fasta. This info is stored in sbjct_start and sbjct_end.
    #
    # query_length      int     qblast output value, refers to length of the input fasta sequence.
    # alignments        List    Objects with values, (including hsps, also a list of objects with values) from qblast.
    #
    # Returns a list of identical alignment hits. The hit may still be longer.
    @staticmethod
    def _make_list_of_dicts_of_hsps_with_0gaps_and_queryLen_equal_alignLen(query_length, alignments):
        identical_aligns_list = []
        is_identical = False
        for alignment in alignments:
            alignment_dict = {'accession_num': 0, 'length': 0, 'hit_def': '', 'hsp_dict': {}}
            for hsp in alignment.hsps:
                is_identical = hsp.expect < 1e-50 and hsp.gaps == 0 and query_length == hsp.align_length == hsp.identities
                if is_identical:
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
            if is_identical:
                identical_aligns_list.append(alignment_dict)
        return identical_aligns_list

    # qblast_result_dict    Dictionary  Built of relevant elements from qblast return object.
    # query_length          int         qblast output value, refers to length of the input fasta sequence.
    # query_letters         int         qblast letters value, refers to number of input fasta sequence characters.
    # path_fastafile        String      FASTA file (with extension) and absolute path.
    # len_hit               int         length of protein that blast hit
    @staticmethod
    def _warn_if_discrepancies_in_query_seq_length(qblast_dict, query_length, query_letters, path_fastafile):
        for alignment in qblast_dict['identical_aligns_list']:
            len_sbjct_prot = alignment['hsp_dict']['sbjct_end'] - alignment['hsp_dict']['sbjct_start'] + 1

            if not query_length == query_letters:
                warnings.warn_explicit(message="Query_length of blast_record is not same as number of query letters in "
                                               "blast_record.",
                                       category=RuntimeWarning, filename="Biopy", lineno=65)
            if not query_length == len(SeqIO.read(path_fastafile, "fasta").seq):
                warnings.warn_explicit(message="Length of fasta input sequence not same as query length of blast_record.",
                                       category=RuntimeWarning, filename="Biopy", lineno=65)
            if query_length < len_sbjct_prot:
                warnings.warn_explicit(message="Length of fasta input (" + str(query_length) + ") is less than length of "
                                               "hit(" + str(len_sbjct_prot) + ").",
                                       category=RuntimeWarning, filename="Biopy", lineno=65)

    from enum import Enum

    # Created to hold constant values.
    # BLST_P        is qblast for protein sequences
    # MAX_ALIGN_..  Currently I have given a choice of using a max number of alignments to return from blast at either
    #               5 or 20
    # MAX_HIT_..    Currently I have given a choice of using a max number of hits to return from blast at either 5 or 20
    class BlastParam(Enum):
        BLST_P = "blastp"
        SWSPRT = "swissprot"
        HOMSAP_ORG = "Homo Sapiens[Organism]"
        MAX_ALIGN_5 = 5
        MAX_HIT_5 = 5
        MAX_ALIGN_20 = 20
        MAX_HIT_20 = 20
