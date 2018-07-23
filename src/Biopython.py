import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from src.GeneralUtilityMethods import GUM


class Biopy(object):

    PATH_OUTPUT_LEAF = 'blastp_swissprot'
    MILLION_BYTES = 1000000

    @staticmethod
    def _read_input_fasta_seq(path_fasta_input_query_seq):
        with open(path_fasta_input_query_seq) as fasta_io:
            fasta_str = fasta_io.read()
        return fasta_str

    @staticmethod
    def _write_qblast_xml_result(qblast_result, path_output_root, path_output_leaf, output_filename):
        path_output = GUM.create_dir_tree(path_output_root, path_output_leaf)
        path_result_file = path_output + '/' + output_filename + '.xml'
        with open(path_result_file, 'w') as out_handle:
            out_handle.write(qblast_result.read())
            qblast_result.close()
        return path_result_file

    # NOTE: time taken for this remote qblast for sequence 1_A.fasta was about 20 seconds.
    @staticmethod
    def find_identical_blastp_hit_swissprot_for_fasta_input(path_fasta_input_query_seq, path_output_root):
        fasta_input_seq = Biopy._read_input_fasta_seq(path_fasta_input_query_seq)
        qblast_result = Biopy._run_blastp(fasta_input_seq)
        output_filename = path_fasta_input_query_seq.split('/')[-1].split('.')[0]
        path_qblast_result = Biopy._write_qblast_xml_result(qblast_result, path_output_root, Biopy.PATH_OUTPUT_LEAF,
                                                              output_filename)
        minimum_filesize_for_parsing = Biopy.MILLION_BYTES

        qblast_xml_size = os.stat(path_qblast_result).st_size

        if qblast_xml_size < minimum_filesize_for_parsing:
            with open(path_qblast_result) as f:
                blast_record = NCBIXML.read(f)
        else:
            with open(path_qblast_result) as f:
                blast_records = NCBIXML.parse(f)

        if not blast_record.BlastOutput.BlastOutput_query-len == len(SeqIO.read(fasta_input_seq, 'fasta').seq):
            raise ValueError('The length of the fasta input sequence is different from the sequence length output by Blast - Something has gone wrong')


        for hit in blast_record.Hit:
            for hsp in hit.Hsp:
                print('Hsp_identity: ' + hsp.Hsp_identity)
                print('Hsp_gaps: ' + hsp.Hsp_gaps)

        blast_record.close()
        blast_records.close()
        return 'bla'

    @staticmethod
    def _run_blastp(fasta_input_seq):
        return NCBIWWW.qblast(program=Biopy.BlastParam.BLST_P.value,
                              database=Biopy.BlastParam.SWSPRT.value,
                              sequence=fasta_input_seq,
                              entrez_query=Biopy.BlastParam.HOMSAP_ORG.value)

    from enum import Enum

    class BlastParam(Enum):
        BLST_P = 'blastp'
        SWSPRT = 'swissprot'
        HOMSAP_ORG = 'Homo Sapiens[Organism]'

