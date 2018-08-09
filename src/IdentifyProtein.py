import os
import json
from src.GeneralUtilityMethods import GUM
from src.Biopython import Biopy


# Class to take any protein identifier, pdb file, FASTA, and identify the protein, potentially via Blast and/or then
# mapping whatever identifier the input has with the recognised identifiers (1_A to 32431_A from RvdK's repaired pdbs)
class IdProt(object):

    DIR_BLASTP = 'blastp'
    PATH_BLASTP = '/' + DIR_BLASTP
    TWO_HUNDRED_KB = 200000

    # (Manually transfer fasta files to an input folder at the path set as a constant global variable method.)
    #
    #
    # path_fastafile        String      Absolute path to FASTA input file (with .fasta extension).
    # write_blastp_json     Boolean     True to write the blastp result dictionary to a json output file.
    # build_idmap_csv       Boolean     True to build a csv file that maps ids (accession and RvdK ids).
    # path_output           String      Absolute path to output files. Typically output_data/<pdbname>/blastp.
    @staticmethod
    def start(path_fastafile, write_blastp_json, build_idmap_csv, path_output):
        with open(path_fastafile) as fasta_io:
            fasta_str = fasta_io.read()
        blastp_result = Biopy.run_blastp(fasta_str)
        print('******************#########***********************' + str(type(blastp_result)))
        fastafile_name = path_fastafile.split('/')[-1].split('.')[0]
        path_blstp_resultXml = IdProt._write_blastp_xml_result(path_output, fastafile_name, blastp_result)
        blastp_result_dict = Biopy.parse_filter_blastp_xml_to_dict(path_blstp_resultXml, fastafile_name, path_fastafile)
        if write_blastp_json:
            IdProt._write_dict_to_json_file(path_output, fastafile_name, blastp_result_dict)
        if build_idmap_csv:
            id_map = IdProt._build_idmap(blastp_result_dict, fastafile_name)
            IdProt._write_idmap_csv(path_output, fastafile_name, id_map)
        return blastp_result_dict

    # path_output           String
    # filename              String
    # blastp_result         ?
    @staticmethod
    def _write_blastp_xml_result(path_output, filename, blastp_result):
        path_output = GUM.create_dir_tree(path_output, filename, IdProt.DIR_BLASTP)
        path_result_file = path_output + "/" + filename + ".xml"
        with open(path_result_file, 'w') as out_handle:
            out_handle.write(blastp_result.read())
            blastp_result.close()
        if os.stat(path_result_file).st_size > IdProt.TWO_HUNDRED_KB:
            raise ValueError("blast xml output size is over 200 KB in size. Something may be wrong.")
        return path_result_file

    @staticmethod
    def _write_dict_to_json_file(path_output, filename, blastp_result_dict):
        path_json_file = path_output + "/" + filename + '/' + IdProt.DIR_BLAST_SWISPROT + "/" + filename + '.json'
        with open(path_json_file, 'w') as f:
            f.write(json.dumps(blastp_result_dict))

    @staticmethod
    def _build_idmap(blastp_result_dict, filename):
        #  TODO
        # lots to do here to parse the dict or xml into a csv that can be uploaded to mysql table
        dict = {}
        return dict

    @staticmethod
    def _write_idmap_csv(path_output, filename, id_map):
        path_output_filename_blstpswp_csvfile = path_output + "/" + filename + '/' + IdProt.DIR_BLAST_SWISPROT + "/" + filename + '.csv'
        with open(path_output_filename_blstpswp_csvfile, 'w') as f:
            f.write(json.dumps(id_map))
        return path_output_filename_blstpswp_csvfile
