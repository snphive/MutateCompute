import os
import json
from src.GeneralUtilityMethods import GUM
from src.Biopython import Biopy
from src.Paths import Paths


# Class to take any protein identifier, pdb file, FASTA, and identify the protein, potentially via Blast and/or then
# mapping whatever identifier the input has with the recognised identifiers (1_A to 32431_A from RvdK's repaired pdbs)
class IdProt(object):

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
        fastafile_name = path_fastafile.split('/')[-1].split('.')[0]
        path_blstp_resultXml = IdProt._write_blastp_xml_result(path_output, fastafile_name, blastp_result)
        blastp_result_dict = Biopy.parse_filter_blastp_xml_to_dict(path_blstp_resultXml, fastafile_name, path_fastafile)

        if write_blastp_json:
            IdProt._write_dict_to_json_file(path_output, fastafile_name, blastp_result_dict)
        if build_idmap_csv:
            id_map = IdProt._build_idmap_for_mysqldb(blastp_result_dict, fastafile_name)
            IdProt._write_idmap_csv(path_output, fastafile_name, id_map)
        return blastp_result_dict

    # Writes the blastp result file to an xml file in a output_data subdir called blastp.
    #
    # path_output       String              Abs path to output root dir, typically /output_data
    # fastafile_name    String              Used for naming the xml file.
    # blastp_result     io.TextIOWrapper    Biopython's blastp output. (Is a buffered text stream).
    #
    # Returns path to the newly-written xml file.
    @staticmethod
    def _write_blastp_xml_result(path_output, fastafile_name, blastp_result):
        path_output_blastp = IdProt._build_output_blastp_dir(path_output)
        path_output_blastp_xmlfile = os.path.join(path_output_blastp, fastafile_name + ".xml")
        with open(path_output_blastp_xmlfile, 'w') as out_handle:
            out_handle.write(blastp_result.read())
            blastp_result.close()

        if os.stat(path_output_blastp_xmlfile).st_size > IdProt.TWO_HUNDRED_KB:
            raise ValueError("blast xml output size is over 200 KB in size. Something may have gone wrong with the "
                             "blastp run.")

        return path_output_blastp_xmlfile

    @staticmethod
    def _write_dict_to_json_file(path_output, fastafile_name, blastp_result_dict):
        path_output_blastp = IdProt._build_output_blastp_dir(path_output)
        path_output_blastp_jsonfile = os.path.join(path_output_blastp, fastafile_name + '.json')
        with open(path_output_blastp_jsonfile, 'w') as f:
            f.write(json.dumps(blastp_result_dict))

    @staticmethod
    def _write_idmap_csv(path_output, filename, id_map):
        path_output_blastp = IdProt._build_output_blastp_dir(path_output)
        path_output_filename_blstpswp_csvfile = os.path.join(path_output_blastp, filename + '.csv')
        with open(path_output_filename_blstpswp_csvfile, 'w') as f:
            f.write(json.dumps(id_map))
        return path_output_filename_blstpswp_csvfile

    @staticmethod
    def _build_idmap_for_mysqldb(blastp_result_dict, filename):
        #  TODO
        # lots to do here to parse the dict or xml into a csv that can be uploaded to mysql table
        dict = {}
        return dict

    @staticmethod
    def _build_output_blastp_dir(path_output):
        path_output_blastp = os.path.join(path_output, Paths.DIR_BLASTP.value)
        try:
            os.makedirs(path_output_blastp)
        except FileExistsError:
            print('Part of all of path already exists. This is absolutely fine.')
        return path_output_blastp
