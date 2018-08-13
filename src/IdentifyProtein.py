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
    def map_seq_to_swsprt_acc_id_and_write_csv(path_fastafile, write_blastp_json, write_idmap_csv, path_output):
        with open(path_fastafile) as fasta_io:
            fasta_str = fasta_io.read()
        blastp_result = Biopy.run_blastp(fasta_str)
        fastafile_name = path_fastafile.split('/')[-1].split('.')[0]
        path_blstp_resultXml = IdProt._write_blast_xml(path_output, fastafile_name, blastp_result)
        blastp_dict = Biopy.parse_filter_blastp_xml_to_dict(path_blstp_resultXml, fastafile_name, path_fastafile)
        if write_blastp_json:
            IdProt._write_dict_to_json_file(path_output, fastafile_name, blastp_dict)
        if write_idmap_csv:
            IdProt._write_idmaps_for_mysqldb(blastp_dict)
        return blastp_dict


    # Writes the blastp result file to an xml file in a output_data subdir called blastp.
    #
    # path_output       String              Abs path to output root dir, typically /output_data
    # fastafile_name    String              Used for naming the xml file.
    # blastp_result     io.TextIOWrapper    Biopython's blastp output. (Is a buffered text stream).
    #
    # Returns path to the newly-written xml file.
    @staticmethod
    def _write_blast_xml(path_output, fastafile_name, blastp_result):
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
    def _write_idmaps_for_mysqldb(blast_dict):
        idmap = {'sequence_id': blast_dict['query_seq_id']}
        fullname_search_str = 'RecName: Full='
        altname_search_str = 'AltName: Full='
        flags_search_str = 'Flags: '
        for alignment in blast_dict['identical_aligns_list']:
            idmap['accession_num'] = alignment['accession_num']
            idmap['length'] = str(alignment['length'])
            full_alt_names_flags = IdProt.extract_names_flags(alignment['hit_def'], fullname_search_str,
                                                              altname_search_str, flags_search_str)
            idmap['full_name'] = full_alt_names_flags[fullname_search_str]
            idmap['altname'] = full_alt_names_flags[altname_search_str]
            idmap['flags'] = full_alt_names_flags[flags_search_str]
            idmap['sbjct_start'] = str(alignment['hsp_dict']['sbjct_start'])
            idmap['sbjct_end'] = str(alignment['hsp_dict']['sbjct_end'])
        IdProt._write_xml(idmap)
        IdProt._write_csv(idmap)

    @staticmethod
    def _write_xml(idmap):
        xml_list = ['<?xml version="1.0" encoding="UTF-8"?>']
        for tagname, value in idmap.items():
            if tagname == 'sequence_id':
                xml_list.append('<' + tagname + '=' + '"' + value + '">' + '</' + tagname + '>')
            else:
                xml_list.append('<' + tagname + '>' + value + '</' + tagname + '>')
        with open(Paths.MC_OUTPUT_BLASTP.value + '/idmap_swissprot.xml', 'w') as idmap_file:
            idmap_file.write(''.join(xml_list))
        return xml_list

    @staticmethod
    def _write_csv(idmap):
        csv_list = []
        for title in idmap.keys():
            csv_list.append(title + ',')
        csv_list = csv_list[:-1]
        for value in idmap.values():
            csv_list.append(value + ',')
        csv_list = csv_list[:-1]
        with open(Paths.MC_OUTPUT_BLASTP.value + '/idmap_swissprot.csv', 'w') as idmap_file:
            idmap_file.write(''.join(csv_list))
        return csv_list

    @staticmethod
    def extract_names_flags(hit_def, *args):
        full_alt_names_flags = {}
        reasonable_max_len_for_names_flags = 25
        for search_str in args:
            try:
                fullname_start_index = hit_def.index(search_str) + len(search_str)
                for i, c in enumerate(hit_def, start=fullname_start_index):
                    if c == ';':
                        semi_colon_index = i
                        full_alt_names_flags[search_str] = hit_def[fullname_start_index:semi_colon_index]
                        if semi_colon_index == fullname_start_index | semi_colon_index > fullname_start_index + \
                                reasonable_max_len_for_names_flags:
                            raise ValueError('The delimiting semi colon is either missing or is beyong the expected '
                                             'max length of: ' + search_str + ' - currently limited to 25 characters.')
            except ValueError:
                print('Did not find ' + search_str + ' in this hit_def: ' + hit_def)
        return full_alt_names_flags


    # Reduces code duplication. Despite very simple code, it is made slightly more verbose with the need for a try
    @staticmethod
    def _write_dict_to_json_file(path_output, fastafile_name, blastp_result_dict):
        path_output_blastp = IdProt._build_output_blastp_dir(path_output)
        path_output_blastp_jsonfile = os.path.join(path_output_blastp, fastafile_name + '.json')
        with open(path_output_blastp_jsonfile, 'w') as f:
            f.write(json.dumps(blastp_result_dict))

    # except clause, making its encapsulation for reuse worthwhile.
    @staticmethod
    def _build_output_blastp_dir(path_output):
        path_output_blastp = os.path.join(path_output, Paths.DIR_BLASTP.value)
        try:
            os.makedirs(path_output_blastp)
        except FileExistsError:
            print('Part of all of path already exists. This is absolutely fine.')
        return path_output_blastp