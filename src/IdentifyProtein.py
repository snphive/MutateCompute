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
    # path_fastafile            String      Absolute path to FASTA input file (with .fasta extension).
    # write_blastp_json         Boolean     True to write the blastp result dictionary to a json output file.
    # path_output               String      Absolute path to output files. Typically output_data/<pdbname>/blastp.
    # write_idmaps_for_mysqldb  Boolean     True to build & write csv/xml files that map ids (accession and RvdK ids).
    @staticmethod
    def map_seq_to_swsprt_acc_id_and_write_files(path_fastafile, path_output, write_idmaps_for_mysldb, make_csv=True,
                                                 make_xml=True, make_json=False):
        with open(path_fastafile) as fasta_io:
            fasta_str = fasta_io.read()
        blastp_result = Biopy.run_blastp(fasta_str)
        fastafile_name = path_fastafile.split('/')[-1].split('.')[0]
        path_blstp_xml = IdProt._write_blast_xml(path_output, fastafile_name, blastp_result)
        blastp_dict = Biopy.parse_and_filter_blastp_xml_to_dict(path_blstp_xml, fastafile_name, path_fastafile)
        if write_idmaps_for_mysldb:
            IdProt._write_idmaps_for_mysqldb(path_output, blastp_dict, make_csv=make_csv, make_xml=make_xml,
                                             make_json=make_json)
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
    def _write_idmaps_for_mysqldb(path_output, blastp_dict, make_xml=True, make_csv=True, make_json=False):
        idmap = {IdProt.Strng.SEQ_ID.value: blastp_dict[IdProt.Strng.QRY_SEQ_ID.value]}
        fullname_search_str = IdProt.Strng.NAME_SRCH_STR.value
        altname_search_str = IdProt.Strng.ALTNAME_SRCH_STR.value
        flags_search_str = IdProt.Strng.FLGS_SRCH_STR.value
        for alignment in blastp_dict[IdProt.Strng.IDENT_ALIGN_LST.value]:
            idmap[IdProt.Strng.ACC_NUM.value] = alignment[IdProt.Strng.ACC_NUM.value]
            idmap[IdProt.Strng.LEN.value] = alignment[IdProt.Strng.LEN.value]
            full_alt_names_flags = IdProt._extract_names_flags(alignment[IdProt.Strng.HIT_DEF.value],
                                                            fullname_search_str, altname_search_str, flags_search_str)
            idmap[IdProt.Strng.NAME.value] = full_alt_names_flags[fullname_search_str]
            idmap[IdProt.Strng.ALTNAME.value] = full_alt_names_flags[altname_search_str]
            idmap[IdProt.Strng.FLGS.value] = full_alt_names_flags[flags_search_str]
            idmap[IdProt.Strng.STRT_POS.value] = alignment[IdProt.Strng.HSP_DICT.value][IdProt.Strng.SBJCT_STRT.value]
            idmap[IdProt.Strng.END_POS.value] = alignment[IdProt.Strng.HSP_DICT.value][IdProt.Strng.SBJCT_END.value]
        if make_csv:
            IdProt._write_idmap_csv(path_output, idmap)
        if make_xml:
            IdProt._write_idmap_xml(path_output, idmap)
        if make_json:
            IdProt._write_idmap_jsonfile(path_output, idmap)

    # Converts data passed in a dictionary data structure to csv format and writes it to the blastp directory.
    #
    # path_output   String          Absolute path to the output root directory for the destination dir of the csv file.
    # idmap         Dictionary      Data structure of values taken from blastp run that are parsed to a csv format.
    #
    # Returns the xml elements (with one attribute value) and newlines in their own element in a list.
    @staticmethod
    def _write_idmap_xml(path_output, idmap):
        xml_list = [IdProt.Strng.XML_PROLOG.value, '\n']
        for tagname, value in idmap.items():
            if tagname == IdProt.Strng.SEQ_ID.value:
                xml_list.append('<' + tagname + ' id=' + '"' + str(value) + '">')
                xml_list.append('\n')
            else:
                xml_list.append('\t<' + tagname + '>' + str(value) + '</' + tagname + '>')
                xml_list.append('\n')
        xml_list.append('</' + IdProt.Strng.SEQ_ID.value + '>')
        xml_list.append('\n')
        path_output_blastp_fastaname = IdProt._build_output_blastp_dir(path_output, idmap[IdProt.Strng.SEQ_ID.value])
        path_output_blastp_fastaname_xmlfile = os.path.join(path_output_blastp_fastaname, 'idmap_swsprt.xml')
        with open(path_output_blastp_fastaname_xmlfile, 'w') as idmap_file:
            xml_str = ''.join(xml_list)
            idmap_file.write(''.join(xml_str))
        return xml_list

    # Converts data passed in a dictionary data structure to csv format and writes it to the blastp directory.
    #
    # path_output   String          Absolute path to the output root directory for the destination dir of the csv file.
    # idmap         Dictionary      Data structure of values taken from blastp run that are parsed to a csv format.
    #
    # Returns the values and commas and newlines in their own element in a list.
    @staticmethod
    def _write_idmap_csv(path_output, idmap):
        csv_list = []
        for title in idmap.keys():
            csv_list.append(title)
            csv_list.append(',')
        csv_list = csv_list[:-1]
        csv_list.append('\n')
        for value in idmap.values():
            csv_list.append(str(value))
            csv_list.append(',')
        csv_list = csv_list[:-1]
        csv_list.append('\n')
        path_output_blastp_fastaname = IdProt._build_output_blastp_dir(path_output, idmap[IdProt.Strng.SEQ_ID.value])
        path_output_blastp_fastaname_csvfile = os.path.join(path_output_blastp_fastaname, 'idmap_swsprt.csv')
        with open(path_output_blastp_fastaname_csvfile, 'w') as idmap_file:
            csv_str = ''.join(csv_list)
            idmap_file.write(''.join(csv_str))
        return csv_list

    # Reduces code duplication. Despite very simple code, it is made slightly more verbose with the need for a try
    @staticmethod
    def _write_idmap_jsonfile(path_output, idmap):
        path_output_blastp_fastaname = IdProt._build_output_blastp_dir(path_output, idmap[IdProt.Strng.SEQ_ID.value])
        path_output_blastp_fastaname_jsonfile = os.path.join(path_output_blastp_fastaname, 'idmap_swsprt.json')
        with open(path_output_blastp_fastaname_jsonfile, 'w') as f:
            f.write(json.dumps(idmap))

    # except clause, making its encapsulation for reuse worthwhile.
    @staticmethod
    def _build_output_blastp_dir(path_output, fastadirname):
        path_output_blastp_fastadirname = os.path.join(path_output, Paths.DIR_BLASTP.value, fastadirname +
                                                       IdProt.Strng.DIR_SUFIX.value)
        try:
            os.makedirs(path_output_blastp_fastadirname)
        except FileExistsError:
            print('Part or all of path already exists. This is absolutely fine.')
        return path_output_blastp_fastadirname

    # NOTE: This method assumes that the 'Recname', 'Altname', 'Flags' are always written in this order. It also
    # assumes that the Flag is the last string of any form in the line.
    # I don't know yet if these assumptions are correct.
    @staticmethod
    def _extract_names_flags(hit_def, *names_flag):
        full_alt_names_flag = {}
        for search_str in names_flag:
            try:
                start_index = hit_def.index(search_str) + len(search_str)
                if search_str == IdProt.Strng.FLGS_SRCH_STR.value:
                    full_alt_names_flag[search_str] = hit_def[start_index:]
                    IdProt.__is_reasonable_length(len(full_alt_names_flag[search_str]))
                    break
                for i, char in enumerate(hit_def[start_index:], start=start_index):
                    if char == ';':
                        end_index = i
                        full_alt_names_flag[search_str] = hit_def[start_index:end_index]
                        IdProt.__is_reasonable_length(len(full_alt_names_flag[search_str]))
                        break
            except ValueError:
                print('Did not find ' + search_str + ' in this hit_def: ' + hit_def)
        return full_alt_names_flag

    @staticmethod
    def __is_reasonable_length(name_or_flag_length):
        reasonable_max_name_or_flag_len = 30
        if name_or_flag_length == 0 | name_or_flag_length > reasonable_max_name_or_flag_len:
            raise ValueError('Full name, altname or flag is either not there at all or seems too long.')
        return True

    from enum import Enum

    class Strng(Enum):
        SEQ_ID = 'sequence_id'
        QRY_SEQ_ID = 'query_seq_id'
        NAME_SRCH_STR = 'RecName: Full='
        ALTNAME_SRCH_STR = 'AltName: Full='
        FLGS_SRCH_STR = 'Flags: '
        IDENT_ALIGN_LST = 'identical_aligns_list'
        ACC_NUM = 'accession_num'
        LEN = 'length'
        HIT_DEF = 'hit_def'
        NAME = 'full_name'
        ALTNAME = 'altname'
        FLGS = 'flags'
        SBJCT_STRT = 'sbjct_start'
        SBJCT_END = 'sbjct_end'
        STRT_POS = 'start_pos'
        END_POS = 'end_pos'
        HSP_DICT = 'hsp_dict'
        XML_PROLOG = '<?xml version="1.0" encoding="UTF-8"?>'
        DIR_SUFIX = '_idmaps'
