import os
import json
from src.Biopython import Biopy
from src.Paths import Paths


# Class to take any protein identifier, pdb file, FASTA, and identify the protein, potentially via Blast and/or then
# mapping whatever identifier the input has with the recognised identifiers (1_A to 32431_A from RvdK's repaired pdbs)
class IdProt(object):

    # (This method relies on the presence of fastafiles in the specified folder, in order for them to be run on Blastp,
    # This transfer is currently done manually)
    #
    # path_fastafile            String      Absolute path to FASTA input file (with .fasta extension).
    # path_output               String      Absolute path to output files. Typically output_data/<pdbname>/blastp.
    # write_idmaps_for_mysqldb  Boolean     True to build dict mapping RvdK ids to swsprt acc & write files.
    # write_csv                 Boolean     True to write csvfiles.
    # write_xml                 Boolean     True to write xmlfiles.
    # write_json                Boolean     True to write jsonfiles.
    #
    # Returns the dictionary data structure representation of the parsed & filtered Blastp run result.
    @staticmethod
    def map_seq_to_swsprt_acc_id_and_write_files(path_fastafile, path_output, write_idmaps_for_mysldb, write_csv=True,
                                                 write_xml=True, write_json=False):
        with open(path_fastafile) as fasta_io:
            fasta_str = fasta_io.read()
        blastp_result = Biopy.run_blastp(fasta_str)
        fastafile_name = path_fastafile.split('/')[-1].split('.')[0]
        path_blstp_xml = IdProt._write_blast_xml(path_output, fastafile_name, blastp_result)
        blastp_dict = Biopy.parse_and_filter_blastp_xml_to_dict(path_blstp_xml, fastafile_name, path_fastafile)
        if write_idmaps_for_mysldb:
            IdProt._write_idmaps_for_mysqldb(path_output, blastp_dict, write_csv=write_csv, write_xml=write_xml,
                                             write_json=write_json)
        return blastp_dict

    # Writes the blastp result file to an xml file in a output_data subdir called blastp.
    #
    # path_output       String              Absolute path to output root dir, typically /output_data.
    # fastafile_name    String              Used for naming the xml file.
    # blastp_result     io.TextIOWrapper    Biopython's raw blastp output. (Is a buffered text stream).
    #
    # Returns path to the newly-written xml file.
    @staticmethod
    def _write_blast_xml(path_output, fastafile_name, blastp_result):
        path_output_blastp = IdProt._build_output_blastp_dir(path_output)
        path_output_blastp_xmlfile = os.path.join(path_output_blastp, fastafile_name + ".xml")
        with open(path_output_blastp_xmlfile, 'w') as out_handle:
            out_handle.write(blastp_result.read())
            blastp_result.close()

        if os.stat(path_output_blastp_xmlfile).st_size > IdProt.Strng.EST_MAX_SZE_BLSTP_RUN_200_KB:
            raise ValueError("blast xml output size is over 200 KB in size. Something may have gone wrong with the "
                             "blastp run.")

        return path_output_blastp_xmlfile

    # Write to file the result of blastp run (that has already been parsed & filtered from the raw blastp result).
    #
    # path_output   String          Absolute path to the output-data dir.
    # blastp_dict   Dictionary      Data structure parsed and filtered from Blastp run.
    # write_xml     Boolean         True to write xmlfile represetation of blastp_dict, True by default.
    # write_csv     Boolean         True to write csvfile represetation of blastp_dict, True by default.
    # write_json    Boolean         True to write jsonfile represetation of blastp_dict, False by default.
    @staticmethod
    def _write_idmaps_for_mysqldb(path_output, blastp_dict, write_xml=True, write_csv=True, write_json=False):
        idmap = {IdProt.Strng.SEQ_ID.value: blastp_dict[IdProt.Strng.QRY_SEQ_ID.value]}
        name_search_str = IdProt.Strng.NAME_SRCH_STR.value
        altname_search_str = IdProt.Strng.ALTNAME_SRCH_STR.value
        flags_search_str = IdProt.Strng.FLGS_SRCH_STR.value
        for alignment in blastp_dict[IdProt.Strng.IDENT_ALIGN_LST.value]:
            idmap[IdProt.Strng.ACC_NUM.value] = alignment[IdProt.Strng.ACC_NUM.value]
            idmap[IdProt.Strng.LEN.value] = alignment[IdProt.Strng.LEN.value]
            full_alt_names_flags = IdProt._extract_names_flags(alignment[IdProt.Strng.HIT_DEF.value], name_search_str,
                                                               altname_search_str, flags_search_str)
            idmap[IdProt.Strng.NAME.value] = full_alt_names_flags[name_search_str]
            idmap[IdProt.Strng.ALTNAME.value] = full_alt_names_flags[altname_search_str]
            idmap[IdProt.Strng.FLGS.value] = full_alt_names_flags[flags_search_str]
            idmap[IdProt.Strng.STRT_POS.value] = alignment[IdProt.Strng.HSP_DICT.value][IdProt.Strng.SBJCT_STRT.value]
            idmap[IdProt.Strng.END_POS.value] = alignment[IdProt.Strng.HSP_DICT.value][IdProt.Strng.SBJCT_END.value]
        if write_csv:
            IdProt._write_idmap_csv(path_output, idmap)
        if write_xml:
            IdProt._write_idmap_xml(path_output, idmap)
        if write_json:
            IdProt._write_idmap_jsonfile(path_output, idmap)

    # Converts data passed in a dictionary data structure to csv format and writes it to the blastp directory.
    #
    # path_output   String          Abs path of output_data dir, destination for xml file.
    # idmap         Dictionary      Data structure of values taken from blastp run that are parsed to a csv format.
    #
    # Returns the xml elements (with one attribute value for id) and newlines in their own element in a list.
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
    # path_output   String          Abs path of output_data dir, destination for csv file.
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

    # Writes out a json file of the idmap (dictionary) using an imported library method from 'json' module.
    #
    # path_output   String      Abs path of output_data dir, destination for json file.
    # idmap         Dictionary  Data structure representation of parsed & filtered Blastp runnresult.
    @staticmethod
    def _write_idmap_jsonfile(path_output, idmap):
        path_output_blastp_fastaname = IdProt._build_output_blastp_dir(path_output, idmap[IdProt.Strng.SEQ_ID.value])
        path_output_blastp_fastaname_jsonfile = os.path.join(path_output_blastp_fastaname, 'idmap_swsprt.json')
        with open(path_output_blastp_fastaname_jsonfile, 'w') as f:
            f.write(json.dumps(idmap))

    # To make directory and all intermediate directories, encapsulated in try/except clause as pre-existing result in
    # exceptions being thrown.
    #
    # path_output   String      Abs path of the output_data root directory of tree to be built here.
    # fastadirname  String      Name of fastafile (without extension)
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
    #
    # hit_def           String      A value from Blastp result, typically includes semi-colon separated text regarding
    #                               names and flags.
    # *search_strngs    List        Search strings, typically name, altname, flags.
    @staticmethod
    def _extract_names_flags(hit_def, *search_strngs):
        full_alt_names_flag = {}
        for search_str in search_strngs:
            try:
                start_index = hit_def.index(search_str) + len(search_str)
                if search_str == IdProt.Strng.FLGS_SRCH_STR.value:
                    full_alt_names_flag[search_str] = hit_def[start_index:]
                    IdProt.__is_reasonable_length(full_alt_names_flag[search_str])
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

    # name_or_flag      String      Text to be tested for length
    #
    # Returns True if text passed is between 2 and 30 characters long. I cannot say for sure that there are no names
    # or flags that are 1 character long or more than 30 characters long - this is running assumption.
    @staticmethod
    def __is_reasonable_length(name_or_flag):
        reasonable_max_name_or_flag_len = 30
        if len(name_or_flag) < 2 | len(name_or_flag) > reasonable_max_name_or_flag_len:
            raise ValueError('Full name, altname or flag is either not there at all or seems too long.')
        return True

    from enum import Enum

    class Strng(Enum):
        EST_MAX_SZE_BLSTP_RUN_200_KB = 200000
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
