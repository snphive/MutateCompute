import os
import json
from src.Biopython import Biopy
from src.Paths import Paths
from src.GeneralUtilityMethods import GUM
import glob
from src.Cluster import Cluster


# Class to take any protein identifier, pdb file, FASTA, and identify the protein, potentially via Blast and/or then
# mapping whatever identifier the input has with the recognised identifiers (1_A to 32431_A from RvdK's repaired pdbs)
class IdProt(object):

    # Expects a directory location of fastafiles (not a fastafile itself).
    # (This method relies on the presence of fastafiles in the specified dir, in order for them to be run on Blastp.
    # This transfer is currently done manually)
    #
    # path_input_fastas         String      Abs path of root dir for input fastafiles (e.g. /input_data/fastas_10).
    # path_output               String      Abs path of root dir for output blastp files (..../output_data/).
    # write_idmaps_for_mysqldb  Boolean     True to build dict mapping RvdK ids to swsprt acc & write files.
    # write_csv                 Boolean     True to write csvfiles.
    # write_xml                 Boolean     True to write xmlfiles.
    # write_json                Boolean     True to write jsonfiles.
    #
    # Returns a list of dictionary data structure representations of each parsed & filtered Blastp run result.
    @staticmethod
    def map_seq_to_swsprt_acc_id_and_write_files(path_input_fastas_dir, use_cluster, path_output,
                                                 write_idmaps_for_mysldb, write_csv=True, write_xml=True,
                                                 write_json=False):
        try:
            path_input_fastas_dir_list = path_input_fastas_dir.split('/')
            if '.fasta' in path_input_fastas_dir_list[-1]:
                path_input_fastas_dir_list = path_input_fastas_dir_list[:-1]
                path_input_fastas_dir = '/'.join(path_input_fastas_dir_list)
            if not os.listdir(path_input_fastas_dir):
                raise ValueError('There are no subdirectories in your specified /input_data/fastas/ directory. '
                                 'Typically expecting /input_data/fastas/<fastafilenames>/, each containing a '
                                 'fastafile. It seems that no fastafiles were copied over from the repository '
                                 'directory. Therefore cannot proceed.')
        except NotADirectoryError:
            raise NotADirectoryError('The input path received: ' + path_input_fastas_dir + ' is is not a directory.')
        path_input_fastafile_list = glob.glob(path_input_fastas_dir + '/*.fasta')
        blastp_dict_list = []
        for path_fastafile in path_input_fastafile_list:
            with open(path_fastafile) as fastafile_opened:
                fasta_str = fastafile_opened.read()
                fastafile_name = path_fastafile.split('/')[-1].split('.')[0]
                if use_cluster:
                    Cluster.write_job_q_bash(job_name=fastafile_name, path_job_q_dir=Paths.CONFIG_JOBQ)
                    Cluster.run_job_q()
                blastp_result = Biopy.run_blastp(fasta_str)
                xml = IdProt._write_raw_blast_xml(path_output, fastafile_name, blastp_result)
                path_blstp_xml = xml
                blastp_dict = Biopy.parse_filter_blastp_xml_to_dict(path_blstp_xml, fastafile_name, path_fastafile)
                blastp_dict_list.append(blastp_dict)
                if write_idmaps_for_mysldb:
                    IdProt._write_idmaps_for_mysqldb(path_output, blastp_dict, write_csv=write_csv,
                                                     write_xml=write_xml, write_json=write_json)
        return blastp_dict_list

    # Writes the blastp result file to an xml file in a output_data subdir called blastp.
    #
    # path_output       String              Absolute path to output root dir, typically /output_data.
    # fastafile_name    String              Used for naming the xml file.
    # blastp_result     io.TextIOWrapper    Biopython's raw blastp output. (Is a buffered text stream).
    #
    # Returns path to the newly-written xml file.
    @staticmethod
    def _write_raw_blast_xml(path_output, fastafile_name, blastp_result):
        path_output_blastp = IdProt._build_dir_tree_with_intermed_dir(path_root=path_output,
                                                                      intermed_dir=Paths.DIR_BLASTP.value,
                                                                      fastadir=fastafile_name)
        path_output_blastp_xmlfile = os.path.join(path_output_blastp, fastafile_name + ".xml")
        with open(path_output_blastp_xmlfile, 'w') as out_handle:
            out_handle.write(blastp_result.read())
            blastp_result.close()

        if os.stat(path_output_blastp_xmlfile).st_size > IdProt.Strng.EST_MAX_SZE_BLSTP_RUN_200_KB.value:
            raise ValueError("blast xml output size is over 200 KB in size. Something may have gone wrong with the "
                             "blastp run.")

        return path_output_blastp_xmlfile

    # Write to file the result of blastp run (that has already been parsed & filtered from the raw blastp result).
    #
    # path_output   String          Absolute path to the output_data dir.
    # blastp_dict   Dictionary      Data structure parsed and filtered from Blastp run.
    # write_xml     Boolean         True to write xmlfile represetation of blastp_dict, True by default.
    # write_csv     Boolean         True to write csvfile represetation of blastp_dict, True by default.
    # write_json    Boolean         True to write jsonfile represetation of blastp_dict, False by default.
    @staticmethod
    def _write_idmaps_for_mysqldb(path_output, blastp_dict, write_xml=True, write_csv=True, write_json=False):
        idmap = {IdProt.Strng.SEQ_ID.value: blastp_dict[IdProt.Strng.QRY_SEQ_ID.value]}
        for alignment in blastp_dict[IdProt.Strng.IDENT_ALIGN_LST.value]:
            idmap[IdProt.Strng.ACC_NUM.value] = alignment[IdProt.Strng.ACC_NUM.value]
            idmap[IdProt.Strng.LEN.value] = alignment[IdProt.Strng.LEN.value]
            full_alt_names_flags = IdProt._extract_names_flags(alignment[IdProt.Strng.HIT_DEF.value],
                                                               IdProt.Strng.NAME_SRCH_STR.value,
                                                               IdProt.Strng.ALTNAME_SRCH_STR.value,
                                                               IdProt.Strng.FLGS_SRCH_STR.value)
            idmap[IdProt.Strng.NAME.value] = full_alt_names_flags[IdProt.Strng.NAME_SRCH_STR.value]
            idmap[IdProt.Strng.ALTNAME.value] = full_alt_names_flags[IdProt.Strng.ALTNAME_SRCH_STR.value]
            try:
                idmap[IdProt.Strng.FLGS.value] = full_alt_names_flags[IdProt.Strng.FLGS_SRCH_STR.value]
            except KeyError:
                print('This blastp hit has no flags.')
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
        path_output_blastp_fastaname = IdProt._build_dir_tree_with_intermed_dir(path_root=path_output,
                                                                            intermed_dir=Paths.DIR_BLASTP.value,
                                                                            fastadir=idmap[IdProt.Strng.SEQ_ID.value])
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
        path_output_blastp_fastaname = IdProt._build_dir_tree_with_intermed_dir(path_root=path_output,
                                                                                intermed_dir=Paths.DIR_BLASTP.value,
                                                                            fastadir=idmap[IdProt.Strng.SEQ_ID.value])
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
        path_output_blastp_fastaname = IdProt._build_dir_tree_with_intermed_dir(path_root=path_output,
                                                                            intermed_dir=Paths.DIR_BLASTP.value,
                                                                            fastadir=idmap[IdProt.Strng.SEQ_ID.value])
        path_output_blastp_fastaname_jsonfile = os.path.join(path_output_blastp_fastaname, 'idmap_swsprt.json')
        with open(path_output_blastp_fastaname_jsonfile, 'w') as f:
            f.write(json.dumps(idmap))

    # To make directory and all intermediate directories, encapsulated in try/except clause as pre-existing result in
    # exceptions being thrown.
    #
    # path_root     String      Abs path of root dir of tree to be built here (.../input_data/ or .../output/data/).
    # intermed_dir  String      Abs path of intermediate dir of tree, subdir to root, in which filename subdir will be.
    # fastadirname  String      Name of fastafile (without extension).
    @staticmethod
    def _build_dir_tree_with_intermed_dir(path_root, intermed_dir, fastadir):
        if fastadir is None:
            path_root_subdirs = os.path.join(path_root, intermed_dir)
        elif intermed_dir is Paths.DIR_BLASTP.value:
            fastadir = fastadir + IdProt.Strng.DIR_SUFIX.value
            path_root_subdirs = os.path.join(path_root, intermed_dir, fastadir)
        else:
            path_root_subdirs = os.path.join(path_root, intermed_dir, fastadir)
        try:
            os.makedirs(path_root_subdirs)
        except FileExistsError:
            print('Part or all of path already exists. This is absolutely fine.')
        return path_root_subdirs

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
                        IdProt.__is_reasonable_length(full_alt_names_flag[search_str])
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
