import os
import json
import glob
import natsort
import threading
import time
from src.Str import Str
from src.Biopython import Biopy
from src.Paths import Paths
from src.Cluster import Cluster
from src.GeneralUtilityMethods import GUM

# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)


# Class to take any protein identifier, pdb file, FASTA, and identify the protein, potentially via Blast and/or then
# mapping whatever identifier the input has with the recognised identifiers (1_A to 32431_A from RvdK's repaired pdbs)
class IdProt(object):

    @staticmethod
    def map_seq_to_swsprt_acc_id_and_write_files(path_input_fastafiles: list, path_output: str, write_idmaps_for_mysqldb: bool,
                                                 write_csv=True, write_xml=True, write_json=False):
        """
        Expects a directory location of fastafiles (not a fastafile itself).
        (This method relies on the presence of fastafiles in the specified dir, in order for them to be run on Blastp.
        This transfer is currently done manually)
        :param path_input_fastafiles: Absolute path of root directory for input fastafiles (e.g. /input_data/fastas_10).
        :param path_output: Absolute path of root directory for blastp output files (..../output_data/).
        :param write_idmaps_for_mysqldb: True (by default) builds dictionary mapping RvdK's ids to swsprt accession nos & write files.
        :param write_csv: True to write csvfiles.
        :param write_xml: True to write xmlfiles.
        :param write_json: True to write jsonfiles.
        :return: List of dictionary data structure representations of each parsed & filtered Blastp run result.
        """
        if isinstance(path_input_fastafiles, str):
            path_input_fastafiles = [path_input_fastafiles]
        blastp_dict_list = []
        # There are problems with using Biopython.Blast on the cluster that I have not yet solved. I may use the
        # blast module that is loaded on the cluster (v 2.5.0+) instead of via Biopython.
        if GUM.using_cluster():
            # for path_fastafile in path_input_fastafile_list:
            #     with open(path_fastafile) as fastafile_opened:
            #         fastafile_name = path_fastafile.split('/')[-1].split('.')[0]
            #         jobname = 'BLSTP_' + fastafile_name
            #         Cluster.write_job_q_bash(job_name=jobname, path_job_q_dir=Paths.CONFIG_JOBQ)
            #         path_output_blastp_fastaname = GUM._os_makedirs(Paths.OUTPUT_BLASTP, fastafile_name)
            #         os.chdir(path_output_blastp_fastaname)
            #         Cluster.run_job_q(path_job_q_dir=Paths.CONFIG_JOBQ)
            #         Cluster.wait_for_grid_engine_job_to_complete(grid_engine_jobname=jobname)
            #         path_blstp_xml = IdProt._write_raw_blast_xml(path_output, fastafile_name,
            #                                                 blastp_result=Biopy.run_blastp(fastafile_opened.read()))
            #         blastp_dict = Biopy.parse_filter_blastp_xml_to_dict(path_blstp_xml, fastafile_name, path_fastafile)
            #         # blastp_dict_list.append(blastp_dict)
            #         if write_idmaps_for_mysqldb:
            #             IdProt._write_idmaps_for_mysqldb(path_output, blastp_dict, write_csv=write_csv,
            #                                              write_xml=write_xml,
            #                                              write_json=write_json)

            python_script_w_paths = os.path.join(Paths.SRC, 'run_BlstpZeus.py') + ' ' + path_input_fastafiles + ' ' \
                                    + path_output + ' ' + Paths.CONFIG_BLST_JOBQ + ' ' + Paths.OUTPUT_BLASTP + ' ' + \
                                    str(write_idmaps_for_mysqldb) + ' ' + str(write_csv) + ' ' + str(write_xml) + \
                                    ' ' + str(write_json)
            Cluster.write_job_q_bash(job_name='IdProtJobs', path_job_q_dir=Paths.CONFIG_BLST_JOBQ,
                                     python_script_with_paths=python_script_w_paths)
            Cluster.run_job_q(path_job_q_dir=Paths.CONFIG_BLST_JOBQ)
        else:
            for path_fastafile in path_input_fastafiles:
                with open(path_fastafile) as f:
                    fasta_str = f.read()
                    fastafile_name = path_fastafile.split('/')[-1].split('.')[0]
                if IdProt._has_all_A_sequence(path_fastafile):
                    print('This sequence has all As, BLAST would think it is a nucleotide sequence and fail. So it is '
                          'not being run: ' + path_fastafile)
                    continue
                blastp_result = Biopy.run_blastp(fasta_str)
                path_blstp_xml = IdProt._write_raw_blast_xml(path_output, fastafile_name, blastp_result)
                blastp_dict = Biopy.parse_filter_blastp_xml_to_dict(path_blstp_xml, fastafile_name, path_fastafile)
                blastp_dict_list.append(blastp_dict)
                if write_idmaps_for_mysqldb:
                    IdProt._write_idmaps_for_mysqldb(path_output, blastp_dict, write_csv=write_csv, write_xml=write_xml,
                                                     write_json=write_json)
        return blastp_dict_list

    @staticmethod
    def _has_all_A_sequence(path_fastafile):
        is_all_A = True
        fasta_str = GUM.get_sequenceOnly_from_fastafile(path_fastafile)
        for aa in fasta_str:
            if aa != 'A':
                is_all_A = False
                break
        return is_all_A

    @staticmethod
    def _write_raw_blast_xml(path_output: str, fastafile_name: str, blastp_result):
        """
        Writes the blastp result file to an xml file in a output_data subdir called blastp.
        :param path_output: Absolute path of output_data root directory, typically ../output_data.
        :param fastafile_name: Used for naming the xml file.
        :param blastp_result: Biopython's raw blastp output. (Buffered text stream (io.StringIO).
        :return: Path to the new raw blastp xml file written by this method.
        """
        if blastp_result is None:
            raise ValueError('Blastp result is none. Something went wrong!')
        path_output_blastp = IdProt._build_dir_tree_with_intermed_dir(path_root=path_output,
                                                                      intermed_dir=Paths.DIR_BLASTP.value,
                                                                      fastadir=fastafile_name)
        path_output_blastp_xmlfile = os.path.join(path_output_blastp, fastafile_name + ".xml")
        with open(path_output_blastp_xmlfile, 'w') as f:
            blstp_read = blastp_result.read()
            f.write(blstp_read)
            blastp_result.close()
        if os.stat(path_output_blastp_xmlfile).st_size > IdProt.Strng.EST_MAX_SZE_BLSTP_RUN_300_KB.value:
            # Not sure this needs to raise ValueError, rather than exception so that it can continue rather than halt execution.
            raise ValueError("blast xml output size is over 300 KB in size. Something may have gone wrong with the "
                             "blastp run.")
        return path_output_blastp_xmlfile

    @staticmethod
    def _write_idmaps_for_mysqldb(path_output: str, blastp_dict: dict, write_xml=True, write_csv=True, write_json=False):
        """
        Write to file the result of blastp run (that has already been parsed & filtered from the raw blastp result).
        :param path_output: Absolute path of output_data root directory.
        :param blastp_dict: Parsed and selected attributes of Blastp result.
        :param write_xml: True (by default) writes xml file representation of blastp_dict.
        :param write_csv: True (by default) writes csv file representation of blastp_dict.
        :param write_json: True (by default) writes json file representation of blastp_dict.
        """
        idmap = {IdProt.Strng.SEQ_ID.value: blastp_dict[IdProt.Strng.QRY_SEQ_ID.value]}
        # time.sleep(1)
        for alignment in blastp_dict[IdProt.Strng.IDENT_ALIGN_LST.value]:
            # time.sleep(1)
            idmap[IdProt.Strng.ACC_NUM.value] = alignment[IdProt.Strng.ACC_NUM.value]
            idmap[IdProt.Strng.LEN.value] = alignment[IdProt.Strng.LEN.value]
            full_alt_names_flags = IdProt._extract_names_flags(alignment[IdProt.Strng.HIT_DEF.value],
                                                               IdProt.Strng.NAME_SRCH_STR.value,
                                                               IdProt.Strng.ALTNAME_SRCH_STR.value,
                                                               IdProt.Strng.FLGS_SRCH_STR.value)
            try:
                idmap[IdProt.Strng.NAME.value] = full_alt_names_flags[IdProt.Strng.NAME_SRCH_STR.value]
            except KeyError:
                print('This blastp hit has no name (??). So none will be added to the result idmap.')
            try:
                idmap[IdProt.Strng.ALTNAME.value] = full_alt_names_flags[IdProt.Strng.ALTNAME_SRCH_STR.value]
            except KeyError:
                print('This blastp hit has no alt name. So none will be added to the result idmap.')
            try:
                idmap[IdProt.Strng.FLGS.value] = full_alt_names_flags[IdProt.Strng.FLGS_SRCH_STR.value]
            except KeyError:
                print('This blastp hit has no flags, so no flag written to result idmap.')
            idmap[IdProt.Strng.STRT_POS.value] = alignment[IdProt.Strng.HSP_DICT.value][IdProt.Strng.SBJCT_STRT.value]
            idmap[IdProt.Strng.END_POS.value] = alignment[IdProt.Strng.HSP_DICT.value][IdProt.Strng.SBJCT_END.value]
        if write_csv:
            IdProt._write_idmap_csv(path_output, idmap)
        if write_xml:
            IdProt._write_idmap_xml(path_output, idmap)
        if write_json:
            IdProt._write_idmap_jsonfile(path_output, idmap)

    @staticmethod
    def _write_idmap_xml(path_output: str, idmap: dict):
        """
        Converts data passed in a dictionary data structure to csv format and writes it to the blastp directory.
        :param path_output: Absolute path of output_data dir, destination for xml file.
        :param idmap: Data structure of values taken from blastp run that are parsed to a csv format.
        :return: Xml elements (with one attribute value for id) and newlines in their own element in a list.
        """
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

    @staticmethod
    def _write_idmap_csv(path_output: str, idmap: dict):
        """
        Converts data passed in a dictionary data structure to csv format and writes it to the blastp directory.
        :param path_output: Abs path of output_data dir, destination for csv file.
        :param idmap: Data structure of values taken from blastp run that are parsed to a csv format.
        :return: The values and commas and newlines in their own element in a list.
        """
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

    @staticmethod
    def _write_idmap_jsonfile(path_output: str, idmap: dict):
        """
        Writes out a json file of the idmap (dictionary) using an imported library method from 'json' module.
        :param path_output: Abs path of the output dir.
        :param idmap:
        """
        path_output_blastp_fastaname = IdProt._build_dir_tree_with_intermed_dir(path_root=path_output,
                                                                            intermed_dir=Paths.DIR_BLASTP.value,
                                                                            fastadir=idmap[IdProt.Strng.SEQ_ID.value])
        path_output_blastp_fastaname_jsonfile = os.path.join(path_output_blastp_fastaname, 'idmap_swsprt.json')
        with open(path_output_blastp_fastaname_jsonfile, 'w') as f:
            f.write(json.dumps(idmap))

    @staticmethod
    def _build_dir_tree_with_intermed_dir(path_root: str, intermed_dir: str, fastadir: str):
        """
        To make directory and all intermediate directories, encapsulated in try/except clause as pre-existing result in
        exceptions being thrown.
        :param path_root: Abs path of the root of this dir tree.
        :param intermed_dir: Name of intermediate directory.
        :param fastadir: Name of Fasta directory.
        :return: Path of new directory and subdirectories.
        """
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
            print(Str.PARTALLPATHEXISTS_MSG.value)
        return path_root_subdirs

    @staticmethod
    def _extract_names_flags(hit_def: str, *search_strngs: str):
        """
        NOTE: This method assumes that the 'Recname', 'Altname', 'Flags' are always written in this order. It also
        assumes that the Flag is the last string of any form in the line.
        I don't know yet if these assumptions are correct.
        :param hit_def: A value from Blastp result, typically includes semi-colon separated text regarding names & flags.
        :param search_strngs: Search strings, typically name, altname, flags.
        :return:
        """
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

    @staticmethod
    def __is_reasonable_length(name_or_flag: str):
        """
        :param name_or_flag: Text to be tested for length.
        :return: True if text passed is between 2 and 30 characters long. I cannot say for sure that there are no names
        or flags that are 1 character long or more than 30 characters long - this is running assumption.
        """
        reasonable_max_name_or_flag_len = 70
        if len(name_or_flag) < 2 | len(name_or_flag) > reasonable_max_name_or_flag_len:
            raise ValueError('Full name, altname or flag is either not there at all or seems too long.')
        return True

    from enum import Enum

    class Strng(Enum):
        EST_MAX_SZE_BLSTP_RUN_300_KB = 300000
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

# pydevd.stoptrace()
