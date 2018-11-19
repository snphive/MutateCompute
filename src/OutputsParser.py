#!/usr/bin/env python3

import os
import gzip
import shutil
import tarfile

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class OutputsParser(object):

    def make_tarfile(self, path_dir_txt_files_to_pack):
        """

        :param path_dir_txt_files_to_pack: Absolute path of directory for source files
        :return:
        """
        if not os.path.isdir(path_dir_txt_files_to_pack):
            raise FileNotFoundError('Directory to pack does not exist with this path.')

        with tarfile.open(path_dir_txt_files_to_pack, 'w:gz') as tar:
            files_only = [f for f in os.listdir(path_dir_txt_files_to_pack) if
                          os.path.isfile(os.path.join(path_dir_txt_files_to_pack, f))]
            for file_to_pack in files_only:
                tar.add(file_to_pack, arcname=os.path.basename(file_to_pack))

    def compress_txt_files_to_gzip(self, path_txt_files_to_compress):
        """
        Compress text file(s) to gzip.
        :param path_txt_files_to_compress: Absolute path of text file or list of text files.
        :return:
        """
        if path_txt_files_to_compress.isinstance(str):
            path_txt_files_to_compress = [path_txt_files_to_compress]
        for path_txt_file_to_compress in path_txt_files_to_compress:
            with open(path_txt_file_to_compress, 'rb') as f_in:
                with gzip.open(path_txt_file_to_compress + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

    def compress_txt_files_to_tar(self, path_txt_files_to_compress):
        """
        Compress text file(s) to tar file.
        :param path_txt_files_to_compress: Absolute path of text file or list of text files.
        :return:
        """
        if path_txt_files_to_compress.isinstance(str):
            path_txt_files_to_compress = [path_txt_files_to_compress]
        for path_txt_file_to_compress in path_txt_files_to_compress:
            with open(path_txt_file_to_compress, 'rb') as f_in:
                with gzip.open(path_txt_file_to_compress + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)


    def parse_agadir_outputs_to_1_file(self, path_output_agadir_3dotsdir):
        """
        Parse output of Agadir to retain wanted information and delete unwanted in order to reduce memory usage,
        and therefore improve data transfer speed.
        :param path_output_agadir_3dotsdir:
        :return:
        """

    def parse_biopython_outputs_to(self, path_):
        """

        :param path_:
        :return:
        """

    def parse_foldx_outputs_to_ddG(self, path_):
        """

        :param path_:
        :return:
        """

    def parse_compress_files(self):
        """

        :return:
        """


