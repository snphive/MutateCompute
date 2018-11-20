#!/usr/bin/env python3

import os
import gzip
import shutil
import tarfile
from src.Str import Str

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class Parser(object):

    def make_tarfile(self, path_dir_txt_files_to_pack):
        """
        Recursively pack contents of specified directory to improve file transfer rates and reduce disk space usage.
        (The tar.gz was observed to reduce disk space usage about 5-fold.)
        :param path_dir_txt_files_to_pack: Absolute path of directory for source files to pack.
        """
        if not os.path.isdir(path_dir_txt_files_to_pack):
            raise FileNotFoundError('Directory to pack does not exist with this path.')

        with tarfile.open(path_dir_txt_files_to_pack + Str.TARGZEXT.value, 'w:gz') as tar:
            tar.add(path_dir_txt_files_to_pack, arcname=os.path.basename(path_dir_txt_files_to_pack))

    # def parse_agadir_outputs_to_1_file(self, path_output_agadir_3dotsdir):
    #     """
    #     Parse output of Agadir to retain wanted information and delete unwanted in order to reduce memory usage,
    #     and therefore improve data transfer speed.
    #     :param path_output_agadir_3dotsdir:
    #     :return:
    #     """
    #
    # def parse_biopython_outputs_to(self, path_):
    #     """
    #
    #     :param path_:
    #     :return:
    #     """
    #
    # def parse_foldx_outputs_to_ddG(self, path_):
    #     """
    #
    #     :param path_:
    #     :return:
    #     """
    #
    # def parse_compress_files(self):
    #     """
    #
    #     :return:
    #     """


