#!/usr/bin/env python3

import os
import tarfile
from src.enums.Str import Str

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


class Parser(object):

    def make_tarfile(self, path_files_to_pack_dir: str):
        """
        Recursively pack contents of specified directory to improve file transfer rates and reduce disk space usage.
        (The tar.gz was observed to reduce disk space usage about 5-fold.)
        :param path_files_to_pack_dir: Absolute path of directory for source files to pack.
        """
        if not os.path.isdir(path_files_to_pack_dir):
            raise FileNotFoundError('Directory to pack does not exist with this path.')

        path_files_to_pack_tar_gz = path_files_to_pack_dir + Str.TARGZEXT.value

        if os.path.exists(path_files_to_pack_tar_gz):
            print(path_files_to_pack_tar_gz + ' already exists. This will be overwritten by new one.')
        with tarfile.open(path_files_to_pack_tar_gz, 'w:gz') as tar:
            tar.add(path_files_to_pack_dir, arcname=os.path.basename(path_files_to_pack_dir))


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


