#!/usr/bin/env python3
"""
Script for running any part of the codebase on relatively small numbers of files typically housed in a single folder.
Those to be run on multiple folders can be done through the KickOff.py module.

PackCompressManual.py can be run locally or on cluster.
"""
import sys
import os
from src.enums.Paths import Paths
from src.tools.OutputsParser import Parser

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

"""
Paths for the entire codebase are set accordingly. "use_cluster" is set to False by default.  
"""
Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))

"""
Select directory for packing into tar.
"""
path_dir_txt_files_to_pack = os.path.join(Paths.OUTPUT_AGADIR, '02064053-3f32-32e6-9660-aaaffc30db87')

"""
Pack files in directory into tar.
"""
Parser().make_tarfile(path_dir_txt_files_to_pack)
