#!/usr/bin/env python3
"""
Script for launching the packing and compressing files into tar files in a specified directory.

PackCompress.py can be run locally or on cluster (the latter via a bash script).
"""
import sys
import os
import glob
from src.enums.Paths import Paths
from src.tools.GeneralUtilityMethods import GUM

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

"""
1. Set up paths. ("use_cluster" is set to False by default.)  
"""
Paths.set_up_paths(use_cluster=(len(sys.argv) > 1 and sys.argv[1].strip(' ') == 'use_cluster=True'))

"""
2. Select directory of files for packing & compressing: Agadir-related
"""
# path_dir_txt_files_to_pack = os.path.join(Paths.OUTPUT_AGADIR, '02064053-3f32-32e6-9660-aaaffc30db87')
path_output_agad_3dots_dir = os.path.join(Paths.OUTPUT_AGADIR, '1...250')
path_files_to_pack_dirs = glob.glob(path_output_agad_3dots_dir + '/*')

"""
3. Select pdbname(s), outputs of which are to be packed & compressed. 
"""
pdbname = ['Repair_14']

"""
4. Select directory of files for packing & compressing: FoldX-related
"""
path_files_to_pack_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname))
path_files_to_pack_dirs.append(os.path.join(Paths.OUTPUT_AC, pdbname))

"""
5. Pack files in directory into tar.
"""
for path_files_to_pack_dir in path_files_to_pack_dirs:
    GUM.make_tarfile(path_files_to_pack_dir)
