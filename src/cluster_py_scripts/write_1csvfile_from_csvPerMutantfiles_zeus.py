#!/usr/bin/env python3
"""
Script for writing one csv file, reading output data from csv files that were written 1 per mutant, from within cluster,
specified in job.q.
"""
import sys
import os
from src.enums.Paths import Paths
from src.enums.Str import Str

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

path_csvfile = sys.argv[1]
path_output_ac_or_bm_dir = sys.argv[2]
pdbname = sys.argv[3]
ddG = ''
buildmodel_csvfile_header = 'pdb,fxmutant,ddG'
analysecomplex_csvfile_header = 'pdb,fxmutant,interaction ddG'
csvfile_header = ''

if Paths.DIR_BM.value in path_output_ac_or_bm_dir:
    ddG = '_ddG'
    csvfile_header = buildmodel_csvfile_header
elif Paths.DIR_AC.value in path_output_ac_or_bm_dir:
    ddG = '_interaction_ddG'
    csvfile_header = analysecomplex_csvfile_header
else:
    print('The output directory path does not contain build_model or analyse_complex. If you want to read/write Agadir-related '
          'csv output files, you need to complete the logic for it here in the run_write_1csvfile_.. script')


path_output_ac_or_bm_ddG_csv_dumpfile = os.path.join(path_output_ac_or_bm_dir, pdbname + ddG + Str.CSVEXT.value)


def _has_fx_bm_header(header_line: str):
    return csvfile_header in header_line


def _has_fx_ac_header(header_line: str):
    return csvfile_header in header_line


if not os.path.exists(path_output_ac_or_bm_ddG_csv_dumpfile):
    with open(path_output_ac_or_bm_ddG_csv_dumpfile, 'w') as ddG_csv:
        ddG_csv.write(csvfile_header + Str.NEWLN.value)

else:
    with open(path_output_ac_or_bm_ddG_csv_dumpfile, 'w') as ddG_csv:
        if len(ddG_csv.readlines()) == 0:
            print('Single csv output file found but file is empty (not even a header). Header is being written now.')
            ddG_csv.write(csvfile_header + Str.NEWLN.value)
        elif _has_fx_bm_header(ddG_csv.readlines()[0]):
            print('Detected header for build_model csv dumpfile in this csvfile. Is proceeding to append csv bm output data.')
        elif _has_fx_ac_header(ddG_csv.readlines()[0]):
            print('Detected header for analyse_complex csv dumpfile in this csvfile. Is proceeding to append csv ac output data.')
        else:
            print('Not detected headers for either build_model or analyse_complex. There is a bug.')


def _is_at_header_line(this_line: str):
    return csvfile_header == this_line


def _is_at_data_line(this_line: str):
    return 'PDB,' in this_line


with open(path_csvfile, 'r') as f:
    for line in f.readlines():
        if _is_at_header_line(line):
            continue
        if _is_at_data_line(line):
                with open(path_output_ac_or_bm_ddG_csv_dumpfile, 'a+') as csv_f:
                    csv_f.write(line + Str.NEWLN.value)



"""
path_output_bm_pdb_avg_csvfiles = []
written_bm_to_csv = False
if write_bm_to_csv:
    fx = FoldX()
    path_output_bm_pdb_fxmutant_dirs = []
    for path_pdbfile in path_pdbfiles:
        pdbname = os.path.basename(path_pdbfile).split('.')[0]
        for specific_fxmutant in specific_fxmutants:
            path_output_bm_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_BM, pdbname, specific_fxmutant))
        if not specific_fxmutants:
            path_output_bm_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_BM, pdbname, '*'))
        for path_output_bm_pdb_fxmutant_dir in path_output_bm_pdb_fxmutant_dirs:
            bm = fx.BuildModel(Cond.INCELL_MAML_FX.value)
            path_output_bm_pdb_avg_csvfile = bm.write_bm_avg_fxout_to_1csvfile_up_2dirlevels(path_output_bm_pdb_fxmutant_dir)
            path_output_bm_pdb_avg_csvfiles.append(path_output_bm_pdb_avg_csvfile)
    written_bm_to_csv = True

# NEED TO KNOW WHAT VALUES ARE SIGNIFICANT, E.G. ANYTHING < 0.01 KCAL/MOL SIGINIFICANT ??
path_output_ac_pdb_sumry_csvfiles = []
if write_ac_to_csv:
    fx = FoldX()
    path_output_ac_pdb_fxmutant_dirs = []
    for path_pdbfile in path_pdbfiles:
        pdbname = os.path.basename(path_pdbfile).split('.')[0]
        for specific_fxmutant in specific_fxmutants:
            path_output_ac_pdb_fxmutant_dirs.append(os.path.join(Paths.OUTPUT_AC, pdbname, specific_fxmutant))
        if not specific_fxmutants:
            path_output_ac_pdb_fxmutant_dirs = glob.glob(os.path.join(Paths.OUTPUT_AC, pdbname, '*'))
        for path_output_ac_pdb_fxmutant_dir in path_output_ac_pdb_fxmutant_dirs:
            ac = fx.AnalyseComplex(Cond.INCELL_MAML_FX.value)
            path_output_ac_pdb_sumry_csvfile = ac.write_ac_sumry_fxout_to_1csvfile_up_2dirlevels(path_output_ac_pdb_fxmutant_dir)
            path_output_ac_pdb_sumry_csvfiles.append(path_output_ac_pdb_sumry_csvfile)
            GUM.linux_remove_dir(path_output_ac_pdb_fxmutant_dir)
            if written_bm_to_csv:
                path_output_bm_pdb_fxmutant_dir = path_output_ac_pdb_fxmutant_dir.split('/')
                path_output_bm_pdb_fxmutant_dir[-3] = Paths.DIR_BM.value
                path_output_bm_pdb_fxmutant_dir = '/'.join(path_output_bm_pdb_fxmutant_dir)
                GUM.linux_remove_dir(path_output_bm_pdb_fxmutant_dir)
"""
