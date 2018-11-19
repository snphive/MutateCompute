import sys
import os
import glob
from src.Biopython import Biopy
from Bio.Blast import NCBIWWW
from src.Cluster import Cluster
from src.IdentifyProtein import IdProt
from src.GeneralUtilityMethods import GUM

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"

# path_input_fastafile_list = sys.argv[1].split(';')
path_input_fastafile_list = glob.glob(sys.argv[1] + '/*.fasta')
path_output = sys.argv[2]
path_config_job = sys.argv[3]
path_output_blastp = sys.argv[4]
write_idmaps_for_mysldb = sys.argv[5] == 'True'
write_csv = sys.argv[6] == 'True'
write_xml = sys.argv[7] == 'True'
write_json = sys.argv[8] == 'True'


for path_fastafile in path_input_fastafile_list:
    with open(path_fastafile) as fastafile_opened:
        fastafile_name = path_fastafile.split('/')[-1].split('.')[0]
        jobname = 'BLSTP_' + fastafile_name
        Cluster.write_job_q_bash(jobname=jobname, path_job_q_dir=path_config_job, queue='all.q', memory_limit_GB='3',
                                 cluster_node='hodor1.vib')
        path_output_blastp_fastaname = GUM._os_makedirs(path_output_blastp, fastafile_name)
        os.chdir(path_output_blastp_fastaname)
        Cluster.run_job_q(path_job_q_dir=path_config_job)
        Cluster.wait_for_grid_engine_job_to_complete(grid_engine_job_prefix=jobname)
        path_raw_blstp_xml = IdProt._write_raw_blast_xml(path_output, fastafile_name,
                                                         blastp_result=NCBIWWW.qblast(
                                                             program=Biopy.BlastParam.BLST_P.value,
                                                             database=Biopy.BlastParam.SWSPRT.value,
                                                             sequence=fastafile_opened.read(),
                                                             entrez_query=Biopy.BlastParam.HOMSAP_ORG.value,
                                                             alignments=Biopy.BlastParam.MAX_ALIGN_20.value,
                                                             hitlist_size=Biopy.BlastParam.MAX_HIT_20.value))
        blastp_dict = Biopy.parse_filter_blastp_xml_to_dict(path_raw_blstp_xml, fastafile_name, path_fastafile)
        # blastp_dict_list.append(blastp_dict)
        if write_idmaps_for_mysldb:
            IdProt._write_idmaps_for_mysqldb(path_output, blastp_dict, write_csv=write_csv, write_xml=write_xml,
                                             write_json=write_json)
