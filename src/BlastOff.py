from src.IdentifyProtein import IdProt
from src.Paths import Paths

IdProt.map_seq_to_swsprt_acc_id_and_write_files(path_input=Paths.INPUT, path_output=Paths.OUTPUT,
                                                write_idmaps_for_mysldb=True, write_csv=True, write_xml=True,
                                                write_json=False)
# As with the MutateFasta, Agadir and FoldX programs, the start of operations begins with identifying which or how many
# fasta files are to be analysed. These must then be copied over from the repository directory (REPO_PDB_FASTA/FASTAs_),
# to the input_data directory.
#
# (And as with the other programs, the necessary destination directory tree for the input files being copied over will
# be created if not already present. This involves a subdirectory bearing the name of the file, which for fastafiles is
# typically '1_A' (unchanged from RvdK's naming scheme).)
#
# E.g. '1_A.fasta' copied from ~/REPO_PDB_FASTA/fastas_100/ to ~/PycharmProjects/MutateCompute/input_data/1_A/
#
#
# PROGRAM           INPUT FILES & DIRECTORIES                           OUTPUT FILES & DIRECTORIES
# -------------------------------------------------------------------------------------------------------------------
# IdentifyProtein   input_data/<fastafilename>/                         output_data/blastp/
#                                                                       output_data/blastp/<fastafilename>_idmaps/
