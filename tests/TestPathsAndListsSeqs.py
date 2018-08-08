from enum import Enum


# All paths are constants and stored in the configuration directory in 'pathsAndDictionaries' yaml file.
# They are read into the Paths enum from which all other classes source all path strings.
# However for tests, the absolute paths are directed towards the "tests" directory, so tests should never need to
# access the yaml configuration file.
#
# All paths are absolute paths unless prefix with an "REL" which indicates they are relative.
# All paths are prefixed with "PATH"
class TPLS(Enum):

    # ABSOLUTE PATHS FOR TESTS ROOT
    MC_TESTS = "/Users/u0120577/PycharmProjects/MutateCompute/tests"

    # CONFIGURATION FILES FROM MAIN DIRECTORY (i.e. real data)
    CONFIG_FOR_READ_ONLY = "/Users/u0120577/PycharmProjects/MutateCompute/configuration"
    # INPUT FILES FROM MAIN DIRECTORY (i.e. real data)
    INPUT_FOR_READ_ONLY = "/Users/u0120577/PycharmProjects/MutateCompute/input_data"

    # REPOSITORY OF PDBs and FASTAs
    REPO_PDB_FASTA = '/Users/u0120577/REPO_PDB_FASTA'

    # These relative paths might not be used - entire config dir gets copied as a whole into /tests/ at start of tests.
    REL_CONFIG = "/configuration"
    REL_CONFIG_AGAD = "/configuration/agadir_config"
    REL_CONFIG_JOBQ = "/configuration/cluster_jobq"
    REL_CONFIG_FX = "/configuration/foldx_config"
    REL_CONFIG_GLOBAL_OPTIONS = "/configuration/global_options"
    REL_CONFIG_FX_ACRUNSCRIPT = "/configuration/foldx_config/ac_runscript"
    REL_CONFIG_FX_BMRUNSCRIPT = "/configuration/foldx_config/bm_runscript"

    # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
    MC_TESTS_CONFIG = MC_TESTS + REL_CONFIG
    MC_TESTS_CONFIG_AGAD = MC_TESTS + REL_CONFIG_AGAD
    MC_TESTS_CONFIG_JOBQ = MC_TESTS + REL_CONFIG_JOBQ
    MC_TESTS_CONFIG_FX = MC_TESTS + REL_CONFIG_FX
    MC_TESTS_CONFIG_GLOBAL_OPTIONS = MC_TESTS + REL_CONFIG_GLOBAL_OPTIONS
    MC_TESTS_CONFIG_FX_ACRUNSCRIPT = MC_TESTS + REL_CONFIG_FX_ACRUNSCRIPT
    MC_TESTS_CONFIG_FX_BMRUNSCRIPT = MC_TESTS + REL_CONFIG_FX_BMRUNSCRIPT

    # REFERENCE FILES - NOTE: THESE ARE ONLY USED FOR TESTS
    REL_REFFILES = "/reference_files"
    # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
    MC_TESTS_REFFILES = MC_TESTS + REL_REFFILES

    # INPUT PATHS ONLY
    # e.g. /tests/input_data/<pdbname>/all_mutants/<fxmutantchainname>
    REL_INPUT = "/input_data"
    REL_ALL_MUTANTS = "/all_mutants"
    # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
    MC_TESTS_INPUT = MC_TESTS + REL_INPUT

    # OUTPUT PATHS ONLY
    # e.g.#1: tests/output_data/<pdbname>/foldx/build_model/<fxmutantchainname>
    # e.g.#2: tests/output_data/<pdbname>/mutate_fasta
    REL_OUTPUT = "/output_data"
    REL_FX_AC = "/foldx/analyse_complex"
    REL_FX_BM = "/foldx/build_model"
    REL_MUTATE_FASTA = "/mutate_fasta"
    # ABSOLUTE PATH BUILT FROM ROOT AND RELATIVE PATHS
    MC_TESTS_OUTPUT = MC_TESTS + REL_OUTPUT

    DICT_AA_1TO3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
                    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}
    DICT_AA_3TO1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'H1S': 'H',
                    'H2S': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
                    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
    LIST_ALL_20_AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                      'Y']

    ZEUS_R_EXE = "/software/shared/apps/general/R/3.1.2/bin/Rscript"
    ZEUS_FOLDX_EXE = "/switchlab/group/tools/FoldX_2015/FoldX"
    ZEUS_AGADIR_EXE = "/switchlab/group/tools/agadir_10042012/agadirwrapper"
    ZEUS_QSUB_EXE = "/opt/sge/bin/lx-amd64/"

    LOCAL_R_EXE = "/usr/local/bin/R"
    LOCAL_FOLDX_EXE = "/Users/u0120577/SNPEFFECT/executables/FoldX"
    LOCAL_AGADIR_EXE = "/switchlab/group/tools/agadir_10042012/agadirwrapper"

    # FASTA SEQUENCES
    FASTA_SEQ_1_A = 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGNFVRVIQTFN' \
                    'RTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRSLTKRNAVRTDQHNSKWLS' \
                    'EPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARLVCSVTDEDGPETHFDELEDVFLLETDN' \
                    'PRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRPGTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSI' \
                    'YPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVLPTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVS' \
                    'QVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRRSRRQDVRHGNPLTQCR'
    FASTA_SEQ_1_B = 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGNFVRVIQTFN' \
                    'RTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRSLTKRNAVRTDQHNSKWLS' \
                    'EPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARLVCSVTDEDGPETHFDELEDVFLLETDN' \
                    'PRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRPGTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSI' \
                    'YPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVLPTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVS' \
                    'QVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRRSRRQDVRHGNPLTQCR'
    FASTA_SEQ_2_A = 'EIVQYGVKNNTTFLECAPKSPQASIKWLLQKDKDRRKEVKLNERIIATSQGLLIRSVQGSDQGLYHCIATENSFKQTIAKINFKVLD'
    FASTA_SEQ_3_A = 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGNFVRVIQTFN' \
                    'RTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRSLTKRNAVRTDQHNSKWLS' \
                    'EPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARLVCSVTDEDGPETHFDELEDVFLLETDN' \
                    'PRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRPGTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSI' \
                    'YPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVLPTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVS' \
                    'QVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRRSRRQDVRHGNPLTQCR'
    FASTA_SEQ_3_B = 'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVGSKDHILSLNINNISQEALSVFWPASTIKVEECKMAGKDPTHGCGNFVRVIQTFN' \
                    'RTHLYVCGSGAFSPVCTYLNRGRRSEDQVFMIDSKCESGKGRCSFNPNVNTVSVMINEELFSGMYIDFMGTDAAIFRSLTKRNAVRTDQHNSKWLS' \
                    'EPMFVDAHVIPDGTDPNDAKVYFFFKEKLTDNNRSTKQIHSMIARICPNDTGGLRSLVNKWTTFLKARLVCSVTDEDGPETHFDELEDVFLLETDN' \
                    'PRTTLVYGIFTTSSSVFKGSAVCVYHLSDIQTVFNGPFAHKEGPNHQLISYQGRIPYPRPGTCPGGAFTPNMRTTKEFPDDVVTFIRNHPLMYNSI' \
                    'YPIHKRPLIVRIGTDYKYTKIAVDRVNAADGRYHVLFLGTDRGTVQKVVVLPTNNSVSGELILEELEVFKNHAPITTMKISSKKQQLYVSSNEGVS' \
                    'QVSLHRCHIYGTACADCCLARDPYCAWDGHSCSRFYPTGKRRSRRQDVRHGNPLTQCR'
    FASTA_SEQ_10_B = 'SCIQFTRHASDVLLNLNRLRSRDILTDVVIVVSREQFRAHKTVLMACSGLFYSIFTDQLKCNLSVINLDPEINPEGFCILLDFMYTSRLNLREGN' \
                     'IMAVMATAMYLQMEHVVDTCRKFIKAS'
    MINI_FASTA_SEQ_WILDTYPE = 'SCI'
    MINI_FASTA_TITLE_SEQ_DICT_ALL_20 = {'S1A': 'ACI', 'S1C': 'CCI', 'S1D': 'DCI', 'S1E': 'ECI', 'S1F': 'FCI',
                                        'S1G': 'GCI', 'S1H': 'HCI', 'S1I': 'ICI', 'S1K': 'KCI', 'S1L': 'LCI',
                                        'S1M': 'MCI', 'S1N': 'NCI', 'S1P': 'PCI', 'S1Q': 'QCI', 'S1R': 'RCI',
                                        'S1S': 'SCI', 'S1T': 'TCI', 'S1V': 'VCI', 'S1W': 'WCI', 'S1Y': 'YCI',
                                        'C2A': 'SAI', 'C2C': 'SCI', 'C2D': 'SDI', 'C2E': 'SEI', 'C2F': 'SFI',
                                        'C2G': 'SGI', 'C2H': 'SHI', 'C2K': 'SKI', 'C2L': 'SLI', 'C2M': 'SMI',
                                        'C2N': 'SNI', 'C2P': 'SPI', 'C2Q': 'SQI', 'C2R': 'SRI', 'C2S': 'SSI',
                                        'C2T': 'STI', 'C2V': 'SVI', 'C2W': 'SWI', 'C2Y': 'SYI', 'I3A': 'SCA',
                                        'I3C': 'SCC', 'I3D': 'SCD', 'I3E': 'SCE', 'I3F': 'SCF', 'I3G': 'SCG',
                                        'I3H': 'SCH', 'I3I': 'SCI', 'I3K': 'SCK', 'I3L': 'SCL', 'I3M': 'SCM',
                                        'I3N': 'SCN', 'I3P': 'SCP', 'I3Q': 'SCQ', 'I3R': 'SCR', 'I3S': 'SCS',
                                        'I3T': 'SCT', 'I3V': 'SCV', 'I3W': 'SCW', 'I3Y': 'SCY'}
