import os
from src.GeneralUtilityMethods import GUM
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC


class MutateFasta(object):

    # def __init__(self, path_input, fastafile_list, mutant_aa_list):
    #     mutants = self.mutate_every_residue_in_fasta_list(path_input, fastafile_list, mutant_aa_list)
    #     self._write_mutants_to_file(fastamutants, path_input, make_fastafile_for_all_mutants,
    #                            make_fastafile_per_mutant,  path_output, make_csv_file=True, make_txt_file=True)

    def __init__(self):
        print('you be calling MutateFasta constructor, so you are.')

    # Mutates FASTA (typically wild-type) sequences at every position, to every residue specified (typically all other
    # 19 residues.)
    #
    # path_input            String      Abs path to input_data directory to which mutant fastafiles will be written.
    # fastafile_list        List        Where concurrency is not being used, the Schedular will pass the whole list of
    #                                   fastafiles here.
    # fastafile             String      Where concurrency is being used, the Scheduler will pass one Fastafile per
    #                                   thread/process.
    # mutant_aa_list        List        All residues that the sequence will be mutated to.
    # write_1_fasta_only    Boolean     Write 1 fastafile for all mutants.
    # write_fasta_per_mut   Boolean     Write 1 fastafile for every mutant.
    # path_output           String      Abs path of output_data dir (to which csv and/or txt files would be written)
    # write_csv             Boolean     True to write 1 csv file contain all mutants/. False by default.
    # write_txt             Boolean     True to write 1 txt file contain all mutants. False by default.
    #
    # Returns a dictionary of the wt-title as key, dict as value which has evert mutant title as key, sequence as value.
    def mutate_every_residue(self, path_input, fastafile_list, fastafile, mutant_aa_list, write_1_fasta_only,
                             write_fasta_per_mut, path_output, write_csv=False, write_txt=False):
        path_fastafile = self.build_complete_paths_for_fastafiles(path_input, fastafile)
        titleSeq = self.make_titleSeqDict_from_fastafile(path_fastafile)
        title_titleSeq = self.convert_titleSeqDict_to_titleTitleSeqDictDict(titleSeq)
        title_titleSeq_w_mutants = self._populate_title_titleSeq_with_mutants(title_titleSeq, mutant_aa_list)
        self._write_mutants(title_titleSeq_w_mutants, path_input, write_1_fasta_only, write_fasta_per_mut,
                            path_output, write_csv, write_txt)
        return title_titleSeq_w_mutants

    # title_titleSeq       Dictionary   Fasta title is key. Value is mutant-title: mutated-sequence dictionary.
    #                                   This already has the dict of wt-title: (wt-title:wt-sequence) dict.
    # mutant_aa_list       List         All residues that the sequence(s) should be mutated to.
    #
    # Returns wt-title: mutant-title:mutant-sequence
    def _populate_title_titleSeq_with_mutants(self, title_titleSeq, mutant_aa_list):
        title_titleSeq_w_mutants = title_titleSeq
        for wt_title in list(title_titleSeq.keys()):
            wt_seq = title_titleSeq[wt_title][wt_title]
            title_titleSeq_w_mutants[wt_title] = self._add_mutantTitle_mutatedSeq_to_dict(
                title_titleSeq[wt_title], wt_title, wt_seq, mutant_aa_list)
        return title_titleSeq_w_mutants

    # A sequence is mutated at every position to every residue that is included in the list (mutant_aa_list).
    #
    # E.g.
    # input as {'230498_A': 'MDVFM'}
    # output as {'230498_A': 'MDVFM','M1A': 'ADVFM','M1C': 'CDVFM', etc.. }
    #
    # titleSeq          Dictionary      Already should have wt-title: wt-title:wt-sequence dictionary.
    # wt_title          String          Name of fasta sequence (same as the fastafile name)
    # wt_seq            String          Amino acid sequence in the fastafile (in FASTA format).
    # mutant_aa_list    List            Specified list of residues, to mutate every residue in protein to.
    #
    # Returns a dictionary of sequence name as key, FASTA sequence itself as value, starting with wild-type and
    # followed by every mutant according to a specified list of residues.
    def _add_mutantTitle_mutatedSeq_to_dict(self, titleSeq, wt_title, wt_seq, mutant_aa_list):
        titleSeqDict_w_mutants = {wt_title: titleSeq[wt_title]}
        mutable_seq = MutableSeq(wt_seq, IUPAC.protein)
        for i, wt_aa in enumerate(mutable_seq):
            for mutant_aa in mutant_aa_list:
                mutant_title = wt_aa + str(i + 1) + mutant_aa
                wt_aa_at_i = mutable_seq[i]
                if not mutant_aa == wt_aa_at_i:
                    mutable_seq[i] = mutant_aa
                    titleSeqDict_w_mutants[mutant_title] = str(mutable_seq)
                    mutable_seq[i] = wt_aa_at_i
        return titleSeqDict_w_mutants

    # Writes the mutants out to fastafiles and/or csv files and/or txt files.
    # These can be written in one file containing all mutants or one file per mutant.
    # The fastafiles are written to /input_data/<fastafilename>/mutants/.
    # The csvfiles and txtfiles are written to /output_data/<fastafilename/sequences/.
    # The reason for fastafile mutants being written to an input folder is that these sequences are generated as direct
    # inputs for the mutation operation.
    #
    # title_titleSeq_w_mutants      Dictionary  Title of wild-type associated to every mutant title:sequence
    #                                           starting with the wild-type title:sequence pair.
    # path_input                    String      Absolute path of input_data dir (where fastafiles will be written).
    # write_1_fasta_only            Boolean     True to write one fastafile containing all mutants, separated by \n.
    # write_fasta_per_mut           Boolean     True to write one fastafile per mutant.
    # path_output                   String      Absolute path of input_data dir (where txt & csv seqs will be written).
    # write_csv                     Boolean     True to write 1 csv file for wt & mutants. False by default.
    # write_txt                     Boolean     True to write 1 txt file for wt & mutants. False by default.
    def _write_mutants(self, title_titleSeq_w_mutants, path_input, write_1_fasta_only, write_fasta_per_mut,
                       path_output, write_csv=False, write_txt=False):
        path_all_mutants_fastafile = None
        path_seqs_w_mutants = None
        for wt_title, title_seq in title_titleSeq_w_mutants.items():
            path_input_fastas_title_mutants = os.path.join(path_input, 'fastas', wt_title, 'mutants')
            try:
                os.makedirs(path_input_fastas_title_mutants)
            except FileExistsError:
                print('Part of all of path already exists. This is absolutely fine.')
            if write_1_fasta_only:
                path_all_mutants_fastafile = open(os.path.join(path_input_fastas_title_mutants,
                                                               wt_title + '_mutants.fasta'), 'w')
            if write_csv or write_txt:
                path_output_title_seqs = os.path.join(path_output, wt_title, 'sequences')
                try:
                    os.makedirs(path_output_title_seqs)
                except FileExistsError:
                    print('Part of all of path already exists. This is absolutely fine.')
                if write_csv:
                    path_seqs_w_mutants = open(os.path.join(path_output_title_seqs, wt_title + '_mutants.csv'), 'w')
                    path_seqs_w_mutants.write(wt_title + ':' + title_seq[wt_title] + ',')
                if write_txt:
                    path_seqs_w_mutants = open(os.path.join(path_output_title_seqs, wt_title + '_mutants.txt'), 'w')
                    path_seqs_w_mutants.write(wt_title + ':' + title_seq[wt_title] + '\n')
            for mut_title, mut_seq in title_seq.items():
                if mut_title is not wt_title:
                    if path_all_mutants_fastafile is not None:
                        path_all_mutants_fastafile.write('>' + mut_title + '\n' + mut_seq + '\n')
                    if path_seqs_w_mutants is not None:
                        if write_csv:
                            path_seqs_w_mutants.write(mut_title + ':' + mut_seq + ',')
                        if write_txt:
                            path_seqs_w_mutants.write(mut_title + ':' + mut_seq + '\n')
                    if write_fasta_per_mut:
                        path_mutant_fastafile = open(os.path.join(path_input_fastas_title_mutants,
                                                                  mut_title + '.fasta'), 'w')
                        path_mutant_fastafile.write('>' + mut_title + '\n' + mut_seq)
                        path_mutant_fastafile.close()
            if path_all_mutants_fastafile is not None:
                path_all_mutants_fastafile.close()
            if path_seqs_w_mutants is not None:
                path_seqs_w_mutants.close()

    # Splits a fastafile's text contents into a dictionary. The key is fastafile's title, value is (amino acid)
    # sequence.
    # E.g.
    # ">1_A\n'RVYLTFDELRETKTSEYFSLSHHPLDYRILLMDEDQDRIYVG...' etc"
    #
    # is saved as
    #
    # {
    #   '1_A': <seq of 1_A>,
    #   '2_A': <seq of 2_A>,
    #    etc..
    # }
    #
    #
    # The fastafile is assumed to be either two lines: 1st line ">title", 2nd line sequence, but the method also allows
    # there to be either no title (i.e. sequence only) or an empty title (">" with no title after it). In these both
    # scenarios, the title is created based on the name of the file itself.
    #
    # path_fastafiles   String OR List of Strings   Absolute path to the fasta file, including the file (incl. .fasta)
    #                                               e.g. ~/../input_data/fastas/1_A/1_A.fasta
    #
    # Returns title:sequence dictionary.
    def make_titleSeqDict_from_fastafile(self, path_fastafiles):
        title_sequence_dict = {}
        if isinstance(path_fastafiles, str):
            path_fastafiles = [path_fastafiles]
        for path_fastafile in path_fastafiles:
            with open(path_fastafile, 'r') as path_fastafile_opened:
                fastafile_lines = path_fastafile_opened.readlines()
                if len(fastafile_lines) > 2 or len(fastafile_lines) < 1:
                    raise ValueError('There is either an extra unexpected carriage return or no sequence data at all')
                elif len(fastafile_lines) == 2:
                    title = GUM._remove_prefix_and_suffix(fastafile_lines[0],  ">", "\n")
                if title == "" or len(fastafile_lines) == 1:
                    title = os.path.splitext(path_fastafiles)[0]
                fasta_seq = fastafile_lines[1] if len(fastafile_lines) == 2 else fastafile_lines[0]
                title_sequence_dict[title] = fasta_seq
        return title_sequence_dict

    # Converts a title:sequence dictionary to a title: {title:seq dictionary} dictionary. One use of this is to add
    # 19 more titles per residue for different titles which represent point mutants.
    # E.g.
    #
    # {
    #   '1_A': <seq of 1_A>,
    #   '2_A': <seq of 2_A>,
    #    etc..
    # }
    #
    # is converted to
    #
    # {
    #   '1_A':
    #           {'1_A'      : <seq of 1_A>,
    #            '1_A_R1A'  : <seq of 1_A_R1A>, ... etc },
    #   '2_A':
    #           {'2_A'      : <seq of 2_A>,
    #            '2_A_P1A'  : <seq of 2_A_P1A>, ... etc },
    #    etc..
    # }
    #
    # title_sequence_dict   Dictionary   Fasta title (i.e. sequence id) is key. Sequence in FASTA format is value.
    #
    # Returns dictionary of title: title:sequence dictionaries.
    def convert_titleSeqDict_to_titleTitleSeqDictDict(self, title_sequence_dict):
        title_titleSeqDict_dict = {}
        for title in title_sequence_dict.keys():
            title_titleSeqDict_dict[title] = title_sequence_dict
        return title_titleSeqDict_dict

    # Combines the absolute path given to every fasta file in the given list. Builds ayn missing subdirectories such as
    # /fastas/<fastaname>
    # This separation of path to input directory and fasta files was done to provide more flexibility (it makes it
    # more testable too).
    #
    # path_input        String              Absolute path of directory holding list of fasta files.
    # fastafiles        String OR List      Target fasta files (including extensions).
    #
    # Returns list of fasta files with absolute path to them, e.g. ~/.../input_data/fastas/1_A/1_A.fasta
    def build_complete_paths_for_fastafiles(self, path_input, fastafiles):
        if isinstance(fastafiles, str):
            fastafiles = [fastafiles]
        return [path_input + '/fastas/' + fastafile.split('.')[0] + '/' + fastafile for fastafile in fastafiles]
