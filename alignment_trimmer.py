#!/usr/bin/env python3


from __future__ import division, print_function
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import argparse
import os


def gap_proportion(msa_column):
    gap_counter = 0
    for aa_or_codon in msa_column:
        if aa_or_codon in ('-', '---'):
            gap_counter += 1
    return gap_counter / len(msa_column)


class Alignment():
    def __init__(self, alignment_fasta, maxgap,
                 trim_positions, mask, is_cds=False):
        self._is_cds = is_cds
        self._letter_n = 3 if is_cds else 1
        self._maxgap = maxgap
        self._fasta = []
        self._trimmed_fasta = []
        with open(alignment_fasta, 'r') as input_handle:
            for fasta_record in SeqIO.parse(input_handle, 'fasta'):
                self._fasta.append(fasta_record)
                self._trimmed_fasta.append([])
        self._first_seq = self._fasta[0]
        self._assert_fasta_is_aligned()
        self._set_mask_char(mask)
        self._parse_trim_positions(trim_positions)
        fname, fextension = os.path.splitext(alignment_fasta)
        self._output_filepath = '{}_{}{}'.format(
            fname,
            "masked" if self._mask_char else "trimmed",
            fextension,
        )
        return

    def _parse_trim_positions(self, trim_pos_list):
        self._trim_positions = set()
        if trim_pos_list:
            for trim_pos in trim_pos_list:
                if '-' in trim_pos:  # it's a range!
                    start_raw, stop_raw = trim_pos.split('-')
                    start = int(start_raw) - 1
                    stop = int(stop_raw)
                    if start >= stop:
                        raise Exception('Invalid range:', trim_pos)
                    self._trim_positions |= set(range(start, stop))
                else:
                    pos = self._str_to_pos(trim_pos)
                    self._trim_positions.add(pos)
        return

    def _str_to_pos(self, pos_str):
        try:
            pos = int(pos_str) - 1  # account for zero-indexing
        except ValueError:
            raise Exception(pos_str, 'is not a valid position!')
        self[pos]  # raises IndexError when pos is out of range of the MSA
        return pos

    def _set_mask_char(self, mask_param):
        if mask_param:
            if len(mask_param) != 1:
                raise Exception('Invalid --mask character:', mask_param)
            self._mask_char = mask_param * self._letter_n
        else:
            self._mask_char = False
        return

    def _assert_fasta_is_aligned(self):
        for fa_entry in self._fasta[1:]:
            if len(fa_entry) != len(self._first_seq):
                raise Exception('The sequences {} and {} have different '
                                'lengths. Is your fasta really aligned'
                                '?'.format(fa_entry.id, self._first_seq.id))
        return

    def __getitem__(self, i):
        ''' returns the n-th column of the multifasta, works for
        CDS and protein fastas '''
        n = self._letter_n  # the number of characters per position
        # three for CDS sequences (codons), one for protein sequences
        start = i * n
        end = start + n
        if start < 0 or end > len(self._first_seq):
            raise IndexError('Cannot access MSA position {} - index out '
                             'of range'.format(i))
        return tuple(str(fa.seq[start:end]) for fa in self._fasta)

    def __iter__(self):
        ''' column-based iteration over MSA fasta, works for
        CDS and protein fastas '''
        for i in range(int(len(self._first_seq) / self._letter_n)):
            yield self[i]

    def _trim_column(self, msa_column, col_index):
        ''' trim an MSA column if the number of gaps exceeds the user-specified
        threshold, or if it is within the user-specified trim-range '''
        col_str = ' '.join(msa_column)
        gap_prop = gap_proportion(msa_column)

        if (col_index in self._trim_positions) or (gap_prop > self._maxgap):
            if self._mask_char:
                action = ' - masked'
            else:
                action = ' - trimmed'
            trim_column = True
        else:
            action = ''
            trim_column = False

        for row_i, aa_or_codon in enumerate(msa_column):
            if trim_column:  # trim (or mask) the column
                if self._mask_char:  # mask it
                    self._trimmed_fasta[row_i].append(self._mask_char)
                else:  # don't add the column at all
                    break
            else:  # keep the column as it was
                self._trimmed_fasta[row_i].append(aa_or_codon)
        print('{: >3} {} ({:.1%} gaps){}'.format(
              col_index + 1, col_str, gap_prop, action))
        return

    def _trim_fasta(self):
        for col_index, msa_column in enumerate(self):
            self._trim_column(msa_column, col_index)
        self._trimmed_seqs = []
        for i, trimmed_seq_list in enumerate(self._trimmed_fasta):
            trimmed_seq = ''.join(trimmed_seq_list)
            input_seq = self._fasta[i]
            out_record = SeqRecord(
                Seq(trimmed_seq, SingleLetterAlphabet),
                id=input_seq.id,
                description=input_seq.description
            )
            self._trimmed_seqs.append(out_record)

    def dump_trimmed_fasta(self):
        self._trim_fasta()
        SeqIO.write(self._trimmed_seqs, self._output_filepath, 'fasta')
        if self._mask_char:
            print('  Wrote "{}"-masked fasta to {}\n'.format(
                self._mask_char, self._output_filepath))
        else:
            print('  Wrote trimmed fasta to {}\n'.format(self._output_filepath))
        return

    def __eq__(self, other):
        """ defines an equality test with the == operator
        returns False when a position of one alignment is a gap
        when it isn't a gap in the other alignment """
        for col_i, self_column in enumerate(self):
            other_column = other[col_i]
            for row_i, aa_or_codon in enumerate(self_column):
                if aa_or_codon in ('-', '---'):
                    if other_column[row_i] not in ('-', '---'):
                        print('\nCodon mismatch in row {} column {}'
                              '\n"{}" does not match "{}"\n'.format(
                                row_i, col_i, aa_or_codon, other_column[row_i]))
                        return False
        return True

    def __ne__(self, other):
        """ defines a non-equality test with the != operator """
        return not self.__eq__(other)


def main(pep_alignment_fasta, cds_alignment_fasta,
         trim_gappy, trim_positions, mask):
    if pep_alignment_fasta:
        pep_al = Alignment(pep_alignment_fasta, maxgap=trim_gappy, mask=mask,
                           trim_positions=trim_positions)
        pep_al.dump_trimmed_fasta()
    if cds_alignment_fasta:
        cds_al = Alignment(cds_alignment_fasta, maxgap=trim_gappy, mask=mask,
                           trim_positions=trim_positions, is_cds=True)
        if pep_alignment_fasta and cds_alignment_fasta:
            if pep_al != cds_al:
                raise Exception('The two alignments have a mismatch'
                                '- see printouts above.')
            else:
                print('Peptide and CDS alignment are in agreement!\n')
        cds_al.dump_trimmed_fasta()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='A program to remove (trim) or mask columns of fasta '
        'multiple sequence alignments. Columns are removed if they exceed a '
        'user-specified "gappyness" treshold, and/or if they are specified by '
        'column index. Also works for codon alignments. The trimmed alignment '
        'is written to [your_alignment]_trimmed.fa')
    parser.add_argument('-p', '--pep_alignment_fasta')
    parser.add_argument('-c', '--cds_alignment_fasta')
    parser.add_argument('--trim_gappy', type=float, default=9999,
                        help='the maximum allowed gap proportion for an MSA '
                        'column, all columns with a higher proportion of gaps '
                        'will be trimmed (default: off)')
    parser.add_argument('--trim_positions', help='specify columns to be '
                        'trimmed from the alignment, e.g. "1-10 42" to trim the '
                        'first ten columns and the 42th column', nargs='+')
    parser.add_argument('--mask', default=False, const='-',
                        action='store', nargs='?', help='mask the alignment '
                        'instead of trimming. if a character is specified after'
                        ' --mask, that character will be used for masking '
                        '(default char: "-")')
    args = parser.parse_args()
    if args.pep_alignment_fasta or args.cds_alignment_fasta:
        main(**vars(args))
    else:
        parser.print_help()
