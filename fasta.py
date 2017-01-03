# Gianluca Croso, Heather Han, Nitin Kumar, Eric Walker
# EN.600.439 Computational Genomics
# Final Project: fasta.py
# Using Ben Langmead's code from:
# https://github.com/BenLangmead/comp-genomics-class/blob/master/notebooks/FASTA.ipynb

import sys
import re
import os


def index_fasta(fh):
    index = []
    current_short_name = None
    current_byte_offset, running_seq_length, running_byte_offset = 0, 0, 0
    line_length_including_ws, line_length_excluding_ws = 0, 0
    for ln in fh:
        ln_stripped = ln.rstrip()
        running_byte_offset += len(ln)
        if ln[0] == '>':
            if current_short_name is not None:
                index.append((current_short_name, running_seq_length,
                              current_byte_offset, line_length_excluding_ws,
                              line_length_including_ws))
            long_name = ln_stripped[1:]
            current_short_name = long_name.split()[0]
            current_byte_offset = running_byte_offset
            running_seq_length = 0
        else:
            line_length_including_ws = max(line_length_including_ws, len(ln))
            line_length_excluding_ws = max(line_length_excluding_ws,
                                           len(ln_stripped))
            running_seq_length += len(ln_stripped)
    if current_short_name is not None:
        index.append((current_short_name, running_seq_length,
                      current_byte_offset, line_length_excluding_ws,
                      line_length_including_ws))
    return index


class FastaIndexed(object):
    """ Encapsulates a set of indexed FASTA files.  Does not load the FASTA
        files into memory but still allows the user to extract arbitrary
        substrings, with the help of the index. """

    __removeWs = re.compile(r'\s+')

    def __init__(self, fafns):
        self.fafhs = {}
        self.faidxs = {}
        self.chr2fh = {}
        self.offset = {}
        self.lens = {}
        self.charsPerLine = {}
        self.bytesPerLine = {}

        for fafn in fafns:
            # Open FASTA file
            self.fafhs[fafn] = fh = open(fafn, 'r')
            # Parse corresponding .fai file
            with open(fafn + '.fai') as idxfh:
                for ln in idxfh:
                    toks = ln.rstrip().split()
                    if len(toks) == 0:
                        continue
                    assert len(toks) == 5
                    # Parse and save the index line
                    chr, ln, offset, charsPerLine, bytesPerLine = toks
                    self.chr2fh[chr] = fh
                    self.offset[chr] = int(offset)  # 0-based
                    self.lens[chr] = int(ln)
                    self.charsPerLine[chr] = int(charsPerLine)
                    self.bytesPerLine[chr] = int(bytesPerLine)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        # Close all the open FASTA files
        for fafh in self.fafhs.values():
            fafh.close()

    def has_name(self, refid):
        return refid in self.offset

    def name_iter(self):
        return self.offset.iterkeys()

    def length_of_ref(self, refid):
        return self.lens[refid]

    def get(self, refid, start, ln):
        ''' Return the specified substring of the reference. '''
        assert refid in self.offset
        if start + ln > self.lens[refid]:
            raise ReferenceOOB('"%s" has length %d; tried to get [%d, %d)' %
                               (refid, self.lens[refid], start, start + ln))
        fh, offset, charsPerLine, bytesPerLine = \
            self.chr2fh[refid], self.offset[refid], \
            self.charsPerLine[refid], self.bytesPerLine[refid]
        byteOff = offset
        byteOff += (start // charsPerLine) * bytesPerLine
        into = start % charsPerLine
        byteOff += into
        fh.seek(byteOff)
        left = charsPerLine - into
        # Count the number of line breaks interrupting the rest of the
        # string we're trying to read
        if ln < left:
            return fh.read(ln)
        else:
            nbreaks = 1 + (ln - left) // charsPerLine
            res = fh.read(ln + nbreaks * (bytesPerLine - charsPerLine))
            res = re.sub(self.__removeWs, '', res)
        return res


def fasta_driver(fn, chrom, offset, length):
    with FastaIndexed([fn]) as fa_idx:
        return fa_idx.get(chrom, offset, length)


def gen_fai(fn):
    if os.path.isfile('./' + fn + '.fai') is False:
        with open(fn) as fh:
            idx = index_fasta(fh)
        with open(fn + '.fai', 'w') as fh:
            fh.write('\n'.join(['\t'.join(map(str, x)) for x in idx]))
    return
