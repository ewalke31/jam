# Gianluca Croso, Heather Han, Nitin Kumar, Eric Walker
# EN.600.439 Computational Genomics
# Final Project: jam.py

import sys
import math
import fasta
import var_decode
import peak_decode
import var_encode
import peak_encode


def main(argv):
    # Check that the number of arguments in the line matches the expected
    if len(argv) != 5:
        sys.stderr.write("Usage: python3 jam.py <commands>" +
                         " input_file ref_genome\n")
        sys.stderr.write("commands:\n" +
                         "    e or d for encode or decode (action)\n" +
                         "    p or v for peaks or variants (application)\n")
        exit(1)

    # Grab parameters
    action = argv[1]
    application = argv[2]
    fin = argv[3]
    ref = argv[4]
    # pre-process the fasta file to get list of the chromosomes & their lengths
    # list printed in .fai file
    fasta.gen_fai(ref)
    chrdict = {}
    chrbits = 1
    posbits = 1
    # use pre-processed fasta .fai file gather info
    with open(ref + '.fai') as fai:
        i = 0
        # For each chromosome in .fai file
        for line in fai:
            # grab name and length
            line = line.strip().split()
            seq = line[0]
            ln = line[1]
            # if encoding associate read name with a number
            if action == 'e':
                chrdict[seq] = i
            # if decoding associate a number with read name
            elif action == 'd':
                chrdict[i] = seq
            # make sure using the least # of bits to handle the #
            # of chromosomes
            if i == (2**chrbits):
                chrbits += 1
            i += 1
            # also make sure the # of position bits can handle
            # the lengths of the reads
            lgln = int(math.ceil(math.log(float(ln), 2)))
            if lgln > posbits:
                posbits = lgln
    # check if encoding
    if action == 'e':
        # if so, see whether to run peak or variant (or if error)
        if application == 'p':
            peak_encode.encode(fin, chrbits, posbits, chrdict)
        elif application == 'v':
            var_encode.encode(fin, ref, chrbits, posbits, chrdict)
        else:
            sys.stderr.write("Invalid application (not implemented)." +
                             " See Usage.\n")
    # else if decoding
    elif action == 'd':
        # check if peak or variant decoding
        if application == 'p':
            peak_decode.decode(fin, ref, chrbits, posbits, chrdict)
        elif application == 'v':
            var_decode.decode(fin, ref, chrbits, posbits, chrdict)
        else:
            sys.stderr.write("Invalid application (not implemented)." +
                             " See Usage.\n")
    else:
        sys.stderr.write("Invalid action (must be e or d). See Usage.\n")
    return


if __name__ == "__main__":
    main(sys.argv)
