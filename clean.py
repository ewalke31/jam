# Gianluca Croso, Heather Han, Nitin Kumar, Eric Walker
# EN.600.439 Computational Genomics
# Final Project: clean.py

import sys


def main(argv):
    valid = ["chrX", "chrY", "chrM", "chr4", "chr2L", "chr2R",
             "chr3L", "chr3R"]
    valid1 = ["SN:chrX", "SN:chrY", "SN:chrM", "SN:chr4", "SN:chr2L",
              "SN:chr2R", "SN:chr3L", "SN:chr3R"]
    if len(argv) < 3:
        print("Usage: python clean.py input.sam output.sam")
        return
    orig = argv[1]
    altered = argv[2]
    with open(orig, 'r') as inp:
        with open(altered, 'w') as out:
            for line in inp:
                linex = line.strip().split()
                if linex[1] in valid1 or linex[2] in valid:
                    out.write(line)
    return


if __name__ == "__main__":
    main(sys.argv)
