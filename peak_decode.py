# Gianluca Croso, Heather Han, Nitin Kumar, Eric Walker
# EN.600.439 Computational Genomics
# Final Project: peak_decode.py

import sys
import fasta


def decode(fin, gn, chr_b, pos_b, dic):
    enc = open(fin, 'rb')
    encoded = enc.read(1000)
    bitrep = bits(encoded)
    with open(fin[:-4] + 'p.sam', 'w') as out:
        head = "@HD\tVN:1.5\tSO:unsorted\n"
        fai = open(gn + ".fai")
        for line in fai:
            line = line.strip().split()
            seq = line[0]
            ln = line[1]
            head += "@SQ\tSN:" + seq + "\tLN:" + ln + "\n"
        out.write(head)
        i = 0
        j = 0
        while i < len(bitrep) - 8:
            flag = flag_convert(int(bitrep[i: i + 6], 2))
            i += 6
            chrom = dic[int(bitrep[i: i + chr_b], 2)]
            i += chr_b
            pos = int(bitrep[i: i + pos_b], 2)
            i += pos_b
            ln = int(bitrep[i: i + 8], 2)
            i += 8
            write_read(out, gn, flag, chrom, pos, ln, j)
            j += 1
            rebuf = int(i/8)
            encoded = enc.read(rebuf)
            bitrep = bitrep[8*rebuf:] + bits(encoded)
            i -= 8*rebuf
    enc.close()
    return


def flag_convert(enc_flag):
    flag = 2*(enc_flag & 1)/1
    flag += 16*(enc_flag & 2)/2
    flag += 32*(enc_flag & 4)/4
    flag += 64*(enc_flag & 8)/8
    flag += 128*(enc_flag & 16)/16
    flag += 1024*(enc_flag & 32)/32
    return int(flag)


def write_read(out, gn, flag, chrom, pos, ln, id):
    read = "r" + str(id) + '\t' + str(flag) + '\t' + chrom + '\t' + str(pos)
    read += "\t255\t" + str(ln) + "M\t" + "*\t0\t0\t"
    read += fasta.fasta_driver(gn, chrom, pos-1, ln) + '\t'
    read += 'D'*ln + '\n'
    out.write(read)
    return


# Rotate left function
# http://www.falatic.com/index.php/108/python-and-bitwise-rotation
rol = lambda val, r_bits, max_bits: \
    (val << r_bits % max_bits) & (2**max_bits-1) | \
    ((val & (2**max_bits-1)) >> (max_bits-(r_bits % max_bits)))


def bits(string):
    bitrep = ""
    for i in range(0, len(string)):
        bitstring = string[i]
        mask = 1 << 7
        for j in range(0, 8):
            if bitstring & mask != 0:
                bitrep += '1'
            else:
                bitrep += '0'
            bitstring = rol(bitstring, 1, 8)
    return bitrep
