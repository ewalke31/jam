# Gianluca Croso, Heather Han, Nitin Kumar, Eric Walker
# EN.600.439 Computational Genomics
# Final Project: peak_encode.py

import sys


def encode(sam, chr_b, pos_b, dic):
    buff = ""  # buffer
    # output file
    fout = sam[:-3] + "pjam"
    out = open(fout, 'wb')
    with open(sam) as f:
        for line in f:
            if line[0] != '@':
                line = line.strip().split()
                # convert sam fields to bitstring
                buff += read(int(line[1]), line[2], int(line[3]),
                             len(line[9]), chr_b, pos_b, dic)
                # write complete bytes to output
                wlen = int(len(buff)/8)
                for i in range(wlen):
                    b = int(buff[i*8:i*8+8], base=2)
                    out.write(bytes([b]))
                # clear buffer
                buff = buff[wlen*8:]
        # pad with 0s and write remaining bytes
        pad = 8 - (len(buff) % 8)
        if pad != 8:
            buff += "0"*pad
    for i in range(int(len(buff)/8)):
        s = buff[i*8: (i+1)*8]
        b = int(s, base=2)
        out.write(bytes([b]))
    out.close()
    return


def read(flg, ch, pos, readlen, cb, pb, dic):
    ''' Takes flags, chromosome number, position of read, read length, the
    number of chromosome bits to write, the number of position bits to write,
    and a dictionary that maps the chromosome names to binary representations.
    Returns the compressed read as a binary sting. '''
    bitstr = ""  # output
    if (flg & 4) != 0:
        return bitstr
    flg_new = int((flg & 2)/2*1 + (flg & 16)/16*2 + (flg & 32)/32*4 +
                  (flg & 64)/64*8 + (flg & 128)/128*16 + (flg & 1024)/1024*32)
    bitstr += "{0:06b}".format(flg_new)  # write saved flags
    chm = dic[ch]
    sc = "{0:0" + str(cb) + "b}"
    bitstr += sc.format(chm)  # chromosome
    sp = "{0:0" + str(pb) + "b}"
    bitstr += sp.format(pos)  # position
    bitstr += "{0:08b}".format(readlen)  # read leangth
    return bitstr
