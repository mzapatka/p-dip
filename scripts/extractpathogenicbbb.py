#!/usr/bin/env python

# !/usr/bin/env python
# Marc Zapatka
# extract bam reads of potential pathogenic origin and write as fastq

import pysam
import Bio.Seq
import warnings
import time
import re
import gzip
import zlib
import logging
import gc


def strRevComp(s):
    # transforms a DNA sequence given in a string to its reverse complement using Biopython functions
    return (Bio.Seq.reverse_complement(s))


def deleteread(alignedRead, outfile1, outfile2):
    global r1
    global r2
    global dr1
    global dr2

    if alignedRead.is_read1:
        n = alignedRead.qname
        if n in r2:
            if alignedRead.is_reverse:
                s= "@" + n  + '\n' + strRevComp(alignedRead.seq) + '\n' + "+" +'\n' + pysam.qualities_to_qualitystring( alignedRead.query_qualities[::-1]) + '\n'
            else: s= "@" + n  + '\n' + alignedRead.seq + '\n' + "+" + '\n' +pysam.qualities_to_qualitystring( alignedRead.query_qualities)+ '\n'
            outfile1.write(s.encode())
            outfile2.write(zlib.decompress(r2[n]))
            del r2[n]
        elif n in dr2:
            del dr2[n]
        else:
            if alignedRead.is_reverse:
                s= "@" + n  + '\n' + strRevComp(alignedRead.seq) + '\n' + "+" +'\n' +   pysam.qualities_to_qualitystring( alignedRead.query_qualities[::-1]) + '\n'
            else: s= "@" + n  + '\n' + alignedRead.seq + '\n' + "+" + '\n' +pysam.qualities_to_qualitystring( alignedRead.query_qualities) + '\n'
            dr1[n]= zlib.compress(s.encode(),1)
    elif alignedRead.is_read2:
        n = alignedRead.qname
        if alignedRead.is_reverse:
             s=  "@" + n + '\n' + strRevComp(alignedRead.seq) + '\n' + "+" +'\n' +  pysam.qualities_to_qualitystring( alignedRead.query_qualities[::-1]) + '\n'
        else: s=  "@" + n  + '\n' + alignedRead.seq + '\n' + "+" + '\n' +pysam.qualities_to_qualitystring( alignedRead.query_qualities) + '\n'
        if n in r1:
            if alignedRead.is_reverse:
                s=  "@" + n + '\n' + strRevComp(alignedRead.seq) + '\n' + "+" +'\n' +   pysam.qualities_to_qualitystring( alignedRead.query_qualities[::-1]) + '\n'
            else: s=  "@" + n  + '\n' + alignedRead.seq + '\n' + "+" + '\n' +pysam.qualities_to_qualitystring( alignedRead.query_qualities) + '\n'
            outfile1.write(zlib.decompress(r1[n]))
            outfile2.write(s.encode())
            del r1[n]
        elif n in dr1:
            del dr1[n]
        else:
            if alignedRead.is_reverse:
                s=  "@" + n + '\n' + strRevComp(alignedRead.seq) + '\n' + "+" +'\n' +  pysam.qualities_to_qualitystring( alignedRead.query_qualities[::-1]) + '\n'
            else: s=  "@" + n  + '\n' + alignedRead.seq + '\n' + "+" + '\n' +pysam.qualities_to_qualitystring( alignedRead.query_qualities)  + '\n'
            dr2[n]= zlib.compress(s.encode(),1)


def saveread(alignedRead, outfile1, outfile2):
    global r1
    global r2
    global dr1
    global dr2

    if alignedRead.is_read1:
        n = alignedRead.qname
        if alignedRead.is_reverse: s=  "@" + n  + '\n' + strRevComp(alignedRead.seq) + '\n' + "+" +'\n' +   pysam.qualities_to_qualitystring( alignedRead.query_qualities[::-1]) + '\n'
        else: s= "@" + n + '\n' + alignedRead.seq + '\n' + "+" + '\n' +pysam.qualities_to_qualitystring( alignedRead.query_qualities) + '\n'
        if n in r2:
            outfile1.write(s.encode())
            outfile2.write(zlib.decompress(r2[n]))
            del r2[n]
        elif n in dr2:
            outfile1.write(s.encode())
            outfile2.write(zlib.decompress(dr2[n]))
            del dr2[n]
        else:
            r1[n] = zlib.compress(s.encode(), 1)

    elif alignedRead.is_read2:
        n = alignedRead.qname
        if alignedRead.is_reverse: s=  "@" + n + '\n' + strRevComp(alignedRead.seq) + '\n' + "+" +'\n' +  pysam.qualities_to_qualitystring( alignedRead.query_qualities[::-1]) + '\n'
        else: s= "@" + n  + '\n' + alignedRead.seq + '\n' + "+" + '\n' +pysam.qualities_to_qualitystring( alignedRead.query_qualities) + '\n'
        if n in r1:
            outfile1.write(zlib.decompress(r1[n]))
            outfile2.write(s.encode())
            del r1[n]
        elif n in dr1:
            outfile1.write(zlib.decompress(dr1[n]))
            outfile2.write(s.encode())
            del dr1[n]
        else:
            r2[n] = zlib.compress(s.encode(), 1)


def writepair(ar1, ar2, outfile1, outfile2,n):
    if ar1.is_reverse:
        s = "@" + n + '\n' + strRevComp(ar1.seq) + '\n' + "+" + '\n' + pysam.qualities_to_qualitystring(ar1.query_qualities[::-1]) + '\n'
        outfile1.write(s.encode())
        s = "@" + n + '\n' + ar2.seq + '\n' + "+" + '\n' + pysam.qualities_to_qualitystring(ar2.query_qualities) + '\n'
        outfile2.write(s.encode())
    else:
        s = "@" + n + '\n' + ar1.seq + '\n' + "+" + '\n' + pysam.qualities_to_qualitystring(ar1.query_qualities) + '\n'
        outfile1.write(s.encode())
        s = "@" + n + '\n' + strRevComp(ar2.seq) + '\n' + "+" + '\n' + pysam.qualities_to_qualitystring(ar2.query_qualities[::-1]) + '\n'
        outfile2.write(s.encode())


def evalreadpair(ar1, ar2, outfile1, outfile2,n,r,M):
    if ar1.is_unmapped or ar2.is_unmapped:
        writepair(ar1, ar2, outfile1, outfile2,n)
        return
    else:
        t=ar1.tid
        if (t != -1):
            if (r[t] == "NC_007605"):
            #maps to herpesvirus4
                writepair(ar1, ar2, outfile1, outfile2,n)
                return
            else:
            #CIGAR match
            #identify largest M string below   31M    # in 20M-30M
                match = list(map(int, M.findall(ar1.cigarstring)))
                if (match!=[]):
                    m = max(match)
                    if m <= 30: # and m >= 20: write out all reads with matches shorter than 31
                        writepair(ar1, ar2, outfile1, outfile2,n)
                        return

    t=ar2.tid
    if (t != -1):
        if (r[t] == "NC_007605"):
            #maps to herpesvirus4
            writepair(ar1, ar2, outfile1, outfile2,n)
            return
        else:
            #CIGAR match
            #identify largest M string in 20M-30M
            match = list(map(int, M.findall(ar2.cigarstring)))
            if match!=[]:
                m = max(match)
                if m <= 30: # and m >= 20: write out all reads with matches shorter than 31
                    writepair(ar1, ar2, outfile1, outfile2,n)
                    return


def evalsingle(samfile,alignedRead, outfile1, outfile2,M):
    #skip supplementary alignments
    #if alignedRead.flag > 2048:
    #    continue
    # DONE by bamcollate
    #save unmapped reads
    if alignedRead.is_unmapped:
        saveread(alignedRead, outfile1, outfile2)
    elif alignedRead.tid != -1:
        #mapped read with reference sequence
        if samfile.getrname(alignedRead.tid) == "NC_007605":
            #maps to herpesvirus4
            saveread(alignedRead, outfile1, outfile2)
        else:
            #CIGAR match
            #identify largest M string in 20M-30M
            match = list(map(int, M.findall(alignedRead.cigarstring)))
            if len(match) > 0:
                m = max(match)
                if m <= 30: # and m >= 20: write out all reads with matches shorter than 31
                    saveread(alignedRead, outfile1, outfile2)
                else:
                    deleteread(alignedRead, outfile1, outfile2)
            else:
                deleteread(alignedRead, outfile1, outfile2)

    else:
        deleteread(alignedRead, outfile1, outfile2)


def performAnalysis(samfile, outfile1, outfile2, unpaired):
    global r1
    global r2
    nread = 0

    r=samfile.references
    M=re.compile(r'(\d+)M')
    while 1:
        try:
            alignedRead1 = next(samfile)
        except StopIteration:
            break
        try:
            alignedRead2 = next(samfile)
        except StopIteration:
            evalsingle(samfile,alignedRead1, outfile1, outfile2,M)
            break
        nread = nread + 2
        if nread % 100000000 == 0:
            logging.info("Processed " + str(nread) + " reads")
            gc.collect()
        n=alignedRead1.qname
        if  n == alignedRead2.qname:
            evalreadpair(alignedRead1, alignedRead2, outfile1, outfile2,n,r,M)
            continue
        #decide whether we would like to keep read
        evalsingle(samfile,alignedRead1, outfile1, outfile2,M)
        evalsingle(samfile,alignedRead2, outfile1, outfile2,M)

        #write remaining reads without pairs in
    for key in r1:
        unpaired.write(zlib.decompress(r1[key]))
    for key in r2:
        unpaired.write(zlib.decompress(r2[key]))

    samfile.close()
    outfile1.close()
    outfile2.close()
    unpaired.close()


if __name__ == '__main__':
    # constructing command line parser
    import optparse

    parser = optparse.OptionParser()
    parser.add_option('--fastqFilePrefix', action='store', type='string', dest='fastqFilePrefix',
                      help='Specify the name of the fastqFile prefix - the _R[12]_sequence.txt extension will be added',
                      default='')
    parser.add_option('--log', action='store', type='string', dest='loglevel', help='loglevel', default='WARN')
    (options, args) = parser.parse_args()
    if (options.fastqFilePrefix == ''):
        print
        "Mandatory parameters missing or wrong. Program will terminate now."
        print
        "\nYour parameter settings:"
        print
        options
        raise SystemExit
    loglevel = options.loglevel
    #read sam file from stdin
    samfile = pysam.AlignmentFile("-", "r")
    outfile1 = gzip.open(options.fastqFilePrefix + '_R1_sequence.gz', "w")
    outfile2 = gzip.open(options.fastqFilePrefix + '_R2_sequence.gz', "w")
    unpaired = gzip.open(options.fastqFilePrefix + '_sequence_unpaired.gz', "w")
    logging.basicConfig(level=loglevel.upper(), format='%(asctime)s %(levelname)s %(message)s')
    logging = logging.getLogger(__name__)
    #store reads until pair is identified
    r1 = {}
    r2 = {}
    dr1 = {}
    dr2 = {}
    performAnalysis(samfile, outfile1, outfile2, unpaired)

