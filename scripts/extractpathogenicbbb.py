#!/usr/bin/env python

# !/usr/bin/env python
# Marc Zapatka
# extract bam reads of potential pathogenic origin and write as fastq

# #bamcollate2 --help
#This is biobambam2 version 2.0.8.
#biobambam2 is distributed under version 3 of the GNU General Public License.

#Key=Value pairs:

#collate=<[1]>                                  : collate pairs
#reset=<>                                       : reset alignments and header like bamreset (for collate=0,1 or 3 only, default enabled for 3)
#level=<[-1]>                                   : compression settings for output bam file (-1=zlib default,0=uncompressed,1=fast,9=best)
#filename=<[stdin]>                             : input filename (default: read file from standard input)
#inputformat=<[bam]>                            : input format: bam
#ranges=<[]>                                    : input ranges (bam input only, collate<2 only, default: read complete file)
#exclude=<[SECONDARY,SUPPLEMENTARY]>            : exclude alignments matching any of the given flags
#disablevalidation=<[0]>                        : disable validation of input data
#colhlog=<[18]>                                 : base 2 logarithm of hash table size used for collation
#colsbs=<[134217728]>                           : size of hash table overflow list in bytes
#T=<[bamcollate2_tbi-pcawg01_12868_1432644314]> : temporary file name
#md5=<[0]>                                      : create md5 check sum (default: 0)
#md5filename=<filename>                         : file name for md5 check sum (default: extend output file name)
#index=<[0]>                                    : create BAM index (default: 0)
#indexfilename=<filename>                       : file name for BAM index file (default: extend output file name)
#readgroups=[<>]                                : read group filter (default: keep all)
#mapqthres=<[-1]>                               : mapping quality threshold (collate=1 only, default: keep all)
#classes=[F,F2,O,O2,S]                          : class filter (collate=1 only, default: keep all)
#resetheadertext=[<>]                           : replacement SAM header text file for reset=1 (default: filter header in source BAM file)
#resetaux=<[1]>                                 : reset auxiliary fields (collate=0,1 only with reset=1)
#auxfilter=[<>]                                 : comma separated list of aux tags to keep if reset=1 and resetaux=0 (default: keep all)
#outputformat=<[bam]>                           : output format (bam)
#O=<[stdout]>                                   : output filename (standard output if unset)
#outputthreads=<[1]>                            : output helper threads (for outputformat=bam only, default: 1)
#inputbuffersize=<[65536]>                      : input buffer size
#replacereadgroupnames=<[]>                     : name of file containing list of read group identifier replacements (default: no replacements)

#Alignment flags: PAIRED,PROPER_PAIR,UNMAP,MUNMAP,REVERSE,MREVERSE,READ1,READ2,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
# bamcollate2 filename=/icgc/pcawg/analysis/PCAWG16/test/test_chr1-1-300000.bam outputformat=sam |less
#$BAMCOLLATE exclude=SUPPLEMENTARY outputformat=sam filename=$BAMFILE |python $EXTRACT --log=INFO --fastqFilePrefix=test12 -
#BAMCOLLATE=/icgc/pcawg/analysis/PCAWG16/scripts/biobambam2/bin/bamcollate2
#BAMFILE=/icgc/pcawg/analysis/PCAWG16/test/test_chr1-1-300000.bam
#EXTRACT=/icgc/pcawg/analysis/PCAWG16/scripts/tools/extractpathogenicbbb.py
#testfile
#time $BAMCOLLATE exclude=SUPPLEMENTARY outputformat=sam filename=$BAMFILE |python $EXTRACT --log=INFO --fastqFilePrefix=test2 -



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
            outfile1.write(s)
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
            outfile2.write(s)
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
            outfile1.write(s)
            outfile2.write(zlib.decompress(r2[n]))
            del r2[n]
        elif n in dr2:
            outfile1.write(s)
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
            match = map(int, M.findall(alignedRead.cigarstring))
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

