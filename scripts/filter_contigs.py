#read in inchworm fasta file and report only contigs with length above limit mincontiglength

import argparse

parser = argparse.ArgumentParser(description='filter contigs based on size')
parser.add_argument('-o', '--outfile', help='outfile')
parser.add_argument('-i','--infile',   help='input fasta file')
parser.add_argument('-l','--length', help='minimun contiglength')
args = parser.parse_args()


#open outfile
outfile = open(args.outfile, 'w')


from Bio import SeqIO
for seq_record in SeqIO.parse(args.infile, "fasta"):
    #print(seq_record.id)
    if len(seq_record) > int(args.length):
        r=SeqIO.write(seq_record, outfile, 'fasta')

outfile.close()
