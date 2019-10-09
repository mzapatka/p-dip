#!/usr/bin/env Rscript
#filtering blast hits
#Rscript --vanilla blast_filter.R

library("optparse")
library("data.table")
library("Rsamtools")
library("parallel")
library("plyr")

option_list = list(
  make_option(
    c("-b", "--blast"),
    type = "character",
    default = NULL,
    help = "blasthits file name",
    metavar = "character"
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    default = "out.txt",
    help = "output file name [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "default number of threads [default= %default]",
    metavar = "number"
  ),
  make_option(
    c("-s", "--sequencebam"),
    type = "character",
    default = NULL,
    help = "bamfile with reads mapped to assembled contigs [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-c", "--contaminants"),
    type = "character",
    default = NULL,
    help = "tab separated file with information on contaminants/other sequences [default= %default]\nFormat: contaminanttype\\tNCBItaxid",
    metavar = "character"
  ),
  make_option(
    c("-i", "--inputreads"),
    type = "integer",
    default = 1,
    help = "number of inputreads for PMER calculation [default= %default]",
    metavar = "number"
  )
  
  
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$blast) | is.null(opt$sequencebam)) {
  print_help(opt_parser)
  stop(
    "At least two argument must be supplied (input blast hits file, aligned contigs BAM file)",
    call. = FALSE
  )
}

if (!file.exists(opt$sequencebam) | !file.exists(opt$blast)) {
  print_help(opt_parser)
  stop("Input BAM file is not existing",
       call. = FALSE)
}

no_of_cores = detectCores()
if (opt$threads > no_of_cores)
  opt$threads = no_of_cores

summaryFunction <- function(region, bamFile, ...)
{
  param <- ScanBamParam(
    what = c('pos', 'qwidth'),
    which = region,
    flag = scanBamFlag(isUnmappedQuery = FALSE)
  )
  x <- scanBam(bamFile, ..., param = param)[[1]]
  #count reads overlapping with at least one base
  return(nreads = length(x$pos))
}

my.identifyhits = function(tile, hl) {
  m = findOverlaps(tile, hl)
  mh = as.data.table(mcols(hl))
  mh = mh[subjectHits(m), ]
  rm(hl)
  mh[, pos := queryHits(m)]
  rm(m)
  #.Sd has a large overhead using .I instead
  setkey(mh, pos)
  #identify and keep only hits with max bitscore per position
  h = mh[mh[, bitscore == max(bitscore), by = pos]$V1, ]
  #remove pos column
  h[, pos := NULL]
  h = unique(h)
  return(h)
}

my.annotatetax = function(taxid,
                          sep = ";",
                          annolist,
                          mc.cores) {
  ts=strsplit(taxid,sep)
  u=unique(unlist(ts)) 
  u=join(data.frame(taxid=u),annolist)  
  
  res=lapply(ts,function(x,u,sep){
    m=u[match(x,u$taxid),"type"]
    m=m[!is.na(m)]
    if(length(m)==0)return("")
    else return(paste0(m,sep=sep))
  },u=u,sep=sep)
  return(unlist(res))
}




X = fread(cmd = paste("gunzip -c", opt$blast))
colnames(X) = c(
  "qseqid",
  "sseqid",
  "sacc",
  "saccver",
  "pident",
  "length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore",
  "staxids",
  "sscinames",
  "scomnames",
  "sblastnames",
  "sskingdoms",
  "stitle",
  "qcovs",
  "qcovhsp"
)
X[, hitid := c(1:nrow(X))]
#x = X[, c("qstart", "qend", "bitscore", "hitid", "qseqid", "length")]
#scolnames(x) = gsub("qseqid", "id", colnames(x))

hl = GRanges(
  seqnames = as.character(X$qseqid),
  ranges = IRanges(
    start = pmin(X$qstart, X$qend),
    end = pmax(X$qstart, X$qend)
  ),
  mcols = X[, c("qseqid", "qstart", "qend", "bitscore", "hitid", "length")]
)
# rm(X)
# gc()
colnames(mcols(hl)) = c("id", "qstart", "qend", "bitscore", "hitid", "length")

sl = !duplicated(as.vector(seqnames(hl)))
seqlengths = mcols(hl)$length[sl]
names(seqlengths) = seqnames(hl)[sl]
rm(sl)
tile = tileGenome(seqlengths, tilewidth = 1)

h = mclapply(unique(seqnames(hl)), function(s, tile, hl) {
  rm(X)
  gc()
  my.identifyhits(tile[seqnames(tile) == s], hl[seqnames(hl) == s])
}, tile = tile, hl = hl, mc.cores = opt$threads)
h = rbindlist(h)
h = unique(h)

#abundance from trinity = average kmer abundance
h[, abundance := unlist(mclapply(h$id, function(y) {
  as.numeric(sub("\\S+;(\\S+)$", "\\1", y))
}, mc.cores = opt$threads))]

hr = GRanges(seqnames = as.character(h$id),
             ranges = IRanges(
               start = pmin(h$qstart, h$qend),
               end = pmax(h$qstart, h$qend)
             ))
targetreads = unlist(
  mclapply(seq_along(hr), function(i, sequencebam, hr) {
    summaryFunction(hr[i], sequencebam)
  }, sequencebam = opt$sequencebam, hr = hr, mc.cores = opt$threads)
)

h[, targetreads := targetreads]
h[, pmer := targetreads / opt$inputreads * 10 ^ 6]



#add annotation based on hitid
res = merge(h, X[, c(colnames(X)[!colnames(X) %in% c(colnames(h), "qseqid")], "hitid"), with =
                   FALSE], by.x = "hitid", by.y = "hitid")
rm(X)
gc()
res = res[, hitid := NULL]
#annotate coverage from bam file of reads aligned to contigs

#filter contaminants based on list with ncbi taxids
otherseq=fread(opt$contaminants,sep="\t",header = TRUE)
#merge multiple annotations
otherseq=otherseq[,.(type=paste0(type,collapse=";")),by=taxid]

res[, other := my.annotatetax(taxid = as.character(res$staxids),
                              sep = ";",
                              annolist = otherseq,
                              mc.cores=opt$threads)]
res[, no := 1:nrow(res)]

options(warn = -1) #max seems to be calculate even if there are no others ?
foo = res[, .(
  no = no,
  contaminant = ifelse(any(
    other %in% c(
      "artificial_sequences",
      "insertion_sequences",
      "plasmids",
      "transposons"
    )
  ), TRUE, FALSE),
  synvirus = ifelse(any(
    other %in% c(
      "plasmids;transposons;midivariant_sequence;insertion_sequences;artificial_sequences;synthetic_viruses",
      "artificial_sequences;synthetic_viruses"
    )
  ), TRUE, FALSE),
  species = staxids,
  max_pmer_contaminant = ifelse(any(
    other %in% c(
      "artificial_sequences",
      "insertion_sequences",
      "plasmids",
      "transposons"
    )
  ),
  max(pmer[other %in% c("artificial_sequences",
                        "insertion_sequences",
                        "plasmids",
                        "transposons")], na.rm = TRUE),
  0),
  max_contaminant_coverage = ifelse(any(
    other %in% c(
      "artificial_sequences",
      "insertion_sequences",
      "plasmids",
      "transposons"
    )
  ),
  as.double(max(
    targetreads[other %in% c("artificial_sequences",
                             "insertion_sequences",
                             "plasmids",
                             "transposons")], na.rm = TRUE
  )),
  0)
), by = id]

options(warn = 0)
setorder(foo, no)
res = cbind(res, foo[, c("contaminant",
                         "max_contaminant_coverage",
                         "max_pmer_contaminant")])
res[,nreads := opt$inputreads]
message(Sys.time(), " Identified ", nrow(res), " hits.")


fwrite(
  res,
  file = opt$out,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)
