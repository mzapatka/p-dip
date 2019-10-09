import os
import re
from snakemake.io import _wildcard_regex, regex
from collections import  namedtuple
from itertools import chain


# Problem io.walk does by default not follow symlinks therefore files are not identified if provided with full path
#snakemake / io.py
# fixed by adding followlinks
os.walk.__defaults__=(True, None, True)

def input_all():
  f=config["infolder"] + "/{pid}/alignment/{id}" + config["bamsuffix"]
  (PIDs, IDs) = glob_wildcards(f)
  print(PIDs)
  print(IDs)
  for i in range(len(PIDs)):
    f= config["logfolder"] + "/" + PIDs[i] + "/" + IDs[i] + ".done"
    print(f)
    yield(f)

rule all:
    input: lambda wc: list(input_all())

    
rule pdip_extract_fastq:
  input:
    expand("{infolder}/{{ID1}}/alignment/{{ID2}}{bamsuffix}",infolder=config["infolder"],bamsuffix=config["bamsuffix"])

  params:
    name="pdip_extract_{ID1}_,_{ID2}",
    scriptfolder=config["scriptfolder"],
    R1=config["R1"]
  threads: 3
  log:
    expand("{logfolder}/{{ID1}}/{{ID2}}_pdip_extract.log",logfolder=config["logfolder"])
  output:
    done=expand("{logfolder}/{{ID1}}/{{ID2}}_fastqextract.done",logfolder=config["logfolder"]),
    R1=expand("{outfolder}/{{ID1}}/{{ID2}}{R1}",outfolder=config["outfolder"],R1=config["R1"]),
    R2=expand("{outfolder}/{{ID1}}/{{ID2}}{R2}",outfolder=config["outfolder"],R2=config["R2"])
  shell:
    """
    set +o nounset
    touch {output.done}
    touch {output.R1}
    touch {output.R2}
    {config[envir]}
    PREFIX={output.R1}
    PREFIX=${{PREFIX/{params.R1}/}}
    bamcollate2 exclude=SUPPLEMENTARY outputformat=bam filename={input} level=1 |mbuffer -m 400M -q |samtools view -h - |mbuffer -m 400M -q |python   {params.scriptfolder}/extractpathogenicbbb.py --log=INFO --fastqFilePrefix=${{PREFIX}} -  > {log} 2>&1
    
    """
    
rule pdip_trim:
  input:
    R1=expand("{outfolder}/{{ID1}}/{{ID2}}{R1}",outfolder=config["outfolder"],R1=config["R1"]),
    R2=expand("{outfolder}/{{ID1}}/{{ID2}}{R2}",outfolder=config["outfolder"],R2=config["R2"])
  params:
    name="pdip_trim_{ID1},_{ID2}"
  threads: 3
  log:
    expand("{logfolder}/{{ID1}}/{{ID2}}_pdip_trim.log",logfolder=config["logfolder"])
  output:
    done=expand("{logfolder}/{{ID1}}/{{ID2}}_trim.done",logfolder=config["logfolder"]),
    R1p=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_R1_trim_1p.gz",outfolder=config["outfolder"])),
    R1u=expand("{outfolder}/{{ID1}}/{{ID2}}_R1_trim_1u.gz",outfolder=config["outfolder"]),
    R2p=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_R2_trim_2p.gz",outfolder=config["outfolder"])),
    R2u=expand("{outfolder}/{{ID1}}/{{ID2}}_R2_trim_2u.gz",outfolder=config["outfolder"]),
    inputreads=expand("{outfolder}/{{ID1}}/{{ID2}}_inputreads",outfolder=config["outfolder"])
  shell:
    """
    touch {output.done}
    touch {output.R1p}
    touch {output.R1u}
    touch {output.R2p}
    touch {output.R2u}
    set +o nounset
    {config[envir]}
    trimmomatic PE -threads {threads} -phred33  {input.R1} {input.R2} {output.R1p} {output.R1u} {output.R2p} {output.R2u} TRAILING:5 MINLEN:50  > {log} 2>&1
    #fastqlines /4 = inputreads
    echo $(( $(gunzip -c {input.R1} {input.R2} |wc -l) /4))  > {output.inputreads}
    """

rule pdip_alignrate:
#using only paired reads
  input:
    R1=expand("{outfolder}/{{ID1}}/{{ID2}}_R1_trim_1p.gz",outfolder=config["outfolder"]),
    R2=expand("{outfolder}/{{ID1}}/{{ID2}}_R2_trim_2p.gz",outfolder=config["outfolder"])
  params:
    name="pdip_hostalign_{ID1}_,_{ID2}",
    hostreference=config["hostreference"]
  threads: 12
  log:
    expand("{logfolder}/{{ID1}}/{{ID2}}_pdip_hostalign.log",logfolder=config["logfolder"])
  output:
    done=expand("{logfolder}/{{ID1}}/{{ID2}}_alignrate.done",logfolder=config["logfolder"])
  shell:
    """
    touch {output.done}
    set +o nounset
    set +o pipefail
     {config[envir]}
bowtie2 -x {params.hostreference}  --threads {threads} -k 1 --very-sensitive-local --skip 1000 --omit-sec-seq --upto 1000000 -1 {input.R1} -2 {input.R2} |grep "% overall alignment rate" 2>&1 |tee {log} 

#Calculating alignment statistics
    """


rule pdip_hostsubtract:
  input:
    R1=expand("{outfolder}/{{ID1}}/{{ID2}}_R1_trim_1p.gz",outfolder=config["outfolder"]),
    R2=expand("{outfolder}/{{ID1}}/{{ID2}}_R2_trim_2p.gz",outfolder=config["outfolder"])
  params:
    name="pdip_hostsubtract_{ID1}_,_{ID2}",
    hostreference=config["hostreference"]
  threads: 12
  log:
    expand("{logfolder}/{{ID1}}/{{ID2}}_pdip_hostsubtract.log",logfolder=config["logfolder"])
  output:
    done=temp(expand("{logfolder}/{{ID1}}/{{ID2}}_hostsubtrat.done",logfolder=config["logfolder"])),
    R1=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_hostsubtracted_R1.gz",outfolder=config["outfolder"])),
    R2=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_hostsubtracted_R2.gz",outfolder=config["outfolder"]))
  shell:
    """
    touch {output.done}
    touch {output.R1}
    touch {output.R2}
    set +o nounset
    set +o pipefail
    {config[envir]}

    UNALIGNED={output.R1}
    UNALIGNED=${{UNALIGNED/_R1.gz/}}
    bowtie2 -x {params.hostreference}  --threads {threads} -k 1 --un-conc-gz ${{UNALIGNED}}% -1 {input.R1} -2 {input.R2} -S /dev/null > {log} 2>&1
    #rename files needed for trinity
    mv ${{UNALIGNED}}1 ${{UNALIGNED}}_R1.gz
    mv ${{UNALIGNED}}2 ${{UNALIGNED}}_R2.gz
    """

rule pdip_assemble:
  input:
    R1=expand("{outfolder}/{{ID1}}/{{ID2}}_hostsubtracted_R1.gz",outfolder=config["outfolder"]),
    R2=expand("{outfolder}/{{ID1}}/{{ID2}}_hostsubtracted_R2.gz",outfolder=config["outfolder"])
  params:
    name="pdip_assemble_{ID1}_,_{ID2}"
  resources:
    mem_mb=100000
  threads: 20
  log:
    expand("{logfolder}/{{ID1}}/{{ID2}}_pdip_assemble.log",logfolder=config["logfolder"])
  output:
    done=temp(expand("{logfolder}/{{ID1}}/{{ID2}}_trinity.done",logfolder=config["logfolder"])),
    R1=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_hostsubtracted_R1.fq",outfolder=config["outfolder"])),
    R2=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_hostsubtracted_R2.fq",outfolder=config["outfolder"])),
    FA=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_trinity/both.fa",outfolder=config["outfolder"])),
    KMER=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_trinity/jellyfish.kmers.fa",outfolder=config["outfolder"])),
    contigs=expand("{outfolder}/{{ID1}}/{{ID2}}_trinity/inchworm.K25.L25.DS.fa",outfolder=config["outfolder"])
  shell:
    """
    touch {output.done}
    touch {output.FA}
    touch {output.KMER}
    mkdir -p `dirname {output.contigs}`
    touch {output.contigs}
    mkfifo {output.R1}
    mkfifo {output.R2}
    set +o nounset
    set +o pipefail
    set -x
    {config[envir]}
    
    OUTPUT=`dirname {output.contigs}`
    #pipe fastq files and add /1 and /2 at the end of the readname if not present
    zcat {input.R1} | awk '{{if(NR%4==1) print gensub(/(\/1)?$/,"\/1","g"); else if(NR%4==3) print "+"; else print $0}}'> {output.R1} &
    zcat {input.R2} | awk '{{if(NR%4==1) print gensub(/(\/2)?$/,"\/2","g"); else if(NR%4==3) print "+"; else print $0}}'> {output.R2} &

    #Problem with reads not ending with /1 and /2 seqtk a component of Trinity fails 
    #Error, not recognizing read name formatting: [HS4_10034:1:2101:5750:30926#8]

    # 18-11-06 added --min_kmer_cov 2 to remove all kmers occuringonly once from the analysis -> speed up but also loss in sensitivity -> compare
    Trinity --min_contig_length 300 --left {output.R1} --right {output.R2} --min_kmer_cov 2 --no_normalize_reads --normalize_max_read_cov 1000 --CPU {threads} --no_run_chrysalis --seqType fq --output ${{OUTPUT}}  --max_memory 100G  --inchworm_cpu {threads} > {log} 2>&1  
    """

rule pdip_select_contigs:
#cutoff based on contig size to blast only larger ones
  input:
    contigs=expand("{outfolder}/{{ID1}}/{{ID2}}_trinity/inchworm.K25.L25.DS.fa",outfolder=config["outfolder"])
  params:
    name="pdip_filtercontigs_{ID1}_,_{ID2}",
    scriptfolder=config["scriptfolder"],
    mincontiglength=config["mincontiglength"]
  resources:
    mem_mb=100000
  threads: 1
  log:
    expand("{logfolder}/{{ID1}}/{{ID2}}_pdip_filtercontigs.log",logfolder=config["logfolder"])
  output:
    done=temp(expand("{logfolder}/{{ID1}}/{{ID2}}_contigsizefilter.done",logfolder=config["logfolder"])),
    contigs=expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.fa",outfolder=config["outfolder"],length=config["mincontiglength"])
  shell:
    """
    set +o nounset
    touch {output.done}
    touch {output.contigs} 
    {config[envir]}
    python {params.scriptfolder}/filter_contigs.py --outfile {output.contigs} --infile {input.contigs} --length {params.mincontiglength}
    """

rule pdip_align_contigs:
  input:
    contigs=expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.fa",outfolder=config["outfolder"],length=config["mincontiglength"]),
    R1=expand("{outfolder}/{{ID1}}/{{ID2}}_hostsubtracted_R1.gz",outfolder=config["outfolder"]),
    R2=expand("{outfolder}/{{ID1}}/{{ID2}}_hostsubtracted_R2.gz",outfolder=config["outfolder"])
  params:
    name="pdip_aligncontigs_{ID1}_,_{ID2}"
  resources:
    mem_mb=100000
  threads: 16
  log:
    expand("{logfolder}/{{ID1}}/{{ID2}}_pdip_aligncontigs.log",logfolder=config["logfolder"])
  output:
    done=temp(expand("{logfolder}/{{ID1}}/{{ID2}}_aligncontigs.done",logfolder=config["logfolder"])),
    contigindex=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}",outfolder=config["outfolder"],length=config["mincontiglength"])),
    bam=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.bam",outfolder=config["outfolder"],length=config["mincontiglength"])),
    srtbam=expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}_srt.bam",outfolder=config["outfolder"],length=config["mincontiglength"]),
    bt1=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.1.bt2",outfolder=config["outfolder"],length=config["mincontiglength"])),
    bt2=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.2.bt2",outfolder=config["outfolder"],length=config["mincontiglength"])),
    bt3=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.3.bt2",outfolder=config["outfolder"],length=config["mincontiglength"])),
    bt4=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.4.bt2",outfolder=config["outfolder"],length=config["mincontiglength"])),
    rev1=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.rev.1.bt2",outfolder=config["outfolder"],length=config["mincontiglength"])),
    rev2=temp(expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.rev.2.bt2",outfolder=config["outfolder"],length=config["mincontiglength"]))
  shell:
    """
    set +o nounset
    touch {output.done}
    touch {output.bam}
    touch {output.srtbam}
    touch {output.contigindex}
    touch {output.bt1}
    touch {output.bt2}
    touch {output.bt3}
    touch {output.bt4}
    touch {output.rev1}
    touch {output.rev2}

    {config[envir]}
 
    #build index for contigs
    bowtie2-build {input.contigs} {output.contigindex} > {log}  2>&1

    #align reads to contigs

    bowtie2 -x {output.contigindex} --threads {threads} -1 {input.R1} -2 {input.R2} --very-sensitive --omit-sec-seq | samtools view -bS - > {output.bam} 2> {log}

    #sort bam file
    samtools sort {output.bam} > {output.srtbam}
    samtools index {output.srtbam} 
    """

                 
rule pdip_blast:
  input:
    contigs=expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}.fa",outfolder=config["outfolder"],length=config["mincontiglength"])
  params:
    name="pdip_blast_{ID1}_,_{ID2}",
    blastdb=config["blastdb"]
  resources:
    mem_mb=100000
  threads: 16
  log:
    expand("{logfolder}/{{ID1}}/{{ID2}}_pdip_blast.log",logfolder=config["logfolder"])
  output:
    done=temp(expand("{logfolder}/{{ID1}}/{{ID2}}_blast.done",logfolder=config["logfolder"])),
    hits=expand("{outfolder}/{{ID1}}/{{ID2}}_contighits_{length}.gz",outfolder=config["outfolder"],length=config["mincontiglength"])
  shell:
    """
    set +o nounset
    export BLASTDB={params.blastdb}
    touch {output.done}
    touch {output.hits}
    {config[envir]}
    blastn -task blastn -evalue 10E-2 -word_size 11 -num_threads {threads} -query {input.contigs}  -db nt \
-outfmt '6 qseqid sseqid sacc saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames sblastnames sskingdoms stitle qcovs qcovhsp' |gzip > {output.hits}
    """


rule pdip_filter_blast:
  input:
    srtbam=expand("{outfolder}/{{ID1}}/{{ID2}}_contigs_{length}_srt.bam",outfolder=config["outfolder"],length=config["mincontiglength"]),
    hits=expand("{outfolder}/{{ID1}}/{{ID2}}_contighits_{length}.gz",outfolder=config["outfolder"],length=config["mincontiglength"]),
    inputreads=expand("{outfolder}/{{ID1}}/{{ID2}}_inputreads",outfolder=config["outfolder"])
  params:
    name="pdip_filter_blast_{ID1}_,_{ID2}",
    scriptsfolder=config["scriptfolder"],
    contaminants=config["contaminants"]
  resources:
    mem_mb=100000
  threads: 10
  log:
    expand("{logfolder}/{{ID1}}/{{ID2}}_pdip_filter_blast.log",logfolder=config["logfolder"])
  output:
    done=expand("{logfolder}/{{ID1}}/{{ID2}}.done",logfolder=config["logfolder"]),
    filterhits=expand("{outfolder}/{{ID1}}/{{ID2}}_contighits_{length}_filtered.csv",outfolder=config["outfolder"],length=config["mincontiglength"])
  shell:
    """
    set +o nounset
    touch {output.done}
    touch {output.filterhits}
    {config[envir]}
    N=`cat {input.inputreads}`
    Rscript --vanilla --verbose {params.scriptsfolder}/blast_filter.R \
    --blast {input.hits} \
    --out {output.filterhits} \
    --threads {threads} \
    --sequencebam {input.srtbam} \
    --contaminants {params.contaminants} \
    --inputreads ${{N}} > {log}  2>&1
    find `dirname {output.filterhits}` -type d |xargs chmod g+rx
    find `dirname {output.filterhits}` -type f |xargs chmod g+r
    """    

  
#rule alternative counting ????
#count reads aligning to full blast hit sequence

#rule format_results:
