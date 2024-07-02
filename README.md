# P-DiP

Pathogen Discovery Pipeline

# clone workflow into working directory
```
git clone path/to/workdir
cd path/to/workdir
```
# edit config as needed
```
vi config.yaml
```
change at least out and log folder

# install dependencies into an isolated conda environment
install miniconda from https://conda.io/miniconda.html
```
conda env create -n pdip --file pdip.yaml 
```
2024-07-02 miniconda environment seems not to build (user issue was raised)
build the environment with micromamba 
(https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
```
micromamba env create -n pdip --file pdip.yaml 
```

# get databases
## blast
```
source activate pdip
cd database
mkdir blast_`date +%Y-%m-%d`
cd blast_`date +%Y-%m-%d`
perl [conda bin folder,pdip/bin]/update_blastdb.pl --decompress nt nr
```

# activate environment
```
source activate pdip
```

# execute workflow
```
snakemake -n -s Snakefile.py --configfile=config.yaml
```

# example run for PBS cluster
```
snakemake  -p --jobs 100 \
--configfile=config.yaml \
-s Snakefile.py \
--cluster-config cluster.json \
--cluster "qsub -j eo -e  cluster.log -N {params.name} -V -l nodes=1:ppn={threads},mem={cluster.mem},walltime={cluster.walltime}" \
--keep-going --latency-wait 100 --rerun-incomplete
```
