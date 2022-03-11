##### STEP 2  - mapping #####
# in config file set simulpath
# Run from /projects/racimolab/people/gsd818/aHLA-vg 
import pandas as pd
import os
configfile: "config.yaml"

## --------------------------------------------------------------------------------

wildcard_constraints:
    vcf="[^-+\.?$]+",
    chrom="[^-+\.?$]+",
    ext="[^.]+",
    fullname="[^-]+",
    dat="[^-]+"

# simulations wrap wildcards
SAMPLE='HG00733'
FULLSPATH = config['simulpath'] + SAMPLE +'/' + SAMPLE + '.filelist'
FULLS = pd.read_table(FULLSPATH,names=['fname'])['fname'].tolist()

df=pd.read_table(config['pangenomes_list'])
PANGEN=pd.unique(df['assemblies_names'])
COVS = config['coverage']
DAM = config['damage']
RLEN=config['readlen']

# graph wrap wildcards
MAPPER=config['mapper']
if MAPPER == 'gaffe':
    MAPPER = 'gaffe{}k{}w{}N'.format(config['mink'], config['minw'], config['covern'])

# choose vcf
VCFV=config['vcf_version']
VCFPATH=config['vcf'][VCFV]['minMAF']

# choose reference
REFV=config['ref_version']

GRAPH=REFV+ '-'+VCFV


## --------------------------------------------------------------------------------
rule all:
    input:
        expand('vg/mapping-stats/{sample}/{fullname}-{graph}-{map}.bam.stats',
            sample=SAMPLE,
            dat=VCFV,
            fullname=FULLS,
            graph=GRAPH,
            map=MAPPER,
            minq=5),
        expand("/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{graph}-{map}.filter-chr6.bam",
            sample=SAMPLE,
            dat=VCFV,
            fullname=FULLS,
            graph=GRAPH,
            map=MAPPER,
            minq=5),

rule bamidx: # f2ix threads
    input:
        bam='/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{graph}-{map}.bam'
    output: 
        bam='/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{graph}-{map}.bam.bai'
    resources:
        mem_mb=50000
    benchmark: 'vg/benchmarks/{sample}/{fullname}-{graph}-{map}-index.benchmark.txt'
    log: 'vg/logs/{sample}/{fullname}-{graph}-{map}-index.log.txt'
    shell:
        "samtools index "
        "{input.bam} "
        "2> {log}"

rule stats:
    input:
        "/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{graph}-{map}.bam",
    output:
        "vg/mapping-stats/{sample}/{fullname}-{graph}-{map}.bam.stats"
    shell:
        "samtools stats {input} > {output}"

rule plots:
    input:
        "vg/mapping-stats/{sample}/{fullname}-{graph}-{map}.bam.stats",
    output:
        "vg/mapping-stats/{sample}/{fullname}-{graph}-{map}.bam.plot"
    shell:
        "samtools plot-bamstats -p {output.plot} {input} "

rule filter_chr6:
    input:
        bam="/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{graph}-{map}.bam",
        bai="/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{graph}-{map}.bam.bai",
    output:
        "/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{graph}-{map}.filter-chr6.bam"
    shell:
        "samtools view -c {input} chr6:25726063-33400644 >> {output}"
