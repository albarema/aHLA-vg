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
    MAPPER = 'gaffe{}k{}w{}N.aug'.format(config['mink'], config['minw'], config['covern'])

# choose vcf
VCFV=config['vcf_version']
VCFPATH=config['vcf'][VCFV]['minMAF']

# choose reference
REFV=config['ref_version']

GRAPH=REFV+ '-'+VCFV


## --------------------------------------------------------------------------------
rule all:
    input:
       expand('/projects/racimolab/data/MHC/vg/simulations/{sample}/{dat}/{fullname}-{graph}-{map}.bam',
              sample=SAMPLE,
              dat=VCFV,
              fullname=FULLS,
              graph=GRAPH,
              map=MAPPER,
              minq=5
       )
## --------------------------------------------------------------------------------

# map reads to the graph using giraffe
rule map_gaffe:
    input:
        r1=config['simulpath'] + "{sample}/{fullname}.adRm.fastq.gz",
        xg='vg/graphs/{dat}/{genome}-{vcf}.aug.xg',
        min='vg/graphs/{dat}/{genome}-{vcf}.k{k}.w{w}.N{n}.aug.min',
        dist='vg/graphs/{dat}/{genome}-{vcf}.aug.dist',
        gbwt='vg/graphs/{dat}/{genome}-{vcf}-N{n}.aug.gbwt',
    output: '/projects/racimolab/data/MHC/vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}-gaffe{k}k{w}w{n}N.aug.bam'
    threads: 6
    resources:
        mem_mb=70000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}.{genome}.{vcf}.gaffe{k}k{w}w{n}N.aug.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{genome}-{vcf}-gaffe{k}k{w}w{n}N.aug.log.txt'
    shell:
        "vg giraffe --progress "
        "--threads {threads} "
        "--minimizer-name {input.min} "
        "--dist-name {input.dist} "
        "--gbwt-name {input.gbwt} "
        "--xg-name {input.xg} "
        "--fastq-in {input.r1} "
        "--output-format BAM > {output} "
        "2> {log}"


# map reads to the graph using mpmap in single-path mode
rule map_map:
    input:
        r1=config['simulpath'] + "{sample}/{fullname}.adRm.fastq.gz",
        xg='vg/graphs/{dat}/{genome}-{vcf}.aug.xg',
        gcsa="vg/graphs/{dat}/{genome}-{vcf}.aug.gcsa",
        gcsalcp="vg/graphs/{dat}/{genome}-{vcf}.aug.gcsa.lcp"
    output: '/projects/racimolab/data/MHC/vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}-map.aug.bam'
    threads: 8
    resources:
        mem_mb=100000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}.{genome}.{vcf}.map.aug.benchmark.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{genome}-{vcf}-map.aug.log.txt'
    shell:
        "vg map --threads {threads} "
        "--xg-name {input.xg} "
        "--min-mem 15 "
        "--band-width 1024 "
        "--gcsa-name {input.gcsa} "
        "--fastq {input.r1} "
        "--surject-to bam "
        "--log-time > {output} 2> {log}"

# compute packed coverage from aligned reads
rule pack:
    input:
        gam='/projects/racimolab/data/MHC/vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.gam',
        xg='vg/graphs/{dat}/{graph}.xg'
    output: '/projects/racimolab/data/MHC/vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.pack'
    threads: 16
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}-{graph}-{map}-q{minq}-pack.benchmark.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{graph}-{map}-q{minq}-pack.log.txt'
    shell:
        "vg pack --xg {input.xg} "
        "--gam {input.gam} "
        "--min-mapq {wildcards.minq} "
        "--threads {threads} "
        "--packs-out {output} 2> {log}"


# call variants from the packed read coverage
rule call_novcf:
    input:
        pack='/projects/racimolab/data/MHC/vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.pack',
        xg='vg/graphs/{dat}/{graph}.xg',
        snarls='vg/graphs/{dat}/{graph}.snarls'
    output:
        vcf='/projects/racimolab/data/MHC/vg/simulation/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.call.vcf.gz',
    params:
        tmp_raw_vcf="{sample}-{dat}-{fullname}-{graph}-q{minq}_calltemp_raw.vcf"
    threads: 20
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}-{graph}-{map}-q{minq}-call.benchmark.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{graph}-{map}-q{minq}-call.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools sort {params.tmp_raw_vcf} | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')

# genotype variants from the packed read coverage
rule call_vcf:
    input:
        pack='/projects/racimolab/data/MHC/vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}.{map}.q{minq}.pack',
        xg='vg/graphs/{dat}/{genome}-{svs}.xg',
        snarls='vg/graphs/{dat}/{genome}-{svs}.snarls',
        vcf='/projects/racimolab/data/MHC/vg/vcf/{vcf}.vcf.gz',
        vcftbi='/projects/racimolab/data/MHC/vg/vcf/{vcf}.vcf.gz.tbi'
    output:
        vcf='{sample}/{dat}/{fullname}-{genome}-{vcf}.{map}.q{minq}.gt.vcf.gz',
        idx='{sample}/{dat}/{fullname}-{genome}-{vcf}.{map}.q{minq}.gt.vcf.gz.tbi'
    params:
        tmp_raw_vcf="{sample}-{dat}-{fullname}-{genome}-{vcf}-q{minq}_genotemp_raw.vcf"
    threads: 20
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}-{genome}-{vcf}-{map}-q{minq}-genotype.benchmark.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{genome}-{vcf}-{map}-q{minq}-genotype.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} -v {input.vcf} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools view -e 'GT=\"0/0\" || GT=\"./.\"' {params.tmp_raw_vcf} | bcftools sort | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')

