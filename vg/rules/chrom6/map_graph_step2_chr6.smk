##### STEP 2  - mapping #####
# run from hla
configfile: "config.yaml"
# in config file set simulpath
import pandas as pd

## --------------------------------------------------------------------------------

wildcard_constraints:
    vcf="[^-+\.?$]+",
    chrom="[^-+\.?$]+",
    ext="[^.]+",
    fullname="[^-]+"

# simulations wrap wildcards
SAMPLE='HG02080.paternal' # 'HG02080.paternal'
# A. map all coverage for giraffe
FULLSPATH = config['simulpath'] + SAMPLE +'/' + SAMPLE + '.filelist'
# B. map only 1c for vg map
#FULLSPATH = config['simulpath'] + SAMPLE +'/' + SAMPLE + '.1c.filelist'

FULLS = pd.read_table(FULLSPATH,names=['fname'])['fname'].tolist()


# graph wrap wildcards
MAPPER=config['mapper']
if MAPPER == 'gaffe':
    MAPPER = 'gaffe{}k{}w{}N'.format(config['mink'], config['minw'], config['covern'])

# choose vcf
VCFV='1kGP' #-1kGP
VCFPATH=config['vcf'][VCFV]['minMAF']
# choose reference
REFV=config['ref_version']
REFPATH=config['ref'][REFV]

GRAPH=REFV+ '-'+VCFV #+'-6'

## --------------------------------------------------------------------------------

rule all:
    input:
        expand('vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.pack',
            sample=SAMPLE,
            dat=VCFV, #1kGP
            fullname=FULLS,
            graph=GRAPH,
            map=MAPPER,
            minq=5
        ),
        expand('vg/simulations/bam/{sample}/{dat}/{fullname}-{graph}.{map}.sorted.bam',
            sample=SAMPLE,
            dat=VCFV, #1kGP
            fullname=FULLS,
            graph=GRAPH,
            map=MAPPER,
            minq=5
)

## --------------------------------------------------------------------------------
#
# Map reads from a sample and call variants

# map reads to the graph using mpmap in single-path mode
rule map_map:
    input:
        r1=config['simulpath'] + "{sample}/{fullname}.adRm.fastq.gz",
        xg='vg/graphs/{dat}/{genome}-{vcf}.xg',
        gcsa="vg/graphs/{dat}/{genome}-{vcf}.gcsa",
        gcsalcp="vg/graphs/{dat}/{genome}-{vcf}.gcsa.lcp"
    output: 'vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}.map.gam'
    threads: 16
    resources:
        mem_mb=100000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}.{genome}.{vcf}.map.benchmark.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{genome}-{vcf}-map.log.txt'
    shell:
        "vg map "
        "--threads {threads} "
        "--xg-name {input.xg} "
        "--gcsa-name {input.gcsa} "
        " --exclude-unaligned"
        "--fastq {input.r1} "
        " --band-width 1024"
        "--log-time "
        "1> {output} 2> {log}"
#-f {input.r2}

# map reads to the graph using giraffe
rule map_gaffe:
    input:
        r1=config['simulpath'] + "{sample}/{fullname}.adRm.fastq.gz",
        xg='vg/graphs/{dat}/{genome}-{vcf}.xg',
        min='vg/graphs/{dat}/{genome}-{vcf}.k{k}.w{w}.N{n}.min',
        dist='vg/graphs/{dat}/{genome}-{vcf}.dist',
        gbwt='vg/graphs/{dat}/{genome}-{vcf}-N{n}.gbwt',
    output: 'vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}.gaffe{k}k{w}w{n}N.gam'
    threads: 16
    resources:
        mem_mb=100000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}.{genome}.{vcf}.gaffe{k}k{w}w{n}N.benchmark.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{genome}-{vcf}-gaffe{k}k{w}w{n}N.log.txt'
    shell:
        "vg giraffe "
        "--progress "
        "--threads {threads} "
        "--minimizer-name {input.min} "
        "--dist-name {input.dist} "
        "--gbwt-name {input.gbwt} "
        "--xg-name {input.xg} "
        "--fastq-in {input.r1} "
        "1> {output} 2> {log}"


rule surject:
    """
    Surject the alignments back into the reference space, yielding a regular BAM file
    """
    input:
        xg='vg/graphs/{dat}/{graph}.xg',
        gam='vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.gam',
    output:
        bam="vg/simulations/bam/{sample}/{dat}/{fullname}-{graph}.{map}.sorted.bam",
        csi="vg/simulations/bam/{sample}/{dat}/{fullname}-{graph}.{map}.sorted.bam.csi",
    log:
        log="data/samples/bam/{sample}/{dat}/{fullname}-{graph}.{map}.sorted.bam.log",
    threads: 8
    shell:
        "( vg surject --bam-output --xg-name {input.xg} --threads {threads} {input.gam} | "
        "  samtools sort --write-index -@ {threads} -O bam -o {output.bam} ) 2> {log}"

# compute packed coverage from aligned reads
rule pack:
    input:
        gam='vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.gam',
        xg='vg/graphs/{dat}/{graph}.xg'
    output: 'vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.pack'
    threads: 16
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}-{graph}-{map}-q{minq}-pack.benchmark.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{graph}-{map}-q{minq}-pack.log.txt'
    shell:
        "vg pack -x {input.xg} -g {input.gam} -Q {wildcards.minq} -t {threads} -o {output} 2> {log}"

# call variants from the packed read coverage
rule call_novcf:
    input:
        pack='vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.pack',
        xg='vg/graphs/{dat}/{graph}.xg',
        snarls='vg/graphs/{dat}/{graph}.snarls'
    output:
        vcf='vg/simulations/vcf/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.call.vcf.gz',
    params:
        tmp_raw_vcf="vg/simulations/vcf/{sample}-{dat}-{fullname}-{graph}-q{minq}_calltemp_raw.vcf",
        bcftools=config['bcftools']
    threads: 20
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}-{graph}-{map}-q{minq}-call.benchmark.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{graph}-{map}-q{minq}-call.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} --ploidy 1 -s {wildcards.sample} --snarls {input.snarls} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("{params.bcftools}/bcftools sort {params.tmp_raw_vcf} | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')

# genotype variants from the packed read coverage
rule call_vcf:
    input:
        pack='vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}.{map}.q{minq}.pack',
        xg='vg/graphs/{dat}/{genome}-{vcf}.xg',
        snarls='vg/graphs/{dat}/{genome}-{vcf}.snarls',
        vcf=config['vcf'][VCFV]['minMAF'],
        vcftbi='{vcf}.vcf.gz.tbi'
    output:
        vcf='vg/simulations/vcf/{sample}/{dat}/{fullname}-{genome}-{vcf}.{map}.q{minq}.gt.vcf.gz',
        idx='vg/simulations/vcf/{sample}/{dat}/{fullname}-{genome}-{vcf}.{map}.q{minq}.gt.vcf.gz.tbi'
    params:
        tmp_raw_vcf="vcf/{sample}-{dat}-{fullname}-{genome}-{vcf}-q{minq}_genotemp_raw.vcf"
    threads: 20
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}-{genome}-{vcf}-{map}-q{minq}-genotype.benchmark.txt'
    log: 'vg/logs/{sample}/{dat}/{fullname}-{genome}-{vcf}-{map}-q{minq}-genotype.log.txt'
    run:
        shell("vg call -k {input.pack} --ploidy 1 -t {threads} -s {wildcards.sample} --snarls {input.snarls} -v {input.vcf} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools view -e 'GT=\"0/0\" || GT=\"./.\"' {params.tmp_raw_vcf} | bcftools sort | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')

