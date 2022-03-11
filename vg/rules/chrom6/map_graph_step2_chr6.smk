##### STEP 2  - mapping #####
# run from hla
configfile: "config.yaml"
# in config file set simulpath
import pandas as pd

## --------------------------------------------------------------------------------------------------

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
    MAPPER = 'gaffe{}k{}w{}N.6.aug'.format(config['mink'], config['minw'], config['covern'])

if MAPPER == "map":
    MAPPER = "map.6.aug"

# choose vcf
VCFV=config['vcf_version']
VCFPATH=config['vcf'][VCFV]['minMAF']

# choose reference
REFV=config['ref_version']



## --------------------------------------------------------------------------------
rule all:
    input:
       expand('/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}-{map}.bam',
              sample='HG00733',
              fullname=FULLS,
              genome=REFV,
              vcf=VCFV,
              map='map.6.aug',
              minq=5
       ),
## --------------------------------------------------------------------------------
#
# Map reads from a sample and call variants

# map reads to the graph using mpmap in single-path mode
# map reads to the graph using giraffe
rule map_gaffe:
    input:
        r1=config['simulpath'] + "{sample}/{fullname}.adRm.fastq.gz",
        min='vg/graphs/{vcf}/{genome}-{vcf}-6.k{k}.w{w}.N{n}.aug.min',
        xg='vg/graphs/{vcf}/{genome}-{vcf}-6.aug.xg',
        dist='vg/graphs/{vcf}/{genome}-{vcf}-6.aug.dist',
        gbwt='vg/graphs/{vcf}/{genome}-{vcf}-chr6-N{n}.aug.gbwt',
        gg='vg/graphs/{vcf}/{genome}-{vcf}-chr6-N{n}.aug.gg',
#        gbz='vg/graphs/{vcf}/GRCh38_no_alts-1kGP-N16.giraffe.gbz'
    output: '/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}-gaffe{k}k{w}w{n}N.6.aug.bam'
    benchmark: 'vg/benchmarks/{sample}/testgbz-{fullname}.{genome}.{vcf}.gaffe{k}k{w}w{n}N.6.aug.benchmark.txt'
    log: 'vg/logs/{sample}/testgbz-{fullname}-{genome}-{vcf}-gaffe{k}k{w}w{n}N.6.aug.log.txt'
    threads: 8
    resources:
        mem_mb=70000
    shell:
        "vg giraffe --progress "
        "--threads {threads} "
        "--minimizer-name {input.min} "
        "--xg-name {input.xg} "
        "--dist-name {input.dist} "
        "--fastq-in {input.r1} "
        "--output-format BAM "
        " > {output} "
        "2> {log}"

# we can algo prodvide gbwt-name, gbwt-graph and xg graph instead of gbz
# xg='vg/graphs/{vcf}/{genome}-{vcf}.xg',
#        "--xg-name {input.xg} "
#        "--gbz-name {input.gbz} "
# map reads to the graph using mpmap in single-path mode
rule map_map:
    input:
        r1=config['simulpath'] + "{sample}/{fullname}.adRm.fastq.gz",
        xg='vg/graphs/{vcf}/{genome}-{vcf}-6.aug.xg',
        gcsa="vg/graphs/{vcf}/{genome}-{vcf}-6.aug.gcsa",
        gcsalcp="vg/graphs/{vcf}/{genome}-{vcf}-6.aug.gcsa.lcp"
    output: '/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}-map.6.aug.bam'
    threads: 10
    resources:
        mem_mb=100000
    benchmark: 'vg/benchmarks/{sample}/{fullname}.{genome}.{vcf}.map.6.aug.benchmark.txt'
    log: 'vg/logs/{sample}/{fullname}-{genome}-{vcf}-map.6.aug.log.txt'
    shell:
        "vg map --threads {threads} "
        "--xg-name {input.xg} "
        "--gcsa-name {input.gcsa} "
        "--fastq {input.r1} "
        "--min-mem 15 "
        "--band-width 1024 "
        "--surject-to bam "
        "--log-time > {output} 2> {log}"

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

