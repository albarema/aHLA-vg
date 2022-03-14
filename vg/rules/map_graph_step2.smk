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

## --------------------------------------------------------------------------------
rule all:
    input:
       expand('/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}-{map}.bam',
              sample='HG00733',
              fullname=FULLS,
              genome=REFV,
              vcf=VCFV,
              map=MAPPER, #map
              minq=5
       ),
## --------------------------------------------------------------------------------
# map reads to the graph using giraffe

rule vg_downsample_graph:
    """
    Sample the the most common local haplotypes present in the input VCF.
    Mapping to a graph with lots of low frequency haplotypes reduces mapping accuracy.
    https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand#sampling-local-haplotypes
    """
    input:
        gbz="vg/graphs/{vcf}/GRCh38_no_alts-1kGP-N16.giraffe.gbz",
        min="vg/graphs/{vcf}/{genome}-{vcf}.k{k}.w{w}.N{n}.min",
        dst="vg/graphs/{vcf}/{genome}-{vcf}.dist",
    output:
        gbz="vg/graphs/{vcf}/GRCh38_no_alts-1kGP-N{n}.sampled.{N}.giraffe.gbz",
        min="vg/graphs/{vcf}/{genome}-{vcf}.k{k}.w{w}.N{n}.sampled.{N}.min",
        dst="vg/graphs/{vcf}/{genome}-{vcf}.sampled.{N}.dist",
    log:
        log="vg/logs/{genome}-{vcf}.k{k}.w{w}.N{n}.sampled.{N}.log",
    benchmark:
        "vg/benchmarks/{genome}-{vcf}.k{k}.w{w}.N{n}.sampled.{N}.tsv"
    shell:
        # TODO ??
        # vg gbwt --xg-name graph.xg --output sampled.gbwt --local-haplotypes haplotypes.gbwt
        "vg gbwt"
        " --gbz-input {input.gbz}"
        " --output vg/graphs/sampled.gbwt"
        " --local-haplotypes vg/graphs/haplotypes.gbwt"
        " --num-paths {wildcards.N}"

rule map_gaffe:
    input:
        r1=config['simulpath'] + "{sample}/{fullname}.adRm.fastq.gz",
        min='vg/graphs/{vcf}/{genome}-{vcf}.k{k}.w{w}.N{n}.min',
        xg='vg/graphs/{vcf}/{genome}-{vcf}.xg',
        dist='vg/graphs/{vcf}/{genome}-{vcf}.dist',
        gbwt='vg/graphs/{vcf}/{genome}-{vcf}-N{n}.gbwt',
        gg='vg/graphs/{vcf}/{genome}-{vcf}-N{n}.gg',
        gbz='vg/graphs/{vcf}/GRCh38_no_alts-1kGP-N16.giraffe.gbz'
    output: '/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}-gaffe{k}k{w}w{n}N.bam'
    benchmark: 'vg/benchmarks/{sample}/testgbz-{fullname}.{genome}.{vcf}.gaffe{k}k{w}w{n}N.benchmark.txt'
    log: 'vg/logs/{sample}/testgbz-{fullname}-{genome}-{vcf}-gaffe{k}k{w}w{n}N.log.txt'
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
        xg='vg/graphs/{vcf}/{genome}-{vcf}.xg',
        gcsa="vg/graphs/{vcf}/{genome}-{vcf}.gcsa",
        gcsalcp="vg/graphs/{vcf}/{genome}-{vcf}.gcsa.lcp"
    output: '/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}-map.bam'
    threads: 10
    resources:
        mem_mb=100000
    benchmark: 'vg/benchmarks/{sample}/{fullname}.{genome}.{vcf}.map.benchmark.txt'
    log: 'vg/logs/{sample}/{fullname}-{genome}-{vcf}-map.log.txt'
    shell:
        "vg map --threads {threads} "
        "--xg-name {input.xg} "
        "--gcsa-name {input.gcsa} "
        "--fastq {input.r1} "
        "--min-mem 15 "
        "--band-width 1024 "
        "--surject-to bam "
        "--log-time > {output} 2> {log}"


# rule surject: # f2ix threads
#     input:
#         gam='/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}.{map}.gam',
#         xg='vg/graphs/{vcf}/{genome}-{vcf}.xg',
#     output: 
#         bam='/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}.{map}.bam'
#     resources:
#         mem_mb=100000
#     benchmark: 'vg/benchmarks/{sample}/{fullname}-{genome}-{vcf}-{map}-2bam-sort.benchmark.txt'
#     log: 'vg/logs/{sample}/{fullname}-{genome}-{vcf}-{map}-2bam-sort.log.txt'
#     shell:
#         "( vg surject --bam-output --xg-name {input.xg} --threads {threads} {input.gam} | "
#         " samtools sort -@ {threads} -O bam -o {output.bam} ) 2> {log}"
    
rule bamidx: # f2ix threads
    input:
        bam='/projects/racimolab/data/MHC/benchmarking/alignment/vg/{sample}/tmp-{fullname}-{graph}.{map}.bam'
    output: 
        bam='/projects/racimolab/data/MHCbenchmarking/alignment/vg/{sample}/tmp-{fullname}-{graph}.{map}.bai'
    resources:
        mem_mb=50000
    benchmark: 'vg/benchmarks/{sample}/{fullname}-{graph}-{map}-2bam-sort.benchmark.txt'
    log: 'vg/logs/{sample}/{fullname}-{graph}-{map}-2bam-sort.log.txt'
    shell:
        "samtools index "
        "{input.bam} "
        "2> {log}"

# compute packed coverage from aligned reads
rule pack:
    input:
        gam='/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/tmp-{fullname}-{genome}-{vcf}.{map}.gam',
        xg='vg/graphs/{vcf}/{genome}-{vcf}.xg'
    output: '/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}.{map}.q{minq}.pack'
    threads: 16
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/{sample}/{fullname}-{genome}-{vcf}-{map}-q{minq}-pack.benchmark.txt'
    log: 'vg/logs/{sample}/{fullname}-{genome}-{vcf}-{map}-q{minq}-pack.log.txt'
    shell:
        "vg pack --xg {input.xg} "
        "--gam {input.gam} "
        "--min-mapq {wildcards.minq} "
        "--threads {threads} "
        "--packs-out {output} 2> {log}"


# call variants from the packed read coverage
rule call_novcf:
    input:
        pack='/projects/racimolab/data/MHC/vg/simulations/{sample}/{vcf}/{fullname}-{graph}.{map}.q{minq}.pack',
        xg='vg/graphs/{vcf}/{graph}.xg',
        snarls='vg/graphs/{vcf}/{graph}.snarls'
    output:
        vcf='/projects/racimolab/data/MHC/vg/simulation/{sample}/{vcf}/{fullname}-{graph}.{map}.q{minq}.call.vcf.gz',
    params:
        tmp_raw_vcf="{sample}-{dat}-{fullname}-{graph}-q{minq}_calltemp_raw.vcf"
    threads: 20
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{sample}/{vcf}/{fullname}-{graph}-{map}-q{minq}-call.benchmark.txt'
    log: 'vg/logs/{sample}/{vcf}/{fullname}-{graph}-{map}-q{minq}-call.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools sort {params.tmp_raw_vcf} | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')

# genotype variants from the packed read coverage
rule call_vcf:
    input:
        pack='/projects/racimolab/data/MHC/vg/simulations/{sample}/{vcf}/{fullname}-{genome}-{vcf}.{map}.q{minq}.pack',
        xg='vg/graphs/{vcf}/{genome}-{svs}.xg',
        snarls='vg/graphs/{vcf}/{genome}-{svs}.snarls',
        vcf='/projects/racimolab/data/MHC/vg/vcf/{vcf}.vcf.gz',
        vcftbi='/projects/racimolab/data/MHC/vg/vcf/{vcf}.vcf.gz.tbi'
    output:
        vcf='{sample}/{vcf}/{fullname}-{genome}-{vcf}.{map}.q{minq}.gt.vcf.gz',
        idx='{sample}/{vcf}/{fullname}-{genome}-{vcf}.{map}.q{minq}.gt.vcf.gz.tbi'
    params:
        tmp_raw_vcf="{sample}-{vcf}-{fullname}-{genome}-{vcf}-q{minq}_genotemp_raw.vcf"
    threads: 20
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{sample}/{vcf}/{fullname}-{genome}-{vcf}-{map}-q{minq}-genotype.benchmark.txt'
    log: 'vg/logs/{sample}/{vcf}/{fullname}-{genome}-{vcf}-{map}-q{minq}-genotype.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} -v {input.vcf} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools view -e 'GT=\"0/0\" || GT=\"./.\"' {params.tmp_raw_vcf} | bcftools sort | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')

# use an index graph
rule map_gaffe_auto:
    input:
        r1=config['simulpath'] + "{sample}/{fullname}.adRm.fastq.gz",
        min='vg/graphs/{vcf}/autoindex/{genome}-{vcf}.allchr_filtered.autoindex.min',
        dist='vg/graphs/{vcf}/autoindex/{genome}-{vcf}.allchr_filtered.autoindex.dist',
        gbz='vg/graphs/{vcf}/autoindex/{genome}-{vcf}.allchr_filtered.autoindex.giraffe.gbz',
    output: '/projects/racimolab/data/MHC/benchmarking/alignment/vg/simulations/{sample}/{fullname}-{genome}-{vcf}.gaffe.autoindex.gam'
    benchmark: 'vg/benchmarks/{sample}/{fullname}.{genome}.{vcf}.gaffe.autoindex.benchmark.txt'
    log: 'vg/logs/{sample}/{fullname}-{genome}-{vcf}-gaffe.autoindex.log.txt'
    threads: 8
    shell:
        "vg giraffe --progress "
        "--threads {threads} "
        "--minimizer-name {input.min} "
        "--dist-name {input.dist} "
        "--gbz-name {input.gbz} "
        "--fastq-in {input.r1} "
        " > {output} "
        "2> {log}"