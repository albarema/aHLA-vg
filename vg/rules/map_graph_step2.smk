##### STEP 2  - mapping #####
# in config file set simulpath
import pandas as pd

## --------------------------------------------------------------------------------

wildcard_constraints:
    vcf="[^-+\.?$]+",
    chrom="[^-+\.?$]+",
    ext="[^.]+",
    fullname="[^-]+"

# simulations wrap wildcards
SAMPLE='HG02080.maternal' # 'HG02080.paternal'
FULLSPATH = config['simulpath'] + SAMPLE +'/' + SAMPLE + '.filelist'
FULLS = pd.read_table(FULLSPATH,names=['fname'])['fname'].tolist()


# graph wrap wildcards
MAPPER=config['mapper']
if MAPPER == 'gaffe':
    MAPPER = 'gaffe{}k{}w{}N'.format(config['mink'], config['minw'], config['covern'])

# choose vcf
VCFV=config['vcf_version']
VCFPATH=config['vcf'][VCFV]['minMAF']
# choose reference
REFV=config['ref_version']
REFPATH=config['ref'][REFV]

GRAPH=REFV+ '-'+VCFV

## --------------------------------------------------------------------------------

rule all:
    input:
        expand('vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.pack',
               sample=SAMPLE,
               dat=VCFV,
               fullname=FULLS,
               graph=GRAPH,
               map=MAPPER,
               minq=5
        )

## --------------------------------------------------------------------------------
#
# Map reads from a sample and call variants
#

# Eventually split the reads into chunk to map in parallel
if config['nb_split_reads'] > 0:
    # split fastq into chunks
    checkpoint split_reads_1:
        input: config['simulpath'] + "{sample}/{fullname}.adRm.fastq.gz"
        output:
            dir=directory('vg/simulations/{sample}/{fullname}/read_chunks')
        threads: 8
        benchmark: 'vg/benchmarks/mapping/{sample}/{fullname}-1-splitreads.benchmark.txt'
        params:
            NLINES=config['nb_split_reads'] * 4
        shell:
            "gzip -cd {input} | split -d -l {params.NLINES} --filter='pigz -p {threads} > ${{FILE}}.fastq.gz' - \"{output.dir}/{wildcards.fullname}.part\""
    # set the input paths for the mapping rules
    read1_in = 'vg/simulations/{sample}/{fullname}/read_chunks/{fullname}.part{part}.fastq.gz'
    map_lab = '{sample}/{dat}/{fullname}.part{part}'
    map_out = 'vg/simulations/{sample}/{dat}/{fullname}/read_chunks/{fullname}-{genome}-{vcf}.part{part}'
else:
    # set the input paths for the mapping rules (no chunks)
    read1_in = config['simulpath'] + "{sample}/{fullname}.fastq.gz"
    map_lab = '{sample}/{dat}/{fullname}'
    map_out = 'vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}'

#shell("gzip -cd {input} | split -d -l {NLINES} --filter='pigz -p {threads} > ${{FILE}}.fastq.gz' - \"{output.dir}/{wildcards.sample}_1.part\"")

# map reads to the graph using mpmap in single-path mode
rule map_map:
    input:
        r1=read1_in,
        xg='vg/graphs/{dat}/{genome}-{vcf}.xg',
        gcsa="vg/graphs/{dat}/{genome}-{vcf}.gcsa",
        gcsalcp="vg/graphs/{dat}/{genome}-{vcf}.gcsa.lcp"
    output: map_out + '.map.gam'
    threads: 16
    resources:
        mem_mb=100000
    benchmark: 'vg/benchmarks/' + map_lab + '.{genome}.{vcf}.map.benchmark.txt'
    log: 'vg/logs/' + map_lab + '-{genome}-{vcf}-map.log.txt'
    shell:
        "vg map -t {threads} -x {input.xg} -g {input.gcsa} -f {input.r1} --log-time > {output} 2> {log}"
#-f {input.r2}

# map reads to the graph using giraffe
rule map_gaffe:
    input:
        r1=read1_in,
        #r2=read2_in,
        xg='vg/graphs/{dat}/{genome}-{vcf}.xg',
        min='vg/graphs/{dat}/{genome}-{vcf}.k{k}.w{w}.N{n}.min',
        dist='vg/graphs/{dat}/{genome}-{vcf}.dist',
        gbwt='vg/graphs/{dat}/{genome}-{vcf}.N{n}.gbwt',
    output: map_out + '.gaffe{k}k{w}w{n}N.gam'
    threads: 16
    resources:
        mem_mb=100000
    benchmark: 'vg/benchmarks/' + map_lab + '.{genome}.{vcf}.gaffe{k}k{w}w{n}N.benchmark.txt'
    log: 'vg/logs/' + map_lab + '-{genome}-{vcf}-gaffe{k}k{w}w{n}N.log.txt'
    shell:
        "vg giraffe -p -t {threads} -m {input.min} -d {input.dist} --gbwt-name {input.gbwt} -x {input.xg} -f {input.r1} > {output} 2> {log}"
#-N {wildcards.sample}
#-f {input.r2}

# If the reads were split, they need to be merged back after alignment
if config['nb_split_reads'] > 0:
    # set variables back to not chunks
    map_lab = '{sample}/{dat}/{fullname}'
    map_out = 'vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}'
    # merge aligned reads
    def aggregate_reads(wildcards):
        checkpoint_output = checkpoints.split_reads_1.get(**wildcards).output[0]
        return expand("{sample}/read_chunks/{fullname}-{graph}.part{part}.{map}.gam",
                      sample=wildcards.sample, fullname=wildcards.fullname, graph=wildcards.graph, map=wildcards.map,
                      part=glob_wildcards(os.path.join(checkpoint_output, wildcards.fullname + ".part{part}.fastq.gz")).part)
    rule merge_gam:
        input: aggregate_reads
        output: 'vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.gam'
        benchmark: 'vg/benchmarks/{sample}/{dat}/{fullname}-{graph}-{map}-mergegam.benchmark.txt'
        shell:
            "cat {input} > {output}"

#'vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}
#{sample}/{dat}/{fullname}

# compute packed coverage from aligned reads
rule pack:
    input:
        gam='vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.gam',
        xg='vg/graphs/{dat}/{graph}.xg'
    output: 'vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.pack'
    threads: 16
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/' + map_lab + '-{graph}-{map}-q{minq}-pack.benchmark.txt'
    log: 'vg/logs/' + map_lab + '-{graph}-{map}-q{minq}-pack.log.txt'
    shell:
        "vg pack -x {input.xg} -g {input.gam} -Q {wildcards.minq} -t {threads} -o {output} 2> {log}"


# call variants from the packed read coverage
rule call_novcf:
    input:
        pack='vg/simulations/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.pack',
        xg='vg/graphs/{dat}/{graph}.xg',
        snarls='vg/graphs/{dat}/{graph}.snarls'
    output:
        vcf='vg/simulation/{sample}/{dat}/{fullname}-{graph}.{map}.q{minq}.call.vcf.gz',
    params:
        tmp_raw_vcf="{sample}-{dat}-{fullname}-{graph}-q{minq}_calltemp_raw.vcf"
    threads: 20
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/' + map_lab + '-{graph}-{map}-q{minq}-call.benchmark.txt'
    log: 'vg/logs/' + map_lab + '-{graph}-{map}-q{minq}-call.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools sort {params.tmp_raw_vcf} | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')

# genotype variants from the packed read coverage
rule call_vcf:
    input:
        pack='vg/simulations/{sample}/{dat}/{fullname}-{genome}-{vcf}.{map}.q{minq}.pack',
        xg='vg/graphs/{dat}/{genome}-{svs}.xg',
        snarls='vg/graphs/{dat}/{genome}-{svs}.snarls',
        vcf='{vcf}.vcf.gz',
        vcftbi='{vcf}.vcf.gz.tbi'
    output:
        vcf='{sample}/{dat}/{fullname}-{genome}-{vcf}.{map}.q{minq}.gt.vcf.gz',
        idx='{sample}/{dat}/{fullname}-{genome}-{vcf}.{map}.q{minq}.gt.vcf.gz.tbi'
    params:
        tmp_raw_vcf="{sample}-{dat}-{fullname}-{genome}-{vcf}-q{minq}_genotemp_raw.vcf"
    threads: 20
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/' + map_lab + '-{genome}-{vcf}-{map}-q{minq}-genotype.benchmark.txt'
    log: 'vg/logs/' + map_lab + '-{genome}-{vcf}-{map}-q{minq}-genotype.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} -v {input.vcf} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools view -e 'GT=\"0/0\" || GT=\"./.\"' {params.tmp_raw_vcf} | bcftools sort | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')

