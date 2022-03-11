#
# Augment graph and call variants
#

# convert a XG index to a PG index
rule pgconvert:
    input: '{graph}.xg'
    output: '{graph}.pg'
    threads: 1
    resources:
        mem_mb=config['mem_pgconvert']
    benchmark: 'benchmarks/{graph}-pgconvert.benchmark.txt'
    log: 'logs/{graph}-pgconvert.log.txt'
    run:
        shell('vg convert {input} -p > {output} 2> {log}')

# augment a graph with aligned reads
rule augment:
    input:
        pg='{graph}.pg',
        gam='{sample}/{sample}-{graph}.{map}.gam'
    output:
        gam='{sample}/{sample}-{graph}.{map}.aug.gam',
        pg='{sample}/{sample}-{graph}.{map}.aug.pg'
    threads: config['cores_augment']
    resources:
        mem_mb=config['mem_augment']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-augment.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-augment.log.txt'             
    run:
        shell('vg augment {input.pg} {input.gam} -t {threads} -m 4 -q 5 -Q 5 -A {output.gam} > {output.pg} 2> {log}')


# prepare the snarls index for the augmented graph
rule index_snarls_aug:
    input: '{sample}/{sample}-{graph}.{map}.aug.pg'
    output: '{sample}/{sample}-{graph}.{map}.aug.snarls'
    threads: config['cores_snarls']
    resources:
        mem_mb=config['mem_snarls']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-aug-snarls.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-aug-snarls.log.txt'
    run:
        shell('vg snarls -t {threads} {input} > {output} 2> {log}')


# compute packed coverage on augmented graph
rule pack_aug:
    input:
        gam='{sample}/{sample}-{graph}.{map}.gam',
        pg='{sample}/{sample}-{graph}.{map}.aug.pg'
    output: '{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack'
    threads: config['cores_pack']
    resources:
        mem_mb=config['mem_pack']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-aug-q{minq}-pack.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-aug-q{minq}-pack.log.txt'
    run:
        shell("vg pack -x {input.pg} -g {input.gam} -Q {wildcards.minq} -t {threads} -o {output} 2> {log}")

# call variants from the packed read coverage
rule call_aug:
    input:
        pack='{sample}/{sample}-{graph}.{map}.aug.q{minq}.pack',
        pg='{sample}/{sample}-{graph}.{map}.aug.pg',
        snarls='{sample}/{sample}-{graph}.{map}.aug.snarls'
    output:
        vcf='{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.vcf.gz',
        idx='{sample}/{sample}-{graph}.{map}.aug.q{minq}.call.vcf.gz.tbi'
    params:
        tmp_raw_vcf="{sample}-{graph}-aug.q{minq}_calltemp_raw.vcf"
    threads: config['cores_call']
    resources:
        mem_mb=config['mem_call']
    benchmark: 'benchmarks/{sample}-{graph}-{map}-aug-q{minq}-call.benchmark.txt'
    log: 'logs/{sample}-{graph}-{map}-aug-q{minq}-call.log.txt'
    run:
        shell("vg call -k {input.pack} -t {threads} -s {wildcards.sample} --snarls {input.snarls} {input.xg} > {params.tmp_raw_vcf} 2> {log}")
        shell("bcftools sort {params.tmp_raw_vcf} | bgzip > {output.vcf}")
        shell("tabix -f -p vcf {output.vcf}")
        shell('rm {params.tmp_raw_vcf}')