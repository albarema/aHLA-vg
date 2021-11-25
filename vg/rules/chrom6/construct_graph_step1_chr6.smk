##### step 1 #####
# NB: already run it for 1kG, now running it for the pangenomes ones. Version without haplo info 
# try to use haplotypes with mapping file - maybe add to dat wildcard 1kGP-mapping or pangenome-mapping
configfile: "config.yaml"

CHROMS = 6
# choose prune options
PRUNEOP=config['prune_version']
OPTS=config['prune_options'][PRUNEOP]

# choose vcf 
VCFV=config['vcf_version']
# choose reference
REFV=config['ref_version']
REFPATH=config['ref'][REFV]

wildcard_constraints:
    vcf="[^-]+",
    chrom="[^-+\.?$]+",
    ext="[^.]+",

## ------------------------------------------------------------------------------------

rule construct_all_gaffe:
    input:
        expand('vg/graphs/{dat}/{genome}-{vcf}-{chrom}.{ext}', dat=VCFV, chrom=6, genome=REFV, vcf=VCFV, ext=['xg', 'trivial.snarls', 'dist']),
        expand('vg/graphs/{dat}/{genome}-{vcf}-{chrom}.k{k}.w{w}.N{n}.min',dat=VCFV, chrom=6, genome=REFV, vcf=VCFV, k=config['mink'],
               w=config['minw'], n=config['covern'])

# indexes used by the default mapper and variant caller
rule construct_all:
    input:
        expand('vg/graphs/{dat}/{genome}-{vcf}-{chrom}.{ext}', dat=VCFV, genome=REFV, vcf=VCFV, chrom=6,ext=['xg', 'snarls', 'gcsa'])


## ------------------------------------------------------------------------------------
rule construct_chr:
    """
    constructing graph
    """
    input:
        ref=config['ref'][REFV],
        vcf=config['vcf'][VCFV]['minMAF'],
        tbi=config['vcf'][VCFV]['minMAF']+ ".tbi"
    output:
        "vg/graphs/{dat}/{genome}-{vcf}-{chrom}.vg"
    log:
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}.vg_construct.log.txt"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.vg_construct.benchmark.log.txt"
    threads: 4
    shell:
        "vg construct "
        "-r {input.ref} "
        "-v {input.vcf} "
        "-C -R chr{wildcards.chrom} "
        "--node-max 32 "
        "--alt-paths "
        "--handle-sv "
        "--progress 1> {output} "
        "2> {log} "

# maybe add --flat-alts: don't chop up alternate alleles from input VCF
# -R 

rule index_xg:
    """
    Make xg index containing the alts paths. 
    """
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.vg",
    output:
        xg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.xg",
    threads: 16
    resources:
        mem_xg=60000
    log:
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}.index.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.index.log"
    shell:
        "vg index "
        "--temp-dir vg/tmp "
        "--xg-name {output.xg} "
        "--xg-alts {input} "
        "--progress 2> {log}"

rule index_snarls:
    """
    prepare the snarls index
    """
    input: 'vg/graphs/{dat}/{genome}-{vcf}-{chrom}.xg'
    output: 'vg/graphs/{dat}/{genome}-{vcf}-{chrom}.snarls'
    threads: 8
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.snarls.benchmark.log'
    log: 'vg/logs/{dat}/{genome}-{vcf}-{chrom}.snarls.log'
    shell:
        "vg snarls "
        "--threads {threads} "
        "{input} > {output} 2> {log}"

rule gbwt_haplo:
    input:
        vcf=config['vcf'][VCFV]['minMAF'],
        tbi=config['vcf'][VCFV]['minMAF']+ ".tbi",
        xg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.xg"
    output:
        gbwt="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gbwt",
    log:
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}-gbwt.log"
    threads: 32
    resources:
        mem_gc=200000
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}-gbwt.log"
    shell:
        "vg gbwt "
        "--temp-dir vg/tmp "
        "--xg-name {input.xg} "
        "--output {output.gbwt} "
        "--vcf-input {input.vcf} "
        "--progress 2> {log} "


rule prune_vg:
    """
    Pruning is used before constructing the GSCA index
    """
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.vg",
    output:
        "vg/graphs/{dat}/{genome}-{vcf}-{chrom}.pruned.vg"
    threads: 4
    benchmark: 
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.prune.benchmark.log"
    log: 
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}.prune.log"
    params:
        opt=OPTS
    shell:
        "vg prune "
        "--threads {threads} "
        "--max-degree 32 "
        "--restore-paths {input.vg} "
        "--progress 1> {output} 2> {log}"

rule gbwt_greedy:
    input:
        vcf=config['vcf'][VCFV]['minMAF'],
        xg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.xg"
    output:
        gbwt="vg/graphs/{dat}/{genome}-{vcf}-{chrom}-N{n}.gbwt",
        gg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}-N{n}.gg"
    threads: 64 # all threads
    resources:
        mem_mb=200000
    log:
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}-gbwt-N{n}.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}-gbwt-N{n}.log"
    shell:
        "vg gbwt "
        "--num-paths {wildcards.n} "
        "--temp-dir vg/tmp "
        "--xg-name {input.xg} "
        "--graph-name {output.gg} "
        "--output {output.gbwt} "
        "--path-cover "
        "--progress 2> {log} "

rule index_gcsa:
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.pruned.vg",
        #mapping='vg/graphs/{dat}/{genome}-{vcf}.ids.mapping', # -f
        gbwt="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gbwt"
    output:
        gcsa="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gcsa",
        gcsalcp="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gcsa.lcp"
    log:
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}.gcsa.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.gcsa.log"
    threads: 32
    resources:
        mem_gc=80000
    shell: 
        "vg index "
        "--temp-dir vg/tmp "
        "--gcsa-out {output.gcsa} "
        "--threads {threads} "
        "{input.vg} "
        "--progress 2> {log}"

rule index_minimizer:
    input:
        xg='vg/graphs/{dat}/{genome}-{vcf}-{chrom}.xg',
        gbwt='vg/graphs/{dat}/{genome}-{vcf}-{chrom}-N{n}.gbwt'
    output: 'vg/graphs/{dat}/{genome}-{vcf}-{chrom}.k{k}.w{w}.N{n}.min'
    threads: 64
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}-minimizer-k{k}-w{w}-N{n}.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}-{chrom}-minimizer-k{k}-w{w}-N{n}.log.txt'
    shell:
        "vg minimizer "
        "--kmer-length {wildcards.k} "
        "--window-length {wildcards.w} "
        "--threads {threads} "
        "--output-name {output} "
        "--gbwt-name {input.gbwt} "
        "{input.xg} 2> {log}"

rule index_trivial_snarls:
    input: 'vg/graphs/{dat}/{genome}-{vcf}-{chrom}.xg'
    output: 'vg/graphs/{dat}/{genome}-{vcf}-{chrom}.trivial.snarls'
    threads: 64 # use all threads
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}-trivialsnarls.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}-{chrom}-trivialsnarls.log.txt'
    shell:
        "vg snarls "
        "--threads {threads} "
        "--include-trivial "
        "{input} "
        "1> {output} 2> {log}"

rule index_distance:
    """
    We might need to remove -x flag - depracted 
    """
    input:
        xg='vg/graphs/{dat}/{genome}-{vcf}-{chrom}.xg',
        snarls='vg/graphs/{dat}/{genome}-{vcf}-{chrom}.trivial.snarls'
    output: 'vg/graphs/{dat}/{genome}-{vcf}-{chrom}.dist'
    threads: 64
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}-distance.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}-{chrom}-distance.log.txt'
    shell:
        "vg index "
        "--threads {threads} "
        "--dist-name {output} "
        "--xg-name {input.xg} "
        "--snarl-name {input.snarls} "
        "2> {log}"