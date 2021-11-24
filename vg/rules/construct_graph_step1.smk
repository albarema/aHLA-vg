##### step 1 #####
configfile: "config.yaml"
CHROMS = list(range(1, 23)) # + ['X']
# choose prune options - NOT WORKING - tab paste as well with the options 
PRUNEOP=config['prune_version']
OPTS=config['prune_options'][PRUNEOP]
# choose vcf 
VCFV=config['vcf_version']
VCFPATH=config['vcf'][VCFV]['allminMAF']
# choose reference
REFV=config['ref_version']
REFPATH=config['ref'][REFV]

wildcard_constraints:
    vcf="[^-+\.?$]+",
    chrom="[^-+\.?$]+",
    ext="[^.]+",

rule construct_chr:
    """
    constructing graph
    """
    input:
        ref=config['ref'][REFV],
        vcf=config['vcf'][VCFV]['allminMAF'],
        tbi=config['vcf'][VCFV]['allminMAF']+ ".tbi"
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

rule node_coord:
    """
    Generate a joint id space across each chromosome graph
    """
    input:
        expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.vg", chrom=CHROMS, allow_missing=True)
    output:
        mapping="vg/graphs/{dat}/{genome}-{vcf}.ids.mapping",
        ba="vg/graphs/{dat}/{genome}-{vcf}.ids.mapping.backup"
    log:
        "vg/logs/{dat}/{genome}-{vcf}.node_coord.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}.node_coord.log"
    shell:
        "touch {output.mapping} ; "
        "vg ids --join --mapping {output.mapping} {input} 2> {log} ; "
        "cp {output.mapping} {output.ba}"

rule index_xg:
    """
    Make xg index containing the alts paths
    """
    input:
        vg=expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.vg", chrom=CHROMS, allow_missing=True),
        map="vg/graphs/{dat}/{genome}-{vcf}.ids.mapping",
    output:
        "vg/graphs/{dat}/{genome}-{vcf}.xg"
    log:
        "vg/logs/{dat}/{genome}-{vcf}.xg.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}.xg.log"
    shell:
        "vg index --temp-dir vg/tmp -x {output} -L {input.vg} --progress 2> {log}"

rule index_snarls:
    """
    prepare the snarls index
    """
    input: 'vg/graphs/{dat}/{genome}-{vcf}.xg'
    output: 'vg/graphs/{dat}/{genome}-{vcf}.snarls'
    threads: 8
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}.snarls.benchmark.log'
    log: 'vg/logs/{dat}/{genome}-{vcf}.snarls.log'
    shell:
        "vg snarls -t {threads} {input} > {output} 2> {log}"

rule gbwt_haplo:
    input:
        vcf="vcf/{dat}/{vcf}.vcf.gz",
        xg="vg/graphs/{dat}/{genome}-{vcf}.xg"
    output:
        gbwt="vg/graphs/{dat}/{genome}-{vcf}.gbwt",
    log:
        "vg/logs/{dat}/{genome}-{vcf}-gbwt.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-gbwt.log"
    shell:
        "vg gbwt -d vg/tmp -x {input.xg} -o {output.gbwt} -v {input.vcf} --progress 2> {log} "

rule prune_vg:
    """
    Pruning is used before constructing the GSCA index
    """
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.vg",
        mapping="vg/graphs/{dat}/{genome}-{vcf}.ids.mapping"
    output: 
        "vg/graphs/{dat}/{genome}-{vcf}-{chrom}.pruned.vg"
    threads: 4
    benchmark: 
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.prune.benchmark.log"
    log: 
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}.prune.log"
    params:
        opt=OPTS #not used cause it was not working
    shell:
        "vg prune -t {threads} -M 32 --restore-paths {input.vg} --progress 1> {output} 2> {log}"

rule gbwt_greedy:
    input:
        vcf=config['vcf'][VCFV]['allminMAF'],
        xg="vg/graphs/{dat}/{genome}-{vcf}.xg"
    output:
        gbwt="vg/graphs/{dat}/{genome}-{vcf}-N{n}.gbwt",
        gg="vg/graphs/{dat}/{genome}-{vcf}-N{n}.gg"
    threads: 64 # all threads
    resources:
        mem_mb=200000
    log:
        "vg/logs/{dat}/{genome}-{vcf}-gbwt-N{n}.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-gbwt-N{n}.log"
    shell:
        "vg gbwt -n {wildcards.n} -d vg/tmp -x {input.xg} -g {output.gg} -o {output.gbwt} -P --progress 2> {log} "


rule index_gcsa:
    input:
        vg=expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.pruned.vg", chrom=CHROMS, allow_missing=True),
        mapping='vg/graphs/{dat}/{genome}-{vcf}.ids.mapping' # -f
    output:
        gcsa="vg/graphs/{dat}/{genome}-{vcf}.gcsa",
        gcsalcp="vg/graphs/{dat}/{genome}-{vcf}.gcsa.lcp"
    params:
        opt=config['gcsa_options']['highdegree']
    log:
        "vg/logs/{dat}/{genome}-{vcf}.gcsa.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}.gcsa.log"
    threads: 64
    resources:
        mem_gc=100000
    shell: 
        "vg index "
        "--temp-dir vg/tmp "
        "-g {output.gcsa} "
        "-t {threads} "
        "{input.vg} --progress 2> {log}"

rule index_minimizer:
    input:
        xg='vg/graphs/{dat}/{genome}-{vcf}.xg',
        gbwt='vg/graphs/{dat}/{genome}-{vcf}-N{n}.gbwt'
    output: 'vg/graphs/{dat}/{genome}-{vcf}.k{k}.w{w}.N{n}.min'
    threads: 64
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-minimizer-k{k}-w{w}-N{n}.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}-minimizer-k{k}-w{w}-N{n}.log.txt'
    shell:
        "vg minimizer -k {wildcards.k} -w {wildcards.w} -t {threads} -i {output} -g {input.gbwt} {input.xg} 2> {log}"

rule index_trivial_snarls:
    input: 'vg/graphs/{dat}/{genome}-{vcf}.xg'
    output: 'vg/graphs/{dat}/{genome}-{vcf}.trivial.snarls'
    threads: 64 # use all threads
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-trivialsnarls.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}-trivialsnarls.log.txt'
    shell:
        "vg snarls -t {threads} --include-trivial {input} > {output} 2> {log}"

rule index_distance:
    """
    We might need to remove -x flag - depracted 
    """
    input:
        xg='vg/graphs/{dat}/{genome}-{vcf}.xg',
        snarls='vg/graphs/{dat}/{genome}-{vcf}.trivial.snarls'
    output: 'vg/graphs/{dat}/{genome}-{vcf}.dist'
    threads: 64
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-distance.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}-distance.log.txt'
    shell:
        "vg index -t {threads} -j {output} -x {input.xg} -s {input.snarls} 2> {log}"