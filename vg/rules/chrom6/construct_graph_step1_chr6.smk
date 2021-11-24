##### step 1 #####
# NB: already run it for 1kG, now running it for the pangenomes ones. Version without haplo info 
# try to use haplotypes with mapping file - maybe add to dat wildcard 1kGP-mapping or pangenome-mapping
configfile: "config.yaml"

CHROMS = 6
# choose prune options
PRUNEOP=config['prune_version']
OPTS=config['prune_options'][PRUNEOP]
# choose vcf 
VCFV='pangenomes'
VCFPATH=config['vcf'][VCFV]['minMAF']
# choose reference
REFV=config['ref_version']
REFPATH=config['ref'][REFV]

wildcard_constraints:
    vcf="[^-]+",
    chrom="[^-+\.?$]+",

rule all:
    input:
        #expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gcsa", chrom=6, dat='pangenomes', genome=REFV, vcf=VCFV),
        expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.snarls", chrom=6, dat='pangenomes', genome=REFV, vcf=VCFV)
        # expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gcsa", chrom=6, dat='1kGP', genome=REFV, vcf=VCFV),


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
        "vg index --temp-dir vg/tmp -x {output.xg} -L {input} --progress 2> {log}"

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
        "vg gbwt -d vg/tmp -x {input.xg} -o {output.gbwt} -v {input.vcf} --progress 2> {log} "


rule prune_vg:
    """
    Pruning is used before constructing the GSCA index
    """
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.vg",
        #mapping="vg/graphs/{dat}/{genome}-{vcf}.ids.mapping" -m {input.mapping}
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
        "vg prune -t {threads} -M 32 --restore-paths {input.vg} --progress 1> {output} 2> {log}"


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
        "-g {output.gcsa} "
        "-t {threads} "
        "{input.vg} --progress 2> {log}"

#    params:
#        opt=config['gcsa_options']['highdegree']