# conda activate aligns - run from home aHLA directory
#################### file to run vg per chr - whole-genome graph alignment #################### 
# INFO from github #
# augment the graph with all variation from the GAM, saving each mapping as a path in the graph.
# softclips of alignment paths are preserved (`-S`).
# Note, this can be much less efficient than the above example if there are many alignments in the GAM
# vg augment x.vg aln.gam -i -S > aug_with_paths.vg


import pandas as pd
import os
configfile: "config.yaml"

# snakemake --snakefile vg/rules/aug_graph_chr6.smk --keep-going --rerun-incomplete --cores 40
## --------------------------------------------------------------------------------
## global parameters 
# CHROMS = list(range(1, 23)) + ['X']
# choose vcf
VCFV='pangenomes'
VCFPATH=config['vcf'][VCFV]['minMAF']
# choose reference 
REFV=config['ref_version']
REFPATH=config['ref'][REFV]

wildcard_constraints:
    vcf="[^-]+",
    chrom="[^-+\.?$]+",
## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        #expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.vg", chrom=6, dat='pangenomes', genome=REFV, vcf=VCFV),
        #expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.vg", chrom=6, dat='1kGP', genome=REFV, vcf='1kGP'),
        expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.snarls",chrom=6, dat='1kGP', genome=REFV, vcf='1kGP'),
        expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.xg",chrom=6, dat='1kGP', genome=REFV, vcf='1kGP'),
        expand("vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.gcsa",chrom=6, dat='1kGP', genome=REFV, vcf='1kGP'),

## --------------------------------------------------------------------------------
## rules

rule vg_map_aug:
    input:
        gcsa="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gcsa", 
        xg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.xg",
        gbwt="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gbwt",  # can be used for potentially more accuarte alignments 
        hla="reference/IMGTHLA/hla_gen.fasta"
    output:
        "vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gam"
    benchmark: 
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.map-aug.benchmark.log"
    threads: 32
    log: 
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}.map-aug.log"
    shell:
        "vg map -F {input.hla} -x {input.xg} -g {input.gcsa} -1 {input.gbwt} --debug --log-time 1> {output} 2> {log}"

rule vg_aug:
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.vg",
        gam="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.gam",
    output:
        vg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.vg",
        gam="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.gam",
    benchmark: 
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.aug.benchmark.log"
    log: 
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}.aug.log"
    shell:
        "vg augment --progress -i {input.vg} {input.gam} -A {output.gam} 1> {output.vg} 2> {log}"
#--min-coverage 1

# prepare the snarls index for the augmented graph
rule index_snarls_aug:
    input: 'vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.vg'
    output: 'vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.snarls'
    threads: 16
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}-aug-snarls.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}-{chrom}-aug-snarls.log.txt'
    run:
        shell('vg snarls -t {threads} {input} > {output} 2> {log}')

rule index_xg:
    """
    Make xg index containing the alts paths. 
    """
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.vg",
    output:
        xg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.xg",
    threads: 16
    resources:
        mem_xg=80000
    log:
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}.index.aug.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.index.aug.log"
    shell:
        "vg index --temp-dir vg/tmp -x {output.xg} -L {input} --progress 2> {log}"

rule index_gcsa:
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.vg",
    output:
        gcsa="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.gcsa",
        gcsalcp="vg/graphs/{dat}/{genome}-{vcf}-{chrom}.aug.gcsa.lcp"
    log:
        "vg/logs/{dat}/{genome}-{vcf}-{chrom}.aug.gcsa.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-{chrom}.aug.gcsa.log"
    threads: 40
    resources:
        mem_gc=200000
    shell: 
        "vg index "
        "--temp-dir vg/tmp "
        "-g {output.gcsa} "
        "-t {threads} "
        "{input.vg} --progress 2> {log}"


