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
# VCFV='pangenomes'
# VCFPATH=config['vcf'][VCFV]['minMAF']
# choose reference 
REFV=config['ref_version']
REFPATH=config['ref'][REFV]
PRUNEOP=config['prune_version']
OPTS=config['prune_options'][PRUNEOP]
wildcard_constraints:
    vcf="[^-]+",
## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        expand('vg/graphs/{dat}/{genome}-{vcf}.aug.{ext}', dat='1kGP', chrom=6, genome=REFV, vcf='1kGP', ext=['xg', 'trivial.snarls', 'dist']),
        expand("vg/graphs/{dat}/{genome}-{vcf}.aug.{ext}", dat='1kGP', genome=REFV, vcf='1kGP', ext=['snarls', 'xg', 'gcsa']),
        expand('vg/graphs/{dat}/{genome}-{vcf}.k{k}.w{w}.N{n}.aug.min',dat='1kGP', chrom=6, genome=REFV, vcf='1kGP', k=config['mink'],
               w=config['minw'], n=config['covern']),

## --------------------------------------------------------------------------------
## rules

rule vg_map_aug:
    input:
        gcsa="vg/graphs/{dat}/{genome}-{vcf}.gcsa", 
        xg="vg/graphs/{dat}/{genome}-{vcf}.xg",
        gbwt="vg/graphs/{dat}/{genome}-{vcf}-N16.gbwt",  # can be used for potentially more accuarte alignments 
        hla="/maps/projects/racimolab/data/MHC/reference/IMGT-HLA/IMGTHLA/hla_gen.fasta"
    output:
        "vg/graphs/{dat}/{genome}-{vcf}.gam"
    benchmark: 
        "vg/benchmarks/{dat}/{genome}-{vcf}.map-aug.benchmark.log"
    threads: 32
    log: 
        "vg/logs/{dat}/{genome}-{vcf}.map-aug.log"
    shell:
        "vg map --fasta {input.hla} "
        "--xg-name {input.xg} "
        "--gcsa-name {input.gcsa} "
        "--gbwt-name {input.gbwt} "
        "--debug "
        "--log-time "
        "1> {output} "
        "2> {log}"

rule vg_aug:
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}-6.vg",
        gam="vg/graphs/{dat}/{genome}-{vcf}.gam",
    output:
        vg="vg/graphs/{dat}/{genome}-{vcf}.aug.vg",
        gam="vg/graphs/{dat}/{genome}-{vcf}.aug.gam",
    benchmark: 
        "vg/benchmarks/{dat}/{genome}-{vcf}.aug.benchmark.log"
    threads: 32
    resources:
        mem_mb=250000
    log: 
        "vg/logs/{dat}/{genome}-{vcf}.aug.log"
    shell:
        "vg augment "
        "--progress "
        "--include-paths "
        "--subgraph "
        "{input.vg} "
        "{input.gam} "
        "--alignment-out {output.gam} "
        "1> {output.vg} "
        "2> {log}"
#--min-coverage 1

# prepare the snarls index for the augmented graph
rule index_snarls_aug:
    input: 'vg/graphs/{dat}/{genome}-{vcf}.aug.vg'
    output: 'vg/graphs/{dat}/{genome}-{vcf}.aug.snarls'
    threads: 16
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-aug-snarls.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}-aug-snarls.log.txt'
    run:
        shell('vg snarls -t {threads} {input} > {output} 2> {log}')

rule index_xg:
    """
    Make xg index containing the alts paths. 
    """
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}.aug.vg",
    output:
        xg="vg/graphs/{dat}/{genome}-{vcf}.aug.xg",
    threads: 16
    resources:
        mem_xg=80000
    log:
        "vg/logs/{dat}/{genome}-{vcf}.index.aug.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}.index.aug.log"
    shell:
        "vg index "
        "--temp-dir vg/tmp "
        "--xg-name {output.xg} "
        "--xg-alts {input} "
        "--progress 2> {log}"


rule prune_vg:
    """
    Pruning is used before constructing the GSCA index
    """
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}.aug.vg",
    output:
        "vg/graphs/{dat}/{genome}-{vcf}.aug.pruned.vg"
    threads: 4
    benchmark: 
        "vg/benchmarks/{dat}/{genome}-{vcf}.aug.prune.benchmark.log"
    log: 
        "vg/logs/{dat}/{genome}-{vcf}.aug.prune.log"
    params:
        opt=OPTS
    shell:
        "vg prune "
        "--threads {threads} "
        "--max-degree 32 "
        "--restore-paths {input.vg} "
        "--progress 1> {output} 2> {log}"

rule index_gcsa:
    input:
        vg="vg/graphs/{dat}/{genome}-{vcf}.aug.pruned.vg",
    output:
        gcsa="vg/graphs/{dat}/{genome}-{vcf}.aug.gcsa",
        gcsalcp="vg/graphs/{dat}/{genome}-{vcf}.aug.gcsa.lcp"
    log:
        "vg/logs/{dat}/{genome}-{vcf}.aug.gcsa.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}.aug.gcsa.log"
    threads: 86
    resources:
        mem_gc=700000
    shell: 
        "vg index "
        "--temp-dir vg/tmp "
        "--gcsa-out {output.gcsa} "
        "--threads {threads} "
        "--size-limit 10000 "
        "{input.vg} "
        "--progress 2> {log}"



rule gbwt_greedy:
    input:
        xg="vg/graphs/{dat}/{genome}-{vcf}.aug.xg"
    output:
        gbwt="vg/graphs/{dat}/{genome}-{vcf}-N{n}.aug.gbwt",
        gg="vg/graphs/{dat}/{genome}-{vcf}-N{n}.aug.gg"
    threads: 64 # all threads
    resources:
        mem_mb=200000
    log:
        "vg/logs/{dat}/{genome}-{vcf}-gbwt-N{n}.aug.log"
    benchmark:
        "vg/benchmarks/{dat}/{genome}-{vcf}-gbwt-N{n}.aug.log"
    shell:
        "vg gbwt "
        "--num-paths {wildcards.n} "
        "--temp-dir vg/tmp "
        "--xg-name {input.xg} "
        "--graph-name {output.gg} "
        "--output {output.gbwt} "
        "--path-cover "
        "--progress 2> {log} "

rule index_minimizer:
    input:
        xg='vg/graphs/{dat}/{genome}-{vcf}.aug.xg',
        gbwt='vg/graphs/{dat}/{genome}-{vcf}-N{n}.aug.gbwt'
    output: 'vg/graphs/{dat}/{genome}-{vcf}.k{k}.w{w}.N{n}.aug.min'
    threads: 64
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}-minimizer-k{k}-w{w}-N{n}.aug.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}-minimizer-k{k}-w{w}-N{n}.aug.log.txt'
    shell:
        "vg minimizer "
        "--kmer-length {wildcards.k} "
        "--window-length {wildcards.w} "
        "--threads {threads} "
        "--output-name {output} "
        "--gbwt-name {input.gbwt} "
        "{input.xg} 2> {log}"

rule index_trivial_snarls:
    input: 'vg/graphs/{dat}/{genome}-{vcf}.aug.xg'
    output: 'vg/graphs/{dat}/{genome}-{vcf}.aug.trivial.snarls'
    threads: 64 # use all threads
    resources:
        mem_mb=120000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}.aug-trivialsnarls.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}.aug-trivialsnarls.log.txt'
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
        xg='vg/graphs/{dat}/{genome}-{vcf}.aug.xg',
        snarls='vg/graphs/{dat}/{genome}-{vcf}.aug.trivial.snarls'
    output: 'vg/graphs/{dat}/{genome}-{vcf}.aug.dist'
    threads: 64
    resources:
        mem_mb=500000
    benchmark: 'vg/benchmarks/{dat}/{genome}-{vcf}.aug-distance.benchmark.txt'
    log: 'vg/logs/{dat}/{genome}-{vcf}.aug-distance.log.txt'
    shell:
        "vg index "
        "--threads {threads} "
        "--dist-name {output} "
        "--xg-name {input.xg} "
        "--snarl-name {input.snarls} "
        "2> {log}"