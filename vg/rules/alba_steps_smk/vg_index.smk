# conda activate aligns - run from home aHLA directory
#  top -d 1 -b | grep vg >> vg.log
#################### file to run vg per chr - whole-genome graph alignment #################### 

import pandas as pd
import os
configfile: "config.yaml"

# snakemake --snakefile vg/rules/vg.smk --keep-going --rerun-incomplete --cores 40
## --------------------------------------------------------------------------------
## global parameters 
CHROMS = list(range(1, 23)) + ['X']

## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        #expand("vg/graphs/GRCh38_no_alt_hap.{dat}.{chrom}.vg", chrom=CHROMS, dat=['1kGP'])
        expand("vg/graphs/GRCh38_no_alt_hap.{dat}.gw.gcsa", dat=['1kGP'])



# a. simulations/GRCh38_reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# b. simulations/GRCh38_contigs/GRCh38_full_analysis_set_plus_decoy_hla.fa
## --------------------------------------------------------------------------------
## rules

"""
rule run_vg:
    """
    constructing graph
    """
    input:
        fa=config['ref_v']['only_chrs'], # get from config file the right fasta version 
        vcf="{dat}/vcf/{dat}.allchr.filtered.vcf.gz", # "{dat}/vcf/{dat}.chr{chrom}.filtered.vcf.gz",
        tbi="{dat}/vcf/{dat}.allchr.filtered.vcf.gz.tbi", #"{dat}/vcf/{dat}.chr{chrom}.filtered.vcf.gz.tbi",
    output:
        "vg/graphs/GRCh38_no_alt_hap.{dat}.{chrom}.vg"
    log:
        "vg/logs/GRCh38_no_alt_hap.{dat}.{chrom}.vg_construct.log"
    benchmark:
        "vg/benchmarks/GRCh38_no_alt_hap.{dat}.{chrom}.vg_construct.benchmark.log"
    shell:
        "vg construct "
        "-C -R chr{wildcards.chrom} "
        "-r {input.fa} "
        "-v {input.vcf} "
        "-t 4 "
        "-m 32 "
        "-a "
        "--progress "
        "1> {output} 2> {log} "
"""

rule node_coord:
    """
    Generate a joint id space across each chromosome graph
    """
    input:
        expand("vg/graphs/GRCh38_no_alt_hap.{dat}.{chrom}.vg", chrom=CHROMS, allow_missing=True)
    output:
        map="vg/graphs/GRCh38_no_alt_hap.{dat}.mapping",
        ba="vg/graphs/GRCh38_no_alt_hap.{dat}.mapping.backup"
    log:
        "vg/logs/GRCh38_no_alt_hap.{dat}.node_coord.log"
    benchmark:
        "vg/benchmarks/GRCh38_no_alt_hap.{dat}.node_coord.log"
    shell:
        "touch {output} ; "
        "vg ids -j -m {output.map} {input} 2> {log} ; "
        "cp {output.map} {output.ba}"

rule with_alt:
    input:
        expand("vg/graphs/GRCh38_no_alt_hap.{dat}.{chrom}.vg", chrom=CHROMS, allow_missing=True)
    output:
        "vg/graphs/GRCh38_no_alt_hap.{dat}.graph-with-alts.xg"
    log:
        "vg/logs/GRCh38_no_alt_hap.{dat}.with_alt.log"
    benchmark:
        "vg/benchmarks/GRCh38_no_alt_hap.{dat}.with_alt.log"
    shell:
        "vg index -x {output} -L {input} --progress 2> {log}"

rule with_gbwt:
    input:
        vcf=expand("{dat}/vcf/{dat}.allchr.filtered.vcf.gz", chrom=CHROMS, allow_missing=True),
        alts="vg/graphs/GRCh38_no_alt_hap.{dat}.graph-with-alts.xg"
    output:
        "vg/graphs/GRCh38_no_alt_hap.{dat}.all.gbwt"
    log:
        "vg/logs/GRCh38_no_alt_hap.{dat}.with_gbwt.log"
    benchmark:
        "vg/benchmarks/GRCh38_no_alt_hap.{dat}.with_gbwt.log"
    shell:
        "vg gbwt -d vg/tmp -x {input.alts} -o {output} -v {input.vcf} --progress 2> {log}"

rule gcsa_index:
    input:
        vg="vg/graphs/GRCh38_no_alt_hap.{dat}.{chrom}.vg",
        gbwt="vg/graphs/GRCh38_no_alt_hap.{dat}.all.gbwt",
        map="vg/graphs/GRCh38_no_alt_hap.{dat}.mapping.backup"
    output:
        vg=temp("vg/graphs/GRCh38_no_alt_hap.{dat}.{chrom}.pruned.vg"),
    threads: 30
    log:
        "vg/logs/GRCh38_no_alt_hap.{dat}.{chrom}.gcsa_index.log"
    benchmark:
        "vg/benchmarks/GRCh38_no_alt_hap.{dat}.{chrom}.gcsa_index.log"
    shell:
        "cp {input.map} vg/graphs/GRCh38_no_alt_hap.{wildcards.dat}.mapping ; "
        "vg prune -u -g {input.gbwt} -t {threads} -a -m vg/graphs/GRCh38_no_alt_hap.{wildcards.dat}.mapping {input.vg} --progress 1> {output.vg} 2> {log}"

rule wg_gcsa:
    input:
        vg=expand("vg/graphs/GRCh38_no_alt_hap.{dat}.{chrom}.pruned.vg", chrom=CHROMS, allow_missing=True),
        map="vg/graphs/GRCh38_no_alt_hap.{dat}.mapping"
    output:
        "vg/graphs/GRCh38_no_alt_hap.{dat}.gw.gcsa"
    log:
        "vg/logs/GRCh38_no_alt_hap.{dat}.gw.gcsa.log"
    benchmark:
        "vg/benchmarks/GRCh38_no_alt_hap.{dat}.gw.gcsa.log"
    threads: 60
    shell: 
        "vg index -b vg/tmp -g {output} -t {threads} -f {input.map} {input.vg} --progress 2> {log}"

