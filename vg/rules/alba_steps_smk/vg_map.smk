# conda activate aligns - run from home aHLA directory
#################### file to run vg per chr - whole-genome graph alignment #################### 

import pandas as pd
import os
configfile: "config.yaml"

# snakemake --snakefile vg/rules/vg_aug.smk --keep-going --rerun-incomplete --cores 40
## --------------------------------------------------------------------------------
## global parameters 
CHROMS = list(range(1, 23)) + ['X']

## --------------------------------------------------------------------------------
## targets

rule all:
    input:

## --------------------------------------------------------------------------------
## rules

rule vg_map_single_end:
    """
    Align a single-end fastq file to the indexed version of the graph

    Using `-k 15 -w 1024` from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02160-7#Sec6
    """
    input:
        xg="data/vg/index/{panel}.xg",
        gcsa="data/vg/index/{panel}.gcsa",
        fastq="data/samples/fastq/adapter_rm/{accession}.collapsed.fastq.gz",
    output:
        gam="data/samples/bam/{panel}.vg_map.{accession}.gam",
    log:
        log="data/samples/bam/{panel}.vg_map.{accession}.gam.log",
    params:
        id=lambda wildcards, input: get_read_group_identifier(input.fastq),
        sm=lambda wildcards, input: get_read_group_sample(input.fastq),
        lb=lambda wildcards, input: get_read_group_library(input.fastq),
    threads: workflow.cores / 4
    shell:
        # TODO add read group headers
        "vg map"
        " --xg-name {input.xg}"
        " --gcsa-name {input.gcsa}"
        " --fastq {input.fastq}"
        " --read-group '@RG\\tID:{params.id}\\tSM:{params.sm}\\tCN:{RG_CENTER_NAME}\\tPL:{RG_PLATFORM}\\tLB:{params.lb}\\tDS:{RG_DESCRIPTION}'"
        " --exclude-unaligned"
        " --min-mem 15"
        " --band-width 1024"
        " --threads {threads}"
        " --log-time"
        " 1> {output.gam} 2> {log}"