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
        expand("vg/graphs/GRCh38_no_alt_hap.{dat}.aug.vg", dat=['1kGP'])

## --------------------------------------------------------------------------------
## rules

rule vg_map_aug:
    input:
        gcsa="vg/graphs/GRCh38_no_alt_hap.{dat}.gw.gcsa", 
        xg="vg/graphs/{dat}.graph-with-alts.xg",
        gbwt="vg/graphs/{dat}.all.gbwt",  # can be used for potentially more accuarte alignments 
        hla="reference/IMGTHLA/hla_gen.fasta"
    output:
        "aln.{dat}.gw.gam"
    log:
        "vg/logs/GRCh38_no_alt_hap.{dat}.vg.map.aug.log",
    shell:
        "vg/vg map -F {input.hla} -x {input.xg} -g {input.gcsa} -1 {input.gbwt} --debug --log-time 1> {output} 2> {log}"

rule vg_aug:
    input:
        vg=expand("vg/graphs/GRCh38_no_alt_hap.{dat}.{chrom}.vg", chrom=CHROMS, allow_missing=True),
        gam="aln.{dat}.gw.gam",
    output:
        vg="vg/graphs/GRCh38_no_alt_hap.{dat}.aug.vg",
        gam="aug.{dat}.gw.gam",
    log:
        log="vg/logs/GRCh38_no_alt_hap.{dat}.vg.aug.log",
    shell:
        "vg/vg augment --min-converage 1  -i {input.vg} {input.gam} -A {output.gam} --progress 1> {output.vg} 2> {log}"


# augment the graph with all variation from the GAM, saving each mapping as a path in the graph.
# softclips of alignment paths are preserved (`-S`).
# Note, this can be much less efficient than the above example if there are many alignments in the GAM
# vg augment x.vg aln.gam -i -S > aug_with_paths.vg