# prepare files to run vg - run from aHLA directory (where config file is). 

import pandas as pd
import os, csv
configfile: "config.yaml"


## --------------------------------------------------------------------------------
## global parameters 
CHROMS = list(range(1, 23)) + ['X']

## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        # by chr, download
        # expand("1kGP/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz", chrom=CHROMS),
        # by chr, filter
        expand("1kGP/vcf/1kGP.chr{chrom}.filtered.vcf.gz", chrom=CHROMS),
        # all chr
        #expand("{dat}/vcf/{dat}.allchr.filtered.vcf.gz", dat=['1kGP']), # hgdp

## --------------------------------------------------------------------------------
## rules

# VCF 1000 Genomes Project
GP_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"

rule get_mds5:
    input:
        md5="1kGP/phased-manifest_July2021.tsv"
    output:
        "1kGP/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz.md5"
    run:
        ## -----------------------------------------------------------------------------
        d = open(input.md5, "rt")
        reader= csv.reader(d, delimiter='\t')
        for i in reader:
            filename = str("1kGP/vcf/")+i[0]
            mdname = filename+str(".md5")
            outfile=open(mdname, "w")
            finalvec= i[2] + "\t" + filename
            outfile.write(''.join(finalvec) + '\n')
            outfile.close()
        ## -----------------------------------------------------------------------------

# Need to download chr X - different path 
rule download_vcf:
    """
    Download the 1000 Genomes Project phase 3 VCFs
    """
    input:
        md5="1kGP/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz.md5"
    output:
        bgz="1kGP/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz",
    threads: 1
    shell:
        "wget --quiet -O {output.bgz} {GP_URL}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz && "
        "md5sum --status --check {input.md5}"

rule download_tbi:
    input:
        md5="1kGP/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi.md5"
    output:
        tbi="1kGP/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"
    threads: 1
    shell:
        "wget --quiet -O {output.tbi} {GP_URL}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi && "
        "md5sum --status --check {input.md5}"

rule filter_vcf_1000gp:
    """
    Filter sites (maf 1%) + only 2504 samples - remove relate individuals
    """
    input:
        vcf="1kGP/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz",
        tbi="1kGP/vcf/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi",
        sampl="1kGP/samples_1kGP_2504.txt"
    output:
        vcf="1kGP/vcf/1kGP.chr{chrom}.filtered.vcf.gz",
        tbi="1kGP/vcf/1kGP.chr{chrom}.filtered.vcf.gz.tbi"
    threads: 4
    shell:
        "bcftools view -S {input.sampl} "
        "--exclude 'MAF<0.01' " 
        " --output-type z"
        " --output-file {output.vcf}"
        " {input.vcf} && "
        "bcftools index --tbi {output.vcf} "


rule merge:
    input:
        expand("{dat}/vcf/{dat}.chr{chrom}.filtered.vcf.gz", chrom=CHROMS, allow_missing=True)
    output:
        "{dat}/vcf/{dat}.allchr.filtered.vcf.gz"
    shell:
        "bcftools concat -Oz -o {output} {input} ; "
        "bcftools index --tbi {output}"

rule get_vcf:
    """
    Filter hgdp vars 
    """
    input:
        os.path.join(config['vcf_hgdp'], "vcf/hgdp_wgs.20190516.full.chr{chrom}.vcf.gz")
    output:
        temp("hgdp/vcf/hgdp.chr{chrom}.filtered.vcf.gz") 
    shell:
        "bcftools view -e 'MAF<0.01' "
        "{input} " 
        "-Oz > {output} "

# GRCH38# GRCh38
# a. simulations/GRCh38_reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# b. simulations/GRCh38_contigs/GRCh38_full_analysis_set_plus_decoy_hla.fa
