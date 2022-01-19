###### conda activate aligns ######
# /projects/racimolab/data/MHC
#  snakemake --snakefile ../../people/gsd818/aHLA-vg/vg/rules/prep_pangenome_vcf_erik.smk --cores 35 --keep-going
import pandas as pd
import os, csv

wildcard_constraints:
    filters="[^.]+",
    maf="\d+"

CHROMS = list(range(1, 23)) + ['X']

CHR6_URL="http://hypervolu.me/~erik/HPRC/wgg.79"
SAMP=pd.read_table('pangenomes/vcf/samples_LOOC_all.txt', header=None)[0].tolist()

## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        # 1. Download
        # expand("vcf/pangenomes/chr{chrom}.pan.fa.a2fb268.4030258.bc221f9.smooth.vcf.gz", chrom=6)
        # 2. Get samples list to keep for vg
        # expand("pangenomes/vcf/samples_lists/samples_LOOC_wo_{sample}.txt",sample='HG00733'), # run first with one individual - EAS chosen
        # 3. Get filtered vcf
        # a. remove ind + conflict sites
        expand("pangenomes/vcf/chr6.pan.fa.smooth_wo_{sample}.{filters}.vcf.gz", sample='HG00733', filters='allfiltered'),
        # b. remvoe large sites and split multiallelics
        # expand("pangenomes/vcf/chr6.pan.fa.smooth_wo_{sample}.sites{len}kb.{filters}.vcf.gz", sample='HG02080', filters='allfiltered', len=5000),
        # c. remove maf filter
        expand("pangenomes/vcf/chr6.pan.fa.smooth_wo_{sample}.{filters}.maf{maf}.vcf.gz", sample='HG00733', filters='allfiltered', maf='01'),


## --------------------------------------------------------------------------------
## rules

rule download_chr6_pan:
    output:
        "pangenomes/vcf/chr{chrom}.pan.fa.a2fb268.4030258.bc221f9.smooth.vcf.gz"
    shell:
        "wget -O {output} {CHR6_URL}/chr{wildcards.chrom}.pan.fa.a2fb268.4030258.bc221f9.smooth.vcf.gz "


rule loov_txt:
    input:
        vcf="pangenomes/vcf/chr6.pan.fa.a2fb268.4030258.bc221f9.smooth.vcf.gz",
    output:
        all="pangenomes/vcf/samples_LOOC_all.txt",
    shell:
        "bcftools query -l {input.vcf} | head -n+44 > {output.all} " 
# head command to remove chm13 and grch38 columns

rule loov_samples:
    input:
        all="pangenomes/vcf/samples_LOOC_all.txt",
    output:
        sampl="pangenomes/vcf/samples_lists/samples_LOOC_wo_{sample}.txt"
    run:
        ## -----------------------------------------------------------------------------
        #d = open(input.all, "rt")
        d = open("pangenomes/vcf/samples_LOOC_all.txt", "rt")
        df= pd.read_csv(d, header=None)
        for i in range(len(df)):
            finalvec = df[~df[0].isin([df[0][i]])]
            finalvec.to_csv(output.sampl, header=None, index=False)
        ## -----------------------------------------------------------------------------

# create a file with all samples minus the one you will align it to

rule filter_chr6_pan:
    """
    Remove individual to test + conflict sites
    """
    input: 
        vcf="pangenomes/vcf/chr6.pan.fa.a2fb268.4030258.bc221f9.smooth.vcf.gz",
        sampl="pangenomes/vcf/samples_lists/samples_LOOC_wo_{sample}.txt"
    output:
        temp("pangenomes/vcf/chr6.pan.tmp.smooth_wo_{sample}.{filters}.vcf.gz")
    log:
        "pangenomes/vcf/chr6.pan.tmp.smooth_wo_{sample}_{filters}.log"
    threads: 4
    run:
        if wildcards.filters == "allfiltered":
            print("allfilters done")
            shell(""" bcftools view -S {input.sampl} --include ' INFO/CONFLICT = "" && INFO/LV == 0 ' --output-type z --output-file {output} {input.vcf} 2> {log} """)
        if wildcards.filters == "confiltered":
            shell(""" bcftools view -S {input.sampl} --include ' INFO/CONFLICT == "" ' --output-type z --output-file {output} {input.vcf} """)
        if wildcards.filters == "lvfiltered":
            shell(""" bcftools view -S {input.sampl} --include ' INFO/LV == 0 ' --output-type z --output-file {output} {input.vcf} """)

rule mod_chr6_name:
    input: 
        "pangenomes/vcf/chr6.pan.tmp.smooth_wo_{sample}.{filters}.vcf.gz"
    output:
        temp("pangenomes/vcf/chr6.pan.smooth_wo_{sample}.{filters}.vcf.gz"),
    shell:
        """
        zcat {input} | awk '{{gsub(/^grch38_grch38#1#/, ""); print}}' | bgzip -c > {output} ; bcftools index -t {output}
        """

rule split_multiallelic:
    input:
        "pangenomes/vcf/chr6.pan.smooth_wo_{sample}.{filters}.vcf.gz"
    output:
        "pangenomes/vcf/chr6.pan.fa.smooth_wo_{sample}.{filters}.vcf.gz"
    shell:
        "bcftools norm --multiallelics -any --output-type z --output {output} {input} ; "
        "bcftools index -t {output}"

rule maf_01:
    input:
        "pangenomes/vcf/chr6.pan.fa.smooth_wo_{sample}.{filters}.vcf.gz"
    output:
        "pangenomes/vcf/chr6.pan.fa.smooth_wo_{sample}.{filters}.{maf}maf.vcf.gz"
    threads: 4
    shell:
        "bcftools view "
        "--exclude 'MAF<0.01' " 
        " --output-type z"
        " --output-file {output}"
        " {input} ; "
        "bcftools index -t {output}"

rule rm_largeSites:
    """
    Remove large sites with more than 10kb alt/ref alleles. Use the file that w/w maf filter. Change {input.X}
    """
    input:
        vcff="pangenomes/vcf/chr6.pan.fa.tmp.smooth_wo_{sample}.{filters}.vcf.gz",
        vcfmaf="pangenomes/vcf/chr6.pan.fa.smooth_wo_{sample}.{filters}.{maf}maf.vcf.gz"
    output:
        "pangenomes/vcf/chr6.pan.fa.smooth_wo_{sample}.sites{len}kb.{filters}.vcf.gz"
    shell:
        "bcftools filter -e '(STRLEN(REF)>={wildcards.len}) | (STRLEN(ALT)>={wildcards.len})' --output-type z --output {output} {input.vcfmaf} ; "
        "bcftools index -t {output}"

# --max-alleles is number of atl alelles (from multiple alignment)
# read length:  bcftools filter -i '(STRLEN(REF)>=1000) | (STRLEN(ALT)>=1000)' 

# bcftools filter -e '(STRLEN(REF)>= 1000) |  (STRLEN(ALT)>= 1000)'  --output-type z --output vcf/pangenomes/chr6.pan.fa.smooth_wo_HG02080.sites1000kb.allfiltered.maf01.vcf.gz vcf/pangenomes/chr6.pan.fa.smooth_wo_HG02080.allfiltered.maf01.vcf.gz ; bcftools index -t vcf/pangenomes/chr6.pan.fa.smooth_wo_HG02080.sites1000kb.allfiltered.maf01.vcf.gz