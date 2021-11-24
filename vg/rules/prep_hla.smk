import pandas as pd
import os, csv
configfile: "config.yaml"

## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        "reference/IMGTHLA/hla_gen.fasta"


## --------------------------------------------------------------------------------
## rules

HLA_URL="https://github.com/ANHIG/IMGTHLA/raw/Latest/fasta/hla_gen.fasta"


rule download_hla_imgtdatabase:
    input:
        md5="reference/IMGTHLA/hla_gen.fasta.md5"
    output:
        fa="reference/IMGTHLA/hla_gen.fasta"
    shell:
        "wget -O {output.fa} {HLA_URL} && "
        "md5sum --status --check {input.md5}"

