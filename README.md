# aHLA-vg
Generate a genome variation graph to map ancient DNA. The graph can be augmented using HLA variation to achieve better accuracy for mapping hypervariable regions (e.g.: HLA). The pipeline is written in snakemake. 

vg pipeline:
- preparation of files
- construction of graphs
- augmentation of graph (optional)
- mapping (output format: bam files)
