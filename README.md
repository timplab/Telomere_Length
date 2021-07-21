# Telomere_Length
This repository contains the pipeline for assessing telomere length from "Chromosome specific telomere lengths and the minimal functional telomere revealed by nanopore sequencing"

# Dependencies 
- Minimap2 (https://github.com/lh3/minimap2)
- TideHunter (https://github.com/Xinglab/TideHunter)
- Porechop (https://github.com/rrwick/Porechop)
- samtools (http://www.htslib.org/download/)
- bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html)
- snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

# Install Telomere_Length
```
git clone git@github.com:timplab/Telomere_Length.git
cd Telomere_Length
vim config.yaml
```
Edit paths in `congfig.yaml` to the location of software installs and input files 

# Run Telomere_Length

```
snakemake --snakefile run_telo.snk --cores 8

```
Additionally, all scripts can run as stand alone executeables (filter.sh and telo_length.R) that take user input from TideHunter and minimap2. Run with -h for instructions.
