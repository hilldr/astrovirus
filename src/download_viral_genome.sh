#! /bin/bash 
## download and index astrovirus genome
## David R. Hill 2018-09-21
## -----------------------------------------------------------------------------

## make directory for genomes
mkdir -p ../data/genomes

## Astrovirus VA1
## GenBank: 4731478
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/885/815/GCF_000885815.1_ViralProj39811/GCF_000885815.1_ViralProj39811_cds_from_genomic.fna.gz ../data/genomes

## index genomes with kallisto
for file in ../data/genomes/*.fna.gz
do
    kallisto index -i $file\.idx $file
done
