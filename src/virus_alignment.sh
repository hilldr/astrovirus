#! /bin/bash
## KALLISTO PROCESSING SCRIPT (astrovirus)
## David R. Hill
## -----------------------------------------------------------------------------
## Setup variables
## start time
DT1=$(date '+%d/%m/%Y %H:%M:%S')
## location of genome index file on Spence lab server
INDEX=../data/genomes/GCF_000885815.1_ViralProj39811_cds_from_genomic.fna.gz.idx
## email address for notifications
EMAIL=d2.david.hill@gmail.com
## This is the directory where the kallisto results will be deposited
RESULTDIR=../results/astrovirus_va1/DESeq2/
## make the folder to deposit results
mkdir -p $RESULTDIR

## this is the directory that contains the fastq directories
for dir in ../data/Run_2127/wobus/*
## for loop will iterate through each directory and find fastq files and run
## kallisto with specified arguments
do
    for file in $dir/*.fastq*
    do
    SHORTNAME=$(basename "$file")
    NAME2="${SHORTNAME##*/}"	
    DIRNAME="${NAME2%.*}"  
    # These settings are for single-end 50 bp reads
    kallisto quant -i $INDEX --output-dir=$RESULTDIR/$DIRNAME --threads=8 \
	     --bootstrap-samples=100 --single --fragment-length=50 --sd=1 $file
    done
done

## Send email notification of script completion
DT2=$(date '+%d/%m/%Y %H:%M:%S')
echo "Kalliso run initiated at $DT1 complete at $DT2" | mail -s "Kallisto complete" $EMAIL
