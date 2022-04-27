#!/bin/bash

###### CHECK USAGE #####
if [ $# != 1 ]; then echo "$0: Usage: Please include a sample ID."; exit; fi
########################

##### MODULES #####
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load gatk/4.1.0.0
###################

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
DATADIR="${PROJDIR}/data/chrX_bam"
INDEXDIR="${PROJDIR}/analysis/0_index_genome"
OUTDIR="${PROJDIR}/analysis/4_duplicates"
#######################

##### VARIABLES #####
sample_id=$1
sample="${DATADIR}/${sample_id}.sorted.bam"
output_prefix="${DATADIR}/${sample_id}.sorted"
reference="${INDEXDIR}/canFam6.masked.fa"
#####################

#check for files
if [ ! -e "${sample}" ]; then echo -e "[4_duplicates]\t$sample does not exist."; exit; fi

if [ ! -e "${reference}" ]; then
  if [ -e "${reference}.gz" ]; then
    echo -e "[4_duplicates]\tReference is gzipped. Please unzip it and try again."
    exit
  else
    echo -e "[4_duplicates]\t${reference}(.gz) does not exist."
    exit
  fi
fi

#GATK command line
gatk --java-options "-Xmx1G" BuildBamIndex -I $sample -R $reference
gatk --java-options "-Xmx1G" MarkDuplicates -R $reference -I $sample -M "${output_prefix}.dup_metrics" -O "${output_prefix}.markedup.bam"
