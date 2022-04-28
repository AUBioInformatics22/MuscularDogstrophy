#!/bin/sh

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bwa/0.7.12
module load samtools/1.2
module load picard/1.79


#path variables
path="/scratch/MDPracticeData/analysis/1_index_genome/"

cd ${path}/

gunzip canFam6.masked.fa.gz

#index genome file
bwa index -p canFam6 canFam6.masked.fa
samtools faidx canFam6.masked.fa
java -Xms2g -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/CreateSequenceDictionary.jar REFERENCE=canFam6.masked.fa OUTPUT=canFam6.masked.dict

#recompress file
gzip canFam6.masked.fa
