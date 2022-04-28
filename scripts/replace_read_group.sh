#!/bin/bash

##### MODULES #####
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load gatk/4.1.0.0
module load samtools/1.3.1
###################

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
DATADIR="${PROJDIR}/data/fastq"
#BAMDIR="${PROJDIR}/data/bam"
#INDEXDIR="${PROJDIR}/analysis/0_index_genome"
#CHRXDIR="${PROJDIR}/data/chrX_bam"
#ALIGNDIR="${PROJDIR}/analysis/3_alignment"
#SCRATCHDIR="/scratch/MD_data"
#######################

##### VARIABLES #####
#sample="0001"
#ref="canFam6"
flowcell="H2LCVDSX3"
lane_id="3"
#flowcell=$(zcat ${DATADIR}/H2LCVDSX3_s1_1_IDT8_UDI_*_i7-IDT8_UDI_*_i5_6797-RN-${sample}.fastq.gz | head -1 | awk 'BEGIN {FS = ":"} {print $3}')
#lane_id=$(zcat ${DATADIR}/H2LCVDSX3_s1_1_IDT8_UDI_*_i7-IDT8_UDI_*_i5_6797-RN-${sample}.fastq.gz | head -1 | awk 'BEGIN {FS = ":"} {print $4}')
#genome_size="$(awk '{sum+=$2} END {print sum}' ${INDEXDIR}/${ref}.masked.fa.fai)"
#####################

###################
# Read Group Header
# ID = $sample.$flowcell.$lane_id
# SM (sample) = $sample
# PU (platform unit) = $flowcell.$lane_id
# PL (platform) = illumina
# LB (library) = $flowcell.$lane_id
###################

#echo "Sample: ${sample}"
#echo "Flowcell: ${flowcell}"
#echo "Lane: ${lane_id}"

gatk AddOrReplaceReadGroups \
    -I 0001.chrX.sorted.markdup.oldRG.bam \
    -O 0001.chrX.sorted.markdup.bam \
    --RGID "0001.$flowcell.$lane_id" \
    --RGSM "Buddy" \
    --RGPU "$flowcell.$lane_id" \
    --RGPL "Illumina" \
    --RGLB "SL522590" 

gatk AddOrReplaceReadGroups \
    -I 0002.chrX.sorted.markdup.oldRG.bam \
    -O 0002.chrX.sorted.markdup.bam \
    --RGID "0002.$flowcell.$lane_id" \
    --RGSM "Dandelion" \
    --RGPU "$flowcell.$lane_id" \
    --RGPL "Illumina" \
    --RGLB "SL522591" 

gatk AddOrReplaceReadGroups \
    -I 0005.chrX.sorted.markdup.oldRG.bam \
    -O 0005.chrX.sorted.markdup.bam \
    --RGID "0005.$flowcell.$lane_id" \
    --RGSM "Camelia" \
    --RGPU "$flowcell.$lane_id" \
    --RGPL "Illumina" \
    --RGLB "SL522594" 

gatk AddOrReplaceReadGroups \
    -I 0006.chrX.sorted.markdup.oldRG.bam \
    -O 0006.chrX.sorted.markdup.bam \
    --RGID "0006.$flowcell.$lane_id" \
    --RGSM "Dottie" \
    --RGPU "$flowcell.$lane_id" \
    --RGPL "Illumina" \
    --RGLB "SL522595"

samtools index 0001.chrX.sorted.markdup.bam
samtools index 0002.chrX.sorted.markdup.bam
samtools index 0005.chrX.sorted.markdup.bam
samtools index 0006.chrX.sorted.markdup.bam
