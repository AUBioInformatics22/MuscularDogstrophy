#!/bin/bash

### Replace the read groups for each sample based on unique metadata.

##### MODULES #####
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load gatk/4.1.0.0
module load samtools/1.3.1
###################

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
MARKDIR="${PROJDIR}/data/marked_chrX_bam"
#######################

##### VARIABLES #####
flowcell="H2LCVDSX3"
lane_id="3"
#####################

###################
# Read Group Header
# ID = sample_name.$flowcell.$lane_id
# SM (sample) = sample_name
# PU (platform unit) = $flowcell.$lane_id
# PL (platform) = illumina
# LB (library) = $flowcell.$lane_id
###################

# Move to the directory the files are in.
cd "${MARKDIR}"
echo "Now in $(pwd)"

# Replace read groups for each sample individually.
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

# index the files with the new read groups
samtools index 0001.chrX.sorted.markdup.bam
samtools index 0002.chrX.sorted.markdup.bam
samtools index 0005.chrX.sorted.markdup.bam
samtools index 0006.chrX.sorted.markdup.bam
