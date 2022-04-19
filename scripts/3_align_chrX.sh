#!/bin/bash

#Use bwa mem to create read groups and align reads to a reference genome
#Subset chrX.
#Collect quality metrics.

##### QUEUE PARAMETERS #####
# medium
# 8 cores
# 16gb
############################


###### CHECK USAGE #####
if [ $# != 1 ]; then echo "$0: Usage: Please include a sample ID."; exit; fi
########################

##### MODULES #####
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bwa/0.7.12
module load samtools/1.3.1
###################

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
DATADIR="${PROJDIR}/data/fastq"
BAMDIR="${PROJDIR}/data/bam"
INDEXDIR="${PROJDIR}/analysis/0_index_genome"
CHRXDIR="${PROJDIR}/data/chrX_bam"
ALIGNDIR="${PROJDIR}/analysis/3_alignment"
SCRATCHDIR="/scratch/MD_data"
#######################

##### VARIABLES #####
sample=$1
ref="canFam6"
forward="${DATADIR}/H2LCVDSX3_s1_1_IDT8_UDI_*_6797-RN-${sample}.fastq.gz"
reverse="${DATADIR}/H2LCVDSX3_s1_2_IDT8_UDI_*_6797-RN-${sample}.fastq.gz"
flowcell=$(zcat ${DATADIR}/H2LCVDSX3_s1_1_IDT8_UDI_*_i7-IDT8_UDI_*_i5_6797-RN-${sample}.fastq.gz | head -1 | awk 'BEGIN {FS = ":"} {print $3}')
lane_id=$(zcat ${DATADIR}/H2LCVDSX3_s1_1_IDT8_UDI_*_i7-IDT8_UDI_*_i5_6797-RN-${sample}.fastq.gz | head -1 | awk 'BEGIN {FS = ":"} {print $4}')
genome_size="$(awk '{sum+=$2} END {print sum}' ${INDEXDIR}/${ref}.masked.fa.fai)"
#####################

#checks
echo "sample: ${sample}"
echo "forward file:" ${forward}
if [ ! -e ${forward} ]; then echo "$0: Error: Forward file ${forward} does not exist."; exit; fi
echo "reverse file:" ${reverse}
if [ ! -e ${reverse} ]; then echo "$0: Error: Forward file ${reverse} does not exist."; exit; fi
echo "flowcell: ${flowcell}"
echo "lane: ${lane_id}"
echo "genome size: ${genome_size}"

#move to bam data directory
cd "${BAMDIR}"
echo "now in $(pwd)"

#TODO: Add a check for index files.
#ls "${INDEXDIR}/*"
#echo "copy index files to working directory"
#cp -v "${INDEXDIR}/${ref}."* .

echo
echo "Running alignment..."

###################
# Read Group Header
# ID = $sample.$flowcell.$lane_id
# SM (sample) = $sample
# PU (platform unit) = $flowcell.$lane_id
# PL (platform) = illumina
# LB (library) = $flowcell.$lane_id
###################

echo $(date +%F%t%T) "Running bwa mem to create sam..."
bwa mem -M -v 2 -t 8 -R "@RG\tID:$sample.$flowcell.$lane_id\tSM:$sample\tPU:$flowcell.$lane_id\tPL:Illumina\tLB:$flowcell.$lane_id" $ref $forward $reverse \
  >${SCRATCHDIR}/$sample.sam

cd "${SCRATCHDIR}"
echo "now in $(pwd)"

echo $(date +%F%t%T) "Running samtools view to convert sam to bam..."
samtools view -Sb -@ 8 $sample.sam >$sample.bam

echo $(date +%F%t%T) "Running samtools sort to sort bam..."
samtools sort -@ 8 -m 1500MB $sample.bam >$sample.sorted.bam

echo $(date +%F%t%T) "Running samtools index on sorted bam..."
samtools index $sample.sorted.bam

#move to bam data directory <-DO NOT UNCOMMENT THIS BLOCK
#cd "${BAMDIR}"
#echo "now in $(pwd)"

##################
#region: chrX
#input file: 0001.sorted.bam
#output file: 0001.chrX.sorted.bam
##################

echo 
echo "Extracting chrX from sample ${sample}..."

echo $(date +%F%t%T) "Running samtools view..."
samtools view -b -@ 8 $sample.sorted.bam chrX >"${CHRXDIR}/${sample}.chrX.sorted.bam"

cd "${CHRXDIR}"
echo "now in $(pwd)"

echo $(date +%F%t%T) "Running samtools index on chrX..."
samtools index ${sample}.chrX.sorted.bam


#collect QC metrics of files
echo $(date +%F%t%T) "Running samtools depth..."
samtools depth -a ${sample}.chrX.sorted.bam \
  | awk "{sum+=\$3; sumsq+=\$3*\$3} END { print \"Average = \",sum/$chrX_size; print \"Stdev = \",sqrt(sumsq/$chrX_size - (sum/$chrX_size)**2)}" \
  >${ALIGNDIR}/${sample}.chrX.coverage_summary.txt

echo $(date +%F%t%T) "Running samtools flagstat..."
samtools flagstat ${sample}.chrX.sorted.bam \
  >${ALIGNDIR}/${sample}.chrX.flagstat_output.txt

echo "all done!"

### NOTE: after all sequences are done, run stat_csv.sh to create csv files for graphing in R. 
