#!/bin/bash

#Compare the coverage when divided by genome size vs. by chrX size.

###### CHECK USAGE #####
#if [ $# != 1 ]; then echo "$0: Usage: Please include a sample ID."; exit; fi
########################

#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load samtools/1.3.1

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
DATADIR="${PROJDIR}/data/fastq"
BAMDIR="${PROJDIR}/data/bam"
CHRXDIR="${PROJDIR}/data/chrX_bam"
MARKDIR="${PROJDIR}/data/marked_chrX_bam"
INDEXDIR="${PROJDIR}/analysis/0_index_genome"
ALIGNDIR="${PROJDIR}/analysis/3_alignment"
DUPDIR="${PROJDIR}/analysis/4_duplicates"
SCRATCHDIR="/scratch/MD_data"
#######################

##### VARIABLES #####
#sample="$1"
mapfile -t sample_list < <(echo -e "0001\n0002\n0005\n0006")
ref="canFam6"
chrX_size="$(awk '/chrX/ {print $2}' ${INDEXDIR}/${ref}.masked.fa.fai)"
#genome_size="$(awk '{sum+=$2} END {print sum}' ${INDEXDIR}/${ref}.masked.fa.fai)"
#####################

#checks
#echo "sample: ${sample}"
#echo "genome size: ${genome_size}"
#echo "chrX size: ${chrX_size}"


cd "${MARKDIR}"
echo "now in $(pwd)"


#if [ ! -s "${sample}.chrX.depth.tmp" ]; then
#  source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#  module load samtools/1.3.1
#  echo $(date +%F%t%T) "Running samtools depth..."
#  samtools depth -a ${sample}.chrX.sorted.markdup.bam >${sample}.chrX.marked.depth.tmp
#fi


#collect QC metrics of files

#if [ ! -s "${sample}.chrX.sums.tmp" ]; then
#  echo "Calculating sums..."
#  mapfile -t sum_sumsq < <(awk '{sum+=$3; sumsq+=$3*$3} END {print sum; print sumsq}' ${sample}.chrX.marked.depth.tmp)
#  echo -e "${sum_sumsq[0]}\n${sum_sumsq[1]}" >"${sample}.chrX.marked.sums.tmp"
#else
#  mapfile -t sum_sumsq <"${sample}.chrX.marked.sums.tmp"
#fi

#checks
#echo "${sum_sumsq[@]}"
#sum="${sum_sumsq[0]}"
#sumsq="${sum_sumsq[1]}"
#echo $sum
#echo $sumsq

#awk -v c="$chrX_size" -v sum="${sum_sumsq[0]}" -v sumsq="${sum_sumsq[1]}" \
#  'END { print "Average = ",sum/c; print "Stdev = ",sqrt(sumsq/c - (sum/c)**2) }' \
#  >"${DUPDIR}/${sample}.chrX.marked.coverage_summary.txt" < <(echo) 

for sample in "${sample_list[@]}"; do
  echo $(date +%F%t%T) "Running samtools flagstat for sample $sample..."
  samtools flagstat ${sample}.chrX.sorted.markdup.bam \
    >${DUPDIR}/${sample}.chrX.marked.flagstat_output.txt
done

echo "all done!"

### NOTE: after all sequences are done, run stat_csv.sh to create csv files for graphing in R. 
