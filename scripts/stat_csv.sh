#!/bin/bash

### Determine coverage (raw, aligned, with duplicates marked, and variants),
### percent mapped, and percent duplicates for given sequences.
### Create csv files that can be used to make graphs in R (create_figures.R)
### Coverage  = (number of reads * read length) / genome size

### The following CSV files will be output in the directory containing
#### the data they are gathering and copied to the designated CSVDIR.
# all_metrics.csv -> Input for create_figures.R
# raw_coverage_data.csv
# aligned_chrX_coverage_data.csv
# marked_chrX_coverage_data.csv
# variants_coverage_data.csv
# percent_mapped_data.csv
# percent_duplicates_data.csv

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
INDEXDIR="${PROJDIR}/analysis/0_index_genome"
RAWDIR="${PROJDIR}/analysis/1_fastqc_raw"
ALIGNDIR="${PROJDIR}/analysis/3_alignment"
DUPDIR="${PROJDIR}/analysis/4_duplicates"
VARDIR="${PROJDIR}/analysis/5_variants"
CSVDIR="${PROJDIR}/analysis/csv"
VCFDIR="${PROJDIR}/data/vcf"
#######################

ref="canFam6"
mapfile -t sample_list < <(echo -e "0001\n0002\n0005\n0006")
output_file="coverage_data.csv"

genome_size="$(awk '{sum+=$2} END {print sum}' ${INDEXDIR}/${ref}.masked.fa.fai)"
chrX_size="$(awk '/chrX/ {print $2}' ${INDEXDIR}/${ref}.masked.fa.fai)"

# Initialize arrays for each metric.
declare -a raw_coverage_list=()
declare -a est_raw_chrX_coverage_list=()
declare -a aligned_chrX_coverage_list=()
declare -a marked_chrX_coverage_list=()
declare -a perc_map_wg_list=()
declare -a perc_map_chrX_list=()
declare -a perc_dup_list=()
declare -a var_SNP_coverage_list=()
declare -a var_INDEL_coverage_list=()

# checks
echo "Samples: ${sample_list[@]}"
echo "Genome size: $genome_size"
echo "chrX size: $chrX_size"


### RAW

cd "${RAWDIR}"
echo "id,raw,est_raw_chrX" >"raw_${output_file}"

#########
# in fastqc_data.txt :
# Total Sequences num_reads
# Sequence length seq_len
#########

for sample in ${sample_list[@]}; do
  unzip -p H2LCVDSX3_s1_1_*${sample}_fastqc.zip H2LCVDSX3_s1_1_*${sample}_fastqc/fastqc_data.txt \
    | head >forward.tmp
  f_num_reads=$(awk '/Total Sequences/ {print $3}' forward.tmp)
  f_seq_len=$(awk '/Sequence length/ {print $3}' forward.tmp)
  f_coverage=$(awk '{printf "%.3f", $1*$2/$3}' <<<"$f_num_reads $f_seq_len $genome_size")

  unzip -p H2LCVDSX3_s1_2_*${sample}_fastqc.zip H2LCVDSX3_s1_2_*${sample}_fastqc/fastqc_data.txt \
    | head >reverse.tmp
  r_num_reads=$(awk '/Total Sequences/ {print $3}' reverse.tmp)
  r_seq_len=$(awk '/Sequence length/ {print $3}' reverse.tmp)
  r_coverage=$(awk '{printf "%.3f", $1*$2/$3}' <<<"$r_num_reads $r_seq_len $genome_size")

  est_chrX_coverage=$(awk '{printf "%.3f", $1*0.04705*$2/$3}' <<<"$f_num_reads $f_seq_len $chrX_size")
  #est_chrX_coverage=$(echo "${f_coverage}" | awk '{print $1*0.04705}')
  est_raw_chrX_coverage_list+=(${est_chrX_coverage})

  if [[ $f_coverage == $r_coverage ]]; then
    raw_coverage_list+=(${f_coverage})
    echo "${sample},${f_coverage},${est_chrX_coverage}" >>"raw_${output_file}"
  else
    echo "Raw coverage: Sample ${sample} forward and reverse coverage do not match. Exiting..."
    exit
  fi
done

# Remove temporary files and copy output file to CSV directory.
rm *.tmp
cp "raw_${output_file}" "${CSVDIR}"


### ALIGNED CHRX AND PERCENT MAPPED

cd "${ALIGNDIR}"
echo "id,aligned_chrX" >"aligned_chrX_${output_file}"
echo "id,percent_mapped_wg,percent_mapped_chrX" >percent_mapped_data.csv

#########
#in coverage_summary.txt :
#  Average = coverage
#in flagstat_output.txt :
#  /\d+/ + 0 mapped (perc_map : N/A)
#########

for sample in ${sample_list[@]}; do
  coverage=$(awk '/Average/ {print $3}' "${sample}.chrX.coverage_summary.txt")
  aligned_chrX_coverage_list+=(${coverage})
  perc_map_wg=$(awk -F "[(|%"] '/[[:digit:]]+ \+ [[:digit:]]+ mapped / {print $2}' "${sample}.flagstat_output.txt")
  perc_map_wg_list+=($perc_map_wg)
  perc_map_chrX=$(awk -F "[(|%"] '/[[:digit:]]+ \+ [[:digit:]]+ mapped / {print $2}' "${sample}.chrX.flagstat_output.txt")
  perc_map_chrX_list+=($perc_map_chrX)
  echo "${sample},${coverage}" >>"aligned_chrX_${output_file}"
  echo "${sample},${perc_map_wg},${perc_map_chrX}" >>"percent_mapped_data.csv"
done

# Copy output files to CSV directory.
cp "aligned_chrX_${output_file}" "${CSVDIR}"
cp percent_mapped_data.csv "${CSVDIR}"


### MARKED CHRX AND PERCENT DUPLICATED

cd "${DUPDIR}"
echo "id,marked_chrX" >"marked_chrX_${output_file}"
echo "id,percent_duplicates" >"percent_duplicates_data.csv"

for sample in "${sample_list[@]}"; do
  coverage=$(awk '/Average/ {print $3}' "${sample}.chrX.marked.coverage_summary.txt")
  marked_chrX_coverage_list+=(${coverage})
  perc_dup=$(awk 'NR == 8 {print $9*100}' "${sample}.chrX.sorted.dup_metrics")
  perc_dup_list+=($perc_dup)
  echo "${sample},${coverage}" >>"marked_chrX_${output_file}"
  echo "${sample},${perc_dup}" >>"percent_duplicates_data.csv"
done

# Copy output files to CSV directory.
cp "marked_chrX_${output_file}" "${CSVDIR}"
cp percent_duplicates_data.csv "${CSVDIR}"


### VARIANTS

cd "${VCFDIR}"
echo "id,var_SNP,var_INDEL" >"variants_${output_file}"

for sample in "${sample_list[@]}"; do
  SNP_coverage=$(awk 'NR>1 {print $3}' "${sample}.chrX.sorted.SNPs.idepth")
  var_SNP_coverage_list+=(${SNP_coverage})
  INDEL_coverage=$(awk 'NR>1 {print $3}' "${sample}.chrX.sorted.INDELs.idepth")
  var_INDEL_coverage_list+=(${INDEL_coverage})
  echo "${sample},${SNP_coverage},${INDEL_coverage}" >>"variants_${output_file}"
done

# Copy output file to CSV directory.
cp "variants_${output_file}" "${CSVDIR}"


### ALL

cd "${CSVDIR}"
echo "id,raw,est_raw_chrX,aligned_chrX,marked_chrX,var_SNP,var_INDEL,percent_mapped_wg,percent_mapped_chrX,percent_duplicates" >"all_metrics.csv"

for (( i=0; i<=$((${#sample_list[@]}-1)); i++ )); do
  echo "${sample_list[i]},${raw_coverage_list[i]},${est_raw_chrX_coverage_list[i]},${aligned_chrX_coverage_list[i]},${marked_chrX_coverage_list[i]},${var_SNP_coverage_list[i]},${var_INDEL_coverage_list[i]},${perc_map_wg_list[i]},${perc_map_chrX_list[i]},${perc_dup_list[i]}" \
    >>"${CSVDIR}/all_metrics.csv"
done

