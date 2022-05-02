#!/bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bcftools/1.13
module load tabix/2013-12-16

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
VCFDIR="${PROJDIR}/data/vcf"
VARDIR="${PROJDIR}/analysis/5_variants"
#######################

sample_list=(0001 0002 0005 0006)

cd "${VCFDIR}"
#echo "Now in $(pwd)"

mkdir -p stats

echo -e "#sample\tunfiltered\tfiltered\tselected_chrX\tselected_DMD\tintersect_chrX\tintersect_DMD" >SNP_numbers.tsv

for sample in "${sample_list[@]}"; do
  #stats for sites that pass filtering
  bcftools stats -f '.,PASS' "$sample.chrX.sorted.SNPs.filtered.vcf.gz" >"stats/$sample.chrX.filtered.filt.vcf_stats.txt"
  
  #stats for selected sites
  bcftools stats "$sample.chrX.sorted.SNPs.filtered.select.filt.vcf.gz" >"stats/$sample.chrX.filtered.select.filt.vcf_stats.txt"

  #stats for selected sites of gene region
  bcftools stats -r chrX:31264482-33330008 "$sample.chrX.sorted.SNPs.filtered.select.filt.vcf.gz" >"stats/$sample.chrX.filtered.select.filt.DMD.vcf_stats.txt"

  #unfiltered
  u=$(zcat "$sample.chrX.sorted.SNPs.vcf.gz" | grep -v "^#" | wc -l)

  #filtered
  ffb=$(grep -v "^#" "stats/$sample.chrX.filtered.filt.vcf_stats.txt" | awk '/number of SNPs/ {print $6}')

  #selected
  fsf=$(zcat "${sample}.chrX.sorted.SNPs.filtered.select.filt.vcf.gz" | grep -v "^#" | wc -l)

  #selected in DMD
  fsfd=$(grep -v "^#" "stats/${sample}.chrX.filtered.select.filt.DMD.vcf_stats.txt" | awk '/number of SNPs/ {print $6}')

  #intersection in DMD
  fsfidb=$(grep -v "^#" "stats/${sample}_isec_DMD.vcf_stats.txt" | awk '/number of SNPs/ {print $6}')
  
  #intersection in chrX
  fsficb=$(grep -v "^#" "stats/${sample}_isec_chrX.vcf_stats.txt" | awk '/number of SNPs/ {print $6}')
  
  #print values for sample to output file
  printf "%s\t" "$sample" "$u" "$ffb" "$fsf" "$fsfd" "$fsficb" "$fsfidb" >>SNP_numbers.tsv
  printf "\n" >>SNP_numbers.tsv
done

cd intersect_SNP_chrX/
#get number of common SNPs (chrX)
all_fsficb=$(zcat all_samples_SNPs_chrX.vcf.gz | grep -v "^#" | wc -l)

cd ../intersect_SNP_DMD/
#get number of common SNPs (DMD)
all_fsfidb=$(zcat all_samples_SNPs_DMD.vcf.gz | grep -v "^#" | wc -l)

cd ..
#print overall values to output file as a check
printf "%s\t" "all" "---" "---" "---" "---" "$all_fsficb" "$all_fsfidb" >>SNP_numbers.tsv
