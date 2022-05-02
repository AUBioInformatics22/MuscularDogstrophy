#!/bin/bash

#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load bcftools/1.13
#module load tabix/2013-12-16

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
VCFDIR="${PROJDIR}/data/vcf"
VARDIR="${PROJDIR}/analysis/5_variants"
#######################

sample_list=(0001 0002 0005 0006)

cd "${VCFDIR}"
#echo "Now in $(pwd)"

mkdir -p stats

#echo -e "Sample ID\tUnfiltered\tFiltered\tSelected\tSelectedBcftoolsStats\tSelectedBcftoolsStatsDMD\tIntersect\tIntersectBcftoolsStats\tIntersectDMD"
echo -e "#sample\tunfiltered\tfiltered\tselected_chrX\tselected_DMD\tintersect_chrX\tintersect_DMD" >SNP_numbers.tsv

for sample in "${sample_list[@]}"; do

#bcftools stats -f '.,PASS' "$sample.chrX.sorted.SNPs.filtered.vcf.gz" >"stats/$sample.chrX.filtered.filt.vcf_stats.txt"

#bcftools stats "$sample.chrX.sorted.SNPs.filtered.select.filt.vcf.gz" >"stats/$sample.chrX.filtered.select.filt.vcf_stats.txt"

#stats for gene region of sample fsf
bcftools stats -r chrX:31264482-33330008 "$sample.chrX.sorted.SNPs.filtered.select.filt.vcf.gz" >"stats/$sample.chrX.filtered.select.filt.DMD.vcf_stats.txt"


u=$(zcat "$sample.chrX.sorted.SNPs.vcf.gz" | grep -v "^#" | wc -l)

ffb=$(grep -v "^#" "stats/$sample.chrX.filtered.filt.vcf_stats.txt" | awk '/number of SNPs/ {print $6}')

fsf=$(zcat "${sample}.chrX.sorted.SNPs.filtered.select.filt.vcf.gz" | grep -v "^#" | wc -l)
#fsfb=$(grep -v "^#" "stats/$sample.chrX.filtered.select.filt.vcf_stats.txt" | awk '/number of SNPs/ {print $6}')

fsfd=$(grep -v "^#" "stats/${sample}.chrX.filtered.select.filt.DMD.vcf_stats.txt" | awk '/number of SNPs/ {print $6}')

fsfidb=$(grep -v "^#" "stats/${sample}_isec_DMD.vcf_stats.txt" | awk '/number of SNPs/ {print $6}')
fsficb=$(grep -v "^#" "stats/${sample}_isec_chrX.vcf_stats.txt" | awk '/number of SNPs/ {print $6}')


printf "%s\t" "$sample" "$u" "$ffb" "$fsf" "$fsfd" "$fsficb" "$fsfidb" >>SNP_numbers.tsv
printf "\n" >>SNP_numbers.tsv
done

cd intersect_SNP_chrX/
#echo "Now in $(pwd)"

all_fsficb=$(zcat all_samples_SNPs_chrX.vcf.gz | grep -v "^#" | wc -l)

cd ../intersect_SNP_DMD/
#echo "Now in $(pwd)"

all_fsfidb=$(zcat all_samples_SNPs_DMD.vcf.gz | grep -v "^#" | wc -l)

cd ..
#echo "Now in $(pwd)"

printf "%s\t" "all" "---" "---" "---" "---" "$all_fsficb" "$all_fsfidb" >>SNP_numbers.tsv
#echo


#my_vcfs=()
#for i in *SNPs.filtered.select.filt.vcf.gz; do
#  my_vcfs+=(${i%.vcf.gz})
#done

#echo "VCF files: ${my_vcfs[@]}"
#echo "Done. Now getting intersection of sites using vcf-compare"
#bcftools stats "${my_vcfs[0]}.vcf.gz" "${my_vcfs[1]}.vcf.gz"
#"${my_vcfs[2]}.vcf.gz" "${my_vcfs[3]}.vcf.gz" 
#>intersection.4upsetR.txt

#vcf-compare ${my_vcfs[0]}.vcf.gz ${my_vcfs[1]}.vcf.gz ${my_vcfs[2]}.vcf.gz ${my_vcfs[3]}.vcf.gz >intersection.4upsetR.txt
#grep ^VN intersection.4upsetR.txt | cut -f 2- >intersection.4upsetR.venn

#bcftools stats -r chrX:31264482-33330008 intersect_SNP/all_samples_SNPs.vcf.gz >stats/all_SNPs.DMD.vcf_stats.txt

#echo "| Sample ID | Filtered | Selected | Intersect | bcftools stat |"
#echo "|:---------:|:---------:|:---------:|:---------:|:---------:|"
