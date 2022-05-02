#!/bin/sh

#LOAD MODULES HERE#######################################
#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load tabix
#module load vcftools
#########################################################

PROJDIR="/home/shared/stevison_group2"
VCFDIR="${PROJDIR}/data/vcf"
cd "${VCFDIR}"
echo "Now in $(pwd)"

#array of VCF files in this directory
my_vcfs=()
for i in *SNPs.filtered.select.filt.DMD.vcf.gz; do
  my_vcfs+=(${i%.vcf.gz})
done

echo "VCF files: ${my_vcfs[@]}"

#first need to remove non-PASS sites from your VCFs
#then, you need to compress and create a tabix index for each compressed VCF file

#let's use a loop!
#for i in ${my_vcfs[@]}; do

#    echo Removing non-PASS sites from $i
#    vcftools --vcf $i.vcf --remove-filtered-all --recode --recode-INFO-all --out $i.filtered

#    echo "Now compressing $i with Bgzip"
#    bgzip -f $i.filtered.recode.vcf

#    echo Now building a tabix index for $i
#    tabix -p vcf $i.filtered.recode.vcf.gz
#done

#echo "Done. Now getting intersection of sites using vcf-compare"
vcf-compare ${my_vcfs[0]}.vcf.gz ${my_vcfs[1]}.vcf.gz ${my_vcfs[2]}.vcf.gz ${my_vcfs[3]}.vcf.gz >intersection_DMD.4upsetR.txt

grep ^VN intersection_DMD.4upsetR.txt | cut -f 2- >intersection_DMD.4upsetR.venn

#paste numbers from resulting file into R script on Canvas. Run R script on your own computer.
