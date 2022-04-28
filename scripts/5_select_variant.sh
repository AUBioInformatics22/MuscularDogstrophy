#!/bin/bash

### Select variants that are hemizygous for an alternate allele in the males.
### Select variants that are heterozygous for an alternate allele in the females.
### Intersect the selected variants and find positions that are in all 4 samples.

# Note: Code handling INDELs is commented out.

##### MODULES #####
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bcftools/1.13
module load tabix/2013-12-16
###################

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
VCFDIR="${PROJDIR}/data/vcf"
VARDIR="${PROJDIR}/analysis/5_variants"
#######################

male_samples=(0001 0002)
female_samples=(0005 0006)
samples=(0001 0002 0005 0006)

cd "${VCFDIR}"
echo "Now in $(pwd)"

echo "Male samples: ${male_samples[@]}"
echo "Female samples: ${female_samples[@]}"


for sample in "${male_samples[@]}"; do
  bcftools view -i 'GT="alt"' \
  -o "${sample}.chrX.sorted.SNPs.filtered.select.vcf" \
  "${sample}.chrX.sorted.SNPs.filtered.vcf.gz" 
  
  bgzip -f "${sample}.chrX.sorted.SNPs.filtered.select.vcf"
  tabix -p vcf "${sample}.chrX.sorted.SNPs.filtered.select.vcf.gz"
  
  echo "Filtered SNPs for sample: ${sample} (M)"

#  bcftools view -i 'GT="alt"' \
#  -o "${sample}.chrX.sorted.INDELs.filtered.select.vcf" \
#  "${sample}.chrX.sorted.INDELs.filtered.vcf.gz" 

#  bgzip -f "${sample}.chrX.sorted.INDELs.filtered.select.vcf"
#  tabix -p vcf "${sample}.chrX.sorted.INDELs.filtered.select.vcf.gz"

# echo "Filtered INDELs for sample: ${sample} (M)"
done

echo "Done with males"


for sample in "${female_samples[@]}"; do
  bcftools view -i 'GT="het"' \
  -o "${sample}.chrX.sorted.SNPs.filtered.select.vcf" \
  "${sample}.chrX.sorted.SNPs.filtered.vcf.gz" 
  
  bgzip -f "${sample}.chrX.sorted.SNPs.filtered.select.vcf"
  tabix -p vcf "${sample}.chrX.sorted.SNPs.filtered.select.vcf.gz"
  
  echo "Filtered SNPs for sample: ${sample} (F)"

#  bcftools view -i 'GT="het"' \
#  -o "${sample}.chrX.sorted.INDELs.filtered.select.vcf" \
#  "${sample}.chrX.sorted.INDELs.filtered.vcf.gz" 

#  bgzip -f "${sample}.chrX.sorted.INDELs.filtered.select.vcf"
#  tabix -p vcf "${sample}.chrX.sorted.INDELs.filtered.select.vcf.gz"

#  echo "Filtered INDELs for sample: ${sample} (F)"
done

echo "Done with females"


echo "Comparing SNPs..."
bcftools isec \
  -n =4 -p intersect_SNP \
  0001.chrX.sorted.SNPs.filtered.select.vcf.gz \
  0002.chrX.sorted.SNPs.filtered.select.vcf.gz \
  0005.chrX.sorted.SNPs.filtered.select.vcf.gz \
  0006.chrX.sorted.SNPs.filtered.select.vcf.gz 

#echo "Comparing INDELs..."
#bcftools isec \
#  -n =4 -p intersect_INDEL \
#  0001.chrX.sorted.INDELs.filtered.select.vcf.gz \
#  0002.chrX.sorted.INDELs.filtered.select.vcf.gz \
#  0005.chrX.sorted.INDELs.filtered.select.vcf.gz \
#  0006.chrX.sorted.INDELs.filtered.select.vcf.gz 

cd intersect_SNP
echo "Now in $(pwd)"  

echo "Zipping SNP intersection files."
bgzip -f 0000.vcf
bgzip -f 0001.vcf
bgzip -f 0002.vcf
bgzip -f 0003.vcf

echo "Indexing SNP intersection files."
tabix -p vcf "0000.vcf.gz"
tabix -p vcf "0001.vcf.gz"
tabix -p vcf "0002.vcf.gz"
tabix -p vcf "0003.vcf.gz"

echo "Merging SNP intersection files."
bcftools merge -f .,PASS 0000.vcf.gz 0001.vcf.gz 0002.vcf.gz 0003.vcf.gz >all_samples_SNPs.vcf
bgzip -f all_samples_SNPs.vcf
tabix -p vcf "all_samples_SNPs.vcf.gz"


#cd ../intersect_INDEL
#echo "Now in $(pwd)"  

#bcftools merge -f .,PASS 0000.vcf.gz 0001.vcf.gz 0002.vcf.gz 0003.vcf.gz >all_samples_INDELs.vcf
#bgzip -f all_samples_INDELs.vcf
#tabix -p vcf "all_samples_INDELs.vcf.gz"

#bgzip -f 0000.vcf
#bgzip -f 0001.vcf
#bgzip -f 0002.vcf
#bgzip -f 0003.vcf

#tabix -p vcf "0000.vcf.gz"
#tabix -p vcf "0001.vcf.gz"
#tabix -p vcf "0002.vcf.gz"
#tabix -p vcf "0003.vcf.gz"

# Comparisons not needed for current analysis, but could be interesting to look at later.
#echo "Comparing male SNPs..."
#bcftools isec \
#  0001.chrX.sorted.SNPs.filtered.select.vcf.gz \
#  0002.chrX.sorted.SNPs.filtered.select.vcf.gz \
#  -n =2 -p intersect

#echo "Comparing male INDELs..."
#bcftools isec \
#  0001.chrX.sorted.INDELs.filtered.select.vcf.gz \
#  0002.chrX.sorted.INDELs.filtered.select.vcf.gz \
#  -n =2 -p intersect

#echo "Comparing female SNPs..."
#bcftools isec \
#  0005.chrX.sorted.SNPs.filtered.select.vcf.gz \
#  0006.chrX.sorted.SNPs.filtered.select.vcf.gz \
#  -n =2 -p intersect
  
#echo "Comparing female INDELs..."
#bcftools isec \
#  0005.chrX.sorted.INDELs.filtered.select.vcf.gz \
#  0006.chrX.sorted.INDELs.filtered.select.vcf.gz \
#  -n =2 -p intersect

echo "all done!"
