#!/bin/bash

##### MODULES #####
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bcftools/1.13
module load tabix/2013-12-16
#module load gatk/4.1.0.0
#source activate gatk
###################

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
VCFDIR="${PROJDIR}/data/vcf"
VARDIR="${PROJDIR}/analysis/5_variants"
SCRATCHDIR="/scratch/MD_data"
#######################

#you can use -L on HaplotypeCaller to exclude regions

#goals:

#filter for alternate alleles autozygous in all male samples (1-4)
#--sample-name 0001 \
#--sample-name 0002 \
#--selectExpressions 

#filter for alternate alleles heterozygous in all female samples (5-8)
#--sample-name 0005 \
#--sample-name 0006 \

#gatk VariantFiltration \     #MALE
#  -V 0001.chrX.sorted.SNPs.filtered.vcf \
#  -O 0001.chrX.sorted.SNPs.filtered.select.vcf \
#  --genotype-filter-expression "isHomVar==1" \
#  --genotype-filter-name "isHomVarFilter" 


#gatk VariantFiltration \     #FEMALE
#  -V 0005...vcf \
#  -O 0005...filtered.vcf \
#  --genotype-filter-expression "isHet==1" \
#  --genotype-filter-name "isHetFilter"


#gatk SelectVariant \
#  -V 0001.chrX.sorted.SNPs.filtered.select.vcf \
#  -O 0001.chrX.sorted.SNPs.filtered.select.only.vcf
#  --set-filtered-gt-to-nocall


#goal:
#find sites that are shared between ALL male samples
#find sites that are shared between ALL female samples
#find sites shared between male and female samples

### VCFtools
### BCFtools

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
  echo "after m snp"

#  bcftools view -i 'GT="alt"' \
#  -o "${sample}.chrX.sorted.INDELs.filtered.select.vcf" \
#  "${sample}.chrX.sorted.INDELs.filtered.vcf.gz" 
#  echo "after m indel"

  bgzip -f "${sample}.chrX.sorted.SNPs.filtered.select.vcf"
#  bgzip -f "${sample}.chrX.sorted.INDELs.filtered.select.vcf"

  tabix -p vcf "${sample}.chrX.sorted.SNPs.filtered.select.vcf.gz"
#  tabix -p vcf "${sample}.chrX.sorted.INDELs.filtered.select.vcf.gz"
done

echo "Done with males"


for sample in "${female_samples[@]}"; do
  bcftools view -i 'GT="het"' \
  -o "${sample}.chrX.sorted.SNPs.filtered.select.vcf" \
  "${sample}.chrX.sorted.SNPs.filtered.vcf.gz" 
  echo "after f snp"

#  bcftools view -i 'GT="het"' \
#  -o "${sample}.chrX.sorted.INDELs.filtered.select.vcf" \
#  "${sample}.chrX.sorted.INDELs.filtered.vcf.gz" 
#  echo "after f indel"

  bgzip -f "${sample}.chrX.sorted.SNPs.filtered.select.vcf"
#  bgzip -f "${sample}.chrX.sorted.INDELs.filtered.select.vcf"

  tabix -p vcf "${sample}.chrX.sorted.SNPs.filtered.select.vcf.gz"
#  tabix -p vcf "${sample}.chrX.sorted.INDELs.filtered.select.vcf.gz"
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

bgzip -f 0000.vcf
bgzip -f 0001.vcf
bgzip -f 0002.vcf
bgzip -f 0003.vcf

tabix -p vcf "0000.vcf.gz"
tabix -p vcf "0001.vcf.gz"
tabix -p vcf "0002.vcf.gz"
tabix -p vcf "0003.vcf.gz"

bcftools merge -f .,PASS 0000.vcf.gz 0001.vcf.gz 0002.vcf.gz 0003.vcf.gz >all_samples_SNPs.vcf
bgzip -f all_samples_SNPs.vcf
tabix -p vcf "all_samples_SNPs.vcf.gz"


#cd ../intersect_INDEL
#echo "Now in $(pwd)"  

#bcftools merge -f .,PASS 0000.vcf.gz 0001.vcf.gz 0002.vcf.gz 0003.vcf.gz >all_samples_INDELs.vcf
#bgzip -f all_samples_INDELs.vcf
#tabix -p vcf "all_samples_INDELs.vcf.gz"


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


#goal: combine vcfs

#goal: analyze variants

### SnpEff





#LOAD MODULES HERE#######################################
#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load tabix
#module load vcftools
#########################################################

##ASSIGN VARIABLES#######################################
#array of VCF files in this directory, there should be 3 - the RAW SNPs, the GATK filtered SNPs, and the Adjusted filtered SNPs
#my_vcfs=()
#for i in *.vcf; do
#  my_vcfs+=(${i%.vcf})
#done

#echo VCF files: ${my_vcfs[@]}
##SAMN03152093.merged.sorted.SNPs.filtered2 SAMN03152093.merged.sorted.SNPs.filtered SAMN03152093.merged.sorted.SNPs


#CODE:#######################################

#first need to remove non-PASS sites from your VCFs
#then, you need to compress and create a tabix index for each compressed VCF file

#let's use a loop!
#for i in ${my_vcfs[@]}
#do

#    echo Removing non-PASS sites from $i
#    vcftools --vcf $i.vcf --remove-filtered-all --recode --recode-INFO-all --out $i.filtered

#    echo "Now compressing $i with Bgzip"
#    bgzip -f $i.filtered.recode.vcf

#    echo Now building a tabix index for $i
#    tabix -p vcf $i.filtered.recode.vcf.gz
#done

#echo "Done. Now getting intersection of sites using vcf-compare"
#vcf-compare ${my_vcfs[0]}.filtered.recode.vcf.gz ${my_vcfs[1]}.filtered.recode.vcf.gz ${my_vcfs[2]}.filtered.recode.vcf.gz >intersection.4upsetR.txt

#grep ^VN intersection.4upsetR.txt | cut -f 2- >intersection.4upsetR.venn

#paste numbers from resulting file into R script on Canvas. Run R script on your own computer.
