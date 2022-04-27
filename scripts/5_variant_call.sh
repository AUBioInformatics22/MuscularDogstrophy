#!/bin/bash 

### Generate VCF files of SNPs and INDELs for each sample.
### Perform hard filtering based on standard quality parameters.
### Assess coverage.

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load gatk/4.1.0.0
module load java/1.8.0_192
source activate gatk

##### DIRECTORIES #####
PROJDIR="/home/shared/stevison_group2"
MARKDIR="${PROJDIR}/data/marked_chrX_bam"
VCFDIR="${PROJDIR}/data/vcf"
INDEXDIR="${PROJDIR}/analysis/0_index_genome"
VARDIR="${PROJDIR}/analysis/5_variant"
SCRATCHDIR="/scratch/MD_data"
#######################

mapfile -t male_samples < <(echo -e "0001\n0002")
mapfile -t female_samples < <(echo -e "0005\n0006")
mapfile -t sample_list < <(echo -e "0001\n0002\n0005\n0006")
ref="${INDEXDIR}/canFam6.masked.fa"

#### COMMAND LINE FOR GATK BELOW ####

for samplenum in "${sample_list[@]}"; do
  sample="$samplenum.chrX.sorted"
  
  cd "${MARKDIR}"
  echo "now in $(pwd)"
  
  module load gcc/7.2.0

  #Build Bam Index
  #gatk --java-options "-Xmx1G" BuildBamIndex -I "${sample}.markdup.bam" -R "$ref"

  #Short Variant Discovery; this step can be multithreaded if you request more cores (try 4 first, then up to 8 if needed)
  #male
  if [[ "$samplenum" == "0001" || "$samplenum" == "0002" ]]; then
    echo "Sample $samplenum is a male"
    gatk --java-options "-Xmx8G" HaplotypeCaller \
      -R "$ref" \
      -I "${sample}.markdup.bam" \
      --sample-ploidy 1 \
      -L chrX \
      -O "${VCFDIR}/$sample.g.vcf.gz" \
      -ERC GVCF
  fi

  #female
  if [[ "$samplenum" == "0005" || "$samplenum" == "0006" ]]; then
    echo "Sample $samplenum is a female"
    gatk --java-options "-Xmx8G" HaplotypeCaller \
      -R "$ref" \
      -I "${sample}.markdup.bam" \
      --sample-ploidy 2 \
      -L chrX \
      -O "${VCFDIR}/$sample.g.vcf.gz" \
      -ERC GVCF
  fi

  cd "${VCFDIR}"
  echo "now in $(pwd)"

  #Joint Genotyping
  gatk --java-options "-Xmx8G" GenotypeGVCFs -R "$ref" -V "$sample.g.vcf.gz" -O "$sample.vcf.gz"

  #Extract only SNPs
  gatk SelectVariants -R "$ref" --variant "$sample.vcf.gz" --select-type-to-include SNP --output "${VCFDIR}/$sample.SNPs.vcf"

  gatk SelectVariants -R "$ref" --variant "$sample.vcf.gz" --select-type-to-include INDEL --output "${VCFDIR}/$sample.INDELs.vcf"

  #Perform variant filtering
  gatk VariantFiltration -R "$ref" --variant "$sample.SNPs.vcf" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --output "$sample.SNPs.filtered.vcf"

#  gatk VariantFiltration -R "$ref" --variant "$sample.INDELs.vcf" \
#    --filter-expression "QD < 2.0" --filter-name "QD2" \
#    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
#    --filter-expression "FS > 60.0" --filter-name "FS60" \
#    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
#    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#    --output "$sample.INDELs.filtered.vcf"

  gatk VariantFiltration \
    -R "$ref" \
    -V "$sample.INDELs.vcf.gz" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "$sample.INDELs.filtered.vcf.gz"

  gzip "$sample.SNPs.vcf" "$sample.SNPs.filtered.vcf"
  gzip "$sample.INDELs.vcf" "$sample.INDELs.filtered.vcf"

  #get depth statistics
  module load vcftools
  vcftools --gzvcf "$sample.SNPs.filtered.vcf.gz" --depth --out "$sample.SNPs" #output will be called $sample.SNPs.idepth
  vcftools --gzvcf "$sample.INDELSs.filtered.vcf.gz" --depth --out "$sample.INDELs" #output will be called $sample.INDELs.idepth

  #print QC metrics to a file
  snp_depth=$(awk 'NR>1 {print $3}' "$sample.SNPs.idepth")
  indel_depth=$(awk 'NR>1 {print $3}' "$sample.INDELs.idepth")
  
  echo SNP depth for $sample is $snp_depth
  echo INDEL depth for $sample is $indel_depth
done
