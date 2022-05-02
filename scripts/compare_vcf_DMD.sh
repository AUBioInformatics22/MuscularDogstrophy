#!/bin/sh

### Compare VCFs to find common SNPs.

##### MODULES #####
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load tabix
module load vcftools
###################

PROJDIR="/home/shared/stevison_group2"
VCFDIR="${PROJDIR}/data/vcf"
cd "${VCFDIR}"
echo "Now in $(pwd)"

#array of VCF files directory
my_vcfs=()
for i in *SNPs.filtered.select.filt.DMD.vcf.gz; do
  my_vcfs+=(${i%.vcf.gz})
done

echo "VCF files: ${my_vcfs[@]}"

# Compare sites in all 4 VCFs.
vcf-compare ${my_vcfs[0]}.vcf.gz ${my_vcfs[1]}.vcf.gz ${my_vcfs[2]}.vcf.gz ${my_vcfs[3]}.vcf.gz >intersection_DMD.4upsetR.txt

# Get lines with intersection info
grep ^VN intersection_DMD.4upsetR.txt | cut -f 2- >intersection_DMD.4upsetR.venn

# Numbers from output file used in create_figures.R to generate upset plot
