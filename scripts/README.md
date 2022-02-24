## scripts

```
# This script requires a list of line delimited SRR IDs called "SRR_Acc_List.txt". It calls the following scripts:
# 1_index_genome.sh
# 2_sra_fastq.sh
# 2a_qc.sh
# 3_split_lanes.sh
# 4_BWA_alignment.sh


##### INITIALIZING #####
#making directories if they don't already exist
mkdir -p analysis/index_genome
mkdir -p analysis/fastqc
mkdir -p analysis/alignment
mkdir -p reports/validation_outputs
mkdir -p reports/validation_errors
mkdir -p reports/fastqdump_outputs

#creating an array of sample SRR IDs
mapfile -t Acc_List < SRR_Acc_List.txt

########################

#./1a_sra_fastq09.sh
#./1b_sra_fastq10.sh
#./1c_sra_fastq14.sh
#./1d_sra_fastq34.sh

#./1_index_genome

for SRR in ${Acc_List[@]}
	do
#		./2_sra_fastq.sh $SRR
		echo "Splitting $SRR"
		./3_split_lanes.sh $SRR
	done
  ```
