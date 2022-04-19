#!/bin/sh


source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc
module load gnu_parallel
module load trimmomatic

DATADIR=/home/aubclsb0103/MD_data/data/fastq
OUTDIR=/home/aubclsb0103/MD_data/data/fastq/quality_trimmed
mkdir -p $OUTDIR

cd $DATADIR

ls *.fastq.gz | cut -d "_" -f 1 | sort | uniq > list


############ Trimmomatic #############
############  Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming parameters
## MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
## SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across  
## requiredQuality: specifies the average quality required.
                
### while loop to process through the names in the list
while read i
do
java -jar /mnt/beegfs/home/aubmxa/.conda/envs/BioInfo_Tools/share/trimmomatic-0.39-1/trimmomatic.jar  PE -threads 8 -phred33 \
        $DATADIR/"$i"_1.fastq.gz $DATADIR/"$i"_2.fastq.gz     \
        "$i"_1_paired.fastq.gz "$i"_1_unpaired.fastq.gz       \
        "$i"_2_paired.fastq.gz "$i"_2_unpaired.fastq.gz       \
        ILLUMINACLIP:Adapters.fa:2:30:10
done<list


############### Now assess Quality with FASTQC again    ##############

ls *_1_paired.fastq.gz | parallel -j+0 --eta 'fastqc {}'
ls *_2_paired.fastq.gz | parallel -j+0 --eta 'fastqc {}'


### copy the results output to my directory for safe keeping
mv *fastqc* $OUTDIR
