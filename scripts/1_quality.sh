#!/bin/sh

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load gnu_parallel
module load fastqc

DATADIR=/home/aubclsb0103/MD_data/data/raw_data2/Springer_DMD
cd $DATADIR

#fastqc -t 8 *_2_*.fastq.gz
#fastqc -t 2 H2LCVDSX3_s1_2_IDT8_UDI_024_i7-IDT8_UDI_024_i5_6797-RN-0001.fastq.gz H2LCVDSX3_s1_2_IDT8_UDI_027_i7-IDT8_UDI_027_i5_6797-RN-0004.fastq.gz

parallel -j+0 --eta 'fastqc {}' ::: H2LCVDSX3_s1_2_IDT8_UDI_024_i7-IDT8_UDI_024_i5_6797-RN-0001.fastq.gz H2LCVDSX3_s1_2_IDT8_UDI_027_i7-IDT8_UDI_027_i5_6797-RN-0004.fastq.gz
