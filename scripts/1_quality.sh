#!/bin/sh

### This script uses FASTQC in parallel to generate quality reports for
### each sequence file. The output is a fastqc.html file for each sequence
### file that can be opened in a browser to view the graphical results.

#### QUEUE PARAMETERS ####
# class
# 4 threads 
# 2gb memory  
# 12:00:00  

#### MODULES #### 
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load gnu_parallel
module load fastqc

#### DIRECTORIES #### 
DATADIR=/home/aubclsb0103/MD_data/data/raw_data2/Springer_DMD

#####################

#move into directory where the sequence files are located
cd $DATADIR

#run fastqc in parallel on each sequence file (forward and reverse)
ls *.fastq.gz | parallel -j+0 --eta 'fastqc {}'
