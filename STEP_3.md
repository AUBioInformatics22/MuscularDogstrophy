## Step 3: Post-alignment processing

We are currently in the process of cutting down to four sequences: two female and two male. With those four sequences they are currently running through the alignment code and then will be running them through the duplicates sequence which has been edited to prepare for the sequences once they have been run through the alignment.  

### 1. Mark the duplicates:

a. Mark duplicates using `4_duplicates.sh` script. Output is post-processed bam file with duplicates marked.  
b. Sort and index the resulting bam file for each sample.  

### 2. Quality Control:

a. Compare quality of BAM files before and after duplicates are marked.  
b. Compare the quality of different samples.  
c. Determine the quality of various aspects of BAM files. (samtools depth for coverage)  
d. Create bar chart with the coverage estimates for each sample, comparing the newly calculated coverage to raw coverage from Step 1 and mapped coverage from Step 2.  
e. Calculate percentage of duplicated reads per sample. (samtools flagstat for number of duplicate reads)  
f. Create histogram of percent duplication.  

### 3. IGV Genome Viewer:

a. Take screen shots from IGV of before and after to show quality improvement during Step 3 for each individual BAM file.
b. Select screenshots that best display the information of sample improvement, highlights, and important differences between specific sets of sequences.

### Discussion:

WHAT we did

### Figures

#### QC

Coverage Estimates:
<img src="analysis/0_figures/3_coverage.png"  alt="Coverage of Sequences">

Percent Duplication:
<img src="analysis/0_figures/percent_duplicates.png"  alt="Percent Duplication of Sequences">


#### IGV

Screenshot Highlights:
_pictures here_
