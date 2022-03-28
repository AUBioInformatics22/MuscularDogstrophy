## Step 3: Post-alignment processing

#### Discussion:
This step of the analysis picks up at the end of step 2. The output from the alignment script from Step 2 was run through the duplicate script: `4_duplicates.sh`, which sorted, indexed, and marked duplicates for the BAM files from each of our 4 samples. `samtools depth` was used to determine the overall coverage of the BAM files (Figure 1). The quality of BAM files before and after marking duplicates was determined using `samtools flagstat` (Figure 2). The plots were generated in `R studio` using the script `xxxx.R`. The `dup_metrics.sh` script was used to determine the preceding metrics. 

3. IGV Genome Viewer:
a. Take screen shots from IGV of before and after to show quality improvement during Step 3 for each individual BAM file. b. Select screenshots that best display the information of sample improvement, highlights, and important differences between specific sets of sequences.


### Figures

<img src="analysis/0_figures/3_coverage.png"  alt="Coverage of Sequences">

__Figure 1.__ Bar plot of coverage across several stages of processing.

<img src="analysis/0_figures/percent_duplicates.png"  alt="Percent Duplication of Sequences">

__Figure 2.__ Bar plot of the found percentage of duplications.

#### IGV

Screenshot Highlights:
_pictures here_

#### Contributions
Jacqueline Barry: graphical analysis and discussion  
Rebecca Nance: IGV analysis  
Cassidy Schnieder: gathered quality metrics  
Kyndall Skelton: sorted, indexed, and marked duplicates  
