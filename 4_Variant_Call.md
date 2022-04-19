## Step 4: Variant calling

### Process

The script `5_variants.sh` was used to perform variant calling, utilizing `GATK HaplotypeCaller`. The script also extracts only SNP variants using `GATK SelectVariants` and filters variants using `GATK VariantFiltration`. Coverage was determined for the resulting VCF files using `vcftools`. These coverage values were plotted alongside data from previous processing steps to facilitate comparison (Figure 1). The VCF files were compared to the BAM files using IGV.

### Discussion

The coverages from Step 1 were calculated for the Raw Whole Genome. When comparing the coverages calculated for the Chromosome X Variants, the coverage is significantly greater. Compared to what was observed in the raw data coverage, the aligned chromosome coverage decreased, then was further decreased after marking duplicates. However, variant coverage increased significantly. Between samples, the first 2 samples (0001 and 0002) had coverages that were fairly equivalent. Sample 0005 had higher coverage, and Sample 0006 had the highest coverage (Table 1).

#### Statistical analysis

We plan to compare the resulting statistical analysis and summary of the datasets. We will include an in depth summary of our data with comments on various quality aspects of our VCF files. We plan to use IGV to show examples of regions where a SNP is considered high quality versus low quality. Our pipeline has generated a high volume of SNP’s and we are in the process of determining those that are of high quality before isolating the area of interest for our final data analysis.

### Figures

<img src="analysis/0_figures/4_coverage.png"  alt="Bar Graph Comparing Coverage at Different Stages of the Pipeline">  

__Figure 1.__ A bar graph showing coverage for each sample at different stages of processing (values in Table 1).

| Sample ID | Raw Whole Genome | Aligned Chromosome X | Marked Chromosome X | Variants Chromosome X |
|:---------:|:----------------:|:--------------------:|:-------------------:|:---------------------:|
|   0001    |      23.701      |       11.3141        |       10.2846       |        29.5533        |
|   0002    |      23.454      |       11.2097        |       10.1347       |        29.5656        |
|   0005    |      17.964      |       15.0552        |       13.7334       |        33.4082        |
|   0006    |      22.089      |       18.4987        |       16.6656       |        40.0527        |

__Table 1.__ Coverage values.

#### Contributions

Jacqueline Barry: graphical analysis and discussion  
Rebecca Nance: command line data assessment/IGV assessment    
Cassidy Schnieder: command line data assessment  
Kyndall Skelton: graphical analysis and discussion  
