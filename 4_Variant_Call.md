## Step 4: Variant calling

### Process

The script `5_variants.sh` was used to perform variant calling, utilizing `GATK HaplotypeCaller`. The script also extracts only SNP variants using `GATK SelectVariants` and filters variants using `GATK VariantFiltration`. Coverage was determined for the resulting VCF files using `vcftools`. These coverage values were plotted alongside data from previous processing steps to facilitate comparison (Figure 1). The VCF files were compared to the BAM files using IGV.

### Discussion

#### Statistical analysis

We plan to compare the resulting statistical analysis and summary of the datasets. We will include an in depth summary of our data with comments on various quality aspects of our VCF files. We plan to use IGV to show examples of regions where a SNP is considered high quality versus low quality. Our pipeline has generated a high volume of SNP’s and we are in the process of determining those that are of high quality before isolating the area of interest for our final data analysis.

#### Coverage

The coverages from Step 1 were calculated for the raw whole genome. Compared to what was observed in the raw data coverage, the aligned chromosome coverage decreased, then further decreased after marking duplicates. However, when comparing the coverages calculated for all Chromosome X variants, the coverage increased significantly compared with the raw whole genome. For each sample, SNP coverage was higher than indels. Between samples, the first two samples (0001 and 0002) had coverages that were very similar. Sample 0005 had markedly higher coverage, and Sample 0006 had the highest coverage overall at all stages of processing (Table 1).

### Figures

<img src="analysis/0_figures/4_coverage.png"  alt="Bar Graph Comparing Coverage at Different Stages of the Pipeline">  

__Figure 1.__ A bar graph showing coverage for each sample at different stages of processing (values in Table 1).

| Sample ID | Raw Whole Genome | Aligned Chromosome X | Marked Chromosome X |Variants Chromosome X SNPs   | Variants Chromosome X Indels|
|:---------:|:----------------:|:--------------------:|:-------------------:|:---------------------------:| :--------------------------:|
|   0001    |      23.701      |       11.3141        |       10.2846       |           29.5533           |            24.6290          |
|   0002    |      23.454      |       11.2097        |       10.1347       |           29.5656           |            25.1147          |
|   0005    |      17.964      |       15.0552        |       13.7334       |           41.4852           |            35.4306          |
|   0006    |      22.089      |       18.4987        |       16.6656       |           47.5243           |            40.5858          |

__Table 1.__ Coverage values.

#### Contributions

Jacqueline Barry: graphical analysis and discussion  
Rebecca Nance: command line data assessment/IGV assessment    
Cassidy Schnieder: command line data assessment  
Kyndall Skelton: graphical analysis and discussion  
