## Step 4: Variant calling

### Process

1. Variant discovery using `GATK HaplotypeCaller`.
2. Variant filtering using `GATK SelectVariants` and `GATK VariantFiltration`.
3. Quality assessment of VCF.
4. Compare quality using `VCFtools`.
5. Compare VCF to BAM using IGV.

### Discussion

Step 4 has been completed we are currently figuring out the best way to present the data and if the filtering was done correctly.

#### Comparing to FastQC of raw data

How did the quality of your file change from Step 1 when you did the FastQC?

#### Comparing to aligned BAM

How did the coverage of your BAM file change from Step 2 when you first generated an alignment to VCFtools depth calculation?

#### Comparing samples

How do the coverage of the samples in my project compare?

#### Statistical analysis

Compare the resulting statistical analysis and summary of the datasets within your group. As before, you should include in your report an in depth summary of your data and comment on various quality aspects of your VCF files. You may use screen shots of IGV to show examples of regions where a SNP is considered high quality versus low quality.

While you should likely use all of these tools for each of your individual VCF files, I do not want a long report with each of these summaries. Instead, you should choose the approach that best displays the information you want to convey. Please only include the highlights and important differences between specific sets of sequences in your group report. 

#### Figures

<img src="analysis/0_figures/4_coverage.png"  alt="Bar Graph Comparing Coverage at Different Stages of the Pipeline">  

__Figure 1.__ A bar graph showing coverage for each sample at different stages of processing (values in Table 1?).
