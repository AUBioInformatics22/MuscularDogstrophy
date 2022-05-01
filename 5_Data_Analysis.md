## Step 5: Data Analysis

### Overview
We manually inspected these variants by uploading the vcf file containing the common SNPs among all dogs into IGV (https://software.broadinstitute.org/software/igv/). The annotated reference file (canFam6.chrX.json) was loaded into IGV along with the VCF file containing the SNPs of interest (all_samples_SNPs.vcf.gz). This variant file contains SNPs on the X chromosome that are hemizygous in all affected males and heterozygous in all carrier females. We then zoomed into the DMD gene to visualize and analyze the SNPs present (Figure 

### Discussion
We first looked for SNPs that were present in exons on the DMD gene. Though we expected to find none, since prior cDNA sequencing indicated no mutations, one SNP was identified within an exon (at position 32505025). This SNP, located in exon 34 (position 32504947-32505117), was hemizygous for C in the affected males, and heterozygous for A/C in the carrier females. Therefore, the mutation causes a change from C to A at position 32505025 on the X chromosome. This is a silent mutation because it does not affect the encoded amino acid (Alanine). Exon 34 is translated from reading frame 1.

<br> 

<img src="analysis/0_figures/DMD_gene_SNPs.png"  alt="SNPs in DMD Gene">  

__Figure 1.__ A screenshot of IGV showing all SNPs of interest within the dystrophin (DMD) gene. Light blue blocks correspond to homozygous (hemizygous) SNPs within the males while dark blue corresponds to heterozygous SNPs within the females.

<img src="analysis/0_figures/DMD_exon_mx.png"  alt="Exon SNP">
<img src="analysis/0_figures/DMD_exon_mx_zoom.png" >

__Figure 2.__ Screenshot of the SNP present within a DMD exon (bottom=zoomed in). 
 
### Discrepancy in IGV Reading Frame

<img src="analysis/0_figures/IGV_vs_UCSC_.png"  alt="IGV vs. UCSC Reading Frame">  

 __Figure 3.__ IGV showed a discrepancy in how the reading frames were displayed. 
 
### Future Directions  
We plan to continue this work by adding more samples (2 additional affected males and 2 additional carrier females) and analyze for INDELs as well. We hope to implement some additional tools (such as Snpeff) to narrow our variants to identify functionally relevant SNPs/INDELs. In this way, we can look for mutations that are more likely to produce an effect on the DMD protein. To validate the final results, we will sequence the PCR product of an unaffected/non-carrier female Springer Spaniel. 
