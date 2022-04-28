## Scripts

### Pipeline scripts

- **[0_index_canFam6.sh](scripts/0_index_canFam6.sh):**
  - Index the canine reference genome (canFam6) using `bwa index`.
  - Index canFam6 using `samtools faidx`.
  - Index canFam6 using `picard CreateSequenceDictionary`.
- **[1_quality.sh](scripts/1_quality.sh):**
  - Generate quality assessment reports using `FastQC`.
- **[2_align_chrX.sh](scripts/2_align_chrX.sh):**
  - Align the samples to the reference genome (canFam6).  
  - Subset the X chromosome using `samtools view`.  
  - Summarize the alignment quality of sequences using `samtools flagstat` and `samtools depth`.
- **[3_duplicates.sh](scripts/3_duplicates.sh):**
  - Sort, index, and mark duplicates for the BAM files.
- **[replace_read_group.sh](scripts/replace_read_group.sh):**
  - Replace the read groups for each sample based on unique metadata.
- **[4_variant_call.sh](scripts/4_variant_call.sh):**
  - Call variants using `GATK HaplotypeCaller`.
  - Extract SNP variants using `GATK SelectVariants`.
  - Hard filter variants using `GATK VariantFiltration`.
  - Determine coverage for the resulting VCF files using `vcftools`.
- **[5_select_variant.sh](scripts/5_select_variant.sh)**
  - Select variants that are hemizygous for an alternate allele in the males using `bcftools view` with the `-i 'GT="alt"'` option.
  - Select variants that are heterozygous for an alternate allele in the females using `bcftools view` with the `-i 'GT="het"'` option.
  - Intersect the selected variants and find positions that are in all 4 samples using `bcftools isec`.
  - Merge files containing selected variants using `bcftools merge`.

### Analysis scripts

- **[dup_metrics.sh](scripts/dup_metrics.sh):**
  - Determine coverage of the marked BAM files using `samtools depth`.
  - Get metrics using `samtools flagstat`.
- **[stat_csv.sh](scripts/stat_csv.sh):** 
  - Create csv files that can be used to make graphs in R (create_figures.R).
  - Calculate raw coverage metrics from fastqc_data.txt files output by `1_quality.sh`.
  - Get aligned chrX coverage data from .coverage_summary.txt files output by `2_align_chrX.sh`.
  - Get marked chrX coverage data from .flagstat_output.txt files output by `dup_metrics.sh`.
  - Get variants coverage data from .idepth files output by `4_variant_call.sh`.
  - Get percent mapped data from .flagstat_output.txt files output by `2_align_chrX.sh`.
  - Get percent percent duplicates data from .dup_metrics files output by `3_duplicates.sh`.
- **[create_figures.R](scripts/create_figures.R):**
  - Generate a bar plot of coverage, including raw and aligned data.  
  - Generate a bar plot of coverage, including raw, aligned, and marked data.
  - Generate a bar plot of coverage, including raw, aligned, marked, and variant data.
  - Generate a bar plot of percent mapped.  
  - Generate a bar plot of percent duplicates.
