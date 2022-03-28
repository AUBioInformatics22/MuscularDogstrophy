## Step 1: Initial Quality Assessment of Raw NGS Data

### Move forward or continue filtering?

__UPDATE THIS!!__ Based on the results, our group agrees that filtering to exclude poor quality reads would be appropriate. We agree that trimming of the last 10 bp should be performed for all samples due to poor quality as evidenced by Phred scores below 20.  Additional trimming may need to be performed on individual samples.  

__TODO!!__ Pick representative graphs from the selection. Add discussion about quality and selected graphs.  

<br>

### Graphical Analysis

#### Coverage

<img src="analysis/0_figures/1_coverage.png"  alt="Raw Data Whole Genome Coverage Bar Graph">  

__Figure 1.__ A bar graph showing coverage for each sample at different stages of processing (values in Table 1).  

| Sample ID | Raw Whole Genome |
|:---------:|:----------------:|
|   0001    |      23.701      |
|   0002    |      23.454      |
|   0005    |      17.964      |
|   0006    |      22.089      |

__Table 1.__ Calculated coverage of the raw data for the whole genome by sequence.  

<br>

#### FastQC Graphs

| __Sample 0001__ |  __Forward__ | __Reverse__ |
| :------: | :------: | :------: |
|Per base sequence quality|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0001">|
| __Sample 0002__ |  __Forward__ | __Reverse__ |
|Per base sequence quality |<img src="analysis/1_fastqc_reports/0002_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0002">|<img src="analysis/1_fastqc_reports/0002_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0002">| 
| __Sample 0005__ |  __Forward__ | __Reverse__ |
|Per base sequence quality|<img src="analysis/1_fastqc_reports/0005_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0005">|<img src="analysis/1_fastqc_reports/0005_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0005">|
| __Sample 0006__ |  __Forward__ | __Reverse__ |
|Per base sequence quality|<img src="analysis/1_fastqc_reports/0006_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0006">|<img src="analysis/1_fastqc_reports/0006_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0006">|
