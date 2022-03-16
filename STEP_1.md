## Step 1: Initial Quality Assessment of Raw NGS Data

### Move forward or continue filtering?
Based on the results, our group agrees that filtering to exclude poor quality reads would be appropriate. Despite inconsistencies in the FASTQC quality reports between Galaxy and ASC, we agree that trimming of the last 10 bp should be performed for all samples due to poor quality as evidenced by Phred scores below 20. Additional trimming may need to be performed on individual samples.  
<br>

### Graphical Analysis
__Per Base Sequence Quality__

| Sample | Per Base Quality Forward | Per Base Quality Reverse |
| :------: | :------: | :------: |
|0001|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0001">|
|0002|<img src="analysis/1_fastqc_reports/0002_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0002">|<img src="analysis/1_fastqc_reports/0002_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0002">|
|0003|<img src="analysis/1_fastqc_reports/0003_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0003">|<img src="analysis/1_fastqc_reports/0003_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0003">|
|0004|<img src="analysis/1_fastqc_reports/0004_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0004">|<img src="analysis/1_fastqc_reports/0004_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0004">|
|0005|<img src="analysis/1_fastqc_reports/0005_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0005">|<img src="analysis/1_fastqc_reports/0005_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0005">|
|0006|<img src="analysis/1_fastqc_reports/0006_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0006">|<img src="analysis/1_fastqc_reports/0006_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0006">|
|0007|<img src="analysis/1_fastqc_reports/0007_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0007">|<img src="analysis/1_fastqc_reports/0007_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0007">|
|0008|<img src="analysis/1_fastqc_reports/0008_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0008">|<img src="analysis/1_fastqc_reports/0008_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0008">|

### FastQC Graphs for Sample 0001

| Quality Check | Forward | Reverse |
| :------: | :------: | :------: |
|[PASS] Per base sequence quality|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/per_base_quality.png"  alt="Per Base Quality of Reverse Reads from Sample 0001">|
|[WARNING] Per tile sequence quality|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/per_tile_quality.png"  alt="Per Tile Quality of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/per_tile_quality.png"  alt="Per Tile Quality of Reverse Reads from Sample 0001">|
|[PASS] Per sequence quality scores|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/per_sequence_quality.png"  alt="Per Sequence Quality of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/per_sequence_quality.png"  alt="Per Sequence Quality of Reverse Reads from Sample 0001">|
|[WARNING] Per base sequence content|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/per_base_sequence_content.png"  alt="Per Base Sequence Content of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/per_base_sequence_content.png"  alt="Per Base Sequence Content of Reverse Reads from Sample 0001">|
|[PASS] Per base N content|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/per_base_n_content.png"  alt="Per Base N Content of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/per_base_n_content.png"  alt="Per Base N Content of Reverse Reads from Sample 0001">|
|[PASS] Sequence Length Distribution|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/sequence_length_distribution.png"  alt="Sequence Length Distribution of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/sequence_length_distribution.png"  alt="Sequence Length Distribution of Reverse Reads from Sample 0001">|
|[PASS] Sequence Duplication Levels|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/duplication_levels.png"  alt="Sequence Duplication Levels of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/duplication_levels.png"  alt="Sequence Duplication Levels of Reverse Reads from Sample 0001">|
|[PASS] Adapter Content|<img src="analysis/1_fastqc_reports/0001_1_fastqc/Images/adapter_content.png"  alt="Adapter Content of Foward Reads from Sample 0001">|<img src="analysis/1_fastqc_reports/0001_2_fastqc/Images/adapter_content.png"  alt="Adapter Content of Reverse Reads from Sample 0001">|
