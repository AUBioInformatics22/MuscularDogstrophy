## Step 2: Alignment of sequencing reads to reference genome

### Summary 

For Step 2 of our project, we indexed the reference genome using BWA, split our fastq files into read groups, and aligned our samples to the reference genome. To index the reference genome, Cassidy and Kyndall first upped the gigabytes in the 1_index_canFam6.sh script to 8, then ran the script with 1 core, 4.5 gigs, and a 2 hour time limit. To split the fastq into read groups, Allison, Cassidy, and Kyndall tried to use the 3_sample_lane_split.sh but were met with a number of errors. The script worked perfectly after Cassidy redownloaded the fastq files from SRAtoolkit, indicating that the error was due to the low quality of the original files. The last step of this project was to align the read groups to the reference genome using the 4_BWA_example.sh script, and this step was by far the most difficult for our group. We encountered a number of hiccups while running this script, including kill events, out of memory errors, and having trouble deciding what cores to use. After days of troubleshooting, we finally got the script to work by __PUT WHAT MADE IT WORK HERE__.


### Analysis
| SRA ID | Percent Mapped | Coverage |
| :-----: | :-----: | :-----: |
| SRR8541909 |  |  |
| SRR8541910 |  |  |
| SRR8541914 |  |  |

Genome Size: 2,312,802,198  
<br>

Mapped coverage  
Percent reads  
Bam file  


### Contributions 

Kyndall and Cassidy did the majority of the work on the command line with Becca and Allison helping to troubleshoot any errors, Allison wrote the report, Jackie made the histograms, and Cassidy formatted the report into markdown and put it into GitHub. 
