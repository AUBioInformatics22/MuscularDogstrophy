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

### Side Note

We haven't been able to get the alignment code to work because the code would run for upwards of three hours and then get an out of memory error, so there is a long delay with troubleshooting. We have tried a variety of combinations of cores and memory, but have not yet been successful. There were no errors within the alignment code other than the out of memory errors. The code has been broken up by taking out the pipes and putting them on seperate lines to avoid the memory error but continued to fail. After making this change, we have run the code with 8 cores with 2gb each as well as 10 cores with 4 gb each, with neither of those being successful. We are still troubleshooting the scripts so this readme will be updated with the data and histograms when the real data has been used and hopefully the out of memory error avoided.
