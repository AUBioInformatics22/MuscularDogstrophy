# Muscular Dogstrophy
### Jacqueline Barry, Rebecca Nance, Cassidy Schneider, Kyndall Skelton


## Overview

Muscular dystrophy is a term that includes many inherited diseases that affects muscle function. Duchenne Muscular Dystrophy (DMD), the most common subtype, occurs in both humans and dogs. This disease is characterized by a defect in the gene that encodes the protein dystrophin, which plays a critical role in muscle structure and function. Consequently, movement and mobility are impaired and affected individuals progressively become weaker. Eventually, the muscles in the heart and/or diaphragm are compromised, leading to premature death or euthanasia. The dystrophin gene is located on the X chromosome and therefore, obeys traditional X-chromosome inheritance pattern. Our lab maintains one of the few canine colonies for studying DMD. Breed-specific mutations have been identified in a few breeds, including Golden Retriever, Labrador Retriever, and Corgi. However, the mutation responsible for DMD in the Springer Spaniel has yet to be elucidated. Previous attempts to identify the mutation using sequencing of cDNA from affected dogs yielded no conclusive mutation. To address this question, we used genomic sequencing data from 2 affected males and 2 carrier females from our Springer Spaniel colony. We aligned our sequences to a canine reference genome `canFam6` and specifically isolated the X chromosome from our sequences. We then marked any duplicates, calculated coverages, and searched for variants at the DMD gene on the X chromosome in order to attempt to identify the mutation associated with muscular dystrophy in canines.


## Samples
Genomic DNA was isolated from 4 affected males and 4 carrier females from the Springer Spaniel Duchenne-Like Muscular Dystrophy dog colony in Dr. Bruce F. Smith's lab at the Auburn University College of Veterinary Medicine. Samples were sequenced as 150bp paired-end reads using the NovaSeq 6000 platform with ~30x coverage at HudsonAlpha Discovery in Huntsville, AL.

| Dog | Description | Library ID | File ID |
| --- | --- | --- | --- |
| Buddy | Affected Male | SL522590 | 0001
| Dandelion | Affected Male | SL522591 | 0002
| Camelia | Carrier Female | SL522594 | 0005
| Dottie | Carrier Female | SL522595 | 0006


<br>

## Reports

#### [Step 1 Report: Initial quality assessment of raw NGS data](1_Quality_Assessment.md)
#### [Step 2 Report: Alignment of sequencing reads to reference genome](2_Alignment.md)
#### [Step 3 Report: Post-alignment processing](3_Duplicates.md)
#### [Step 4 Report: Variant calling](4_Variant_Call.md)
#### [Data Analysis](Data_Analysis.md)
<br>

## Practice Data Reports

Before the sequencing data was received, example data was used to develop the pipeline.

#### [Step 1 Practice Report: Initial quality assessment of raw NGS data](practice_reports/STEP_1.md)
#### [Step 2 Practice Report: Alignment of sequencing reads to reference genome](practice_reports/STEP_2.md)
