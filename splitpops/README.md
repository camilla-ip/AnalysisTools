# splitpops

Split a file of Illumina paired-end reads for a Hepatitis C Virus (HCV) sample into one BAM file per sub-genotype population. Only sub-genotypes that attract at least M% of the read pairs (default 1%) are reported as a population.

__Input__  
- Sample paired-end reads (BAM format)
- A set of HCV target references (FASTA format, containing one FASTA entry for each HCV reference genome of known subgenotype, where file name and reference ID of format SUBGENOTYPE_ID, e.g., 1a_ABCDEFG)
- Minimum percentage of sample reads to a target reference sub-genotype matches for the target references (range 0 to 100, default 1.0)

__Method__  

- Map each of the sample reads to the target references of known sub-genotype, discarding any non-primary or supplementary mappings.
- If  at least one read of the pair were able to be mapped to a target reference genome, the pair was allocated the sub-genotype of the highest-scoring match. If both mapped to different sub-genotypes, the allocated sub-genotype read pair was marked as "ambiguous". If neither read mapped, the sub-genotype was designated "unclassified".

__Output__  
- A summary of the populations found as a tab-separated file containing "dataid, popN, readclass, readcnt, readpct, bampath", where popN is "1" for the populatio with the highest percentage, readclass is a sub-genotype (e.g., 1a, 2b), and bampath is the absolute path name to a BAM file (e.g., DATAID.splitpops.pop1.bam) (DATAID.splitpops.stats.txt)
- A set of BAM files containing read pairs from each population (DATAID.splitpops.pop1.bam, DATAID.splitpops.pop2.bam, ...)

This programme was originally written by me for the Snork HCV analysis pipeline (unpublished) at the [STOP-HCV Consortium](https://www.stop-hcv.ox.ac.uk/home), parts of which were incorporated in to the Public Health England HCV analysis pipeline published in "Technical Validation of a Hepatitis C Virus Whole Genome Sequencing Assay for Detection of Genotype and Antiviral Resistance in the Clinical Pathway", Frontiers in Microbiology, 2020.

Camilla Ip  
2017, The Wellcome Centre for Human Genetics, University of Oxford, UK
