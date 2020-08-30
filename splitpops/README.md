# splitpops

Split a file of Illumina paired-end reads for an Hepatitis C Virus (HCV) sample into one file per sub-genotype population. Only sub-genotypes that attract at least M% of the read pairs (default 1%) are reported as a population.

## Input
* Sample paired-end reads [BAM format]
* HCV target references [FASTA format, where file name and reference ID of format SUBGENOTYPE_ID, e.g., 1a_ABCDEFG] 
* Minimum percentage of sample reads to a target reference subgenotype matches for the  Target References [range 0 to 100, default 1.0]

## Method
* Map each of the sample reads to the target references of known subgenotype, discarding any non-primary or supplementary mappings.
* If  at least one read of the pair were able to be mapped to a target reference genome, the pair was allocated the sub-genotype of the highest-scoring match. If both mapped to different sub-genotypes, the allocated sub-genotype read pair was marked as "ambigous". If neither read mapped, the sub-genotype was designated "unclassified".

## Output
* A summary of the populations found as a tab-separated file containing "dataid, popN, readclass, readcnt, readpct, bampath", where popN is "1" for the populatio with the highest percentage, readclass is a sub-genotype (e.g., 1a, 2b), and bampath is the absolute path name to a BAM file (e.g., DATAID.splitpops.pop1.bam) [DATAID.splitpops.stats.txt]
* A set of BAM files containing read pairs from each population [DATAID.splitpops.pop1.bam, DATAID.splitpops.pop2.bam, ...]

Camilla Ip, 2017  
The Wellcome Centre for Human Genetics, University of Oxford, UK
