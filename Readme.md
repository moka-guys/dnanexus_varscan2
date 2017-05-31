# DNAnexus Varscan2.4.3 v1.1

## What does this app do?

This app is a variant and somatic mutation/CNV caller. Generates a vcf.

The app generates an mpileup from a bam file then uses Varscan2 to call variants. 

This app uses Varscan2 mpileup2cns command to make consensus calls on SNP/Indel/Reference from a mpileup file. 

## What are typical use cases for this app?

This app is typically used after sequencing alignment, to generate a vcf of variants identified. 

## What inputs are required for this app to run?

This app requires the following data:

- Compressed reference genome including `*.fa` and `*.fa.fai` (`*.tar.gz`)
- BAM files (`*.bam`)
- BED file of regions of intrest, for filtering output vcf (`*.bed`) (optional)

This app requires the following inputs [default]:

-	min-coverage:	Minimum read depth at a position to make a call [10]
-	min-reads2:	Minimum supporting reads at a position to call variants [5]
-	min-avg-qual:	Minimum base quality at a position to count a read [15]
-	min-var-freq:	Minimum variant allele frequency threshold [0.01]
-	min-freq-for-hom:	Minimum frequency to call homozygote [0.75]
-	p-value	Default: p-value threshold for calling variants [0.05]
-	strand-filter:	Ignore variants with >90% support on one strand [False (0)]
-	output-vcf:	Outputs in VCF format [True (1)]
-	variants:	Report only variant (SNP/indel) positions [True (1)]


## What does this app output?

This app outputs a vcf (or.txt) per sample detaling all variants identified. For detailed information about the analysis, consult the Varscan manual at:

http://varscan.sourceforge.net/using-varscan.html


## How does this app work?

This app runs Samtools to generate an mpileup of each bam file.
Varscan2, then makes consensus calls (SNP/Indel/Reference) from a mpileup file based on user-defined parameters. By defult outputs in vcf format.


## Custom modifications
If a bed file is provided, the varscan output (`*.vcf`) are filtered to regions of interest denoted in the bed file. Bedfiltered vcfs are moved to a seperate directory# dnanexus_varscan2
