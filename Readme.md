# DNAnexus Varscan2.4.3 v1.2

## What does this app do?

This app applies Varscan2, a variant caller well suited for somatic samples.

Varscan2 variant calling is performed after alignment and the output vcf will be uploaded into Ingenuity for annotation and filtering.

This app will perform variant calling on multiple samples. 


## What are typical use cases for this app?

Varscan2 is used to detect variation (CNV and SNV) in somatic cancer samples being tested using the SWIFT amplicon panels. 


## What inputs are required for this app to run?

This app requires the following data:

- Compressed reference genome including `*.fa` and `*.fa.fai` (`*.tar.gz`)
- BAM file(s) (`*.bam`)
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
This app will output a vcf file for each sample detailing CNV or SNV variants called. Output vcfs will be passed to Ingenuity for annotation and filtering.
For detailed information about the analysis, consult the Varscan manual at:
http://varscan.sourceforge.net/using-varscan.html


## How does this app work?

- The app loops through the list of input BAM files
- Samtools converts each bamfile to a mpileup file 
- Variant calling is performed using Varscan2 mpileup2cns function
- If a BED file is supplied bedtools is used to filter the vcf (removing variants from off-target alignment)


## Custom modifications
If a bed file is provided, the varscan output (`*.vcf`) are filtered to regions of interest denoted in the bed file. Bedfiltered vcfs are moved to a seperate directory# dnanexus_varscan2
