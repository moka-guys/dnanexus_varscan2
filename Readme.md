# DNAnexus Varscan2 v1.8.0
## What does this app do?
This app applies Varscan2 ([v2.4.3](https://dkoboldt.github.io/varscan/)), a variant caller well suited for somatic samples.

Varscan2 variant calling is performed after alignment and produces a vcf. If a BED file is also given a second, bed filtered VCF is produced.

## What are typical use cases for this app?
Varscan2 is used to detect variation (CNV and SNV) in somatic cancer samples being tested using the SWIFT amplicon panels. 

The output vcf(s) can be uploaded into Ingenuity for variant interpretation.

## What inputs are required for this app to run?
This app requires the following data:
- Compressed reference genome including `*.fa` and `*.fa.fai` (`*.tar.gz`)
- BAM file(s) (`*.bam`). If multiple BAM files are given a seperate analysis will be performed on each BAM file.
- BED file of regions of interest, for filtering output vcf (`*.bed`) (optional) - NB if provided, ensure there are no overlapping regions, as any variants found within multiple regions will be called multiple times.

This app requires the following inputs for mpileupcns [app default]:
-	min-coverage: Minimum read depth at a position to make a call [10]
-	min-reads2: Minimum supporting reads at a position to call variants [5]
-	min-var-freq: Minimum variant allele frequency threshold [0.01]
-	min-freq-for-hom: Minimum frequency to call homozygote [0.75]
-	p-value	Default: p-value threshold for calling variants [0.05]
-	strand-filter: Ignore variants with >90% support on one strand [True (1)]
-	output-vcf: Outputs in VCF format [True (1)]
-	variants: Report only variant (SNP/indel) positions [True (1)]
-	samplename: samplename to be included in the VCF headerline (string). If not provided the bamfile prefix will be used.

The following samtools mpileup can be specified. if not given the app defaults are applied (stated in square brackets):
-	min-MQ: Minimum mapping quality for an alignment to be used (-q) [20]
-   extra arguments (string)

Minimum base call quality
-	min-BQ: Minimum base quality for a base to be considered. Used by varscan (--min-avg-qual) and mpileup(-Q) [15]

## How does this app work?
- The app loops through the list of input BAM files
- It creates a flagstat file using `samtools flagstat` (samtools v1.9)
- It checks if the BAM file is empty using `samtools view -c` - if it is the BAM file is skipped 
- If there are aligned reads `samtools mpileup` creates an mpileup file.
- Variant calling is performed on the mpileup file using `varscan mpileup2cns` (Varscan2 v2.4.3)
- If a BED file is supplied `bedtools intersect` is used to filter the vcf (removing variants from off-target alignment) (bedtools v2.25.0)


## What does this app output?
This app will output:
- A vcf file for each sample detailing CNV or SNV variants called. As of V1.7.1 the reference genome will be stated in the VCF header
- If a BED file is provided a bed filtered vcf file for each sample, detailing CNV or SNV variants called.
- samtools mpileup file. This is the input to varscan and can be used to indicate the read depth available to Varscan.
- samtools flagstat file. Summarises the bitwise flags present in the BAM file. This can be used to explain an empty mpileup file which does not get processed by Varscan. This output can be displayed by MultiQC.

For detailed information about the analysis, consult the [Varscan manual](https://dkoboldt.github.io/varscan/using-varscan.html) and [Samtools flagstat](http://www.htslib.org/doc/samtools.html).

VCF files are output to `/output`, mpileup to `/coverage/mpileup` and flagstat to `/QC`.
