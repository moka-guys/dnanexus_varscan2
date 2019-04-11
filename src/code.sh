#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#Grab inputs
dx-download-all-inputs --except ref_genome --parallel

# make output folders
mkdir -p ~/out/varscan_vcf/output ~/out/varscan_vcf_bed/output ~/out/flagstat/QC/
echo "echo arguments:"
echo "min_coverage" $min_coverage
echo "min_reads2" $min_reads2
echo "min_avg_qual" $min_avg_qual
echo "min_var_freq" $min_var_freq
echo "min_freq_for_hom" $min_freq_for_hom
echo "p_value" $p_value
echo "strand_filter" $strand_filter
echo "output_vcf" $output_vcf
echo "variants" $variants

# compile user specified options/inputs required to run Varscan. Append optional inputs, if specified.
opts=" --min-coverage $min_coverage --min-reads2 $min_reads2"

if [ "$min_avg_qual" != "" ]; then
  opts="$opts --min_avg_qual $min_avg_qual"
fi

if [ "$min_var_freq" != "" ]; then
	opts="$opts --min-var-freq $min_var_freq" 
fi

if [ "$min_freq_for_hom" != "" ]; then
	opts="$opts --min-freq-for-hom $min_freq_for_hom" 
fi

if [ "$p_value" != "" ]; then
	opts="$opts --p-value $p_value" 
fi 

if [ "$strand_filter" == "true" ]; then
	opts="$opts --strand-filter 1"
else
	opts="$opts --strand-filter 0"
fi

if [ "$output_vcf" == "true" ]; then
	opts="$opts --output-vcf 1"
else
	opts="$opts --output-vcf 0"
fi

if [ "$variants" == "true" ]; then
	opts="$opts --variants"
fi

# make directory for reference genome and unpackage the reference genome
mkdir genome
dx cat "$ref_genome" | tar zxvf - -C genome  
# => genome/<ref>, genome/<ref>.ann, genome/<ref>.bwt, etc.

# rename genome files to grch37 so that the VCF header states the reference to be grch37.fa, which then allows Ingenuity to accept the VCFs (otherwise VCF header would have reference as genome.fa which Ingenuity won't accept)
mv  genome/*.fa  genome/grch37.fa
mv  genome/*.fa.fai  genome/grch37.fa.fai
# mv genome.dict grch37.dict
genome_file=`ls genome/*.fa`

# Show the java version the worker is using
# Note: Have only specified java8 in json 
echo $(java -version)


# Calculate 80% of memory size, for java
head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}' >.mem_in_mb.txt
java="java -Xmx$(<.mem_in_mb.txt)m"

# Run variant annotator
mark-section "Run Varscan VariantAnnotator"
# loop through array of all bam files input, run varscan for each bam file. 
for (( i=0; i<${#bam_file_path[@]}; i++ )); 
# show name of current bam file be run
do echo ${bam_file_prefix[i]}
#if BAM is empty
if [ $(samtools view -c ${bam_file_path[i]}) -eq 0 ]; then
	# skip and write to stdout
	echo "empty BAM. skipping...."
# if not empty perform variant calling
else
	# generate an mpileup from bam file
	samtools mpileup -f $genome_file -B -d 500000 -q 1 ${bam_file_path[i]} > ${bam_file_prefix[i]}.mpileup
	# generate a flagstat output 
	samtools flagstat  ${bam_file_path[i]} > ~/out/flagstat/QC/${bam_file_prefix[i]}.flagstat
	# test if the mpileup file is empty - if it is skip varscan variant calling
	if [ $(cat ${bam_file_prefix[i]}.mpileup | wc -l ) -eq 0 ]; then
		# skip and write to stdout
		echo "empty mpileup file. skipping...."
	else
		#Call varscan on mpileup file using mpileupcns function. write vcf direct to output folder
		cat ${bam_file_prefix[i]}.mpileup | $java -jar /usr/bin/VarScan.v2.4.3.jar mpileup2cns $opts > ~/out/varscan_vcf/output/${bam_file_prefix[i]}.varscan.vcf
		# Rename sample in vcf to corrospond to bam file name (aka sample name). Varscan defult is to name samples 'Sample1'
		sed -i 's/Sample1/'"${bam_file_prefix[i]}"'/' ~/out/varscan_vcf/output/${bam_file_prefix[i]}.varscan.vcf
		# if bedfile provided filter vcf to contain only variants within genomic regions specified by the bed file.
		if [ "$bed_file" != "" ]; then
			# use sed to remove chr from chromosome in bedfile. write to temp vcf (not output from this app)
			sed 's/chr//' ~/out/varscan_vcf/output/${bam_file_prefix[i]}.varscan.vcf > ${bam_file_prefix[i]}.temp.vcf
			# use bedtools to instersect BED file and VCF file. Write VCF direct to the output folder
			/usr/bin/bedtools2/bin/bedtools intersect -header -a ${bam_file_prefix[i]}.temp.vcf -b ${bed_file_path} > ~/out/varscan_vcf_bed/output/${bam_file_prefix[i]}.varscan.bedfiltered.vcf
		fi
	fi
fi
done 

# upload output vcfs
mark-section "Upload output"
dx-upload-all-outputs --parallel
mark-success
