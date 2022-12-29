#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#Grab inputs
dx-download-all-inputs --except ref_genome --parallel

# make output folders
mkdir -p ~/out/varscan_vcf/output ~/out/varscan_vcf_bed/output ~/out/flagstat/QC/ ~/out/mpileup_file/coverage/mpileup/

# compile user specified options/inputs required to run Varscan. Append optional inputs, if specified.
opts=" --min-coverage $min_coverage --min-reads2 $min_reads2"

if [ "$min_BQ" != "" ]; then
  opts="$opts --min_avg_qual $min_BQ"
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

# Parse the reference genome file name for the genome build used and set $genomebuild appropriately.
# This reference build will be added to the VCF header.
genomebuild="unknown"
if [[ $ref_genome_name =~ .*37.* ]]
then
	genomebuild="grch37"
elif [[ $ref_genome_name =~ .*19.* ]]
then
	genomebuild="hg19"
elif [[ $ref_genome_name =~ .*38.* ]]
then
	genomebuild="grch38"
elif [[ $ref_genome_name =~ .*20.* ]]
then
	genomebuild="hg20"
else
	echo "$ref_genome_name does not contain a parsable reference genome name"
fi

# rename reference genoms
mv genome/*.fa  genome/$genomebuild.fa
mv genome/*.fa.fai  genome/$genomebuild.fa.fai
# capture the fasta file as a variable for mpileup
genome_file=`ls genome/*.fa`

# Show the java version the worker is using
# Note: Have only specified java8 in json 
echo $(java -version)


# Calculate 80% of memory size, for java
head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}' >.mem_in_mb.txt
java="java -Xmx$(<.mem_in_mb.txt)m"

mark-section "Run Varscan"
# loop through array of all bam files input, run varscan for each bam file. 
for (( i=0; i<${#bam_file_path[@]}; i++ )); 
# show name of current bam file
do echo ${bam_file_prefix[i]}
# if $samplename is not given set it to the bamfile prefix
if [ "$samplename" == "" ]; then
	samplename=${bam_file_prefix[i]}
fi
# generate a flagstat output 
samtools flagstat  ${bam_file_path[i]} > ~/out/flagstat/QC/${bam_file_prefix[i]}.flagstat

#check if BAM is empty
if [ $(samtools view -c ${bam_file_path[i]}) -eq 0 ]; then
	# skip and write to stdout
	echo "empty BAM. skipping...."

# if not empty perform variant calling
else
	# build the argument string, including the optional inputs if required 
	# -a outputs all abses, even is 0 coverage
	# -B disables BAQ
	# -d max number of reads to count (saves memory - set very high to override default)
	mpileup_opts="-a -B -d 500000"
	if [ "$min_MQ" != "" ]; then
	mpileup_opts="$mpileup_opts -q $min_MQ"
	fi
	if [ "$min_BQ" != "" ]; then
	mpileup_opts="$mpileup_opts -Q $min_BQ"
	fi
	if [ "$bed_file" != "" ]; then
	mpileup_opts="$mpileup_opts -l $bed_file_path"
	fi
	if [ "$mpileup_extra_opts" != "" ]; then
	mpileup_opts="$mpileup_opts $mpileup_extra_opts"
	fi
	# generate an mpileup from bam file
	samtools mpileup -f $genome_file $mpileup_opts ${bam_file_path[i]} > out/mpileup_file/coverage/mpileup/${bam_file_prefix[i]}.mpileup
	
	# test if the mpileup file is empty - if it is skip varscan variant calling
	if [ $(cat out/mpileup_file/coverage/mpileup/${bam_file_prefix[i]}.mpileup | wc -l ) -eq 0 ]; then
		# skip and write to stdout
		echo "empty mpileup file. skipping...."
	else
		#Call varscan on mpileup file using mpileupcns function. write vcf direct to output folder
		cat out/mpileup_file/coverage/mpileup/${bam_file_prefix[i]}.mpileup | $java -jar /usr/bin/VarScan.v2.4.3.jar mpileup2cns $opts > out/varscan_vcf/output/${bam_file_prefix[i]}.varscan.vcf
		# Rename sample in vcf to given samplename. Varscan defult is to name samples 'Sample1'
		sed -i 's/Sample1/'"$samplename"'/' ~/out/varscan_vcf/output/${bam_file_prefix[i]}.varscan.vcf
		# add the reference genome into the last line of the header
		sed -i "s/#CHROM/##REFERENCE=$genomebuild\n#CHROM/" ~/out/varscan_vcf/output/${bam_file_prefix[i]}.varscan.vcf

		# if bedfile provided filter vcf to contain only variants within genomic regions specified by the bed file.
		if [ "$bed_file" != "" ]; then
			# use sed to remove chr from chromosome in bedfile. write to temp vcf (not output from this app)
			sed 's/chr//' ~/out/varscan_vcf/output/${bam_file_prefix[i]}.varscan.vcf > ${bam_file_prefix[i]}.temp.vcf
			# use bedtools to instersect BED file and VCF file. Write VCF direct to the output folder
			bedtools intersect -header -a ${bam_file_prefix[i]}.temp.vcf -b ${bed_file_path} > ~/out/varscan_vcf_bed/output/${bam_file_prefix[i]}.varscan.bedfiltered.vcf
		fi
	fi
fi
done 

# upload output vcfs
mark-section "Upload output"
dx-upload-all-outputs --parallel
mark-success
