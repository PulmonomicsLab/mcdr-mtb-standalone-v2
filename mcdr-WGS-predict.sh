#sudo apt-get install trim-galore bwa samtools freebayes

if [ $# -ne 2 ] ; then
	echo "Error !!! Requires 2 arguments. $# argument(s) provided. Exiting ... (ERR_CODE: 1000)"
	exit 1000
fi

declare freebayes_path
declare samtools_path
declare bwa_path
declare trim_galore_path
declare vcflib_path
declare bgzip_path
declare bcftools_path

declare trim_galore_cores
declare bwa_mem_cores
declare samtools_cores

declare script_path
declare input_folder
declare input_id
declare input_count
declare output_folder
declare output_path
declare -a files
declare -a input_ids

script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )/"
source $script_path"config.sh"

input_folder=$1
output_folder=$2

if [ ! -d "$input_folder" ] ; then
	echo "Error !!! Invalid 'input_folder' given. Exiting ... (ERR_CODE: 1001)"
	exit 1001
fi
if [ ! -d "$output_folder" ] ; then
	echo "Error !!! Invalid 'output_folder' given. Exiting ... (ERR_CODE: 1002)"
	exit 1002
fi

if [ "${input_folder[@]: -1: 1}" != "/" ] ; then
	input_folder=$input_folder"/"
fi

if [ "${output_folder[@]: -1: 1}" != "/" ] ; then
	output_folder=$output_folder"/"
fi

# files=(`ls $input_folder*.fastq.gz`)
files=(`find $input_folder -maxdepth 1 -type f -name "*.fastq.gz"`)
for ((i=0,j=0; i<${#files[@]}; ++i)) ; do
	declare temp
	temp=${files[$i]##*/}
	temp=${temp%%.*}
	temp=${temp%_*}
	if [[ ! " ${input_ids[*]} " =~ " $temp " ]]; then
		input_ids[$j]=$temp
		((j=j+1))
	fi
# 	echo $j" => "$temp" => "${input_ids[j-1]}
	unset temp
done

input_count=${#input_ids[@]}

if [[ $input_count -eq 0 ]] ; then
	echo "Error !!! No input .fastq files. Exiting ... (ERR_CODE: 1001)"
	exit 1001
fi

echo "Running with following arguments:"
echo "script_path = "$script_path
echo "input_folder = "$input_folder
echo "Sample IDs = {"${input_ids[@]}"}"
echo "output_folder = "$output_folder
echo "reference_genome = "$script_path"data/GCF_000195955.2_ASM19595v2_genomic.fna"
echo "freebayes_path = "$freebayes_path
echo "samtools_path = "$samtools_path
echo "bwa_path = "$bwa_path
echo "trim_galore_path = "$trim_galore_path
echo "vcflib_path = "$vcflib_path
echo "bgzip_path = "$bgzip_path
echo "bcftools_path = "$bcftools_path
echo "trim_galore_cores = "$trim_galore_cores
echo "bwa_mem_cores = "$bwa_mem_cores
echo "samtools_cores = "$samtools_cores

echo ""

((k=1))
for input_id in ${input_ids[@]} ; do
	output_path=$output_folder$input_id"/"
	echo "($k / $input_count) Running for $input_id (output_path = $output_path) ..."
	echo ""

	echo "Doing mkdir(s) ..."
	mkdir $output_path
	mkdir $output_path"reference"
	cp $script_path"data/GCF_000195955.2_ASM19595v2_genomic.fna" $output_path"reference"
	echo "Done mkdir(s)"
	echo ""

	#quality preprocess
	echo "Doing trim-galore ..."
	$trim_galore_path --cores "$trim_galore_cores" --paired  $input_folder$input_id\_1.fastq.gz $input_folder$input_id\_2.fastq.gz --output_dir $output_path >/dev/null 2>&1
	echo "Done trim-galore. Exit status $?"
	echo ""

	#indexing reference
	echo "Doing bwa index ..."
	$bwa_path index $output_path"reference/GCF_000195955.2_ASM19595v2_genomic.fna" >/dev/null 2>&1
	echo "Done bwa index. Exit status $?"
	echo ""

	#alignment
	echo "Doing bwa mem ..."
	$bwa_path mem -t "$bwa_mem_cores" $output_path"reference/GCF_000195955.2_ASM19595v2_genomic.fna" -R "@RG\tID:$input_id\_\tSM:$input_id\_\tPL:ILLUMINA\tLB:$input_id\_" $output_path$input_id\_1_val_1.fq.gz $output_path$input_id\_2_val_2.fq.gz | samtools sort "-@"$samtools_cores -o $output_path$input_id.bam
	echo "Done bwa mem. Exit status $?"
	echo ""

	#indexing the aligned bam
	echo "Doing samtools index ..."
	$samtools_path index $output_path$input_id.bam >/dev/null 2>&1
	echo "Done samtools index. Exit status $?"
	echo ""

	#sorting the bam on names of the reads
	echo "Doing samtools sort ..."
	$samtools_path sort -n "-@"$samtools_cores -o $output_path$input_id\_namesort.bam $output_path$input_id.bam >/dev/null 2>&1
	echo "Done samtools sort. Exit status $?"
	echo ""

	#fix the paired mates
	echo "Doing samtools fixmate ..."
	$samtools_path fixmate -m "-@"$samtools_cores $output_path$input_id\_namesort.bam $output_path$input_id\_fix.bam >/dev/null 2>&1
	echo "Done samtools fixmate. Exit status $?"
	echo ""

	#sorting the bam on positions
	echo "Doing samtools sort ..."
	$samtools_path sort "-@"$samtools_cores -o $output_path$input_id\_positionsort.bam $output_path$input_id\_fix.bam >/dev/null 2>&1
	echo "Done samtools sort. Exit status $?"
	echo ""

	#marking the duplicates
	echo "Doing samtools markdup ..."
	$samtools_path markdup "-@"$samtools_cores $output_path$input_id\_positionsort.bam $output_path$input_id\_markdup.bam >/dev/null 2>&1
	echo "Done markdup. Exit status $?"
	echo ""

	#variant calling (ploidy 1)
	echo "Doing freebayes ..."
	$freebayes_path -f $output_path"reference/GCF_000195955.2_ASM19595v2_genomic.fna" -p 1 $output_path$input_id\_markdup.bam > $output_path$input_id.vcf
	echo "Done freebayes. Exit status $?"
	echo ""

	if [[ $input_count -gt 1 ]] ; then
		echo "Doing zip and index ..."
		$bgzip_path -c $output_path$input_id.vcf > $output_path$input_id.vcf.gz && $bcftools_path index $output_path$input_id.vcf.gz
		echo "Done zip and index. Exit status $?"
		echo ""
	fi

	echo "($k / $input_count) Done for $input_id ..."
	echo ""

	((k=$k+1))

done

if [[ $input_count -gt 1 ]] ; then
	echo "Doing VCF merge ..."
	declare -a zipped_vcf_files
	for ((i=0; i<${#input_ids[@]}; ++i)) ; do
		zipped_vcf_files[$i]=$output_folder${input_ids[$i]}"/"${input_ids[$i]}".vcf.gz"
	done
	# echo ${zipped_vcf_files[@]}
	$bcftools_path merge -m all ${zipped_vcf_files[@]} > $output_folder"merged.vcf"
	echo "Done VCF merge. Exit status $?"
	echo ""

	echo "Doing vcf2tsv ..."
	$vcflib_path vcfbreakmulti $output_folder"merged.vcf" | $vcflib_path vcf2tsv -g > $output_folder"merged.tsv"
	echo "Done vcf2tsv. Exit status $?"
	echo ""

	echo "Doing predictions ..."
	Rscript $script_path"predict_batch.R" $script_path"data/" $output_folder"merged.tsv" $output_folder
	echo "Done predictions"
else
	echo "Doing vcf2tsv ..."
	$vcflib_path vcfbreakmulti $output_path${input_ids[0]}".vcf" | $vcflib_path vcf2tsv -g > $output_folder${input_ids[0]}".tsv"
	echo "Done vcf2tsv. Exit status $?"
	echo ""

	echo "Doing predictions ..."
	Rscript $script_path"predict.R" $script_path"data/" $output_folder${input_ids[0]}".tsv" $output_folder
	echo "Done predictions"
fi
