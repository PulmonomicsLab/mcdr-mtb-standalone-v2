if [ $# -ne 2 ] ; then
	echo "Error !!! Requires 2 arguments. $# argument(s) provided. Exiting ... (ERR_CODE: 1000)"
	exit 1000
fi

declare vcflib_path

declare script_path
declare input_folder
declare input_count
declare output_folder
declare -a files
declare -a input_ids

script_path="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )/"
source $script_path"config.sh"

input_folder=$1
output_folder=$2

source $script_path"config.sh"

if [ ! -d "$input_folder" ] ; then
	echo "Error !!! Invalid 'input_folder' given. Exiting ... (ERR_CODE: 1001)"
	exit 1001
fi
if [ ! -d "$output_folder" ] ; then
	echo "Error !!! Invalid 'output_folder' given. Exiting ... (ERR_CODE: 1002)"
	exit 1002
fi

if [ "${input_folder[@]: -1: 1}" != "/" ] ; then
	input_folder="$input_folder/"
fi

if [ "${output_folder[@]: -1: 1}" != "/" ] ; then
	output_folder="$output_folder/"
fi

# files=(`ls $input_folder*.vcf`)
files=(`find $input_folder -maxdepth 1 -type f -name "*.vcf"`)
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
	echo "Error !!! No input VCF file. Exiting ... (ERR_CODE: 1001)"
	exit 1001
elif [[ $input_count -eq 1 ]] ; then
	echo "Error !!! Single VCF file found in input folder. Exiting ... (ERR_CODE: 1001)"
	echo "Use 'predict_single.sh' for prediction from single files."
	exit 1001
fi

echo "Running with following arguments:"
echo "script_path = $script_path"
echo "input_folder = $input_folder"
echo "Sample IDs = {"${input_ids[@]}"}"
echo "output_folder = $output_folder"
echo "vcflib_path = $vcflib_path"

echo ""

echo "Doing zip and index ..."
for ((i=0; i<${#input_ids[@]}; ++i)) ; do
	$bgzip_path -c "$input_folder${input_ids[$i]}.vcf" > $output_folder${input_ids[$i]}.vcf.gz && $bcftools_path index $output_folder${input_ids[$i]}.vcf.gz
done
echo "Done zip and index. Exit status $?"
echo ""

echo "Doing VCF merge ..."
declare -a zipped_vcf_files
for ((i=0; i<${#input_ids[@]}; ++i)) ; do
	zipped_vcf_files[$i]=$output_folder${input_ids[$i]}".vcf.gz"
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

