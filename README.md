# mcdr-mtb-standalone-v2
## Introduction
**mcdr-mtb-standalone-v2** is a pipeline for prediction of drug resistance class **(MDR, Pre-XDR, XDR or Susceptible)** of *Mycobacterium tuberculosis* isolates from whole genome sequencing **(WGS)** data. The users can start with paired-end FASTQ files or variant calling output VCF files.

## Workflows

There are 4 different scripts for prediction:
* ***mcdr-WGS-predict.sh:*** Prediction of drug resistance class of single or multiple MTB WGS data (paired-end *fastq* files).
* ***mcdr-VCF-predict.sh:*** Prediction of drug resistance class of single MTB VCF file.
* ***mcdr-multi-VCF-predict.sh:*** Prediction of drug resistance class of multiple MTB VCF files.
* ***mcdr-merge-VCF-predict.sh:*** Prediction of drug resistance class of MTB merged VCF file.

## Pre-requisites
* trim-galore (version `0.6.7`) - quality check and trimming of read sequences
* bwa (`0.7.17-r1188`) - reference based alignment
* samtools (`1.13`)- processing the BAM files
* freebayes (`v1.3.6`) - variant calling
* libvcflib-tools (`1.0.7`) - processing the VCF files
* libvcflib-dev (`1.0.7`)	- processing the VCF files
* bgzip (`1.13+ds`) - zipping files
* R (`4.2.2 Patched (2022-11-10 r83330)`) - data operation and machine learning model usage


## Installation
### Step 1: Install dependent packages/tools
*(For ubuntu)*

    sudo apt-get install trim-galore bwa samtools freebayes libvcflib-tools libvcflib-dev bgzip

*(For other distributions)*

The installation steps for the different packages/tools are given in the following links:

* trim-galore - https://github.com/FelixKrueger/TrimGalore
* bwa - https://github.com/lh3/bwa
* samtools, bcftools, bgzip(htstools) - http://www.htslib.org/download/
* freebayes - https://github.com/freebayes/freebayes
* vcflib - https://github.com/vcflib/vcflib

R should be installed in the user system/PC. R installation steps are given in https://cran.r-project.org/.

### Step 2: Install mcdr-mtb-standalone-v2
  #### I. Download the software from GitHub repository
   Create a clone of the repository:

      git clone https://github.com/AbhirupaGhosh/mcdr-mtb-standalone-v2

   **Note:** Creating a clone of the repository requires `git` to be installed.

   The `git` can be installed using

      sudo apt-get install git
  **OR**

   Download using `wget`:

      wget https://github.com/AbhirupaGhosh/mcdr-mtb-standalone-v2/archive/refs/heads/main.zip
      unzip main.zip
  **Note:** wget can be installed using

      sudo apt-get install wget
  #### II. Make the shell scripts executable

    chmod +x INSTALLATION_DIR/mcdr-mtb-standalone-v2 config.sh mcdr-WGS-predict.sh mcdr-VCF-predict.sh mcdr-multi-VCF-predict.sh mcdr-merge-VCF-predict.sh

  `INSTALLATION_DIR` = Directory where mcdr-mtb-standalone-v2 is installed
  #### III. update the paths in config.sh (optional)

  The `config.sh` looks like

	freebayes_path=/usr/bin/freebayes
	samtools_path=/usr/bin/samtools
	bwa_path=/usr/bin/bwa
	trim_galore_path=/usr/bin/trim_galore
	vcflib_path=/usr/bin/vcflib
	bgzip_path=/usr/bin/bgzip
	bcftools_path=/usr/bin/bcftools
	trim_galore_cores=4
	bwa_mem_cores=4
	samtools_cores=4

Note: It shows the default paths of the executables files for `freebayes`, `samtools`, `bwa`, `trim galore!`, `vcflib`, `bgzip` and `bcftools`. The users need to update the paths of the executables, in case these tools were installed in ways other than the `apt-get install` command. 

## Usage

Initially change the directory to the directory where mcdr-mtb-standalone-v2 is installed.

    cd INSTALLATION_DIR/mcdr-mtb-standalone-v2

Different operations can be performed by calling the appropriate scripts with two command-line arguments: `INPUT_DIR` and `OUTPUT_DIR`.

`INPUT_DIR` = the path (absolute or relative) of the folder containing the input files.

`OUTPUT_DIR` = the path (absolute or relative) of the folder in which mcdr-mtb-standalone-v2 will store the outputs.

The executable script, and contents of `INPUT_DIR` and `OUTPUT_DIR` depends on the choice of operations. The different operations are explained below.

### 1) Prediction of drug resistance class of one or many MTB paired-end FASTQ files.
**(*.fastq* to drug resistance class (S, M, P, X))**

    ./mcdr-WGS-predict.sh INPUT_DIR OUTPUT_DIR

`INPUT_DIR` must contain paired end FASTQ (*<sample_id>*_1.fastq.gz & *<sample_id>*_2.fastq.gz) of 1 or more isolates.

`OUTPUT_DIR` will contain
* one folder for each ISOLATE ID. Each folder will contain
  * the VCF file (*<sample_id>*.vcf)
  * the intermediate BAM files (*<sample_id>*.bam, *<sample_id>*_sorted.bam)
* the MERGED.vcf file (only in case of multiple isolates)
* the intermediate TSV file (*<sample_id>*.tsv/MERGED.tsv)
* the prediction result performed with the full model (prediction.tsv)
* the SHAP result for each isolate with the 37-feature model (shap_result_37_features_*<sample_id>*.tsv)
* the SHAP result for each isolate with the 100-feature model (shap_result_100_features_*<sample_id>*.tsv)
* the SHAP result plot for each isolate with the 37-feature model (shap_plot_37_features_*<sample_id>*.tsv)
* the SHAP result plot for each isolate with the 100-feature model (shap_plot_100_features_*<sample_id>*.tsv)

**Note:** The merged.vcf file will not be present if there was only one isolate in the INPUT_DIR. <br/>
The prediction and the SHAP results output will also be displayed on the terminal.

----

### 2) Prediction of drug resistance class of a single MTB VCF file.
**(.vcf to drug resistance class (S, M, P, X))**

    ./mcdr-VCF-predict.sh INPUT_DIR OUTPUT_DIR

`INPUT_DIR` must contain VCF file (*<sample_id>*.vcf)

`OUTPUT_DIR` will contain
* intermediate .tsv file (*<sample_id>*.tsv)
* the prediction result performed with the full model (prediction.tsv)
* the SHAP result for the isolate with the 37-feature model (shap_result_37_features_*<sample_id>*.tsv)
* the SHAP result for the isolate with the 100-feature model (shap_result_100_features_*<sample_id>*.tsv)
* the SHAP result plot for the isolate with the 37-feature model (shap_plot_37_features_*<sample_id>*.tsv)
* the SHAP result plot for the isolate with the 100-feature model (shap_plot_100_features_*<sample_id>*.tsv)

The prediction and the SHAP results output will also be displayed on the terminal.

----

### 3) Prediction of drug resistance class of multiple MTB VCF files.
**(multiple .vcf to drug resistance class (S, M, P, X))**

    ./mcdr-multi-VCF-predict.sh INPUT_DIR OUTPUT_DIR

`INPUT_DIR` must contain more than 1 VCF file (*<sample_id_1>*.vcf, *<sample_id_2>*.vcf, … , *<sample_id_n>*.vcf)

`OUTPUT_DIR` will contain
* the MERGED.vcf file
* the intermediate .tsv file (MERGED.tsv)
* the prediction result performed with the full model (prediction.tsv)
* the SHAP result for each isolate with the 37-feature model (shap_result_37_features_*<sample_id>*.tsv)
* the SHAP result for each isolate with the 100-feature model (shap_result_100_features_*<sample_id>*.tsv)
* the SHAP result plot for each isolate with the 37-feature model (shap_plot_37_features_*<sample_id>*.tsv)
* the SHAP result plot for each isolate with the 100-feature model (shap_plot_100_features_*<sample_id>*.tsv)

The prediction and the SHAP results output will also be displayed on the terminal.

----

### 4) Prediction of drug resistance class of merged MTB VCF file.
**(merged.vcf to drug resistance class (S, M, P, X))**

    ./mcdr-merge-VCF-predict.sh INPUT_DIR OUTPUT_DIR

`INPUT_DIR` must contain a merged VCF file (MERGED.vcf) with information of one or more isolates.

`OUTPUT_DIR` will contain
* the intermediate .tsv file (MERGED.tsv)
* the prediction result performed with the full model (prediction.tsv)
* the SHAP result for each isolate with the 37-feature model (shap_result_37_features_*<sample_id>*.tsv)
* the SHAP result for each isolate with the 100-feature model (shap_result_100_features_*<sample_id>*.tsv)
* the SHAP result plot for each isolate with the 37-feature model (shap_plot_37_features_*<sample_id>*.tsv)
* the SHAP result plot for each isolate with the 100-feature model (shap_plot_100_features_*<sample_id>*.tsv)

The prediction and the SHAP results output will also be displayed on the terminal.
 
----


## Demo runs

### 1. Prediction from FASTQ
#### a. Single isolate

  Step 1. Create an Input directory

   	   mkdir /home/username/Input_Dir1

  Step 2. Get Data

  Download the whole genome sequencing FASTQ files of a MTB isolate run, ERR137249 (ERR137249_1.fastq & ERR137249_2.fastq) from https://www.ebi.ac.uk/ena/browser/view/ERR137249

  Step 3. Store these files in `Input_Dir1`
 
  Step 4. Create an Output directory

   	   mkdir /home/username/Output_Dir1

  Step 5. Go to the mcdr-mtb-standalone-v2 installation directory

   	   cd /home/username/Documents/mcdr-mtb-standalone-v2/

  Step 6. Run mcdr-WGS-predict.sh

   	   ./mcdr-WGS-predict.sh /home/username/Input_Dir1/ /home/username/Output_Dir1/

`Input_Dir1` contains ERR137249_1.fastq, ERR137249_2.fastq

`Output_Dir1` contains -
* Folder - ERR137249
* ERR137249.tsv

The ERR137249 folder contains -

* reference folder - reference genome and index files
* Trim galore outputs - ERR137249_1_val_1.fq.gz, ERR137249_2_val_2.fq.gz, ERR137249_1_trimming_report.txt, ERR137249_2_trimming_report.txt
* Bwa-mem output - ERR137249.bam
* Intermediate BAM files - ERR137249_fix.bam, ERR137249_namesort.bam, ERR137249_positionsort.bam, ERR137249_markdup.bam
* BAM index - ERR137249.bam.bai
* Freebayes output - ERR137249.vcf
* VCF compressed - ERR137249.vcf.gz
* VCF index - ERR137249.vcf.gz.csi

Printed Output -

```
ISOLATE PREDICTION   	DIFF RI
1  ERR137249  		S 0.04866350  0
```
#### b. Multiple isolates

  Step 1. Create an Input directory

    mkdir /home/username/Input_Dir2

  Step 2. Get Data

Download the whole genome sequencing FASTQ files of MTB ISOLATE runs, ERR137249 (ERR137249_1.fastq & ERR137249_2.fastq) and SRR1103491 (SRR1103491_1.fastq & SRR1103491_2.fastq) from https://www.ebi.ac.uk/ena/browser/view/ERR137249 and https://www.ebi.ac.uk/ena/browser/view/SRR1103491

  Step 3. Store these files in `Input_Dir2`
 
  Step 4. Create an Output directory

    mkdir /home/username/Output_Dir2

  Step 5. Go to the mcdr-mtb-standalone-v2 installation directory

    cd /home/username/Documents/mcdr-mtb-standalone-v2/

  Step 6. Run mcdr-WGS-predict.sh

    ./mcdr-WGS-predict.sh /home/username/Input_Dir1/ /home/username/Output_Dir2/

`Input_Dir2` contains ERR137249_1.fastq, ERR137249_2.fastq, SRR1103491_1.fastq, SRR1103491_2.fastq

`Output_Dir2` contains -
* Two folders - ERR137249, SRR1103491
* merged.vcf
* merged.tsv

Each of the ERR137249 and SRR1103491 named folder contains -

* reference folder - reference genome and index files
* Trim galore outputs - ISOLATENAME_1_val_1.fq.gz, ISOLATENAME_2_val_2.fq.gz, ISOLATENAME_1_trimming_report.txt, ISOLATENAME_2_trimming_report.txt
* Bwa-mem output - ISOLATENAME.bam
* Intermediate BAM files - ISOLATENAME_fix.bam, ISOLATENAME_namesort.bam, ISOLATENAME_positionsort.bam, ISOLATENAME_markdup.bam
* BAM index - ISOLATENAME.bam.bai
* Freebayes output - ISOLATENAME.vcf
* VCF compressed - ISOLATENAME.vcf.gz
* VCF index - ISOLATENAME.vcf.gz.csi

Printed Output -
```
ISOLATE PREDICTION   	DIFF RI
1 ERR137249  		S 0.04866350  0
2 SRR1103491  		M 0.02846038  0
```
### 2. Prediction from single VCF

  Step 1. Create an Input directory

    mkdir /home/username/Input_Dir3

  Step 2. Store a VCF file generated from variant calling of a MTB isolate based on MTB H37Rv reference genome  (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz)

  Step 3. Create an Output directory

    mkdir /home/username/Output_Dir3

  Step 4. Go to the mcdr-mtb-standalone-v2 installation directory

    cd /home/username/Documents/mcdr-mtb-standalone-v2/

  Step 5. Run mcdr-VCF-predict.sh

    ./mcdr-VCF-predict.sh /home/username/Input_Dir3/ /home/username/Output_Dir3/

`Input_Dir3` contains ERR137249.vcf

`Output_Dir3` contains - ERR137249.tsv

 Printed Output -
```
ISOLATE PREDICTION   	DIFF RI
1  ERR137249  		S 0.04866350  0
```
### 3. Prediction from multiple VCFs

  Step 1. Create an Input directory

    mkdir /home/username/Input_Dir4

  Step 2. Store multiple VCF files generated from variant calling of  MTB isolates based on MTB H37Rv reference genome  (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz)

  Step 3. Create an Output directory

    mkdir /home/username/Output_Dir4

  Step 4. Go to the mcdr-mtb-standalone-v2 installation directory

    cd /home/username/Documents/mcdr-mtb-standalone-v2/

  Step 5. Run mcdr-multi-VCF-predict.sh

    ./mcdr-multi-VCF-predict.sh /home/username/Input_Dir4/ /home/username/Output_Dir4/

`Input_Dir4` contains ERR137249.vcf, SRR1103491.vcf

`Output_Dir4` contains -
  * Compressed VCFs - ERR137249.vcf.gz and SRR1103491.vcf.gz
  * VCF indices - ERR137249.vcf.gz.csi and SRR1103491.vcf.gz.csi
  * Merged.vcf
  * Merged.tsv

Printed Output -
```
ISOLATE PREDICTION   	DIFF RI
1 ERR137249  		S 0.04866350  0
2 SRR1103491  		M 0.02846038  0
```

### 4. Prediction from merged VCFs

  Step 1. Create an Input directory

    mkdir /home/username/Input_Dir5

  Step 2. Store a merged VCF file generated by merging multiple VCF files from variant calling of  MTB isolates based on MTB H37Rv reference genome  (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz)

  Step 3. Create an Output directory

    mkdir /home/username/Output_Dir5

  Step 4. Go to the mcdr-mtb-standalone-v2 installation directory

    cd /home/username/Documents/mcdr-mtb-standalone-v2/

  Step 5. Run mcdr-merge-predict.sh

    ./mcdr-merge-predict.sh /home/username/Input_Dir5/ /home/username/Output_Dir5/

`Input_Dir5` contains merged.vcf
`Output_Dir5` contains Merged.tsv

Printed Output -
```
ISOLATE PREDICTION   	DIFF RI
1 ERR137249  		S 0.04866350  0
2 SRR1103491  		M 0.02846038  0
```
## Team

**Abhirupa Ghosh, Sudipto Bhattacharjee and Sudipto Saha**

## Disclaimer

The scripts were tried and tested on the Ubuntu Operating system.

This tool is strictly for Research Use Only. By using this tool the user acknowledges no intended medical purpose such as patient diagnosis.
