#!/usr/bin/env bash

# script: Kinnex_16S_decat_demux.sh
# run skera and lima on a Revio multiplexed Kinnex 16S RUN
#
# Stephane Plaisance - VIB-NC 2024-06-03 v1.0
# small edits in the header below: v1.0.1
# rewritten to improve file structure: v1.1.0
# added samplesheet validation: v1.2.2
#
# visit our Git: https://github.com/Nucleomics-VIB

# requirements:
# conda env: Kinnex_16S_decat_demux_env (installed from conda_env_setup.yaml)
# => conda env create -f conda_env_setup.yaml
# config.yaml (edited and pointing to existing files)
# yq (version 4+ to read the yaml config file into bash variables)
# PacBio barcode files copied to barcode_files)
# a BAM file resulting from a Kinnex 16S run
# a samplesheet linking barcode pairs to sample names (custom made)
# All parameters have been externalised from the code and are listed in config.yaml

# improve error handling and report
set -euo pipefail
IFS=$'\n\t'

version="2025-06-24; 1.2.2"

# script basedir
BASEDIR=$(dirname "$(readlink -f "$0")")

myenv="Kinnex_16S_decat_demux_env"

# prefix used to extract barcode pair from BAM filenames
bam_prefix="HiFi."

################################
########## FUNCTIONS ###########
################################

function parse_yaml() {
    local yaml_file=$1
    local prefix="${2:-}"
    local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
    sed -ne "s|^\($s\):|\1|" \
         -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
         -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p" "$yaml_file" |
    awk -F$fs '{
        indent = length($1)/2;
        vname[indent] = $2;
        for (i in vname) {if (i > indent) {delete vname[i]}}
        if (length($3) > 0) {
            vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
            printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
        }
    }'
}

function PrintConfig() {
echo "####################"
echo "# run configuration"
echo "#"
echo "- runfolder: $runfolder"
echo "- movie: $movie"
echo "- hififolder: $hififolder"
echo "- barcode_number: $bcnum"
for i in "${barcode_indices[@]}"; do
  bam_var="bcM000${i}_bam"
  samplesheet_var="bcM000${i}_samplesheet"
  echo "- bcM000${i} BAM: ${!bam_var}"
  echo "- bcM000${i} Samplesheet: ${!samplesheet_var}"
done
echo "- outfolder: $outfolder"
echo "- inputs: $inputs"
echo "- skera_results: $skera_results"
echo "- lima_results: $lima_results"
echo "- fastq_results: $fastq_results"
echo "- lima_min_len: $lima_min_len"
echo "- log_skera: $log_skera"
echo "- log_lima: $log_lima"
echo "- nthr_skera: $nthr_skera"
echo "- nthr_lima: $nthr_lima"
echo "- par_bam2fastq: $par_bam2fastq"
echo "- nthr_bam2fastq: $nthr_bam2fastq"
echo "- mincnt: $mincnt"
echo "- qc_format: $qc_format"
echo "- final_results: $final_results"
echo "reference files:"
echo "- ${BASEDIR}/barcode_files/MAS-Seq_Adapter_v2/mas12_primers.fasta"
echo "- ${BASEDIR}/barcode_files/Kinnex16S_384plex_primers/Kinnex16S_384plex_primers.fasta"
}

function check_csv() {
local csv_file="$1"

echo -e "\n# Validating CSV samplesheet: $(basename $csv_file)"

# Check if file exists
if [ ! -f "$csv_file" ]; then
    echo "Error: CSV file not found: $csv_file"
    return 1
fi

# Convert to Unix format and remove any BOM
local temp_csv=$(mktemp)
dos2unix < "$csv_file" > "$temp_csv" 2>/dev/null

# Check for Windows line endings (CR+LF) or Mac line endings (CR only) in original file
if grep -q $'\r' "$csv_file"; then
    echo "Error: CSV file contains carriage return characters (^M). Please use Unix line endings (LF only)"
    rm "$temp_csv"
    return 1
fi

# Check header (first line)
local header=$(head -n 1 "$temp_csv")
if [ "$header" != "Barcode,Bio Sample" ]; then
    echo "Error: CSV header must be exactly 'Barcode,Bio Sample', found: '$header'"
    rm "$temp_csv"
    return 1
fi

# Arrays to track duplicates
local -a barcodes_seen=()
local -a biosamples_seen=()

# Check each data row (skip header)
local line_num=1
while IFS= read -r line; do
    ((line_num++))
    
    # Skip empty lines
    [ -z "$line" ] && continue
    
    # Count columns (commas + 1)
    local col_count=$(echo "$line" | tr -cd ',' | wc -c)
    ((col_count++))
    
    if [ "$col_count" -ne 2 ]; then
        echo "Error: Line $line_num must have exactly 2 columns, found $col_count: '$line'"
        rm "$temp_csv"
        return 1
    fi
    
    # Extract both columns
    local barcode=$(echo "$line" | cut -d, -f1)
    local bio_sample=$(echo "$line" | cut -d, -f2)
    
    # Check if Bio Sample is empty
    if [ -z "$bio_sample" ]; then
        echo "Error: Line $line_num has empty Bio Sample column: '$line'"
        rm "$temp_csv"
        return 1
    fi
    
    # Check if Bio Sample contains only valid characters (a-z, A-Z, 0-9, -, _, .)
    if [[ ! "$bio_sample" =~ ^[a-zA-Z0-9_.-]+$ ]]; then
        echo "Error: Line $line_num Bio Sample '$bio_sample' contains invalid characters. Only a-z, A-Z, 0-9, -, _, . are allowed"
        rm "$temp_csv"
        return 1
    fi
    
    # Check for duplicate barcodes
    for seen_barcode in "${barcodes_seen[@]}"; do
        if [ "$barcode" = "$seen_barcode" ]; then
            echo "Error: Line $line_num has duplicate barcode '$barcode'"
            rm "$temp_csv"
            return 1
        fi
    done
    barcodes_seen+=("$barcode")
    
    # Check for duplicate bio samples
    for seen_biosample in "${biosamples_seen[@]}"; do
        if [ "$bio_sample" = "$seen_biosample" ]; then
            echo "Error: Line $line_num has duplicate Bio Sample name '$bio_sample'"
            rm "$temp_csv"
            return 1
        fi
    done
    biosamples_seen+=("$bio_sample")
    
done < <(tail -n +2 "$temp_csv")

rm "$temp_csv"
echo "CSV validation passed: $(basename $csv_file)"
return 0
}

function CopyRunData() {
local flag_file="${outfolder}/${inputs}/CopyRunData_ok"

echo -e "\n# Copying RUN data locally"

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo -e "\nCopyRunData: already done."
    return 0 # Exit the function successfully
fi

mkdir -p "${outfolder}/${inputs}"

# Check if the adapter folder exists
if [ -d "${runfolder}/${hififolder}" ]; then
    cp "${runfolder}/${hififolder}"/*bc*.bam* "${outfolder}/${inputs}/"
else
    echo "HiFi folder not found: ${runfolder}/${hififolder}"
    return 1 # Exit the function with an error status
fi

# Check if the sample sheet exists
for i in "${barcode_indices[@]}"; do
  samplesheet_var="bcM000${i}_samplesheet"
  samplesheet_file="${!samplesheet_var}"
  
  if [ -f "$runfolder/$samplesheet_file" ]; then
    # Validate CSV format before copying
    if ! check_csv "$runfolder/$samplesheet_file"; then
      echo "CSV validation failed for $samplesheet_file"
      return 1
    fi
    
    cp "${runfolder}/$samplesheet_file" "${outfolder}/${inputs}/"
    echo "Copied $samplesheet_file to ${outfolder}/${inputs}/"
  else
    echo "File $samplesheet_file does not exist in $runfolder"
    return 1 # Exit the function with an error status
  fi
done

# Write the flag file upon successful completion
touch "${flag_file}"
}

function SkeraSplit {
local flag_file="${outfolder}/${skera_results}/SkeraSplit_ok"

echo -e "\n# Running Skera de-concatenation"

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo -e "\nSkeraSplit: already done."
    return 0 # Exit the function successfully
fi

for i in "${barcode_indices[@]}"; do
mkdir -p "${outfolder}/${skera_results}/bc0${i}"
done

# run for each bc bam file and samplesheet
for i in "${barcode_indices[@]}"; do
  bam_var="${movie}.hifi_reads.bcM000${i}.bam"
  bam_file="${outfolder}/${inputs}/${bam_var}"

  if [[ -f "$bam_file" ]]; then
    cmd="$(which skera) split \
      ${bam_file} \
      ${BASEDIR}/barcode_files/MAS-Seq_Adapter_v2/mas12_primers.fasta \
      ${outfolder}/${skera_results}/bc0${i}/skera.bam \
      --num-threads ${nthr_skera} \
      --log-level ${log_skera} \
      --log-file ${outfolder}/${skera_results}/skera_run_bc0${i}-log.txt"

    echo "# ${cmd}"
    eval ${cmd}

  else
    echo "One or both files are missing"
    return 1 # Exit the function with an error status
  fi
  
done

# Write the flag file upon successful completion
touch "${flag_file}"
}

function Lima() {
local flag_file="${outfolder}/${lima_results}/Lima_ok"

echo -e "\n# Running Lima demultiplexing"

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo -e "\nLima: already done."
    return 0 # Exit the function successfully
fi

for i in "${barcode_indices[@]}"; do
mkdir -p "${outfolder}/${lima_results}/bc0${i}"
done

for i in "${barcode_indices[@]}"; do
  bam_file="${outfolder}/${skera_results}/bc0${i}/skera.bam"
  samplesheet_var="bcM000${i}_samplesheet"
  samplesheet_file="${outfolder}/${inputs}/${!samplesheet_var}"
  
  if [[ -f "${samplesheet_file}" && -f "${bam_file}" ]]; then

    cmd="$(which lima) \
      ${bam_file} \
      ${BASEDIR}/barcode_files/Kinnex16S_384plex_primers/Kinnex16S_384plex_primers.fasta \
      ${outfolder}/${lima_results}/bc0${i}/HiFi.bam \
      --hifi-preset ASYMMETRIC \
      --min-length ${lima_min_len} \
      --split-named \
      --split-subdirs \
      --biosample-csv ${samplesheet_file} \
      --num-threads ${nthr_lima} \
      --log-level ${log_lima} \
      --log-file ${outfolder}/${lima_results}/lima_run_bc0${i}-log.txt"

    echo "# ${cmd}"
    eval ${cmd}

    echo -e "# Creating barcode QC report"
    projectnum=$(basename ${samplesheet_file} | cut -d "_" -f 1 | tr -d "\n")
    cmd="${BASEDIR}/scripts/barcode_QC_Kinnex.sh \
      -i ${outfolder}/${lima_results}/bc0${i}/HiFi.lima.counts \
      -r ${BASEDIR}/scripts/barcode_QC_Kinnex.Rmd \
      -m ${mincnt} \
      -f ${qc_format} \
      -p ${projectnum} \
      -s ${samplesheet_file}"
    
    echo "# ${cmd}"
    if eval ${cmd}; then
      # Move QC files if they exist
      if ls barcode_QC_Kinnex.* 1> /dev/null 2>&1; then
        mv barcode_QC_Kinnex.* "${outfolder}/${lima_results}/bc0${i}/"
      else
        echo "Warning: No barcode_QC_Kinnex.* files found to move"
      fi
    else
      echo "Error: barcode_QC_Kinnex.sh failed"
      return 1
    fi

  else
    echo "One or both files are missing"
    return 1 # Exit the function with an error status
  fi

done

# Write the flag file upon successful completion
touch "${flag_file}"
}

function bam2fastq() {
local flag_file="${outfolder}/${fastq_results}/bam2fastq_ok"

echo -e "\n# Converting BAM data to FastQ and renaming samples"

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo -e "\nbam2fastq: already done."
    return 0 # Exit the function successfully
fi

for i in "${barcode_indices[@]}"; do
mkdir -p "${outfolder}/${fastq_results}/bc0${i}"
done

for i in "${barcode_indices[@]}"; do
  lima_folder="${outfolder}/${lima_results}/bc0${i}"
  samplesheet_var="bcM000${i}_samplesheet"
  samplesheet_file="${outfolder}/${inputs}/${!samplesheet_var}"
  joblist="${outfolder}/${fastq_results}/bc0${i}/job.list"
  
  if [[ -f "${samplesheet_file}" && -d "${lima_folder}" ]]; then

    # initialize
    cat /dev/null > ${joblist}
    
    echo -e "\n# Preparing job list from all bc0${i} lima BAM files"
    
    for bam in $(find ${lima_folder} -name "*.bam"); do
    # rename sample from samplesheet 'Bio Sample'
    pfx=$(basename ${bam%.bam})
    bcpair=${pfx#${bam_prefix}}
    
    # Check if barcode pair exists in samplesheet
    if ! grep -q "${bcpair}" "${samplesheet_file}"; then
      echo "Error: Sample name not found in the samplesheet for barcode pair ${bcpair}"
      return 1
    fi
    
    biosample=$(grep ${bcpair} ${samplesheet_file} | \
      dos2unix | cut -d, -f 2 | tr -d "\n" | tr ' ' '_')
    
    # Verify that biosample is not empty and contains valid characters
    if [[ -z "${biosample}" ]]; then
      echo "Error: Sample name not found in the samplesheet for barcode pair ${bcpair}"
      return 1
    fi
    
    # Check if biosample contains only valid filename characters (alphanumeric, underscore, hyphen, dot)
    if [[ ! "${biosample}" =~ ^[a-zA-Z0-9_.-]+$ ]]; then
      echo "Error: Sample name not found in the samplesheet for barcode pair ${bcpair}"
      return 1
    fi

    echo "$(which bam2fastq) \
        ${bam} \
        --output ${outfolder}/${fastq_results}/bc0${i}/${biosample} \
        --num-threads ${nthr_bam2fastq}" >> ${joblist}
    done

    # execute job list in batches of \${nthr_bam2fastq_par}
    cmd="parallel -j ${par_bam2fastq} --joblog my_job_log.log \
      < ${joblist} && (rm ${joblist} my_job_log.log)"

    echo -e "\n# Executing bc0${i} ${joblist} in parallel batches"

    echo "# ${cmd}"
    eval ${cmd}
  else
    echo "One or both inputs are missing"
    return 1 # Exit the function with an error status
  fi
done

# Write the flag file upon successful completion
touch "${flag_file}"
}



##################################
##########  MAIN SCRIPT ##########
##################################

usage() {
  cat << EOF
Usage: $(basename "$0") -c <config.yaml>

Options:
  -c <config.yaml>   Path to the YAML configuration file (required).

Description:
  This script demultiplexes PacBio Kinnex 16S sequencing data using skera and lima.
  All parameters are externalized in the config.yaml file.

Example:
  $(basename "$0") -c my_config.yaml

EOF
}

# Process command-line arguments
while getopts "c:h" opt; do
  case "$opt" in
    c)
      CONFIG_FILE="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
    h)
      usage
      exit 0
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      exit 1
      ;;
  esac
done

# Check if the config file argument was provided
if [ -z "$CONFIG_FILE" ]; then
  echo "Error: Configuration file is required. Please use the -c option."
  usage
  exit 1
fi

# Check if the specified config file exists
if [ ! -f "$CONFIG_FILE" ]; then
  echo "Error: Configuration file '$CONFIG_FILE' not found."
  exit 1
fi

# Activate conda environment
source /etc/profile.d/conda.sh
conda activate "${myenv}" || {
  echo "# the conda environment ${myenv} was not found on this machine"
  echo "# please read the top part of the script!"
  exit 1
}

# Check if required executables are present
declare -a arr=( "yq" "skera" "lima" "bam2fastq" "pigz" "scripts/barcode_QC_Kinnex.sh" )
for prog in "${arr[@]}"; do
  hash "$prog" 2>/dev/null || { echo "# required ${prog} not found in PATH"; exit 1; }
done

# Read config in and prepare data and folders
eval "$(parse_yaml "$CONFIG_FILE")"

# --- Barcode detection logic ---
declare -a barcode_indices=()
for var in $(compgen -A variable | grep -E '^bcM[0-9]{4}_bam$'); do
  bam_val="${!var}"
  idx=$(echo "$var" | sed -E 's/^bcM0*([0-9]+)_bam$/\1/')
  base="bcM$(printf "%04d" "$idx")"
  samplesheet_var="${base}_samplesheet"
  samplesheet_val="${!samplesheet_var}"
  if [[ -n "$bam_val" || -n "$samplesheet_val" ]]; then
    barcode_indices+=("$idx")
  fi
done
bcnum=${#barcode_indices[@]}

echo "# found ${bcnum} barcoded bam HiFi files"

# create output folder
mkdir -p "${outfolder}"

# redirect all outputs to a log file
cat /dev/null > "${outfolder}/runlog.txt"
exec &> >(tee -a "${outfolder}/runlog.txt")

# Print configuration
PrintConfig

# Run the pipeline steps
time CopyRunData || { echo "CopyRunData failed"; exit 1; }
time SkeraSplit   || { echo "SkeraSplit failed"; exit 1; }
time Lima         || { echo "Lima failed"; exit 1; }
time bam2fastq    || { echo "bam2fastq failed"; exit 1; }

exit 0

########## future additions

copy the lima simmaries to data_transfer

#!/bin/bash

# List of barcodes
barcodes=(bc01 bc02 bc03 bc04)

# Output directory
outdir="data_transfer/demux_results"
mkdir -p "$outdir"

# Loop over each barcode
for bc in "${barcodes[@]}"; do
    for file in HiFi.lima.counts HiFi.lima.summary barcode_QC_Kinnex.html; do
        src="lima_results/${bc}/${file}"
        if [[ -f "$src" ]]; then
            # Insert barcode before extension
            base="${file%.*}"
            ext="${file##*.}"
            # Handle files with no extension
            if [[ "$base" == "$file" ]]; then
                newname="${base}_${bc}"
            else
                newname="${base}_${bc}.${ext}"
            fi
            cp "$src" "$outdir/$newname"
        else
            echo "Warning: $src does not exist."
        fi
    done
done
