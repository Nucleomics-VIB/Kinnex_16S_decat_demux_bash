#!/usr/bin/env bash

# script: Kinnex_16S_decat_demux.sh
# run skera and lima on a Revio multiplexed Kinnex 16S RUN
#
# Stephane Plaisance - VIB-NC 2024-06-03 v1.0
# small edits in the header below: v1.0.1
# rewritten to improve file structure: v1.1.0
# added samplesheet validation: v1.2.2
# added post_processing: v2.0.0
# added archive creation + md5sum;: v2.1.0
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

version="2025-06-24; 2.1.0"

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

echo -e "\n"

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

function post_processing() {
local flag_file="${outfolder}/post_processing_ok"

echo -e "\n# Post-processing: preparing outputs for delivery"

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo -e "\npost_processing: already done."
    return 0 # Exit the function successfully
fi

# Create run_QC directory for delivery outputs
local qc_dir="${outfolder}/run_QC"
mkdir -p "${qc_dir}"

echo "# Collecting QC files from lima results"

# Find all relevant files in lima_results (excluding flag files and report files)
find "${outfolder}/${lima_results}" -maxdepth 2 -type f \
    -not -name "Lima_ok" \
    -not -name "HiFi.lima.report" \
    -not -name "*-log.txt" \
    -not -name "*.xml" \
    -not -name "*.json" | while read -r file; do
    
    # Extract barcode pattern (bc0X) from the file path
    barcode_pattern=$(echo "$file" | grep -o 'bc0[0-4]')
    
    # If barcode pattern is found, copy file with renamed format
    if [ -n "$barcode_pattern" ]; then
        filename=$(basename "$file")
        destination="${qc_dir}/${barcode_pattern}_${filename}"
        
        echo "  Copying $(basename "$file") -> ${barcode_pattern}_${filename}"
        cp "$file" "$destination"
    fi
done

echo "# Collecting CSV samplesheet files"

# Copy CSV samplesheet files used for lima
for i in "${barcode_indices[@]}"; do
    samplesheet_var="bcM000${i}_samplesheet"
    samplesheet_file="${outfolder}/${inputs}/${!samplesheet_var}"
    
    if [ -f "$samplesheet_file" ]; then
        filename=$(basename "$samplesheet_file")
        destination="${qc_dir}/bc0${i}_${filename}"
        
        echo "  Copying samplesheet $(basename "$samplesheet_file") -> bc0${i}_${filename}"
        cp "$samplesheet_file" "$destination"
    else
        echo "  Warning: Samplesheet file not found: $samplesheet_file"
    fi
done

echo "# QC files prepared in: ${qc_dir}"

# Check for and copy movie report PDF if it exists
local movie_report="${outfolder}/${movie}.report.pdf"
if [ -f "$movie_report" ]; then
    echo "# Copying movie report PDF"
    cp "$movie_report" "${qc_dir}/"
    echo "  Copied ${movie}.report.pdf to run_QC"
else
    echo "# Movie report PDF not found: ${movie}.report.pdf"
fi

# Merge FASTQ results if multiple barcode subfolders exist
merge_fastq_results

# Create delivery archive
create_archive

echo "# Post-processing completed successfully"

# Write the flag file upon successful completion
touch "${flag_file}"
}

function merge_fastq_results() {
local fastq_dir="${outfolder}/${fastq_results}"
local final_dir="${outfolder}/fastq_final"
local flag_file="${final_dir}/merge_fastq_results_ok"

echo -e "\n# Checking barcode subfolders in fastq results"

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo -e "\nmerge_fastq_results: already done."
    return 0 # Exit the function successfully
fi

# Count the number of bc subfolders using mapfile
mapfile -t bc_folders < <(find "$fastq_dir" -maxdepth 1 -type d -name "bc0*" 2>/dev/null)
local bc_count=${#bc_folders[@]}

echo "# Found $bc_count barcode subfolder(s) in fastq results"

if [ "$bc_count" -eq 0 ]; then
    echo "# No barcode subfolders found, skipping processing"
    return 0
elif [ "$bc_count" -eq 1 ]; then
    echo "# Single barcode subfolder found, copying files directly to fastq_final"
    mkdir -p "$final_dir"
    
    # Copy all FASTQ files from the single bc folder directly to fastq_final
    local single_bc_folder="${bc_folders[0]}"
    while IFS= read -r -d '' file; do
        filename=$(basename "$file")
        echo "  Copying $filename"
        cp "$file" "$final_dir/$filename"
    done < <(find "$single_bc_folder" -type f \( -name "*.fastq*" -o -name "*.fq*" \) -print0 2>/dev/null)
    
    echo "# FASTQ files copied to: $final_dir"
    
    # Write the flag file upon successful completion
    touch "${flag_file}"
    return 0
fi

echo "# Multiple barcode subfolders detected, merging FASTQ files"
mkdir -p "$final_dir"

# Create associative array to track files by name
declare -A file_list

# First pass: collect all unique filenames
for bc_folder in "${bc_folders[@]}"; do
    if [ -d "$bc_folder" ]; then
        while IFS= read -r -d '' file; do
            filename=$(basename "$file")
            if [[ "$filename" =~ \.(fastq|fq)(\.gz)?$ ]]; then
                file_list["$filename"]=1
            fi
        done < <(find "$bc_folder" -type f \( -name "*.fastq*" -o -name "*.fq*" \) -print0 2>/dev/null)
    fi
done

# Second pass: merge files with same names
for filename in "${!file_list[@]}"; do
    local output_file="$final_dir/$filename"
    local temp_files=()
    
    # Collect all files with this name from different bc folders
    for bc_folder in "${bc_folders[@]}"; do
        local source_file="$bc_folder/$filename"
        if [ -f "$source_file" ]; then
            temp_files+=("$source_file")
        fi
    done
    
    if [ ${#temp_files[@]} -eq 0 ]; then
        continue
    elif [ ${#temp_files[@]} -eq 1 ]; then
        # Only one file with this name, just copy it
        echo "  Copying $filename (single file)"
        cp "${temp_files[0]}" "$output_file"
    else
        # Multiple files with same name, concatenate them
        echo "  Merging $filename (${#temp_files[@]} files)"
        
        # Handle compressed vs uncompressed files
        if [[ "$filename" =~ \.gz$ ]]; then
            # Compressed files - concatenate compressed data
            cat "${temp_files[@]}" > "$output_file"
        else
            # Uncompressed files - simple concatenation
            cat "${temp_files[@]}" > "$output_file"
        fi
    fi
done

echo "# FASTQ files processed successfully in: $final_dir"
echo "# Total unique files processed: ${#file_list[@]}"

# Write the flag file upon successful completion
touch "${flag_file}"
}

function create_archive() {
local flag_file="${outfolder}/create_archive_ok"
local archive_name
local archive_path

archive_name="${movie}.tar.gz"
archive_path="${outfolder}/${archive_name}"

echo -e "\n# Creating delivery archive"

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo -e "\ncreate_archive: already done."
    return 0 # Exit the function successfully
fi

# Check if required directories exist
local qc_dir="${outfolder}/run_QC"
local final_dir="${outfolder}/fastq_final"

if [ ! -d "$qc_dir" ] && [ ! -d "$final_dir" ]; then
    echo "# Neither run_QC nor fastq_final directories exist, skipping archive creation"
    return 0
fi

echo "# Creating archive: $archive_name"

# Change to output directory to create relative paths in the archive
cd "${outfolder}" || {
    echo "Error: Cannot change to output directory: ${outfolder}"
    return 1
}

# Build tar command with existing directories
local tar_args=()
if [ -d "run_QC" ]; then
    tar_args+=("run_QC")
    echo "# Including run_QC directory in archive"
fi

if [ -d "fastq_final" ]; then
    tar_args+=("fastq_final")
    echo "# Including fastq_final directory in archive"
fi

# Create the archive
if [ ${#tar_args[@]} -gt 0 ]; then
    # Join array elements with commas for display
    local dirs_list
    dirs_list=$(printf "%s, " "${tar_args[@]}")
    dirs_list=${dirs_list%, }  # Remove trailing comma and space
    echo "# Creating archive with directories: ${dirs_list}"
    echo "# Creating archive and computing md5sum simultaneously"
    
    # Create archive and compute md5sum in one pass using tee and md5sum
    if tar -cz "${tar_args[@]}" | tee "${archive_name}" | md5sum | sed "s/-$/${archive_name}/" > "${archive_name}.md5"; then
        echo "# Archive created successfully: ${archive_path}"
        echo "# Archive size: $(du -h "${archive_name}" | cut -f1)"
        echo "# MD5 checksum saved to: ${archive_name}.md5"
        echo "# MD5: $(cat "${archive_name}.md5")"
    else
        echo "Error: Failed to create archive or compute md5sum"
        return 1
    fi
else
    echo "# No directories to archive"
    return 0
fi

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
time post_processing || { echo "post_processing failed"; exit 1; }
time create_archive || { echo "create_archive failed"; exit 1; }

exit 0
