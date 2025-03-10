#!/usr/bin/env bash

# script: Kinnex_16S_decat_demux.sh
# run skera and lima on a Revio multiplexed Kinnex 16S RUN
#
# Stephane Plaisance - VIB-NC 2024-06-03 v1.0
# small edits in the header below: v1.0.1
# rewritten to improve file structure: v1.1.0
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

version="2025-03-07; 1.2.0"

# script basedir
BASEDIR=$(dirname "$(readlink -f "$0")")

myenv="Kinnex_16S_decat_demux_env"

################################
########## FUNCTIONS ###########
################################

function parse_yaml() {
    local prefix=$2
    local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
    sed -ne "s|^\($s\):|\1|" \
         -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
         -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p" $1 |
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
echo "bcM0001 BAM: $bcM0001_bam"
echo "bcM0001 Samplesheet: $bcM0001_samplesheet"
echo "bcM0002 BAM: $bcM0002_bam"
echo "bcM0002 Samplesheet: $bcM0002_samplesheet"
echo "bcM0003 BAM: $bcM0003_bam"
echo "bcM0003 Samplesheet: $bcM0003_samplesheet"
echo "bcM0004 BAM: $bcM0004_bam"
echo "bcM0004 Samplesheet: $bcM0004_samplesheet"
echo "- outfolder: $outfolder"
echo "- inputs: $inputs"
echo "- skera_results: $skera_results"
echo "- lima_results: $lima_results"
echo "- fastq_results: $fastq_results"
echo "- log_skera: $log_skera"
echo "- log_lima: $log_lima"
echo "- nthr_skera: $nthr_skera"
echo "- nthr_lima: $nthr_lima"
echo "- par_bam2fastq: $par_bam2fastq"
echo "- nthr_bam2fastq: $nthr_bam2fastq"
echo "- mincnt: $mincnt"
echo "- qc_format: $qc_format"
echo "- final_results: $final_results"
}

function CopyRunData() {
local flag_file="${outfolder}/${inputs}/CopyRunData_ok"

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo "CopyRunData: already done."
    return 0 # Exit the function successfully
fi

mkdir -p "${outfolder}/${inputs}"

echo -e "\n# Copying RUN data locally"

# Check if the adapter folder exists
if [ -d "${runfolder}/${hififolder}" ]; then
    cp ${runfolder}/${hififolder}/*bc*.bam* ${outfolder}/${inputs}/
else
    echo "HiFi folder not found: ${runfolder}/${hififolder}"
    return 1 # Exit the function with an error status
fi

# Check if the sample sheet exists
for i in {1..4}; do
  samplesheet_var="bcM000${i}_samplesheet"
  samplesheet_file="${!samplesheet_var}"
  
  if [ -f "$runfolder/$samplesheet_file" ]; then
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

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo "SkeraSplit: already done."
    return 0 # Exit the function successfully
fi

mkdir -p "${outfolder}/${skera_results}/bc0"{1..4}

echo -e "\n# Running Skera de-concatenation"

# run for each bc bam file and samplesheet
for i in {1..4}; do
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

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo "Lima: already done."
    return 0 # Exit the function successfully
fi

mkdir -p "${outfolder}/${lima_results}/bc0"{1..4}

echo -e "\n# Running Lima demultiplexing"

for i in {1..4}; do
  bam_file="${outfolder}/${skera_results}/bc0${i}/skera.bam"
  samplesheet_var="bcM000${i}_samplesheet"
  samplesheet_file="${outfolder}/${inputs}/${!samplesheet_var}"
  
  if [[ -f "${samplesheet_file}" && -f "${bam_file}" ]]; then

    cmd="$(which lima) \
      ${bam_file} \
      ${BASEDIR}/barcode_files/Kinnex16S_384plex_primers/Kinnex16S_384plex_primers.fasta \
      ${outfolder}/${lima_results}/bc0${i}/HiFi.bam \
      --hifi-preset ASYMMETRIC \
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
    eval ${cmd} && mv barcode_QC_Kinnex.* ${outfolder}/${lima_results}/bc0${i}/

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

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo "bam2fastq: already done."
    return 0 # Exit the function successfully
fi

mkdir -p "${outfolder}/${fastq_results}/bc0"{1..4}

for i in {1..4}; do
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
    bcpair=${pfx#HiFi.}
    biosample=$(grep ${bcpair} ${samplesheet_file} | \
      dos2unix | cut -d, -f 2 | tr -d "\n" | tr ' ' '_')

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

# Process command-line arguments
while getopts "c:" opt; do
  case "$opt" in
    c)
      CONFIG_FILE="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
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

# create output folder
mkdir -p "${outfolder}"

# redirect all outputs to a log file
cat /dev/null > "${outfolder}/runlog.txt"
exec &> >(tee -a "${outfolder}/runlog.txt")

# Print configuration
PrintConfig

# Run the pipeline steps
time CopyRunData

time SkeraSplit

time Lima

time bam2fastq

exit 0









#########################
####### leftovers #######
#########################

#time BundleResults

#time createArchive

function BundleResults() {
local flag_file="${outfolder}/BundleResults_ok"

# Check if the flag file exists and echo "already done" if it does
if [ -f "${flag_file}" ]; then
    echo "BundleResults: already done."
    return 0 # Exit the function successfully
fi

mkdir -p "${outfolder}/${final_results}"

echo -e "\n# Copying files to final_result folder for transfer"

# copy run_QC from RUN folder
cp ${runfolder}/*.pdf ${outfolder}/${final_results}/

# copy Zymo control PDF anad README.txt
cp -r info ${outfolder}/${final_results}/

cp ${outfolder}/${inputs}/${samplesheet} ${outfolder}/${final_results}/
cp ${outfolder}/${skera_results}/${movie}.skera.summary.csv ${outfolder}/${final_results}/skera.summary.csv
cp ${outfolder}/${lima_results}/HiFi.lima.* ${outfolder}/${final_results}/

# move fastq to save room
mv ${outfolder}/${fastq_results} ${outfolder}/${final_results}/

echo -e "# Creating barcode QC report"
projectnum=$(echo ${samplesheet} | cut -d "_" -f 1 | tr -d "\n")
cmd="${BASEDIR}/scripts/barcode_QC_Kinnex.sh \
  -i ${outfolder}/${final_results}/HiFi.lima.counts \
  -r scripts/barcode_QC_Kinnex.Rmd \
  -m ${mincnt} \
  -f ${qc_format} \
  -p ${projectnum} \
  -s ${outfolder}/${inputs}/${samplesheet}"

echo "# ${cmd}"
eval ${cmd} && mv barcode_QC_Kinnex.${qc_format} ${outfolder}/${final_results}/

# Write the flag file upon successful completion
touch "$flag_file"
}

function createArchive() {
local flag_file="Archive_ok"

# Check if the flag file exists and echo "already done" if it does
if [ -f "# create barcode plots
${flag_file}" ]; then
    echo "createArchive: already done."
    return 0 # Exit the function successfully
fi

echo -e "\n# Creating TGZ archive of ${final_results} and its md5sum"

thr=8
pfx="$(echo ${samplesheet} | cut -d '_' -f 1 | tr -d '\n')_archive"

cd ${outfolder}
{ tar cvf - "${final_results}" \
  | pigz -p ${thr} \
  | tee >(md5sum > ${pfx}.tgz_md5.txt) > ${pfx}.tgz; \
  } 2> ${pfx}_content.log

echo -e "# archive ${pfx}.tgz and md5sum ${pfx}.tgz_md5.txt were created"
# fix file path in md5sum
sed -i "s|-|${pfx}.tgz|g" ${pfx}.tgz_md5.txt

echo -e "\n# Checking md5sum"
md5sum -c ${pfx}.tgz_md5.txt | tee -a ${pfx}.tgz_md5-check.txt

# write flag if checksum is OK
if grep -q "OK" "${pfx}.tgz_md5-check.txt"; then
    # Write the flag file upon successful completion
    touch "$flag_file"
else
    echo "Flag file not created. Verification failed."
fi
}