
Kinnex_16S_decat_demux.sh






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
    echo -e  "\ncreateArchive: already done."
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