#!/bin/bash

DEST=/hdfs/store/user/mithakor/2024_Categorized_MTT

# List of input dirs
datasets=(
JetMET0_Run2024C-MINIv6NANOv15-v1_NANOAOD_18Sep25_1713
JetMET0_Run2024D-MINIv6NANOv15-v1_NANOAOD_18Sep25_1715
JetMET0_Run2024E-MINIv6NANOv15-v1_NANOAOD_18Sep25_1716
JetMET0_Run2024F-MINIv6NANOv15-v2_NANOAOD_18Sep25_1718
JetMET0_Run2024G-MINIv6NANOv15-v2_NANOAOD_18Sep25_1720
JetMET0_Run2024H-MINIv6NANOv15-v2_NANOAOD_18Sep25_1723
JetMET0_Run2024I-MINIv6NANOv15-v2_NANOAOD_18Sep25_1724
JetMET0_Run2024I-MINIv6NANOv15_v2-v1_NANOAOD_18Sep25_1726
JetMET1_Run2024C-MINIv6NANOv15-v1_NANOAOD_18Sep25_1728
JetMET1_Run2024D-MINIv6NANOv15-v1_NANOAOD_18Sep25_1730
JetMET1_Run2024E-MINIv6NANOv15-v1_NANOAOD_18Sep25_1732
JetMET1_Run2024F-MINIv6NANOv15-v2_NANOAOD_18Sep25_1734
JetMET1_Run2024G-MINIv6NANOv15-v2_NANOAOD_18Sep25_1739
JetMET1_Run2024H-MINIv6NANOv15-v2_NANOAOD_18Sep25_1743
JetMET1_Run2024I-MINIv6NANOv15-v1_NANOAOD_18Sep25_1744
JetMET1_Run2024I-MINIv6NANOv15_v2-v2_NANOAOD_18Sep25_1745
)

# Loop over datasets
for ds in "${datasets[@]}"; do
  # Strip the trailing _<date>_<time>_<id>
  base_jobname=$(echo "$ds" | sed -E 's/_18Sep25_[0-9]+$//')

  echo ">>> Submitting job for $ds (jobName=$base_jobname)"
  python3 submit_jobs.py \
    --inputDir /hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/$ds \
    --destination $DEST     \
    --jobName $base_jobname  \
    --submitDirPath /nfs_scratch/mithakor/Corrections_MTT   \
    --year 2024 \

done

echo "Submitted all."