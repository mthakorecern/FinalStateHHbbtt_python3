#!/bin/bash

# Common options
DEST=/hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/WithFatJetID_FatJetVetoMaps
JSON_BASE=/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024

# List of input dirs
datasets=(
JetMET0_Run2024C-MINIv6NANOv15-v1_NANOAOD_10Sep25_1352_000
JetMET0_Run2024D-MINIv6NANOv15-v1_NANOAOD_10Sep25_1518_002
JetMET0_Run2024E-MINIv6NANOv15-v1_NANOAOD_10Sep25_1521_004
JetMET0_Run2024F-MINIv6NANOv15-v2_NANOAOD_10Sep25_1524_006
JetMET0_Run2024G-MINIv6NANOv15-v2_NANOAOD_10Sep25_1528_008
JetMET0_Run2024H-MINIv6NANOv15-v2_NANOAOD_10Sep25_1534_010
JetMET0_Run2024I-MINIv6NANOv15-v2_NANOAOD_10Sep25_1536_012
JetMET0_Run2024I-MINIv6NANOv15_v2-v1_NANOAOD_10Sep25_1537_013
JetMET1_Run2024C-MINIv6NANOv15-v1_NANOAOD_10Sep25_1516_001
JetMET1_Run2024D-MINIv6NANOv15-v1_NANOAOD_10Sep25_1520_003
JetMET1_Run2024E-MINIv6NANOv15-v1_NANOAOD_10Sep25_1522_005
JetMET1_Run2024F-MINIv6NANOv15-v2_NANOAOD_10Sep25_1526_007
JetMET1_Run2024G-MINIv6NANOv15-v2_NANOAOD_10Sep25_1531_009
JetMET1_Run2024H-MINIv6NANOv15-v2_NANOAOD_10Sep25_1535_011
JetMET1_Run2024I-MINIv6NANOv15-v1_NANOAOD_10Sep25_1538_014
JetMET1_Run2024I-MINIv6NANOv15_v2-v2_NANOAOD_10Sep25_1539_015
)

# Loop over datasets
for ds in "${datasets[@]}"; do
  # Strip the trailing _<date>_<time>_<id>
  base_jobname=$(echo "$ds" | sed -E 's/_10Sep25_[0-9]+_[0-9]+$//')

  echo ">>> Submitting job for $ds (jobName=$base_jobname)"
  python3 submit_jobs.py \
    --inputDir /hdfs/store/user/mithakor/2024_skimmed/$ds \
    --destination $DEST \
    --jobName $base_jobname \
    --goldenjson $JSON_BASE/GoldenJSON_2024.json \
    --jetidjson $JSON_BASE/jetid.json \
    --vetomapjson $JSON_BASE/jetvetomaps.json \
    --jercjson $JSON_BASE/jet_jerc.json.gz
done

echo "Submitted all."