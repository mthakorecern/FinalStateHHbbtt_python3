#!/bin/bash
set -euo pipefail

# Common JSON inputs
GOLDENJSON="/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/GoldenJSON_2024.json"
JETIDJSON="/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/jetid.json"
VETOMAPJSON="/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/jetvetomaps.json"
JERCJSON="/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/jet_jerc.json.gz"

# Directories
INPUT_BASE="/hdfs/store/user/mithakor/2024_skimmed"
OUTPUT_BASE="/hdfs/store/user/mithakor/2024_JETcorrected"
TEMP_BASE="/nfs_scratch/mithakor/temp"

# List of input dirs
INPUT_DIRS=(
  # JetMET0_Run2024C-MINIv6NANOv15-v1_NANOAOD_10Sep25_1352_000
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

# Loop
for d in "${INPUT_DIRS[@]}"; do
  INPUT_DIR="$INPUT_BASE/$d"

  OUTNAME=$(echo "$d" | sed -E 's/(_NANOAOD.*)//')

  echo ">>> Processing $d → $OUTNAME.root"

  python3 addJetIDbranch_JetVetoMaps_JECapplied_branches.py \
    --goldenjson "$GOLDENJSON" \
    --jetidjson "$JETIDJSON" \
    --vetomapjson "$VETOMAPJSON" \
    --jercjson "$JERCJSON" \
    --inputDir "$INPUT_DIR" \
    --outputDir "$OUTPUT_BASE" \
    --tempDir "$TEMP_BASE/$OUTNAME" \
    --finalHadd \
    --finalOut "$OUTNAME.root" \
    &> "log_${OUTNAME}.txt"

done
