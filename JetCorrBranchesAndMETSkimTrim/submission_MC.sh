#!/bin/bash
set -euo pipefail

# Common JSON inputs
GOLDENJSON="/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/GoldenJSON_2024.json"
JETIDJSON="/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/jetid.json"
VETOMAPJSON="/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/jetvetomaps.json"
JERCJSON="/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/jet_jerc.json.gz"

# Directories
INPUT_BASE="/hdfs/store/user/mithakor/2024_skimmed_hadded"
OUTPUT_BASE="/hdfs/store/user/mithakor/2024_JETcorrected"
TEMP_BASE="/nfs_scratch/mithakor/temp"


python3 addJetIDbranch_JetVetoMaps_JECapplied_branches.py \
  --goldenjson "$GOLDENJSON"  \
  --jetidjson "$JETIDJSON" \
  --vetomapjson "$VETOMAPJSON" \
  --jercjson "$JERCJSON" \
  --inputDir "$INPUT_BASE" \
  --outputDir "$OUTPUT_BASE" \
  --tempDir "$TEMP_BASE" \
  --isMC \
  &> log_MC.txt &
