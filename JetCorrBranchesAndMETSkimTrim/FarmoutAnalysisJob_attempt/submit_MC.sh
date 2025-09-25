DEST=/hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/WithFatJetID_FatJetVetoMaps
JSON_BASE=/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024

python3 submit_jobs.py \
    --inputList Without_Top.txt \
    --destination $DEST \
    --jobName MC_Without_Top \
    --jetidjson $JSON_BASE/jetid.json \
    --vetomapjson $JSON_BASE/jetvetomaps.json \
    --jercjson $JSON_BASE/jet_jerc.json.gz  \
    --isMC

# python3 submit_jobs.py \
#     --inputList Top.txt \
#     --destination $DEST \
#     --jobName Top \
#     --jetidjson $JSON_BASE/jetid.json \
#     --vetomapjson $JSON_BASE/jetvetomaps.json \
#     --jercjson $JSON_BASE/jet_jerc.json.gz  \
#     --isMC