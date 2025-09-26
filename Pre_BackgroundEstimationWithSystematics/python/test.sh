python3 CommonAnalysisWithSystematics_6.py  \
    --inputFile /hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/WithFatJetID_FatJetVetoMaps/MC_Without_Top_24Sep25_1529/nanoPostProc_helper-GluGlutoRadiontoHHto2B2Tau_M-1000.root \
    --outputFile GluGlutoRadiontoHHto2B2Tau_M-1000.root   \
    --year 2024  \
    --isMC  \
    --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python  \
    &> log_test.txt &

python3 CommonAnalysisWithSystematics_6.py  \
    --inputFile /hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/WithFatJetID_FatJetVetoMaps/JETMET_Hadded/JetMET0_Run2024C-MINIv6NANOv15-v1_NANOAOD.root \
    --outputFile JetMET0_Run2024C-MINIv6NANOv15-v1_NANOAOD_Output.root   \
    --year 2024  \
    --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python  \
    &> log_data_test.txt &


