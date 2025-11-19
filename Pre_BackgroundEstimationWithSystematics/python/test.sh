python3 CommonAnalysisWithSystematics_6.py  \
    --inputFile /hdfs/store/user/mithakor/2024_skimmed_hadded/WWto4Q_TuneCP5_13p6TeV_powheg-pythia8_0.root \
    --outputFile test_output_WWto4Q_TuneCP5_13p6TeV_powheg-pythia8_0.root   \
    --year 2024  \
    --isMC  \
    --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python  \
    --runNominal    \
    &> log_test_WWto4Q_TuneCP5_13p6TeV_powheg-pythia8_0.txt &

# python3 CommonAnalysisWithSystematics_6.py  \
#     --inputFile /hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/WithFatJetID_FatJetVetoMaps/JETMET_Hadded/JetMET0_Run2024C-MINIv6NANOv15-v1_NANOAOD.root \
#     --outputFile JetMET0_Run2024C-MINIv6NANOv15-v1_NANOAOD_Output.root   \
#     --year 2024  \
#     --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python  \
#     &> log_data_test.txt &


