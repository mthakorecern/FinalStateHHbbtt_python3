python3 CommonAnalysisWithSystematics_6.py  \
    --inputFile /hdfs/store/user/mithakor/2024_skimmed_hadded/GluGlutoRadiontoHHto2B2Tau_M-1000.root \
    --outputFile output_Signal.root   \
    --year 2024  \
    --isMC  \
    --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python  \
    --runNominal    \
    &> log_test.txt &

python3 CommonAnalysisWithSystematics_6.py  \
    --inputFile /hdfs/store/user/mithakor/2024_skimmed/JetMET0_Run2024C-MINIv6NANOv15-v1_NANOAOD_10Sep25_1352_000/singleFileSkimForSubmission-00cd90c5-134d-4fb8-9ae9-2f78749e80df.root \
    --outputFile output_data.root   \
    --year 2024  \
    --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python  \
    &> log_data_test.txt &


