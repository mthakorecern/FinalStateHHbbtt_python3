# python3 submit_jobs.py  \
#     --inputDir /hdfs/store/user/mithakor/Tester \
#     --destination /hdfs/store/user/mithakor/2024_Categorized_MTT    \
#     --jobName Test  \
#     --submitDirPath /nfs_scratch/mithakor/Corrections_MTT   \
#     --year 2024 \
#     --isMC  \


# python3 CommonAnalysisWithSystematics_6.py  \
#     --inputFile /hdfs/store/user/mithakor/2024_skimmed_hadded/TbarBQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8.root \
#     --outputFile TbarBQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8.root \
#     --year 2024 \
#     --isMC  \
#     --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python \
#     --runNominal    \
#     &> log_TbarBQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8.txt &  

# python3 CommonAnalysisWithSystematics_6.py  \
#     --inputFile /hdfs/store/user/mithakor/2024_skimmed_hadded/TBbarQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8.root \
#     --outputFile TBbarQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8.root \
#     --year 2024 \
#     --isMC  \
#     --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python \
#     --runNominal    \
#     &> log_TBbarQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8.txt &  

# python3 CommonAnalysisWithSystematics_6.py  \
#     --inputFile /hdfs/store/user/mithakor/2024_skimmed_hadded/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8.root \
#     --outputFile ZZto4L_TuneCP5_13p6TeV_powheg-pythia8.root \
#     --year 2024 \
#     --isMC  \
#     --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python \
#     --runNominal    \
    # &> log_ZZto4L_TuneCP5_13p6TeV_powheg-pythia8.txt &  


    

python3 CommonAnalysisWithSystematics_6.py  \
    --inputFile /hdfs/store/user/mithakor/2024_skimmed_hadded/GluGlutoRadiontoHHto2B2Tau_M-1000.root \
    --outputFile Radion_output.root \
    --year 2024 \
    --isMC  \
    --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python \
    --runNominal    \
    &> Radion.txt &  