# python3 submit_jobs.py  \
#     --inputDir /hdfs/store/user/mithakor/Tester \
#     --destination /hdfs/store/user/mithakor/2024_Categorized_MTT    \
#     --jobName Test  \
#     --submitDirPath /nfs_scratch/mithakor/Corrections_MTT   \
#     --year 2024 \
#     --isMC  \


python3 CommonAnalysisWithSystematics_6.py  \
    --inputFile /hdfs/store/user/mithakor/Tester/nanoPostProc_helper-GluGlutoRadiontoHHto2B2Tau_M-1000.root \
    --outputFile GluGlutoRadiontoHHto2B2Tau_M-1000.root \
    --year 2024 \
    --isMC  \
    --cutflowDir /afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python \
    --runNominal  
