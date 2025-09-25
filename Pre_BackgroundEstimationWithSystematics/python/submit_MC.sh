SOURCE=/hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/WithFatJetID_FatJetVetoMaps/MC_Without_Top_24Sep25_1529
DEST=/hdfs/store/user/mithakor/2024_Categorized_MTT/MC



# python3 submit_jobs.py  \
#     --inputDir $SOURCE  \
#     --destination $DEST     \
#     --jobName MC_With_Top  \
#     --submitDirPath /nfs_scratch/mithakor/Corrections_MTT/MC   \
#     --year 2024 \
#     --isMC  \


python3 submit_jobs.py  \
    --inputDir $SOURCE  \
    --destination $DEST     \
    --jobName MC_Without_top  \
    --submitDirPath /nfs_scratch/mithakor/Corrections_MTT/MC   \
    --year 2024 \
    --isMC  \