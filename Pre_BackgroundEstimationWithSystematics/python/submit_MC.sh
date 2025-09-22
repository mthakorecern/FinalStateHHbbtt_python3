SOURCE=/hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/MC_18Sep25_2312
DEST=/hdfs/store/user/mithakor/2024_Categorized_MTT



# python3 submit_jobs.py  \
#     --inputDir $SOURCE  \
#     --destination $DEST     \
#     --jobName MC_Without_Top  \
#     --submitDirPath /nfs_scratch/mithakor/Corrections_MTT   \
#     --year 2024 \
#     --isMC  \


python3 submit_jobs.py  \
    --inputDir $SOURCE  \
    --destination $DEST     \
    --jobName MC_With_Top  \
    --submitDirPath /nfs_scratch/mithakor/Corrections_MTT   \
    --year 2024 \
    --isMC  \