# SOURCE_Top=/hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/JME_JetFatJetMET_JES_JER_Updated/Top_05Oct25_1728

# SOURCE_Without_Top=/hdfs/store/user/mithakor/2024_JETID_JetVetoMaps_JES_JER_GlobalparT3mass_softdrop/JME_JetFatJetMET_JES_JER_Updated/MC_Without_Top_05Oct25_1550

SOURCE=/hdfs/store/user/mithakor/2024_skimmed_hadded


DEST=/hdfs/store/user/mithakor/2024_Categorized_MTT/MC

# python3 submit_jobs.py  \
#     --inputDir $SOURCE_Top  \
#     --destination $DEST     \
#     --jobName MC_With_Top_MET_120_Nominal  \
#     --submitDirPath /nfs_scratch/mithakor/Corrections_MTT/MC   \
#     --year 2024 \
#     --isMC  \
#     --runNominal    \
#     &> log_MC_With_Top.txt & 


# python3 submit_jobs.py  \
#     --inputDir $SOURCE_Without_Top  \
#     --destination $DEST     \
#     --jobName MC_Without_top_MET_120_Nominal  \
#     --submitDirPath /nfs_scratch/mithakor/Corrections_MTT/MC   \
#     --year 2024 \
#     --isMC  \
#     --runNominal    \
#     &> log_MC_Without_Top.txt &

python3 submit_jobs.py  \
    --inputDir $SOURCE  \
    --destination $DEST     \
    --jobName MC_MET_120_Nocorrections_veto_final  \
    --submitDirPath /nfs_scratch/mithakor/Corrections_MTT/MC   \
    --year 2024 \
    --isMC  \
    --runNominal    \
    &> log_MC_MET_120_Nocorrections.txt &