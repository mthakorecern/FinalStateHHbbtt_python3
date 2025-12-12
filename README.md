# FinalStateHHbbtt

---

## Purpose
This setup is for aimed at performing the tasks of:
- event and object selection: Filters events based on muons, electrons, jets, and MET thresholds defined in config files. 
- channel categorization: Assigns events to different τ decay channels (e.g. eμ, μτ, ττ).
- Addition of additional branches: Computes and appends new high-level physics variables (mass, isolation, b-tag flags), then drops unused branches to reduce file size.


It was originally written by Ganesh Parida for Run-2 UL analysis. It will be modified for Run-3 analysis (with Python 3).

## Setup and Installation

```
cmsrel CMSSW_15_0_10
cd CMSSW_15_0_10/src/
cmsenv
git clone https://github.com/mthakorecern/FinalStateHHbbtt_python3.git
cd $CMSSW_BASE/src/
scram b -j4
```
## Exeuction of code 

Inside the folder ```Pre_BackgroundEstimationWithSystematics```, one can try the command 
```
python3 CommonAnalysisWithSystematics_6.py  \
    --inputFile /hdfs/store/user/mithakor/Branch_addition_categorization/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_1.root  \
    --outputFile Test_Outputroot    \
    --year 2024 \
    --isMC \
    --cutflowDir $pwd \
    --runNominal    \
```

The folder ```/hdfs/store/user/mithakor/Branch_addition_categorization``` has some root files corresponding to Data, Signal and Background, which can be used for testing this code. 


