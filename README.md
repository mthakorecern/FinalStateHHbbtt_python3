# FinalStateHHbbtt

---

## Purpose
This setup is for aimed at performing the tasks of:
- event and object selection: Filters events based on muons, electrons, jets, and MET thresholds defined in config files. 
- channel categorization: Assigns events to different τ decay channels (e.g. eμ, μτ, ττ).
- Addition of additional branches: Computes and appends new high-level physics variables (mass, isolation, b-tag flags), then drops unused branches to reduce file size.

It was originally written by Ganesh Parida for Run-2 UL analysis. It will be modified for Run-3 analysis (with Python 3).

### Setup the environment/required dependencies:

```
cmsrel CMSSW_15_0_2
cd CMSSW_15_0_2/src/
cmsenv
```

### Obtain Grid proxy
```
voms-proxy-init -voms cms -out /tmp/x509up_u10196 -valid 192:00
```
### Setting up NanoAOD Tools

By default NanoAOD Tools is already integrated into CMSSW. However, there are some additional modules helpful for the analysis, which were written by Common Analysis Tools Group and present externally in a repository. We will set it up properly. In the src area, we do

```
git clone git@github.com:mthakorecern/nanoAOD-tools-modules.git  PhysicsTools/NATModules
cd PhysicsTools/NATModules
git switch softdropmass_corrections
cd $CMSSW_BASE/src 
scram b -j 8
``` 
### Setting up Fast MTT 

This tool is essential for reconstruction of the Higgs four vector using Taus and its decay products. In the src area, we set it up as 

```
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit -b fastMTT_19_02_2019
git clone https://github.com/SVfit/SVfitTF TauAnalysis/SVfitTF
cd $CMSSW_BASE/src
scram b -j 8
```

### Setting up Tau POG’s SF module
Although this module is not required right now we set it up as there might be some calculations that might require this module

```
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs
cmsenv
scram b -j8
```

### Main code
Now that these three modules are installed, then one can directly go onto setting up the main code within the src area

```
cmsenv
git clone https://github.com/mthakorecern/FinalStateHHbbtt_python3.git
cd FinalStateHHbbtt
git switch NanoAODdefaultbranches
cd $CMSSW_BASE/src
scram b -j 8
```

This last scram should set up everything. The main script for channel categorization is in FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python - CommonAnalysisWithSystematics_6.py

## Running the code

For testing, To run this code for simply one file, do

``` 
python3 CommonAnalysisWithSystematics_6.py  \
    --inputFile /hdfs/store/user/mithakor/Branch_addition_categorization/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_1.root  \
    --outputFile Test_Output.root    \
    --year 2024 \
    --isMC \
    --cutflowDir $pwd \
    --runNominal  
```

This folder also contains some root files that can be tested ```/hdfs/store/user/mithakor/Branch_addition_categorization```
