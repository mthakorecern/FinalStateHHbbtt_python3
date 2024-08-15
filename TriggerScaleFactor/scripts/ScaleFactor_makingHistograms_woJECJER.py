#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
#Import the JEC and JER NanoAODTools correction Modules
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
import glob
import numpy as np
import json

import argparse
#Import and load the MET-XY correction header file for UL
#MET phi correction: https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection_withUL17andUL18andUL16.h, https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections#xy_Shift_Correction_MET_phi_modu
print 'gRoot ProcessLine XYMETCorrection = ', ROOT.gROOT.ProcessLine(".L /afs/hep.wisc.edu/home/parida/public/HHbbtt_Analysis_Scripts/TriggerScaleFactor/scripts/ULMETXY_Correction.h")

class dataEfficiencty(Module):
    def __init__(self):
        self.writeHistFile=True
        #self.year = year

    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)
        #self.binEdges = np.array([80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 700.0, 800.0, 900.0, 1000.0])
        self.binEdges = np.array([80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 350.0, 400.0, 450.0, 500.0])
        self.nbin = len(self.binEdges) - 1
        #self.h_data_denominator =  ROOT.TH1F("h_data_denominator" , "; E_{T}^{miss} [GeV]" , 46, 80., 1000.)
        #self.h_data_numerator =  ROOT.TH1F("h_data_numerator" , "; E_{T}^{miss} [GeV]" , 46, 80., 1000.)
        self.h_data_denominator =  ROOT.TH1F("h_data_denominator" , "; E_{T}^{miss} [GeV]" , self.nbin, self.binEdges)
        self.h_data_numerator =  ROOT.TH1F("h_data_numerator" , "; E_{T}^{miss} [GeV]" , self.nbin, self.binEdges)
        self.delR_Muon_fatJet =  ROOT.TH1F("delR" , "; delR" , 50, 0., 5. )
        self.addObject(self.h_data_denominator)
        self.addObject(self.h_data_numerator)
        self.addObject(self.delR_Muon_fatJet)

    def analyze(self, event):
        #met       = Object(event, "MET")
        #hlt       = Object(event, "HLT")
        Muons = Collection(event, 'Muon', 'nMuon')
        FatJet = Collection(event,'FatJet','nFatJet')
        #goodMuons = filter(lambda j : ((j.pt > 30) and (j.pfRelIso04_all<0.25) and (j.tightId) and (abs(j.eta)<2.4)),Muons)
        goodMuons = filter(lambda j : ((j.pt > 30) and (j.pfRelIso04_all<0.15) and (j.tightId) and (abs(j.eta)<2.4)),Muons)
        goodFatJet = filter(lambda x : ((x.pt > 170) and (abs(x.eta) < 2.5) and (x.jetId>1)),FatJet)
        #goodFatJet = filter(lambda x : (abs(x.eta) < 2.5),FatJet)

        #METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, TString year, bool isMC, int npv, bool isUL =false,bool ispuppi=false)
        #if(isMC && year == "2016" && !isUL) runera = y2016MC;
        #else if(isMC && year == "2017" && !isUL) {runera = y2017MC; usemetv2 =true;}
        #else if(isMC && year == "2018" && !isUL) runera = y2018MC;
        #else if(isMC && year == "2016APV" && isUL) runera = yUL2016MCAPV;
        #else if(isMC && year == "2016nonAPV" && isUL) runera = yUL2016MCnonAPV;
        #else if(isMC && year == "2017" && isUL) runera = yUL2017MC;
        #else if(isMC && year == "2018" && isUL) runera = yUL2018MC;
        if (args.year=="2016"):
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi, event.run,"2016nonAPV",False,event.PV_npvs,True,False)
        elif (args.year=="2016APV"):
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi, event.run,"2016APV",False,event.PV_npvs,True,False)
        elif (args.year=="2017"):
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi, event.run,"2017",False,event.PV_npvs,True,False)
        elif (args.year=="2018"):
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi, event.run,"2018",False,event.PV_npvs,True,False)


        METcorrected_pt = METcorrected_pt_phi[0]
        METcorrected_phi = METcorrected_pt_phi[1]

        #Throw away events that have corrected MET < 80 GeV and have zero good Muons
        if ((len(goodMuons)!=1) or (METcorrected_pt<80)):
            return False

        if (args.year=="2016" or args.year=="2016APV"):
            refTrigger = ((event.HLT_IsoMu24) or (event.HLT_IsoTkMu24))
        else:
            refTrigger = ((event.HLT_IsoMu24))

        #self.h_datapassreftrig.Fill(refTrigger)



        if not refTrigger:
            return False

        #for jet in goodFatJet:
        if (len(goodFatJet)==0):
            return False

        if (goodMuons[0].p4().DeltaR(goodFatJet[0].p4()) <= 0.8):
            return False

        self.delR_Muon_fatJet.Fill(goodMuons[0].p4().DeltaR(FatJet[0].p4()))
        #if (goodMuons[0].p4().DeltaR(FatJet[0].p4()) > 0.8):
        self.h_data_denominator.Fill(METcorrected_pt)
        
        
        

        #if(refTrigger):
        #try:
        if (args.year=="2016"):
            metTrigger = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned or event.HLT_PFMET170_HBHE_BeamHaloCleaned)
            #if(metTrigger and ((goodMuons[0].p4().DeltaR(FatJet[0].p4())) > 0.8)):
            if(metTrigger):
                self.h_data_numerator.Fill(METcorrected_pt)
        #except:
        elif (args.year=="2016APV"):
            metTriggerSansBeamHaloCleaned = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned)
            #if(metTriggerSansBeamHaloCleaned and ((goodMuons[0].p4().DeltaR(FatJet[0].p4())) > 0.8)):
            if(metTriggerSansBeamHaloCleaned):
                self.h_data_numerator.Fill(METcorrected_pt)

        elif (args.year=="2017"):
            metTrigger2017 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
            if(metTrigger2017):
                self.h_data_numerator.Fill(METcorrected_pt)

        elif (args.year=="2018"):
            metTrigger2018 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
            if(metTrigger2018):
                self.h_data_numerator.Fill(METcorrected_pt)
            



        return True          

class mcEfficiencty(Module):
    def __init__(self,year):
        self.year = year
        if self.year == '2016': 
            self.LHCLumi =  16.81e15
        elif self.year == '2016APV':
            self.LHCLumi = 19.52e15
        elif self.year == '2017':
            self.LHCLumi = 41.48e15
        elif self.year == '2018':
            self.LHCLumi = 59.83e15       
        self.writeHistFile=True
    
    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)
        #float self.binEdges[22] = [80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 700.0, 800.0, 900.0, 1000.0]
        #self.binEdges = np.array([80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 700.0, 800.0, 900.0, 1000.0])
        self.binEdges = np.array([80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 350.0, 400.0, 450.0, 500.0])
        self.nbin = len(self.binEdges) - 1
        #self.h_mc_denominator =  ROOT.TH1F("h_mc_denominator" , "; E_{T}^{miss} [GeV]" , 46, 80., 1000.)
        #self.h_mc_numerator =  ROOT.TH1F("h_mc_numerator" , "; E_{T}^{miss} [GeV]" , 46, 80., 1000.)
        self.h_mc_denominator =  ROOT.TH1F("h_mc_denominator" , "; E_{T}^{miss} [GeV]" , self.nbin, self.binEdges)
        self.h_mc_numerator =  ROOT.TH1F("h_mc_numerator" , "; E_{T}^{miss} [GeV]" , self.nbin, self.binEdges)
        self.delR_Muon_fatJet =  ROOT.TH1F("delR" , "; delR" , 50, 0., 5. )
        self.addObject(self.h_mc_denominator )
        self.addObject(self.h_mc_numerator )

        self.xsJSON = os.environ['CMSSW_BASE']+"/src/MetaData/xsJSON/%s_Samples.json"%(self.year) #Add a customization to different yrs
        with open(self.xsJSON,'r') as xsjson:
            self.xsInfo = json.load(xsjson)


        self.countsJSON = os.environ['CMSSW_BASE']+'/src/MetaData/SumOfGenWeights/%s_weight.json'%(self.year) #Add a customization to different yrs
        with open(self.countsJSON,'r') as countsjson:
            self.countsInfo = json.load(countsjson)



    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        print(" The inputFile name gotten from the beginFile fn : ",type(inputFile),"   ",inputFile.GetName())
        filenameForJson = (inputFile.GetName().strip().split('/')[-1]).split('.')[-2]
        self.filename = filenameForJson
        self.XS = self.xsInfo[self.filename]['XS'] * 1e-12
        self.totalNumberOfEvents = self.countsInfo[self.filename]
        print(" Filename for JSON : ",self.filename," XS = ",self.XS," sumofgenwts = ",self.totalNumberOfEvents)



    def analyze(self, event):

        eventXsWt = (self.XS * self.LHCLumi * event.genWeight) / self.totalNumberOfEvents
        #met       = Object(event, "MET")
        #hlt       = Object(event, "HLT")
        Muons = Collection(event, 'Muon', 'nMuon')
        FatJet = Collection(event,'FatJet','nFatJet')
        #goodMuons = filter(lambda j : ((j.pt > 30) and (j.pfRelIso04_all<0.25) and (j.tightId) and (abs(j.eta)<2.4)),Muons)
        goodMuons = filter(lambda j : ((j.pt > 30) and (j.pfRelIso04_all<0.15) and (j.tightId) and (abs(j.eta)<2.4)),Muons)
        goodFatJet = filter(lambda x : ((x.pt > 170) and (abs(x.eta) < 2.5) and (x.jetId>1)),FatJet)
        #goodFatJet = filter(lambda x : (abs(x.eta) < 2.5),FatJet)


        #METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, TString year, bool isMC, int npv, bool isUL =false,bool ispuppi=false)
        #if(isMC && year == "2016" && !isUL) runera = y2016MC;
        #else if(isMC && year == "2017" && !isUL) {runera = y2017MC; usemetv2 =true;}
        #else if(isMC && year == "2018" && !isUL) runera = y2018MC;
        #else if(isMC && year == "2016APV" && isUL) runera = yUL2016MCAPV;
        #else if(isMC && year == "2016nonAPV" && isUL) runera = yUL2016MCnonAPV;
        #else if(isMC && year == "2017" && isUL) runera = yUL2017MC;
        #else if(isMC && year == "2018" && isUL) runera = yUL2018MC;
        if (args.year=="2016"):
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi,event.run,"2016nonAPV",True,event.PV_npvs,True,False)
        elif (args.year=="2016APV"):
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi,event.run,"2016APV",True,event.PV_npvs,True,False)
        elif (args.year=="2017"):
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi,event.run,"2017",True,event.PV_npvs,True,False)
        elif (args.year=="2018"):
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi,event.run,"2018",True,event.PV_npvs,True,False)
        
        METcorrected_pt = METcorrected_pt_phi[0]
        METcorrected_phi = METcorrected_pt_phi[1]


        #Throw away events that have corrected MET < 80 GeV and have zero good Muons
        if ((len(goodMuons)!=1) or (METcorrected_pt<80)):
            return False

        if (args.year=="2016" or args.year=="2016APV"):
            refTrigger = ((event.HLT_IsoMu24) or (event.HLT_IsoTkMu24))
        else:
            refTrigger = ((event.HLT_IsoMu24))


        #self.h_mcpassreftrig.Fill(refTrigger)

        if not refTrigger:
            return False

        #for jet in goodFatJet:
        if (len(goodFatJet)==0):
            return False

        if (goodMuons[0].p4().DeltaR(goodFatJet[0].p4()) <= 0.8):
            return False
        
        #if (event.FinalWeighting<=0):
        #    return False

        self.delR_Muon_fatJet.Fill(goodMuons[0].p4().DeltaR(FatJet[0].p4()))
        
        #self.h_mc_denominator.Fill(METcorrected_pt,event.FinalWeighting)
        self.h_mc_denominator.Fill(METcorrected_pt,event.pileupWeighting*eventXsWt)
        #self.h_mc_denominator.Fill(METcorrected_pt)
        #try:
        if (args.year=="2016"):
            metTrigger = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned or event.HLT_PFMET170_HBHE_BeamHaloCleaned)
            if(metTrigger):
                #self.h_mc_numerator.Fill(METcorrected_pt,event.FinalWeighting)
                self.h_mc_numerator.Fill(METcorrected_pt,event.pileupWeighting*eventXsWt)
        #except:
        elif (args.year=="2016APV"): 
            metTriggerSansBeamHaloCleaned = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned)
            if(metTriggerSansBeamHaloCleaned):
                #self.h_mc_numerator.Fill(METcorrected_pt,event.FinalWeighting)
                self.h_mc_numerator.Fill(METcorrected_pt,event.pileupWeighting*eventXsWt)

        elif (args.year=="2017"):
            metTrigger2017 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
            if(metTrigger2017):
                #self.h_mc_numerator.Fill(METcorrected_pt,event.FinalWeighting)
                self.h_mc_numerator.Fill(METcorrected_pt,event.pileupWeighting*eventXsWt)                            
        
        elif (args.year=="2018"):
            metTrigger2018 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
            if(metTrigger2018):
                #self.h_mc_numerator.Fill(METcorrected_pt,event.FinalWeighting)
                self.h_mc_numerator.Fill(METcorrected_pt,event.pileupWeighting*eventXsWt)      

        return True


    

parser = argparse.ArgumentParser(description='Script create efficiency histograms for data and mc to be used to compute scalefactors')
#parser.add_argument('--dataPath',help="enter the path to data files",default="")
parser.add_argument('--mcPath',help="enter the path to mc files",default="")
parser.add_argument('--year','-y',help="Enter the year of data taking",choices=["2016","2016APV","2017","2018"], required=True)
args = parser.parse_args()

#fnames_data = glob.glob(args.dataPath + "/*.root")

data = lambda: dataEfficiencty()
mc = lambda: mcEfficiencty(args.year)

#Consider all the MC samples except the signal samples
fnames_mc = list(set(glob.glob(args.mcPath + "/*.root"))-set(glob.glob(args.mcPath + "/Radion*.root")))
#fnames_mc = glob.glob(args.mcPath + "/WJetsToLNu_HT-400To600*.root")
print (fnames_mc)

eventSelection2016 =["PV_ndof > 4",
					"abs(PV_z) < 24",
					"sqrt(PV_x*PV_x+PV_y*PV_y) < 2",
                    "Flag_goodVertices",
                    "Flag_globalSuperTightHalo2016Filter",
                    "Flag_HBHENoiseFilter",
                    "Flag_HBHENoiseIsoFilter",
                    "Flag_EcalDeadCellTriggerPrimitiveFilter",
                    "Flag_BadPFMuonFilter",
                    "Flag_BadPFMuonDzFilter",
                    "Flag_eeBadScFilter",
                    "Flag_hfNoisyHitsFilter"]

eventSelection2017_18 =["PV_ndof > 4",
					"abs(PV_z) < 24",
					"sqrt(PV_x*PV_x+PV_y*PV_y) < 2",
                    "Flag_goodVertices",
                    "Flag_globalSuperTightHalo2016Filter",
                    "Flag_HBHENoiseFilter",
                    "Flag_HBHENoiseIsoFilter",
                    "Flag_EcalDeadCellTriggerPrimitiveFilter",
                    "Flag_BadPFMuonFilter",
                    "Flag_BadPFMuonDzFilter",
                    "Flag_hfNoisyHitsFilter",
                    "Flag_eeBadScFilter",
                    "Flag_ecalBadCalibFilter"]
if(args.year=="2016" or args.year=="2016APV"):
    preselection = "&&".join(eventSelection2016)
elif(args.year=="2017" or args.year=="2018"):
    preselection = "&&".join(eventSelection2017_18)


#Create the JME modules for year 2016 MC and 2016 Data
#jetmetCorrectorUL2016MC = createJMECorrector(isMC=True, dataYear='UL2016', jesUncert="Total",jetType="AK4PFchs")
#jetmetCorrectorUL2016RunF = createJMECorrector(isMC=False, dataYear='UL2016', runPeriod="F",jesUncert="Total",jetType="AK4PFchs")
#jetmetCorrectorUL2016RunG = createJMECorrector(isMC=False, dataYear='UL2016', runPeriod="G",jesUncert="Total",jetType="AK4PFchs")
#jetmetCorrectorUL2016RunH = createJMECorrector(isMC=False, dataYear='UL2016', runPeriod="H",jesUncert="Total",jetType="AK4PFchs")

if (args.year=="2016"):
    #Create and run post processor for 2016 MC
    print("Creating/Running post processor for 2016 MC")
    ##post2016_MC = PostProcessor("2016/.",fnames_mc,cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016MC(),mc()],histFileName="2016MC_MCHistos.root",histDirName="MCHistos")
    #post2016_MC = PostProcessor("2016/.",fnames_mc,cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016MC(),mc()],noOut=False,histFileName="2016/2016MC_MCHistos.root",histDirName="MCHistos",outputbranchsel="drop.txt")
    post2016_MC = PostProcessor("2016/.",fnames_mc,cut=preselection,branchsel=None,modules=[mc()],noOut=True,histFileName="2016/2016MC_MCHistos.root",histDirName="MCHistos")    
    #post2016_MC = PostProcessor(".",["/hdfs/store/user/gparida/HHbbtt/Hadded_Skimmed_Files/Full_Production_August19_23_EleBitMapFix/2016/RadionTohhTohtatahbb_narrow_M-4500.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016MC(),mc()],histFileName="2016MC_MCHistos.root",histDirName="MCHistos")
    post2016_MC.run()


    #Creating and running Post processor of 2016 Data
    print("Creating/Running post processor for 2016 Data - SingleMu_2016F")
    #post2016_RunF=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016F.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016RunF(),data()],noOut=True,histFileName="2016RunF_DataHistos.root",histDirName="DataHistos")
    post2016_RunF=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016F.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016/2016RunF_DataHistos.root",histDirName="DataHistos")

    post2016_RunF.run()

    print("Creating/Running post processor for 2016 Data - SingleMu_2016G")
    #post2016_RunG=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016G.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016RunF(),data()],noOut=True,histFileName="2016RunF_DataHistos.root",histDirName="DataHistos")
    post2016_RunG=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016G.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016/2016RunG_DataHistos.root",histDirName="DataHistos")
    post2016_RunG.run()

    print("Creating/Running post processor for 2016 Data - SingleMu_2016H")
    #post2016_RunH=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016H.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016RunH(),data()],noOut=True,histFileName="2016RunH_DataHistos.root",histDirName="DataHistos")
    post2016_RunH=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016H.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016/2016RunH_DataHistos.root",histDirName="DataHistos")
    post2016_RunH.run()

elif (args.year=="2016APV"):
    #Create and run post processor for 2016APV MC
    print("Creating/Running post processor for 2016APV MC")
    APV2016_MC = PostProcessor("2016APV/.",fnames_mc,cut=preselection,branchsel=None,modules=[mc()],noOut=True,histFileName="2016APV/2016APVMC_MCHistos.root",histDirName="MCHistos")
    APV2016_MC.run()

    #Creating and running Post processor of 2016APV Run2016B-ver1_HIPM  Data
    print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016B-ver1_HIPM")
    Run2016B_ver1_HIPM=PostProcessor("2016APV/.",[args.mcPath + "/SingleMu/SingleMu_Run2016B-ver1_HIPM.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016APV/Run2016B_ver1_HIPM_DataHistos.root",histDirName="DataHistos")
    Run2016B_ver1_HIPM.run()

    #Creating and running Post processor of 2016APV Run2016B-ver2_HIPM  Data
    print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016B-ver2_HIPM")
    Run2016B_ver2_HIPM=PostProcessor("2016APV/.",[args.mcPath + "/SingleMu/SingleMu_Run2016B-ver2_HIPM.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016APV/Run2016B_ver2_HIPM_DataHistos.root",histDirName="DataHistos")
    Run2016B_ver2_HIPM.run()

    #Creating and running Post processor of 2016APV Run2016C-HIPM  Data
    print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016C-HIPM")
    Run2016C_HIPM=PostProcessor("2016APV/.",[args.mcPath + "/SingleMu/SingleMu_Run2016C-HIPM.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016APV/Run2016C_HIPM_DataHistos.root",histDirName="DataHistos")
    Run2016C_HIPM.run()

    #Creating and running Post processor of 2016APV Run2016D-HIPM  Data
    print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016D-HIPM")
    Run2016D_HIPM=PostProcessor("2016APV/.",[args.mcPath + "/SingleMu/SingleMu_Run2016D-HIPM.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016APV/Run2016D_HIPM_DataHistos.root",histDirName="DataHistos")
    Run2016D_HIPM.run()

    #Creating and running Post processor of 2016APV Run2016E-HIPM  Data
    print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016E-HIPM")
    Run2016E_HIPM=PostProcessor("2016APV/.",[args.mcPath + "/SingleMu/SingleMu_Run2016E-HIPM.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016APV/Run2016E_HIPM_DataHistos.root",histDirName="DataHistos")
    Run2016E_HIPM.run()

    #Creating and running Post processor of 2016APV Run2016F-HIPM  Data
    print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016F-HIPM")
    Run2016F_HIPM=PostProcessor("2016APV/.",[args.mcPath + "/SingleMu/SingleMu_Run2016F-HIPM.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016APV/Run2016F_HIPM_DataHistos.root",histDirName="DataHistos")
    Run2016F_HIPM.run()

elif (args.year=="2018"):
    print("Creating/Running post processor for 2018 MC")
    post2018_MC = PostProcessor("2018/.",fnames_mc,cut=preselection,branchsel=None,modules=[mc()],noOut=True,histFileName="2018/2018MC_MCHistos.root",histDirName="MCHistos")
    post2018_MC.run()  

    print("Creating/Running post processor for 2018 Data - SingleMu_Run2018A")
    Run2018A = PostProcessor("2018/.",[args.mcPath + "/SingleMu/SingleMu_Run2018A.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2018_Leagcy.json",noOut=True,histFileName="2018/Run2018A_DataHistos.root",histDirName="DataHistos")
    Run2018A.run()        

    print("Creating/Running post processor for 2018 Data - SingleMu_Run2018B")
    Run2018B = PostProcessor("2018/.",[args.mcPath + "/SingleMu/SingleMu_Run2018B.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2018_Leagcy.json",noOut=True,histFileName="2018/Run2018B_DataHistos.root",histDirName="DataHistos")
    Run2018B.run() 

    print("Creating/Running post processor for 2018 Data - SingleMu_Run2018C")
    Run2018C = PostProcessor("2018/.",[args.mcPath + "/SingleMu/SingleMu_Run2018C.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2018_Leagcy.json",noOut=True,histFileName="2018/Run2018C_DataHistos.root",histDirName="DataHistos")
    Run2018C.run()     

    print("Creating/Running post processor for 2018 Data - SingleMu_Run2018D")
    #data2018D_files = glob.glob(args.mcPath + "/SingleMu/SingleMu_Run2018D*.root")
    Run2018D = PostProcessor("2018/.",glob.glob(args.mcPath + "/SingleMu/SingleMu_Run2018D*.root"),cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2018_Leagcy.json",noOut=True,histFileName="2018/Run2018D_DataHistos.root",histDirName="DataHistos")
    Run2018D.run()     


elif (args.year=="2017"):
    print("Creating/Running post processor for 2017 MC")
    post2017_MC = PostProcessor("2017/.",fnames_mc,cut=preselection,branchsel=None,modules=[mc()],noOut=True,histFileName="2017/2017MC_MCHistos.root",histDirName="MCHistos")
    post2017_MC.run()  

    print("Creating/Running post processor for 2017 Data - SingleMu_Run2017B")
    Run2017B = PostProcessor("2017/.",[args.mcPath + "/SingleMu/SingleMu_Run2017B.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2017_Leagcy.json",noOut=True,histFileName="2017/Run2017B_DataHistos.root",histDirName="DataHistos")
    Run2017B.run() 

    print("Creating/Running post processor for 2017 Data - SingleMu_Run2017C")
    Run2017C = PostProcessor("2017/.",[args.mcPath + "/SingleMu/SingleMu_Run2017C.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2017_Leagcy.json",noOut=True,histFileName="2017/Run2017C_DataHistos.root",histDirName="DataHistos")
    Run2017C.run() 

    print("Creating/Running post processor for 2017 Data - SingleMu_Run2017D")
    Run2017D = PostProcessor("2017/.",[args.mcPath + "/SingleMu/SingleMu_Run2017D.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2017_Leagcy.json",noOut=True,histFileName="2017/Run2017D_DataHistos.root",histDirName="DataHistos")
    Run2017D.run() 

    print("Creating/Running post processor for 2017 Data - SingleMu_Run2017E")
    Run2017E = PostProcessor("2017/.",[args.mcPath + "/SingleMu/SingleMu_Run2017E.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2017_Leagcy.json",noOut=True,histFileName="2017/Run2017E_DataHistos.root",histDirName="DataHistos")
    Run2017E.run() 

    print("Creating/Running post processor for 2017 Data - SingleMu_Run2017F")
    Run2017F = PostProcessor("2017/.",[args.mcPath + "/SingleMu/SingleMu_Run2017F.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2017_Leagcy.json",noOut=True,histFileName="2017/Run2017F_DataHistos.root",histDirName="DataHistos")
    Run2017F.run() 