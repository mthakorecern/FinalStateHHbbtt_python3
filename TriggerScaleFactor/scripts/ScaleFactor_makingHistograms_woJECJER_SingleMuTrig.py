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

import argparse
#Import and load the MET-XY correction header file for UL
print 'gRoot ProcessLine XYMETCorrection = ', ROOT.gROOT.ProcessLine(".L /afs/hep.wisc.edu/home/parida/public/HHbbtt_Analysis_Scripts/TriggerScaleFactor/scripts/ULMETXY_Correction.h")

class dataEfficiencty(Module):
    def __init__(self):
        self.writeHistFile=True
    
    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)
        self.h_datapassreftrig  = ROOT.TH1F("h_datapassreftrig" , "; passed ref trigger" , 2 , 0. , 2. )
        self.h_datapassmettrig  = ROOT.TH1F("h_datapassmettrig" ,  "; passed met trigger",2 ,0. ,2. )
        self.h_data_denominator =  ROOT.TH1F("h_data_denominator" , "; E_{T}^{miss} [GeV]" , 21, 80., 500. )
        self.h_data_numerator =  ROOT.TH1F("h_data_numerator" , "; E_{T}^{miss} [GeV]" , 21, 80., 500. )
        
        self.h_data_denominator_MuPt =  ROOT.TH1F("h_data_denominator_MuPt" , "; E_{T}^{miss} [GeV]" , 20, 20., 420. )
        self.h_data_numerator_MuPt =  ROOT.TH1F("h_data_numerator_MuPt" , "; E_{T}^{miss} [GeV]" , 20, 20., 420. )

        self.addObject(self.h_datapassreftrig )
        self.addObject(self.h_data_denominator )
        self.addObject(self.h_data_numerator )
        self.addObject(self.h_datapassmettrig)
        self.addObject(self.h_data_denominator_MuPt)
        self.addObject(self.h_data_numerator_MuPt)

    def analyze(self, event):
        #met       = Object(event, "MET")
        #hlt       = Object(event, "HLT")
        Muons = Collection(event, 'Muon', 'nMuon')
        goodMuons = filter(lambda j : (j.pt > 30) and (j.pfRelIso04_all<0.25) and (j.tightId),Muons)

        METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi, event.run,"2016nonAPV",False,event.PV_npvs,True,False)
        METcorrected_pt = METcorrected_pt_phi[0]
        METcorrected_phi = METcorrected_pt_phi[1]

        #Throw away events that have corrected MET < 80 GeV and have zero good Muons
        if ((len(goodMuons)==0) or (METcorrected_pt<80)):
            return False

        refTrigger = ((event.HLT_IsoMu24) or (event.HLT_IsoTkMu24))


        self.h_datapassreftrig.Fill(refTrigger)

        #if not refTrigger:
        #    return False
        #self.h_mc_denominator.Fill(METcorrected_pt,event.FinalWeighting)
        self.h_data_denominator.Fill(METcorrected_pt)
        self.h_data_denominator_MuPt.Fill(goodMuons[0].pt)


        #if(event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight 
        #  or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
        #  or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight
        #  or event.HLT_PFMET110_PFMHT110_IDTight
        #  or event.HLT_PFMET120_PFMHT120_IDTight
        #  or event.HLT_PFMET170_HBHECleaned
        #  or event.HLT_PFMET170_HBHE_BeamHaloCleaned):
        if (refTrigger):
            self.h_data_numerator.Fill(METcorrected_pt)
            self.h_data_numerator_MuPt.Fill(goodMuons[0].pt)
            self.h_datapassmettrig.Fill(1)
        
        return True

class mcEfficiencty(Module):
    def __init__(self):
        #self.year = year
        #self.isMC isMC
        #print ("FileName =",filename)
        self.writeHistFile=True
    
    def beginJob(self,histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)
        self.h_mcpassreftrig  = ROOT.TH1F("h_mcpassreftrig" , "; passed ref trigger" , 2 , 0. , 2. )
        self.h_mcpassmettrig  = ROOT.TH1F("h_mcpassmettrig" ,  "; passed met trigger",2 ,0. ,2. )
        self.h_mc_denominator =  ROOT.TH1F("h_mc_denominator" , "; E_{T}^{miss} [GeV]" , 21, 80., 500. )
        self.h_mc_numerator =  ROOT.TH1F("h_mc_numerator" , "; E_{T}^{miss} [GeV]" , 21, 80., 500. )
        self.addObject(self.h_mcpassreftrig )
        self.addObject(self.h_mc_denominator )
        self.addObject(self.h_mc_numerator )
        self.addObject(self.h_mcpassmettrig)

    def analyze(self, event):
        #met       = Object(event, "MET")
        #hlt       = Object(event, "HLT")
        Muons = Collection(event, 'Muon', 'nMuon')
        #goodMuons=[]
        goodMuons = filter(lambda j : (j.pt > 30) and (j.pfRelIso04_all<0.25) and (j.tightId),Muons)

        METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(event.MET_pt, event.MET_phi,event.run,"2016nonAPV",True,event.PV_npvs,True,False)
        METcorrected_pt = METcorrected_pt_phi[0]
        METcorrected_phi = METcorrected_pt_phi[1]

        #Throw away events that have corrected MET < 80 GeV and have zero good Muons
        if ((len(goodMuons)==0) or (METcorrected_pt<80)):
            return False

        #refTrigger = hlt.IsoMu24

        refTrigger = ((event.HLT_IsoMu24) or (event.HLT_IsoTkMu24))


        #self.h_mcpassreftrig.Fill(refTrigger,event.FinalWeighting)
        #self.h_mcpassreftrig.Fill(refTrigger)

        #if not refTrigger:
            #return False
        
        self.h_mc_denominator.Fill(METcorrected_pt,event.FinalWeighting)
        #self.h_mc_denominator.Fill(METcorrected_pt)

#        if(event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight 
#          or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight
#          or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight
#          or event.HLT_PFMET110_PFMHT110_IDTight
#          or event.HLT_PFMET120_PFMHT120_IDTight
#          or event.HLT_PFMET170_HBHECleaned
#          or event.HLT_PFMET170_HBHE_BeamHaloCleaned):
        if(refTrigger):
            self.h_mc_numerator.Fill(METcorrected_pt,event.FinalWeighting)
            self.h_mcpassmettrig.Fill(1,event.FinalWeighting)

        
        return True

data = lambda: dataEfficiencty()
mc = lambda: mcEfficiencty()
    

parser = argparse.ArgumentParser(description='Script create efficiency histograms for data and mc to be used to compute scalefactors')
#parser.add_argument('--dataPath',help="enter the path to data files",default="")
parser.add_argument('--mcPath',help="enter the path to mc files",default="")
args = parser.parse_args()



eventSelectionAND = [	"PV_ndof > 4",
						"abs(PV_z) < 24",
						"sqrt(PV_x*PV_x+PV_y*PV_y) < 2",
						"Flag_goodVertices",
						"Flag_globalSuperTightHalo2016Filter", 
						"Flag_HBHENoiseIsoFilter",
						"Flag_HBHENoiseFilter",
						"Flag_EcalDeadCellTriggerPrimitiveFilter",
						"Flag_BadPFMuonFilter",
						"Flag_eeBadScFilter"]

preselection = "&&".join(eventSelectionAND)

#fnames_data = glob.glob(args.dataPath + "/*.root")

fnames_mc = list(set(glob.glob(args.mcPath + "/*.root"))-set(glob.glob(args.mcPath + "/Radion*.root")))
print (fnames_mc)


###Crate and run post processor for 2016 MC
##print("Creating/Running post processor for 2016 MC")
###post2016_MC = PostProcessor("2016/.",fnames_mc,cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016MC(),mc()],histFileName="2016MC_MCHistos.root",histDirName="MCHistos")
##post2016_MC = PostProcessor("2016/.",fnames_mc,cut=preselection,branchsel=None,modules=[mc()],noOut=True,histFileName="2016MC_MCHistos.root",histDirName="MCHistos")
###post2016_MC = PostProcessor(".",["/hdfs/store/user/gparida/HHbbtt/Hadded_Skimmed_Files/Full_Production_August19_23_EleBitMapFix/2016/RadionTohhTohtatahbb_narrow_M-4500.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016MC(),mc()],histFileName="2016MC_MCHistos.root",histDirName="MCHistos")
##post2016_MC.run()

#Creating and running Post processor of 2016 Data
print("Creating/Running post processor for 2016 Data - SingleMu_2016F")
#post2016_RunF=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016F.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016RunF(),data()],noOut=True,histFileName="2016RunF_DataHistos.root",histDirName="DataHistos")
post2016_RunF=PostProcessor("2016/TestWithIdIso/.",[args.mcPath + "/SingleMu/SingleMu_Run2016F.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016/TestWithIdIso/2016RunF_DataHistos.root",histDirName="DataHistos")
post2016_RunF.run()

print("Creating/Running post processor for 2016 Data - SingleMu_2016G")
#post2016_RunG=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016G.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016RunF(),data()],noOut=True,histFileName="2016RunF_DataHistos.root",histDirName="DataHistos")
post2016_RunG=PostProcessor("2016/TestWithIdIso/.",[args.mcPath + "/SingleMu/SingleMu_Run2016G.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016/TestWithIdIso/2016RunG_DataHistos.root",histDirName="DataHistos")
post2016_RunG.run()

print("Creating/Running post processor for 2016 Data - SingleMu_2016H")
#post2016_RunH=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016H.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016RunH(),data()],noOut=True,histFileName="2016RunH_DataHistos.root",histDirName="DataHistos")
post2016_RunH=PostProcessor("2016/TestWithIdIso/.",[args.mcPath + "/SingleMu/SingleMu_Run2016H.root"],cut=preselection,branchsel=None,modules=[data()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json",noOut=True,histFileName="2016/TestWithIdIso/2016RunH_DataHistos.root",histDirName="DataHistos")
post2016_RunH.run()



