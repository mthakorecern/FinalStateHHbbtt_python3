#!/usr/bin/env python
#This script is used for files that are pre-processed to have a good Data-MC agreement in the SingleMuon Dataset
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

	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.out = wrappedOutputTree
		#self.out.branch("EventMass", "F")

	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		pass

	def analyze(self, event):
		FatJet = Collection(event,'FatJet','nFatJet')

		if ((FatJet[event.fatjet_index].msoftdrop<30) or (FatJet[event.fatjet_index].pt_nom<=200) or (event.ngood_MediumJets!=0)):
			return False


		self.h_data_denominator.Fill(event.METcorrected_pt)
		
		
		

		#if(refTrigger):
		#try:
		if (args.year=="2016"):
			#metTrigger = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned or event.HLT_PFMET170_HBHE_BeamHaloCleaned)
			#if(metTrigger and ((goodMuons[0].p4().DeltaR(FatJet[0].p4())) > 0.8)):
			if(event.met_trig_bool==1):
				self.h_data_numerator.Fill(event.METcorrected_pt)
		#except:
		elif (args.year=="2016APV"):
			#metTriggerSansBeamHaloCleaned = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned)
			#if(metTriggerSansBeamHaloCleaned and ((goodMuons[0].p4().DeltaR(FatJet[0].p4())) > 0.8)):
			if(event.met_trig_bool==1):
				self.h_data_numerator.Fill(event.METcorrected_pt)

		elif (args.year=="2017"):
			#metTrigger2017 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
			if(event.met_trig_bool==1):
				self.h_data_numerator.Fill(event.METcorrected_pt)

		elif (args.year=="2018"):
			#metTrigger2018 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
			if(event.met_trig_bool==1):
				self.h_data_numerator.Fill(event.METcorrected_pt)
			



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
		FatJet = Collection(event,'FatJet','nFatJet')        

		#print ("Index of fatjet = ",event.fatjet_index)
		if ((FatJet[event.fatjet_index].msoftdrop<30) or (FatJet[event.fatjet_index].pt_nom<=200) or (event.ngood_MediumJets!=0)):
		#if ((FatJet[event.fatjet_index[0]].msoftdrop<30) or (FatJet[event.fatjet_index[0]].pt_nom<=200)):
			return False


		self.h_mc_denominator.Fill(event.METcorrected_pt,event.FinalWeighting*event.pileupWeighting*event.combinedWZgenPtDeborahWeight*event.Top_pt_rwt)
		#self.h_mc_denominator.Fill(METcorrected_pt)
		#try:
		if (args.year=="2016"):
			#metTrigger = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned or event.HLT_PFMET170_HBHE_BeamHaloCleaned)
			if(event.met_trig_bool==1):
				#self.h_mc_numerator.Fill(METcorrected_pt,event.FinalWeighting)
				self.h_mc_numerator.Fill(event.METcorrected_pt,event.FinalWeighting*event.pileupWeighting*event.combinedWZgenPtDeborahWeight*event.Top_pt_rwt)
		#except:
		elif (args.year=="2016APV"): 
			#metTriggerSansBeamHaloCleaned = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned)
			if(event.met_trig_bool==1):
				#self.h_mc_numerator.Fill(METcorrected_pt,event.FinalWeighting)
				self.h_mc_numerator.Fill(event.METcorrected_pt,event.FinalWeighting*event.pileupWeighting*event.combinedWZgenPtDeborahWeight*event.Top_pt_rwt)

		elif (args.year=="2017"):
			#metTrigger2017 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
			if(event.met_trig_bool==1):
				#self.h_mc_numerator.Fill(METcorrected_pt,event.FinalWeighting)
				self.h_mc_numerator.Fill(event.METcorrected_pt,event.FinalWeighting*event.pileupWeighting*event.combinedWZgenPtDeborahWeight*event.Top_pt_rwt)                            
		
		elif (args.year=="2018"):
			#metTrigger2018 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
			if(event.met_trig_bool==1):
				#self.h_mc_numerator.Fill(METcorrected_pt,event.FinalWeighting)
				self.h_mc_numerator.Fill(event.METcorrected_pt,event.FinalWeighting*event.pileupWeighting*event.combinedWZgenPtDeborahWeight*event.Top_pt_rwt)      

		return True


	

parser = argparse.ArgumentParser(description='Script create efficiency histograms for data and mc to be used to compute scalefactors (currently wotj )')
#parser.add_argument('--dataPath',help="enter the path to data files",default="")
parser.add_argument('--mcPath',help="enter the path to mc files",default="")
parser.add_argument('--year','-y',help="Enter the year of data taking",choices=["2016","2016APV","2017","2018"], required=True)
args = parser.parse_args()

#fnames_data = glob.glob(args.dataPath + "/*.root")

data = lambda: dataEfficiencty()
mc = lambda: mcEfficiencty(args.year)

#Consider all the MC samples except the signal samples
fnames_mc = list(set(glob.glob(args.mcPath + "/*.root"))-set(glob.glob(args.mcPath + "/*Radion*.root"))-set(glob.glob(args.mcPath + "/*Graviton*.root"))-set(glob.glob(args.mcPath + "/*SingleMu_*.root"))-set(glob.glob(args.mcPath + "/*MET*.root")))
#fnames_mc = glob.glob(args.mcPath + "/WJetsToLNu_HT-400To600*.root")
print ([(nameStrip.strip().split('/')[-1]).split('.')[-2] for nameStrip in fnames_mc])



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
	post2016_MC = PostProcessor("2016/.",fnames_mc,cut=None,branchsel=None,modules=[mc()],noOut=True,histFileName="2016/2016MC_MCHistos.root",histDirName="MCHistos")    
	#post2016_MC = PostProcessor(".",["/hdfs/store/user/gparida/HHbbtt/Hadded_Skimmed_Files/Full_Production_August19_23_EleBitMapFix/2016/RadionTohhTohtatahbb_narrow_M-4500.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016MC(),mc()],histFileName="2016MC_MCHistos.root",histDirName="MCHistos")
	post2016_MC.run()


	#Creating and running Post processor of 2016 Data
	print("Creating/Running post processor for 2016 Data - SingleMu_2016F")
	#post2016_RunF=PostProcessor("2016/.",[args.mcPath + "SingleMu/SingleMu/SingleMu_Run2016F.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016RunF(),data()],noOut=True,histFileName="2016RunF_DataHistos.root",histDirName="DataHistos")
	post2016_RunF=PostProcessor("2016/.",glob.glob(args.mcPath + "/SingleMu_Run2016F*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2016/2016RunF_DataHistos.root",histDirName="DataHistos")

	post2016_RunF.run()

	print("Creating/Running post processor for 2016 Data - SingleMu_2016G")
	#post2016_RunG=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu_Run2016G.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016RunF(),data()],noOut=True,histFileName="2016RunF_DataHistos.root",histDirName="DataHistos")
	post2016_RunG=PostProcessor("2016/.",glob.glob(args.mcPath + "/SingleMu_Run2016G*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2016/2016RunG_DataHistos.root",histDirName="DataHistos")
	post2016_RunG.run()

	print("Creating/Running post processor for 2016 Data - SingleMu_2016H")
	#post2016_RunH=PostProcessor("2016/.",[args.mcPath + "/SingleMu/SingleMu/SingleMu_Run2016H.root"],cut=preselection,branchsel=None,modules=[jetmetCorrectorUL2016RunH(),data()],noOut=True,histFileName="2016RunH_DataHistos.root",histDirName="DataHistos")
	post2016_RunH=PostProcessor("2016/.",glob.glob(args.mcPath + "/SingleMu_Run2016H*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2016/2016RunH_DataHistos.root",histDirName="DataHistos")
	post2016_RunH.run()

elif (args.year=="2016APV"):
	#Create and run post processor for 2016APV MC
	print("Creating/Running post processor for 2016APV MC")
	APV2016_MC = PostProcessor("2016APV/.",fnames_mc,cut=None,branchsel=None,modules=[mc()],noOut=True,histFileName="2016APV/2016APVMC_MCHistos.root",histDirName="MCHistos")
	APV2016_MC.run()

	#Creating and running Post processor of 2016APV Run2016B-ver1_HIPM  Data
	print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016B-ver1_HIPM")
	Run2016B_ver1_HIPM=PostProcessor("2016APV/.",glob.glob(args.mcPath + "/SingleMu_Run2016B-ver1_HIPM*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2016APV/Run2016B_ver1_HIPM_DataHistos.root",histDirName="DataHistos")
	Run2016B_ver1_HIPM.run()

	#Creating and running Post processor of 2016APV Run2016B-ver2_HIPM  Data
	print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016B-ver2_HIPM")
	Run2016B_ver2_HIPM=PostProcessor("2016APV/.",glob.glob(args.mcPath + "/SingleMu_Run2016B-ver2_HIPM*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2016APV/Run2016B_ver2_HIPM_DataHistos.root",histDirName="DataHistos")
	Run2016B_ver2_HIPM.run()

	#Creating and running Post processor of 2016APV Run2016C-HIPM  Data
	print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016C-HIPM")
	Run2016C_HIPM=PostProcessor("2016APV/.",glob.glob(args.mcPath + "/SingleMu_Run2016C-HIPM*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2016APV/Run2016C_HIPM_DataHistos.root",histDirName="DataHistos")
	Run2016C_HIPM.run()

	#Creating and running Post processor of 2016APV Run2016D-HIPM  Data
	print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016D-HIPM")
	Run2016D_HIPM=PostProcessor("2016APV/.",glob.glob(args.mcPath + "/SingleMu_Run2016D-HIPM*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2016APV/Run2016D_HIPM_DataHistos.root",histDirName="DataHistos")
	Run2016D_HIPM.run()

	#Creating and running Post processor of 2016APV Run2016E-HIPM  Data
	print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016E-HIPM")
	Run2016E_HIPM=PostProcessor("2016APV/.",glob.glob(args.mcPath + "/SingleMu_Run2016E-HIPM*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2016APV/Run2016E_HIPM_DataHistos.root",histDirName="DataHistos")
	Run2016E_HIPM.run()

	#Creating and running Post processor of 2016APV Run2016F-HIPM  Data
	print("Creating/Running post processor for 2016APV Data - SingleMu_Run2016F-HIPM")
	Run2016F_HIPM=PostProcessor("2016APV/.",glob.glob(args.mcPath + "/SingleMu_Run2016F-HIPM*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2016APV/Run2016F_HIPM_DataHistos.root",histDirName="DataHistos")
	Run2016F_HIPM.run()

elif (args.year=="2018"):
	print("Creating/Running post processor for 2018 MC")
	post2018_MC = PostProcessor("2018/.",fnames_mc,cut=None,branchsel=None,modules=[mc()],noOut=True,histFileName="2018/2018MC_MCHistos.root",histDirName="MCHistos")
	post2018_MC.run()  

	print("Creating/Running post processor for 2018 Data - SingleMu_Run2018A")
	Run2018A = PostProcessor("2018/.",glob.glob(args.mcPath + "/SingleMu_Run2018A*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2018/Run2018A_DataHistos.root",histDirName="DataHistos")
	Run2018A.run()        

	print("Creating/Running post processor for 2018 Data - SingleMu_Run2018B")
	Run2018B = PostProcessor("2018/.",glob.glob(args.mcPath + "/SingleMu_Run2018B*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2018/Run2018B_DataHistos.root",histDirName="DataHistos")
	Run2018B.run() 

	print("Creating/Running post processor for 2018 Data - SingleMu_Run2018C")
	Run2018C = PostProcessor("2018/.",glob.glob(args.mcPath + "/SingleMu_Run2018C*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2018/Run2018C_DataHistos.root",histDirName="DataHistos")
	Run2018C.run()     

	print("Creating/Running post processor for 2018 Data - SingleMu_Run2018D")
	#data2018D_files = glob.glob(args.mcPath + "/SingleMu/SingleMu/SingleMu_Run2018D*.root")
	Run2018D = PostProcessor("2018/.",glob.glob(args.mcPath + "/SingleMu_Run2018D*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2018/Run2018D_DataHistos.root",histDirName="DataHistos")
	Run2018D.run()     


elif (args.year=="2017"):
	print("Creating/Running post processor for 2017 MC")
	post2017_MC = PostProcessor("2017/.",fnames_mc,cut=None,branchsel=None,modules=[mc()],noOut=True,histFileName="2017/2017MC_MCHistos.root",histDirName="MCHistos")
	post2017_MC.run()  

	print("Creating/Running post processor for 2017 Data - SingleMu_Run2017B")
	Run2017B = PostProcessor("2017/.",glob.glob(args.mcPath + "/SingleMu_Run2017B*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2017/Run2017B_DataHistos.root",histDirName="DataHistos")
	Run2017B.run() 

	print("Creating/Running post processor for 2017 Data - SingleMu_Run2017C")
	Run2017C = PostProcessor("2017/.",glob.glob(args.mcPath + "/SingleMu_Run2017C*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2017/Run2017C_DataHistos.root",histDirName="DataHistos")
	Run2017C.run() 

	print("Creating/Running post processor for 2017 Data - SingleMu_Run2017D")
	Run2017D = PostProcessor("2017/.",glob.glob(args.mcPath + "/SingleMu_Run2017D*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2017/Run2017D_DataHistos.root",histDirName="DataHistos")
	Run2017D.run() 

	print("Creating/Running post processor for 2017 Data - SingleMu_Run2017E")
	Run2017E = PostProcessor("2017/.",glob.glob(args.mcPath + "/SingleMu_Run2017E*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2017/Run2017E_DataHistos.root",histDirName="DataHistos")
	Run2017E.run() 

	print("Creating/Running post processor for 2017 Data - SingleMu_Run2017F")
	Run2017F = PostProcessor("2017/.",glob.glob(args.mcPath + "/SingleMu_Run2017F*.root"),cut=None,branchsel=None,modules=[data()],noOut=True,histFileName="2017/Run2017F_DataHistos.root",histDirName="DataHistos")
	Run2017F.run() 