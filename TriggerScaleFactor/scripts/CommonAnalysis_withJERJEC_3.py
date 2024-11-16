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
from re import search
import multiprocessing as  mp
from FinalStateHHbbtt.Pre_BackgroundEstimation.Corrections.WandZptWeightingModule.wandzgenptweightProcessor import wandzgenptweight

import argparse
#Import and load the MET-XY correction header file for UL
#MET phi correction: https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection_withUL17andUL18andUL16.h, https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections#xy_Shift_Correction_MET_phi_modu


class mettrigsf(Module):
	def __init__(self,filename,year=2016, isData=False):
		self.filename = filename
		self.year = year
		self.isData = isData
		self.isMC = not self.isData
		if self.isMC:
			if self.year == '2016': 
				self.LHCLumi =  16.81e15
			elif self.year == '2016APV':
				self.LHCLumi = 19.52e15
			elif self.year == '2017':
				self.LHCLumi = 41.48e15
			elif self.year == '2018':
				self.LHCLumi = 59.83e15 
		self.jetFV = ROOT.TLorentzVector(0.0,0.0,0.0,0.0)   
		self.higgsBBFV = ROOT.TLorentzVector(0.0,0.0,0.0,0.0)
		self.muonFV = ROOT.TLorentzVector(0.0,0.0,0.0,0.0)

		if self.year=="2016":
			self.LooseJet = 0.0480
			self.MediumJet = 0.2489
			self.TightJet = 0.6377
		elif self.year=="2016APV":
			self.LooseJet = 0.0508
			self.MediumJet = 0.2598
			self.TightJet = 0.6502
		elif self.year=="2017":
			self.LooseJet = 0.0532
			self.MediumJet = 0.3040
			self.TightJet = 0.7476
		elif self.year=="2018":
			self.LooseJet = 0.0490
			self.MediumJet = 0.2783
			self.TightJet =  0.7100

	def beginJob(self):
		if (self.isMC == True):
			self.xsJSON = os.environ['CMSSW_BASE']+"/src/MetaData/xsJSON/%s_Samples.json"%(self.year) #Add a customization to different yrs
			with open(self.xsJSON,'r') as xsjson:
				self.xsInfo = json.load(xsjson)


			self.countsJSON = os.environ['CMSSW_BASE']+'/src/MetaData/SumOfGenWeights/%s_weight.json'%(self.year) #Add a customization to different yrs
			with open(self.countsJSON,'r') as countsjson:
				self.countsInfo = json.load(countsjson)

	def endJob(self):
		pass

	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		#Define output branches to check if preselection is working
		self.out = wrappedOutputTree 
		self.out.branch("met_trig_bool","I")
		#self.out.branch("METcorrected_pt","F")
		#self.out.branch("METcorrected_phi","F")
		self.out.branch("delR_muon_fatjet","F")
		self.out.branch("fatjet_index","I")
		self.out.branch("muon_index","I")
		self.out.branch("ngood_Jets","I")
		self.out.branch("ngood_LooseJets","I")
		self.out.branch("ngood_MediumJets","I")
		self.out.branch("ngood_TightJets","I")
		self.out.branch("index_gJets","I",lenVar="ngood_Jets")
		self.out.branch("index_gLooseJets","I",lenVar="ngood_LooseJets")
		self.out.branch("index_gMediumJets","I",lenVar="ngood_MediumJets")
		self.out.branch("index_gTightJets","I",lenVar="ngood_TightJets")
		self.out.branch("FinalWeighting","F")
		self.out.branch("MT","F")
		self.out.branch("delphi_MET_L","F")
		self.out.branch("Top_pt_rwt","F")
		self.out.branch("Top_wt","F")
		self.out.branch("genTop_pt","F")
		self.out.branch("antiTop_wt","F")
		self.out.branch("genantiTop_pt","F")

		if (self.isMC==True):
			print(" The inputFile name gotten from the beginFile fn : ",type(inputFile),"   ",inputFile.GetName())
			filenameForJson = (inputFile.GetName().strip().split('/')[-1]).split('.')[-2]
			self.filename = filenameForJson
			self.XS = self.xsInfo[self.filename]['XS'] * 1e-12
			self.totalNumberOfEvents = self.countsInfo[self.filename]
			print(" Filename for JSON : ",self.filename," XS = ",self.XS," sumofgenwts = ",self.totalNumberOfEvents)


	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		pass

	def analyze(self, event):

		Jet = Collection(event,'Jet','nJet')
		FatJet = Collection(event,'FatJet','nFatJet')

 		def applyPOGselectionToAK4_for2018Veto(ak4Object_enu):
			if (((ak4Object_enu[1].pt_nom>15) and (ak4Object_enu[1].pt_nom<50)) and (ak4Object_enu[1].eta < -1.3) and (ak4Object_enu[1].eta > -3.2) and (ak4Object_enu[1].phi < -0.87) and (ak4Object_enu[1].phi > -1.57)):
				if ((ak4Object_enu[1].jetId>1) and (ak4Object_enu[1].puId>=4)):
					return True
			if ((ak4Object_enu[1].pt_nom>=50) and (ak4Object_enu[1].eta < -1.3) and (ak4Object_enu[1].eta > -3.2) and (ak4Object_enu[1].phi < -0.87) and (ak4Object_enu[1].phi > -1.57)):
				if (ak4Object_enu[1].jetId>1):
					return True
			return False
 		def applyPOGselectionToAK4(ak4Object_enu):
			if (((ak4Object_enu[1].pt_nom>30) and (ak4Object_enu[1].pt_nom<50)) and (abs(ak4Object_enu[1].eta)<2.5)):
				if ((ak4Object_enu[1].jetId>1) and (ak4Object_enu[1].puId>=4)):
					return True
			if ((ak4Object_enu[1].pt_nom>=50) and (abs(ak4Object_enu[1].eta)<2.5)):
				if (ak4Object_enu[1].jetId>1):
					return True
			return False

		if (args.year == "2018"):
			#print ("It comes to the year part",event.run," Data = ",self.isData," MC = ",self.isMC)
			if ((self.isData==True) and (event.run>=319077)):
				#print ("Entering the Data section")
				Jet_Vetoenu = filter(applyPOGselectionToAK4_for2018Veto,enumerate(Jet))
				FatJet_Vetoenu = filter(lambda x: (x[1].pt_nom > 170) and (x[1].eta < -1.1) and (x[1].eta > -2.7) and (x[1].phi < -0.67) and (x[1].phi > -1.77) and (x[1].jetId>1), enumerate(FatJet))				
				if ((len(FatJet_Vetoenu)>0) or (len(Jet_Vetoenu)>0)):
					return False    
			if ((self.isMC==True)):
				#print ("Entering the MC section")
				frac = np.random.uniform(0.0, 1.0)
				if (frac>0.352275515):
					Jet = Collection(event,'Jet','nJet')
					Jet_Vetoenu = filter(applyPOGselectionToAK4_for2018Veto,enumerate(Jet))
					FatJet_Vetoenu = filter(lambda x: (x[1].pt_nom > 170) and (x[1].eta < -1.1) and (x[1].eta > -2.7) and (x[1].phi < -0.67) and (x[1].phi > -1.77) and (x[1].jetId>1), enumerate(FatJet))
					if ((len(FatJet_Vetoenu)>0) or (len(Jet_Vetoenu)>0)):
						return False            			

		def JetFatJetOverlap(jetObject_enu): #FatJet and Jet overlap, separation > 1.2 (0.8 + 0.4)
			self.jetFV.SetPtEtaPhiM(jetObject_enu[1].pt_nom,jetObject_enu[1].eta,jetObject_enu[1].phi,jetObject_enu[1].mass_nom)
			deltaR = self.jetFV.DeltaR(self.higgsBBFV)
			if deltaR > 1.2:
				return True
			else:
				return False

		def FatJetMuonOverlap(muonObject_enu): #FatJet and Jet overlap, separation > 1.2 (0.8 + 0.4)
			self.muonFV.SetPtEtaPhiM(muonObject_enu[1].pt,muonObject_enu[1].eta,muonObject_enu[1].phi,muonObject_enu[1].mass)
			deltaR = self.muonFV.DeltaR(self.higgsBBFV)
			if deltaR > 0.8:
				return True
			else:
				return False


		Muons = Collection(event, 'Muon', 'nMuon')

		Muons_enu = filter(lambda j : ((j[1].pt > 30) and (j[1].pfRelIso04_all<0.15) and (j[1].tightId) and (abs(j[1].eta)<2.4)),enumerate(Muons))
		FatJet_enu = filter(lambda x : ((x[1].pt > 170) and (abs(x[1].eta) < 2.5) and (x[1].jetId>1)),enumerate(FatJet))
		if (len(FatJet_enu)==0):
			return False
		self.higgsBBFV.SetPtEtaPhiM(FatJet_enu[0][1].pt_nom,FatJet_enu[0][1].eta,FatJet_enu[0][1].phi,FatJet_enu[0][1].mass_nom)
		Jet_enu =  filter(applyPOGselectionToAK4,enumerate(Jet))
		Jet_enu = filter(JetFatJetOverlap,Jet_enu)

		Muons_enu = filter(FatJetMuonOverlap,Muons_enu)

		if ((len(Muons_enu)!=1)):
			return False

		#METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, TString year, bool isMC, int npv, bool isUL =false,bool ispuppi=false)
		#if(isMC && year == "2016" && !isUL) runera = y2016MC;
		#else if(isMC && year == "2017" && !isUL) {runera = y2017MC; usemetv2 =true;}
		#else if(isMC && year == "2018" && !isUL) runera = y2018MC;
		#else if(isMC && year == "2016APV" && isUL) runera = yUL2016MCAPV;
		#else if(isMC && year == "2016nonAPV" && isUL) runera = yUL2016MCnonAPV;
		#else if(isMC && year == "2017" && isUL) runera = yUL2017MC;
		#else if(isMC && year == "2018" && isUL) runera = yUL2018MC;

		#Throw away events that have corrected MET < 80 GeV and have zero good Muons
		#if ((len(goodMuons)!=1) or (METcorrected_pt<150)):



		#if (args.year=="2016" or args.year=="2016APV"):
		#	refTrigger = ((event.HLT_IsoMu24) or (event.HLT_IsoTkMu24))
		#else:
		#	refTrigger = ((event.HLT_IsoMu24))

		#self.h_datapassreftrig.Fill(refTrigger)



		#if not refTrigger:
		#	return False

		#for jet in goodFatJet:



		if (self.isMC==True):
			self.out.fillBranch("FinalWeighting",(self.XS * self.LHCLumi * event.genWeight) / self.totalNumberOfEvents)

		else:
			self.out.fillBranch("FinalWeighting",1.0)

		delphi_MET_L = Muons_enu[0][1].phi - event.METcorrected_phi 

		if (delphi_MET_L > math.pi):
			delphi_MET_L -= 2 * math.pi
		if (delphi_MET_L < (-math.pi)):
			delphi_MET_L += 2 * math.pi

		#print ("delphi_MET_L = ",delphi_MET_L)
		MT = math.sqrt(2 * Muons_enu[0][1].pt * event.METcorrected_pt * (1 - math.cos(delphi_MET_L)))

		self.out.fillBranch("MT",MT)
		self.out.fillBranch("delphi_MET_L",delphi_MET_L)
		self.out.fillBranch("delR_muon_fatjet",Muons_enu[0][1].p4().DeltaR(self.higgsBBFV))
		self.out.fillBranch("fatjet_index",FatJet_enu[0][0])
		self.out.fillBranch("muon_index",Muons_enu[0][0])
		
		self.out.fillBranch("ngood_Jets",len(Jet_enu))
		self.out.fillBranch("index_gJets",[x[0] for x in Jet_enu])

		Jet_enu_Loose = filter(lambda x:x[1].btagDeepFlavB >= self.LooseJet,Jet_enu)
		self.out.fillBranch("ngood_LooseJets",len(Jet_enu_Loose))
		self.out.fillBranch("index_gLooseJets",[x[0] for x in Jet_enu_Loose])

		Jet_enu_Medium = filter(lambda x:x[1].btagDeepFlavB >= self.MediumJet,Jet_enu)
		self.out.fillBranch("ngood_MediumJets",len(Jet_enu_Medium))
		self.out.fillBranch("index_gMediumJets",[x[0] for x in Jet_enu_Medium])
	
		Jet_enu_Tight = filter(lambda x:x[1].btagDeepFlavB >= self.TightJet,Jet_enu)
		self.out.fillBranch("ngood_TightJets",len(Jet_enu_Tight))
		self.out.fillBranch("index_gTightJets",[x[0] for x in Jet_enu_Tight])

		SF_top = 1.0
		SF_antitop = 1.0
		#Compute the top reweighting
		if search("TTTo",self.filename):
			#print ("This is a TTbar sample : ",self.filename)
			genParticles = Collection(event, "GenPart")
			#Filter the genPart collection to pick out top (pdgid is 6) and make sure top quark in MC is defined at parton level, after radiation and before decay, using the 'isLastCopy' i.e status 62
			top = list(filter(lambda gen : (gen.pdgId == 6) and (gen.status == 62), genParticles))
			if len(top) > 0:
				SF_top = 0.103*np.exp(-0.0118*top[0].pt) - 0.000134*top[0].pt + 0.973
			
			#Repeat the same for antitop
			antitop = list(filter(lambda gen : (gen.pdgId == -6) and (gen.status == 62), genParticles))
			if len(antitop) > 0:
				SF_antitop = 0.103*np.exp(-0.0118*antitop[0].pt) - 0.000134*antitop[0].pt + 0.973
		
			#This branch is used for plotting "Top_pt_rwt"
			self.out.fillBranch("Top_pt_rwt",np.sqrt(SF_top*SF_antitop))
			self.out.fillBranch("Top_wt",SF_top)
			self.out.fillBranch("antiTop_wt",SF_antitop)
			self.out.fillBranch("genTop_pt",top[0].pt)
			self.out.fillBranch("genantiTop_pt",antitop[0].pt)

		else:
			self.out.fillBranch("Top_pt_rwt",1.0)
			self.out.fillBranch("Top_wt",1.0)
			self.out.fillBranch("antiTop_wt",1.0)
			self.out.fillBranch("genTop_pt",-99.99)
			self.out.fillBranch("genantiTop_pt",-99.99)		
		

		#if(refTrigger):
		#try:
		if (args.year=="2016"):
			#metTrigger = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned or event.HLT_PFMET170_HBHE_BeamHaloCleaned)
			#Following the OR (Trig):https://cms-pub-talk.web.cern.ch/t/trigger-object-review/30610
			#We are ONLY keeping METNoMu triggers
			metTrigger = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight) 
			if(metTrigger):
				self.out.fillBranch("met_trig_bool",1)
			else:
				self.out.fillBranch("met_trig_bool",0)
			return True


				
		#except:
		elif (args.year=="2016APV"):
			#metTriggerSansBeamHaloCleaned = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned)
			#metTrigger = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned or event.HLT_PFMET170_HBHE_BeamHaloCleaned)
			#Following the OR (Trig):https://cms-pub-talk.web.cern.ch/t/trigger-object-review/30610
			#We are ONLY keeping METNoMu triggers
			#if(metTriggerSansBeamHaloCleaned and ((goodMuons[0].p4().DeltaR(FatJet[0].p4())) > 0.8)):
			metTriggerSansBeamHaloCleaned = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight)
			if(metTriggerSansBeamHaloCleaned):
				self.out.fillBranch("met_trig_bool",1)
			else:
				self.out.fillBranch("met_trig_bool",0)
			return True


		elif (args.year=="2017"):
			#metTrigger2017 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
			#metTrigger = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned or event.HLT_PFMET170_HBHE_BeamHaloCleaned)
			#Following the OR (Trig):https://cms-pub-talk.web.cern.ch/t/trigger-object-review/30610
			#We are ONLY keeping METNoMu triggers
			metTrigger2017 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight)
			if(metTrigger2017):
				self.out.fillBranch("met_trig_bool",1)
			else:
				self.out.fillBranch("met_trig_bool",0)
			return True

		elif (args.year=="2018"):
			#metTrigger2018 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET120_PFMHT120_IDTight)
			#metTrigger = (event.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight or event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_PFMET110_PFMHT110_IDTight or event.HLT_PFMET120_PFMHT120_IDTight or event.HLT_PFMET170_HBHECleaned or event.HLT_PFMET170_HBHE_BeamHaloCleaned)
			#Following the OR (Trig):https://cms-pub-talk.web.cern.ch/t/trigger-object-review/30610
			#We are ONLY keeping METNoMu triggers
			metTrigger2018 = (event.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight or event.HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight)
			if(metTrigger2018):
				self.out.fillBranch("met_trig_bool",1)
			else:
				self.out.fillBranch("met_trig_bool",0)
			return True

def call_postpoc(files):
	nameStrip=files.strip()
	filename = (nameStrip.split('/')[-1]).split('.')[-2]
	print(filename)
	
	if (search("Run", filename)):
		print ("This is a ",args.year," Data file = ",filename)
		mainModule = lambda: mettrigsf(filename,args.year,True)
		wandzWtModule = lambda: wandzgenptweight(filename,args.year,True)
	else:
		print ("This is a ",args.year," MC file = ", filename)
		mainModule = lambda: mettrigsf(filename,args.year,False)
		wandzWtModule = lambda: wandzgenptweight(filename,args.year,False)

	if (args.year=="2018"):
		if (search("Run", filename)):
			print ("This is Data - Use Golden Json file")	
			p = PostProcessor(args.tempStore,[files], cut="(HLT_IsoMu24) && (METcorrected_pt>=150)", branchsel="input_drop.txt",modules=[mainModule(),wandzWtModule()],jsonInput="GOLDEN_JSONS/2018_Leagcy.json", postfix="",noOut=False,outputbranchsel="output_drop.txt")
		else:
			p = PostProcessor(args.tempStore,[files], cut="(HLT_IsoMu24) && (METcorrected_pt>=150)", branchsel="input_drop.txt",modules=[mainModule(),wandzWtModule()], postfix="",noOut=False,outputbranchsel="output_drop.txt")

	if (args.year=="2017"):
		if (search("Run", filename)):
			print ("This is Data - Use Golden Json file")	
			p = PostProcessor(args.tempStore,[files], cut="(HLT_IsoMu24) && (METcorrected_pt>=150)", branchsel="input_drop.txt",modules=[mainModule(),wandzWtModule()],jsonInput="GOLDEN_JSONS/2017_Leagcy.json", postfix="",noOut=False,outputbranchsel="output_drop.txt")
		else:
			p = PostProcessor(args.tempStore,[files], cut="(HLT_IsoMu24) && (METcorrected_pt>=150)", branchsel="input_drop.txt",modules=[mainModule(),wandzWtModule()], postfix="",noOut=False,outputbranchsel="output_drop.txt")

	if (args.year=="2016APV"):
		if (search("Run", filename)):
			print ("This is Data - Use Golden Json file")	
			p = PostProcessor(args.tempStore,[files], cut="((HLT_IsoMu24) || (HLT_IsoTkMu24)) && (METcorrected_pt>=150)", branchsel="input_drop.txt",modules=[mainModule(),wandzWtModule()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json", postfix="",noOut=False,outputbranchsel="output_drop.txt")
		else:
			p = PostProcessor(args.tempStore,[files], cut="((HLT_IsoMu24) || (HLT_IsoTkMu24)) && (METcorrected_pt>=150)", branchsel="input_drop.txt",modules=[mainModule(),wandzWtModule()], postfix="",noOut=False,outputbranchsel="output_drop.txt")

	if (args.year=="2016"):
		if (search("Run", filename)):
			print ("This is Data - Use Golden Json file")	
			p = PostProcessor(args.tempStore,[files], cut="((HLT_IsoMu24) || (HLT_IsoTkMu24)) && (METcorrected_pt>=150)", branchsel="input_drop.txt",modules=[mainModule(),wandzWtModule()],jsonInput="GOLDEN_JSONS/2016_Leagcy.json", postfix="",noOut=False,outputbranchsel="output_drop.txt")
		else:
			p = PostProcessor(args.tempStore,[files], cut="((HLT_IsoMu24) || (HLT_IsoTkMu24)) && (METcorrected_pt>=150)", branchsel="input_drop.txt",modules=[mainModule(),wandzWtModule()], postfix="",noOut=False,outputbranchsel="output_drop.txt")


	p.run()	
	print("###############MOVING THE OUTPUT FILE BACK TO HDFS#######################")
	#os.system("hadoop fs -moveFromLocal -f "+args.tempStore+"/"+filename+".root"+" "+outputDir+"/.") #This currently doesnot work due to username differences - it takes parida by default
	os.system("mv "+args.tempStore+"/"+filename+".root"+" "+outputDir+"/.")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Output of this script will be used for computing MET SFs. These outfiles are stored under (/hdfs/store/user/gparida/HHbbtt/Hadded_Skimmed_Files/Full_Production_CMSSW_13_0_13_Nov24_23/JetCorrBranchesAndMETSkimTrim_METSFs)')	
	parser.add_argument('--inputLocation','-i',help="enter the path to the location of input file set",default="")
	parser.add_argument('--outputLocation','-o',help="enter the path where yu want the output files to be stored",default =".")
	parser.add_argument('--ncores','-n',help ="number of cores for parallel processing", default=1)
	parser.add_argument('--year','-y',help='specify the run - to make sure right triggers are used',choices=['2016','2016APV','2017','2018'])
	parser.add_argument('--tempStore','-t',help='Temporary staging area for files before moving out to hdfs', required=True)

	args = parser.parse_args()

	#Input file location typically:
	#fnames = glob.glob(args.inputLocation + "/RadionTohhTohtatahbb_narrow_M-1000*.root")  #making a list of input MC/DATA files
	#if (args.data):
	#	fnames = glob.glob(args.inputLocation + "/SingleMu/*.root")  #making a list of input MC/DATA files
	#else:
	#Remove the Radion and Graviton Signal files
	fnames = list(set(glob.glob(args.inputLocation + "/*.root"))-set(glob.glob(args.inputLocation + "/*Radion*.root"))-set(glob.glob(args.inputLocation + "/*Graviton*.root")))

	outputDir = args.outputLocation

	argList = list()
	for file in fnames:
		argList.append(file)

	if int(args.ncores) == 1:
		for arr in argList:
			print ("Using a single thread ")
			call_postpoc(arr)
	
	else:
		print ("Using Multithreading")
		pool = mp.Pool(int(args.ncores))
		print ("list", argList)
		res=pool.map(call_postpoc, argList)



		











		






		
			
		






			



					





		
	




		
		
		




