

# Copied from the parent script #CommonAnalysis_4.py. I will list the changes that goes in this
# I need to implement:
# 200 GeV threshold for FatJets
# Electron transition veto
# Tau Energy Scale (nominal and variations if possible)
import math
import argparse
import multiprocessing as mp
from re import search
import numpy as np
import glob
from collections import OrderedDict
import json
from FinalStateHHbbtt.Pre_BackgroundEstimationWithSystematics.TauEnergyScaleModule.TauEnergyScaleForHPSandBoosted import TauEnergyScaleForHPSandBoosted
from FinalStateHHbbtt.fastMTTPython.fastMTTtool import *
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import time
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
# from FinalStateHHbbtt.Pre_BackgroundEstimation.filterSDmass_HbbtaggingAndWeights import filterSDmass_HbbtaggingAndWeights
# from FinalStateHHbbtt.Pre_BackgroundEstimation.Corrections.WandZptWeightingModule.wandzgenptweightProcessor import wandzgenptweight
# from FinalStateHHbbtt.Pre_BackgroundEstimation.AddEventHTandWtransversemass import HT_MT_vars
# from FinalStateHHbbtt.Pre_BackgroundEstimation.AddBunchofDelRs import AddBunchofDelRs
# from FinalStateHHbbtt.Pre_BackgroundEstimation.Corrections.topPtWeighingModule.topPtprocessor import topPtweight
# from FinalStateHHbbtt.Pre_BackgroundEstimation.TauEnergyScaleForHPS import TauCorrectionProducer
# from FinalStateHHbbtt.Pre_BackgroundEstimation.OrthogonalityChecks.HadronicTauToAK8Jet.matchingAK8ToTaus_semileptonic import AddAK8MatchingToTaus
# from FinalStateHHbbtt.Pre_BackgroundEstimation.OrthogonalityChecks.Standalone_AdditionalLepVetoFlag.Standalone_AdditionalLepVetoFlag import Standalone_AdditionalLepVetoFlag


class cutsAndcategories(Module):
    def __init__(self, filename, year, isData, runNominal=False):
        # self.writeHistFile=True
        self.runNominal = runNominal
        print("ACTIVATE: cutsAndcategories")
        self.year = year
        if ((self.year == "2016APV") or (self.year == "2016")):
            self.year_unc = "2016"
        elif ((self.year == "2017") or (self.year == "2018")):
            self.year_unc = year
        self.isMC = not isData
        self.isData = isData
        self.filename = filename
        # self. jesUnc = ["","jesUp","jesDown","jerUp","jerDown","tesUp","tesDown"]

        if ((self.isData) or (self.runNominal)):
            self.jesUnc = [""]
        else:
            self.jesUnc = [
                "",
                "jesTotalUp",
                "jesTotalDown",
                "jerUp",
                "jerDown",
                "jesAbsoluteUp",
                "jesAbsoluteDown",
                "jesAbsolute_%sUp" %
                (self.year_unc),
                "jesAbsolute_%sDown" %
                (self.year_unc),
                "jesBBEC1Up",
                "jesBBEC1Down",
                "jesBBEC1_%sUp" %
                (self.year_unc),
                "jesBBEC1_%sDown" %
                (self.year_unc),
                "jesEC2Up",
                "jesEC2Down",
                "jesEC2_%sUp" %
                (self.year_unc),
                "jesEC2_%sDown" %
                (self.year_unc),
                "jesFlavorQCDUp",
                "jesFlavorQCDDown",
                "jesHFUp",
                "jesHFDown",
                "jesHF_%sUp" %
                (self.year_unc),
                "jesHF_%sDown" %
                (self.year_unc),
                "jesRelativeBalUp",
                "jesRelativeBalDown",
                "jesRelativeSample_%sUp" %
                (self.year_unc),
                "jesRelativeSample_%sDown" %
                (self.year_unc),
                "UnclustUp",
                "UnclustDown",
                "tesUp",
                "tesDown"]

        # Instantiating the fastMTT tool once (to be used multiple times later)
        self.theFastMTTtool = fastMTTtool()

        # Define the TLorentz vectors as memeber variables to avoid multiple
        # reinitializations within loops
        self.lepFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.eleFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.muFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.tauFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.jetFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.leadingMatch = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.subleadingMatch = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.subsubleadingMatch = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)

        self.pair1FV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.pair2FV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)

        self.higgsTTFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.higgsTTvisFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.higgsBBFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.RadionFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.RadionvisFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)

        self.met = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)

        # Thresholds for different WPs for jets for various years
        # 2016pre:https://btv-wiki.docs.cern.ch/ScaleFactors/UL2016preVFP/
        # 2016post:https://btv-wiki.docs.cern.ch/ScaleFactors/UL2016postVFP/
        # 2017:https://btv-wiki.docs.cern.ch/ScaleFactors/UL2017/
        # 2018:https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
        if self.year == "2016":
            self.LooseJet = 0.0480
            self.MediumJet = 0.2489
            self.TightJet = 0.6377
        elif self.year == "2016APV":
            self.LooseJet = 0.0508
            self.MediumJet = 0.2598
            self.TightJet = 0.6502
        elif self.year == "2017":
            self.LooseJet = 0.0532
            self.MediumJet = 0.3040
            self.TightJet = 0.7476
        elif self.year == "2018":
            self.LooseJet = 0.0490
            self.MediumJet = 0.2783
            self.TightJet = 0.7100

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # self.cutflow_hist = inputFile.cutflow
        self.cutflow_dict = OrderedDict([
            ("Events Generated", int(inputFile.cutflow.GetBinContent(2))),
            ("met_80", int(inputFile.cutflow.GetBinContent(3))),
            ("nFatJet>0", int(inputFile.cutflow.GetBinContent(4))),
            ("METFilters,Triggers,PVCond,met_150", inputTree.GetEntries()),
            ("met_180", 0),
            ("AK8_sel_NoAK8bbTag", 0),
            (" ...breakdown..> Atleast_one_Tau (Reco + ID + cleaning)", 0),
            (" ...breakdown..> Atleast_one_Electron (Reco + IDnoIso + cleaning)", 0),
            (" ...breakdown..> Atleast_one_Muon (Reco + IDnoIso + cleaning)", 0),
            ("Atleast_2leptons_anykind (no Iso-cut for e/mu)", 0),
            ("Atleast_one_pair_anykind (Iso-cut applied for e/mu)", 0),
            (" ...breakdown..> Atleast_one_TauTau_pair", 0),
            (" ...breakdown..> Atleast_one_TauElectron_pair", 0),
            (" ...breakdown..> Atleast_one_TauMuon_pair", 0),
            (" ...breakdown..> TT_channel (max pt pair)", 0),
            (" ...breakdown..> ET_channel (max pt pair)", 0),
            (" ...breakdown..> MT_channel (max pt pair)", 0),
            ("DeltaR_LL<1.5 and abs(Hbb_met_phi) > 1 cut", 0),
            ("Visible Mass HTT > 20 cut", 0),
            ("Medium AK4 btag veto", 0)])
        if ((self.year == "2018") and (self.isMC)):
            self.totalEvents = int(inputTree.GetEntries())
            print(("Total Events in MC file = ", self.totalEvents))
        self.out = wrappedOutputTree
        self.out.branch("eventnominal", "I")
        for sys in self.jesUnc:
            if sys == "":
                self.out.branch("Hbb_lep1_deltaR%s" % (sys), "F")
                self.out.branch("Hbb_lep2_deltaR%s" % (sys), "F")
                self.out.branch("softdropmass%s" % (sys), "F")
                self.out.branch("softdropmassnom%s" % (sys), "F")
                self.out.branch("pnetmass%s" % (sys), "F")
            self.out.branch("channel%s" % (sys), "I")
            self.out.branch("boost%s" % (sys), "I")

            if sys == "":
                self.out.branch("HTT_m%s" % (sys), "F")
                self.out.branch("HTT_eta%s" % (sys), "F")
                self.out.branch("HTT_phi%s" % (sys), "F")
                self.out.branch("HTT_pt%s" % (sys), "F")
                self.out.branch("HTTvis_m%s" % (sys), "F")
                self.out.branch("HTTvis_eta%s" % (sys), "F")
                self.out.branch("HTTvis_phi%s" % (sys), "F")
                self.out.branch("HTTvis_pt%s" % (sys), "F")
                self.out.branch("HTTvis_deltaR%s" % (sys), "F")
                self.out.branch("nallTaus%s" % (sys), "I")
                self.out.branch("allTaus_pt%s" % (sys), "F",
                                lenVar="nallTaus%s" % (sys))
                self.out.branch("allTaus_eta%s" % (sys), "F",
                                lenVar="nallTaus%s" % (sys))
                self.out.branch("allTaus_phi%s" % (sys), "F",
                                lenVar="nallTaus%s" % (sys))
                self.out.branch("allTaus_mass%s" %
                                (sys), "F", lenVar="nallTaus%s" % (sys))
                self.out.branch("allTaus_decayMode%s" %
                                (sys), "F", lenVar="nallTaus%s" % (sys))
                self.out.branch("Xvis_m%s" % (sys), "F")
                self.out.branch("Xvis_eta%s" % (sys), "F")
                self.out.branch("Xvis_phi%s" % (sys), "F")
                self.out.branch("Xvis_pt%s" % (sys), "F")

            self.out.branch("X_m%s" % (sys), "F")
            self.out.branch("X_eta%s" % (sys), "F")
            self.out.branch("X_phi%s" % (sys), "F")
            self.out.branch("X_pt%s" % (sys), "F")

            self.out.branch("ngood_Taus%s" % (sys), "I")
            self.out.branch("ngood_boostedTaus%s" % (sys), "I")
            self.out.branch("ngood_Electrons%s" % (sys), "I")
            self.out.branch("ngood_Muons%s" % (sys), "I")
            self.out.branch("ngood_FatJets%s" % (sys), "I")
            self.out.branch("ngood_Jets%s" % (sys), "I")
            # self.out.branch("ngood_LooseJets%s"%(sys),"I")
            self.out.branch("ngood_MediumJets%s" % (sys), "I")
            self.out.branch("ngood_TightJets%s" % (sys), "I")

            self.out.branch("index_gElectrons%s" % (sys), "I",
                            lenVar="ngood_Electrons%s" % (sys))
            self.out.branch("index_gMuons%s" % (sys), "I",
                            lenVar="ngood_Muons%s" % (sys))
            self.out.branch("index_gTaus%s" % (sys), "I",
                            lenVar="ngood_Taus%s" % (sys))
            self.out.branch("index_gboostedTaus%s" %
                            (sys), "I", lenVar="ngood_boostedTaus%s" % (sys))
            self.out.branch("index_gFatJets%s" % (sys), "I",
                            lenVar="ngood_FatJets%s" % (sys))
            self.out.branch("index_gJets%s" % (sys), "I",
                            lenVar="ngood_Jets%s" % (sys))
            # self.out.branch("index_gLooseJets%s"%(sys),"I",lenVar="ngood_LooseJets%s"%(sys))
            self.out.branch("index_gMediumJets%s" % (sys), "I",
                            lenVar="ngood_MediumJets%s" % (sys))
            self.out.branch("index_gTightJets%s" % (sys), "I",
                            lenVar="ngood_TightJets%s" % (sys))
            self.out.branch("Hbb_met_phi%s" % (sys), "F")

        # keep track of diff wps used
        # self.out.branch("eleID","F")
        # self.out.branch("muID","F")
        # self.out.branch("tauID","F")
        # self.out.branch("btauID","F")

        # Flag for addtional lepton veto
        # self.out.branch("addlepton_vetoflag_all","I")
        # self.out.branch("addlepton_vetoflag_semi","I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        with open(args.tempStore + "/" + "cutflow2_" + self.filename + ".json", "w") as f:
            json.dump(self.cutflow_dict, f, indent=4)
        # pass

    def analyze(self, event):
        # Fill the cutflow_hist for "Total Postproc"
        # self.cutflow2_hist.Fill(0.5)
        # self.cutflow_dict["METFilters,Triggers,PVCond,met_150"] += 1

        # self.out.fillBranch("eleID",self.eleID)
        # self.out.fillBranch("muID",self.muID)
        # self.out.fillBranch("tauID",self.tauID)
        # self.out.fillBranch("btauID",self.btauID)

        def gettaupt(tau, sys):
            if (sys == "tesUp"):
                return tau.pt_tesUp
            elif (sys == "tesDown"):
                return tau.pt_tesDown
            else:
                return tau.pt_nom

        def gettaumass(tau, sys):
            if (sys == "tesUp"):
                return tau.mass_tesUp
            elif (sys == "tesDown"):
                return tau.mass_tesDown
            else:
                return tau.mass_nom

        def getjetpt(jet, sys):
            if ((sys == "") or (sys == "UnclustUp") or (
                    sys == "UnclustDown") or (sys == "tesUp") or (sys == "tesDown")):
                return jet.pt_nom
            elif sys == "jesTotalUp":
                return jet.pt_jesTotalUp
            elif sys == "jesTotalDown":
                return jet.pt_jesTotalDown
            elif sys == "jerUp":
                return jet.pt_jerUp
            elif sys == "jerDown":
                return jet.pt_jerDown
            elif sys == "jesAbsoluteUp":
                return jet.pt_jesAbsoluteUp
            elif sys == "jesAbsoluteDown":
                return jet.pt_jesAbsoluteDown
            elif sys == "jesAbsolute_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesAbsolute_2016Up
                elif (self.year_unc == "2017"):
                    return jet.pt_jesAbsolute_2017Up
                elif (self.year_unc == "2018"):
                    return jet.pt_jesAbsolute_2018Up
            elif sys == "jesAbsolute_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesAbsolute_2016Down
                elif (self.year_unc == "2017"):
                    return jet.pt_jesAbsolute_2017Down
                elif (self.year_unc == "2018"):
                    return jet.pt_jesAbsolute_2018Down
            elif sys == "jesBBEC1Up":
                return jet.pt_jesBBEC1Up
            elif sys == "jesBBEC1Down":
                return jet.pt_jesBBEC1Down
            elif sys == "jesBBEC1_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesBBEC1_2016Up
                elif (self.year_unc == "2017"):
                    return jet.pt_jesBBEC1_2017Up
                elif (self.year_unc == "2018"):
                    return jet.pt_jesBBEC1_2018Up
            elif sys == "jesBBEC1_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesBBEC1_2016Down
                elif (self.year_unc == "2017"):
                    return jet.pt_jesBBEC1_2017Down
                elif (self.year_unc == "2018"):
                    return jet.pt_jesBBEC1_2018Down
            elif sys == "jesEC2Up":
                return jet.pt_jesEC2Up
            elif sys == "jesEC2Down":
                return jet.pt_jesEC2Down
            elif sys == "jesEC2_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesEC2_2016Up
                elif (self.year_unc == "2017"):
                    return jet.pt_jesEC2_2017Up
                elif (self.year_unc == "2018"):
                    return jet.pt_jesEC2_2018Up
            elif sys == "jesEC2_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesEC2_2016Down
                elif (self.year_unc == "2017"):
                    return jet.pt_jesEC2_2017Down
                elif (self.year_unc == "2018"):
                    return jet.pt_jesEC2_2018Down
            elif sys == "jesFlavorQCDUp":
                return jet.pt_jesFlavorQCDUp
            elif sys == "jesFlavorQCDDown":
                return jet.pt_jesFlavorQCDDown
            elif sys == "jesHFUp":
                return jet.pt_jesHFUp
            elif sys == "jesHFDown":
                return jet.pt_jesHFDown
            elif sys == "jesHF_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesHF_2016Up
                elif (self.year_unc == "2017"):
                    return jet.pt_jesHF_2017Up
                elif (self.year_unc == "2018"):
                    return jet.pt_jesHF_2018Up
            elif sys == "jesHF_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesHF_2016Down
                elif (self.year_unc == "2017"):
                    return jet.pt_jesHF_2017Down
                elif (self.year_unc == "2018"):
                    return jet.pt_jesHF_2018Down
            elif sys == "jesRelativeBalUp":
                return jet.pt_jesRelativeBalUp
            elif sys == "jesRelativeBalDown":
                return jet.pt_jesRelativeBalDown
            elif sys == "jesRelativeSample_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesRelativeSample_2016Up
                elif (self.year_unc == "2017"):
                    return jet.pt_jesRelativeSample_2017Up
                elif (self.year_unc == "2018"):
                    return jet.pt_jesRelativeSample_2018Up
            elif sys == "jesRelativeSample_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.pt_jesRelativeSample_2016Down
                elif (self.year_unc == "2017"):
                    return jet.pt_jesRelativeSample_2017Down
                elif (self.year_unc == "2018"):
                    return jet.pt_jesRelativeSample_2018Down

        def getjetmass(jet, sys):
            if ((sys == "") or (sys == "UnclustUp") or (
                    sys == "UnclustDown") or (sys == "tesUp") or (sys == "tesDown")):
                return jet.mass_nom
            elif sys == "jesTotalUp":
                return jet.mass_jesTotalUp
            elif sys == "jesTotalDown":
                return jet.mass_jesTotalDown
            elif sys == "jerUp":
                return jet.mass_jerUp
            elif sys == "jerDown":
                return jet.mass_jerDown
            elif sys == "jesAbsoluteUp":
                return jet.mass_jesAbsoluteUp
            elif sys == "jesAbsoluteDown":
                return jet.mass_jesAbsoluteDown
            elif sys == "jesAbsolute_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesAbsolute_2016Up
                elif (self.year_unc == "2017"):
                    return jet.mass_jesAbsolute_2017Up
                elif (self.year_unc == "2018"):
                    return jet.mass_jesAbsolute_2018Up
            elif sys == "jesAbsolute_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesAbsolute_2016Down
                elif (self.year_unc == "2017"):
                    return jet.mass_jesAbsolute_2017Down
                elif (self.year_unc == "2018"):
                    return jet.mass_jesAbsolute_2018Down
            elif sys == "jesBBEC1Up":
                return jet.mass_jesBBEC1Up
            elif sys == "jesBBEC1Down":
                return jet.mass_jesBBEC1Down
            elif sys == "jesBBEC1_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesBBEC1_2016Up
                elif (self.year_unc == "2017"):
                    return jet.mass_jesBBEC1_2017Up
                elif (self.year_unc == "2018"):
                    return jet.mass_jesBBEC1_2018Up
            elif sys == "jesBBEC1_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesBBEC1_2016Down
                elif (self.year_unc == "2017"):
                    return jet.mass_jesBBEC1_2017Down
                elif (self.year_unc == "2018"):
                    return jet.mass_jesBBEC1_2018Down
            elif sys == "jesEC2Up":
                return jet.mass_jesEC2Up
            elif sys == "jesEC2Down":
                return jet.mass_jesEC2Down
            elif sys == "jesEC2_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesEC2_2016Up
                elif (self.year_unc == "2017"):
                    return jet.mass_jesEC2_2017Up
                elif (self.year_unc == "2018"):
                    return jet.mass_jesEC2_2018Up
            elif sys == "jesEC2_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesEC2_2016Down
                elif (self.year_unc == "2017"):
                    return jet.mass_jesEC2_2017Down
                elif (self.year_unc == "2018"):
                    return jet.mass_jesEC2_2018Down
            elif sys == "jesFlavorQCDUp":
                return jet.mass_jesFlavorQCDUp
            elif sys == "jesFlavorQCDDown":
                return jet.mass_jesFlavorQCDDown
            elif sys == "jesHFUp":
                return jet.mass_jesHFUp
            elif sys == "jesHFDown":
                return jet.mass_jesHFDown
            elif sys == "jesHF_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesHF_2016Up
                elif (self.year_unc == "2017"):
                    return jet.mass_jesHF_2017Up
                elif (self.year_unc == "2018"):
                    return jet.mass_jesHF_2018Up
            elif sys == "jesHF_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesHF_2016Down
                elif (self.year_unc == "2017"):
                    return jet.mass_jesHF_2017Down
                elif (self.year_unc == "2018"):
                    return jet.mass_jesHF_2018Down
            elif sys == "jesRelativeBalUp":
                return jet.mass_jesRelativeBalUp
            elif sys == "jesRelativeBalDown":
                return jet.mass_jesRelativeBalDown
            elif sys == "jesRelativeSample_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesRelativeSample_2016Up
                elif (self.year_unc == "2017"):
                    return jet.mass_jesRelativeSample_2017Up
                elif (self.year_unc == "2018"):
                    return jet.mass_jesRelativeSample_2018Up
            elif sys == "jesRelativeSample_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return jet.mass_jesRelativeSample_2016Down
                elif (self.year_unc == "2017"):
                    return jet.mass_jesRelativeSample_2017Down
                elif (self.year_unc == "2018"):
                    return jet.mass_jesRelativeSample_2018Down

        def getMETpt(sys):
            if ((sys == "") or (sys == "tesUp") or (sys == "tesDown")):
                return event.METcorrected_pt
            elif sys == "jesTotalUp":
                return event.METcorrected_ptScaleUp
            elif sys == "jesTotalDown":
                return event.METcorrected_ptScaleDown
            elif sys == "jerUp":
                return event.METcorrected_ptResUp
            elif sys == "jerDown":
                return event.METcorrected_ptResDown
            elif sys == "jesAbsoluteUp":
                return event.METcorrected_ptScaleAbsoluteUp
            elif sys == "jesAbsoluteDown":
                return event.METcorrected_ptScaleAbsoluteDown
            elif sys == "jesAbsolute_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleAbsolute_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleAbsolute_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleAbsolute_2018Up
            elif sys == "jesAbsolute_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleAbsolute_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleAbsolute_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleAbsolute_2018Down
            elif sys == "jesBBEC1Up":
                return event.METcorrected_ptScaleBBEC1Up
            elif sys == "jesBBEC1Down":
                return event.METcorrected_ptScaleBBEC1Down
            elif sys == "jesBBEC1_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleBBEC1_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleBBEC1_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleBBEC1_2018Up
            elif sys == "jesBBEC1_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleBBEC1_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleBBEC1_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleBBEC1_2018Down
            elif sys == "jesEC2Up":
                return event.METcorrected_ptScaleEC2Up
            elif sys == "jesEC2Down":
                return event.METcorrected_ptScaleEC2Down
            elif sys == "jesEC2_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleEC2_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleEC2_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleEC2_2018Up
            elif sys == "jesEC2_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleEC2_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleEC2_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleEC2_2018Down
            elif sys == "jesFlavorQCDUp":
                return event.METcorrected_ptScaleFlavorQCDUp
            elif sys == "jesFlavorQCDDown":
                return event.METcorrected_ptScaleFlavorQCDDown
            elif sys == "jesHFUp":
                return event.METcorrected_ptScaleHFUp
            elif sys == "jesHFDown":
                return event.METcorrected_ptScaleHFDown
            elif sys == "jesHF_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleHF_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleHF_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleHF_2018Up
            elif sys == "jesHF_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleHF_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleHF_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleHF_2018Down
            elif sys == "jesRelativeBalUp":
                return event.METcorrected_ptScaleRelativeBalUp
            elif sys == "jesRelativeBalDown":
                return event.METcorrected_ptScaleRelativeBalDown
            elif sys == "jesRelativeSample_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleRelativeSample_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleRelativeSample_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleRelativeSample_2018Up
            elif sys == "jesRelativeSample_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_ptScaleRelativeSample_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_ptScaleRelativeSample_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_ptScaleRelativeSample_2018Down
            elif sys == "UnclustUp":
                return event.METcorrected_ptUnclustUp
            elif sys == "UnclustDown":
                return event.METcorrected_ptUnclustDown

        def getMETphi(sys):
            if ((sys == "") or (sys == "UnclustUp") or (
                    sys == "UnclustDown") or (sys == "tesUp") or (sys == "tesDown")):
                return event.METcorrected_phi
            elif sys == "jesTotalUp":
                return event.METcorrected_phiScaleUp
            elif sys == "jesTotalDown":
                return event.METcorrected_phiScaleDown
            elif sys == "jerUp":
                return event.METcorrected_phiResUp
            elif sys == "jerDown":
                return event.METcorrected_phiResDown
            elif sys == "jesAbsoluteUp":
                return event.METcorrected_phiScaleAbsoluteUp
            elif sys == "jesAbsoluteDown":
                return event.METcorrected_phiScaleAbsoluteDown
            elif sys == "jesAbsolute_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleAbsolute_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleAbsolute_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleAbsolute_2018Up
            elif sys == "jesAbsolute_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleAbsolute_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleAbsolute_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleAbsolute_2018Down
            elif sys == "jesBBEC1Up":
                return event.METcorrected_phiScaleBBEC1Up
            elif sys == "jesBBEC1Down":
                return event.METcorrected_phiScaleBBEC1Down
            elif sys == "jesBBEC1_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleBBEC1_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleBBEC1_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleBBEC1_2018Up
            elif sys == "jesBBEC1_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleBBEC1_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleBBEC1_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleBBEC1_2018Down
            elif sys == "jesEC2Up":
                return event.METcorrected_phiScaleEC2Up
            elif sys == "jesEC2Down":
                return event.METcorrected_phiScaleEC2Down
            elif sys == "jesEC2_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleEC2_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleEC2_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleEC2_2018Up
            elif sys == "jesEC2_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleEC2_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleEC2_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleEC2_2018Down
            elif sys == "jesFlavorQCDUp":
                return event.METcorrected_phiScaleFlavorQCDUp
            elif sys == "jesFlavorQCDDown":
                return event.METcorrected_phiScaleFlavorQCDDown
            elif sys == "jesHFUp":
                return event.METcorrected_phiScaleHFUp
            elif sys == "jesHFDown":
                return event.METcorrected_phiScaleHFDown
            elif sys == "jesHF_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleHF_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleHF_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleHF_2018Up
            elif sys == "jesHF_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleHF_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleHF_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleHF_2018Down
            elif sys == "jesRelativeBalUp":
                return event.METcorrected_phiScaleRelativeBalUp
            elif sys == "jesRelativeBalDown":
                return event.METcorrected_phiScaleRelativeBalDown
            elif sys == "jesRelativeSample_%sUp" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleRelativeSample_2016Up
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleRelativeSample_2017Up
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleRelativeSample_2018Up
            elif sys == "jesRelativeSample_%sDown" % (self.year_unc):
                if (self.year_unc == "2016"):
                    return event.METcorrected_phiScaleRelativeSample_2016Down
                elif (self.year_unc == "2017"):
                    return event.METcorrected_phiScaleRelativeSample_2017Down
                elif (self.year_unc == "2018"):
                    return event.METcorrected_phiScaleRelativeSample_2018Down

        def applyPOGselectionToAK4(ak4Object_enu, sys):
            # Based on the Jet MET object contact suggestion - we can remove the pile up ID requirememt.
            # Link: https://cms-pub-talk.web.cern.ch/t/jme-or/29619/4
            # if ((ak4Object_enu[1].pt_nom>30) and
            # (abs(ak4Object_enu[1].eta)<2.5)):
            if ((getjetpt(ak4Object_enu[1], sys) > 30) and (
                    abs(ak4Object_enu[1].eta) < 2.5)):
                if (ak4Object_enu[1].jetId > 1):
                    return True
            return False

        def applyPOGselectionToAK4_for2018Veto(ak4Object_enu):
            # Based on the Jet MET object contact suggestion - we can remove the pile up ID requirememt even for the veto (since it only makes it more stringent)
            # Link: https://cms-pub-talk.web.cern.ch/t/jme-or/29619/4
            if ((ak4Object_enu[1].pt_nom > 15) and (ak4Object_enu[1].eta < -
                                                    1.3) and (ak4Object_enu[1].eta > -
                                                              3.2) and (ak4Object_enu[1].phi < -
                                                                        0.87) and (ak4Object_enu[1].phi > -
               1.57)):
                if (ak4Object_enu[1].jetId > 1):
                    return True
            return False

        def pass_cuts_EleID(electronObject_enu):
            for cutnr in range(0, 10):
                if cutnr == 7:
                    continue
                # if (electronObject_enu[1].vidNestedWPBitmap >> (cutnr*3) &
                # 0x7) < self.eleID:
                if (electronObject_enu[1].vidNestedWPBitmap >> (
                        cutnr * 3) & 0x7) < 2:
                    return False
            return True

        # FatJet and Jet overlap, separation > 1.2 (0.8 + 0.4)
        def JetFatJetOverlap(jetObject_enu, sys):
            # self.jetFV.SetPtEtaPhiM(jetObject_enu[1].pt_nom,jetObject_enu[1].eta,jetObject_enu[1].phi,jetObject_enu[1].mass_nom)
            self.jetFV.SetPtEtaPhiM(
                getjetpt(
                    jetObject_enu[1],
                    sys),
                jetObject_enu[1].eta,
                jetObject_enu[1].phi,
                getjetmass(
                    jetObject_enu[1],
                    sys))
            deltaR = self.jetFV.DeltaR(self.higgsBBFV)
            if deltaR > 1.2:
                return True
            else:
                return False

        # Function used by electron(muon) cleaning vs AK8 > 0.8
        def FatJetConeIsolation(leptonObject_enu):
            self.lepFV.SetPtEtaPhiM(
                leptonObject_enu[1].pt,
                leptonObject_enu[1].eta,
                leptonObject_enu[1].phi,
                leptonObject_enu[1].mass)
            deltaR = self.lepFV.DeltaR(self.higgsBBFV)
            if deltaR > 0.8:
                return True
            else:
                return False

        # Function used by Tau(boostedTau) cleaning vs AK8 > 0.8
        def FatJetTauOverlap(tauObject_enu, boost):
            if boost == 1:
                self.tauFV.SetPtEtaPhiM(
                    tauObject_enu[1].pt,
                    tauObject_enu[1].eta,
                    tauObject_enu[1].phi,
                    tauObject_enu[1].mass)
            else:
                # self.tauFV.SetPtEtaPhiM(tauObject_enu[1].pt_nom,tauObject_enu[1].eta,tauObject_enu[1].phi,tauObject_enu[1].mass_nom)
                self.tauFV.SetPtEtaPhiM(
                    tauObject_enu[1].pt,
                    tauObject_enu[1].eta,
                    tauObject_enu[1].phi,
                    tauObject_enu[1].mass)
            deltaR = self.tauFV.DeltaR(self.higgsBBFV)
            if deltaR > 1.5:
                return True
            else:
                return False

        # Function used by Tau(boostedTau) cleaning vs Electrons > 0.05
        def ElectronTauOverlap(tauObject_enu, elecoll_enu, boost):
            if boost == 1:
                self.tauFV.SetPtEtaPhiM(
                    tauObject_enu[1].pt,
                    tauObject_enu[1].eta,
                    tauObject_enu[1].phi,
                    tauObject_enu[1].mass)
            else:
                # self.tauFV.SetPtEtaPhiM(tauObject_enu[1].pt_nom,tauObject_enu[1].eta,tauObject_enu[1].phi,tauObject_enu[1].mass_nom)
                self.tauFV.SetPtEtaPhiM(
                    tauObject_enu[1].pt,
                    tauObject_enu[1].eta,
                    tauObject_enu[1].phi,
                    tauObject_enu[1].mass)
            for electron in elecoll_enu:
                self.eleFV.SetPtEtaPhiM(
                    electron[1].pt, electron[1].eta, electron[1].phi, 0.0)
                deltaR = self.tauFV.DeltaR(self.eleFV)
                if deltaR <= 0.05:
                    return False
            return True

        # Function used by Tau(boostedTau) cleaning vs Muons > 0.05
        def MuonTauOverlap(tauObject_enu, mucoll_enu, boost):
            if boost == 1:
                self.tauFV.SetPtEtaPhiM(
                    tauObject_enu[1].pt,
                    tauObject_enu[1].eta,
                    tauObject_enu[1].phi,
                    tauObject_enu[1].mass)
            else:
                # self.tauFV.SetPtEtaPhiM(tauObject_enu[1].pt_nom,tauObject_enu[1].eta,tauObject_enu[1].phi,tauObject_enu[1].mass_nom)
                self.tauFV.SetPtEtaPhiM(
                    tauObject_enu[1].pt,
                    tauObject_enu[1].eta,
                    tauObject_enu[1].phi,
                    tauObject_enu[1].mass)
            for muon in mucoll_enu:
                self.muFV.SetPtEtaPhiM(
                    muon[1].pt, muon[1].eta, muon[1].phi, muon[1].mass)
                deltaR = self.tauFV.DeltaR(self.muFV)
                if deltaR <= 0.05:
                    return False
            return True

        def removeOverlapOfAK4WithLightHeavyLeptons(
                ak4Object_enu,
                gTau_index,
                Tau_coll,
                gbTau_index,
                bTau_coll,
                gEle_index,
                Ele_coll,
                gMu_index,
                Mu_coll,
                sys):
            # self.jetFV.SetPtEtaPhiM(ak4Object_enu[1].pt_nom,ak4Object_enu[1].eta,ak4Object_enu[1].phi,ak4Object_enu[1].mass_nom)
            self.jetFV.SetPtEtaPhiM(
                getjetpt(
                    ak4Object_enu[1],
                    sys),
                ak4Object_enu[1].eta,
                ak4Object_enu[1].phi,
                getjetmass(
                    ak4Object_enu[1],
                    sys))

            for index in gEle_index:
                self.eleFV.SetPtEtaPhiM(
                    Ele_coll[index].pt,
                    Ele_coll[index].eta,
                    Ele_coll[index].phi,
                    0.0)
                deltaR = self.jetFV.DeltaR(self.eleFV)
                if deltaR <= 0.4:
                    return False

            for index in gMu_index:
                self.muFV.SetPtEtaPhiM(
                    Mu_coll[index].pt,
                    Mu_coll[index].eta,
                    Mu_coll[index].phi,
                    Mu_coll[index].mass)
                deltaR = self.jetFV.DeltaR(self.muFV)
                if deltaR <= 0.4:
                    return False

            for index in gTau_index:
                self.tauFV.SetPtEtaPhiM(
                    Tau_coll[index].pt,
                    Tau_coll[index].eta,
                    Tau_coll[index].phi,
                    Tau_coll[index].mass)
                deltaR = self.jetFV.DeltaR(self.tauFV)
                if deltaR <= 0.4:
                    return False

            for index in gbTau_index:
                self.tauFV.SetPtEtaPhiM(
                    bTau_coll[index].pt,
                    bTau_coll[index].eta,
                    bTau_coll[index].phi,
                    bTau_coll[index].mass)
                deltaR = self.jetFV.DeltaR(self.tauFV)
                if deltaR <= 0.4:
                    return False

            return True

        def selfPairing(col1, tag):
            combinedPt = -10
            index1 = -1
            index2 = -1
            if len(col1) <= 1:
                return (combinedPt, index1, index2)
            for i in range(len(col1)):
                for j in range(i + 1, len(col1)):
                    self.pair1FV.SetPtEtaPhiM(
                        col1[i][1].pt, col1[i][1].eta, col1[i][1].phi, col1[i][1].mass)
                    self.pair2FV.SetPtEtaPhiM(
                        col1[j][1].pt, col1[j][1].eta, col1[j][1].phi, col1[j][1].mass)

                    deltaR_Pair = self.pair1FV.DeltaR(self.pair2FV)
                    if ((deltaR_Pair <= 0.05)):
                        continue
                    sumFourVector = self.pair1FV + self.pair2FV
                    pt = sumFourVector.Pt()
                    if pt >= combinedPt:
                        combinedPt = pt
                        index1 = col1[i][0]
                        index2 = col1[j][0]
            return (combinedPt, index1, index2)

        def crossPairing(col1, col2, tag):
            # Include a protection for tag string
            tagprotection = ((tag == "be") or (tag == "te")
                             or (tag == "bm") or (tag == "tm"))
            if not tagprotection:
                print("The tags are not correct in cross pairing ##ERROR##------")
                sys.exit()

            combinedPt = -10
            index1 = -1
            index2 = -1
            if (len(col1) == 0 or len(col2) == 0):
                return (combinedPt, index1, index2)
            for i in range(len(col1)):
                for j in range(len(col2)):
                    self.pair1FV.SetPtEtaPhiM(
                        col1[i][1].pt, col1[i][1].eta, col1[i][1].phi, col1[i][1].mass)
                    self.pair2FV.SetPtEtaPhiM(
                        col2[j][1].pt, col2[j][1].eta, col2[j][1].phi, col2[j][1].mass)
                    deltaR_Pair = (self.pair1FV.DeltaR(self.pair2FV))
                    if ((deltaR_Pair <= 0.05)):
                        continue
                    # Found the bug in the code - I was only applying isolation
                    # for boostedTaus but none for HPS
                    if ((tag == "be") or (tag == "te")):
                        if not ElectronIsolationCut(
                                col1[i][1], col2[j][1], tag):
                            continue

                    elif ((tag == "bm") or (tag == "tm")):
                        if not MuonIsolationCut(col1[i][1], col2[j][1], tag):
                            continue

                    sumFourVector = self.pair1FV + self.pair2FV
                    pt = sumFourVector.Pt()
                    if pt >= combinedPt:
                        combinedPt = pt
                        index1 = col1[i][0]
                        index2 = col2[j][0]
            return (combinedPt, index1, index2)

        def MuonIsolationCut(tau, muo, tag):
            tagprotection = ((tag == "bm") or (tag == "tm"))
            if not tagprotection:
                print("The tags are not correct in Electron Isolation ##ERROR##------")
                sys.exit()

            isTau = ""

            self.tauFV.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass)
            self.muFV.SetPtEtaPhiM(muo.pt, muo.eta, muo.phi, muo.mass)
            deltaR = (self.tauFV).DeltaR(self.muFV)
            isolationCut = 0.25

            if (tag == "tm"):
                if ((muo.pfRelIso04_all) < isolationCut):
                    return True
                else:
                    return False

            if deltaR < 0.7:
                isTau = "close"
            if isTau == "close":
                self.leadingMatch.SetPtEtaPhiM(
                    tau.LeadingMuonPt,
                    tau.LeadingMuonEta,
                    tau.LeadingMuonPhi,
                    tau.LeadingMuonM)
                self.subleadingMatch.SetPtEtaPhiM(
                    tau.SubLeadingMuonPt,
                    tau.SubLeadingMuonEta,
                    tau.SubLeadingMuonPhi,
                    tau.SubLeadingMuonM)
                self.subsubleadingMatch.SetPtEtaPhiM(
                    tau.SubSubLeadingMuonPt,
                    tau.SubSubLeadingMuonEta,
                    tau.SubSubLeadingMuonPhi,
                    tau.SubSubLeadingMuonM)
                if (self.muFV.DeltaR(self.leadingMatch) < 0.05):
                    if ((tau.LeadingMuonCorrIso / tau.LeadingMuonPt) < isolationCut):
                        # print ("Leading Muon Matched")
                        return True
                    else:
                        return False
                elif (self.muFV.DeltaR(self.subleadingMatch) < 0.05):
                    if ((tau.SubLeadingMuonCorrIso /
                         tau.SubLeadingMuonPt) < isolationCut):
                        # print ("subLeading Muon Matched")
                        return True
                    else:
                        return False
                elif (self.muFV.DeltaR(self.subsubleadingMatch) < 0.05):
                    if ((tau.SubSubLeadingMuonCorrIso /
                         tau.SubSubLeadingMuonPt) < isolationCut):
                        # print ("subsubLeading Muon Matched")
                        return True
                    else:
                        return False
                else:
                    if ((muo.pfRelIso04_all) < isolationCut):
                        # print ("pf Muon Isolation")
                        return True
                    else:
                        return False

            elif isTau == "":
                if ((muo.pfRelIso04_all) < isolationCut):
                    return True
                else:
                    return False

        # passing the inidivial eletrons from the collection to apply the
        # correction
        def ElectronIsolationCut(tau, ele, tag):
            tagprotection = ((tag == "be") or (tag == "te"))
            if not tagprotection:
                print("The tags are not correct in Electron Isolation ##ERROR##------")
                sys.exit()

            isTau = ""
            isolationCut = 0.0

            if abs(ele.eta) <= 1.479:
                # Barrel values
                # loose =  0.112 + (0.506/ele.pt)
                # medium = 0.0478 + (0.506/ele.pt)
                # tight = 0.0287 + (0.506/ele.pt)
                isolationCut = 0.112 + (0.506 / ele.pt)
            elif (abs(ele.eta) > 1.479) and (abs(ele.eta) <= 2.5):
                # Endcap values
                # loose =  0.108 + (0.963/ele.pt)
                # medium = 0.0658 + (0.963/ele.pt)
                # tight =0.0445 + (0.963/ele.pt)
                isolationCut = 0.108 + (0.963 / ele.pt)
            else:
                return False

            if (tag == "te"):
                if ((ele.pfRelIso03_all) < isolationCut):
                    return True
                else:
                    return False

            self.tauFV.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass)
            self.eleFV.SetPtEtaPhiM(ele.pt, ele.eta, ele.phi, 0.0)
            deltaR = (self.tauFV).DeltaR(self.eleFV)

            if deltaR < 0.6:
                isTau = "close"
            if isTau == "close":
                self.leadingMatch.SetPtEtaPhiM(
                    tau.LeadingElectronPt,
                    tau.LeadingElectronEta,
                    tau.LeadingElectronPhi,
                    0.0)
                self.subleadingMatch.SetPtEtaPhiM(
                    tau.SubLeadingElectronPt,
                    tau.SubLeadingElectronEta,
                    tau.SubLeadingElectronPhi,
                    0.0)
                self.subsubleadingMatch.SetPtEtaPhiM(
                    tau.SubSubLeadingElectronPt,
                    tau.SubSubLeadingElectronEta,
                    tau.SubSubLeadingElectronPhi,
                    0.0)
                if (self.eleFV.DeltaR(self.leadingMatch) < 0.05):
                    if ((tau.LeadingElectronCorrIso /
                         tau.LeadingElectronPt) < isolationCut):
                        # print ("Leading Ele Matched")
                        return True
                    else:
                        return False
                elif (self.eleFV.DeltaR(self.subleadingMatch) < 0.05):
                    if ((tau.SubLeadingElectronCorrIso /
                         tau.SubLeadingElectronPt) < isolationCut):
                        # print ("subLeading Ele Matched")
                        return True
                    else:
                        return False
                elif (self.eleFV.DeltaR(self.subsubleadingMatch) < 0.05):
                    if ((tau.SubSubLeadingElectronCorrIso /
                         tau.SubSubLeadingElectronPt) < isolationCut):
                        # print ("subsubLeading Ele Matched")
                        return True
                    else:
                        return False
                else:
                    if ((ele.pfRelIso03_all) < isolationCut):
                        print("use pf Ele Isolation")
                        return True
                    else:
                        return False

            elif isTau == "":
                if ((ele.pfRelIso03_all) < isolationCut):
                    return True
                else:
                    return False

        def ElectronIsolationCut_addlepveto(eleObject_enu):
            isolationCut = 0.0
            if abs(eleObject_enu[1].eta) <= 1.479:
                # Barrel values
                # loose =  0.112 + (0.506/eleObject_enu[1].pt)
                # medium = 0.0478 + (0.506/eleObject_enu[1].pt)
                # tight = 0.0287 + (0.506/eleObject_enu[1].pt)
                isolationCut = 0.112 + (0.506 / eleObject_enu[1].pt)
            elif (abs(eleObject_enu[1].eta) > 1.479) and (abs(eleObject_enu[1].eta) <= 2.5):
                # Endcap values
                # loose =  0.108 + (0.963/eleObject_enu[1].pt)
                # medium = 0.0658 + (0.963/eleObject_enu[1].pt)
                # tight =0.0445 + (0.963/eleObject_enu[1].pt)
                isolationCut = 0.108 + (0.963 / eleObject_enu[1].pt)
            else:
                return False
            if ((eleObject_enu[1].pfRelIso03_all) < isolationCut):
                return True
            else:
                return False

        def fillBranchesWithDefault(sys):
            if sys == "":
                self.out.fillBranch("Hbb_lep1_deltaR%s" % (sys), -1.00)
                self.out.fillBranch("Hbb_lep2_deltaR%s" % (sys), -1.00)
                self.out.fillBranch("softdropmassnom%s" % (sys), -1.00)
                self.out.fillBranch("softdropmass%s" % (sys), -1.00)
                self.out.fillBranch("pnetmass%s" % (sys), -1.00)

                self.out.fillBranch("HTT_m%s" % (sys), -1.00)
                self.out.fillBranch("HTT_eta%s" % (sys), -1.00)
                self.out.fillBranch("HTT_phi%s" % (sys), -99.99)
                self.out.fillBranch("HTT_pt%s" % (sys), -1.00)
                self.out.fillBranch("HTTvis_m%s" % (sys), -1.00)
                self.out.fillBranch("HTTvis_eta%s" % (sys), -1.00)
                self.out.fillBranch("HTTvis_phi%s" % (sys), -99.99)
                self.out.fillBranch("HTTvis_pt%s" % (sys), -1.00)
                self.out.fillBranch("HTTvis_deltaR%s" % (sys), -1.00)
                self.out.fillBranch("nallTaus%s" % (sys), 0)
                self.out.fillBranch("allTaus_pt%s" % (sys), [])
                self.out.fillBranch("allTaus_eta%s" % (sys), [])
                self.out.fillBranch("allTaus_phi%s" % (sys), [])
                self.out.fillBranch("allTaus_mass%s" % (sys), [])
                self.out.fillBranch("allTaus_decayMode%s" % (sys), [])
                self.out.fillBranch("Xvis_m%s" % (sys), -1.00)
                self.out.fillBranch("Xvis_eta%s" % (sys), -1.00)
                self.out.fillBranch("Xvis_phi%s" % (sys), -99.99)
                self.out.fillBranch("Xvis_pt%s" % (sys), -1.00)

            self.out.fillBranch("channel%s" % (sys), -1)
            self.out.fillBranch("boost%s" % (sys), -1)

            self.out.fillBranch("Hbb_met_phi%s" % (sys), -99.99)
            self.out.fillBranch("X_m%s" % (sys), -1.00)
            self.out.fillBranch("X_eta%s" % (sys), -1.00)
            self.out.fillBranch("X_phi%s" % (sys), -99.99)
            self.out.fillBranch("X_pt%s" % (sys), -1.00)

            self.out.fillBranch("ngood_Taus%s" % (sys), 0)
            self.out.fillBranch("ngood_boostedTaus%s" % (sys), 0)
            self.out.fillBranch("ngood_Electrons%s" % (sys), 0)
            self.out.fillBranch("ngood_Muons%s" % (sys), 0)
            self.out.fillBranch("ngood_FatJets%s" % (sys), 0)
            self.out.fillBranch("ngood_Jets%s" % (sys), 0)
            # self.out.fillBranch("ngood_LooseJets%s"%(sys),0)
            self.out.fillBranch("ngood_MediumJets%s" % (sys), 0)
            self.out.fillBranch("ngood_TightJets%s" % (sys), 0)

            self.out.fillBranch("index_gElectrons%s" % (sys), [])
            self.out.fillBranch("index_gMuons%s" % (sys), [])
            self.out.fillBranch("index_gTaus%s" % (sys), [])
            self.out.fillBranch("index_gboostedTaus%s" % (sys), [])
            self.out.fillBranch("index_gFatJets%s" % (sys), [])
            self.out.fillBranch("index_gJets%s" % (sys), [])
            # self.out.fillBranch("index_gLooseJets%s"%(sys),[])
            self.out.fillBranch("index_gMediumJets%s" % (sys), [])
            self.out.fillBranch("index_gTightJets%s" % (sys), [])

        Jet = Collection(event, 'Jet', 'nJet')
        Electron = Collection(event, 'Electron', 'nElectron')
        FatJet = Collection(event, 'FatJet', 'nFatJet')

        if (self.year == "2018"):
            if ((self.isData) and (event.run >= 319077)):
                Jet_Vetoenu = list(filter(
                    applyPOGselectionToAK4_for2018Veto, enumerate(Jet)))
                # Adding buffer region (in eta-phi) for AK8 Jet veto for HEM:
                # https://mattermost.web.cern.ch/cms-hh-bbtautau/pl/6fk8on61rj8d7jnrwxnwj9uiph
                FatJet_Vetoenu = [
                    x for x in enumerate(FatJet) if (
                        x[1].eta < -
                        1.1) and (
                        x[1].eta > -
                        2.7) and (
                        x[1].phi < -
                        0.67) and (
                        x[1].phi > -
                        1.77) and (
                        x[1].jetId > 1)]
                Electron_Vetoenu = [
                    x for x in enumerate(Electron) if x[1].pt > 10 and (
                        (int(
                            x[1].cutBased)) >= 2) and (
                        x[1].eta < -
                        1.3) and (
                        x[1].eta > -
                        2.5) and (
                        x[1].phi < -
                        0.87) and (
                            x[1].phi > -
                            1.57)]
                if ((len(FatJet_Vetoenu) > 0) or (len(Jet_Vetoenu) > 0)
                        or (len(Electron_Vetoenu) > 0)):
                    return False
            if ((self.isMC)):
                frac = np.random.uniform(0.0, 1.0)
                # (1-0.647724485)
                if (frac > 0.352275515):
                    Jet_Vetoenu = list(filter(
                        applyPOGselectionToAK4_for2018Veto, enumerate(Jet)))
                    # Adding buffer region (in eta-phi) for AK8 Jet veto for
                    # HEM:
                    # https://mattermost.web.cern.ch/cms-hh-bbtautau/pl/6fk8on61rj8d7jnrwxnwj9uiph
                    FatJet_Vetoenu = [
                        x for x in enumerate(FatJet) if (
                            x[1].eta < -
                            1.1) and (
                            x[1].eta > -
                            2.7) and (
                            x[1].phi < -
                            0.67) and (
                            x[1].phi > -
                            1.77) and (
                            x[1].jetId > 1)]
                    Electron_Vetoenu = [
                        x for x in enumerate(Electron) if x[1].pt > 10 and (
                            (int(
                                x[1].cutBased)) >= 2) and (
                            x[1].eta < -
                            1.3) and (
                            x[1].eta > -
                            2.5) and (
                            x[1].phi < -
                            0.87) and (
                            x[1].phi > -
                            1.57)]
                    if ((len(FatJet_Vetoenu) > 0) or (len(Jet_Vetoenu) > 0)
                            or (len(Electron_Vetoenu) > 0)):
                        return False

        if (self.isMC):
            if (self.year_unc == "2018"):
                min_upvar = min([event.METcorrected_ptScaleUp,
                                 event.METcorrected_ptScaleAbsoluteUp,
                                 event.METcorrected_ptScaleAbsolute_2018Up,
                                 event.METcorrected_ptScaleBBEC1Up,
                                 event.METcorrected_ptScaleBBEC1_2018Up,
                                 event.METcorrected_ptScaleEC2Up,
                                 event.METcorrected_ptScaleEC2_2018Up,
                                 event.METcorrected_ptScaleFlavorQCDUp,
                                 event.METcorrected_ptScaleHFUp,
                                 event.METcorrected_ptScaleHF_2018Up,
                                 event.METcorrected_ptScaleRelativeBalUp,
                                 event.METcorrected_ptScaleRelativeSample_2018Up,
                                 event.METcorrected_ptResUp,
                                 event.METcorrected_ptUnclustUp])
                min_downvar = min([event.METcorrected_ptScaleDown,
                                   event.METcorrected_ptScaleAbsoluteDown,
                                   event.METcorrected_ptScaleAbsolute_2018Down,
                                   event.METcorrected_ptScaleBBEC1Down,
                                   event.METcorrected_ptScaleBBEC1_2018Down,
                                   event.METcorrected_ptScaleEC2Down,
                                   event.METcorrected_ptScaleEC2_2018Down,
                                   event.METcorrected_ptScaleFlavorQCDDown,
                                   event.METcorrected_ptScaleHFDown,
                                   event.METcorrected_ptScaleHF_2018Down,
                                   event.METcorrected_ptScaleRelativeBalDown,
                                   event.METcorrected_ptScaleRelativeSample_2018Down,
                                   event.METcorrected_ptResDown,
                                   event.METcorrected_ptUnclustDown])
            elif (self.year_unc == "2017"):
                min_upvar = min([event.METcorrected_ptScaleUp,
                                 event.METcorrected_ptScaleAbsoluteUp,
                                 event.METcorrected_ptScaleAbsolute_2017Up,
                                 event.METcorrected_ptScaleBBEC1Up,
                                 event.METcorrected_ptScaleBBEC1_2017Up,
                                 event.METcorrected_ptScaleEC2Up,
                                 event.METcorrected_ptScaleEC2_2017Up,
                                 event.METcorrected_ptScaleFlavorQCDUp,
                                 event.METcorrected_ptScaleHFUp,
                                 event.METcorrected_ptScaleHF_2017Up,
                                 event.METcorrected_ptScaleRelativeBalUp,
                                 event.METcorrected_ptScaleRelativeSample_2017Up,
                                 event.METcorrected_ptResUp,
                                 event.METcorrected_ptUnclustUp])
                min_downvar = min([event.METcorrected_ptScaleDown,
                                   event.METcorrected_ptScaleAbsoluteDown,
                                   event.METcorrected_ptScaleAbsolute_2017Down,
                                   event.METcorrected_ptScaleBBEC1Down,
                                   event.METcorrected_ptScaleBBEC1_2017Down,
                                   event.METcorrected_ptScaleEC2Down,
                                   event.METcorrected_ptScaleEC2_2017Down,
                                   event.METcorrected_ptScaleFlavorQCDDown,
                                   event.METcorrected_ptScaleHFDown,
                                   event.METcorrected_ptScaleHF_2017Down,
                                   event.METcorrected_ptScaleRelativeBalDown,
                                   event.METcorrected_ptScaleRelativeSample_2017Down,
                                   event.METcorrected_ptResDown,
                                   event.METcorrected_ptUnclustDown])
            elif (self.year_unc == "2016"):
                min_upvar = min([event.METcorrected_ptScaleUp,
                                 event.METcorrected_ptScaleAbsoluteUp,
                                 event.METcorrected_ptScaleAbsolute_2016Up,
                                 event.METcorrected_ptScaleBBEC1Up,
                                 event.METcorrected_ptScaleBBEC1_2016Up,
                                 event.METcorrected_ptScaleEC2Up,
                                 event.METcorrected_ptScaleEC2_2016Up,
                                 event.METcorrected_ptScaleFlavorQCDUp,
                                 event.METcorrected_ptScaleHFUp,
                                 event.METcorrected_ptScaleHF_2016Up,
                                 event.METcorrected_ptScaleRelativeBalUp,
                                 event.METcorrected_ptScaleRelativeSample_2016Up,
                                 event.METcorrected_ptResUp,
                                 event.METcorrected_ptUnclustUp])
                min_downvar = min([event.METcorrected_ptScaleDown,
                                   event.METcorrected_ptScaleAbsoluteDown,
                                   event.METcorrected_ptScaleAbsolute_2016Down,
                                   event.METcorrected_ptScaleBBEC1Down,
                                   event.METcorrected_ptScaleBBEC1_2016Down,
                                   event.METcorrected_ptScaleEC2Down,
                                   event.METcorrected_ptScaleEC2_2016Down,
                                   event.METcorrected_ptScaleFlavorQCDDown,
                                   event.METcorrected_ptScaleHFDown,
                                   event.METcorrected_ptScaleHF_2016Down,
                                   event.METcorrected_ptScaleRelativeBalDown,
                                   event.METcorrected_ptScaleRelativeSample_2016Down,
                                   event.METcorrected_ptResDown,
                                   event.METcorrected_ptUnclustDown])

            # Now reject events (skim further)
            if ((event.METcorrected_pt < 180) and (
                    min_upvar < 180) and (min_downvar < 180)):
                return False
            # Delete the min variables - not used for the rest of the code
            del min_upvar
            del min_downvar
        elif (self.isData):
            if ((event.METcorrected_pt < 180)):
                return False

        # FatJet_skim_enu = filter(lambda x: (x[1].pt_nom >= 180) and (abs(x[1].eta) < 2.5) and (x[1].jetId>1) and (x[1].msoftdrop_nom>=30), enumerate(FatJet))
        FatJet_skim_enu = [x for x in enumerate(FatJet) if (
            x[1].pt_nom >= 180) and (abs(x[1].eta) < 2.5) and (x[1].jetId > 1)]
        if (len(FatJet_skim_enu) == 0):
            return False

        # If the event survives this then no need of this variable
        del FatJet_skim_enu

        passAsingleSystematic = 0

        Muon = Collection(event, 'Muon', 'nMuon')
        Tau = Collection(event, 'Tau', 'nTau')
        boostedTau = Collection(event, 'boostedTau', 'nboostedTau')

        nominal_bool = 0
        for sys in self.jesUnc:
            # Fill the cutflow_hist for "met >= 180"
            # self.cutflow2_hist.Fill(1.5)

            if ((getMETpt(sys) < 180)):
                fillBranchesWithDefault(sys)
                # return False
                continue

            if sys == "":
                self.cutflow_dict["met_180"] += 1

            # FatJet_enu = filter(lambda x: (x[1].pt_nom >= 200) and (abs(x[1].eta) < 2.5) and (x[1].jetId>1) and (x[1].msoftdrop_nom>=30), enumerate(FatJet))
            # FatJet_enu = filter(lambda x: (getjetpt(x[1],sys) >= 200) and (abs(x[1].eta) < 2.5) and (x[1].jetId>1) and (x[1].msoftdrop_nom>=30), enumerate(FatJet))
            FatJet_enu = [x for x in enumerate(FatJet) if (
                getjetpt(x[1], sys) >= 200) and (abs(x[1].eta) < 2.5) and (x[1].jetId > 1)]
            if (len(FatJet_enu) == 0):
                # move on to the next systematics
                fillBranchesWithDefault(sys)
                continue

            # Fill the cutflow_hist for "AK8_sel"
            # self.cutflow2_hist.Fill(2.5)
            if sys == "":
                self.cutflow_dict["AK8_sel_NoAK8bbTag"] += 1

            # HbbScoreList = [(obj_enu[1].particleNetLegacy_Xbb/(obj_enu[1].particleNetLegacy_Xbb + obj_enu[1].particleNetLegacy_QCD)) for obj_enu in FatJet_enu]
            HbbPtList = [(getjetpt(obj_enu[1], sys)) for obj_enu in FatJet_enu]
            # zipPair = zip(FatJet_enu,HbbScoreList)
            zipPair = list(zip(FatJet_enu, HbbPtList))
            FatJet_enu = [
                fatjetobject_enu for fatjetobject_enu,
                _ in sorted(
                    zipPair,
                    key=lambda x: x[1],
                    reverse=True)]

            # remove the HbbPtList, zipPair
            del HbbPtList
            del zipPair

            Jet_enu = [x for x in enumerate(Jet) if applyPOGselectionToAK4(
                x, sys)]

            # Tau_idDeepTau2018v2p5VSe : 1 = VVVLoose, 2 = VVLoose, 3 = VLoose, 4 = Loose, 5 = Medium, 6 = Tight, 7 = VTight, 8 = VVTight* (Recommendation: VVloose, Tight)
            # Tau_idDeepTau2018v2p5VSmu : 1 = VLoose, 2 = Loose, 3 = Medium, 4 = Tight (VLoose)
            # Tau_idDeepTau2018v2p5VSjet : 1 = VVVLoose, 2 = VVLoose, 3 = VLoose, 4 = Loose, 5 = Medium, 6 = Tight, 7 = VTight, 8 = VVTight* (Loose)
            # With TES applied
            # Tau_enu = filter(lambda x: (x[1].pt_nom > 20) and (abs(x[1].eta) < 2.5) and (abs(x[1].dz)<0.2) and (x[1].idDecayModeNewDMs) and (x[1].idDeepTau2018v2p5VSjet >= 4) and (x[1].idDeepTau2018v2p5VSe >= 2) and (x[1].idDeepTau2018v2p5VSmu >= 1), enumerate(Tau)) #The HPS tau id are not bitmps anymore
            # With TAU ES applied
            # Tau_enu = filter(lambda x: (x[1].pt > 20) and (abs(x[1].eta) <
            # 2.5) and (abs(x[1].dz)<0.2) and (x[1].idDecayModeNewDMs) and
            # (x[1].idDeepTau2018v2p5VSjet >= 4) and (x[1].idDeepTau2018v2p5VSe
            # >= 2) and (x[1].idDeepTau2018v2p5VSmu >= 1), enumerate(Tau)) #The
            # HPS tau id are not bitmps anymore
            Tau_enu = [
                x for x in enumerate(Tau) if (
                    gettaupt(
                        x[1], sys) > 20) and (
                    abs(
                        x[1].eta) < 2.5) and (
                    abs(
                        x[1].dz) < 0.2) and (
                    x[1].idDecayModeNewDMs) and (
                    x[1].idDeepTau2018v2p5VSjet >= 4) and (
                    x[1].idDeepTau2018v2p5VSe >= 2) and (
                    x[1].idDeepTau2018v2p5VSmu >= 1)]  # The HPS tau id are not bitmps anymore

            # boostedTau_enu = filter(lambda x: (x[1].pt > 20) and (abs(x[1].eta) < 2.5) and (x[1].rawDeepTau2018v2p7VSjet>=0.85), enumerate(boostedTau))
            boostedTau_enu = [x for x in enumerate(boostedTau) if (gettaupt(x[1], sys) > 20) and (
                abs(x[1].eta) < 2.5) and (x[1].rawDeepTau2018v2p7VSjet >= 0.85)]

            # self.higgsBBFV.SetPtEtaPhiM(FatJet_enu[0][1].pt_nom,FatJet_enu[0][1].eta,FatJet_enu[0][1].phi,FatJet_enu[0][1].mass_nom)
            self.higgsBBFV.SetPtEtaPhiM(
                getjetpt(
                    FatJet_enu[0][1],
                    sys),
                FatJet_enu[0][1].eta,
                FatJet_enu[0][1].phi,
                getjetmass(
                    FatJet_enu[0][1],
                    sys))

            # incorporate ECAL transistion veto : https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaUL2016To2018#General_note_about_ID_SFs
            # https://cms-pub-talk.web.cern.ch/t/egamma-or/30044/2
            Electron_enu = [x for x in enumerate(Electron) if (x[1].pt > 10) and (
                (abs(x[1].eta) <= 1.44) or (abs(x[1].eta) >= 1.57))]
            Electron_enu = list(filter(pass_cuts_EleID, Electron_enu))

            ###################################################################
            # Code for addtional lepton veto
            # Electron_addlep_enu = filter(lambda x: (x[1].pt > 10) and ((abs(x[1].eta) <= 1.44) or (abs(x[1].eta) >= 1.57)) and (x[1].cutBased>=2), enumerate(Electron))
            # Electron_addlep_enu = filter(ElectronIsolationCut_addlepveto,Electron_addlep_enu)
            # Electron_addlep_enu = filter(FatJetConeIsolation,Electron_addlep_enu)

            # Muon_addlep_enu = filter(lambda x: x[1].pt > 10 and (abs(x[1].eta) < 2.5) and x[1].looseId and (x[1].pfRelIso04_all<0.25), enumerate(Muon))
            # Muon_addlep_enu = filter(FatJetConeIsolation,Muon_addlep_enu)

            ###################################################################

            # if self.muID == 2:
            Muon_enu = [x for x in enumerate(Muon) if x[1].pt > 15 and (
                abs(x[1].eta) < 2.4) and x[1].looseId]
            # elif self.muID == 3:
            # Muon_enu = filter(lambda x: x[1].pt > 10 and (abs(x[1].eta) < 2.5) and x[1].mediumId, enumerate(Muon))
            # elif self.muID == 4:
            # Muon_enu = filter(lambda x: x[1].pt > 10 and (abs(x[1].eta) < 2.5) and x[1].tightId, enumerate(Muon))

            # Object cleaning procedures
            Jet_enu = [x for x in Jet_enu if JetFatJetOverlap(x, sys)]

            Electron_enu = list(filter(FatJetConeIsolation, Electron_enu))
            Muon_enu = list(filter(FatJetConeIsolation, Muon_enu))

            Tau_enu = [x for x in Tau_enu if FatJetTauOverlap(x, boost=0)]
            boostedTau_enu = [
                x for x in boostedTau_enu if FatJetTauOverlap(
                    x, boost=1)]

            Tau_enu = [x for x in Tau_enu if ElectronTauOverlap(
                x, Electron_enu, boost=0)]
            boostedTau_enu = [x for x in boostedTau_enu if ElectronTauOverlap(
                x, Electron_enu, boost=1)]

            Tau_enu = [x for x in Tau_enu if MuonTauOverlap(
                x, Muon_enu, boost=0)]
            boostedTau_enu = [x for x in boostedTau_enu if MuonTauOverlap(
                x, Muon_enu, boost=1)]

            # Cutflow counting
            if sys == "":
                if (len(Tau_enu) > 0 or len(boostedTau_enu) > 0):
                    self.cutflow_dict[" ...breakdown..> Atleast_one_Tau (Reco + ID + cleaning)"] += 1
                if (len(Muon_enu) > 0):
                    self.cutflow_dict[" ...breakdown..> Atleast_one_Muon (Reco + IDnoIso + cleaning)"] += 1
                if (len(Electron_enu) > 0):
                    self.cutflow_dict[
                        " ...breakdown..> Atleast_one_Electron (Reco + IDnoIso + cleaning)"] += 1

            enoughleptonstopair = (
                ((len(Tau_enu) +
                  len(Muon_enu) +
                    len(Electron_enu)) >= 2) or (
                    (len(boostedTau_enu) +
                     len(Muon_enu) +
                        len(Electron_enu)) >= 2))

            if not enoughleptonstopair:
                # move to the next systematics
                fillBranchesWithDefault(sys)
                continue

            if (sys == ""):
                if (((len(Tau_enu) + len(Muon_enu) + len(Electron_enu)) >= 2)
                        or ((len(boostedTau_enu) + len(Muon_enu) + len(Electron_enu)) >= 2)):
                    self.cutflow_dict["Atleast_2leptons_anykind (no Iso-cut for e/mu)"] += 1

            list = {}
            list["bb"] = selfPairing(boostedTau_enu, "bb")
            list["tt"] = selfPairing(Tau_enu, "tt")
            list["be"] = crossPairing(boostedTau_enu, Electron_enu, "be")
            list["bm"] = crossPairing(boostedTau_enu, Muon_enu, "bm")
            list["te"] = crossPairing(Tau_enu, Electron_enu, "te")
            list["tm"] = crossPairing(Tau_enu, Muon_enu, "tm")

            Keymax = max(list, key=lambda x: list[x][0])

            if (list[Keymax][0] < 0):
                # return False
                # Electron and Muon isolation failed - move on to the next
                # systematics
                fillBranchesWithDefault(sys)
                continue

            if (sys == ""):
                self.cutflow_dict["Atleast_one_pair_anykind (Iso-cut applied for e/mu)"] += 1
                if ((list["bb"][0] > 0) or (list["tt"][0] > 0)):
                    self.cutflow_dict[" ...breakdown..> Atleast_one_TauTau_pair"] += 1
                if ((list["be"][0] > 0) or (list["te"][0] > 0)):
                    self.cutflow_dict[" ...breakdown..> Atleast_one_TauElectron_pair"] += 1
                if ((list["bm"][0] > 0) or (list["tm"][0] > 0)):
                    self.cutflow_dict[" ...breakdown..> Atleast_one_TauMuon_pair"] += 1

            gFatJet_index = [FatJet_enu[0][0]]
            gTau_index = []
            gboostedTau_index = []
            gElectron_index = []
            gMuon_index = []
            gJet_index = []
            firstLepton = fastMTTlepton()
            secondLepton = fastMTTlepton()
            theMET = fastMTTmet(
                # measuredX = event.METcorrected_pt * math.cos(event.METcorrected_phi),
                # measuredY = event.METcorrected_pt * math.sin(event.METcorrected_phi),
                measuredX=getMETpt(sys) * math.cos(getMETphi(sys)),
                measuredY=getMETpt(sys) * math.sin(getMETphi(sys)),
                xx=event.MET_covXX,
                xy=event.MET_covXY,
                yy=event.MET_covYY)

            if Keymax == "bb":
                self.out.fillBranch("channel%s" % (sys), 0)
                self.out.fillBranch("boost%s" % (sys), 1)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 2)
                gboostedTau_index = [list[Keymax][1], list[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (sys),
                                        [boostedTau[gboostedTau_index[0]].decayMode,
                                         boostedTau[gboostedTau_index[1]].decayMode])
                # firstLepton  =  fastMTTlepton(pt = boostedTau[gboostedTau_index[0]].pt,eta = boostedTau[gboostedTau_index[0]].eta,phi = boostedTau[gboostedTau_index[0]].phi,m = boostedTau[gboostedTau_index[0]].mass,leptonType = 'Tau',tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
                firstLepton = fastMTTlepton(pt=gettaupt(boostedTau[gboostedTau_index[0]],
                                                        sys),
                                            eta=boostedTau[gboostedTau_index[0]].eta,
                                            phi=boostedTau[gboostedTau_index[0]].phi,
                                            m=gettaumass(boostedTau[gboostedTau_index[0]],
                                                         sys),
                                            leptonType='Tau',
                                            tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
                self.pair1FV.SetPtEtaPhiM(gettaupt(boostedTau[gboostedTau_index[0]],
                                                   sys),
                                          boostedTau[gboostedTau_index[0]].eta,
                                          boostedTau[gboostedTau_index[0]].phi,
                                          gettaumass(boostedTau[gboostedTau_index[0]],
                                                     sys))
                secondLepton = fastMTTlepton(pt=gettaupt(boostedTau[gboostedTau_index[1]],
                                                         sys),
                                             eta=boostedTau[gboostedTau_index[1]].eta,
                                             phi=boostedTau[gboostedTau_index[1]].phi,
                                             m=gettaumass(boostedTau[gboostedTau_index[1]],
                                                          sys),
                                             leptonType='Tau',
                                             tauDecayMode=boostedTau[gboostedTau_index[1]].decayMode)
                self.pair2FV.SetPtEtaPhiM(gettaupt(boostedTau[gboostedTau_index[1]],
                                                   sys),
                                          boostedTau[gboostedTau_index[1]].eta,
                                          boostedTau[gboostedTau_index[1]].phi,
                                          gettaumass(boostedTau[gboostedTau_index[1]],
                                                     sys))
                # Fill the cutflow_hist for "FH_channel"
                # self.cutflow2_hist.Fill(4.5)
                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> TT_channel (max pt pair)"] += 1
            elif Keymax == "tt":
                self.out.fillBranch("channel%s" % (sys), 0)
                self.out.fillBranch("boost%s" % (sys), 0)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 2)
                gTau_index = [list[Keymax][1], list[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (
                        sys), [Tau[gTau_index[0]].decayMode, Tau[gTau_index[0]].decayMode])
                # firstLepton  =  fastMTTlepton(pt = Tau[gTau_index[0]].pt,eta = Tau[gTau_index[0]].eta,phi = Tau[gTau_index[0]].phi,m = Tau[gTau_index[0]].mass,leptonType = 'Tau',tauDecayMode=Tau[gTau_index[0]].decayMode)
                firstLepton = fastMTTlepton(pt=gettaupt(Tau[gTau_index[0]],
                                                        sys),
                                            eta=Tau[gTau_index[0]].eta,
                                            phi=Tau[gTau_index[0]].phi,
                                            m=gettaumass(Tau[gTau_index[0]],
                                                         sys),
                                            leptonType='Tau',
                                            tauDecayMode=Tau[gTau_index[0]].decayMode)
                self.pair1FV.SetPtEtaPhiM(gettaupt(
                    Tau[gTau_index[0]], sys), Tau[gTau_index[0]].eta, Tau[gTau_index[0]].phi, gettaumass(Tau[gTau_index[0]], sys))
                secondLepton = fastMTTlepton(pt=gettaupt(Tau[gTau_index[1]],
                                                         sys),
                                             eta=Tau[gTau_index[1]].eta,
                                             phi=Tau[gTau_index[1]].phi,
                                             m=gettaumass(Tau[gTau_index[1]],
                                                          sys),
                                             leptonType='Tau',
                                             tauDecayMode=Tau[gTau_index[1]].decayMode)
                self.pair2FV.SetPtEtaPhiM(gettaupt(
                    Tau[gTau_index[1]], sys), Tau[gTau_index[1]].eta, Tau[gTau_index[1]].phi, gettaumass(Tau[gTau_index[0]], sys))
                # firstLepton  =  fastMTTlepton(pt = Tau[gTau_index[0]].pt_nom,eta = Tau[gTau_index[0]].eta,phi = Tau[gTau_index[0]].phi,m = Tau[gTau_index[0]].mass_nom,leptonType = 'Tau',tauDecayMode=Tau[gTau_index[0]].decayMode)
                # self.pair1FV.SetPtEtaPhiM(Tau[gTau_index[0]].pt_nom, Tau[gTau_index[0]].eta,Tau[gTau_index[0]].phi,Tau[gTau_index[0]].mass_nom)
                # secondLepton =  fastMTTlepton(pt = Tau[gTau_index[1]].pt_nom,eta = Tau[gTau_index[1]].eta,phi = Tau[gTau_index[1]].phi,m = Tau[gTau_index[1]].mass_nom,leptonType = 'Tau',tauDecayMode=Tau[gTau_index[1]].decayMode)
                # self.pair2FV.SetPtEtaPhiM(Tau[gTau_index[1]].pt_nom, Tau[gTau_index[1]].eta,Tau[gTau_index[1]].phi,Tau[gTau_index[1]].mass_nom)
                # Fill the cutflow_hist for "FH_channel"
                # self.cutflow2_hist.Fill(4.5)
                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> TT_channel (max pt pair)"] += 1
            elif Keymax == "be":
                self.out.fillBranch("channel%s" % (sys), 1)
                self.out.fillBranch("boost%s" % (sys), 1)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 1)
                gboostedTau_index = [list[Keymax][1]]
                gElectron_index = [list[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (
                        sys), [boostedTau[gboostedTau_index[0]].decayMode])
                # firstLepton  =  fastMTTlepton(pt = boostedTau[gboostedTau_index[0]].pt,eta = boostedTau[gboostedTau_index[0]].eta,phi = boostedTau[gboostedTau_index[0]].phi,m = boostedTau[gboostedTau_index[0]].mass,leptonType = 'Tau',tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
                firstLepton = fastMTTlepton(pt=gettaupt(boostedTau[gboostedTau_index[0]],
                                                        sys),
                                            eta=boostedTau[gboostedTau_index[0]].eta,
                                            phi=boostedTau[gboostedTau_index[0]].phi,
                                            m=gettaumass(boostedTau[gboostedTau_index[0]],
                                                         sys),
                                            leptonType='Tau',
                                            tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
                self.pair1FV.SetPtEtaPhiM(gettaupt(boostedTau[gboostedTau_index[0]],
                                                   sys),
                                          boostedTau[gboostedTau_index[0]].eta,
                                          boostedTau[gboostedTau_index[0]].phi,
                                          gettaumass(boostedTau[gboostedTau_index[0]],
                                                     sys))
                secondLepton = fastMTTlepton(pt=Electron[gElectron_index[0]].pt,
                                             eta=Electron[gElectron_index[0]].eta,
                                             phi=Electron[gElectron_index[0]].phi,
                                             m=0.51100e-3,
                                             leptonType='Electron',
                                             tauDecayMode=-1)
                self.pair2FV.SetPtEtaPhiM(
                    Electron[gElectron_index[0]].pt, Electron[gElectron_index[0]].eta, Electron[gElectron_index[0]].phi, 0.0)
                # Fill the cutflow_hist for "SL_channel"
                # self.cutflow2_hist.Fill(5.5)
                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> ET_channel (max pt pair)"] += 1
            elif Keymax == "te":
                self.out.fillBranch("channel%s" % (sys), 1)
                self.out.fillBranch("boost%s" % (sys), 0)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 1)
                gTau_index = [list[Keymax][1]]
                gElectron_index = [list[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (sys), [
                                        Tau[gTau_index[0]].decayMode])
                firstLepton = fastMTTlepton(pt=gettaupt(Tau[gTau_index[0]],
                                                        sys),
                                            eta=Tau[gTau_index[0]].eta,
                                            phi=Tau[gTau_index[0]].phi,
                                            m=gettaumass(Tau[gTau_index[0]],
                                                         sys),
                                            leptonType='Tau',
                                            tauDecayMode=Tau[gTau_index[0]].decayMode)
                self.pair1FV.SetPtEtaPhiM(gettaupt(
                    Tau[gTau_index[0]], sys), Tau[gTau_index[0]].eta, Tau[gTau_index[0]].phi, gettaumass(Tau[gTau_index[0]], sys))
                # firstLepton  =  fastMTTlepton(pt = Tau[gTau_index[0]].pt_nom,eta = Tau[gTau_index[0]].eta,phi = Tau[gTau_index[0]].phi,m = Tau[gTau_index[0]].mass_nom,leptonType = 'Tau',tauDecayMode=Tau[gTau_index[0]].decayMode)
                # self.pair1FV.SetPtEtaPhiM(Tau[gTau_index[0]].pt_nom, Tau[gTau_index[0]].eta,Tau[gTau_index[0]].phi,Tau[gTau_index[0]].mass_nom)
                secondLepton = fastMTTlepton(pt=Electron[gElectron_index[0]].pt,
                                             eta=Electron[gElectron_index[0]].eta,
                                             phi=Electron[gElectron_index[0]].phi,
                                             m=0.51100e-3,
                                             leptonType='Electron',
                                             tauDecayMode=-1)
                self.pair2FV.SetPtEtaPhiM(
                    Electron[gElectron_index[0]].pt, Electron[gElectron_index[0]].eta, Electron[gElectron_index[0]].phi, 0.0)
                # Fill the cutflow_hist for "SL_channel"
                # self.cutflow2_hist.Fill(5.5)
                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> ET_channel (max pt pair)"] += 1
            elif Keymax == "bm":
                self.out.fillBranch("channel%s" % (sys), 2)
                self.out.fillBranch("boost%s" % (sys), 1)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 1)
                gboostedTau_index = [list[Keymax][1]]
                gMuon_index = [list[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (
                        sys), [boostedTau[gboostedTau_index[0]].decayMode])
                firstLepton = fastMTTlepton(pt=gettaupt(boostedTau[gboostedTau_index[0]],
                                                        sys),
                                            eta=boostedTau[gboostedTau_index[0]].eta,
                                            phi=boostedTau[gboostedTau_index[0]].phi,
                                            m=gettaumass(boostedTau[gboostedTau_index[0]],
                                                         sys),
                                            leptonType='Tau',
                                            tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
                self.pair1FV.SetPtEtaPhiM(gettaupt(boostedTau[gboostedTau_index[0]],
                                                   sys),
                                          boostedTau[gboostedTau_index[0]].eta,
                                          boostedTau[gboostedTau_index[0]].phi,
                                          gettaumass(boostedTau[gboostedTau_index[0]],
                                                     sys))
                secondLepton = fastMTTlepton(pt=Muon[gMuon_index[0]].pt,
                                             eta=Muon[gMuon_index[0]].eta,
                                             phi=Muon[gMuon_index[0]].phi,
                                             m=Muon[gMuon_index[0]].mass,
                                             leptonType='Muon',
                                             tauDecayMode=-1)
                self.pair2FV.SetPtEtaPhiM(Muon[gMuon_index[0]].pt,
                                          Muon[gMuon_index[0]].eta,
                                          Muon[gMuon_index[0]].phi,
                                          Muon[gMuon_index[0]].mass)
                # Fill the cutflow_hist for "SL_channel"
                # self.cutflow2_hist.Fill(5.5)
                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> MT_channel (max pt pair)"] += 1

            elif Keymax == "tm":
                self.out.fillBranch("channel%s" % (sys), 2)
                self.out.fillBranch("boost%s" % (sys), 0)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 1)
                gTau_index = [list[Keymax][1]]
                gMuon_index = [list[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (sys), [
                                        Tau[gTau_index[0]].decayMode])
                firstLepton = fastMTTlepton(pt=gettaupt(Tau[gTau_index[0]],
                                                        sys),
                                            eta=Tau[gTau_index[0]].eta,
                                            phi=Tau[gTau_index[0]].phi,
                                            m=gettaumass(Tau[gTau_index[0]],
                                                         sys),
                                            leptonType='Tau',
                                            tauDecayMode=Tau[gTau_index[0]].decayMode)
                self.pair1FV.SetPtEtaPhiM(gettaupt(
                    Tau[gTau_index[0]], sys), Tau[gTau_index[0]].eta, Tau[gTau_index[0]].phi, gettaumass(Tau[gTau_index[0]], sys))
                # firstLepton  =  fastMTTlepton(pt = Tau[gTau_index[0]].pt_nom,eta = Tau[gTau_index[0]].eta,phi = Tau[gTau_index[0]].phi,m = Tau[gTau_index[0]].mass_nom,leptonType = 'Tau',tauDecayMode=Tau[gTau_index[0]].decayMode)
                # self.pair1FV.SetPtEtaPhiM(Tau[gTau_index[0]].pt_nom, Tau[gTau_index[0]].eta,Tau[gTau_index[0]].phi,Tau[gTau_index[0]].mass_nom)
                secondLepton = fastMTTlepton(pt=Muon[gMuon_index[0]].pt,
                                             eta=Muon[gMuon_index[0]].eta,
                                             phi=Muon[gMuon_index[0]].phi,
                                             m=Muon[gMuon_index[0]].mass,
                                             leptonType='Muon',
                                             tauDecayMode=-1)
                self.pair2FV.SetPtEtaPhiM(Muon[gMuon_index[0]].pt,
                                          Muon[gMuon_index[0]].eta,
                                          Muon[gMuon_index[0]].phi,
                                          Muon[gMuon_index[0]].mass)
                # Fill the cutflow_hist for "SL_channel"
                # self.cutflow2_hist.Fill(5.5)
                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> MT_channel (max pt pair)"] += 1

            self.met.SetPtEtaPhiM(getMETpt(sys), 0.0, getMETphi(sys), 0.0)
            pass_quality_delcuts = (
                ((abs(
                    self.higgsBBFV.DeltaPhi(
                        self.met))) > 1) and (
                    (self.pair1FV.DeltaR(
                        self.pair2FV)) > 0) and (
                    (self.pair1FV.DeltaR(
                        self.pair2FV)) < 1.5))

            if not pass_quality_delcuts:
                # move on to the next systematic
                fillBranchesWithDefault(sys)
                continue
            # delete the temp variable
            del pass_quality_delcuts

            # Fill cut flow
            if (sys == ""):
                self.cutflow_dict["DeltaR_LL<1.5 and abs(Hbb_met_phi) > 1 cut"] += 1

            self.out.fillBranch("Hbb_met_phi%s" %
                                (sys), self.higgsBBFV.DeltaPhi(self.met))

            self.theFastMTTtool.setFirstLepton(firstLepton)
            self.theFastMTTtool.setSecondLepton(secondLepton)
            self.theFastMTTtool.setTheMET(theMET)
            higgsFV_list = self.theFastMTTtool.getFastMTTfourvector()
            self.higgsTTFV.SetPtEtaPhiM(
                higgsFV_list[0],
                higgsFV_list[1],
                higgsFV_list[2],
                higgsFV_list[3])
            self.RadionFV = self.higgsTTFV + self.higgsBBFV

            self.higgsTTvisFV = self.pair1FV + self.pair2FV
            self.RadionvisFV = self.higgsTTvisFV + self.higgsBBFV

            if (self.higgsTTvisFV.M() <= 20):
                fillBranchesWithDefault(sys)
                # move to next systematic
                continue

            # Fill cut flow
            if (sys == ""):
                self.cutflow_dict["Visible Mass HTT > 20 cut"] += 1

            Jet_enu = [
                x for x in Jet_enu if removeOverlapOfAK4WithLightHeavyLeptons(
                    x,
                    gTau_index,
                    Tau,
                    gboostedTau_index,
                    boostedTau,
                    gElectron_index,
                    Electron,
                    gMuon_index,
                    Muon,
                    sys)]
            gJet_index = [x[0] for x in Jet_enu]
            Jet_enu_Loose = [
                x for x in Jet_enu if x[1].btagDeepFlavB >= self.LooseJet]
            gJet_Looseindex = [x[0] for x in Jet_enu_Loose]
            Jet_enu_Medium = [
                x for x in Jet_enu if x[1].btagDeepFlavB >= self.MediumJet]
            gJet_Mediumindex = [x[0] for x in Jet_enu_Medium]
            Jet_enu_Tight = [
                x for x in Jet_enu if x[1].btagDeepFlavB >= self.TightJet]
            gJet_Tightindex = [x[0] for x in Jet_enu_Tight]

            # compute the fastMTT vector
            # theFastMTTtool = fastMTTtool()

            if sys == "":
                self.out.fillBranch("HTT_m%s" % (sys), self.higgsTTFV.M())
                self.out.fillBranch("HTT_eta%s" % (sys), self.higgsTTFV.Eta())
                self.out.fillBranch("HTT_phi%s" % (sys), self.higgsTTFV.Phi())
                self.out.fillBranch("HTT_pt%s" % (sys), self.higgsTTFV.Pt())
                self.out.fillBranch("HTTvis_m%s" %
                                    (sys), self.higgsTTvisFV.M())
                self.out.fillBranch("HTTvis_eta%s" %
                                    (sys), self.higgsTTvisFV.Eta())
                self.out.fillBranch("HTTvis_phi%s" %
                                    (sys), self.higgsTTvisFV.Phi())
                self.out.fillBranch("HTTvis_pt%s" %
                                    (sys), self.higgsTTvisFV.Pt())
                self.out.fillBranch("HTTvis_deltaR%s" %
                                    (sys), self.pair1FV.DeltaR(self.pair2FV))

            if sys == "":
                self.out.fillBranch("Hbb_lep1_deltaR%s" %
                                    (sys), self.higgsBBFV.DeltaR(self.pair1FV))
                self.out.fillBranch("Hbb_lep2_deltaR%s" %
                                    (sys), self.higgsBBFV.DeltaR(self.pair2FV))
                self.out.fillBranch("softdropmassnom%s" %
                                    (sys), FatJet_enu[0][1].msoftdrop_nom)
                self.out.fillBranch("softdropmass%s" %
                                    (sys), FatJet_enu[0][1].msoftdrop)
                self.out.fillBranch("pnetmass%s" % (
                    sys), FatJet_enu[0][1].particleNetLegacy_mass)

            self.out.fillBranch("X_m%s" % (sys), self.RadionFV.M())
            self.out.fillBranch("X_eta%s" % (sys), self.RadionFV.Eta())
            self.out.fillBranch("X_phi%s" % (sys), self.RadionFV.Phi())
            self.out.fillBranch("X_pt%s" % (sys), self.RadionFV.Pt())
            if sys == "":
                self.out.fillBranch("Xvis_m%s" % (sys), self.RadionvisFV.M())
                self.out.fillBranch("Xvis_eta%s" %
                                    (sys), self.RadionvisFV.Eta())
                self.out.fillBranch("Xvis_phi%s" %
                                    (sys), self.RadionvisFV.Phi())
                self.out.fillBranch("Xvis_pt%s" % (sys), self.RadionvisFV.Pt())

            # Fill cut flow
            # _##if (sys==""):
            # _##	if ((self.pair1FV.DeltaR(self.pair2FV)>0) and (self.pair1FV.DeltaR(self.pair2FV)<1.5) and (self.RadionFV.M()>=750) and (self.RadionFV.M()<=5010)):
            # _##		self.cutflow_dict["750<=X_m<=5010 cut"] += 1

            self.out.fillBranch("ngood_boostedTaus%s" %
                                (sys), len(gboostedTau_index))
            self.out.fillBranch("ngood_Taus%s" % (sys), len(gTau_index))
            self.out.fillBranch("ngood_Electrons%s" %
                                (sys), len(gElectron_index))
            self.out.fillBranch("ngood_Muons%s" % (sys), len(gMuon_index))
            self.out.fillBranch("ngood_FatJets%s" % (sys), len(gFatJet_index))
            self.out.fillBranch("ngood_Jets%s" % (sys), len(gJet_index))
            # self.out.fillBranch("ngood_LooseJets%s"%(sys),len(gJet_Looseindex))
            self.out.fillBranch("ngood_MediumJets%s" %
                                (sys), len(gJet_Mediumindex))
            self.out.fillBranch("ngood_TightJets%s" %
                                (sys), len(gJet_Tightindex))

            # Fill cut flow
            if (sys == ""):
                if ((len(gJet_Mediumindex) == 0)):
                    self.cutflow_dict["Medium AK4 btag veto"] += 1

            self.out.fillBranch("index_gboostedTaus%s" %
                                (sys), gboostedTau_index)
            self.out.fillBranch("index_gTaus%s" % (sys), gTau_index)
            self.out.fillBranch("index_gElectrons%s" % (sys), gElectron_index)
            self.out.fillBranch("index_gMuons%s" % (sys), gMuon_index)
            self.out.fillBranch("index_gFatJets%s" % (sys), gFatJet_index)
            self.out.fillBranch("index_gJets%s" % (sys), gJet_index)
            # self.out.fillBranch("index_gLooseJets%s"%(sys),gJet_Looseindex)
            self.out.fillBranch("index_gMediumJets%s" %
                                (sys), gJet_Mediumindex)
            self.out.fillBranch("index_gTightJets%s" % (sys), gJet_Tightindex)

            if ((Keymax == "bb") or (Keymax == "tt")):
                if (sys == ""):
                    self.out.fillBranch("allTaus_pt%s" % (
                        sys), [self.pair1FV.Pt(), self.pair2FV.Pt()])
                    self.out.fillBranch("allTaus_eta%s" % (
                        sys), [self.pair1FV.Eta(), self.pair2FV.Eta()])
                    self.out.fillBranch("allTaus_phi%s" % (
                        sys), [self.pair1FV.Phi(), self.pair2FV.Phi()])
                    self.out.fillBranch("allTaus_mass%s" % (
                        sys), [self.pair1FV.M(), self.pair2FV.M()])
                # if ((len(Electron_addlep_enu)!=0) or (len(Muon_addlep_enu)!=0)):
                # self.out.fillBranch("addlepton_vetoflag_all",1)
                # self.out.fillBranch("addlepton_vetoflag_semi",0)
                # elif ((len(Electron_addlep_enu)==0) and (len(Muon_addlep_enu)==0)):
                # self.out.fillBranch("addlepton_vetoflag_all",0)
                # self.out.fillBranch("addlepton_vetoflag_semi",0)
            else:
                if (sys == ""):
                    self.out.fillBranch("allTaus_pt%s" %
                                        (sys), [self.pair1FV.Pt()])
                    self.out.fillBranch("allTaus_eta%s" %
                                        (sys), [self.pair1FV.Eta()])
                    self.out.fillBranch("allTaus_phi%s" %
                                        (sys), [self.pair1FV.Phi()])
                    self.out.fillBranch("allTaus_mass%s" %
                                        (sys), [self.pair1FV.M()])
                # if ((len(Electron_addlep_enu)!=0) or (len(Muon_addlep_enu)!=0)):
                # self.out.fillBranch("addlepton_vetoflag_all",1)
                # self.out.fillBranch("addlepton_vetoflag_semi",1)
                # elif ((len(Electron_addlep_enu)==0) and (len(Muon_addlep_enu)==0)):
                # self.out.fillBranch("addlepton_vetoflag_all",0)
                # self.out.fillBranch("addlepton_vetoflag_semi",0)
            if (sys == ""):
                nominal_bool = 1
            passAsingleSystematic = passAsingleSystematic + 1

        if passAsingleSystematic > 0:
            if (nominal_bool == 1):
                self.out.fillBranch("eventnominal", 1)
            else:
                self.out.fillBranch("eventnominal", 0)
            return True
        else:
            return False


def call_postpoc(files):
    nameStrip = files.strip()
    filename = (nameStrip.split('/')[-1]).split('.')[-2]
    print(filename)

    if (search("Run", filename)):
        print(("This is a ", args.year, " Data file = ", filename))

        def tesModule(): return TauEnergyScaleForHPSandBoosted(
            args.year, True, tauID_wp='Loose', ele_wp='VVLoose')
        # mainModule = lambda: cutsAndcategories(filename,args.year,True,args.eleID,args.muID,args.tauID,args.btauID)

        def mainModule(): return cutsAndcategories(
            filename, args.year, True, args.runNominal)
        # fatjetVarModule = lambda: filterSDmass_HbbtaggingAndWeights(filename,args.year,True)
        # wandzWtModule = lambda: wandzgenptweight(filename,args.year,True)
        # HTModule = lambda: HT_MT_vars(args.year,True)
        # delRModule = lambda: AddBunchofDelRs(args.year,True)
        # topPtweightModule = lambda: topPtweight(filename,args.year,True)
        # tauToAK8matchFlag = lambda: AddAK8MatchingToTaus(args.year,True)
        # additionalleptonveto = lambda: Standalone_AdditionalLepVetoFlag(args.year,True)
    else:
        print(("This is a ", args.year, " MC file = ", filename))

        def tesModule(): return TauEnergyScaleForHPSandBoosted(
            args.year, False, tauID_wp='Loose', ele_wp='VVLoose')
        # mainModule = lambda: cutsAndcategories(filename,args.year,False,args.eleID,args.muID,args.tauID,args.btauID)

        def mainModule(): return cutsAndcategories(
            filename, args.year, False, args.runNominal)
        # fatjetVarModule = lambda: filterSDmass_HbbtaggingAndWeights(filename,args.year,False)
        # wandzWtModule = lambda: wandzgenptweight(filename,args.year,False)
        # HTModule = lambda: HT_MT_vars(args.year,False)
        # delRModule = lambda: AddBunchofDelRs(args.year,False)
        # topPtweightModule = lambda: topPtweight(filename,args.year,False)
        # tauToAK8matchFlag = lambda: AddAK8MatchingToTaus(args.year,False)
        # additionalleptonveto = lambda: Standalone_AdditionalLepVetoFlag(args.year,False)

    if (args.Remove):
        if (search("Run", filename)):
            p = PostProcessor(
                args.tempStore,
                [files],
                cut=preselection,
                branchsel='Keep_Drop_txt/Data_remove_oldWtnaming.txt',
                modules=[
                    tesModule(),
                    mainModule()],
                postfix=post,
                noOut=False,
                outputbranchsel='Keep_Drop_txt/Data_remove_oldWtnaming.txt')  # ,histFileName="{}/cutflow2_{}.root".format(args.tempStore, filename),histDirName="cutflow2")
            # p = PostProcessor(args.tempStore,[files], cut=None, branchsel='Keep_Drop_txt/Data_remove_oldWtnaming.txt',modules=[tesModule(),mainModule(),fatjetVarModule(),wandzWtModule(),HTModule(),delRModule(),topPtweightModule(),tauToAK8matchFlag(),additionalleptonveto()], postfix=post,noOut=False,outputbranchsel='Keep_Drop_txt/Data_remove_oldWtnaming.txt')#,histFileName="{}/cutflow2_{}.root".format(args.tempStore, filename),histDirName="cutflow2")
            # p = PostProcessor(args.tempStore,[files], cut=None, branchsel='Keep_Drop_txt/Data_remove_oldWtnaming.txt',modules=[], postfix=post,noOut=False,outputbranchsel='Keep_Drop_txt/Data_remove_oldWtnaming.txt')

        else:
            p = PostProcessor(
                args.tempStore,
                [files],
                cut=preselection,
                branchsel='Keep_Drop_txt/MC_remove_oldWtnaming.txt',
                modules=[
                    tesModule(),
                    mainModule()],
                postfix=post,
                noOut=False,
                outputbranchsel='Keep_Drop_txt/MC_remove_oldWtnaming.txt')  # ,maxEntries=12000)#,histFileName="{}/cutflow2_{}.root".format(args.tempStore, filename),histDirName="cutflow2")
            # p = PostProcessor(args.tempStore,[files], cut=None, branchsel='Keep_Drop_txt/MC_remove_oldWtnaming.txt',modules=[tesModule(),mainModule(),fatjetVarModule(),wandzWtModule(),HTModule(),delRModule(),topPtweightModule(),tauToAK8matchFlag(),additionalleptonveto()], postfix=post,noOut=False,outputbranchsel='Keep_Drop_txt/MC_remove_oldWtnaming.txt')#,histFileName="{}/cutflow2_{}.root".format(args.tempStore, filename),histDirName="cutflow2")
            # p = PostProcessor(args.tempStore,[files], cut=None, branchsel='Keep_Drop_txt/MC_remove_oldWtnaming.txt',modules=[], postfix=post,noOut=False,outputbranchsel='Keep_Drop_txt/MC_remove_oldWtnaming.txt')
    else:
        p = PostProcessor(
            args.tempStore,
            [files],
            cut=None,
            branchsel=None,
            modules=[
                mainModule(),
                fatjetVarModule(),
                wandzWtModule()],
            postfix=post,
            noOut=False,
            outputbranchsel=None)
    p.run()
    print("###############MOVING THE OUTPUT FILE BACK TO HDFS#######################")
    if (args.transferOff):
        print(("Files are in temp area. NOT transferred to hdfs = ", args.tempStore))
    else:
        # This hadoop command doesnot work in singlarity - temporarily use mv
        # os.system("hadoop fs -moveFromLocal -f "+args.tempStore+"/"+filename+".root"+" "+outputDir+"/.") #This currently doesnot work due to username differences - it takes parida by default
        # Ensure the output directory exists
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        # Move the main ROOT file
        os.system("mv " + args.tempStore + "/" +
                  filename + ".root " + outputDir + "/.")
        print("Moved the ROOT file to output directory.")

        # Ensure the cutflow directory exists
        cutflowDir = os.path.join(outputDir, "cutflow_" + args.year)
        if not os.path.exists(cutflowDir):
            os.makedirs(cutflowDir)

        # Move the cutflow rootfile
        os.system("mv " + args.tempStore + "/" + "cutflow2_" +
                  filename + ".json " + cutflowDir + "/.")
        print("Moved the cutflow rootfile to the cutflow directory.")


if __name__ == "__main__":
    start_time = time.time()
    parser = argparse.ArgumentParser(
        description='Analysis framwork code part 2. Apply object selections, and get event categorization. FastMTT branches')
    parser.add_argument(
        '--inputLocation',
        '-i',
        help="enter the path to the location of input file set",
        default="")
    parser.add_argument(
        '--outputLocation',
        '-o',
        help="enter the path where yu want the output files to be stored",
        default=".")
    parser.add_argument(
        '--ncores',
        '-n',
        help="number of cores for parallel processing",
        default=1)
    parser.add_argument(
        '--postfix', help="string at the end of output file names", default="")
    parser.add_argument(
        '--year',
        '-y',
        help='specify the run - to make sure right triggers are used',
        choices=[
            '2016',
            '2016APV',
            '2017',
            '2018'])
    # parser.add_argument('--eleID','-eId',help="Specify the ID for the electrons to be selected 2-Loose, 3-Medium, 4-Tight",choices=['2','3','4'],default='2') #default is loose
    # parser.add_argument('--muID','-mId',help="Specify the ID for the muons to be selected 2-Loose, 3-Medium, 4-Tight (to keep it consistent with the electrons)",choices=['2','3','4'],default='2') #default is loose
    # parser.add_argument('--tauID','-tId',help="DeepTau (NOT A BIT MAP!!!) 1 = VVVLoose, 2 = VVLoose, 3 = VLoose, 4 = Loose, 5 = Medium, 6 = Tight, 7 = VTight, 8 = VVTight",choices=['1','2','3','4','5','6','7','8'],default='4') #default is loose
    # parser.add_argument('--btauID','-btId',help="DeepBoostedTau Raw
    # score",default='0.85')#default is 0.85
    parser.add_argument(
        '--tempStore',
        '-t',
        help='Temporary staging area for files before moving out to hdfs',
        required=True)
    parser.add_argument(
        '--transferOff',
        '-toff',
        help='with this flag the files will be not moved to hdfs from the temp area',
        action='store_true')
    parser.add_argument(
        '--runNominal',
        '-rnom',
        help='Run the script without the systematic variations',
        action='store_true')
    parser.add_argument(
        '--Remove',
        '-r',
        help="Remove [old naming convention branches] crossSectionWeighting, pileupWeighting, pileupWeight_UP, pileupWeight_DOWN, FinalWeighting, FinalWeighting_pileupWeight_UP, FinalWeighting_pileupWeight_DOWN",
        action='store_true')

    args = parser.parse_args()

    trigger_2016 = ["HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
                    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                    "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"]

    trigger_2017_18 = ["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                       "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"]

    if ((args.year == "2016APV") or (args.year == "2016")):
        preselection = "(" + "||".join(trigger_2016) + ")"
    elif ((args.year == "2017") or (args.year == "2018")):
        preselection = "(" + "||".join(trigger_2017_18) + ")"

    # fnames = glob.glob(args.inputLocation + "/*.root")  #making a list of input MC/DATA files
    # fnames = glob.glob(args.inputLocation + "/*TTToSemiLeptonic_2.root")
    fnames = glob.glob(args.inputLocation + "/*.root")
    outputDir = args.outputLocation

    post = args.postfix
    argList = list()
    for file in fnames:
        argList.append(file)

    if int(args.ncores) == 1:
        for arr in argList:
            print("Using a single thread ")
            call_postpoc(arr)

    else:
        try:
            print("Using Multithreading")
            pool = mp.Pool(int(args.ncores))
            print(("list", argList))
            res = pool.map(call_postpoc, argList)
        except Exception as error:
            print("MultiProc error - needs debugging")
            print(("An exception occurred:", type(error).__name__))

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(("Elapsed Time: {:.2f} seconds".format(elapsed_time)))
    elapsed_time_hours = elapsed_time / 3600  # Convert seconds to hours
    print(("Elapsed Time: {:.2f} hours".format(elapsed_time_hours)))
