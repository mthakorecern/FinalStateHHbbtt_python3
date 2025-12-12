#!/usr/bin/env python3

import os
import time
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import shutil
import math
import argparse
import json
import subprocess
from collections import OrderedDict

from FinalStateHHbbtt.fastMTTPython.fastMTTtool import *
from FinalStateHHbbtt.Pre_BackgroundEstimationWithSystematics.TauEnergyScaleModule.TauEnergyScaleForHPSandBoosted import TauEnergyScaleForHPSandBoosted

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NATModules.modules.fatjetId import fatJetId
from PhysicsTools.NATModules.modules.jetId import jetId
from PhysicsTools.NATModules.modules.jetVetoMap import jetVMAP
from PhysicsTools.NATModules.modules.fatJetvetoMap import fatJetVMAP
from PhysicsTools.NATModules.modules.xsWeights import XSWeightOnly

import math

class cutsAndcategories(Module):
    def __init__(self, filename, year, isData, runNominal=False, cutflowDir=None):

        self.runNominal = runNominal
        print("ACTIVATE: cutsAndcategories")
        self.year = year
        if self.year == "2024":
            self.year_unc = "2024"
        
        self.isMC = not isData
        self.isData = isData
        self.filename = filename
        self.cutflowDir = cutflowDir

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

        self.jetLeadFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.subjet1FV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.subjet2FV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)

        self.met = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)

        ## https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer24/#ak4-b-tagging
        if self.year == "2024":
            self.LooseJet = 0.0246
            self.MediumJet = 0.1272
            self.TightJet = 0.4678
        

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self.isMC:
            self.cutflow_dict = OrderedDict([
                ("Skimming Stage: Event-Level Sum of GenWeights", int(inputFile["cutflow"].GetBinContent(1))),                
                ("Skimming Stage: Runs Tree-Level Sum of GenWeights", int(inputFile["cutflow"].GetBinContent(2))),
                ("Skimming Stage: Events before any cuts", int(inputFile["cutflow"].GetBinContent(3))),                
                ("Skimming Stage: FatJet Requirement (nFatJet > 0)", int(inputFile["cutflow"].GetBinContent(4))),
                ("Skimming Stage: PuppiMET_pt Threshold (PuppiMET_pt > 120)", int(inputFile["cutflow"].GetBinContent(5))),
                ("Skimming Stage: Passing Filter: Flag_goodVertices", int(inputFile["cutflow"].GetBinContent(6))),
                ("Skimming Stage: Passing Filter: Flag_globalSuperTightHalo2016Filter", int(inputFile["cutflow"].GetBinContent(7))),
                ("Skimming Stage: Passing Filter: Flag_EcalDeadCellTriggerPrimitiveFilter", int(inputFile["cutflow"].GetBinContent(8))),
                ("Skimming Stage: Passing Filter: Flag_BadPFMuonFilter", int(inputFile["cutflow"].GetBinContent(9))),
                ("Skimming Stage: Passing Filter: Flag_BadPFMuonDzFilter", int(inputFile["cutflow"].GetBinContent(10))),
                ("Skimming Stage: Passing Filter: Flag_hfNoisyHitsFilter", int(inputFile["cutflow"].GetBinContent(11))),
                ("Skimming Stage: Passing Filter: Flag_eeBadScFilter", int(inputFile["cutflow"].GetBinContent(12))),
                ("Skimming Stage: Passing Filter: Flag_ecalBadCalibFilter", int(inputFile["cutflow"].GetBinContent(13))),
                ("Skimming Stage: Good Primary Vertices (PV_ndof > 4) && (abs(PV_z) < 24) && (sqrt(PV_x*PV_x+PV_y*PV_y) < 2)", int(inputFile["cutflow"].GetBinContent(14))),
                # ("Skimming Stage: Tau requirments (nboostedTau > 0) || (nTau > 0)", int(inputFile["cutflow"].GetBinContent(15))),
                ("Skimming Stage: Events available after Skimming:", inputTree.GetEntries()),
                
                ## Pre-selection Cuts
                ("Pre-selection: PuppiMET_pt > 120", 0),
                ("Events surviving the FatJet skim  (pt > 180, |eta| < 2.5, jetId>1)", 0),
                ("Events after FatJet Skimming and PuppiMET_pt >  120",0),
                ("Events after all object-level selections (before overlap cleaning)", 0),


                (" ...breakdown..> Atleast_one_Tau (Reco + ID + cleaning)", 0),
                # (" ...breakdown..> Atleast_one_Electron (Reco + IDnoIso + cleaning)", 0),
                (" ...breakdown..> Atleast_one_Muon (Reco + IDnoIso + cleaning)", 0),
                
                ("Atleast_2leptons_anykind", 0),
                #("Atleast_one_pair_anykind (Iso-cut applied for e/mu)", 0),
                (" ...breakdown..> Atleast_one_TauTau_pair", 0),
                # (" ...breakdown..> Atleast_one_TauElectron_pair", 0),
                (" ...breakdown..> Atleast_one_TauMuon_pair", 0),
                (" ...breakdown..> TT_channel (max pt pair)", 0),
                # (" ...breakdown..> ET_channel (max pt pair)", 0),
                (" ...breakdown..> MT_channel (max pt pair)", 0),
                
                #("DeltaR_LL<1.5 and abs(Hbb_met_phi) > 1 cut", 0),
                
                #("Visible Mass HTT > 20 cut", 0),
                ("Medium AK4 b-tag veto (events with 0 medium b-tagged jets)", 0),
                ("Final surviving events", 0),])

        else:
            self.cutflow_dict = OrderedDict([
                ("Skimming Stage: Events before any cuts", int(inputFile["cutflow"].GetBinContent(1))),                
                ("Skimming Stage: FatJet Requirement (nFatJet > 0)", int(inputFile["cutflow"].GetBinContent(2))),
                ("Skimming Stage: PuppiMET_pt Threshold (PuppiMET_pt > 120)", int(inputFile["cutflow"].GetBinContent(3))),
                ("Skimming Stage: Passing Filter: Flag_goodVertices", int(inputFile["cutflow"].GetBinContent(4))),
                ("Skimming Stage: Passing Filter: Flag_globalSuperTightHalo2016Filter", int(inputFile["cutflow"].GetBinContent(5))),
                ("Skimming Stage: Passing Filter: Flag_EcalDeadCellTriggerPrimitiveFilter", int(inputFile["cutflow"].GetBinContent(6))),
                ("Skimming Stage: Passing Filter: Flag_BadPFMuonFilter", int(inputFile["cutflow"].GetBinContent(7))),
                ("Skimming Stage: Passing Filter: Flag_BadPFMuonDzFilter", int(inputFile["cutflow"].GetBinContent(8))),
                ("Skimming Stage: Passing Filter: Flag_hfNoisyHitsFilter", int(inputFile["cutflow"].GetBinContent(9))),
                ("Skimming Stage: Passing Filter: Flag_eeBadScFilter", int(inputFile["cutflow"].GetBinContent(10))),
                ("Skimming Stage: Passing Filter: Flag_ecalBadCalibFilter", int(inputFile["cutflow"].GetBinContent(11))),
                ("Skimming Stage: Good Primary Vertices (PV_ndof > 4) && (abs(PV_z) < 24) && (sqrt(PV_x*PV_x+PV_y*PV_y) < 2)", int(inputFile["cutflow"].GetBinContent(12))),
                # ("Skimming Stage: Tau requirments (nboostedTau > 0) || (nTau > 0)", int(inputFile["cutflow"].GetBinContent(13))),
                ("Skimming Stage: Events available after Skimming:", inputTree.GetEntries()),
                
                ## Pre-selection Cuts
                ("Pre-selection: PuppiMET_pt > 120", 0),
                ("Events surviving the FatJet skim  (pt > 180, |eta| < 2.5, jetId>1)", 0),
                ("Events after FatJet Skimming and PuppiMET_pt >  120",0),
                ("Events after all object-level selections (before overlap cleaning)", 0),


                (" ...breakdown..> Atleast_one_Tau (Reco + ID + cleaning)", 0),
                # (" ...breakdown..> Atleast_one_Electron (Reco + IDnoIso + cleaning)", 0),
                (" ...breakdown..> Atleast_one_Muon (Reco + IDnoIso + cleaning)", 0),
                
                ("Atleast_2leptons_anykind", 0),
                #("Atleast_one_pair_anykind (Iso-cut applied for e/mu)", 0),
                (" ...breakdown..> Atleast_one_TauTau_pair", 0),
                # (" ...breakdown..> Atleast_one_TauElectron_pair", 0),
                (" ...breakdown..> Atleast_one_TauMuon_pair", 0),
                (" ...breakdown..> TT_channel (max pt pair)", 0),
                # (" ...breakdown..> ET_channel (max pt pair)", 0),
                (" ...breakdown..> MT_channel (max pt pair)", 0),
                
                #("DeltaR_LL<1.5 and abs(Hbb_met_phi) > 1 cut", 0),
                
                #("Visible Mass HTT > 20 cut", 0),
                ("Medium AK4 b-tag veto (events with 0 medium b-tagged jets)", 0),
                ("Final surviving events", 0),])


         
        if ((self.year == "2024") and (self.isMC)):
            self.totalEvents = int(inputTree.GetEntries())
            print(("Total Events in MC file = ", self.totalEvents))
        self.out = wrappedOutputTree
        self.out.branch("eventnominal", "I")
        self.out.branch("FatJet_globalParT3Xbb_mass", "F", lenVar="nFatJet")

        
        for sys in self.jesUnc:
            if sys == "":
                self.out.branch("Hbb_lep1_deltaR%s" % (sys), "F")
                self.out.branch("Hbb_lep2_deltaR%s" % (sys), "F")
                self.out.branch("softdropmass%s" % (sys), "F")
                self.out.branch("pnetmass%s" % (sys), "F")

                self.out.branch("HTT_m%s" % (sys), "F")
                self.out.branch("HTT_eta%s" % (sys), "F")
                self.out.branch("HTT_phi%s" % (sys), "F")
                self.out.branch("HTT_pt%s" % (sys), "F")
                self.out.branch("HTT_HPS_m%s" % (sys), "F")
                self.out.branch("HTT_HPS_eta%s" % (sys), "F")
                self.out.branch("HTT_HPS_phi%s" % (sys), "F")
                self.out.branch("HTT_boosted_m%s" % (sys), "F")
                self.out.branch("HTT_boosted_eta%s" % (sys), "F")
                self.out.branch("HTT_boosted_phi%s" % (sys), "F")

                self.out.branch("HTTvis_m%s" % (sys), "F")
                self.out.branch("HTTvis_eta%s" % (sys), "F")
                self.out.branch("HTTvis_phi%s" % (sys), "F")
                self.out.branch("HTTvis_pt%s" % (sys), "F")
                self.out.branch("HTTvis_deltaR%s" % (sys), "F")
                
                self.out.branch("HTTvis_HPS_m%s" % (sys), "F")
                self.out.branch("HTTvis_HPS_eta%s" % (sys), "F")
                self.out.branch("HTTvis_HPS_phi%s" % (sys), "F")
                self.out.branch("HTTvis_boosted_m%s" % (sys), "F")
                self.out.branch("HTTvis_boosted_eta%s" % (sys), "F")
                self.out.branch("HTTvis_boosted_phi%s" % (sys), "F")

                ### Lepton-specific Fast MTT Variables
                for prefix in ["HPS", "boosted"]:
                    for lep in ["Ele", "Mu"]:
                        self.out.branch(f"HTT_{prefix}_{lep}_m%s" % (sys), "F")
                        self.out.branch(f"HTT_{prefix}_{lep}_eta%s" % (sys), "F")
                        self.out.branch(f"HTT_{prefix}_{lep}_phi%s" % (sys), "F")

                self.out.branch("deltaR_tau_ele", "F")
                self.out.branch("deltaR_tau_mu", "F")
                self.out.branch("deltaR_tau1_tau2", "F")
                self.out.branch("deltaPhi_tau_ele", "F")
                self.out.branch("deltaPhi_tau_mu", "F")
                self.out.branch("deltaPhi_tau1_tau2", "F")
                
                self.out.branch("deltaR_hbb_httvis", "F")
                self.out.branch("deltaR_hbb_ak4", "F")
                self.out.branch("deltaPhi_hbb_ak4", "F")
                self.out.branch("deltaPhi_hbb_httvis", "F")
                self.out.branch("deltaR_hbb_htt", "F")
                self.out.branch("deltaPhi_hbb_htt", "F")
                self.out.branch("deltaPhi_hbb_leadingtau", "F")
                self.out.branch("deltaPhi_hbb_subleadingtau", "F")
                self.out.branch("deltaPhi_hbb_leadingele", "F")
                self.out.branch("deltaPhi_hbb_leadingmu", "F")
                

                self.out.branch("deltaPhi_met_tautau%s" % sys, "F")
                self.out.branch("deltaPhi_met_leadingtau%s" % sys, "F")
                self.out.branch("deltaPhi_met_subleadingtau%s" % sys, "F")
                self.out.branch("deltaPhi_met_leadingele%s" % sys, "F")
                self.out.branch("deltaPhi_met_leadingmu%s" % sys, "F")

                self.out.branch("deltaR_ak4_leadtau", "F")
                self.out.branch("deltaR_ak4_subtau", "F")
                self.out.branch("deltaR_ak4_ele", "F")
                self.out.branch("deltaR_ak4_mu", "F")
                self.out.branch("deltaR_httvis_ak4lead", "F")
                self.out.branch("deltaPhi_httvis_ak4lead", "F")
                self.out.branch("deltaPhi_met_ak4lead", "F")
                self.out.branch("deltaR_hbb_ak4lead", "F")
                self.out.branch("deltaPhi_hbb_ak4lead", "F")

                self.out.branch("deltaPhi_ak4_leadtau", "F")
                self.out.branch("deltaPhi_ak4_subtau", "F")
                self.out.branch("deltaPhi_ak4_ele", "F")
                self.out.branch("deltaPhi_ak4_mu", "F")

                self.out.branch("deltaR_subjets", "F")
                self.out.branch("deltaPhi_subjets", "F")

                self.out.branch("deltaR_subjet1_leadtau", "F")
                self.out.branch("deltaR_subjet1_subtau", "F")
                self.out.branch("deltaR_subjet1_ele", "F")
                self.out.branch("deltaR_subjet1_mu", "F")

                self.out.branch("deltaPhi_subjet1_leadtau", "F")
                self.out.branch("deltaPhi_subjet1_subtau", "F")
                self.out.branch("deltaPhi_subjet1_ele", "F")
                self.out.branch("deltaPhi_subjet1_mu", "F")

                self.out.branch("deltaR_subjet2_leadtau", "F")
                self.out.branch("deltaR_subjet2_subtau", "F")
                self.out.branch("deltaR_subjet2_ele", "F")
                self.out.branch("deltaR_subjet2_mu", "F")

                self.out.branch("deltaPhi_subjet2_leadtau", "F")
                self.out.branch("deltaPhi_subjet2_subtau", "F")
                self.out.branch("deltaPhi_subjet2_ele", "F")
                self.out.branch("deltaPhi_subjet2_mu", "F")
                self.out.branch("pt_balance_hbb_htt_abs",  "F")
                self.out.branch("pt_balance_hbb_htt_signed", "F")
                self.out.branch("fatjet_tau21", "F")
                self.out.branch("fatjet_tau32", "F")

                self.out.branch("subjet1_tau21", "F")
                self.out.branch("subjet1_tau32", "F")
                self.out.branch("subjet2_tau21", "F")
                self.out.branch("subjet2_tau32", "F")


                self.out.branch("Tau_rawDeepTauVSjet_logit", "F", lenVar="ngood_Taus")
                self.out.branch("boostedTau_rawDeepTauVSjet_logit", "F", lenVar="ngood_boostedTaus")
                
                self.out.branch("nallTaus%s" % (sys), "I")
                self.out.branch("allTaus_pt%s" % (sys), "F", lenVar="nallTaus%s" % (sys))
                self.out.branch("allTaus_eta%s" % (sys), "F", lenVar="nallTaus%s" % (sys))
                self.out.branch("allTaus_phi%s" % (sys), "F", lenVar="nallTaus%s" % (sys))
                self.out.branch("allTaus_mass%s" % (sys), "F", lenVar="nallTaus%s" % (sys))
                self.out.branch("allTaus_decayMode%s" % (sys), "F", lenVar="nallTaus%s" % (sys))

                self.out.branch("Xvis_m%s" % (sys), "F")
                self.out.branch("Xvis_eta%s" % (sys), "F")
                self.out.branch("Xvis_phi%s" % (sys), "F")
                self.out.branch("Xvis_pt%s" % (sys), "F")

            self.out.branch("channel%s" % (sys), "I")
            self.out.branch("boost%s" % (sys), "I")

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
            self.out.branch("ngood_LooseJets%s"%(sys),"I")
            self.out.branch("ngood_MediumJets%s" % (sys), "I")
            self.out.branch("ngood_TightJets%s" % (sys), "I")

            self.out.branch("index_gElectrons%s" % (sys), "I", lenVar="ngood_Electrons%s" % (sys))
            self.out.branch("index_gMuons%s" % (sys), "I", lenVar="ngood_Muons%s" % (sys))
            self.out.branch("index_gTaus%s" % (sys), "I", lenVar="ngood_Taus%s" % (sys))
            self.out.branch("index_gboostedTaus%s" % (sys), "I", lenVar="ngood_boostedTaus%s" % (sys))
            self.out.branch("index_gFatJets%s" % (sys), "I", lenVar="ngood_FatJets%s" % (sys))
            self.out.branch("index_gJets%s" % (sys), "I", lenVar="ngood_Jets%s" % (sys))
            self.out.branch("index_gLooseJets%s"%(sys),"I",lenVar="ngood_LooseJets%s"%(sys))
            self.out.branch("index_gMediumJets%s" % (sys), "I", lenVar="ngood_MediumJets%s" % (sys))
            self.out.branch("index_gTightJets%s" % (sys), "I", lenVar="ngood_TightJets%s" % (sys))

            self.out.branch("nDeltaR_boosted_HPS_preclean", "I")
            self.out.branch("nDeltaR_boosted_HPS_postclean", "I")
            self.out.branch("deltaR_boosted_HPS_preclean", "F", lenVar="nDeltaR_boosted_HPS_preclean")
            self.out.branch("deltaR_boosted_HPS_postclean", "F", lenVar="nDeltaR_boosted_HPS_postclean")
 
            
            self.out.branch("Hbb_met_phi%s" % (sys), "F")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        output_dir = os.path.dirname(outputFile.GetName())
        base_name = os.path.splitext(os.path.basename(outputFile.GetName()))[0]

        local_json_path = os.path.join(output_dir, f"{base_name}_cutflow.json")
        with open(local_json_path, "w") as outfile:
            json.dump(self.cutflow_dict, outfile, indent=4)
        print(f"[INFO] Cutflow JSON written locally: {local_json_path}")

        cutflow_dir = getattr(self, "cutflowDir", None)

        if cutflow_dir is None or cutflow_dir.strip() == "":
            print("[INFO] No cutflowDir provided. Defaulting to ROOT output directory.")
            cutflow_dir = output_dir  
        os.makedirs(cutflow_dir, exist_ok=True)

        final_json_path = os.path.join(cutflow_dir, f"{base_name}_cutflow.json")
        try:
            shutil.move(local_json_path, final_json_path)
            print(f"[INFO] Cutflow JSON moved to: {final_json_path}")
        except Exception as e:
            print(f"[ERROR] Failed to move cutflow JSON to {final_json_path}")
            print(f"       Reason: {e}")


    def analyze(self, event):
        def gettaupt(tau, sys):
            if (sys == "tesUp"):
                return tau.pt_tesUp
            elif (sys == "tesDown"):
                return tau.pt_tesDown
            else:
                return tau.pt

        def gettaumass(tau, sys):
            if (sys == "tesUp"):
                return tau.mass_tesUp
            elif (sys == "tesDown"):
                return tau.mass_tesDown
            else:
                return tau.mass

        def getjetpt(jet, sys):
            if ((sys == "") or (sys == "UnclustUp") or (
                    sys == "UnclustDown") or (sys == "tesUp") or (sys == "tesDown")):
                return jet.pt
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
                if (self.year_unc == "2024"):
                    return jet.pt_jesAbsolute_2024Up
            elif sys == "jesAbsolute_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.pt_jesAbsolute_2024Down
            elif sys == "jesBBEC1Up":
                return jet.pt_jesBBEC1Up
            elif sys == "jesBBEC1Down":
                return jet.pt_jesBBEC1Down
            elif sys == "jesBBEC1_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.pt_jesBBEC1_2024Up
            elif sys == "jesBBEC1_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.pt_jesBBEC1_2024Down
            elif sys == "jesEC2Up":
                return jet.pt_jesEC2Up
            elif sys == "jesEC2Down":
                return jet.pt_jesEC2Down
            elif sys == "jesEC2_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.pt_jesEC2_2024Up
            elif sys == "jesEC2_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.pt_jesEC2_2024Down
            elif sys == "jesFlavorQCDUp":
                return jet.pt_jesFlavorQCDUp
            elif sys == "jesFlavorQCDDown":
                return jet.pt_jesFlavorQCDDown
            elif sys == "jesHFUp":
                return jet.pt_jesHFUp
            elif sys == "jesHFDown":
                return jet.pt_jesHFDown
            elif sys == "jesHF_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.pt_jesHF_2024Up
            elif sys == "jesHF_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.pt_jesHF_2024Down
            elif sys == "jesRelativeBalUp":
                return jet.pt_jesRelativeBalUp
            elif sys == "jesRelativeBalDown":
                return jet.pt_jesRelativeBalDown
            elif sys == "jesRelativeSample_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.pt_jesRelativeSample_2024Up
            elif sys == "jesRelativeSample_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.pt_jesRelativeSample_2024Down
                
        def getjetmass(jet, sys):
            if ((sys == "") or (sys == "UnclustUp") or (
                    sys == "UnclustDown") or (sys == "tesUp") or (sys == "tesDown")):
                return jet.mass
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
                if (self.year_unc == "2024"):
                    return jet.mass_jesAbsolute_2024Up
            elif sys == "jesAbsolute_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.mass_jesAbsolute_2024Down
            elif sys == "jesBBEC1Up":
                return jet.mass_jesBBEC1Up
            elif sys == "jesBBEC1Down":
                return jet.mass_jesBBEC1Down
            elif sys == "jesBBEC1_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.mass_jesBBEC1_2024Up
            elif sys == "jesBBEC1_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.mass_jesBBEC1_2024Down
            elif sys == "jesEC2Up":
                return jet.mass_jesEC2Up
            elif sys == "jesEC2Down":
                return jet.mass_jesEC2Down
            elif sys == "jesEC2_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.mass_jesEC2_2024Up
            elif sys == "jesEC2_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.mass_jesEC2_2024Down
            elif sys == "jesFlavorQCDUp":
                return jet.mass_jesFlavorQCDUp
            elif sys == "jesFlavorQCDDown":
                return jet.mass_jesFlavorQCDDown
            elif sys == "jesHFUp":
                return jet.mass_jesHFUp
            elif sys == "jesHFDown":
                return jet.mass_jesHFDown
            elif sys == "jesHF_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.mass_jesHF_2024Up
            elif sys == "jesHF_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.mass_jesHF_2024Down
            elif sys == "jesRelativeBalUp":
                return jet.mass_jesRelativeBalUp
            elif sys == "jesRelativeBalDown":
                return jet.mass_jesRelativeBalDown
            elif sys == "jesRelativeSample_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.mass_jesRelativeSample_2024Up
            elif sys == "jesRelativeSample_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return jet.mass_jesRelativeSample_2024Down

        def getMETpt(sys):
            if ((sys == "") or (sys == "tesUp") or (sys == "tesDown")):
                return event.PuppiMET_pt
            elif sys == "jesTotalUp":
                return event.PuppiMET_pt_jesTotalUp
            elif sys == "jesTotalDown":
                return event.PuppiMET_pt_jesTotalDown
            elif sys == "jerUp":
                return event.PuppiMET_pt_jerUp
            elif sys == "jerDown":
                return event.PuppiMET_pt_jerDown
            elif sys == "jesAbsoluteUp":
                return event.PuppiMET_pt_jesAbsoluteUp
            elif sys == "jesAbsoluteDown":
                return event.PuppiMET_pt_jesAbsoluteDown
            elif sys == "jesAbsolute_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesAbsolute_2024Up
            elif sys == "jesAbsolute_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesAbsolute_2024Down
            elif sys == "jesBBEC1Up":
                return event.PuppiMET_pt_jesBBEC1Up
            elif sys == "jesBBEC1Down":
                return event.PuppiMET_pt_jesBBEC1Down
            elif sys == "jesBBEC1_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesBBEC1_2024Up
            elif sys == "jesBBEC1_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesBBEC1_2024Down
            elif sys == "jesEC2Up":
                return event.PuppiMET_pt_jesEC2Up
            elif sys == "jesEC2Down":
                return event.PuppiMET_pt_jesEC2Down
            elif sys == "jesEC2_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesEC2_2024Up
            elif sys == "jesEC2_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesEC2_2024Down
            elif sys == "jesFlavorQCDUp":
                return event.PuppiMET_pt_jesFlavorQCDUp
            elif sys == "jesFlavorQCDDown":
                return event.PuppiMET_pt_jesFlavorQCDDown
            elif sys == "jesHFUp":
                return event.PuppiMET_pt_jesHFUp
            elif sys == "jesHFDown":
                return event.PuppiMET_pt_jesHFDown
            elif sys == "jesHF_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesHF_2024Up
            elif sys == "jesHF_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesHF_2024Down
            elif sys == "jesRelativeBalUp":
                return event.PuppiMET_pt_jesRelativeBalUp
            elif sys == "jesRelativeBalDown":
                return event.PuppiMET_pt_jesRelativeBalDown
            elif sys == "jesRelativeSample_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesRelativeSample_2024Up
            elif sys == "jesRelativeSample_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_pt_jesRelativeSample_2024Down
            elif sys == "UnclustUp":
                return event.PuppiMET_ptUnclusteredUp
            elif sys == "UnclustDown":
                return event.PuppiMET_ptUnclusteredDown

        def getMETphi(sys):
            if ((sys == "")  or (sys == "UnclustUp") or (
                    sys == "UnclustDown") or (sys == "tesUp") or (sys == "tesDown")):
                return event.PuppiMET_phi
            elif sys == "jesTotalUp":
                return event.PuppiMET_phi_jesTotalUp
            elif sys == "jesTotalDown":
                return event.PuppiMET_phi_jesTotalDown
            elif sys == "jerUp":
                return event.PuppiMET_phi_jerUp
            elif sys == "jerDown":
                return event.PuppiMET_phi_jerDown
            elif sys == "jesAbsoluteUp":
                return event.PuppiMET_phi_jesAbsoluteUp
            elif sys == "jesAbsoluteDown":
                return event.PuppiMET_phi_jesAbsoluteDown
            elif sys == "jesAbsolute_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesAbsolute_2024Up
            elif sys == "jesAbsolute_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesAbsolute_2024Down
            elif sys == "jesBBEC1Up":
                return event.PuppiMET_phi_jesBBEC1Up
            elif sys == "jesBBEC1Down":
                return event.PuppiMET_phi_jesBBEC1Down
            elif sys == "jesBBEC1_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesBBEC1_2024Up
            elif sys == "jesBBEC1_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesBBEC1_2024Down
            elif sys == "jesEC2Up":
                return event.PuppiMET_phi_jesEC2Up
            elif sys == "jesEC2Down":
                return event.PuppiMET_phi_jesEC2Down
            elif sys == "jesEC2_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesEC2_2024Up
            elif sys == "jesEC2_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesEC2_2024Down
            elif sys == "jesFlavorQCDUp":
                return event.PuppiMET_phi_jesFlavorQCDUp
            elif sys == "jesFlavorQCDDown":
                return event.PuppiMET_phi_jesFlavorQCDDown
            elif sys == "jesHFUp":
                return event.PuppiMET_phi_jesHFUp
            elif sys == "jesHFDown":
                return event.PuppiMET_phi_jesHFDown
            elif sys == "jesHF_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesHF_2024Up
            elif sys == "jesHF_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesHF_2024Down
            elif sys == "jesRelativeBalUp":
                return event.PuppiMET_phi_jesRelativeBalUp
            elif sys == "jesRelativeBalDown":
                return event.PuppiMET_phi_jesRelativeBalDown
            elif sys == "jesRelativeSample_%sUp" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesRelativeSample_2024Up
            elif sys == "jesRelativeSample_%sDown" % (self.year_unc):
                if (self.year_unc == "2024"):
                    return event.PuppiMET_phi_jesRelativeSample_2024Down
            elif sys == "UnclustUp":
                return event.PuppiMET_phiUnclusteredUp
            elif sys == "UnclustDown":
                return event.PuppiMET_phiUnclusteredDown

        def applyPOGselectionToAK4(ak4Object_enu, sys):
            """
            Select central jets (2024 baseline):
            - pt > 30 GeV
            - |eta| < 2.5
            - jetId > 1
            """

            if getjetpt(ak4Object_enu[1], sys) > 30 and abs(ak4Object_enu[1].eta) < 2.5 and ak4Object_enu[1].jetId > 1: 
                #and getattr(event, "Flag_JetVetoed", 0) == 0:
                return True
            return False

        def logit(x):
            return math.log(x/(1 - x))

        # def pass_cuts_EleID(electronObject_enu):
        #     for cutnr in range(0, 10):
        #         if cutnr == 7:
        #             continue
        #         # if (electronObject_enu[1].vidNestedWPBitmap >> (cutnr*3) &
        #         # 0x7) < self.eleID:
        #         if (electronObject_enu[1].vidNestedWPBitmap >> (
        #                 cutnr * 3) & 0x7) < 2:
        #             return False
        #     return True

        
        
        # FatJet and Jet overlap, separation > 1.2 
        def JetFatJetOverlap(jetObject_enu, sys):
            self.jetFV.SetPtEtaPhiM(getjetpt(jetObject_enu[1],sys), jetObject_enu[1].eta, jetObject_enu[1].phi, getjetmass(jetObject_enu[1], sys))
            deltaR = self.jetFV.DeltaR(self.higgsBBFV)
            if deltaR > 1.2: #(0.8 + 0.4)
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
            if deltaR > 0.8: #(0.8)
                return True
            else:
                return False

        # Function used by Tau(boostedTau) cleaning vs AK8 > 0.8 
        def FatJetTauOverlap(tauObject_enu, boost):
            if boost == 1:
                self.tauFV.SetPtEtaPhiM(tauObject_enu[1].pt, tauObject_enu[1].eta, tauObject_enu[1].phi, tauObject_enu[1].mass)
            else:
                self.tauFV.SetPtEtaPhiM(tauObject_enu[1].pt, tauObject_enu[1].eta, tauObject_enu[1].phi, tauObject_enu[1].mass)
            deltaR = self.tauFV.DeltaR(self.higgsBBFV)
            if deltaR > 1.5: #(0.8+0.4) but taking slightly conservative approach and going for 1.5
                return True
            else:
                return False

        # Function used by Tau(boostedTau) cleaning vs Electrons > 0.05. This is to retain the Taus, which are close to Electrons for pairing
        def ElectronTauOverlap(tauObject_enu, elecoll_enu, boost):
            if boost == 1:
                self.tauFV.SetPtEtaPhiM(tauObject_enu[1].pt, tauObject_enu[1].eta, tauObject_enu[1].phi,tauObject_enu[1].mass)
            else:
                self.tauFV.SetPtEtaPhiM(tauObject_enu[1].pt, tauObject_enu[1].eta, tauObject_enu[1].phi, tauObject_enu[1].mass)
            for electron in elecoll_enu:
                self.eleFV.SetPtEtaPhiM(electron[1].pt, electron[1].eta, electron[1].phi, 0.0)
                deltaR = self.tauFV.DeltaR(self.eleFV)
                if deltaR <= 0.05: 
                    return False
            return True

        # Function used by Tau(boostedTau) cleaning vs Muons > 0.05. This is to retain the Muons, which are close to Electrons for pairing
        def MuonTauOverlap(tauObject_enu, mucoll_enu, boost):
            if boost == 1:
                self.tauFV.SetPtEtaPhiM(tauObject_enu[1].pt,tauObject_enu[1].eta, tauObject_enu[1].phi, tauObject_enu[1].mass)
            else:
                self.tauFV.SetPtEtaPhiM(tauObject_enu[1].pt, tauObject_enu[1].eta, tauObject_enu[1].phi, tauObject_enu[1].mass)
            for muon in mucoll_enu:
                self.muFV.SetPtEtaPhiM(muon[1].pt, muon[1].eta, muon[1].phi, muon[1].mass)
                deltaR = self.tauFV.DeltaR(self.muFV)
                if deltaR <= 0.05:
                    return False
            return True

        def removeOverlapOfAK4WithLightHeavyLeptons(ak4Object_enu, gTau_index, Tau_coll, gbTau_index, bTau_coll, gEle_index, Ele_coll, gMu_index,   Mu_coll, sys):
            self.jetFV.SetPtEtaPhiM(getjetpt(ak4Object_enu[1], sys), ak4Object_enu[1].eta, ak4Object_enu[1].phi, getjetmass(ak4Object_enu[1], sys))
            # for index in gEle_index:
            #     self.eleFV.SetPtEtaPhiM(Ele_coll[index].pt, Ele_coll[index].eta, Ele_coll[index].phi, 0.0)
            #     deltaR = self.jetFV.DeltaR(self.eleFV)
            #     if deltaR <= 0.4:
            #         return False

            for index in gMu_index:
                self.muFV.SetPtEtaPhiM(Mu_coll[index].pt, Mu_coll[index].eta, Mu_coll[index].phi, Mu_coll[index].mass)
                deltaR = self.jetFV.DeltaR(self.muFV)
                if deltaR <= 0.4:
                    return False

            for index in gTau_index:
                self.tauFV.SetPtEtaPhiM(Tau_coll[index].pt, Tau_coll[index].eta, Tau_coll[index].phi, Tau_coll[index].mass)
                deltaR = self.jetFV.DeltaR(self.tauFV)
                if deltaR <= 0.4:
                    return False

            for index in gbTau_index:
                self.tauFV.SetPtEtaPhiM(bTau_coll[index].pt, bTau_coll[index].eta, bTau_coll[index].phi,bTau_coll[index].mass)
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
                    
                    ## Commenting out the isolation part for now. Might add it in the future.

                    # if ((tag == "be") or (tag == "te")):
                    #     if not ElectronIsolationCut(
                    #             col1[i][1], col2[j][1], tag):
                    #         continue

                    # elif ((tag == "bm") or (tag == "tm")):
                    #     if not MuonIsolationCut(col1[i][1], col2[j][1], tag):
                    #         continue

                    sumFourVector = self.pair1FV + self.pair2FV
                    pt = sumFourVector.Pt()
                    if pt >= combinedPt:
                        combinedPt = pt
                        index1 = col1[i][0]
                        index2 = col2[j][0]
            return (combinedPt, index1, index2)

        # def MuonIsolationCut(tau, muo, tag):
        #     tagprotection = ((tag == "bm") or (tag == "tm"))
        #     if not tagprotection:
        #         print("The tags are not correct in Electron Isolation ##ERROR##------")
        #         sys.exit()

        #     isTau = ""

        #     self.tauFV.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass)
        #     self.muFV.SetPtEtaPhiM(muo.pt, muo.eta, muo.phi, muo.mass)
        #     deltaR = (self.tauFV).DeltaR(self.muFV)
        #     isolationCut = 0.25

        #     if (tag == "tm"):
        #         if ((muo.pfRelIso04_all) < isolationCut):
        #             return True
        #         else:
        #             return False

        #     if deltaR < 0.7:
        #         isTau = "close"
        #     if isTau == "close":
        #         self.leadingMatch.SetPtEtaPhiM(tau.LeadingMuonPt, tau.LeadingMuonEta, tau.LeadingMuonPhi,tau.LeadingMuonM)
        #         self.subleadingMatch.SetPtEtaPhiM(tau.SubLeadingMuonPt,tau.SubLeadingMuonEta,tau.SubLeadingMuonPhi,tau.SubLeadingMuonM)
        #         self.subsubleadingMatch.SetPtEtaPhiM(tau.SubSubLeadingMuonPt, tau.SubSubLeadingMuonEta,tau.SubSubLeadingMuonPhi, tau.SubSubLeadingMuonM)
        #         if (self.muFV.DeltaR(self.leadingMatch) < 0.05):
        #             if ((tau.LeadingMuonCorrIso / tau.LeadingMuonPt) < isolationCut):
        #                 # print ("Leading Muon Matched")
        #                 return True
        #             else:
        #                 return False
        #         elif (self.muFV.DeltaR(self.subleadingMatch) < 0.05):
        #             if ((tau.SubLeadingMuonCorrIso /
        #                  tau.SubLeadingMuonPt) < isolationCut):
        #                 # print ("subLeading Muon Matched")
        #                 return True
        #             else:
        #                 return False
        #         elif (self.muFV.DeltaR(self.subsubleadingMatch) < 0.05):
        #             if ((tau.SubSubLeadingMuonCorrIso /
        #                  tau.SubSubLeadingMuonPt) < isolationCut):
        #                 # print ("subsubLeading Muon Matched")
        #                 return True
        #             else:
        #                 return False
        #         else:
        #             if ((muo.pfRelIso04_all) < isolationCut):
        #                 # print ("pf Muon Isolation")
        #                 return True
        #             else:
        #                 return False

        #     elif isTau == "":
        #         if ((muo.pfRelIso04_all) < isolationCut):
        #             return True
        #         else:
        #             return False

        # passing the inidivial eletrons from the collection to apply the
        # correction
        # def ElectronIsolationCut(tau, ele, tag):
        #     tagprotection = ((tag == "be") or (tag == "te"))
        #     if not tagprotection:
        #         print("The tags are not correct in Electron Isolation ##ERROR##------")
        #         sys.exit()

        #     isTau = ""
        #     isolationCut = 0.0

        #     if abs(ele.eta) <= 1.479:
        #         isolationCut = 0.194 + (0.535 / ele.pt)
        #     elif (abs(ele.eta) > 1.479) and (abs(ele.eta) <= 2.5):
        #         # Endcap values
        #         # loose =  0.108 + (0.963/ele.pt)
        #         # medium = 0.0658 + (0.963/ele.pt)
        #         # tight =0.0445 + (0.963/ele.pt)
        #         isolationCut = 0.184 + (0.519 / ele.pt)
        #     else:
        #         return False

        #     if (tag == "te"):
        #         if ((ele.pfRelIso03_all) < isolationCut):
        #             return True
        #         else:
        #             return False

        #     self.tauFV.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass)
        #     self.eleFV.SetPtEtaPhiM(ele.pt, ele.eta, ele.phi, 0.0)
        #     deltaR = (self.tauFV).DeltaR(self.eleFV)

        #     if deltaR < 0.6:
        #         isTau = "close"
        #     if isTau == "close":
        #         self.leadingMatch.SetPtEtaPhiM(
        #             tau.LeadingElectronPt,
        #             tau.LeadingElectronEta,
        #             tau.LeadingElectronPhi,
        #             0.0)
        #         self.subleadingMatch.SetPtEtaPhiM(
        #             tau.SubLeadingElectronPt,
        #             tau.SubLeadingElectronEta,
        #             tau.SubLeadingElectronPhi,
        #             0.0)
        #         self.subsubleadingMatch.SetPtEtaPhiM(
        #             tau.SubSubLeadingElectronPt,
        #             tau.SubSubLeadingElectronEta,
        #             tau.SubSubLeadingElectronPhi,
        #             0.0)
        #         if (self.eleFV.DeltaR(self.leadingMatch) < 0.05):
        #             if ((tau.LeadingElectronCorrIso /
        #                  tau.LeadingElectronPt) < isolationCut):
        #                 # print ("Leading Ele Matched")
        #                 return True
        #             else:
        #                 return False
        #         elif (self.eleFV.DeltaR(self.subleadingMatch) < 0.05):
        #             if ((tau.SubLeadingElectronCorrIso /
        #                  tau.SubLeadingElectronPt) < isolationCut):
        #                 # print ("subLeading Ele Matched")
        #                 return True
        #             else:
        #                 return False
        #         elif (self.eleFV.DeltaR(self.subsubleadingMatch) < 0.05):
        #             if ((tau.SubSubLeadingElectronCorrIso /
        #                  tau.SubSubLeadingElectronPt) < isolationCut):
        #                 # print ("subsubLeading Ele Matched")
        #                 return True
        #             else:
        #                 return False
        #         else:
        #             if ((ele.pfRelIso03_all) < isolationCut):
        #                 print("use pf Ele Isolation")
        #                 return True
        #             else:
        #                 return False

        #     elif isTau == "":
        #         if ((ele.pfRelIso03_all) < isolationCut):
        #             return True
        #         else:
        #             return False

        # def ElectronIsolationCut_addlepveto(eleObject_enu):
        #     #For 2024: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun3 
        #     isolationCut = 0.0
        #     if abs(eleObject_enu[1].eta) <= 1.479:
        #         isolationCut = 0.194 + (0.535 / eleObject_enu[1].pt)
        #     elif (abs(eleObject_enu[1].eta) > 1.479) and (abs(eleObject_enu[1].eta) <= 2.5):
        #         isolationCut = 0.184 + (0.519 / eleObject_enu[1].pt)
        #     else:
        #         return False
        #     if ((eleObject_enu[1].pfRelIso03_all) < isolationCut):
        #         return True
        #     else:
        #         return False
        
        
        def fillBranchesWithDefault(sys):
            if sys == "":
                self.out.fillBranch("Hbb_lep1_deltaR%s" % (sys), -99.99)
                self.out.fillBranch("Hbb_lep2_deltaR%s" % (sys), -99.99)
                # self.out.fillBranch("softdropmassnom%s" % (sys), -1.00)
                self.out.fillBranch("softdropmass%s" % (sys), -99.99)
                self.out.fillBranch("pnetmass%s" % (sys), -99.99)

                self.out.fillBranch("HTT_m%s" % (sys), -99.99)
                self.out.fillBranch("HTT_eta%s" % (sys), -99.99)
                self.out.fillBranch("HTT_phi%s" % (sys), -99.99)
                self.out.fillBranch("HTT_pt%s" % (sys), -99.99)
                self.out.fillBranch("HTTvis_m%s" % (sys), -99.99)
                self.out.fillBranch("HTTvis_eta%s" % (sys), -99.99)
                self.out.fillBranch("HTTvis_phi%s" % (sys), -99.99)
                self.out.fillBranch("HTTvis_pt%s" % (sys), -99.99)
                self.out.fillBranch("HTTvis_deltaR%s" % (sys), -99.99)

                self.out.fillBranch("HTT_HPS_m", -99.99)
                self.out.fillBranch("HTT_HPS_eta", -99.99)
                self.out.fillBranch("HTT_HPS_phi", -99.99)

                self.out.fillBranch("HTT_boosted_m", -99.99)
                self.out.fillBranch("HTT_boosted_eta", -99.99)
                self.out.fillBranch("HTT_boosted_phi", -99.99)

                self.out.fillBranch("HTTvis_HPS_m", -99.99)
                self.out.fillBranch("HTTvis_HPS_eta", -99.99)
                self.out.fillBranch("HTTvis_HPS_phi", -99.99)
                
                self.out.fillBranch("HTTvis_boosted_m", -99.99)
                self.out.fillBranch("HTTvis_boosted_eta", -99.99)
                self.out.fillBranch("HTTvis_boosted_phi", -99.99)

                for prefix in ["HPS", "boosted"]:
                    for lep in ["Ele", "Mu"]:
                        self.out.fillBranch(f"HTT_{prefix}_{lep}_m", -99.99)
                        self.out.fillBranch(f"HTT_{prefix}_{lep}_eta", -99.99)
                        self.out.fillBranch(f"HTT_{prefix}_{lep}_phi", -99.99)

                self.out.fillBranch("deltaR_tau_ele", -99.99)
                self.out.fillBranch("deltaR_tau_mu", -99.99)              
                self.out.fillBranch("deltaPhi_tau_ele", -99.99)
                self.out.fillBranch("deltaPhi_tau_mu", -99.99)
  
                
                self.out.fillBranch("deltaPhi_met_tautau", -99.99)
                self.out.fillBranch("deltaPhi_met_leadingtau", -99.99)
                self.out.fillBranch("deltaPhi_met_subleadingtau", -99.99)
                self.out.fillBranch("deltaPhi_met_leadingele", -99.99)
                self.out.fillBranch("deltaPhi_met_leadingmu", -99.99)

                self.out.fillBranch("deltaR_hbb_httvis", -99.99)
                self.out.fillBranch("deltaPhi_hbb_httvis", -99.99)
                self.out.fillBranch("deltaPhi_hbb_htt", -99.99)
                self.out.fillBranch("deltaPhi_hbb_leadingtau", -99.99)
                self.out.fillBranch("deltaPhi_hbb_subleadingtau", -99.99)
                self.out.fillBranch("deltaPhi_hbb_leadingele", -99.99)
                self.out.fillBranch("deltaPhi_hbb_leadingmu", -99.99)

                self.out.fillBranch("deltaPhi_tau1_tau2", -99.99)
                self.out.fillBranch("deltaR_tau1_tau2", -99.99)

                self.out.fillBranch("deltaR_httvis_ak4lead", -99.99)
                self.out.fillBranch("deltaPhi_httvis_ak4lead", -99.99)
                self.out.fillBranch("deltaR_hbb_ak4lead", -99.99)
                self.out.fillBranch("deltaPhi_hbb_ak4lead", -99.99)
                self.out.fillBranch("deltaR_ak4_leadtau", -99.99)
                self.out.fillBranch("deltaPhi_ak4_leadtau", -99.99)
                self.out.fillBranch("deltaR_ak4_subtau", -99.99)
                self.out.fillBranch("deltaPhi_ak4_subtau", -99.99)
                self.out.fillBranch("deltaR_ak4_ele", -99.99)
                self.out.fillBranch("deltaPhi_ak4_ele", -99.99)
                self.out.fillBranch("deltaR_ak4_mu", -99.99)
                self.out.fillBranch("deltaPhi_ak4_mu", -99.99)
                self.out.fillBranch("deltaPhi_met_ak4lead", -99.99)


                self.out.fillBranch("deltaR_subjets", -99.99)
                self.out.fillBranch("deltaPhi_subjets", -99.99)
                
                self.out.fillBranch("deltaR_subjet1_leadtau", -99.99)
                self.out.fillBranch("deltaPhi_subjet1_leadtau", -99.99)
                self.out.fillBranch("deltaR_subjet1_subtau", -99.99)
                self.out.fillBranch("deltaPhi_subjet1_subtau", -99.99)
                self.out.fillBranch("deltaR_subjet1_ele", -99.99)
                self.out.fillBranch("deltaPhi_subjet1_ele", -99.99)
                self.out.fillBranch("deltaR_subjet1_mu", -99.99)
                self.out.fillBranch("deltaPhi_subjet1_mu", -99.99)

                self.out.fillBranch("deltaR_subjet2_leadtau", -99.99)
                self.out.fillBranch("deltaPhi_subjet2_leadtau", -99.99)
                self.out.fillBranch("deltaR_subjet2_subtau", -99.99)
                self.out.fillBranch("deltaPhi_subjet2_subtau", -99.99)
                self.out.fillBranch("deltaR_subjet2_ele", -99.99)
                self.out.fillBranch("deltaPhi_subjet2_ele", -99.99)
                self.out.fillBranch("deltaR_subjet2_mu", -99.99)
                self.out.fillBranch("deltaPhi_subjet2_mu", -99.99)

                self.out.fillBranch("fatjet_tau21", -99.99)
                self.out.fillBranch("fatjet_tau32", -99.99)

                self.out.fillBranch("subjet1_tau21", -99.99)
                self.out.fillBranch("subjet1_tau32", -99.99)
                self.out.fillBranch("subjet2_tau21", -99.99)
                self.out.fillBranch("subjet2_tau32", -99.99)

                self.out.fillBranch("Tau_rawDeepTauVSjet_logit", [])
                self.out.fillBranch("boostedTau_rawDeepTauVSjet_logit", [])

                self.out.fillBranch("pt_balance_hbb_htt_abs", -99.99)
                self.out.fillBranch("pt_balance_hbb_htt_signed", -99.99)



                self.out.fillBranch("nallTaus%s" % (sys), 0)
                self.out.fillBranch("allTaus_pt%s" % (sys), [])
                self.out.fillBranch("allTaus_eta%s" % (sys), [])
                self.out.fillBranch("allTaus_phi%s" % (sys), [])
                self.out.fillBranch("allTaus_mass%s" % (sys), [])
                self.out.fillBranch("allTaus_decayMode%s" % (sys), [])
                self.out.fillBranch("Xvis_m%s" % (sys), -99.99)
                self.out.fillBranch("Xvis_eta%s" % (sys), -99.99)
                self.out.fillBranch("Xvis_phi%s" % (sys), -99.99)
                self.out.fillBranch("Xvis_pt%s" % (sys), -99.99)

            self.out.fillBranch("channel%s" % (sys), -1)
            self.out.fillBranch("boost%s" % (sys), -1)

            self.out.fillBranch("Hbb_met_phi%s" % (sys), -99.99)
            self.out.fillBranch("X_m%s" % (sys), -99.99)
            self.out.fillBranch("X_eta%s" % (sys), -99.99)
            self.out.fillBranch("X_phi%s" % (sys), -99.99)
            self.out.fillBranch("X_pt%s" % (sys), -99.99)

            self.out.fillBranch("ngood_Taus%s" % (sys), 0)
            self.out.fillBranch("ngood_boostedTaus%s" % (sys), 0)
            self.out.fillBranch("ngood_Electrons%s" % (sys), 0)
            self.out.fillBranch("ngood_Muons%s" % (sys), 0)
            self.out.fillBranch("ngood_FatJets%s" % (sys), 0)
            self.out.fillBranch("ngood_Jets%s" % (sys), 0)
            self.out.fillBranch("ngood_LooseJets%s"%(sys),0)
            self.out.fillBranch("ngood_MediumJets%s" % (sys), 0)
            self.out.fillBranch("ngood_TightJets%s" % (sys), 0)

            self.out.fillBranch("index_gElectrons%s" % (sys), [])
            self.out.fillBranch("index_gMuons%s" % (sys), [])
            self.out.fillBranch("index_gTaus%s" % (sys), [])
            self.out.fillBranch("index_gboostedTaus%s" % (sys), [])
            self.out.fillBranch("index_gFatJets%s" % (sys), [])
            self.out.fillBranch("index_gJets%s" % (sys), [])
            self.out.fillBranch("index_gLooseJets%s"%(sys),[])
            self.out.fillBranch("index_gMediumJets%s" % (sys), [])
            self.out.fillBranch("index_gTightJets%s" % (sys), [])

        Jet = Collection(event, 'Jet', 'nJet')
        Electron = Collection(event, 'Electron', 'nElectron')
        FatJet = Collection(event, 'FatJet', 'nFatJet')

        if (self.isMC):
            if (self.year_unc == "2024"):          
                if event.PuppiMET_pt < 120:# and (
                    return False
        elif (self.isData):
            if ((event.PuppiMET_pt < 120)):
                return False
        
        self.cutflow_dict["Pre-selection: PuppiMET_pt > 120"] += 1

       
        FatJet_skim_enu = [x for x in enumerate(FatJet) if (x[1].pt > 180) and abs(x[1].eta) < 2.5 and (x[1].jetId > 1)]
        if (len(FatJet_skim_enu) == 0):
            return False

        self.cutflow_dict["Events surviving the FatJet skim  (pt > 180, |eta| < 2.5, jetId>1)"] += 1
        
        del FatJet_skim_enu

        passAsingleSystematic = 0

        Muon = Collection(event, 'Muon', 'nMuon')
        Tau = Collection(event, 'Tau', 'nTau')
        boostedTau = Collection(event, 'boostedTau', 'nboostedTau')

        nominal_bool = 0
        for sys in self.jesUnc:
            fillBranchesWithDefault(sys)

            
            del self.theFastMTTtool
            self.theFastMTTtool = fastMTTtool()

            for vec in [self.lepFV, self.eleFV, self.muFV, self.tauFV,
                        self.jetFV, self.leadingMatch, self.subleadingMatch,
                        self.subsubleadingMatch, self.higgsTTFV, self.higgsTTvisFV,
                        self.higgsBBFV, self.RadionFV, self.RadionvisFV,
                        self.pair1FV, self.pair2FV, self.met,
                        self.jetLeadFV, self.subjet1FV, self.subjet2FV]:
                vec.SetPxPyPzE(0, 0, 0, 0)
        
            if ((getMETpt(sys) < 120)):
                # fillBranchesWithDefault(sys)
                continue

            if sys == "":
                self.cutflow_dict["Events after FatJet Skimming and PuppiMET_pt >  120"] += 1

            FatJet_enu = [x for x in enumerate(FatJet) if (getjetpt(x[1], sys) > 180) and (abs(x[1].eta) < 2.5 and (x[1].jetId > 1))]

            if (len(FatJet_enu) == 0):
                # fillBranchesWithDefault(sys)
                continue

            HbbPtList = [(getjetpt(obj_enu[1], sys)) for obj_enu in FatJet_enu]
            zipPair = list(zip(FatJet_enu, HbbPtList))
            FatJet_enu = [
                fatjetobject_enu for fatjetobject_enu,
                _ in sorted(
                    zipPair,
                    key=lambda x: x[1],
                    reverse=True)]

            del HbbPtList
            del zipPair
                        
            Jet_enu = [x for x in enumerate(Jet) if applyPOGselectionToAK4(x, sys)]
            
            Tau_enu = [x for x in enumerate(Tau) if (gettaupt(x[1], sys) > 20) and (abs(x[1].eta) < 2.5) and (abs(x[1].dz) < 0.2) and (x[1].idDecayModeNewDMs) and (x[1].idDeepTau2018v2p5VSjet >= 4) and (x[1].idDeepTau2018v2p5VSe >= 2) and (x[1].idDeepTau2018v2p5VSmu >= 1)]
            
            boostedTau_enu = [x for x in enumerate(boostedTau) if (gettaupt(x[1], sys) > 25) and (abs(x[1].eta) < 2.5) and (x[1].rawBoostedDeepTauRunIIv2p0VSjet >= 0.85)]

            self.higgsBBFV.SetPtEtaPhiM(getjetpt(FatJet_enu[0][1],sys), FatJet_enu[0][1].eta, FatJet_enu[0][1].phi, getjetmass(FatJet_enu[0][1], sys))

            # Electron_enu = [x for x in enumerate(Electron) if x[1].pt > 10 and (abs(x[1].eta) < 2.5) and x[1].cutBased >= 2]
            
            Muon_enu = [x for x in enumerate(Muon) if x[1].pt > 15 and (abs(x[1].eta) < 2.4) and x[1].looseId]

            deltaR_boosted_HPS_preclean = []
            for b in boostedTau_enu:
                for h in Tau_enu:
                    bfv = ROOT.TLorentzVector()
                    hfv = ROOT.TLorentzVector()
                    bfv.SetPtEtaPhiM(gettaupt(b[1], sys), b[1].eta, b[1].phi, gettaumass(b[1], sys))
                    hfv.SetPtEtaPhiM(gettaupt(h[1], sys), h[1].eta, h[1].phi, gettaumass(h[1], sys))
                    deltaR_boosted_HPS_preclean.append(bfv.DeltaR(hfv))

            self.out.fillBranch("nDeltaR_boosted_HPS_preclean", len(deltaR_boosted_HPS_preclean))
            self.out.fillBranch("deltaR_boosted_HPS_preclean", deltaR_boosted_HPS_preclean)

                        
            # Object cleaning procedures
            if sys == "":  # only count once for nominal
                self.cutflow_dict["Events after all object-level selections (before overlap cleaning)"] += 1

            Jet_enu = [x for x in Jet_enu if JetFatJetOverlap(x, sys)]
            # Electron_enu = list(filter(FatJetConeIsolation, Electron_enu))
            Muon_enu = list(filter(FatJetConeIsolation, Muon_enu))
            Tau_enu = [x for x in Tau_enu if FatJetTauOverlap(x, boost=0)]
            # Tau_enu = [x for x in Tau_enu if ElectronTauOverlap(x, Electron_enu, boost=0)]
            Tau_enu = [x for x in Tau_enu if MuonTauOverlap(x, Muon_enu, boost=0)]
            boostedTau_enu = [x for x in boostedTau_enu if FatJetTauOverlap(x, boost=1)]
            # boostedTau_enu = [x for x in boostedTau_enu if ElectronTauOverlap(x, Electron_enu, boost=1)]
            boostedTau_enu = [x for x in boostedTau_enu if MuonTauOverlap(x, Muon_enu, boost=1)]
          
            

            # Cutflow counting
            if sys == "":
                if (len(Tau_enu) > 0 or len(boostedTau_enu) > 0):
                    self.cutflow_dict[" ...breakdown..> Atleast_one_Tau (Reco + ID + cleaning)"] += 1
                if (len(Muon_enu) > 0):
                    self.cutflow_dict[" ...breakdown..> Atleast_one_Muon (Reco + IDnoIso + cleaning)"] += 1
                # if (len(Electron_enu) > 0):
                #     self.cutflow_dict[
                #         " ...breakdown..> Atleast_one_Electron (Reco + IDnoIso + cleaning)"] += 1

            enoughleptonstopair = (((len(Tau_enu) + len(Muon_enu)) >= 2) or ((len(boostedTau_enu) + len(Muon_enu)) >= 2))
            if not enoughleptonstopair:
                # move to the next systematics
                # fillBranchesWithDefault(sys)
                continue

            if (sys == ""):
                if (((len(Tau_enu) + len(Muon_enu)) >= 2) or ((len(boostedTau_enu) + len(Muon_enu)) >= 2)):
                    self.cutflow_dict["Atleast_2leptons_anykind"] += 1

            pairDict = {}
            pairDict["bb"] = selfPairing(boostedTau_enu, "bb")
            pairDict["tt"] = selfPairing(Tau_enu, "tt")
            # pairDict["be"] = crossPairing(boostedTau_enu, Electron_enu, "be")
            pairDict["bm"] = crossPairing(boostedTau_enu, Muon_enu, "bm")
            # pairDict["te"] = crossPairing(Tau_enu, Electron_enu, "te")
            pairDict["tm"] = crossPairing(Tau_enu, Muon_enu, "tm")

            Keymax = max(pairDict, key=lambda x: pairDict[x][0])

            if (pairDict[Keymax][0] < 0):
                # return False
                # Electron and Muon isolation failed - move on to the next
                # systematics
                # fillBranchesWithDefault(sys)
                continue

            if (sys == ""):
                #self.cutflow_dict["Atleast_one_pair_anykind"] += 1
                if ((pairDict["bb"][0] > 0) or (pairDict["tt"][0] > 0)):
                    self.cutflow_dict[" ...breakdown..> Atleast_one_TauTau_pair"] += 1
                # if ((pairDict["be"][0] > 0) or (pairDict["te"][0] > 0)):
                #     self.cutflow_dict[" ...breakdown..> Atleast_one_TauElectron_pair"] += 1
                if ((pairDict["bm"][0] > 0) or (pairDict["tm"][0] > 0)):
                    self.cutflow_dict[" ...breakdown..> Atleast_one_TauMuon_pair"] += 1

            
            ##### Everything related to subjects anf N-subjettiness of FatJets

            has_sj1 = False
            has_sj2 = False
            sj1 = None
            sj2 = None

            if len(FatJet_enu) > 0:
                fat = FatJet_enu[0][1]
                if fat.tau1 > 0:
                    fat_tau21 = fat.tau2 / fat.tau1
                    self.out.fillBranch("fatjet_tau21", fat_tau21)
                else:
                    self.out.fillBranch("fatjet_tau21", -99.99)

                if fat.tau2 > 0:
                    fat_tau32 = fat.tau3 / fat.tau2
                    self.out.fillBranch("fatjet_tau32", fat_tau32)
                else:
                    self.out.fillBranch("fatjet_tau32", -99.99)

                if hasattr(event, "nSubJet"):
                    SubJet = Collection(event, "SubJet", "nSubJet")

                    sj1_idx = getattr(fat, "subJetIdx1", -1)
                    sj2_idx = getattr(fat, "subJetIdx2", -1)

                    if 0 <= sj1_idx < len(SubJet):
                        sj1 = SubJet[sj1_idx]
                        self.subjet1FV.SetPtEtaPhiM(sj1.pt, sj1.eta, sj1.phi, sj1.mass)
                        has_sj1 = True

                    if 0 <= sj2_idx < len(SubJet):
                        sj2 = SubJet[sj2_idx]
                        self.subjet2FV.SetPtEtaPhiM(sj2.pt, sj2.eta, sj2.phi, sj2.mass)
                        has_sj2 = True

            
            if has_sj1 and has_sj2:

                self.out.fillBranch("deltaR_subjets", self.subjet1FV.DeltaR(self.subjet2FV))
                self.out.fillBranch("deltaPhi_subjets", abs(self.subjet1FV.DeltaPhi(self.subjet2FV)))

                sj1_tau21 = sj1.tau2 / sj1.tau1 if sj1.tau1 > 0 else -99.99
                sj1_tau32 = sj1.tau3 / sj1.tau2 if sj1.tau2 > 0 else -99.99
                self.out.fillBranch("subjet1_tau21", sj1_tau21)
                self.out.fillBranch("subjet1_tau32", sj1_tau32)

                sj2_tau21 = sj2.tau2 / sj2.tau1 if sj2.tau1 > 0 else -99.99
                sj2_tau32 = sj2.tau3 / sj2.tau2 if sj2.tau2 > 0 else -99.99
                self.out.fillBranch("subjet2_tau21", sj2_tau21)
                self.out.fillBranch("subjet2_tau32", sj2_tau32)

            else:
                self.out.fillBranch("deltaR_subjets", -99.99)
                self.out.fillBranch("deltaPhi_subjets", -99.99)
                self.out.fillBranch("subjet1_tau21", -99.99)
                self.out.fillBranch("subjet1_tau32", -99.99)
                self.out.fillBranch("subjet2_tau21", -99.99)
                self.out.fillBranch("subjet2_tau32", -99.99)



            gFatJet_index = [FatJet_enu[0][0]]
            gTau_index = []
            gboostedTau_index = []
            gElectron_index = []
            gMuon_index = []
            gJet_index = []
            firstLepton = fastMTTlepton()
            secondLepton = fastMTTlepton()
            theMET = fastMTTmet(
                measuredX=getMETpt(sys) * math.cos(getMETphi(sys)),
                measuredY=getMETpt(sys) * math.sin(getMETphi(sys)),
                xx=event.PuppiMET_covXX, xy=event.PuppiMET_covXY, yy=event.PuppiMET_covYY)
            
            self.met.SetPtEtaPhiM(getMETpt(sys), 0.0, getMETphi(sys), 0.0)


            if Keymax == "bb":
                self.out.fillBranch("channel%s" % (sys), 0)
                self.out.fillBranch("boost%s" % (sys), 1)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 2)
                gboostedTau_index = [pairDict[Keymax][1], pairDict[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (sys), [boostedTau[gboostedTau_index[0]].decayMode, boostedTau[gboostedTau_index[1]].decayMode])
                
                firstLepton = fastMTTlepton(pt=gettaupt(boostedTau[gboostedTau_index[0]], sys), eta=boostedTau[gboostedTau_index[0]].eta, phi=boostedTau[gboostedTau_index[0]].phi, m=gettaumass(boostedTau[gboostedTau_index[0]], sys),leptonType='Tau', tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
            
                self.pair1FV.SetPtEtaPhiM(gettaupt(boostedTau[gboostedTau_index[0]], sys), boostedTau[gboostedTau_index[0]].eta, boostedTau[gboostedTau_index[0]].phi, gettaumass(boostedTau[gboostedTau_index[0]], sys))
                
                secondLepton = fastMTTlepton(pt=gettaupt(boostedTau[gboostedTau_index[1]], sys), eta=boostedTau[gboostedTau_index[1]].eta, phi=boostedTau[gboostedTau_index[1]].phi, m=gettaumass(boostedTau[gboostedTau_index[1]], sys), leptonType='Tau', tauDecayMode=boostedTau[gboostedTau_index[1]].decayMode)
               
                self.pair2FV.SetPtEtaPhiM(gettaupt(boostedTau[gboostedTau_index[1]], sys), boostedTau[gboostedTau_index[1]].eta, boostedTau[gboostedTau_index[1]].phi, gettaumass(boostedTau[gboostedTau_index[1]], sys))

                self.theFastMTTtool.setFirstLepton(firstLepton)
                self.theFastMTTtool.setSecondLepton(secondLepton)
                self.theFastMTTtool.setTheMET(theMET)
                higgsFV_list = self.theFastMTTtool.getFastMTTfourvector()

                self.higgsTTFV.SetPtEtaPhiM(higgsFV_list[0], higgsFV_list[1], higgsFV_list[2], higgsFV_list[3])
                self.higgsTTvisFV = self.pair1FV + self.pair2FV
                
                self.out.fillBranch("deltaPhi_met_tautau", abs(self.met.DeltaPhi(self.higgsTTvisFV)))
                self.out.fillBranch("deltaPhi_met_leadingtau", abs(self.met.DeltaPhi(self.pair1FV)))
                self.out.fillBranch("deltaPhi_met_subleadingtau", abs(self.met.DeltaPhi(self.pair2FV)))

                self.out.fillBranch("deltaR_hbb_httvis", self.higgsBBFV.DeltaR(self.higgsTTvisFV))
                self.out.fillBranch("deltaPhi_hbb_httvis", abs(self.higgsBBFV.DeltaPhi(self.higgsTTvisFV)))                
                
                self.out.fillBranch("deltaR_hbb_htt", self.higgsBBFV.DeltaR(self.higgsTTFV))
                self.out.fillBranch("deltaPhi_hbb_htt", abs(self.higgsBBFV.DeltaPhi(self.higgsTTFV)))
                
                self.out.fillBranch("deltaPhi_hbb_leadingtau", abs(self.higgsBBFV.DeltaPhi(self.pair1FV)))
                self.out.fillBranch("deltaPhi_hbb_subleadingtau", abs(self.higgsBBFV.DeltaPhi(self.pair2FV)))

                self.out.fillBranch("deltaPhi_tau1_tau2", abs(self.pair1FV.DeltaPhi(self.pair2FV)))
                self.out.fillBranch("deltaR_tau1_tau2", self.pair1FV.DeltaR(self.pair2FV))
                
                self.out.fillBranch("HTT_boosted_m", self.higgsTTFV.M())
                self.out.fillBranch("HTT_boosted_eta", self.higgsTTFV.Eta())
                self.out.fillBranch("HTT_boosted_phi", self.higgsTTFV.Phi())

                self.out.fillBranch("HTTvis_boosted_m", self.higgsTTvisFV.M())
                self.out.fillBranch("HTTvis_boosted_eta", self.higgsTTvisFV.Eta())
                self.out.fillBranch("HTTvis_boosted_phi", self.higgsTTvisFV.Phi())
    
                if has_sj1 and has_sj2:
                    self.out.fillBranch("deltaR_subjet1_leadtau", self.subjet1FV.DeltaR(self.pair1FV))
                    self.out.fillBranch("deltaPhi_subjet1_leadtau", abs(self.subjet1FV.DeltaPhi(self.pair1FV)))
                    self.out.fillBranch("deltaR_subjet1_subtau", self.subjet1FV.DeltaR(self.pair2FV))
                    self.out.fillBranch("deltaPhi_subjet1_subtau", abs(self.subjet1FV.DeltaPhi(self.pair2FV)))

                    self.out.fillBranch("deltaR_subjet2_leadtau", self.subjet2FV.DeltaR(self.pair1FV))
                    self.out.fillBranch("deltaPhi_subjet2_leadtau", abs(self.subjet2FV.DeltaPhi(self.pair1FV)))
                    self.out.fillBranch("deltaR_subjet2_subtau", self.subjet2FV.DeltaR(self.pair2FV))
                    self.out.fillBranch("deltaPhi_subjet2_subtau", abs(self.subjet2FV.DeltaPhi(self.pair2FV)))

                
                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> TT_channel (max pt pair)"] += 1
            
            
            elif Keymax == "tt":
                self.out.fillBranch("channel%s" % (sys), 0)
                self.out.fillBranch("boost%s" % (sys), 0)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 2)
                gTau_index = [pairDict[Keymax][1], pairDict[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (sys), [Tau[gTau_index[0]].decayMode, Tau[gTau_index[1]].decayMode])
                
                firstLepton = fastMTTlepton(pt=gettaupt(Tau[gTau_index[0]], sys),eta=Tau[gTau_index[0]].eta, phi=Tau[gTau_index[0]].phi, m=gettaumass(Tau[gTau_index[0]], sys), leptonType='Tau', tauDecayMode=Tau[gTau_index[0]].decayMode)
                self.pair1FV.SetPtEtaPhiM(gettaupt(Tau[gTau_index[0]], sys), Tau[gTau_index[0]].eta, Tau[gTau_index[0]].phi, gettaumass(Tau[gTau_index[0]], sys))
                
                secondLepton = fastMTTlepton(pt=gettaupt(Tau[gTau_index[1]], sys), eta=Tau[gTau_index[1]].eta, phi=Tau[gTau_index[1]].phi, m=gettaumass(Tau[gTau_index[1]], sys),
                                             leptonType='Tau', tauDecayMode=Tau[gTau_index[1]].decayMode)
                self.pair2FV.SetPtEtaPhiM(gettaupt(Tau[gTau_index[1]], sys), Tau[gTau_index[1]].eta, Tau[gTau_index[1]].phi, gettaumass(Tau[gTau_index[1]], sys))

                self.theFastMTTtool.setFirstLepton(firstLepton)
                self.theFastMTTtool.setSecondLepton(secondLepton)
                self.theFastMTTtool.setTheMET(theMET)
                higgsFV_list = self.theFastMTTtool.getFastMTTfourvector()

                self.higgsTTFV.SetPtEtaPhiM(higgsFV_list[0], higgsFV_list[1], higgsFV_list[2], higgsFV_list[3])
                self.higgsTTvisFV = self.pair1FV + self.pair2FV

                self.out.fillBranch("deltaPhi_met_tautau", abs(self.met.DeltaPhi(self.higgsTTvisFV)))
                self.out.fillBranch("deltaR_hbb_httvis", self.higgsBBFV.DeltaR(self.higgsTTvisFV))
                self.out.fillBranch("deltaPhi_hbb_httvis", abs(self.higgsBBFV.DeltaPhi(self.higgsTTvisFV)))

                self.out.fillBranch("deltaPhi_met_leadingtau", abs(self.met.DeltaPhi(self.pair1FV)))
                self.out.fillBranch("deltaPhi_met_subleadingtau", abs(self.met.DeltaPhi(self.pair2FV)))

                self.out.fillBranch("deltaR_hbb_htt", self.higgsBBFV.DeltaR(self.higgsTTFV))
                self.out.fillBranch("deltaPhi_hbb_htt", abs(self.higgsBBFV.DeltaPhi(self.higgsTTFV)))
                
                self.out.fillBranch("deltaPhi_hbb_leadingtau", abs(self.higgsBBFV.DeltaPhi(self.pair1FV)))
                self.out.fillBranch("deltaPhi_hbb_subleadingtau", abs(self.higgsBBFV.DeltaPhi(self.pair2FV)))

                self.out.fillBranch("deltaPhi_tau1_tau2", abs(self.pair1FV.DeltaPhi(self.pair2FV)))
                self.out.fillBranch("deltaR_tau1_tau2", self.pair1FV.DeltaR(self.pair2FV))

                self.out.fillBranch("HTT_HPS_m", self.higgsTTFV.M())
                self.out.fillBranch("HTT_HPS_eta", self.higgsTTFV.Eta())
                self.out.fillBranch("HTT_HPS_phi", self.higgsTTFV.Phi())

                self.out.fillBranch("HTTvis_HPS_m", self.higgsTTvisFV.M())
                self.out.fillBranch("HTTvis_HPS_eta", self.higgsTTvisFV.Eta())
                self.out.fillBranch("HTTvis_HPS_phi", self.higgsTTvisFV.Phi())
                
                if has_sj1 and has_sj2:
                    self.out.fillBranch("deltaR_subjet1_leadtau", self.subjet1FV.DeltaR(self.pair1FV))
                    self.out.fillBranch("deltaPhi_subjet1_leadtau", abs(self.subjet1FV.DeltaPhi(self.pair1FV)))
                    self.out.fillBranch("deltaR_subjet1_subtau", self.subjet1FV.DeltaR(self.pair2FV))
                    self.out.fillBranch("deltaPhi_subjet1_subtau", abs(self.subjet1FV.DeltaPhi(self.pair2FV)))

                    self.out.fillBranch("deltaR_subjet2_leadtau", self.subjet2FV.DeltaR(self.pair1FV))
                    self.out.fillBranch("deltaPhi_subjet2_leadtau", abs(self.subjet2FV.DeltaPhi(self.pair1FV)))
                    self.out.fillBranch("deltaR_subjet2_subtau", self.subjet2FV.DeltaR(self.pair2FV))
                    self.out.fillBranch("deltaPhi_subjet2_subtau", abs(self.subjet2FV.DeltaPhi(self.pair2FV)))



                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> TT_channel (max pt pair)"] += 1
            
            
            # elif Keymax == "be":
            #     self.out.fillBranch("channel%s" % (sys), 1)
            #     self.out.fillBranch("boost%s" % (sys), 1)
            #     if (sys == ""):
            #         self.out.fillBranch("nallTaus%s" % (sys), 1)
            #     gboostedTau_index = [pairDict[Keymax][1]]
            #     gElectron_index = [pairDict[Keymax][2]]
            #     if (sys == ""):
            #         self.out.fillBranch("allTaus_decayMode%s" % (sys), [boostedTau[gboostedTau_index[0]].decayMode])
                
            #     firstLepton = fastMTTlepton(pt=gettaupt(boostedTau[gboostedTau_index[0]],sys), eta=boostedTau[gboostedTau_index[0]].eta, phi=boostedTau[gboostedTau_index[0]].phi, m=gettaumass(boostedTau[gboostedTau_index[0]], sys), leptonType='Tau', tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
            #     self.pair1FV.SetPtEtaPhiM(gettaupt(boostedTau[gboostedTau_index[0]], sys), boostedTau[gboostedTau_index[0]].eta, boostedTau[gboostedTau_index[0]].phi, gettaumass(boostedTau[gboostedTau_index[0]],sys))
                
            #     secondLepton = fastMTTlepton(pt=Electron[gElectron_index[0]].pt, eta=Electron[gElectron_index[0]].eta, phi=Electron[gElectron_index[0]].phi, m=0.51100e-3, leptonType='Electron', tauDecayMode=-1)
            #     self.pair2FV.SetPtEtaPhiM(Electron[gElectron_index[0]].pt, Electron[gElectron_index[0]].eta, Electron[gElectron_index[0]].phi, 0.0)

            #     self.theFastMTTtool.setFirstLepton(firstLepton)
            #     self.theFastMTTtool.setSecondLepton(secondLepton)
            #     self.theFastMTTtool.setTheMET(theMET)
            #     higgsFV_list = self.theFastMTTtool.getFastMTTfourvector()

            #     self.higgsTTFV.SetPtEtaPhiM(higgsFV_list[0], higgsFV_list[1], higgsFV_list[2], higgsFV_list[3])

            #     self.out.fillBranch("deltaPhi_met_leadingtau", abs(self.met.DeltaPhi(self.pair1FV)))
            #     self.out.fillBranch("deltaPhi_met_leadingele", abs(self.met.DeltaPhi(self.pair2FV)))
            #     self.out.fillBranch("deltaPhi_hbb_htt", abs(self.higgsBBFV.DeltaPhi(self.higgsTTFV)))

            #     self.out.fillBranch("deltaPhi_hbb_leadingtau", abs(self.higgsBBFV.DeltaPhi(self.pair1FV)))
            #     self.out.fillBranch("deltaPhi_hbb_leadingele", abs(self.higgsBBFV.DeltaPhi(self.pair2FV)))

            #     self.out.fillBranch("deltaPhi_tau_ele", abs(self.pair1FV.DeltaPhi(self.pair2FV)))
            #     self.out.fillBranch("deltaR_tau_ele", self.pair1FV.DeltaR(self.pair2FV))


            #     self.out.fillBranch("HTT_boosted_Ele_m",self.higgsTTFV.M())
            #     self.out.fillBranch("HTT_boosted_Ele_eta", self.higgsTTFV.Eta())
            #     self.out.fillBranch("HTT_boosted_Ele_phi", self.higgsTTFV.Phi())

            #     if has_sj1 and has_sj2:
            #         self.out.fillBranch("deltaR_subjet1_leadtau", self.subjet1FV.DeltaR(self.pair1FV))
            #         self.out.fillBranch("deltaPhi_subjet1_leadtau", abs(self.subjet1FV.DeltaPhi(self.pair1FV)))
            #         self.out.fillBranch("deltaR_subjet1_ele", self.subjet1FV.DeltaR(self.pair2FV))
            #         self.out.fillBranch("deltaPhi_subjet1_ele", abs(self.subjet1FV.DeltaPhi(self.pair2FV)))
                    
            #         self.out.fillBranch("deltaR_subjet2_leadtau", self.subjet2FV.DeltaR(self.pair1FV))
            #         self.out.fillBranch("deltaPhi_subjet2_leadtau", abs(self.subjet2FV.DeltaPhi(self.pair1FV)))
            #         self.out.fillBranch("deltaR_subjet2_ele", self.subjet2FV.DeltaR(self.pair2FV))
            #         self.out.fillBranch("deltaPhi_subjet2_ele", abs(self.subjet2FV.DeltaPhi(self.pair2FV)))

            

            #     if (sys == ""):
            #         self.cutflow_dict[" ...breakdown..> ET_channel (max pt pair)"] += 1
            
            
            # elif Keymax == "te":
            #     self.out.fillBranch("channel%s" % (sys), 1)
            #     self.out.fillBranch("boost%s" % (sys), 0)
            #     if (sys == ""):
            #         self.out.fillBranch("nallTaus%s" % (sys), 1)
            #     gTau_index = [pairDict[Keymax][1]]
            #     gElectron_index = [pairDict[Keymax][2]]
            #     if (sys == ""):
            #         self.out.fillBranch("allTaus_decayMode%s" % (sys), [Tau[gTau_index[0]].decayMode])
                
            #     firstLepton = fastMTTlepton(pt=gettaupt(Tau[gTau_index[0]], sys), eta=Tau[gTau_index[0]].eta, phi=Tau[gTau_index[0]].phi, m=gettaumass(Tau[gTau_index[0]], sys), leptonType='Tau', tauDecayMode=Tau[gTau_index[0]].decayMode)
            #     self.pair1FV.SetPtEtaPhiM(gettaupt(Tau[gTau_index[0]], sys), Tau[gTau_index[0]].eta, Tau[gTau_index[0]].phi, gettaumass(Tau[gTau_index[0]], sys))
                
            #     secondLepton = fastMTTlepton(pt=Electron[gElectron_index[0]].pt, eta=Electron[gElectron_index[0]].eta, phi=Electron[gElectron_index[0]].phi, m=0.51100e-3, leptonType='Electron', tauDecayMode=-1)
            #     self.pair2FV.SetPtEtaPhiM(Electron[gElectron_index[0]].pt, Electron[gElectron_index[0]].eta, Electron[gElectron_index[0]].phi, 0.0)

               
            #     self.theFastMTTtool.setFirstLepton(firstLepton)
            #     self.theFastMTTtool.setSecondLepton(secondLepton)
            #     self.theFastMTTtool.setTheMET(theMET)
            #     higgsFV_list = self.theFastMTTtool.getFastMTTfourvector()

            #     self.higgsTTFV.SetPtEtaPhiM(higgsFV_list[0], higgsFV_list[1], higgsFV_list[2], higgsFV_list[3])

            #     self.out.fillBranch("deltaPhi_met_leadingtau", abs(self.met.DeltaPhi(self.pair1FV)))
            #     self.out.fillBranch("deltaPhi_met_leadingele", abs(self.met.DeltaPhi(self.pair2FV)))
            #     self.out.fillBranch("deltaPhi_hbb_htt", abs(self.higgsBBFV.DeltaPhi(self.higgsTTFV)))

            #     self.out.fillBranch("deltaPhi_hbb_leadingtau", abs(self.higgsBBFV.DeltaPhi(self.pair1FV)))
            #     self.out.fillBranch("deltaPhi_hbb_leadingele", abs(self.higgsBBFV.DeltaPhi(self.pair2FV)))
                
            #     self.out.fillBranch("deltaPhi_tau_ele", abs(self.pair1FV.DeltaPhi(self.pair2FV)))
            #     self.out.fillBranch("deltaR_tau_ele", self.pair1FV.DeltaR(self.pair2FV))

                
            #     self.out.fillBranch("HTT_HPS_Ele_m", self.higgsTTFV.M())
            #     self.out.fillBranch("HTT_HPS_Ele_eta", self.higgsTTFV.Eta())
            #     self.out.fillBranch("HTT_HPS_Ele_phi", self.higgsTTFV.Phi())

            #     if has_sj1 and has_sj2:
            #         self.out.fillBranch("deltaR_subjet1_leadtau", self.subjet1FV.DeltaR(self.pair1FV))
            #         self.out.fillBranch("deltaPhi_subjet1_leadtau", abs(self.subjet1FV.DeltaPhi(self.pair1FV)))
            #         self.out.fillBranch("deltaR_subjet1_ele", self.subjet1FV.DeltaR(self.pair2FV))
            #         self.out.fillBranch("deltaPhi_subjet1_ele", abs(self.subjet1FV.DeltaPhi(self.pair2FV)))
                    
            #         self.out.fillBranch("deltaR_subjet2_leadtau", self.subjet2FV.DeltaR(self.pair1FV))
            #         self.out.fillBranch("deltaPhi_subjet2_leadtau", abs(self.subjet2FV.DeltaPhi(self.pair1FV)))
            #         self.out.fillBranch("deltaR_subjet2_ele", self.subjet2FV.DeltaR(self.pair2FV))
            #         self.out.fillBranch("deltaPhi_subjet2_ele", abs(self.subjet2FV.DeltaPhi(self.pair2FV)))

            #     if (sys == ""):
            #         self.cutflow_dict[" ...breakdown..> ET_channel (max pt pair)"] += 1
            
            
            elif Keymax == "bm":
                self.out.fillBranch("channel%s" % (sys), 2)
                self.out.fillBranch("boost%s" % (sys), 1)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 1)
                gboostedTau_index = [pairDict[Keymax][1]]
                gMuon_index = [pairDict[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (
                        sys), [boostedTau[gboostedTau_index[0]].decayMode])
                
                firstLepton = fastMTTlepton(pt=gettaupt(boostedTau[gboostedTau_index[0]], sys), eta=boostedTau[gboostedTau_index[0]].eta, phi=boostedTau[gboostedTau_index[0]].phi, m=gettaumass(boostedTau[gboostedTau_index[0]], sys), leptonType='Tau', tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)

                self.pair1FV.SetPtEtaPhiM(gettaupt(boostedTau[gboostedTau_index[0]], sys), boostedTau[gboostedTau_index[0]].eta, boostedTau[gboostedTau_index[0]].phi, gettaumass(boostedTau[gboostedTau_index[0]],sys))
                
                secondLepton = fastMTTlepton(pt=Muon[gMuon_index[0]].pt, eta=Muon[gMuon_index[0]].eta, phi=Muon[gMuon_index[0]].phi, m=Muon[gMuon_index[0]].mass, leptonType='Muon', tauDecayMode=-1)
                self.pair2FV.SetPtEtaPhiM(Muon[gMuon_index[0]].pt, Muon[gMuon_index[0]].eta,  Muon[gMuon_index[0]].phi, Muon[gMuon_index[0]].mass)

                self.theFastMTTtool.setFirstLepton(firstLepton)
                self.theFastMTTtool.setSecondLepton(secondLepton)
                self.theFastMTTtool.setTheMET(theMET)
                higgsFV_list = self.theFastMTTtool.getFastMTTfourvector()

                self.higgsTTFV.SetPtEtaPhiM(higgsFV_list[0], higgsFV_list[1], higgsFV_list[2], higgsFV_list[3])

                self.out.fillBranch("deltaPhi_met_leadingtau", abs(self.met.DeltaPhi(self.pair1FV)))
                self.out.fillBranch("deltaPhi_met_leadingmu", abs(self.met.DeltaPhi(self.pair2FV)))
                self.out.fillBranch("deltaPhi_hbb_htt", abs(self.higgsBBFV.DeltaPhi(self.higgsTTFV)))
                self.out.fillBranch("deltaPhi_hbb_leadingtau", abs(self.higgsBBFV.DeltaPhi(self.pair1FV)))
                self.out.fillBranch("deltaPhi_hbb_leadingmu", abs(self.higgsBBFV.DeltaPhi(self.pair2FV)))

                self.out.fillBranch("deltaPhi_tau_mu", abs(self.pair1FV.DeltaPhi(self.pair2FV)))
                self.out.fillBranch("deltaR_tau_mu", self.pair1FV.DeltaR(self.pair2FV))

                self.out.fillBranch("HTT_boosted_Mu_m", self.higgsTTFV.M())
                self.out.fillBranch("HTT_boosted_Mu_eta", self.higgsTTFV.Eta())
                self.out.fillBranch("HTT_boosted_Mu_phi", self.higgsTTFV.Phi())

                if has_sj1 and has_sj2:
                    self.out.fillBranch("deltaR_subjet1_leadtau", self.subjet1FV.DeltaR(self.pair1FV))
                    self.out.fillBranch("deltaPhi_subjet1_leadtau", abs(self.subjet1FV.DeltaPhi(self.pair1FV)))
                    self.out.fillBranch("deltaR_subjet1_mu", self.subjet1FV.DeltaR(self.pair2FV))
                    self.out.fillBranch("deltaPhi_subjet1_mu", abs(self.subjet1FV.DeltaPhi(self.pair2FV)))
                    
                    self.out.fillBranch("deltaR_subjet2_leadtau", self.subjet2FV.DeltaR(self.pair1FV))
                    self.out.fillBranch("deltaPhi_subjet2_leadtau", abs(self.subjet2FV.DeltaPhi(self.pair1FV)))
                    self.out.fillBranch("deltaR_subjet2_mu", self.subjet2FV.DeltaR(self.pair2FV))
                    self.out.fillBranch("deltaPhi_subjet2_mu", abs(self.subjet2FV.DeltaPhi(self.pair2FV)))



                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> MT_channel (max pt pair)"] += 1

            elif Keymax == "tm":
                self.out.fillBranch("channel%s" % (sys), 2)
                self.out.fillBranch("boost%s" % (sys), 0)
                if (sys == ""):
                    self.out.fillBranch("nallTaus%s" % (sys), 1)
                gTau_index = [pairDict[Keymax][1]]
                gMuon_index = [pairDict[Keymax][2]]
                if (sys == ""):
                    self.out.fillBranch("allTaus_decayMode%s" % (sys), [
                                        Tau[gTau_index[0]].decayMode])
                
                firstLepton = fastMTTlepton(pt=gettaupt(Tau[gTau_index[0]],  sys),  eta=Tau[gTau_index[0]].eta,  phi=Tau[gTau_index[0]].phi, m=gettaumass(Tau[gTau_index[0]], sys), leptonType='Tau', tauDecayMode=Tau[gTau_index[0]].decayMode)
                
                self.pair1FV.SetPtEtaPhiM(gettaupt(Tau[gTau_index[0]], sys), Tau[gTau_index[0]].eta, Tau[gTau_index[0]].phi, gettaumass(Tau[gTau_index[0]], sys))
                
                secondLepton = fastMTTlepton(pt=Muon[gMuon_index[0]].pt, eta=Muon[gMuon_index[0]].eta, phi=Muon[gMuon_index[0]].phi, m=Muon[gMuon_index[0]].mass, leptonType='Muon', tauDecayMode=-1)

                self.pair2FV.SetPtEtaPhiM(Muon[gMuon_index[0]].pt, Muon[gMuon_index[0]].eta, Muon[gMuon_index[0]].phi, Muon[gMuon_index[0]].mass)

               
                self.theFastMTTtool.setFirstLepton(firstLepton)
                self.theFastMTTtool.setSecondLepton(secondLepton)
                self.theFastMTTtool.setTheMET(theMET)
                higgsFV_list = self.theFastMTTtool.getFastMTTfourvector()

                self.higgsTTFV.SetPtEtaPhiM(higgsFV_list[0], higgsFV_list[1], higgsFV_list[2], higgsFV_list[3])

                self.out.fillBranch("deltaPhi_met_leadingtau", abs(self.met.DeltaPhi(self.pair1FV)))
                self.out.fillBranch("deltaPhi_met_leadingmu", abs(self.met.DeltaPhi(self.pair2FV)))
                self.out.fillBranch("deltaPhi_hbb_htt", abs(self.higgsBBFV.DeltaPhi(self.higgsTTFV)))
                self.out.fillBranch("deltaPhi_hbb_leadingtau", abs(self.higgsBBFV.DeltaPhi(self.pair1FV)))
                self.out.fillBranch("deltaPhi_hbb_leadingmu", abs(self.higgsBBFV.DeltaPhi(self.pair2FV)))

                self.out.fillBranch("deltaPhi_tau_mu", abs(self.pair1FV.DeltaPhi(self.pair2FV)))
                self.out.fillBranch("deltaR_tau_mu", self.pair1FV.DeltaR(self.pair2FV))

                self.out.fillBranch("HTT_HPS_Mu_m", self.higgsTTFV.M())
                self.out.fillBranch("HTT_HPS_Mu_eta", self.higgsTTFV.Eta())
                self.out.fillBranch("HTT_HPS_Mu_phi", self.higgsTTFV.Phi())

                if has_sj1 and has_sj2:
                    self.out.fillBranch("deltaR_subjet1_leadtau", self.subjet1FV.DeltaR(self.pair1FV))
                    self.out.fillBranch("deltaPhi_subjet1_leadtau", abs(self.subjet1FV.DeltaPhi(self.pair1FV)))
                    self.out.fillBranch("deltaR_subjet1_mu", self.subjet1FV.DeltaR(self.pair2FV))
                    self.out.fillBranch("deltaPhi_subjet1_mu", abs(self.subjet1FV.DeltaPhi(self.pair2FV)))
                    
                    self.out.fillBranch("deltaR_subjet2_leadtau", self.subjet2FV.DeltaR(self.pair1FV))
                    self.out.fillBranch("deltaPhi_subjet2_leadtau", abs(self.subjet2FV.DeltaPhi(self.pair1FV)))
                    self.out.fillBranch("deltaR_subjet2_mu", self.subjet2FV.DeltaR(self.pair2FV))
                    self.out.fillBranch("deltaPhi_subjet2_mu", abs(self.subjet2FV.DeltaPhi(self.pair2FV)))
                
                
                if (sys == ""):
                    self.cutflow_dict[" ...breakdown..> MT_channel (max pt pair)"] += 1
            
                       
            # pass_quality_delcuts = (
            #     ((abs(
            #         self.higgsBBFV.DeltaPhi(
            #             self.met))) > 1) and (
            #         (self.pair1FV.DeltaR(
            #             self.pair2FV)) > 0) and (
            #         (self.pair1FV.DeltaR(
            #             self.pair2FV)) < 1.5))

            # if not pass_quality_delcuts:
            #     # move on to the next systematic
            #     fillBranchesWithDefault(sys)
            #     continue
            # # delete the temp variable
            # del pass_quality_delcuts

            # # Fill cut flow
            # if (sys == ""):
            #     self.cutflow_dict["DeltaR_LL<1.5 and abs(Hbb_met_phi) > 1 cut"] += 1

            self.out.fillBranch("Hbb_met_phi%s" %
                                (sys), self.higgsBBFV.DeltaPhi(self.met))

            # self.theFastMTTtool.setFirstLepton(firstLepton)
            # self.theFastMTTtool.setSecondLepton(secondLepton)
            # self.theFastMTTtool.setTheMET(theMET)
            # higgsFV_list = self.theFastMTTtool.getFastMTTfourvector()
            
            # self.higgsTTFV.SetPtEtaPhiM(higgsFV_list[0], higgsFV_list[1], higgsFV_list[2], higgsFV_list[3])
            self.RadionFV = self.higgsTTFV + self.higgsBBFV

            # self.higgsTTvisFV = self.pair1FV + self.pair2FV
            self.RadionvisFV = self.higgsTTvisFV + self.higgsBBFV

            pt_balance_hbb_htt_abs = abs(self.higgsBBFV.Pt() - self.higgsTTFV.Pt()) / (self.higgsBBFV.Pt() + self.higgsTTFV.Pt())
            self.out.fillBranch("pt_balance_hbb_htt_abs", pt_balance_hbb_htt_abs)
            pt_balance_hbb_htt_signed = (self.higgsBBFV.Pt() - self.higgsTTFV.Pt()) / (self.higgsBBFV.Pt() + self.higgsTTFV.Pt())
            self.out.fillBranch("pt_balance_hbb_htt_signed", pt_balance_hbb_htt_signed)

            # if (self.higgsTTvisFV.M() <= 20):
            #     fillBranchesWithDefault(sys)
            #     # move to next systematic
            #     continue

            # # Fill cut flow
            # if (sys == ""):
            #     self.cutflow_dict["Visible Mass HTT > 20 cut"] += 1

            Jet_enu = [x for x in Jet_enu if removeOverlapOfAK4WithLightHeavyLeptons(
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

            # if len(Jet_enu) == 0:
            #     fillBranchesWithDefault(sys)
            #     continue

            Jet_enu = sorted(Jet_enu, key=lambda x: getjetpt(x[1], sys), reverse=True)

            deltaR_boosted_HPS_postclean = []
            for b in boostedTau_enu:
                for h in Tau_enu:
                    bfv = ROOT.TLorentzVector()
                    hfv = ROOT.TLorentzVector()
                    bfv.SetPtEtaPhiM(
                        gettaupt(b[1], sys), b[1].eta, b[1].phi, gettaumass(b[1], sys)
                    )
                    hfv.SetPtEtaPhiM(
                        gettaupt(h[1], sys), h[1].eta, h[1].phi, gettaumass(h[1], sys)
                    )
                    deltaR_boosted_HPS_postclean.append(bfv.DeltaR(hfv))

            self.out.fillBranch("nDeltaR_boosted_HPS_postclean", len(deltaR_boosted_HPS_postclean))
            self.out.fillBranch("deltaR_boosted_HPS_postclean", deltaR_boosted_HPS_postclean)
            
            gJet_index = [x[0] for x in Jet_enu]
            Jet_enu_Loose = [
               x for x in Jet_enu if x[1].btagUParTAK4B >= self.LooseJet]
            gJet_Looseindex = [x[0] for x in Jet_enu_Loose]
            Jet_enu_Medium = [
                x for x in Jet_enu if x[1].btagUParTAK4B >= self.MediumJet]
            gJet_Mediumindex = [x[0] for x in Jet_enu_Medium]
            Jet_enu_Tight = [
                x for x in Jet_enu if x[1].btagUParTAK4B >= self.TightJet]
            gJet_Tightindex = [x[0] for x in Jet_enu_Tight]

            if len(Jet_enu) > 0:
                leadJet = Jet_enu[0][1]
                self.jetLeadFV.SetPtEtaPhiM(getjetpt(leadJet, sys), leadJet.eta, leadJet.phi, getjetmass(leadJet, sys))
            
            self.out.fillBranch("deltaR_hbb_ak4lead", self.higgsBBFV.DeltaR(self.jetLeadFV) if len(Jet_enu) > 0 else -99.99)
            self.out.fillBranch("deltaPhi_hbb_ak4lead", abs(self.higgsBBFV.DeltaPhi(self.jetLeadFV)) if len(Jet_enu) > 0 else -99.99 )
            self.out.fillBranch("deltaPhi_met_ak4lead", abs(self.met.DeltaPhi(self.jetLeadFV)) if len(Jet_enu) > 0 else -99.99)
            
            if Keymax == "bb":
                self.out.fillBranch("deltaR_ak4_leadtau", self.jetLeadFV.DeltaR(self.pair1FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_leadtau", abs(self.jetLeadFV.DeltaPhi(self.pair1FV)) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaR_ak4_subtau", self.jetLeadFV.DeltaR(self.pair2FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_subtau", abs(self.jetLeadFV.DeltaPhi(self.pair2FV)) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaR_httvis_ak4lead", self.higgsTTvisFV.DeltaR(self.jetLeadFV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_httvis_ak4lead", abs(self.higgsTTvisFV.DeltaPhi(self.jetLeadFV)) if len(Jet_enu) > 0 else -99.99)
            elif Keymax == "tt":
                self.out.fillBranch("deltaR_ak4_leadtau", self.jetLeadFV.DeltaR(self.pair1FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_leadtau", abs(self.jetLeadFV.DeltaPhi(self.pair1FV)) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaR_ak4_subtau", self.jetLeadFV.DeltaR(self.pair2FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_subtau", abs(self.jetLeadFV.DeltaPhi(self.pair2FV)) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaR_httvis_ak4lead", self.higgsTTvisFV.DeltaR(self.jetLeadFV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_httvis_ak4lead", abs(self.higgsTTvisFV.DeltaPhi(self.jetLeadFV)) if len(Jet_enu) > 0 else -99.99)
            elif Keymax == "be":
                self.out.fillBranch("deltaR_ak4_leadtau", self.jetLeadFV.DeltaR(self.pair1FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_leadtau", abs(self.jetLeadFV.DeltaPhi(self.pair1FV)) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaR_ak4_ele", self.jetLeadFV.DeltaR(self.pair2FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_ele", abs(self.jetLeadFV.DeltaPhi(self.pair2FV)) if len(Jet_enu) > 0 else -99.99)
            elif Keymax == "te":
                self.out.fillBranch("deltaR_ak4_leadtau", self.jetLeadFV.DeltaR(self.pair1FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_leadtau", abs(self.jetLeadFV.DeltaPhi(self.pair1FV)) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaR_ak4_ele", self.jetLeadFV.DeltaR(self.pair2FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_ele", abs(self.jetLeadFV.DeltaPhi(self.pair2FV)) if len(Jet_enu) > 0 else -99.99)
            elif Keymax == "bm":
                self.out.fillBranch("deltaR_ak4_leadtau", self.jetLeadFV.DeltaR(self.pair1FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_leadtau", abs(self.jetLeadFV.DeltaPhi(self.pair1FV)) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaR_ak4_mu", self.jetLeadFV.DeltaR(self.pair2FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_mu", abs(self.jetLeadFV.DeltaPhi(self.pair2FV)) if len(Jet_enu) > 0 else -99.99)
            elif Keymax == "tm":
                self.out.fillBranch("deltaR_ak4_leadtau", self.jetLeadFV.DeltaR(self.pair1FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_leadtau", abs(self.jetLeadFV.DeltaPhi(self.pair1FV)) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaR_ak4_mu", self.jetLeadFV.DeltaR(self.pair2FV) if len(Jet_enu) > 0 else -99.99)
                self.out.fillBranch("deltaPhi_ak4_mu", abs(self.jetLeadFV.DeltaPhi(self.pair2FV)) if len(Jet_enu) > 0 else -99.99)


            tau_logit_list = []
            for idx in gTau_index:
                tau_obj = Tau[idx]
                score = tau_obj.rawDeepTau2018v2p5VSjet
                tau_logit_list.append(logit(score))

            self.out.fillBranch("Tau_rawDeepTauVSjet_logit", tau_logit_list)


            boosted_logit_list = []
            for idx in gboostedTau_index:
                btau_obj = boostedTau[idx]
                score = btau_obj.rawBoostedDeepTauRunIIv2p0VSjet
                boosted_logit_list.append(logit(score))

            self.out.fillBranch("boostedTau_rawDeepTauVSjet_logit", boosted_logit_list)


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
                # self.out.fillBranch("softdropmassnom%s" %
                #                     (sys), FatJet_enu[0][1].msoftdrop_nom)
                self.out.fillBranch("softdropmass%s" %
                                    (sys), FatJet_enu[0][1].msoftdrop)
                self.out.fillBranch("pnetmass%s" % (
                    sys), FatJet_enu[0][1].particleNetLegacy_mass)
                # self.out.fillBranch("globalparT3mass%s" %(
                #     sys), FatJet_enu[0][1].globalParT3Xbb_mass)
                

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
            self.out.fillBranch("ngood_LooseJets%s"%(sys),len(gJet_Looseindex))
            self.out.fillBranch("ngood_MediumJets%s" %
                                (sys), len(gJet_Mediumindex))
            self.out.fillBranch("ngood_TightJets%s" %
                                (sys), len(gJet_Tightindex))

            # Fill cut flow
            if (sys == ""):
                if ((len(gJet_Mediumindex) == 0)):
                    self.cutflow_dict["Medium AK4 b-tag veto (events with 0 medium b-tagged jets)"] += 1

            self.out.fillBranch("index_gboostedTaus%s" %
                                (sys), gboostedTau_index)
            self.out.fillBranch("index_gTaus%s" % (sys), gTau_index)
            self.out.fillBranch("index_gElectrons%s" % (sys), gElectron_index)
            self.out.fillBranch("index_gMuons%s" % (sys), gMuon_index)
            self.out.fillBranch("index_gFatJets%s" % (sys), gFatJet_index)
            self.out.fillBranch("index_gJets%s" % (sys), gJet_index)
            self.out.fillBranch("index_gLooseJets%s"%(sys),gJet_Looseindex)
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
            passAsingleSystematic += 1

            del Tau_enu, boostedTau_enu, Muon_enu, Jet_enu # Electron_enu,
            del firstLepton, secondLepton, theMET





        gp3_masses = []
        for fj in FatJet:
            if hasattr(fj, "globalParT3_massCorrX2p"):
                raw_mass = fj.mass * (1.0 - fj.rawFactor)
                gp3_masses.append(fj.globalParT3_massCorrX2p * raw_mass)
            else:
                gp3_masses.append(-999.0)
        self.out.fillBranch("FatJet_globalParT3Xbb_mass", gp3_masses)
        
        if passAsingleSystematic > 0:
            if (nominal_bool == 1):
                self.out.fillBranch("eventnominal", 1)
                self.cutflow_dict["Final surviving events"] += 1
            else:
                self.out.fillBranch("eventnominal", 0)
            return True
        else:
            return False


            
def call_postpoc():
    inputFile = args.inputFile
    outputFile = args.outputFile
    isMC = args.isMC
    year = args.year

    filename = os.path.basename(inputFile).replace('nanoPostProc_helper-', '').rsplit('.',1)[0]
    
    print(f"Processing file {inputFile}")
    print(f"Output will be {outputFile}")

    if not isMC:
        print(("This is a ", args.year, " Data file: ", filename))

        # def tesModule():
        #     return TauEnergyScaleForHPSandBoosted(args.year, True, tauID_wp='Loose', ele_wp='VVLoose')
        def mainModule():
            return cutsAndcategories(filename, args.year, True, args.runNominal, args.cutflowDir)
        
    else:
        print(("This is a ", args.year, " MC file: ", filename))

        # def tesModule():
        #     return TauEnergyScaleForHPSandBoosted(args.year, False, tauID_wp='Loose', ele_wp='VVLoose')
        
        def mainModule():
            return cutsAndcategories(filename, args.year, False, args.runNominal, args.cutflowDir)
              
    module_inst = mainModule()   
    
    p = PostProcessor(
        outputDir=os.path.dirname(outputFile),
        inputFiles=[inputFile],
        cut=preselection,
        branchsel=None,
        modules=[#tesModule(),
                jetId("jetid.json", jetType="AK4PUPPI"),
                fatJetId("jetid.json", jetType="AK8PUPPI"),
                jetVMAP("jetvetomaps.json",
                           corrName="Summer24Prompt24_RunBCDEFGHI_V1",
                           veto_map_name="jetvetomap"),
                fatJetVMAP("jetvetomaps.json",
                           corrName="Summer24Prompt24_RunBCDEFGHI_V1",
                           veto_map_name="jetvetomap"),
                module_inst,
                 XSWeightOnly(
                filename=os.path.basename(inputFile),
                year=year,
                isData=not isMC,
            )], 
        postfix="",
        # maxEntries=2000,
        noOut=False,
        outputbranchsel="Datadrop.txt",
        jsonInput=None if isMC else "GoldenJSON_2024.json"
        )
    
    p.run()

    produced = os.path.join(os.path.dirname(args.outputFile), os.path.basename(args.inputFile))
    if os.path.exists(produced) and produced != args.outputFile:
        shutil.move(produced, args.outputFile)
    elif not os.path.exists(args.outputFile):
        raise RuntimeError(f"[ERROR] Expected output {args.outputFile} not produced!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Single-file postprocessing driver")
    parser.add_argument("--inputFile", required=True, help="Path to input NanoAOD file")
    parser.add_argument("--outputFile", required=True, help="Name of the ROOT output file")
    parser.add_argument("--year", required=True, choices=["2024"])
    parser.add_argument("--isMC", action="store_true", help="Set this flag if processing MC")
    parser.add_argument("--cutflowDir", required=True, help="Directory to save cutflow JSON")
    parser.add_argument("--runNominal", action="store_true", help="Disable systematics")
    args = parser.parse_args()

    start_time = time.time()

    met_selection = ["PuppiMET_pt > 120"]
    Tau_selection = ["nTau > 0 || nboostedTau > 0"]


    trigger_2024 = [
        "HLT_IsoMu24", 
        "HLT_Mu50",
        "HLT_CascadeMu100",
        "HLT_HighPtTkMu100"
    ]

    if args.year == "2024":
        preselection = (
        "(" + "&&".join(met_selection) + ")" + "&&(" + "||".join(trigger_2024) + ")" + "&&(" + "&&".join(Tau_selection) + ")"
        )
    
    print(f"Applying pre-selection cuts:",{preselection})

    call_postpoc()
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(("Elapsed Time: {:.2f} seconds".format(elapsed_time)))