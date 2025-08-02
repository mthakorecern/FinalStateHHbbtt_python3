import math
import argparse
import multiprocessing as np
from re import search
import numpy as np
import glob
from FinalStateHHbbtt.fastMTTPython.fastMTTtool import *
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class cutsAndcategories(Module):
    def __init__(self, year, isData, eleID, muID, tauID, btauID):
        self.year = year
        self.isMC = not isData
        self.isData = isData

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

        # Electron, Muon, Tau, boostedTau IDs defined
        self.eleID = int(eleID)
        self.muID = int(muID)
        self.tauID = int(tauID)
        self.btauID = float(btauID)

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
        self.out = wrappedOutputTree

        self.out.branch("channel", "I")
        self.out.branch("boost", "I")

        self.out.branch("HTT_m", "F")
        self.out.branch("HTT_eta", "F")
        self.out.branch("HTT_phi", "F")
        self.out.branch("HTT_pt", "F")
        self.out.branch("HTTvis_m", "F")
        self.out.branch("HTTvis_eta", "F")
        self.out.branch("HTTvis_phi", "F")
        self.out.branch("HTTvis_pt", "F")
        self.out.branch("HTTvis_deltaR", "F")

        self.out.branch("X_m", "F")
        self.out.branch("X_eta", "F")
        self.out.branch("X_phi", "F")
        self.out.branch("X_pt", "F")
        self.out.branch("Xvis_m", "F")
        self.out.branch("Xvis_eta", "F")
        self.out.branch("Xvis_phi", "F")
        self.out.branch("Xvis_pt", "F")

        self.out.branch("nallTaus", "I")
        self.out.branch("ngood_Taus", "I")
        self.out.branch("ngood_boostedTaus", "I")
        self.out.branch("ngood_Electrons", "I")
        self.out.branch("ngood_Muons", "I")
        self.out.branch("ngood_FatJets", "I")
        self.out.branch("ngood_Jets", "I")
        self.out.branch("ngood_LooseJets", "I")
        self.out.branch("ngood_MediumJets", "I")
        self.out.branch("ngood_TightJets", "I")

        self.out.branch("allTaus_pt", "F", lenVar="nallTaus")
        self.out.branch("allTaus_eta", "F", lenVar="nallTaus")
        self.out.branch("allTaus_phi", "F", lenVar="nallTaus")
        self.out.branch("allTaus_mass", "F", lenVar="nallTaus")
        self.out.branch("allTaus_decayMode", "F", lenVar="nallTaus")

        self.out.branch("index_gElectrons", "I", lenVar="ngood_Electrons")
        self.out.branch("index_gMuons", "I", lenVar="ngood_Muons")
        self.out.branch("index_gTaus", "I", lenVar="ngood_Taus")
        self.out.branch("index_gboostedTaus", "I", lenVar="ngood_boostedTaus")
        self.out.branch("index_gFatJets", "I", lenVar="ngood_FatJets")
        self.out.branch("index_gJets", "I", lenVar="ngood_Jets")
        self.out.branch("index_gLooseJets", "I", lenVar="ngood_LooseJets")
        self.out.branch("index_gMediumJets", "I", lenVar="ngood_MediumJets")
        self.out.branch("index_gTightJets", "I", lenVar="ngood_TightJets")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        def applyPOGselectionToAK4(ak4Object_enu):
            if (((ak4Object_enu[1].pt_nom > 30) and (ak4Object_enu[1].pt_nom < 50)) and (
                    abs(ak4Object_enu[1].eta) < 2.5)):
                if ((ak4Object_enu[1].jetId > 1) and (
                        ak4Object_enu[1].puId >= 4)):
                    return True
            if ((ak4Object_enu[1].pt_nom >= 50) and (
                    abs(ak4Object_enu[1].eta) < 2.5)):
                if (ak4Object_enu[1].jetId > 1):
                    return True
            return False

        def pass_cuts_EleID(electronObject_enu):
            for cutnr in range(0, 10):
                if cutnr == 7:
                    continue
                if (electronObject_enu[1].vidNestedWPBitmap >> (
                        cutnr * 3) & 0x7) < self.eleID:
                    return False
            return True

        # FatJet and Jet overlap, separation > 1.2 (0.8 + 0.4)
        def JetFatJetOverlap(jetObject_enu):
            self.jetFV.SetPtEtaPhiM(
                jetObject_enu[1].pt_nom,
                jetObject_enu[1].eta,
                jetObject_enu[1].phi,
                jetObject_enu[1].mass_nom)
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
        def FatJetTauOverlap(tauObject_enu):
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
        def ElectronTauOverlap(tauObject_enu, elecoll_enu):
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
        def MuonTauOverlap(tauObject_enu, mucoll_enu):
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
                Mu_coll):
            self.jetFV.SetPtEtaPhiM(
                ak4Object_enu[1].pt_nom,
                ak4Object_enu[1].eta,
                ak4Object_enu[1].phi,
                ak4Object_enu[1].mass_nom)

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

                    if (tag == "be"):
                        if not ElectronIsolationCut(
                                col1[i][1], col2[j][1], tag):
                            continue

                    elif (tag == "bm"):
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
            isTau = ""
            isolationCut = 0.0

            if abs(ele.eta) <= 1.479:
                # isolationCut = 0.175
                isolationCut = 0.198 + (0.506 / ele.pt)
            elif (abs(ele.eta) > 1.479) and (abs(ele.eta) <= 2.5):
                # isolationCut = 0.159
                isolationCut = 0.203 + (0.963 / ele.pt)
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

        if (event.METcorrected_pt < 180):
            return False

        FatJet = Collection(event, 'FatJet', 'nFatJet')
        Jet = Collection(event, 'Jet', 'nJet')
        Tau = Collection(event, 'Tau', 'nTau')
        boostedTau = Collection(event, 'boostedTau', 'nboostedTau')
        Electron = Collection(event, 'Electron', 'nElectron')
        Muon = Collection(event, 'Muon', 'nMuon')

        Jet_enu = list(filter(applyPOGselectionToAK4, enumerate(Jet)))

        Tau_enu = [x for x in enumerate(Tau) if (x[1].pt > 20) and (abs(x[1].eta) < 2.5) and (abs(
            x[1].dz) < 0.2) and (x[1].idDecayModeNewDMs) and (x[1].idDeepTau2018v2p5VSjet & self.tauID == self.tauID)]

        boostedTau_enu = [x for x in enumerate(boostedTau) if (x[1].pt > 20) and (
            abs(x[1].eta) < 2.5) and (x[1].rawDeepTau2018v2p7VSjet >= self.btauID)]

        FatJet_enu = [x for x in enumerate(FatJet) if (
            x[1].pt_nom > 250) and (abs(x[1].eta) < 2.5) and (x[1].jetId > 1)]
        if (len(FatJet_enu) == 0):
            return False

        self.higgsBBFV.SetPtEtaPhiM(
            FatJet_enu[0][1].pt_nom,
            FatJet_enu[0][1].eta,
            FatJet_enu[0][1].phi,
            FatJet_enu[0][1].mass_nom)

        Electron_enu = [x for x in enumerate(Electron) if x[1].pt > 10]
        Electron_enu = list(filter(pass_cuts_EleID, Electron_enu))

        if self.muID == 2:
            Muon_enu = [x for x in enumerate(Muon) if x[1].pt > 10 and (
                abs(x[1].eta) < 2.5) and x[1].looseId]
        elif self.muID == 3:
            Muon_enu = [x for x in enumerate(Muon) if x[1].pt > 10 and (
                abs(x[1].eta) < 2.5) and x[1].mediumId]
        elif self.muID == 4:
            Muon_enu = [x for x in enumerate(Muon) if x[1].pt > 10 and (
                abs(x[1].eta) < 2.5) and x[1].tightId]

        # Object cleaning procedures
        Jet_enu = list(filter(JetFatJetOverlap, Jet_enu))

        Electron_enu = list(filter(FatJetConeIsolation, Electron_enu))
        Muon_enu = list(filter(FatJetConeIsolation, Muon_enu))

        Tau_enu = list(filter(FatJetTauOverlap, Tau_enu))
        boostedTau_enu = list(filter(FatJetTauOverlap, boostedTau_enu))

        Tau_enu = [x for x in Tau_enu if ElectronTauOverlap(x, Electron_enu)]
        boostedTau_enu = [
            x for x in boostedTau_enu if ElectronTauOverlap(
                x, Electron_enu)]

        Tau_enu = [x for x in Tau_enu if MuonTauOverlap(x, Muon_enu)]
        boostedTau_enu = [
            x for x in boostedTau_enu if MuonTauOverlap(
                x, Muon_enu)]

        list = {}
        list["bb"] = selfPairing(boostedTau_enu, "bb")
        list["tt"] = selfPairing(Tau_enu, "tt")
        list["be"] = crossPairing(boostedTau_enu, Electron_enu, "be")
        list["bm"] = crossPairing(boostedTau_enu, Muon_enu, "bm")
        list["te"] = crossPairing(Tau_enu, Electron_enu, "te")
        list["tm"] = crossPairing(Tau_enu, Muon_enu, "tm")

        Keymax = max(list, key=lambda x: list[x][0])

        if (list[Keymax][0] < 0):
            return False
        gFatJet_index = [FatJet_enu[0][0]]
        gTau_index = []
        gboostedTau_index = []
        gElectron_index = []
        gMuon_index = []
        gJet_index = []
        firstLepton = fastMTTlepton()
        secondLepton = fastMTTlepton()
        theMET = fastMTTmet(
            measuredX=event.METcorrected_pt * math.cos(event.METcorrected_phi),
            measuredY=event.METcorrected_pt * math.sin(event.METcorrected_phi),
            xx=event.MET_covXX,
            xy=event.MET_covXY,
            yy=event.MET_covYY)

        if Keymax == "bb":
            self.out.fillBranch("channel", 0)
            self.out.fillBranch("boost", 1)
            self.out.fillBranch("nallTaus", 2)
            gboostedTau_index = [list[Keymax][1], list[Keymax][2]]
            self.out.fillBranch("allTaus_decayMode",
                                [boostedTau[gboostedTau_index[0]].decayMode,
                                 boostedTau[gboostedTau_index[1]].decayMode])
            firstLepton = fastMTTlepton(pt=boostedTau[gboostedTau_index[0]].pt,
                                        eta=boostedTau[gboostedTau_index[0]].eta,
                                        phi=boostedTau[gboostedTau_index[0]].phi,
                                        m=boostedTau[gboostedTau_index[0]].mass,
                                        leptonType='Tau',
                                        tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
            self.pair1FV.SetPtEtaPhiM(boostedTau[gboostedTau_index[0]].pt,
                                      boostedTau[gboostedTau_index[0]].eta,
                                      boostedTau[gboostedTau_index[0]].phi,
                                      boostedTau[gboostedTau_index[0]].mass)
            secondLepton = fastMTTlepton(pt=boostedTau[gboostedTau_index[1]].pt,
                                         eta=boostedTau[gboostedTau_index[1]].eta,
                                         phi=boostedTau[gboostedTau_index[1]].phi,
                                         m=boostedTau[gboostedTau_index[1]].mass,
                                         leptonType='Tau',
                                         tauDecayMode=boostedTau[gboostedTau_index[1]].decayMode)
            self.pair2FV.SetPtEtaPhiM(boostedTau[gboostedTau_index[1]].pt,
                                      boostedTau[gboostedTau_index[1]].eta,
                                      boostedTau[gboostedTau_index[1]].phi,
                                      boostedTau[gboostedTau_index[1]].mass)

        elif Keymax == "tt":
            self.out.fillBranch("channel", 0)
            self.out.fillBranch("boost", 0)
            self.out.fillBranch("nallTaus", 2)
            gTau_index = [list[Keymax][1], list[Keymax][2]]
            self.out.fillBranch("allTaus_decayMode",
                                [Tau[gTau_index[0]].decayMode,
                                 Tau[gTau_index[0]].decayMode])
            firstLepton = fastMTTlepton(pt=Tau[gTau_index[0]].pt,
                                        eta=Tau[gTau_index[0]].eta,
                                        phi=Tau[gTau_index[0]].phi,
                                        m=Tau[gTau_index[0]].mass,
                                        leptonType='Tau',
                                        tauDecayMode=Tau[gTau_index[0]].decayMode)
            self.pair1FV.SetPtEtaPhiM(Tau[gTau_index[0]].pt,
                                      Tau[gTau_index[0]].eta,
                                      Tau[gTau_index[0]].phi,
                                      Tau[gTau_index[0]].mass)
            secondLepton = fastMTTlepton(pt=Tau[gTau_index[1]].pt,
                                         eta=Tau[gTau_index[1]].eta,
                                         phi=Tau[gTau_index[1]].phi,
                                         m=Tau[gTau_index[1]].mass,
                                         leptonType='Tau',
                                         tauDecayMode=Tau[gTau_index[1]].decayMode)
            self.pair2FV.SetPtEtaPhiM(Tau[gTau_index[1]].pt,
                                      Tau[gTau_index[1]].eta,
                                      Tau[gTau_index[1]].phi,
                                      Tau[gTau_index[1]].mass)

        elif Keymax == "be":
            self.out.fillBranch("channel", 1)
            self.out.fillBranch("boost", 1)
            self.out.fillBranch("nallTaus", 1)
            gboostedTau_index = [list[Keymax][1]]
            gElectron_index = [list[Keymax][2]]
            self.out.fillBranch("allTaus_decayMode",
                                [boostedTau[gboostedTau_index[0]].decayMode])
            firstLepton = fastMTTlepton(pt=boostedTau[gboostedTau_index[0]].pt,
                                        eta=boostedTau[gboostedTau_index[0]].eta,
                                        phi=boostedTau[gboostedTau_index[0]].phi,
                                        m=boostedTau[gboostedTau_index[0]].mass,
                                        leptonType='Tau',
                                        tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
            self.pair1FV.SetPtEtaPhiM(boostedTau[gboostedTau_index[0]].pt,
                                      boostedTau[gboostedTau_index[0]].eta,
                                      boostedTau[gboostedTau_index[0]].phi,
                                      boostedTau[gboostedTau_index[0]].mass)
            secondLepton = fastMTTlepton(pt=Electron[gElectron_index[0]].pt,
                                         eta=Electron[gElectron_index[0]].eta,
                                         phi=Electron[gElectron_index[0]].phi,
                                         m=0.51100e-3,
                                         leptonType='Electron',
                                         tauDecayMode=-1)
            self.pair2FV.SetPtEtaPhiM(Electron[gElectron_index[0]].pt,
                                      Electron[gElectron_index[0]].eta,
                                      Electron[gElectron_index[0]].phi,
                                      0.0)

        elif Keymax == "te":
            self.out.fillBranch("channel", 1)
            self.out.fillBranch("boost", 0)
            self.out.fillBranch("nallTaus", 1)
            gTau_index = [list[Keymax][1]]
            gElectron_index = [list[Keymax][2]]
            self.out.fillBranch("allTaus_decayMode", [
                                Tau[gTau_index[0]].decayMode])
            firstLepton = fastMTTlepton(pt=Tau[gTau_index[0]].pt,
                                        eta=Tau[gTau_index[0]].eta,
                                        phi=Tau[gTau_index[0]].phi,
                                        m=Tau[gTau_index[0]].mass,
                                        leptonType='Tau',
                                        tauDecayMode=Tau[gTau_index[0]].decayMode)
            self.pair1FV.SetPtEtaPhiM(Tau[gTau_index[0]].pt,
                                      Tau[gTau_index[0]].eta,
                                      Tau[gTau_index[0]].phi,
                                      Tau[gTau_index[0]].mass)
            secondLepton = fastMTTlepton(pt=Electron[gElectron_index[0]].pt,
                                         eta=Electron[gElectron_index[0]].eta,
                                         phi=Electron[gElectron_index[0]].phi,
                                         m=0.51100e-3,
                                         leptonType='Electron',
                                         tauDecayMode=-1)
            self.pair2FV.SetPtEtaPhiM(Electron[gElectron_index[0]].pt,
                                      Electron[gElectron_index[0]].eta,
                                      Electron[gElectron_index[0]].phi,
                                      0.0)

        elif Keymax == "bm":
            self.out.fillBranch("channel", 2)
            self.out.fillBranch("boost", 1)
            self.out.fillBranch("nallTaus", 1)
            gboostedTau_index = [list[Keymax][1]]
            gMuon_index = [list[Keymax][2]]
            self.out.fillBranch("allTaus_decayMode",
                                [boostedTau[gboostedTau_index[0]].decayMode])
            firstLepton = fastMTTlepton(pt=boostedTau[gboostedTau_index[0]].pt,
                                        eta=boostedTau[gboostedTau_index[0]].eta,
                                        phi=boostedTau[gboostedTau_index[0]].phi,
                                        m=boostedTau[gboostedTau_index[0]].mass,
                                        leptonType='Tau',
                                        tauDecayMode=boostedTau[gboostedTau_index[0]].decayMode)
            self.pair1FV.SetPtEtaPhiM(boostedTau[gboostedTau_index[0]].pt,
                                      boostedTau[gboostedTau_index[0]].eta,
                                      boostedTau[gboostedTau_index[0]].phi,
                                      boostedTau[gboostedTau_index[0]].mass)
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

        elif Keymax == "tm":
            self.out.fillBranch("channel", 2)
            self.out.fillBranch("boost", 0)
            self.out.fillBranch("nallTaus", 1)
            gTau_index = [list[Keymax][1]]
            gMuon_index = [list[Keymax][2]]
            self.out.fillBranch("allTaus_decayMode", [
                                Tau[gTau_index[0]].decayMode])
            firstLepton = fastMTTlepton(pt=Tau[gTau_index[0]].pt,
                                        eta=Tau[gTau_index[0]].eta,
                                        phi=Tau[gTau_index[0]].phi,
                                        m=Tau[gTau_index[0]].mass,
                                        leptonType='Tau',
                                        tauDecayMode=Tau[gTau_index[0]].decayMode)
            self.pair1FV.SetPtEtaPhiM(Tau[gTau_index[0]].pt,
                                      Tau[gTau_index[0]].eta,
                                      Tau[gTau_index[0]].phi,
                                      Tau[gTau_index[0]].mass)
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
                Muon)]
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
        theFastMTTtool = fastMTTtool()
        theFastMTTtool.setFirstLepton(firstLepton)
        theFastMTTtool.setSecondLepton(secondLepton)
        theFastMTTtool.setTheMET(theMET)
        higgsFV_list = theFastMTTtool.getFastMTTfourvector()
        self.higgsTTFV.SetPtEtaPhiM(
            higgsFV_list[0],
            higgsFV_list[1],
            higgsFV_list[2],
            higgsFV_list[3])
        self.RadionFV = self.higgsTTFV + self.higgsBBFV

        self.higgsTTvisFV = self.pair1FV + self.pair2FV
        self.RadionvisFV = self.higgsTTvisFV + self.higgsBBFV

        self.out.fillBranch("HTT_m", self.higgsTTFV.M())
        self.out.fillBranch("HTT_eta", self.higgsTTFV.Eta())
        self.out.fillBranch("HTT_phi", self.higgsTTFV.Phi())
        self.out.fillBranch("HTT_pt", self.higgsTTFV.Pt())
        self.out.fillBranch("HTTvis_m", self.higgsTTvisFV.M())
        self.out.fillBranch("HTTvis_eta", self.higgsTTvisFV.Eta())
        self.out.fillBranch("HTTvis_phi", self.higgsTTvisFV.Phi())
        self.out.fillBranch("HTTvis_pt", self.higgsTTvisFV.Pt())
        self.out.fillBranch("HTTvis_deltaR", self.pair1FV.DeltaR(self.pair2FV))

        self.out.fillBranch("X_m", self.RadionFV.M())
        self.out.fillBranch("X_eta", self.RadionFV.Eta())
        self.out.fillBranch("X_phi", self.RadionFV.Phi())
        self.out.fillBranch("X_pt", self.RadionFV.Pt())
        self.out.fillBranch("Xvis_m", self.RadionvisFV.M())
        self.out.fillBranch("Xvis_eta", self.RadionvisFV.Eta())
        self.out.fillBranch("Xvis_phi", self.RadionvisFV.Phi())
        self.out.fillBranch("Xvis_pt", self.RadionvisFV.Pt())

        self.out.fillBranch("ngood_boostedTaus", len(gboostedTau_index))
        self.out.fillBranch("ngood_Taus", len(gTau_index))
        self.out.fillBranch("ngood_Electrons", len(gElectron_index))
        self.out.fillBranch("ngood_Muons", len(gMuon_index))
        self.out.fillBranch("ngood_FatJets", len(gFatJet_index))
        self.out.fillBranch("ngood_Jets", len(gJet_index))
        self.out.fillBranch("ngood_LooseJets", len(gJet_Looseindex))
        self.out.fillBranch("ngood_MediumJets", len(gJet_Mediumindex))
        self.out.fillBranch("ngood_TightJets", len(gJet_Tightindex))

        self.out.fillBranch("index_gboostedTaus", gboostedTau_index)
        self.out.fillBranch("index_gTaus", gTau_index)
        self.out.fillBranch("index_gElectrons", gElectron_index)
        self.out.fillBranch("index_gMuons", gMuon_index)
        self.out.fillBranch("index_gFatJets", gFatJet_index)
        self.out.fillBranch("index_gJets", gJet_index)
        self.out.fillBranch("index_gLooseJets", gJet_Looseindex)
        self.out.fillBranch("index_gMediumJets", gJet_Mediumindex)
        self.out.fillBranch("index_gTightJets", gJet_Tightindex)

        if ((Keymax == "bb") or (Keymax == "tt")):
            self.out.fillBranch(
                "allTaus_pt", [
                    self.pair1FV.Pt(), self.pair2FV.Pt()])
            self.out.fillBranch(
                "allTaus_eta", [
                    self.pair1FV.Eta(), self.pair2FV.Eta()])
            self.out.fillBranch(
                "allTaus_phi", [
                    self.pair1FV.Phi(), self.pair2FV.Phi()])
            self.out.fillBranch(
                "allTaus_mass", [
                    self.pair1FV.M(), self.pair2FV.M()])
        else:
            self.out.fillBranch("allTaus_pt", [self.pair1FV.Pt()])
            self.out.fillBranch("allTaus_eta", [self.pair1FV.Eta()])
            self.out.fillBranch("allTaus_phi", [self.pair1FV.Phi()])
            self.out.fillBranch("allTaus_mass", [self.pair1FV.M()])
        return True


def call_postpoc(files):
    nameStrip = files.strip()
    filename = (nameStrip.split('/')[-1]).split('.')[-2]
    print(filename)

    if (search("Run", filename)):
        print(("This is a ", args.year, " Data file = ", filename))

        def mainModule(): return cutsAndcategories(
            args.year, True, args.eleID, args.muID, args.tauID, args.btauID)
    else:
        print(("This is a ", args.year, " MC file = ", filename))

        def mainModule(): return cutsAndcategories(
            args.year, False, args.eleID, args.muID, args.tauID, args.btauID)

    p = PostProcessor(
        args.tempStore,
        [files],
        cut=None,
        branchsel=None,
        modules=[
            mainModule()],
        postfix=post,
        noOut=False,
        outputbranchsel=None)
    p.run()
    print("###############MOVING THE OUTPUT FILE BACK TO HDFS#######################")
    os.system(
        "hadoop fs -moveFromLocal -f " +
        args.tempStore +
        "/" +
        filename +
        ".root" +
        " " +
        outputDir +
        "/.")  # This currently doesnot work due to username differences - it takes parida by default


if __name__ == "__main__":
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
        '--postfix',
        help="string at the end of output file names",
        default="")
    parser.add_argument(
        '--year',
        '-y',
        help='specify the run - to make sure right triggers are used',
        choices=[
            '2016',
            '2016APV',
            '2017',
            '2018'])
    parser.add_argument(
        '--eleID',
        '-eId',
        help="Specify the ID for the electrons to be selected 2-Loose, 3-Medium, 4-Tight",
        choices=[
            '2',
            '3',
            '4'])
    parser.add_argument(
        '--muID',
        '-mId',
        help="Specify the ID for the muons to be selected 2-Loose, 3-Medium, 4-Tight (to keep it consistent with the electrons)",
        choices=[
            '2',
            '3',
            '4'])
    parser.add_argument(
        '--tauID',
        '-tId',
        help="DeepTau 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose, 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight",
        choices=[
            '1',
            '2',
            '4',
            '8',
            '16',
            '32',
            '64',
            '128'])
    parser.add_argument(
        '--btauID',
        '-btId',
        help="DeepBoostedTau Raw score",
        required=True)
    parser.add_argument(
        '--tempStore',
        '-t',
        help='Temporary staging area for files before moving out to hdfs',
        required=True)

    args = parser.parse_args()

    # making a list of input MC/DATA files
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
            pool = np.Pool(int(args.ncores))
            print(("list", argList))
            res = pool.map(call_postpoc, argList)
        except Exception as error:
            print("MultiProc error - needs debugging")
            print(("An exception occurred:", type(error).__name__))
