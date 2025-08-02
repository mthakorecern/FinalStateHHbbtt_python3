import argparse
import multiprocessing as np
from re import search
import numpy as np
import glob
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
# Import the JEC and JER NanoAODTools correction Modules

print(('gRoot ProcessLine XYMETCorrection = ', ROOT.gROOT.ProcessLine(
    ".L /afs/hep.wisc.edu/home/parida/public/HHbbtt_Analysis_Scripts/CMSSW_10_6_27/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/XYMETCorrection/ULMETXY_Correction.h")))


class addJecsMET(Module):
    def __init__(self, year, isData):
        self.year = year
        self.isMC = not isData
        self.isData = isData
        self.jesSys = "Total"

        if self.jesSys == "All":
            self.jesUnc = [
                "",
                "AbsoluteMPFBias",
                "AbsoluteScale",
                "AbsoluteStat",
                "FlavorQCD",
                "Fragmentation",
                "PileUpDataMC",
                "PileUpPtBB",
                "PileUpPtEC1",
                "PileUpPtEC2",
                "PileUpPtHF",
                "PileUpPtRef",
                "RelativeFSR",
                "RelativeJEREC1",
                "RelativeJEREC2",
                "RelativeJERHF",
                "RelativePtBB",
                "RelativePtEC1",
                "RelativePtEC2",
                "RelativePtHF",
                "RelativeBal",
                "RelativeSample",
                "RelativeStatEC",
                "RelativeStatFSR",
                "RelativeStatHF",
                "SinglePionECAL",
                "SinglePionHCAL",
                "TimePtEta"]
        elif self.jesSys == "Total":
            self.jesUnc = [""]

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("METcorrected_pt", "F")
        self.out.branch("METcorrected_phi", "F")
        if self.isMC:
            for sys in self.jesUnc:
                self.out.branch("METcorrected_ptScale" + sys + "Up", "F")
                self.out.branch("METcorrected_ptScale" + sys + "Down", "F")
                self.out.branch("METcorrected_phiScale" + sys + "Up", "F")
                self.out.branch("METcorrected_phiScale" + sys + "Down", "F")

            self.out.branch("METcorrected_ptResUp", "F")
            self.out.branch("METcorrected_ptResDown", "F")
            self.out.branch("METcorrected_phiResUp", "F")
            self.out.branch("METcorrected_phiResDown", "F")
            self.out.branch("METcorrected_ptUnclustUp", "F")
            self.out.branch("METcorrected_ptUnclustDown", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        FatJet = Collection(event, 'FatJet', 'nFatJet')
        goodFatJet = [x for x in FatJet if (
            (x.pt > 170) and (abs(x.eta) < 2.5) and (x.jetId > 1))]
        if (len(goodFatJet) == 0):
            return False
        jesBranches = {}

        def getJESMETpt(sys):
            if sys == "Up":
                return event.MET_T1Smear_pt_jesTotalUp
            elif sys == "Down":
                return event.MET_T1Smear_pt_jesTotalDown

        def getJESMETphi(sys):
            if sys == "Up":
                return event.MET_T1Smear_phi_jesTotalUp
            elif sys == "Down":
                return event.MET_T1Smear_phi_jesTotalDown

        # METXYCorr_Met_MetPhi(double uncormet, double uncormet_phi, int runnb, TString year, bool isMC, int npv, bool isUL =false,bool ispuppi=false)
        # if(isMC && year == "2016" && !isUL) runera = y2016MC;
        # else if(isMC && year == "2017" && !isUL) {runera = y2017MC; usemetv2 =true;}
        # else if(isMC && year == "2018" && !isUL) runera = y2018MC;
        # else if(isMC && year == "2016APV" && isUL) runera = yUL2016MCAPV;
        # else if(isMC && year == "2016nonAPV" && isUL) runera = yUL2016MCnonAPV;
        # else if(isMC && year == "2017" && isUL) runera = yUL2017MC;
        # else if(isMC && year == "2018" && isUL) runera = yUL2018MC;

        if self.isData:
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(
                event.MET_T1_pt,
                event.MET_T1_phi,
                event.run,
                self.year,
                self.isMC,
                event.PV_npvs,
                True,
                False)
        else:
            METcorrected_pt_phi = ROOT.METXYCorr_Met_MetPhi(
                event.MET_T1Smear_pt,
                event.MET_T1Smear_phi,
                event.run,
                self.year,
                self.isMC,
                event.PV_npvs,
                True,
                False)
        METcorrected_pt = METcorrected_pt_phi[0]
        METcorrected_phi = METcorrected_pt_phi[1]

        if self.isMC:
            for sys in self.jesUnc:
                jesBranches["METcorrected_pt_phiScale" +
                            sys +
                            "Up"] = ROOT.METXYCorr_Met_MetPhi(getJESMETpt(sys +
                                                                          "Up"), getJESMETphi(sys +
                                                                                              "Up"), event.run, self.year, self.isMC, event.PV_npvs, True, False)
                jesBranches["METcorrected_pt_phiScale" +
                            sys +
                            "Down"] = ROOT.METXYCorr_Met_MetPhi(getJESMETpt(sys +
                                                                            "Down"), getJESMETphi(sys +
                                                                                                  "Down"), event.run, self.year, self.isMC, event.PV_npvs, True, False)
            METcorrected_pt_phiResUp = ROOT.METXYCorr_Met_MetPhi(
                event.MET_T1Smear_pt_jerUp,
                event.MET_T1Smear_phi_jerUp,
                event.run,
                self.year,
                self.isMC,
                event.PV_npvs,
                True,
                False)
            METcorrected_pt_phiResDown = ROOT.METXYCorr_Met_MetPhi(
                event.MET_T1Smear_pt_jerDown,
                event.MET_T1Smear_phi_jerDown,
                event.run,
                self.year,
                self.isMC,
                event.PV_npvs,
                True,
                False)

            for sys in self.jesUnc:
                jesBranches["METcorrected_ptScale" +
                            sys +
                            "Up"] = jesBranches["METcorrected_pt_phiScale" +
                                                sys +
                                                "Up"][0]
                jesBranches["METcorrected_phiScale" +
                            sys +
                            "Up"] = jesBranches["METcorrected_pt_phiScale" +
                                                sys +
                                                "Up"][1]
                jesBranches["METcorrected_ptScale" +
                            sys +
                            "Down"] = jesBranches["METcorrected_pt_phiScale" +
                                                  sys +
                                                  "Down"][0]
                jesBranches["METcorrected_phiScale" +
                            sys +
                            "Down"] = jesBranches["METcorrected_pt_phiScale" +
                                                  sys +
                                                  "Down"][1]

            METcorrected_ptResUp = METcorrected_pt_phiResUp[0]
            METcorrected_phiResUp = METcorrected_pt_phiResUp[1]
            METcorrected_ptResDown = METcorrected_pt_phiResDown[0]
            METcorrected_phiResDown = METcorrected_pt_phiResDown[1]

            METcorrected_ptX = METcorrected_pt * math.cos(METcorrected_phi)
            METcorrected_ptY = METcorrected_pt * math.sin(METcorrected_phi)

            METcorrected_ptUnclustXUp = METcorrected_ptX + event.MET_MetUnclustEnUpDeltaX
            METcorrected_ptUnclustYUp = METcorrected_ptY + event.MET_MetUnclustEnUpDeltaY
            METcorrected_ptUnclustUp = math.sqrt(
                pow(METcorrected_ptUnclustXUp, 2) + pow(METcorrected_ptUnclustYUp, 2))

            METcorrected_ptUnclustXDown = METcorrected_ptX - event.MET_MetUnclustEnUpDeltaX
            METcorrected_ptUnclustYDown = METcorrected_ptY - event.MET_MetUnclustEnUpDeltaY
            METcorrected_ptUnclustDown = math.sqrt(
                pow(METcorrected_ptUnclustXDown, 2) + pow(METcorrected_ptUnclustYDown, 2))

        self.out.fillBranch("METcorrected_pt", METcorrected_pt)
        self.out.fillBranch("METcorrected_phi", METcorrected_phi)
        if self.isMC:
            for sys in self.jesUnc:
                self.out.fillBranch("METcorrected_ptScale" + sys + "Up",
                                    jesBranches["METcorrected_ptScale" + sys + "Up"])
                self.out.fillBranch("METcorrected_ptScale" +
                                    sys +
                                    "Down", jesBranches["METcorrected_ptScale" +
                                                        sys +
                                                        "Down"])
                self.out.fillBranch("METcorrected_phiScale" + sys + "Up",
                                    jesBranches["METcorrected_phiScale" + sys + "Up"])
                self.out.fillBranch("METcorrected_phiScale" +
                                    sys +
                                    "Down", jesBranches["METcorrected_phiScale" +
                                                        sys +
                                                        "Down"])
            self.out.fillBranch("METcorrected_ptResUp", METcorrected_ptResUp)
            self.out.fillBranch(
                "METcorrected_ptResDown",
                METcorrected_ptResDown)
            self.out.fillBranch("METcorrected_phiResUp", METcorrected_phiResUp)
            self.out.fillBranch(
                "METcorrected_phiResDown",
                METcorrected_phiResDown)
            # Systematics - Unclustered MET
            self.out.fillBranch(
                "METcorrected_ptUnclustUp",
                METcorrected_ptUnclustUp)
            self.out.fillBranch(
                "METcorrected_ptUnclustDown",
                METcorrected_ptUnclustDown)
        return True


def call_postpoc(files):
    nameStrip = files.strip()
    filename = (nameStrip.split('/')[-1]).split('.')[-2]
    print(filename)
    # Create the JME modules 2016 Data and MC
    AK4jetmetCorrectorUL2016MC = createJMECorrector(
        isMC=True, dataYear='UL2016', jesUncert="Total", jetType="AK4PFchs")
    AK4jetmetCorrectorUL2016RunF = createJMECorrector(
        isMC=False,
        dataYear='UL2016',
        runPeriod="F",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2016RunG = createJMECorrector(
        isMC=False,
        dataYear='UL2016',
        runPeriod="G",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2016RunH = createJMECorrector(
        isMC=False,
        dataYear='UL2016',
        runPeriod="H",
        jesUncert="Total",
        jetType="AK4PFchs")

    AK8jetmetCorrectorUL2016MC = createJMECorrector(
        isMC=True, dataYear='UL2016', jesUncert="Total", jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2016RunF = createJMECorrector(
        isMC=False,
        dataYear='UL2016',
        runPeriod="F",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2016RunG = createJMECorrector(
        isMC=False,
        dataYear='UL2016',
        runPeriod="G",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2016RunH = createJMECorrector(
        isMC=False,
        dataYear='UL2016',
        runPeriod="H",
        jesUncert="Total",
        jetType="AK8PFPuppi")

    # Create the JME modules 2016APV Data and MC
    AK4jetmetCorrectorUL2016preVFPMC = createJMECorrector(
        isMC=True, dataYear='UL2016_preVFP', jesUncert="Total", jetType="AK4PFchs")
    AK4jetmetCorrectorUL2016preVFPRunB = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="B",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2016preVFPRunC = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="C",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2016preVFPRunD = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="D",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2016preVFPRunE = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="E",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2016preVFPRunF = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="F",
        jesUncert="Total",
        jetType="AK4PFchs")

    AK8jetmetCorrectorUL2016preVFPMC = createJMECorrector(
        isMC=True, dataYear='UL2016_preVFP', jesUncert="Total", jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2016preVFPRunB = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="B",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2016preVFPRunC = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="C",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2016preVFPRunD = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="D",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2016preVFPRunE = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="E",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2016preVFPRunF = createJMECorrector(
        isMC=False,
        dataYear='UL2016_preVFP',
        runPeriod="F",
        jesUncert="Total",
        jetType="AK8PFPuppi")

    # Create the JME modules 2017 Data and MC
    AK4jetmetCorrectorUL2017MC = createJMECorrector(
        isMC=True, dataYear='UL2017', jesUncert="Total", jetType="AK4PFchs")
    AK4jetmetCorrectorUL2017RunB = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="B",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2017RunC = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="C",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2017RunD = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="D",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2017RunE = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="E",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2017RunF = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="F",
        jesUncert="Total",
        jetType="AK4PFchs")

    AK8jetmetCorrectorUL2017MC = createJMECorrector(
        isMC=True, dataYear='UL2017', jesUncert="Total", jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2017RunB = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="B",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2017RunC = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="C",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2017RunD = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="D",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2017RunE = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="E",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2017RunF = createJMECorrector(
        isMC=False,
        dataYear='UL2017',
        runPeriod="F",
        jesUncert="Total",
        jetType="AK8PFPuppi")

    # Create the JME modules 2018 Data and MC
    AK4jetmetCorrectorUL2018MC = createJMECorrector(
        isMC=True, dataYear='UL2018', jesUncert="Total", jetType="AK4PFchs")
    AK4jetmetCorrectorUL2018RunA = createJMECorrector(
        isMC=False,
        dataYear='UL2018',
        runPeriod="A",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2018RunB = createJMECorrector(
        isMC=False,
        dataYear='UL2018',
        runPeriod="B",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2018RunC = createJMECorrector(
        isMC=False,
        dataYear='UL2018',
        runPeriod="C",
        jesUncert="Total",
        jetType="AK4PFchs")
    AK4jetmetCorrectorUL2018RunD = createJMECorrector(
        isMC=False,
        dataYear='UL2018',
        runPeriod="D",
        jesUncert="Total",
        jetType="AK4PFchs")

    AK8jetmetCorrectorUL2018MC = createJMECorrector(
        isMC=True, dataYear='UL2018', jesUncert="Total", jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2018RunA = createJMECorrector(
        isMC=False,
        dataYear='UL2018',
        runPeriod="A",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2018RunB = createJMECorrector(
        isMC=False,
        dataYear='UL2018',
        runPeriod="B",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2018RunC = createJMECorrector(
        isMC=False,
        dataYear='UL2018',
        runPeriod="C",
        jesUncert="Total",
        jetType="AK8PFPuppi")
    AK8jetmetCorrectorUL2018RunD = createJMECorrector(
        isMC=False,
        dataYear='UL2018',
        runPeriod="D",
        jesUncert="Total",
        jetType="AK8PFPuppi")

    if (args.year == "2016nonAPV"):
        if (search("Run", filename)):
            def addJecsAndMET(): return addJecsMET(args.year, isData=True)
            if (search("SingleMu_Run2016F", filename)):
                print("This is a 2016 SingleMu_Run2016F Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2016RunF(),
                        AK8jetmetCorrectorUL2016RunF(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2016_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2016G", filename)):
                print("This is a 2016 SingleMu_Run2016G Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2016RunG(),
                        AK8jetmetCorrectorUL2016RunG(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2016_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2016H", filename)):
                print("This is a 2016 SingleMu_Run2016H Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2016RunH(),
                        AK8jetmetCorrectorUL2016RunH(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2016_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            else:
                print("Incorrect data file name for 2016")

        else:
            print("This is a 2016 MC file")
            def addJecsAndMET(): return addJecsMET(args.year, isData=False)
            p = PostProcessor(
                args.tempStore,
                [files],
                cut=preselection,
                branchsel="Keep_Drop_txt/MC_keep_and_drop_SFs.txt",
                modules=[
                    AK4jetmetCorrectorUL2016MC(),
                    AK8jetmetCorrectorUL2016MC(),
                    addJecsAndMET()],
                postfix=post,
                noOut=False,
                outputbranchsel="Keep_Drop_txt/MC_keep_and_drop_SFs.txt")

    if (args.year == "2016APV"):
        if (search("Run", filename)):
            def addJecsAndMET(): return addJecsMET(args.year, isData=True)
            if (search("SingleMu_Run2016B-ver1_HIPM", filename)):
                print("This is a 2016APV SingleMu_Run2016B-ver1_HIPM Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2016preVFPRunB(),
                        AK8jetmetCorrectorUL2016preVFPRunB(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2016_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2016B-ver2_HIPM", filename)):
                print("This is a 2016APV SingleMu_Run2016B-ver2_HIPM Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2016preVFPRunB(),
                        AK8jetmetCorrectorUL2016preVFPRunB(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2016_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2016C-HIPM", filename)):
                print("This is a 2016APV SingleMu_Run2016C-HIPM Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2016preVFPRunC(),
                        AK8jetmetCorrectorUL2016preVFPRunC(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2016_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2016D-HIPM", filename)):
                print("This is a 2016APV SingleMu_Run2016D-HIPM Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2016preVFPRunD(),
                        AK8jetmetCorrectorUL2016preVFPRunD(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2016_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2016E-HIPM", filename)):
                print("This is a 2016APV SingleMu_Run2016E-HIPM Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2016preVFPRunE(),
                        AK8jetmetCorrectorUL2016preVFPRunE(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2016_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2016F-HIPM", filename)):
                print("This is a 2016APV SingleMu_Run2016F-HIPM Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2016preVFPRunF(),
                        AK8jetmetCorrectorUL2016preVFPRunF(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2016_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")

        else:
            print("This is a 2016APV MC file")
            def addJecsAndMET(): return addJecsMET(args.year, isData=False)
            p = PostProcessor(
                args.tempStore,
                [files],
                cut=preselection,
                branchsel="Keep_Drop_txt/MC_keep_and_drop_SFs.txt",
                modules=[
                    AK4jetmetCorrectorUL2016preVFPMC(),
                    AK8jetmetCorrectorUL2016preVFPMC(),
                    addJecsAndMET()],
                postfix=post,
                noOut=False,
                outputbranchsel="Keep_Drop_txt/MC_keep_and_drop_SFs.txt")

    if (args.year == "2017"):
        if (search("Run", filename)):
            def addJecsAndMET(): return addJecsMET(args.year, isData=True)
            if (search("SingleMu_Run2017B", filename)):
                print("This is a 2017 SingleMu_Run2017B Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2017RunB(),
                        AK8jetmetCorrectorUL2017RunB(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2017_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2017C", filename)):
                print("This is a 2017 SingleMu_Run2017C Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2017RunC(),
                        AK8jetmetCorrectorUL2017RunC(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2017_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2017D", filename)):
                print("This is a 2017 SingleMu_Run2017D Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2017RunD(),
                        AK8jetmetCorrectorUL2017RunD(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2017_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2017E", filename)):
                print("This is a 2017 SingleMu_Run2017E Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2017RunE(),
                        AK8jetmetCorrectorUL2017RunE(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2017_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2017F", filename)):
                print("This is a 2017 SingleMu_Run2017F Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2017RunF(),
                        AK8jetmetCorrectorUL2017RunF(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2017_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")

        else:
            print("This is a 2017 MC file")
            def addJecsAndMET(): return addJecsMET(args.year, isData=False)
            p = PostProcessor(
                args.tempStore,
                [files],
                cut=preselection,
                branchsel="Keep_Drop_txt/MC_keep_and_drop_SFs.txt",
                modules=[
                    AK4jetmetCorrectorUL2017MC(),
                    AK8jetmetCorrectorUL2017MC(),
                    addJecsAndMET()],
                postfix=post,
                noOut=False,
                outputbranchsel="Keep_Drop_txt/MC_keep_and_drop_SFs.txt")

    if (args.year == "2018"):
        if (search("Run", filename)):
            def addJecsAndMET(): return addJecsMET(args.year, isData=True)
            if (search("SingleMu_Run2018A", filename)):
                print("This is a 2018 SingleMu_Run2018A Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2018RunA(),
                        AK8jetmetCorrectorUL2018RunA(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2018_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2018B", filename)):
                print("This is a 2018 SingleMu_Run2018B Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2018RunB(),
                        AK8jetmetCorrectorUL2018RunB(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2018_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2018C", filename)):
                print("This is a 2018 SingleMu_Run2018C Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2018RunC(),
                        AK8jetmetCorrectorUL2018RunC(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2018_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
            elif (search("SingleMu_Run2018D", filename)):
                print("This is a 2018 SingleMu_Run2018D Data file")
                p = PostProcessor(
                    args.tempStore,
                    [files],
                    cut=preselection,
                    branchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt",
                    modules=[
                        AK4jetmetCorrectorUL2018RunD(),
                        AK8jetmetCorrectorUL2018RunD(),
                        addJecsAndMET()],
                    postfix=post,
                    jsonInput="GOLDEN_JSONS/2018_Leagcy.json",
                    noOut=False,
                    outputbranchsel="Keep_Drop_txt/Data_keep_and_drop_SFs.txt")
        else:
            print("This is a 2018 MC file")
            def addJecsAndMET(): return addJecsMET(args.year, isData=False)
            p = PostProcessor(
                args.tempStore,
                [files],
                cut=preselection,
                branchsel="Keep_Drop_txt/MC_keep_and_drop_SFs.txt",
                modules=[
                    AK4jetmetCorrectorUL2018MC(),
                    AK8jetmetCorrectorUL2018MC(),
                    addJecsAndMET()],
                postfix=post,
                noOut=False,
                outputbranchsel="Keep_Drop_txt/MC_keep_and_drop_SFs.txt")

    p.run()
    print("###############MOVING THE OUTPUT FILE BACK TO HDFS#######################")
    # os.system("hadoop fs -moveFromLocal -f
    # "+args.tempStore+"/"+filename+".root"+" "+outputDir+"/.") #This
    # currently doesnot work due to username differences - it takes parida by
    # default
    os.system(
        "mv " +
        args.tempStore +
        "/" +
        filename +
        ".root" +
        " " +
        outputDir +
        "/.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Script to Handle root file preparation to split into channels. Input should be a singular files for each dataset or data already with some basic selections applied')
    # parser.add_argument('--Channel',help="enter either tt or et or mt. For boostedTau test enter test",required=True,choices=['tt', 'et', 'mt'])
    parser.add_argument(
        '--inputLocation',
        '-i',
        help="enter the path to the location of input file set",
        default="")
    parser.add_argument(
        '--outputLocation',
        '-o',
        help="enter the path where yu want the output files to be stored",
        default="")
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
            '2016nonAPV',
            '2016APV',
            '2017',
            '2018'])
    parser.add_argument(
        '--tempStore',
        '-t',
        help='Temporary staging area for files before moving out to hdfs',
        required=True)

    args = parser.parse_args()

    met_selection = ["MET_pt>150"]

    # Verified with UL Twiki:
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2016_data_and_MC_UL
    eventSelection2016 = ["(PV_ndof > 4)",
                          "(abs(PV_z) < 24)",
                          "(sqrt(PV_x*PV_x+PV_y*PV_y) < 2)",
                          "Flag_goodVertices",
                          "Flag_globalSuperTightHalo2016Filter",
                          "Flag_HBHENoiseFilter",
                          "Flag_HBHENoiseIsoFilter",
                          "Flag_EcalDeadCellTriggerPrimitiveFilter",
                          "Flag_BadPFMuonFilter",
                          "Flag_BadPFMuonDzFilter",
                          "Flag_eeBadScFilter",
                          "Flag_hfNoisyHitsFilter"]

    # Verified with UL Twiki:
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2018_2017_data_and_MC_UL
    eventSelection2017_18 = ["(PV_ndof > 4)",
                             "(abs(PV_z) < 24)",
                             "(sqrt(PV_x*PV_x+PV_y*PV_y) < 2)",
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

    trigger_2016APV = ["HLT_IsoMu24", "HLT_IsoTkMu24"]

    trigger_2016 = ["HLT_IsoMu24", "HLT_IsoTkMu24"]

    # trigger_2017_18=["HLT_PFMETNoMu120_PFMHTNoMu120_IDTight","HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight","HLT_PFMET120_PFMHT120_IDTight"]
    trigger_2017_18 = ["HLT_IsoMu24"]

    if (args.year == "2016nonAPV"):
        preselection = "(" + "&&".join(met_selection) + ")" + "&&" + "(" + "&&".join(
            eventSelection2016) + ")" + "&&" + "(" + "||".join(trigger_2016) + ")"
    elif (args.year == "2016APV"):
        preselection = "(" + "&&".join(met_selection) + ")" + "&&" + "(" + "&&".join(
            eventSelection2016) + ")" + "&&" + "(" + "||".join(trigger_2016APV) + ")"
    elif (args.year == "2017"):
        preselection = "(" + "&&".join(met_selection) + ")" + "&&" + "(" + "&&".join(
            eventSelection2017_18) + ")" + "&&" + "(" + "||".join(trigger_2017_18) + ")"
    elif (args.year == "2018"):
        preselection = "(" + "&&".join(met_selection) + ")" + "&&" + "(" + "&&".join(
            eventSelection2017_18) + ")" + "&&" + "(" + "||".join(trigger_2017_18) + ")"

    # making a list of input MC files
    fnames_mc = glob.glob(args.inputLocation + "/*.root")
    # making a list of input MC files
    fnames_data = glob.glob(args.inputLocation + "/SingleMu/*.root")
    fnames = fnames_mc + fnames_data
    # fnames = fnames_data

    outputDir = args.outputLocation

    # printing the preselection conditions
    print(("############################__YEAR__", args.year,
          "____##############################################"))
    print(("PRESELECTIONS = ", preselection))
    print("###################################################################################################")

    post = args.postfix
    argList = list()
    for file in fnames:
        argList.append(file)

    if int(args.ncores) == 1:
        for arr in argList:
            print("Using a single thread ")
            call_postpoc(arr)

    else:
        # try:
        print("Using Multithreading")
        pool = np.Pool(int(args.ncores))
        print(("list", argList))
        res = pool.map(call_postpoc, argList)
        # except Exception as error:
        # print ("MultiProc error - needs debugging")
        # print("An exception occurred:", type(error).__name__)
