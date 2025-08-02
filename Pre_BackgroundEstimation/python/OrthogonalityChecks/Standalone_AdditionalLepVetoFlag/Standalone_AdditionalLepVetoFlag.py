from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import argparse
import multiprocessing as mp
from re import search
import numpy as np
import glob
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True


class Standalone_AdditionalLepVetoFlag(Module):
    def __init__(self, year=2016, isData=False):
        print("ACTIAVTE: AddAK8MatchingToTaus")
        self.year = year
        self.isData = isData
        self.isMC = not self.isData
        self.higgsBBFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)
        self.lepFV = ROOT.TLorentzVector(0.0, 0.0, 0.0, 0.0)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # Define output branches to check if preselection is working
        self.out = wrappedOutputTree
        # Flag for addtional lepton veto
        self.out.branch("addlepton_vetoflag_all", "I")
        self.out.branch("addlepton_vetoflag_semi", "I")
        self.out.branch("nvetoElectron", "I")
        self.out.branch("index_gVetoElectrons", "I", lenVar="nvetoElectron")
        self.out.branch("nvetoMuon", "I")
        self.out.branch("index_gVetoMuons", "I", lenVar="nvetoMuon")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        # self.eventCount+=1
        FatJet = Collection(event, 'FatJet', 'nFatJet')
        self.higgsBBFV.SetPtEtaPhiM(FatJet[event.index_gFatJets[0]].pt_nom,
                                    FatJet[event.index_gFatJets[0]].eta,
                                    FatJet[event.index_gFatJets[0]].phi,
                                    FatJet[event.index_gFatJets[0]].mass_nom)

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

        def distinctFromExistinglepton(leptonObject_enu):
            # print ("See if it picking up the index = ",leptonObject_enu[0])
            if event.channel == 1:
                print("Picking up the same electron as Tau_Electron!!!!!!!!!!")
                if (event.index_gElectrons[0] == leptonObject_enu[0]):
                    return False
            elif event.channel == 2:
                print("Picking up the same muon as Muon!!!!!!!!!!")
                if (event.index_gMuons[0] == leptonObject_enu[0]):
                    return False
            return True

        # Create the collections
        Electron = Collection(event, 'Electron', 'nElectron')
        Muon = Collection(event, 'Muon', 'nMuon')

        #######################################################################
        # Code for addtional lepton veto
        Electron_addlep_enu = [x for x in enumerate(Electron) if (x[1].pt > 10) and (
            (abs(x[1].eta) <= 1.44) or (abs(x[1].eta) >= 1.57)) and (x[1].cutBased >= 2)]
        Electron_addlep_enu = list(
            filter(
                ElectronIsolationCut_addlepveto,
                Electron_addlep_enu))
        Electron_addlep_enu = list(
            filter(
                FatJetConeIsolation,
                Electron_addlep_enu))
        Electron_addlep_enu = list(
            filter(
                distinctFromExistinglepton,
                Electron_addlep_enu))

        Muon_addlep_enu = [x for x in enumerate(Muon) if x[1].pt > 10 and (
            abs(x[1].eta) < 2.5) and x[1].looseId and (x[1].pfRelIso04_all < 0.25)]
        Muon_addlep_enu = list(filter(FatJetConeIsolation, Muon_addlep_enu))
        Muon_addlep_enu = list(
            filter(
                distinctFromExistinglepton,
                Muon_addlep_enu))

        #######################################################################

        if ((event.channel == 0)):
            if ((len(Electron_addlep_enu) != 0) or (len(Muon_addlep_enu) != 0)):
                self.out.fillBranch("addlepton_vetoflag_all", 1)
                self.out.fillBranch("addlepton_vetoflag_semi", 0)
            elif ((len(Electron_addlep_enu) == 0) and (len(Muon_addlep_enu) == 0)):
                self.out.fillBranch("addlepton_vetoflag_all", 0)
                self.out.fillBranch("addlepton_vetoflag_semi", 0)
        else:
            if ((len(Electron_addlep_enu) != 0) or (len(Muon_addlep_enu) != 0)):
                self.out.fillBranch("addlepton_vetoflag_all", 1)
                self.out.fillBranch("addlepton_vetoflag_semi", 1)
            elif ((len(Electron_addlep_enu) == 0) and (len(Muon_addlep_enu) == 0)):
                self.out.fillBranch("addlepton_vetoflag_all", 0)
                self.out.fillBranch("addlepton_vetoflag_semi", 0)

        self.out.fillBranch("nvetoElectron", len(Electron_addlep_enu))
        self.out.fillBranch(
            "index_gVetoElectrons", [
                x[0] for x in Electron_addlep_enu])
        self.out.fillBranch("nvetoMuon", len(Muon_addlep_enu))
        self.out.fillBranch(
            "index_gVetoMuons", [
                x[0] for x in Muon_addlep_enu])
        return True


def call_postpoc(files):
    nameStrip = files.strip()
    filename = (nameStrip.split('/')[-1]).split('.')[-2]
    print(filename)

    if (search("Run", filename)):
        print(("This is a ", args.year, " Data file = ", filename))
        def mainModule(): return Standalone_AdditionalLepVetoFlag(args.year, True)
    else:
        print(("This is a ", args.year, " MC file = ", filename))
        def mainModule(): return Standalone_AdditionalLepVetoFlag(args.year, False)

    p = PostProcessor(
        args.tempStore,
        [files],
        cut=None,
        branchsel=None,
        modules=[
            mainModule()],
        postfix="",
        noOut=False,
        outputbranchsel=None)
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
        description='I had used VLoose selection for HPS taus - which didnot have SF computed. This script is used to select events after tighenting HPS selection to Loose - after channel assignement. Not Ideal. StopGap measure')
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
        '--year',
        '-y',
        help='specify the run - to make sure right triggers are used',
        choices=[
            '2016',
            '2016APV',
            '2017',
            '2018'])
    parser.add_argument(
        '--tempStore',
        '-t',
        help='Temporary staging area for files before moving out to hdfs',
        required=True)

    args = parser.parse_args()

    # Input file location typically:
    # fnames = glob.glob(args.inputLocation +
    # "/RadionTohhTohtatahbb_narrow_M-1000*.root")  #making a list of input
    # MC/DATA files
    # making a list of input MC/DATA files
    fnames = glob.glob(args.inputLocation + "/Radion*.root")

    outputDir = args.outputLocation

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
