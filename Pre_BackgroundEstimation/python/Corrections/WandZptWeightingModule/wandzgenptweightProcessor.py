from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from FinalStateHHbbtt.Pre_BackgroundEstimation.Corrections.WandZptWeightingModule.kFactorTool import *
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import argparse
import multiprocessing as np
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


# from kFactorTool import *


class wandzgenptweight(Module):
    def __init__(self, filename, year=2016, isData=False):
        print("ACTIAVTE: wandzgenptweight")
        self.filename = filename
        self.year = year
        self.isData = isData
        self.isMC = not self.isData
        self.kFactorTool = KFactorTool(year=self.year)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # Define output branches to check if preselection is working
        self.eventCount = 0
        self.out = wrappedOutputTree
        self.out.branch("ewkWWeight", "F")
        self.out.branch("ewkZWeight", "F")
        self.out.branch("ewkWDeborahsWeight", "F")
        self.out.branch("ewkZDeborahsWeight", "F")
        self.out.branch("qcdWWeight", "F")
        self.out.branch("qcdZTo2LWeight", "F")
        self.out.branch("GenV_pt", "F")
        self.out.branch("WnnloWeight", "F")
        self.out.branch("ZnnloWeight", "F")
        self.out.branch("combinedWZgenPtWeight", "F")
        self.out.branch("combinedWZgenPtDeborahWeight", "F")
        self.out.branch("combinedWZgenPtDylanWeight", "F")

        # Systematics - QCD Scale Factors
        self.out.branch("qcdWWeightRenUp", "F")
        self.out.branch("qcdWWeightRenDown", "F")
        self.out.branch("qcdWWeightFacUp", "F")
        self.out.branch("qcdWWeightFacDown", "F")
        self.out.branch("qcdZTo2LWeightRenUp", "F")
        self.out.branch("qcdZTo2LWeightRenDown", "F")
        self.out.branch("qcdZTo2LWeightFacUp", "F")
        self.out.branch("qcdZTo2LWeightFacDown", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        ewkWWeight = ewkZWeight = qcdWWeight = qcdZTo2LWeight = ewkWDeborahsWeight = ewkZDeborahsWeight = combinedWZgenPtDylanWeight = WnnloWeight = ZnnloWeight = combinedWZgenPtWeight = combinedWZgenPtDeborahWeight = 1
        qcdWWeightRenUp = qcdWWeightRenDown = qcdWWeightFacUp = qcdWWeightFacDown = 1
        qcdZTo2LWeightRenUp = qcdZTo2LWeightRenDown = qcdZTo2LWeightFacUp = qcdZTo2LWeightFacDown = 1
        GenV_pt = -9
        GenV = []
        if (search("DYJetsToLL_M-50", self.filename)):
            genParticles = Collection(event, "GenPart")
            GenV = list([gen for gen in genParticles if (
                gen.pdgId == 23) and (gen.status == 22)])
            if len(GenV) > 0:
                GenV_pt = GenV[0].pt
                ewkZWeight *= self.kFactorTool.getEWKZ(GenV_pt)
                ewkZDeborahsWeight *= self.kFactorTool.getDeborahEWKZ(GenV_pt)
                qcdZTo2LWeight *= self.kFactorTool.getQCDZTo2L(GenV_pt)
                qcdZTo2LWeightRenUp *= self.kFactorTool.getRenUpZTo2L(GenV_pt)
                qcdZTo2LWeightRenDown *= self.kFactorTool.getRenDownZTo2L(
                    GenV_pt)
                qcdZTo2LWeightFacUp *= self.kFactorTool.getFacUpZTo2L(GenV_pt)
                qcdZTo2LWeightFacDown *= self.kFactorTool.getFacDownZTo2L(
                    GenV_pt)
                ZnnloWeight = 0.934

                if (self.year == "2016APV"):
                    combinedWZgenPtWeight = qcdZTo2LWeight * ewkZWeight
                    combinedWZgenPtDeborahWeight = qcdZTo2LWeight * ewkZDeborahsWeight
                elif (self.year == "2016"):
                    combinedWZgenPtWeight = qcdZTo2LWeight * ewkZWeight
                    combinedWZgenPtDeborahWeight = qcdZTo2LWeight * ewkZDeborahsWeight
                elif (self.year == "2017"):
                    combinedWZgenPtWeight = qcdZTo2LWeight * ewkZWeight * 0.934
                    combinedWZgenPtDeborahWeight = qcdZTo2LWeight * ewkZDeborahsWeight * 0.934
                elif (self.year == "2018"):
                    combinedWZgenPtWeight = qcdZTo2LWeight * ewkZWeight * 0.934
                    combinedWZgenPtDeborahWeight = qcdZTo2LWeight * ewkZDeborahsWeight * 0.934

            combinedWZgenPtDylanWeight = 1.23
            # combinedWZgenPtWeight = qcdZTo2LWeight*ewkZWeight*0.934
            # combinedWZgenPtDeborahWeight = qcdZTo2LWeight*ewkZDeborahsWeight*0.934

        elif (search("WJet", self.filename)):
            genParticles = Collection(event, "GenPart")
            GenV = list([gen for gen in genParticles if (
                abs(gen.pdgId) == 24) and gen.status == 22])
            if len(GenV) > 0:
                GenV_pt = GenV[0].pt
                ewkWWeight *= self.kFactorTool.getEWKW(GenV_pt)
                ewkWDeborahsWeight *= self.kFactorTool.getDeborahEWKW(GenV_pt)
                qcdWWeight *= self.kFactorTool.getQCDW(GenV_pt)
                qcdWWeightRenUp *= self.kFactorTool.getRenUpW(GenV_pt)
                qcdWWeightRenDown *= self.kFactorTool.getRenDownW(GenV_pt)
                qcdWWeightFacUp *= self.kFactorTool.getFacUpW(GenV_pt)
                qcdWWeightFacDown *= self.kFactorTool.getFacDownW(GenV_pt)
                WnnloWeight = 0.9135
                if (self.year == "2016APV"):
                    combinedWZgenPtWeight = qcdWWeight * ewkWWeight
                    # VERIFY again..Deborah said QCD to be taken from python
                    # for 2016 but we have all same tunes unlike victor
                    combinedWZgenPtDeborahWeight = qcdWWeight * ewkWDeborahsWeight
                elif (self.year == "2016"):
                    combinedWZgenPtWeight = qcdWWeight * ewkWWeight
                    # VERIFY again..Deborah said QCD to be taken from python
                    # for 2016 but we have all same tunes unlike victor
                    combinedWZgenPtDeborahWeight = qcdWWeight * ewkWDeborahsWeight
                elif (self.year == "2017"):
                    combinedWZgenPtWeight = qcdWWeight * ewkWWeight * 0.9135
                    combinedWZgenPtDeborahWeight = qcdWWeight * ewkWDeborahsWeight * 0.9135
                elif (self.year == "2018"):
                    combinedWZgenPtWeight = qcdWWeight * ewkWWeight * 0.9135
                    combinedWZgenPtDeborahWeight = qcdWWeight * ewkWDeborahsWeight * 0.9135

            combinedWZgenPtDylanWeight = 1.21
            # combinedWZgenPtWeight = qcdWWeight*ewkWWeight*0.9135
            # combinedWZgenPtDeborahWeight = qcdWWeight*ewkWDeborahsWeight*0.9135

        self.out.fillBranch("ewkWWeight", ewkWWeight)
        self.out.fillBranch("ewkZWeight", ewkZWeight)
        self.out.fillBranch("ewkWDeborahsWeight", ewkWDeborahsWeight)
        self.out.fillBranch("ewkZDeborahsWeight", ewkZDeborahsWeight)
        self.out.fillBranch("qcdWWeight", qcdWWeight)
        self.out.fillBranch("qcdZTo2LWeight", qcdZTo2LWeight)
        self.out.fillBranch("GenV_pt", GenV_pt)
        self.out.fillBranch("WnnloWeight", WnnloWeight)
        self.out.fillBranch("ZnnloWeight", ZnnloWeight)
        self.out.fillBranch("combinedWZgenPtWeight", combinedWZgenPtWeight)
        self.out.fillBranch(
            "combinedWZgenPtDeborahWeight",
            combinedWZgenPtDeborahWeight)
        self.out.fillBranch(
            "combinedWZgenPtDylanWeight",
            combinedWZgenPtDylanWeight)
        # Systematics - QCD Scale Factors
        self.out.fillBranch("qcdWWeightRenUp", qcdWWeightRenUp)
        self.out.fillBranch("qcdWWeightRenDown", qcdWWeightRenDown)
        self.out.fillBranch("qcdWWeightFacUp", qcdWWeightFacUp)
        self.out.fillBranch("qcdWWeightFacDown", qcdWWeightFacDown)
        self.out.fillBranch("qcdZTo2LWeightRenUp", qcdZTo2LWeightRenUp)
        self.out.fillBranch("qcdZTo2LWeightRenDown", qcdZTo2LWeightRenDown)
        self.out.fillBranch("qcdZTo2LWeightFacUp", qcdZTo2LWeightFacUp)
        self.out.fillBranch("qcdZTo2LWeightFacDown", qcdZTo2LWeightFacDown)

        return True


def call_postpoc(files):
    nameStrip = files.strip()
    filename = (nameStrip.split('/')[-1]).split('.')[-2]
    print(filename)

    if (search("Run", filename)):
        print(("This is a ", args.year, " Data file = ", filename))
        def mainModule(): return wandzgenptweight(filename, args.year, True)
    else:
        print(("This is a ", args.year, " MC file = ", filename))
        def mainModule(): return wandzgenptweight(filename, args.year, False)

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
    parser = argparse.ArgumentParser(description='Add W and Z genpt weight')
    parser.add_argument(
        '--inputLocation',
        '-i',
        help="enter the path to the location of input file set",
        default="")
    parser.add_argument(
        '--outputLocation',
        '-o',
        help="enter the path where you want the output files to be stored",
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

    # fnames = list(set(glob.glob(args.inputLocation +
    # "/*.root"))-set(glob.glob(args.inputLocation + "/*MET*Run*.root")))
    # #making a list of input MC/DATA files
    fnames = glob.glob(args.inputLocation + "/*MET*Run*.root")
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
            pool = np.Pool(int(args.ncores))
            print(("list", argList))
            res = pool.map(call_postpoc, argList)
        except Exception as error:
            print("MultiProc error - needs debugging")
            print(("An exception occurred:", type(error).__name__))
