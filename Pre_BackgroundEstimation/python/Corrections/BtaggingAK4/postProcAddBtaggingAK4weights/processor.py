import os, sys
import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
import glob
import numpy as np
from re import search
import multiprocessing as  np
import argparse
import math



from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from FinalStateHHbbtt.Pre_BackgroundEstimation.Corrections.BtaggingAK4.BTaggingTool import *

class addbtagwts(Module):
    def __init__(self, year=2016, isData=False):
        self.year = year
        self.isData = isData
        self.isMC = not self.isData
        if self.isMC:
            self.btagTool_loose = BTagWeightTool(tagger="deepjetflavb",wp='loose', sigmabc='central', sigmalight='central',year=self.year)
            self.btagTool_medium = BTagWeightTool(tagger="deepjetflavb",wp='medium', sigmabc='central', sigmalight='central',year=self.year)
            self.btagTool_tight = BTagWeightTool(tagger="deepjetflavb",wp='tight', sigmabc='central', sigmalight='central',year=self.year)
           
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        #Define output branches to check if preselection is working
        self.out = wrappedOutputTree
        self.out.branch("bjetWeight_medium","F")
        self.out.branch("bjetWeight_loose", "F")
        self.out.branch("bjetWeight_tight", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets = Collection(event, "Jet")

        #Filter the jet collection to only have jets interesting for the analysis
        jetsgood_loose = []
        jetsgood_medium = []
        jetsgood_tight = []

        if self.isData:
            self.out.fillBranch("bjetWeight_medium",1.0)
            self.out.fillBranch("bjetWeight_loose",1.0)
            self.out.fillBranch("bjetWeight_tight",1.0)
            return True            


        for i in range(event.ngood_LooseJets):
            jetsgood_loose.append(jets[event.index_gLooseJets[i]])

        for i in range(event.ngood_MediumJets):
            jetsgood_medium.append(jets[event.index_gMediumJets[i]])

        for i in range(event.ngood_TightJets):
            jetsgood_tight.append(jets[event.index_gTightJets[i]])

        
        print("Number of good jets : ", event.ngood_Jets)
        print("Number of Loose jets : ", event.ngood_LooseJets," weight : ",self.btagTool_loose.getWeight(jetsgood_loose))
        print("Number of Medium jets : ", event.ngood_MediumJets," weight : ",self.btagTool_medium.getWeight(jetsgood_medium))
        print("Number of Tight jets : ", event.ngood_TightJets," weight : ",self.btagTool_tight.getWeight(jetsgood_tight))

        self.out.fillBranch("bjetWeight_loose",self.btagTool_loose.getWeight(jetsgood_loose))
        self.out.fillBranch("bjetWeight_medium",self.btagTool_medium.getWeight(jetsgood_medium))
        self.out.fillBranch("bjetWeight_tight",self.btagTool_tight.getWeight(jetsgood_tight))

        return True        


def call_postpoc(files):
	nameStrip=files.strip()
	filename = (nameStrip.split('/')[-1]).split('.')[-2]
	print(filename)
	
	if (search("Run", filename)):
		print ("This is a ",args.year," Data file = ",filename)
		mainModule = lambda: addbtagwts(args.year,True)
	else:
		print ("This is a ",args.year," MC file = ", filename)
		mainModule = lambda: addbtagwts(args.year,False)
	
	p = PostProcessor(args.tempStore,[files], cut=None, branchsel=None,modules=[mainModule()], postfix="",noOut=False,outputbranchsel=None)
	p.run()	
	#print("###############MOVING THE OUTPUT FILE BACK TO HDFS#######################")
	#os.system("hadoop fs -moveFromLocal -f "+args.tempStore+"/"+filename+".root"+" "+outputDir+"/.") #This currently doesnot work due to username differences - it takes parida by default

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Analysis framwork code part 2. Apply object selections, and get event categorization. FastMTT branches')	
	parser.add_argument('--inputLocation','-i',help="enter the path to the location of input file set",default="")
	parser.add_argument('--outputLocation','-o',help="enter the path where yu want the output files to be stored",default =".")
	parser.add_argument('--ncores','-n',help ="number of cores for parallel processing", default=1)
	parser.add_argument('--year','-y',help='specify the run - to make sure right triggers are used',choices=['2016','2016APV','2017','2018'])
	parser.add_argument('--tempStore','-t',help='Temporary staging area for files before moving out to hdfs', required=True)

	args = parser.parse_args()

	fnames = glob.glob(args.inputLocation + "/TTToSemiLeptonic*.root")  #making a list of input MC/DATA files
	outputDir = args.outputLocation

	argList = list()
	for file in fnames:
		argList.append(file)

	if int(args.ncores) == 1:
		for arr in argList:
			print ("Using a single thread ")
			call_postpoc(arr)
	
	else:
		try:
			print ("Using Multithreading")
			pool = np.Pool(int(args.ncores))
			print ("list", argList)
			res=pool.map(call_postpoc, argList)
		except Exception as error:
			print ("MultiProc error - needs debugging")
			print("An exception occurred:", type(error).__name__)