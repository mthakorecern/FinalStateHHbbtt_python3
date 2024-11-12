import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import argparse
import glob
from re import search

class TauToElectronFakeRate(Module):
	def __init__(self,year):
		self.year = year
		self.tauFV= ROOT.TLorentzVector(0.0,0.0,0.0,0.0)
		self.higgsBBFV = ROOT.TLorentzVector(0.0,0.0,0.0,0.0)


	def beginJob(self):
		self.denominator_count = 0.0
		self.numerator_count = 0.0
		#pass

	def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		self.denominator_count_file = 0.0
		self.numerator_count_file = 0.0


	def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
		print("\n\n\n>>>>>>>>>>>...END FILE........<<<<<<<<<<<<")
		print("Numerator = ", self.numerator_count_file," Denomiantor = ",self.denominator_count_file)
		if (self.denominator_count_file==0):
			print("Denominator is zezo")
		else:
			print("Fake rate = %.3f%%" % (float(self.numerator_count_file * 100) / self.denominator_count_file))
		#print("Fake rate = ", float((self.numerator_count_file*100)/self.denominator_count_file))		

	def endJob(self):
		print("\n\n\n>>>>>>>>>>>...Full DY Process........<<<<<<<<<<<<")
		print("Numerator = ", self.numerator_count," Denomiantor = ",self.denominator_count)
		print("Fake rate = %.3f%%" % (float(self.numerator_count * 100) / self.denominator_count))
		#pass


	def analyze(self, event):
		def has_lepton_descendant(genparts, tau, target_pdgId):
			""" Recursively checks if a tau has a specific lepton (e.g., electron or muon) as a descendant. """
			for gen in genparts:
				if gen.genPartIdxMother == tau._index:
					if abs(gen.pdgId) == target_pdgId:
						return gen  # Return the matching lepton
					# Recursively check the next level in the chain
					descendant = has_lepton_descendant(genparts, gen, target_pdgId)
					if descendant:
						return descendant
			return None

		def FatJetTauOverlap(tauObject_enu): # Function used by Tau(boostedTau) cleaning vs AK8 > 0.8
			self.tauFV.SetPtEtaPhiM(tauObject_enu.pt,tauObject_enu.eta,tauObject_enu.phi,tauObject_enu.mass)
			deltaR = self.tauFV.DeltaR(self.higgsBBFV)
			if deltaR > 1.5:
				return True
			else:
				return False

		if (event.METcorrected_pt < 180):
			return False
		FatJet = Collection(event, "FatJet")
		FatJet_enu = filter(lambda x: (x.pt_nom >= 200) and (abs(x.eta) < 2.5) and (x.jetId>1) and (x.msoftdrop_nom>=30),FatJet)
		if (len(FatJet_enu)==0):
			return False
		
		self.higgsBBFV.SetPtEtaPhiM(FatJet_enu[0].pt_nom,FatJet_enu[0].eta,FatJet_enu[0].phi,FatJet_enu[0].mass_nom)

		taus = Collection(event, "Tau")
		boostedtaus = Collection(event, "boostedTau")
		Tau_enu = filter(lambda x: (x.pt > 20) and (abs(x.eta) < 2.5) and (abs(x.dz)<0.2) and (x.idDecayModeNewDMs) and (x.idDeepTau2018v2p5VSjet >= 4) and (x.idDeepTau2018v2p5VSe >= 2) and (x.idDeepTau2018v2p5VSmu >= 1),taus) #The HPS tau id are not bitmps anymore
		boostedTau_enu = filter(lambda x: (x.pt > 20) and (abs(x.eta) < 2.5) and (x.rawDeepTau2018v2p7VSjet>=0.85), boostedtaus)

		Tau_enu = filter(FatJetTauOverlap,Tau_enu)
		boostedTau_enu = filter(FatJetTauOverlap, boostedTau_enu)

		if ((len(Tau_enu)==0) and (len(boostedTau_enu)==0)):
			return False

		# Access gen particles collection
		genparts = Collection(event, "GenPart")

		# Select "original" tau leptons with no tau parent
		taus = []
		for gen in genparts:
		    if abs(gen.pdgId) == 15:
		        # Check if this tau has another tau as its parent
		        has_tau_parent = False
		        if gen.genPartIdxMother >= 0:
		            parent = genparts[gen.genPartIdxMother]
		            if abs(parent.pdgId) == 15:
		                has_tau_parent = True
		        # If no tau parent, consider this an "original" tau
		        if not has_tau_parent:
		            taus.append(gen)
		if len(taus) != 2:
			#print ("Number of taus in the event = ", len(taus))
			return False  # Skip events without exactly two tau leptons

		# Check for tau decays: one to an electron, the other hadronic
		electron_from_tau = None
		hadronic_tau = None

		for tau in taus:
			# Check if tau has an electron descendant
			electron_descendant = has_lepton_descendant(genparts, tau, 11)  # 11 for electron
			if electron_descendant:
				electron_from_tau = electron_descendant
			else:
				# Check for any lepton descendant to classify as hadronic or leptonic
				if not has_lepton_descendant(genparts, tau, 13):  # 13 for muon
					hadronic_tau = tau

		# Count denominator if we have one tau decaying to an electron and the other hadronic
		if electron_from_tau and hadronic_tau:
			# Apply event-level pre-selections both for numerator and Denominator

			#print ("Filling Denominator")
			self.denominator_count += 1*event.crossSectionWeighting*event.pileupWeighting
			self.denominator_count_file += 1*event.crossSectionWeighting*event.pileupWeighting


			# Check if the electron from tau decay matches any of the filtered taus or boosted taus
			matched = False
			for tau in Tau_enu:
				if electron_from_tau.p4().DeltaR(tau.p4()) <= 0.1:
					matched = True
					break
			for btau in boostedTau_enu:
				if electron_from_tau.p4().DeltaR(btau.p4()) <= 0.1:
					matched = True
					break

			if matched:
				#print ("Filling ###Numerator###")
				self.numerator_count += 1*event.crossSectionWeighting*event.pileupWeighting  # Increment numerator count
				self.numerator_count_file += 1*event.crossSectionWeighting*event.pileupWeighting
			return True
		else:
			return False


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Analysis framwork code part 2. Apply object selections, and get event categorization. FastMTT branches')	
	parser.add_argument('--inputLocation','-i',help="enter the path to the location of input file set",default="")
	parser.add_argument('--outputLocation','-o',help="enter the path where yu want the output files to be stored",default =".")
	parser.add_argument('--ncores','-n',help ="number of cores for parallel processing", default=1)
	parser.add_argument('--postfix',help="string at the end of output file names", default="")
	parser.add_argument('--year','-y',help='specify the run - to make sure right triggers are used',choices=['2016','2016APV','2017','2018'])
	parser.add_argument('--tempStore','-t',help='Temporary staging area for files before moving out to hdfs', required=True)
	parser.add_argument('--transferOff','-toff',help='with this flag the files will be not moved to hdfs from the temp area',action='store_true')
	parser.add_argument('--Remove','-r',help="Remove [old naming convention branches] crossSectionWeighting, pileupWeighting, pileupWeight_UP, pileupWeight_DOWN, FinalWeighting, FinalWeighting_pileupWeight_UP, FinalWeighting_pileupWeight_DOWN",action='store_true')

	args = parser.parse_args()

	fnames = glob.glob(args.inputLocation + "/DYJetsToLL_M-50*.root")  #making a list of input MC/DATA files
	outputDir = args.outputLocation

	mainModule = lambda: TauToElectronFakeRate(args.year)
	p = PostProcessor(args.tempStore,fnames, cut=None, branchsel=None,modules=[mainModule()], postfix="",noOut=True,outputbranchsel=None)
	p.run()

