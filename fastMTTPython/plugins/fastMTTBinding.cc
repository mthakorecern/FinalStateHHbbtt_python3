//Andrew Loeliger
//Quick pybind11 binding to the FastMTT mass calculation algorithm
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
//#include "TMatrixD.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <map>

using namespace classic_svFit;


LorentzVector getBestHiggsVector(			   
				 double measuredMETx, //met x component
			   double measuredMETy, //met y component
			   double covMETxx, //met cov 0,0 entry
			   double covMETxy, //met cov 0,1 and 1,0 entry
			   double covMETyy, //met cov 1,1 entry
			   int firstLeptonDecayCode,
			   double firstLeptonPt,
			   double firstLeptonEta,
			   double firstLeptonPhi,
			   double firstLeptonM,
			   int firstLeptonTauDecayMode,
			   int secondLeptonDecayCode,
			   double secondLeptonPt,
			   double secondLeptonEta,
			   double secondLeptonPhi,
			   double secondLeptonM,
			   int secondTauLeptonDecayMode)
{
  //Our input is a bit unorganized, let's spend some time putting that back together
  //first, let's construct the covariance matrix we're going to need
  TMatrixD covMET(2, 2);
  covMET[0][0] = covMETxx;
  covMET[1][0] = covMETxy;
  covMET[0][1] = covMETxy;
  covMET[1][1] = covMETyy;
  
  //Okay, now we need to construct the measured lepton objects to be used by the algorithm
  //let's first make a map between the decay code, and the SVFit decay enum

  std::map <int, MeasuredTauLepton::kDecayType> decayMapping { {0, MeasuredTauLepton::kUndefinedDecayType}, {1, MeasuredTauLepton::kTauToHadDecay}, {2, MeasuredTauLepton::kTauToMuDecay}, {3, MeasuredTauLepton::kTauToElecDecay} };

  std::vector<MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(
			       MeasuredTauLepton(decayMapping[firstLeptonDecayCode], 
						 firstLeptonPt,
						 firstLeptonEta,
						 firstLeptonPhi,
						 firstLeptonM,
						 firstLeptonTauDecayMode)
			       );
  measuredTauLeptons.push_back(
			       MeasuredTauLepton(decayMapping[secondLeptonDecayCode],
						 secondLeptonPt,
						 secondLeptonEta,
						 secondLeptonPhi,
						 secondLeptonM,
						 secondTauLeptonDecayMode)
			       );
  
  //The data is assembled, now we run fastMTT
  FastMTT aFastMTTAlgo;
  aFastMTTAlgo.run(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
  
  return ttP4;

}

double retrieveFastMTTMass(
			   double measuredMETx, //met x component
			   double measuredMETy, //met y component
			   double covMETxx, //met cov 0,0 entry
			   double covMETxy, //met cov 0,1 and 1,0 entry
			   double covMETyy, //met cov 1,1 entry
			   int firstLeptonDecayCode,
			   double firstLeptonPt,
			   double firstLeptonEta,
			   double firstLeptonPhi,
			   double firstLeptonM,
			   int firstLeptonTauDecayMode,
			   int secondLeptonDecayCode,
			   double secondLeptonPt,
			   double secondLeptonEta,
			   double secondLeptonPhi,
			   double secondLeptonM,
			   int secondTauLeptonDecayMode)
{
  return getBestHiggsVector(measuredMETx, measuredMETy, covMETxx, covMETxy, covMETyy, firstLeptonDecayCode, firstLeptonPt, firstLeptonEta, firstLeptonPhi, firstLeptonM, firstLeptonTauDecayMode, secondLeptonDecayCode, secondLeptonPt, secondLeptonEta, secondLeptonPhi, secondLeptonM, secondTauLeptonDecayMode).M();
}

double retrieveFastMTTPt(
			 double measuredMETx, //met x component
			 double measuredMETy, //met y component
			 double covMETxx, //met cov 0,0 entry
			 double covMETxy, //met cov 0,1 and 1,0 entry
			 double covMETyy, //met cov 1,1 entry
			 int firstLeptonDecayCode,
			 double firstLeptonPt,
			 double firstLeptonEta,
			 double firstLeptonPhi,
			 double firstLeptonM,
			 int firstLeptonTauDecayMode,
			 int secondLeptonDecayCode,
			 double secondLeptonPt,
			 double secondLeptonEta,
			 double secondLeptonPhi,
			 double secondLeptonM,
			 int secondTauLeptonDecayMode
			 )
{
  return getBestHiggsVector(measuredMETx, measuredMETy, covMETxx, covMETxy, covMETyy, firstLeptonDecayCode, firstLeptonPt, firstLeptonEta, firstLeptonPhi, firstLeptonM, firstLeptonTauDecayMode, secondLeptonDecayCode, secondLeptonPt, secondLeptonEta, secondLeptonPhi, secondLeptonM, secondTauLeptonDecayMode).Pt();
}

double retrieveFastMTTPhi(
			 double measuredMETx, //met x component
			 double measuredMETy, //met y component
			 double covMETxx, //met cov 0,0 entry
			 double covMETxy, //met cov 0,1 and 1,0 entry
			 double covMETyy, //met cov 1,1 entry
			 int firstLeptonDecayCode,
			 double firstLeptonPt,
			 double firstLeptonEta,
			 double firstLeptonPhi,
			 double firstLeptonM,
			 int firstLeptonTauDecayMode,
			 int secondLeptonDecayCode,
			 double secondLeptonPt,
			 double secondLeptonEta,
			 double secondLeptonPhi,
			 double secondLeptonM,
			 int secondTauLeptonDecayMode
			 )
{
  return getBestHiggsVector(measuredMETx, measuredMETy, covMETxx, covMETxy, covMETyy, firstLeptonDecayCode, firstLeptonPt, firstLeptonEta, firstLeptonPhi, firstLeptonM, firstLeptonTauDecayMode, secondLeptonDecayCode, secondLeptonPt, secondLeptonEta, secondLeptonPhi, secondLeptonM, secondTauLeptonDecayMode).Phi();
}

double retrieveFastMTTEta(
			 double measuredMETx, //met x component
			 double measuredMETy, //met y component
			 double covMETxx, //met cov 0,0 entry
			 double covMETxy, //met cov 0,1 and 1,0 entry
			 double covMETyy, //met cov 1,1 entry
			 int firstLeptonDecayCode,
			 double firstLeptonPt,
			 double firstLeptonEta,
			 double firstLeptonPhi,
			 double firstLeptonM,
			 int firstLeptonTauDecayMode,
			 int secondLeptonDecayCode,
			 double secondLeptonPt,
			 double secondLeptonEta,
			 double secondLeptonPhi,
			 double secondLeptonM,
			 int secondTauLeptonDecayMode
			 )
{
  return getBestHiggsVector(measuredMETx, measuredMETy, covMETxx, covMETxy, covMETyy, firstLeptonDecayCode, firstLeptonPt, firstLeptonEta, firstLeptonPhi, firstLeptonM, firstLeptonTauDecayMode, secondLeptonDecayCode, secondLeptonPt, secondLeptonEta, secondLeptonPhi, secondLeptonM, secondTauLeptonDecayMode).Eta();
}
//std::tuple<double, double, double, double> retrieveFastMTTFourVector(
std::vector<double> retrieveFastMTTFourVector(
			 double measuredMETx, //met x component
			 double measuredMETy, //met y component
			 double covMETxx, //met cov 0,0 entry
			 double covMETxy, //met cov 0,1 and 1,0 entry
			 double covMETyy, //met cov 1,1 entry
			 int firstLeptonDecayCode,
			 double firstLeptonPt,
			 double firstLeptonEta,
			 double firstLeptonPhi,
			 double firstLeptonM,
			 int firstLeptonTauDecayMode,
			 int secondLeptonDecayCode,
			 double secondLeptonPt,
			 double secondLeptonEta,
			 double secondLeptonPhi,
			 double secondLeptonM,
			 int secondTauLeptonDecayMode
			 )
{
  //return {2.,3.0,4.0,6.0};
  LorentzVector FV = getBestHiggsVector(measuredMETx, measuredMETy, covMETxx, covMETxy, covMETyy, firstLeptonDecayCode, firstLeptonPt, firstLeptonEta, firstLeptonPhi, firstLeptonM, firstLeptonTauDecayMode, secondLeptonDecayCode, secondLeptonPt, secondLeptonEta, secondLeptonPhi, secondLeptonM, secondTauLeptonDecayMode);
  return{FV.Pt(),FV.eta(),FV.Phi(),FV.M()};
  //return std::make_tuple(FV.Pt(),FV.eta(),FV.Phi(),FV.M());
  //return getBestHiggsVector(measuredMETx, measuredMETy, covMETxx, covMETxy, covMETyy, firstLeptonDecayCode, firstLeptonPt, firstLeptonEta, firstLeptonPhi, firstLeptonM, firstLeptonTauDecayMode, secondLeptonDecayCode, secondLeptonPt, secondLeptonEta, secondLeptonPhi, secondLeptonM, secondTauLeptonDecayMode);
}

PYBIND11_MODULE(pluginfastMTT_binding, m)
{
  m.doc() = "A binding to call fast mtt evaluations from python";
  m.def("fastMTTmass", &retrieveFastMTTMass, "get the mass from the fast MTT Calculation");
  m.def("fastMTTpt", &retrieveFastMTTPt, "get the pt from the fast MTT Calculation");
  m.def("fastMTTphi", &retrieveFastMTTPhi, "get the phi from the fast MTT Calculation");
  m.def("fastMTTeta", &retrieveFastMTTEta, "get the eta from the fast MTT Calculation");
  m.def("fastMTTfourvector", &retrieveFastMTTFourVector, "get the full Higgs four vector from fast MTT Calculation");
}
