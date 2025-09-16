import correctionlib
evaluator = correctionlib.CorrectionSet.from_file("/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/jet_jerc.json")

# 2) One per line
for k in evaluator.keys():
    print(k)

