#!/usr/bin/env python3
import correctionlib

# 1. Load the correction set
evaluator = correctionlib.CorrectionSet.from_file(
    "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/jetvetomaps.json"
)

# 2. Print all available top-level keys
print("\n🔑 All correction keys found:\n")
for key in evaluator.keys():
    print(" -", key)

# 3. Inspect each correction's internal structure
print("\n📊 Detailed correction structure:\n")
for key in evaluator.keys():
    corr = evaluator[key]
    print(f"Correction: {key}")
    print(f"  - Description: {corr.description}")
    print(f"  - Version: {corr.version}")
    print(f"  - Inputs:")
    for inp in corr.inputs:
        print(f"     • {inp.name} (type: {inp.type})")
    print(f"  - Output type: {corr.output}\n")

# 4. Safely access a specific correction (optional)
target_key = "Summer24Prompt24_RunBCDEFGHI_V1"  # ✅ use a key that exists
if target_key in evaluator.keys():
    corr = evaluator[target_key]
    print(f"\n✅ Accessed correction: {target_key}")
    print(f"   Description: {corr.description}")
    print(f"   Inputs: {[inp.name for inp in corr.inputs]}")
else:
    print(f"\n❌ Correction '{target_key}' not found! Check the list above.")
