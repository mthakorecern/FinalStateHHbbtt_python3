#!/usr/bin/env python3
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import argparse
import os
import shutil
import sys

sys.path.append(os.path.join(os.environ["CMSSW_BASE"], "python"))
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from PhysicsTools.NATModules.modules.jetId import jetId
from PhysicsTools.NATModules.modules.fatjetId import fatJetId
from PhysicsTools.NATModules.modules.jetVetoMap import jetVMAP
from PhysicsTools.NATModules.modules.fatJetvetoMap import fatJetVMAP
from PhysicsTools.NATModules.modules.applyJercFJERC import ApplyJercAll  


# class EventCounter(Module):
#     def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#         print(">>> Events in this file before any skim:", inputTree.GetEntries())
#     def analyze(self, event):
#         return True


def process_file(inputFile, outputFile,
                 goldenjson, jetidjson, vetomapjson, jercjson, isMC):

    met_selection = ["PuppiMET_pt >= 120"]

    eventSelection2024 = [
        "Flag_goodVertices",
        "Flag_globalSuperTightHalo2016Filter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter",
        "Flag_BadPFMuonDzFilter",
        "Flag_hfNoisyHitsFilter",
        "Flag_eeBadScFilter",
        "Flag_ecalBadCalibFilter",
        "PV_npvsGood > 0"
    ]

    trigger_2024 = [
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF"
    ]

    # Build preselection string
    preselection = (
        "(" + "&&".join(met_selection) + ")" +
        "&&(" + "&&".join(eventSelection2024) + ")" +
        "&&(" + "||".join(trigger_2024) + ")"
    )

    print(">>> Preselection string:", preselection)

    modules = []

    modules.append(jetId("jetid.json", jetType="AK4PUPPI"))
    print("Added JetID module")
    
    modules.append(fatJetId("jetid.json", jetType="AK8PUPPI"))
    print("Added FatJetID module")

    # JERC (Jet Energy Corrections + JER smearing + MET)
    modules.append(
        ApplyJercAll(
            year="2024",        # could expose as CLI arg if you want
            isData=not isMC,
            jercjson=jercjson,
            era=None,
            year_unc="2024"
        )
    )
    print("Added ApplyJercAll module")

    modules.append(jetVMAP("jetvetomaps.json",
                           corrName="Summer24Prompt24_RunBCDEFGHI_V1",
                           veto_map_name="jetvetomap"))
    print("Added JetVetoMap module")

    modules.append(fatJetVMAP("jetvetomaps.json",
                           corrName="Summer24Prompt24_RunBCDEFGHI_V1",
                           veto_map_name="jetvetomap"))
    print("Added FatJetVetoMap module")

    # PostProcessor
    if isMC:
        print(">>> [INFO] Running on MC: golden JSON will be ignored")
    else:
        print(f">>> [INFO] Running on Data: applying golden JSON {goldenjson}")

    # PostProcessor
    p = PostProcessor(
        outputDir=os.path.dirname(outputFile),
        inputFiles=[inputFile],
        cut=preselection,
        branchsel=None,
        modules=modules,
        postfix="",
        provenance=True,
        noOut=False,
        justcount=False,
        fwkJobReport=False,
        haddFileName=None,
        jsonInput=None if isMC else goldenjson
    )
    p.run()

    # Ensure output is renamed properly
    produced = os.path.join(os.path.dirname(outputFile), os.path.basename(inputFile))
    if os.path.exists(produced) and produced != outputFile:
        shutil.move(produced, outputFile)
    elif not os.path.exists(outputFile):
        raise RuntimeError(f"[ERROR] Expected output {outputFile} not produced!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NanoAODTools helper for farmout (1 input → 1 output).")
    parser.add_argument("--inputFile", required=True, help="Single input ROOT file")
    parser.add_argument("--outputFile", required=True, help="Output ROOT file")
    parser.add_argument("--goldenjson", required=False, help="Golden JSON (Data only)")
    parser.add_argument("--jetidjson", required=True, help="JetID JSON file")
    parser.add_argument("--vetomapjson", required=True, help="Jet veto map JSON file")
    parser.add_argument("--jercjson", required=True, help="JERC JSON file (jet_jerc.json[.gz])")
    parser.add_argument("--isMC", action="store_true", help="Flag: running on MC")

    args = parser.parse_args()

    print(">>> nanoPostProc_helper.py called with arguments:")
    for k, v in vars(args).items():
        print(f"{k}: {v}")

    # if not os.path.exists(args.inputFile) and not args.inputFile.startswith("root://"):
    #     print(f"[ERROR] Input file {args.inputFile} does not exist or is not accessible!")
    # else:
    #     print(f"[INFO] Opening input file: {args.inputFile}")

    process_file(
        args.inputFile, args.outputFile,
        args.goldenjson, args.jetidjson, args.vetomapjson, args.jercjson,
        args.isMC)
