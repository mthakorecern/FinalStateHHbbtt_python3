#!/usr/bin/env python3
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import argparse
import os
import shutil
import sys

# Ensure local python modules are visible on worker nodes
sys.path.append(os.path.join(os.environ["CMSSW_BASE"], "python"))
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

# Custom modules
from PhysicsTools.NATModules.modules.jetVetoMap import jetVMAP
from PhysicsTools.NATModules.modules.jetId import jetId
from PhysicsTools.NATModules.modules.jetCorr import jetJERC
from PhysicsTools.NATModules.modules.fatjetcorr import fatJetJERC


class EventCounter(Module):
    """Simple module to print # events before processing"""
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        print(">>> Events in this file before any skim:", inputTree.GetEntries())
    def analyze(self, event):
        return True


def process_file(inputFile, outputFile,
                 goldenjson, jetidjson, vetomapjson, jercjson,
                 isMC, jersmearjson=None):

    # ------------------------
    # Define 2024 preselection
    # ------------------------

    met_selection = ["PuppiMET_pt > 150"]

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


        
    """Run NanoAODTools modules on a single file."""
    modules = [EventCounter()]

    # JetID
    modules.append(jetId(jetidjson, jetType="AK4PUPPI"))
    print("Added JetID module")

    # JetVetoMap
    modules.append(jetVMAP(vetomapjson,
                           corrName="Summer24Prompt24_RunBCDEFGHI_V1",
                           veto_map_name="jetvetomap"))
    print("Added JetVetoMap module")

    if isMC:
        # AK4 jets
        modules.append(
            jetJERC(
                json_JERC=jercjson,
                json_JERsmear=None,
                L1Key="Summer24Prompt24_V1_MC_L1FastJet_AK4PFPuppi",
                L2Key="Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi",
                L3Key="Summer24Prompt24_V1_MC_L3Absolute_AK4PFPuppi",
                L2L3Key="Summer24Prompt24_V1_MC_L2L3Residual_AK4PFPuppi", 
                scaleTotalKey="Summer24Prompt24_V1_MC_Total_AK4PFPuppi",
                smearKey=None,
                JERKey="Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi",
                JERsfKey="Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi",
                overwritePt=False,
                usePhiDependentJEC=False,
                useRunDependentJEC=False,
                isMC=True
            )
        )
        # FatJets
        modules.append(
            fatJetJERC(
                json_JERC=jercjson,
                json_JERsmear=None,
                L1Key="Summer24Prompt24_V1_MC_L1FastJet_AK4PFPuppi",  # same JSON keys
                L2Key="Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi",
                L3Key="Summer24Prompt24_V1_MC_L3Absolute_AK4PFPuppi",
                L2L3Key="Summer24Prompt24_V1_MC_L2L3Residual_AK4PFPuppi", 
                scaleTotalKey="Summer24Prompt24_V1_MC_Total_AK4PFPuppi",
                smearKey=None,
                JERKey="Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi",
                JERsfKey="Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi",
                overwritePt=False,
                usePhiDependentJEC=False,
                useRunDependentJEC=False,
                isMC=True
            )
        )
        print("Added Jet/FatJet JERC (MC)")

    else:
        # AK4 jets
        modules.append(
            jetJERC(
                json_JERC=jercjson,
                json_JERsmear=None,
                L1Key="Summer24Prompt24_V1_DATA_L1FastJet_AK4PFPuppi",
                L2Key="Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi",
                L3Key="Summer24Prompt24_V1_DATA_L3Absolute_AK4PFPuppi",
                L2L3Key="Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi",
                scaleTotalKey=None,
                smearKey=None,
                JERKey=None,
                JERsfKey=None,
                overwritePt=False,
                usePhiDependentJEC=False,
                useRunDependentJEC=True,
                isMC=False
            )
        )
        # FatJets
        modules.append(
            fatJetJERC(
                json_JERC=jercjson,
                json_JERsmear=None,
                L1Key="Summer24Prompt24_V1_DATA_L1FastJet_AK4PFPuppi",
                L2Key="Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi",
                L3Key="Summer24Prompt24_V1_DATA_L3Absolute_AK4PFPuppi",
                L2L3Key="Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi",
                scaleTotalKey=None,
                smearKey=None,
                JERKey=None,
                JERsfKey=None,
                overwritePt=False,
                usePhiDependentJEC=False,
                useRunDependentJEC=True,
                isMC=False
            )
        )
        print("Added Jet/FatJet JERC (Data)")

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
    parser.add_argument("--jercjson", required=True, help="JERC JSON file")
    parser.add_argument("--jersmearjson", default=None, help="JER smearing JSON file (MC only)")
    parser.add_argument("--isMC", action="store_true", help="Flag: running on MC")

    args = parser.parse_args()

    print(">>> nanoPostProc_helper.py called with arguments:")
    for k, v in vars(args).items():
        print(f"    {k}: {v}")

    if not os.path.exists(args.inputFile) and not args.inputFile.startswith("root://"):
        print(f"[ERROR] Input file {args.inputFile} does not exist or is not accessible!")
    else:
        print(f"[INFO] Opening input file: {args.inputFile}")

    process_file(
        args.inputFile, args.outputFile,
        args.goldenjson, args.jetidjson, args.vetomapjson,
        args.jercjson, args.isMC, args.jersmearjson
    )
