#!/usr/bin/env python3
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import json
import re
import os
import sys
import glob
import argparse
import multiprocessing as mp
import subprocess
import shutil

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module 


# Importing veto map module
from PhysicsTools.NATModules.modules.jetVetoMap import jetVMAP ## New Module added from https://github.com/cms-cat/nanoAOD-tools-modules/blob/master/python/modules/jetVetoMap.py

from PhysicsTools.NATModules.modules.jetId import jetId ## New Module added from https://github.com/cms-cat/nanoAOD-tools-modules/blob/master/python/modules/jetVetoMap.py

from PhysicsTools.NATModules.modules.jetCorr import jetJERC ## New module added from https://github.com/cms-cat/nanoAOD-tools-modules/blob/master/python/modules/jetCorr.py 

## Custom FatJet correction Module
from PhysicsTools.NATModules.modules.fatjetcorr import fatJetJERC  # AK8 (new)



class EventCounter(Module):
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        print(">>> Events in this file before any skim:", inputTree.GetEntries())
    def analyze(self, event):
        return True

def chunk_processor(args):
    (chunk_files,
     goldenjson,
     jetidjson,
     vetomapjson,
     jercjson,
     jersmearjson,
     isMC,
     outputDir,
     chunk_outfile,
     finalHadd) = args

    modules = []

    if not isMC:
        print(f"Processing Data with the Golden JSON: {goldenjson}")
    else:
        print("Processing MC (no Golden JSON applied).")
    
    modules.append(EventCounter())
    ## JET ID Branch Creation
    modules.append(jetId(jetidjson, jetType="AK4PUPPI"))
    print(f"Added JetID module")
  

    ## Applying JetVeto Maps
    modules.append(jetVMAP(vetomapjson,
                            corrName="Summer24Prompt24_RunBCDEFGHI_V1",
                            veto_map_name="jetvetomap"))
    print("Added JetVetoMap module")

    if isMC:
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
                overwritePt=True,
                usePhiDependentJEC=False,
                useRunDependentJEC=False,
                isMC=True
            )
        )
        print("Added AK4 JetJERC module (MC)")

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
                overwritePt=True,
                usePhiDependentJEC=False,
                useRunDependentJEC=False,
                isMC=True
            )
        )
        print("Added FatJetJERC module (MC)")
    else:
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
                overwritePt=True,
                usePhiDependentJEC=False,
                useRunDependentJEC=True,
                isMC=False
            )
        )
        print("Added AK4 JetJERC module (Data)")

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
                overwritePt=True,
                usePhiDependentJEC=False,
                useRunDependentJEC=True,
                isMC=False
            )
        )
        print("Added FatJetJERC module (Data)")

    # Always run in "normal output" mode
    p = PostProcessor(
        outputDir=os.path.dirname(chunk_outfile),
        inputFiles=chunk_files,
        cut=None,
        branchsel=None,
        modules=modules,
        postfix="",            # don’t add _Skim suffix
        provenance=True,
        noOut=False,           # ensure we actually write output
        justcount=False,
        fwkJobReport=False,
        haddFileName=None,     # avoid merge mode
        jsonInput=goldenjson if not isMC else None
    )
    p.run()

    # --- One-to-one: move immediately to HDFS ---
    if not finalHadd and chunk_outfile:
        hdfs_out = os.path.join(outputDir, os.path.basename(chunk_outfile))
        print(f">>> Moving {chunk_outfile} -> {hdfs_out}")
        os.makedirs(outputDir, exist_ok=True)
        try:
            os.rename(chunk_outfile, hdfs_out)   # atomic if same FS
        except OSError:
            shutil.move(chunk_outfile, hdfs_out) # fallback if cross-FS



def run_postproc(goldenjson, jetidjson, vetomapjson, jercjson, jersmearjson,
                 isMC, inputDir, outputDir, tempDir, finalOut, finalHadd=False):

    inputFiles = glob.glob(os.path.join(inputDir, "*.root"))
    if not inputFiles:
        raise RuntimeError(f"No ROOT files found in {inputDir}")
    print(f"Running with {len(inputFiles)} input files")

    os.makedirs(outputDir, exist_ok=True)

    # --- Case 1: Final hadd requested ---
    if finalHadd:
        os.makedirs(tempDir, exist_ok=True)
        chunks = [inputFiles[i::os.cpu_count()] for i in range(os.cpu_count()) if inputFiles[i::os.cpu_count()]]
        worker_args = []
        for chunk in chunks:
            # each file inside chunk keeps its original basename inside tempDir
            for infile in chunk:
                out_name = os.path.basename(infile)
                local_out = os.path.join(tempDir, out_name)
                worker_args.append(([infile], goldenjson, jetidjson, vetomapjson,
                                    jercjson, jersmearjson, isMC, tempDir, local_out,
                                    finalHadd))

        with mp.Pool(min(len(worker_args), os.cpu_count())) as pool:
            pool.map(chunk_processor, worker_args)

        # now just grab all .root files from tempDir (they have original names)
        chunk_files = sorted(glob.glob(os.path.join(tempDir, "*.root")))

        if not finalOut:
            dataset_name = os.path.basename(os.path.normpath(inputDir))
            finalOut = f"{dataset_name}.root"

        final_outfile = os.path.join(outputDir, os.path.basename(finalOut))
        print(f"Running final haddnano.py {final_outfile}")
        subprocess.check_call(["haddnano.py", final_outfile] + chunk_files)
        print(f"Final merged output available at: {final_outfile}")

        # cleanup
        for cf in chunk_files:
            os.remove(cf)
        os.rmdir(tempDir)


    # --- Case 2: One-to-one mapping (default) ---
    else:  # one-to-one mapping
        worker_args = []
        for infile in inputFiles:
            out_name = os.path.basename(infile)
            # write to temp first
            local_out = os.path.join(tempDir, out_name)
            worker_args.append(([infile], goldenjson, jetidjson, vetomapjson,
                                jercjson, jersmearjson, isMC, outputDir, local_out,
                                finalHadd))


        with mp.Pool(min(len(worker_args), os.cpu_count())) as pool:
            pool.map(chunk_processor, worker_args)

        print(f"Processed {len(inputFiles)} files → outputs in {outputDir}")

    print("All processing finished successfully.")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run NanoAODTools with GoldenJSON + JetVeto + JetCorr (AK4 & AK8).\n\n"
                    "Output modes:\n"
                    "  * Default (no --finalHadd): One-to-one mapping. Each input file produces\n"
                    "    one output file with the same name in --outputDir.\n"
                    "  * With --finalHadd: All input files are split into chunks, each chunk\n"
                    "    processed in parallel, and intermediate outputs merged into a single\n"
                    "    ROOT file (--finalOut). Temporary chunk files are cleaned up."
                    , formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--goldenjson", required=True, help="Golden JSON file (Data only)")
    parser.add_argument("--jetidjson", required=True, help="JetID JSON file")
    parser.add_argument("--vetomapjson", required=True, help="Jet veto map JSON file")
    parser.add_argument("--jercjson", required=True, help="JERC JSON file (AK4/AK8 Puppi)")
    parser.add_argument("--jersmearjson", required=None, help="JER smearing JSON file (MC only)")
    parser.add_argument("--isMC", action="store_true", help="Set this flag for MC")
    parser.add_argument("--inputDir", required=True, help="Directory with input ROOT files")
    parser.add_argument("--outputDir", default="output", help="Directory for final output files")
    parser.add_argument("--tempDir", default="/nfs_scratch/mithakor/temp",
                        help="Temporary directory for intermediate chunk outputs (only used if --finalHadd)")
    parser.add_argument("--finalHadd", action="store_true",
                        help="Merge all chunk outputs into one final file (enables hadd mode)")
    parser.add_argument("--finalOut", default="merged_output.root",
                        help="Final merged ROOT file name (only used if --finalHadd)")

    args = parser.parse_args()

    run_postproc(args.goldenjson, args.jetidjson, args.vetomapjson,
                 args.jercjson, args.jersmearjson, args.isMC,
                 args.inputDir, args.outputDir, args.tempDir,
                 finalOut=args.finalOut, finalHadd=args.finalHadd)