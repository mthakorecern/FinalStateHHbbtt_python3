#!/usr/bin/env python3
import argparse
import glob
import os
import datetime
import sys

# Ensure local python modules are visible on worker nodes
sys.path.append(os.path.join(os.environ["CMSSW_BASE"], "python"))

def normalize_path(p):
    if p.startswith("/hdfs/store/"):
        return p.replace("/hdfs", "", 1)  # drop the leading /hdfs
    return p

def main(args):
    # 1. Collect input files
    if args.inputList:
        with open(args.inputList) as f:
            inputFiles = [line.strip() for line in f if line.strip()]
        print(f"[INFO] Using {len(inputFiles)} files from list {args.inputList}")
    elif args.inputDir:
        inputFiles = sorted(glob.glob(os.path.join(args.inputDir, "*.root")))
        print(f"[INFO] Found {len(inputFiles)} input files in {args.inputDir}")
    else:
        print("[ERROR] You must provide either --inputDir or --inputList")
        sys.exit(1)

    if not inputFiles:
        print("[ERROR] No input files found!")
        sys.exit(1)

    # 2. Unique job name
    timestamp = datetime.datetime.now().strftime('%d%b%y_%H%M')
    job_name = f"{args.jobName}_{timestamp}"
    overallSubmitDir = os.path.join(args.submitDirPath, job_name)
    dagLocation = os.path.join(overallSubmitDir, "dags")
    os.makedirs(os.path.join(dagLocation, "daginputs"), exist_ok=True)

    # 3. Write file list for farmout
    inputFileTextName = os.path.join(dagLocation, "daginputs", job_name + "_input.txt")
    with open(inputFileTextName, "w") as f:
        f.write("\n".join(normalize_path(p) for p in inputFiles))

    # 4. Build farmout command
    helper = "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/FarmoutAnalysisJob_attempt/nanoPostProc_helper.py"

    extra_inputs = ",".join(filter(None, [
        args.goldenjson if not args.isMC else "",
        args.jetidjson,
        args.vetomapjson,
        args.jercjson,
        helper,
    ]))

    commandList = [
        'farmoutAnalysisJobs',
        '--fwklite',
        '--infer-cmssw-path',
        '--input-files-per-job=1',
        '--use-singularity=rhel9',
        f'--input-file-list={inputFileTextName}',
        '--assume-input-files-exist',
        '--max-usercode-size=350',
        '--use-hdfs',
        f'--submit-dir={overallSubmitDir}/submit',
        f'--output-dag-file={dagLocation}/dag',
        f'--output-dir={args.destination}/{job_name}',
        '--opsys=rhel9',
        '--memory-requirements=5000',
        '--disk-requirements=10000',
        '--input-dir=/',
        f'--extra-inputs={extra_inputs}',
        job_name,
        helper,   # must be executable with shebang
        '--',
        '\'--inputFile=$inputFileNames\'',
        '\'--outputFile=$outputFileName\'',
        ] + ([f'--goldenjson={os.path.basename(args.goldenjson)}'] if (args.goldenjson and not args.isMC) else []) \
        + [
            f'--jetidjson={os.path.basename(args.jetidjson)}',
            f'--vetomapjson={os.path.basename(args.vetomapjson)}',
            f'--jercjson={os.path.basename(args.jercjson)}',
        ] + (['--isMC'] if args.isMC else [])

    theCommand = " ".join([c for c in commandList if c.strip()])
    print(f"\n[INFO] Submitting job with command:\n{theCommand}\n")
    retcode = os.system(theCommand)
    print(f"[DEBUG] farmoutAnalysisJobs exited with code {retcode}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Farmout submission for NanoAODTools postprocessing")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--inputDir", help="Directory with input ROOT files")
    group.add_argument("--inputList", help="Text file with one ROOT file per line")

    parser.add_argument("--destination", required=True, help="HDFS or local destination for output")
    parser.add_argument("--jobName", required=True, help="Base name for this submission")
    parser.add_argument("--submitDirPath", default="/nfs_scratch/"+os.environ["USER"]+"/JetFatPuppiMET_JES_JER_VetoMaps_Condor",
                        help="Scratch area for submit files (default =/nfs_scratch/mithakor/JetFatPuppiMET_JES_JER_VetoMaps_Condor)")
    parser.add_argument("--goldenjson", required=False)
    parser.add_argument("--jetidjson", required=True)
    parser.add_argument("--vetomapjson", required=True)
    parser.add_argument("--jercjson", required=True)
    parser.add_argument("--isMC", action="store_true", help="Pass this flag to run on MC")
    args = parser.parse_args()
    main(args)
