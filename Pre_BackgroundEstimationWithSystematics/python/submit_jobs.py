#!/usr/bin/env python3
import argparse
import glob
import os
import datetime
import sys

sys.path.append(os.path.join(os.environ["CMSSW_BASE"], "python"))

def normalize_path(p):
    if p.startswith("/hdfs/"):
        return "file:" + p  
    return p

def main(args):
    
    if args.inputList:
        with open(args.inputList) as f:
            inputFiles = [line.strip() for line in f if line.strip()]
        print(f"Using {len(inputFiles)} files from list {args.inputList}")
    elif args.inputDir:
        inputFiles = sorted(glob.glob(os.path.join(args.inputDir, "*.root")))
        print(f"Found {len(inputFiles)} input files in {args.inputDir}")
    else:
        print("You must provide either --inputDir or --inputList")
        sys.exit(1)

    if not inputFiles:
        print("No input files found!")
        sys.exit(1)

    timestamp = datetime.datetime.now().strftime('%d%b%y_%H%M')
    job_name = f"{args.jobName}_{timestamp}"
    overallSubmitDir = os.path.join(args.submitDirPath, job_name)
    dagLocation = os.path.join(overallSubmitDir, "dags")
    os.makedirs(os.path.join(dagLocation, "daginputs"), exist_ok=True)

    inputFileTextName = os.path.join(dagLocation, "daginputs", job_name + "_input.txt")
    with open(inputFileTextName, "w") as f:
        f.write("\n".join(normalize_path(p) for p in inputFiles))

    helper = "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/CommonAnalysisWithSystematics_6.py"
    extra_inputs = ",".join(filter(None, [
        "/afs/hep.wisc.edu/user/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/GoldenJSON_2024.json" if not args.isMC else "",
        "/afs/hep.wisc.edu/user/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/JetCorrBranchesAndMETSkimTrim/2024/jetid.json",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/jetvetomaps.json",
        helper,
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/Datadrop.txt",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/2024_DIB_samples.json",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/2024_DY_samples.json",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/2024_QCD_samples.json",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/2024_Radion_samples.json",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/2024_STop_samples.json",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/2024_TTbar_samples.json",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/2024_WJets_samples.json",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/2024_new_processed_samples.json",
        "/afs/hep.wisc.edu/home/mithakor/HH_bb_tautau_Analysis/Branch_addition_systematics/CMSSW_15_0_10/src/FinalStateHHbbtt/Pre_BackgroundEstimationWithSystematics/python/2024_old_samples.json.json"
    ]))

    commandList = [
        'farmoutAnalysisJobs',
        '--fwklite',
        '--infer-cmssw-path',
        '--input-files-per-job=1',
        '--job-generates-output-name',
        '--use-singularity=rhel9',
        f'--input-file-list={inputFileTextName}',
        '--assume-input-files-exist',
        '--max-usercode-size=350',
        '--use-hdfs',
        f'--submit-dir={overallSubmitDir}/submit',
        f'--output-dag-file={dagLocation}/dag',
        f'--output-dir={args.destination}/{job_name}',
        '--opsys=rhel9',
        '--memory-requirements=2000',
        '--disk-requirements=10000',
        '--input-dir=/',
        f'--extra-inputs={extra_inputs}',
        job_name,
        helper,   # must be executable with shebang
        '--',
        '\'--inputFile=$inputFileNames\'',
        '\'--outputFile=$outputFileName\'',
        f'\'--year={args.year}\'',
        f"'--cutflowDir={args.destination}/{job_name}/cutflows'"

    ]

    if args.runNominal:
        commandList.append("--runNominal")
    
    if args.isMC:
        commandList.append("--isMC")


    theCommand = " ".join([c for c in commandList if c.strip()])
    print(f"\nSubmitting job with command:\n{theCommand}\n")
    retcode = os.system(theCommand)
    print(f"farmoutAnalysisJobs exited with code {retcode}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Farmout submission for NanoAODTools postprocessing")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--inputDir", help="Directory with input ROOT files")
    group.add_argument("--inputList", help="Text file with one ROOT file per line")

    parser.add_argument("--destination", required=True, help="HDFS or local destination for output")
    parser.add_argument("--jobName", required=True, help="Base name for this submission")
    parser.add_argument("--submitDirPath", default="/nfs_scratch/"+os.environ["USER"]+"/JET_JES_JER_Softdrop_Condor_Jobs", help="Scratch area for submit files")
    parser.add_argument("--year", required=True, choices=["2024"])
    parser.add_argument("--isMC", action="store_true", help="Pass this flag to run on MC")
    parser.add_argument("--runNominal", action="store_true", help="Disable systematics")
    args = parser.parse_args()

    main(args)
