import ROOT
import numpy as np
import csv

# Define the mass points (in GeV)
mass_points = [
    1000,
    1200,
    1400,
    1600,
    1800,
    2000,
    2500,
    3000,
    3500,
    4000,
    4500]
year_list = ['2016', '2016APV', '2017', '2018']
# Full set of HLT paths
full_hlt_paths = {
    "2016": [
        "HLT_PFMET170_HBHECleaned",
        "HLT_PFMET170_HBHE_BeamHaloCleaned",
        "HLT_PFMET110_PFMHT110_IDTight",
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"],
    "2016APV": [
        "HLT_PFMET170_HBHECleaned",
        "HLT_PFMET170_HBHE_BeamHaloCleaned",
        "HLT_PFMET110_PFMHT110_IDTight",
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"],
    "2017": [
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"],
    "2018": [
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"]}
# Add more paths as needed


# Subset of HLT paths
subset_hlt_paths = {
    "2016": [
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"],
    "2016APV": [
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"],
    "2017": [
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"],
    "2018": [
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"]}

# Function to calculate signal efficiency for a set of HLT paths


def calculate_efficiency(input_file, hlt_paths, tree_name='Events'):
    # Open the ROOT file
    file = ROOT.TFile(input_file)
    tree = file.Get(tree_name)
    cutflow = file.Get('cutflow')

    # Total number of events in the tree
    total_events = cutflow.GetBinContent(1)

    # Initialize a counter for events passing any HLT path
    passing_events = 0

    # Loop over events
    for event in tree:
        if ((event.HTTvis_deltaR < 1.5) and (tree.HTTvis_m > 20)
                and (abs(tree.Hbb_met_phi) > 1)):
            for path in hlt_paths:
                # Check if the event passes the HLT path
                if getattr(
                        event, path):  # Assuming HLT decision is stored as a boolean
                    passing_events += 1
                    # Exit the loop when one path passes (assuming OR
                    # condition)
                    break

    # Calculate efficiency as ratio of passing events to total events
    efficiency = passing_events / total_events if total_events > 0 else 0
    file.Close()
    return efficiency

# Function to calculate percentage decrease with a - sign


def calculate_percentage_decrease(full_efficiency, subset_efficiency):
    if full_efficiency > 0:
        decrease = ((full_efficiency - subset_efficiency)
                    * 100) / (full_efficiency)
    else:
        decrease = 0
    return -decrease  # Return negative percentage to represent the decrease


basePath = '/hdfs/store/user/gparida/HHbbtt/Framework_Processed_Files/Full_Production_CMSSW_13_0_13_Nov24_23/CommonAnalysis_4scripts_ImprovementChecks/AK8Ordering_pt200/'
# Loop over each mass point
for year in year_list:
    # Output CSV file name

    print(("Entering the ERA = ", year))
    print("")
    print(("The full HLT Path set : ", full_hlt_paths[year]))
    print("")
    print(("The noMuOnly HLT Path set : ", subset_hlt_paths[year]))
    print("")
    output_csv = 'hlt_efficiency_comparison_{}.csv'.format(year)
    with open(output_csv, mode='w') as file:  # In Python 2, open in binary mode for writing
        writer = csv.writer(file)
        # Write the header row
        writer.writerow(['Mass Point (TeV)',
                         'Efficiency (Full Set)',
                         'Efficiency (NoMuSubset)',
                         'Percentage Decrease (%)'])
        for mass in mass_points:
            print(("......mass...", mass, "...file....",
                  'RadionTohhTohtatahbb_narrow_M-{}.root'.format(int(mass))))
            # Adjust filename format as needed
            input_file = basePath + year + \
                '/RadionTohhTohtatahbb_narrow_M-{}.root'.format(int(mass))
            # Calculate efficiency for the full set of HLT paths

            full_efficiency = calculate_efficiency(
                input_file, full_hlt_paths[year])

            # Calculate efficiency for the subset of HLT paths

            subset_efficiency = calculate_efficiency(
                input_file, subset_hlt_paths[year])
            # Calculate percentage decrease
            percentage_decrease = calculate_percentage_decrease(
                full_efficiency, subset_efficiency)
            # Convert mass to TeV for the CSV output
            mass_in_tev = mass / 1000.0
            # Write row to CSV
            writer.writerow([mass_in_tev, full_efficiency,
                            subset_efficiency, percentage_decrease])

    print(("Efficiencies and percentage decreases saved to for {} year in {}".format(
        year, output_csv)))
