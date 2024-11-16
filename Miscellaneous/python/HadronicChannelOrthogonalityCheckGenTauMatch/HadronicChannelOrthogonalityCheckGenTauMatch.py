import ROOT

# Base path where the NanoAOD files are located
base_path = "/hdfs/store/user/gparida/HHbbtt/Framework_Processed_Files/Full_Production_CMSSW_13_0_13_Nov24_23/CommonAnalysis_4scripts_btagProcessorWeights_Aug26_24/2016/signalregion/TauToAK8Matching_Files/AdditionalLeptonVeto"

# Define the mass points and signals
mass_points = [1000, 3000]  # in TeV
signals = ["HHbtautau", "HHbbWW1L1Nu2Q", "HHbbZZ2L2Nu"]

# Loop over each mass point and signal
txtfile = open("lepveto-miniiso/output_percentages.txt", "w")
for mass in mass_points:
    print("Processing mass...",mass)
    txtfile.write("\n\n\n\n\nProcessing mass...{}".format(mass))
    for signal in signals:
        print("Processing signal...",signal)
        # Construct the full path to the file
        if signal == "HHbtautau":
            file_path = "{}/RadionTohhTohtatahbb_narrow_M-{}.root".format(base_path,str(mass))
        elif signal == "HHbbWW1L1Nu2Q":
            file_path = "{}/RadionHHbbwww_Lnu2Q_narrow_M-{}.root".format(base_path,str(mass))
        elif signal == "HHbbZZ2L2Nu":
            file_path = "{}/RadionHHbbwww_2L2nu_narrow_M-{}.root".format(base_path,str(mass))

        total_events_SR = 0
        total_events_SR_tt = 0
        total_events_SR_et = 0
        total_events_SR_mt = 0

        total_events_SR_tt_rejected = 0
        total_events_SR_et_rejected = 0
        total_events_SR_mt_rejected = 0

        # Open the NanoAOD root file
        file = ROOT.TFile(file_path)
        tree = file.Get("Events")  # Assuming the tree is called 'Events'

        # Create a 2D histogram for genpartFlav of leading and subleading hadronic taus
        n_bins = 8
        hist2d = ROOT.TH2F("hist2d", "2D Normalized (%) {} ({} TeV);;".format(signal,str(mass/1000)), 
                           n_bins, -1, 7, n_bins, -1, 7)

        # Loop over events and fill the histogram for the fully hadronic channel (channel == 0)
        for event in tree:
            if ((event.HTTvis_deltaR>0) and (event.HTTvis_deltaR<1.5) and (event.ngood_MediumJets==0) and (event.Hbb_Loose==1) and (event.X_m>=750) and (event.X_m<=5010) and (abs(event.Hbb_met_phi)>1) and (event.HTTvis_m>20)): 
                total_events_SR += 1
                if (event.channel==0):
                    total_events_SR_tt += 1
                elif (event.channel==1):
                    total_events_SR_et += 1
                    if (event.miniiso_addlepton_vetoflag_all==1):
                        total_events_SR_et_rejected+=1
                elif (event.channel==2):
                    total_events_SR_mt += 1
                    if (event.miniiso_addlepton_vetoflag_all==1):
                        total_events_SR_mt_rejected+=1                                         
            if ((event.channel == 0) and (event.HTTvis_deltaR>0) and (event.HTTvis_deltaR<1.5) and (event.ngood_MediumJets==0) and (event.Hbb_Loose==1) and (event.X_m>=750) and (event.X_m<=5010) and (abs(event.Hbb_met_phi)>1) and (event.HTTvis_m>20)):  # Select fully hadronic channel
                # Extract leading and subleading tau genpartFlav values from AK8toTau_gentruth
                #if len(event.AK8toTau_gentruth) >= 2:
                if (event.miniiso_addlepton_vetoflag_all==0):
                    leading_tau_flav = event.AK8toTau_gentruth[0]
                    subleading_tau_flav = event.AK8toTau_gentruth[1]
                
                    # Fill the 2D histogram with the genpartFlav of leading and subleading taus
                    hist2d.Fill(leading_tau_flav, subleading_tau_flav)
                elif (event.miniiso_addlepton_vetoflag_all==1):
                    total_events_SR_tt_rejected += 1

                    

        # Normalize the histogram to represent percentages
        hist2d.Scale(100.0 / hist2d.Integral())

        print ("Total events in SR for {} is ".format(signal), total_events_SR)
        #print ("Fully Hadronic events in SR for {} is ".format(signal), total_events_SR_tt,"....percentage of total is... ",(float(total_events_SR_tt*100)/total_events_SR))
        #print ("Tau-Electron events in SR for {} is ".format(signal), total_events_SR_et,"....percentage of total is... ",(float(total_events_SR_et*100)/total_events_SR))
        #print ("Tau-Muon events in SR for {} is ".format(signal), total_events_SR_mt,"....percentage of total is... ",(float(total_events_SR_mt*100)/total_events_SR))
        print ("Fully Hadronic vetoed events in SR for {} is ".format(signal), total_events_SR_tt_rejected,"....percentage of total is... ",(float(total_events_SR_tt_rejected*100)/total_events_SR_tt))
        print ("Tau-Electron vetoed events in SR for {} is ".format(signal), total_events_SR_et_rejected,"....percentage of total is... ",(float(total_events_SR_et_rejected*100)/total_events_SR_et))
        print ("Tau-Muon vetoed events in SR for {} is ".format(signal), total_events_SR_mt_rejected,"....percentage of total is... ",(float(total_events_SR_mt_rejected*100)/total_events_SR_mt))

        #with open("lepveto-stdiso/output_percentages.txt", "w") as txtfile:
        # Write total events in SR
        txtfile.write("\n\n\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n\nTotal events in SR for {} is {}\n".format(signal, total_events_SR))
        # Uncomment the following lines to include fully hadronic and tau-lepton categories if needed
        # file.write("Fully Hadronic events in SR for {} is {} ....percentage of total is... {:.2f}%\n".format(signal, total_events_SR_tt, (float(total_events_SR_tt * 100) / total_events_SR)))
        # file.write("Tau-Electron events in SR for {} is {} ....percentage of total is... {:.2f}%\n".format(signal, total_events_SR_et, (float(total_events_SR_et * 100) / total_events_SR)))
        # file.write("Tau-Muon events in SR for {} is {} ....percentage of total is... {:.2f}%\n".format(signal, total_events_SR_mt, (float(total_events_SR_mt * 100) / total_events_SR)))
        # Write vetoed events in SR
        txtfile.write("Fully Hadronic vetoed events in SR for {} is {} ....percentage of total is... {:.2f}%\n".format(signal, total_events_SR_tt_rejected, (float(total_events_SR_tt_rejected * 100) / total_events_SR_tt)))
        txtfile.write("Tau-Electron vetoed events in SR for {} is {} ....percentage of total is... {:.2f}%\n".format(signal, total_events_SR_et_rejected, (float(total_events_SR_et_rejected * 100) / total_events_SR_et)))
        txtfile.write("Tau-Muon vetoed events in SR for {} is {} ....percentage of total is... {:.2f}%\n".format(signal, total_events_SR_mt_rejected, (float(total_events_SR_mt_rejected * 100) / total_events_SR_mt)))

        # Create a canvas to draw the histogram
        canvas = ROOT.TCanvas("canvas", "2D Normalized Tau GenPartFlav - {} TeV {}".format(str(mass),signal), 800, 600)
        ROOT.gStyle.SetOptStat(0)

        canvas.SetLeftMargin(0.15)  # Increase the left margin
        canvas.SetRightMargin(0.1)  # Adjust the right margin if necessary
        canvas.SetBottomMargin(0.10)  # Adjust the bottom margin if x-axis labels need space
        canvas.SetTopMargin(0.10)  # Adjust the top margin for the title

        # Set a single color palette (e.g., shades of blue)
        ROOT.gStyle.SetPalette(30)  # You can switch to ROOT.kGrayScale for grayscale
        ROOT.gStyle.SetPaintTextFormat("4.1f")  # Format to show the percentage with one decimal place

        # Draw the 2D histogram
        hist2d.Draw("COLZ TEXT")  # "TEXT" adds the numbers in each box

        hist2d.GetXaxis().SetTitle("Leading #tau")
        hist2d.GetXaxis().SetLabelSize(0.05)
        hist2d.GetYaxis().SetTitle("SubLeading #tau")
        hist2d.GetYaxis().SetTitleOffset(1.9)
        hist2d.GetYaxis().SetLabelSize(0.05)
        # Set axis labels for the bins (1-6 for genpartFlav)
        labels = ["","jets #rightarrow #tau_{h}", "e #rightarrow #tau_{h}", "#mu #rightarrow #tau_{h}", "#tau_{e} #rightarrow #tau_{h}", "#tau_{#mu} #rightarrow #tau_{h}", "#tau_{h} #rightarrow #tau_{h}",""]
        for i in range(1,8):
            #print(i)
            hist2d.GetXaxis().SetBinLabel(i, labels[i-1])
            hist2d.GetYaxis().SetBinLabel(i, labels[i-1])

        # Draw a color palette to show the normalized percentage
        hist2d.Draw("COLZ TEXT")  # Redraw with the text option to show percentages

        # Update and save the canvas with a filename reflecting the mass point and signal
        output_filename = "lepveto-miniiso/tau_genpartFlav_2D_{}TeV_{}_lepveto.png".format(str(mass),signal)
        canvas.Update()
        canvas.SaveAs(output_filename)

        # Show the canvas (optional, remove if running in batch mode)
        canvas.Draw()

        # Close the file after processing
        file.Close()

txtfile.close()
print("All plots generated and saved.")
