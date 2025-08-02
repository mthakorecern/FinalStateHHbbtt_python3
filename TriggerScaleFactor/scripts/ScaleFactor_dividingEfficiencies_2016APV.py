

import ROOT

dataRootFile_2016B_ver1_HIPM = ROOT.TFile(
    "2016APV/Run2016B_ver1_HIPM_DataHistos.root", "READ")
dataDir__2016B_ver1_HIPM = dataRootFile_2016B_ver1_HIPM.GetDirectory(
    "DataHistos")
dataDeno_2016B_ver1_HIPM = dataDir__2016B_ver1_HIPM.Get("h_data_denominator")
dataDeno_2016B_ver1_HIPM.Sumw2()
dataNum_2016B_ver1_HIPM = dataDir__2016B_ver1_HIPM.Get("h_data_numerator")
dataNum_2016B_ver1_HIPM.Sumw2()

dataRootFile_2016B_ver2_HIPM = ROOT.TFile(
    "2016APV/Run2016B_ver2_HIPM_DataHistos.root", "READ")
dataDir_2016B_ver2_HIPM = dataRootFile_2016B_ver2_HIPM.GetDirectory(
    "DataHistos")
dataDeno_2016B_ver2_HIPM = dataDir_2016B_ver2_HIPM.Get("h_data_denominator")
dataDeno_2016B_ver2_HIPM.Sumw2()
dataNum_2016B_ver2_HIPM = dataDir_2016B_ver2_HIPM.Get("h_data_numerator")
dataNum_2016B_ver2_HIPM.Sumw2()

dataRootFile_2016C_HIPM = ROOT.TFile(
    "2016APV/Run2016C_HIPM_DataHistos.root", "READ")
dataDir_2016C_HIPM = dataRootFile_2016C_HIPM.GetDirectory("DataHistos")
dataDeno_2016C_HIPM = dataDir_2016C_HIPM.Get("h_data_denominator")
dataDeno_2016C_HIPM.Sumw2()
dataNum_2016C_HIPM = dataDir_2016C_HIPM.Get("h_data_numerator")
dataNum_2016C_HIPM.Sumw2()

dataRootFile_2016D_HIPM = ROOT.TFile(
    "2016APV/Run2016D_HIPM_DataHistos.root", "READ")
dataDir_2016D_HIPM = dataRootFile_2016D_HIPM.GetDirectory("DataHistos")
dataDeno_2016D_HIPM = dataDir_2016D_HIPM.Get("h_data_denominator")
dataDeno_2016D_HIPM.Sumw2()
dataNum_2016D_HIPM = dataDir_2016D_HIPM.Get("h_data_numerator")
dataNum_2016D_HIPM.Sumw2()

dataRootFile_2016E_HIPM = ROOT.TFile(
    "2016APV/Run2016E_HIPM_DataHistos.root", "READ")
dataDir_2016E_HIPM = dataRootFile_2016E_HIPM.GetDirectory("DataHistos")
dataDeno_2016E_HIPM = dataDir_2016E_HIPM.Get("h_data_denominator")
dataDeno_2016E_HIPM.Sumw2()
dataNum_2016E_HIPM = dataDir_2016E_HIPM.Get("h_data_numerator")
dataNum_2016E_HIPM.Sumw2()

dataRootFile_2016F_HIPM = ROOT.TFile(
    "2016APV/Run2016F_HIPM_DataHistos.root", "READ")
dataDir_2016F_HIPM = dataRootFile_2016F_HIPM.GetDirectory("DataHistos")
dataDeno_2016F_HIPM = dataDir_2016F_HIPM.Get("h_data_denominator")
dataDeno_2016F_HIPM.Sumw2()
dataNum_2016F_HIPM = dataDir_2016F_HIPM.Get("h_data_numerator")
dataNum_2016F_HIPM.Sumw2()

dataNum_2016B_ver1_HIPM.Add(dataNum_2016B_ver2_HIPM)
dataNum_2016B_ver1_HIPM.Add(dataNum_2016C_HIPM)
dataNum_2016B_ver1_HIPM.Add(dataNum_2016D_HIPM)
dataNum_2016B_ver1_HIPM.Add(dataNum_2016E_HIPM)
dataNum_2016B_ver1_HIPM.Add(dataNum_2016F_HIPM)

dataDeno_2016B_ver1_HIPM.Add(dataDeno_2016B_ver2_HIPM)
dataDeno_2016B_ver1_HIPM.Add(dataDeno_2016C_HIPM)
dataDeno_2016B_ver1_HIPM.Add(dataDeno_2016D_HIPM)
dataDeno_2016B_ver1_HIPM.Add(dataDeno_2016E_HIPM)
dataDeno_2016B_ver1_HIPM.Add(dataDeno_2016F_HIPM)

dataEff = ROOT.TEfficiency(dataNum_2016B_ver1_HIPM, dataDeno_2016B_ver1_HIPM)
dataEff.SetStatisticOption(2)
dataEff_hist = dataNum_2016B_ver1_HIPM.Clone()
dataEff_hist.Sumw2()
dataEff_hist.Divide(dataDeno_2016B_ver1_HIPM)
dataEff_hist.Print("all")

mcRootFile = ROOT.TFile("2016APV/2016APVMC_MCHistos.root", "READ")
mcDir = mcRootFile.GetDirectory("MCHistos")

print("Going to compute MC Eff")

mcDeno = mcDir.Get("h_mc_denominator")
print("Printing the denomiantor")
# mcDeno.Rebin(46)
mcDeno.Print("all")
mcNum = mcDir.Get("h_mc_numerator")
print("Printing the numerator")
mcNum.Print("all")
# mcNum.Rebin(46)

for i in range(1, 47):
    if (mcNum.GetBinContent(i) > mcDeno.GetBinContent(i)):
        print("Error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print((i, " Deno = ", mcDeno.GetBinContent(
            i), " Num = ", mcNum.GetBinContent(i)))

# print("Printing the numerator")
# mcNum.Print("all")
##############################################################
mcNum.SetLineColor(8)
mcDeno.SetLineColor(9)
mcNum.GetXaxis().SetTitle("MET [GeV]")
mcNum.GetYaxis().SetTitle("Entries")
mcDeno.GetXaxis().SetTitle("MET [GeV]")
mcDeno.GetYaxis().SetTitle("Entries")

dataNum_2016B_ver1_HIPM.SetLineColor(8)
dataDeno_2016B_ver1_HIPM.SetLineColor(9)
dataNum_2016B_ver1_HIPM.GetXaxis().SetTitle("MET [GeV]")
dataNum_2016B_ver1_HIPM.GetYaxis().SetTitle("Entries")
dataDeno_2016B_ver1_HIPM.GetXaxis().SetTitle("MET [GeV]")
dataDeno_2016B_ver1_HIPM.GetYaxis().SetTitle("Entries")
#############################################################


# mcEff = ROOT.TEfficiency(mcDeno,mcNum)
mcEff = ROOT.TEfficiency(mcNum, mcDeno)
mcEff.SetStatisticOption(2)
mcEff_hist = mcNum.Clone()
mcEff_hist.Divide(mcDeno)
mcEff_hist.Print("all")


MET_SF = dataEff_hist.Clone()
MET_SF.Divide(mcEff_hist)
MET_SF.SetMarkerColor(2)
MET_SF.SetMarkerStyle(8)
MET_SF.Print("all")


MetTrigSF = ROOT.TFile("2016APV/2016APV_MetTriggerSFs.root", "RECREATE")
METTrigSFdir = MetTrigSF.mkdir("SF")

METTrigSFdir.WriteObject(dataEff, "DataEfficiency")
METTrigSFdir.WriteObject(mcEff, "McEfficiency")
METTrigSFdir.WriteObject(MET_SF, "MET_SF")


# Plot the Data and MC efficiency on the same canvas
# MC
mcEff.SetLineColor(8)
mcEff.SetMarkerStyle(21)
mcEff.SetMarkerColor(8)
mcEff.SetMarkerSize(0.8)
mcEff.SetTitle("Efficiency Plot; MET (GeV) ; Efficiency")
# mcEff.SetTitle("Efficiency Plot")
# mcEff.GetXaxis().SetTitle("MET (GeV)")
# mcEff.GetYaxis().SetTitle("Efficiency")
# mcEff.SetMinimum(-0.05)
# mcEff.SetMaximum(1.1)

# Data
dataEff.SetLineColor(49)
dataEff.SetMarkerStyle(25)
dataEff.SetMarkerColor(49)
dataEff.SetMarkerSize(0.8)
dataEff.SetTitle("Efficiency Plot; MET (GeV) ; Efficiency")
# dataEff.SetMinimum(-0.05)
# dataEff.SetMaximum(1.1)
# dataEff.SetTitle("Efficiency Plot")
# dataEff.GetXaxis().SetTitle("MET (GeV)")
# dataEff.GetYaxis().SetTitle("Efficiency")
# dataEff.GetYaxis().SetRangeUser(-0.05,1.1)

can1 = ROOT.TCanvas("canvas1", "efficiency")
can1.SetGrid()
mcEff.Draw("ap")
dataEff.Draw("same p")
graph = mcEff.GetPaintedGraph()
# graph.SetMinimum(0.0)
# graph.SetMaximum(1.1)
legend = ROOT.TLegend(0.5100287, 0.12, 0.70, 0.22)
legend.SetFillStyle(1001)
legend.AddEntry(mcEff, "MC", "ep")
legend.AddEntry(dataEff, "Data", "ep")
legend.Draw("same")
can1.SaveAs("2016APV/Data_MC_Eff_2016APV.pdf")
can1.SaveAs("2016APV/Data_MC_Eff_2016APV.png")

# Plot the SF
MET_SF.SetLineColor(8)
MET_SF.SetStats(0)
MET_SF.SetMarkerStyle(21)
MET_SF.SetMarkerColor(8)
MET_SF.SetMarkerSize(1)
MET_SF.SetTitle("Scale factors [MET]")
MET_SF.GetXaxis().SetTitle("MET (GeV)")
MET_SF.GetYaxis().SetTitle("Data/MC")
MET_SF.GetYaxis().SetRangeUser(-0.05, 1.1)
can2 = ROOT.TCanvas("canvas2", "SF")
can2.SetGrid()
MET_SF.Draw("same")
can2.SaveAs("2016APV/SF_2016APV.pdf")
can2.SaveAs("2016APV/SF_2016APV.png")


theCanvas = ROOT.TCanvas("theCanvas", "theCanvas")
theCanvas.Divide(1, 2)

plotPad = ROOT.gPad.GetPrimitive('theCanvas_1')
ratioPad = ROOT.gPad.GetPrimitive('theCanvas_2')
plotPad.SetPad("pad1", "plot", 0, 0.25, 1, 1)
plotPad.SetFillColor(0)
plotPad.SetBorderMode(0)
plotPad.SetBorderSize(1)
# plotPad.SetTickx(1)
# plotPad.SetTicky(1)
plotPad.SetGrid()
plotPad.SetLeftMargin(0.15)  # 0.15
plotPad.SetRightMargin(0.15)  # 0.1
plotPad.SetTopMargin(0.14)  # 0.122 if the exponent is not present
plotPad.SetBottomMargin(0.025)
plotPad.SetFrameFillStyle(0)
plotPad.SetFrameLineStyle(0)
plotPad.SetFrameLineWidth(1)  # 1
plotPad.SetFrameBorderMode(0)
plotPad.SetFrameBorderSize(1)
# plotPad.SetErrorX()
# plotPad.SetLogy(1)
# plotPad.SetOptTitle(0)
#
ratioPad.SetPad("pad2", "ratio", 0, 0, 1, 0.25)
ratioPad.SetFillColor(0)
ratioPad.SetTopMargin(0.02)
ratioPad.SetBottomMargin(0.35)
ratioPad.SetLeftMargin(0.15)
ratioPad.SetRightMargin(0.15)
# ratioPad.SetTickx(0)
# ratioPad.SetTicky(0)
ratioPad.SetFrameLineWidth(1)
ratioPad.SetGridy()

ratioPad.cd()
MET_SF.Draw()
MET_SF.GetYaxis().SetNdivisions(5)
MET_SF.GetXaxis().SetTitle("MET[GeV]")
MET_SF.GetYaxis().SetTitle("#epsilon_{Data}/#epsilon_{MC}")
MET_SF.GetYaxis().SetRangeUser(0.5, 1.40)
MET_SF.GetYaxis().SetTitleOffset(0.30)

MET_SF.GetYaxis().SetTitleSize(0.12)
MET_SF.GetYaxis().SetLabelSize(0.085)
MET_SF.GetXaxis().SetLabelSize(0.1)
MET_SF.GetXaxis().SetTitleSize(0.15)
MET_SF.GetXaxis().SetTitleOffset(0.70)
MET_SF.GetXaxis().SetRangeUser(80, 500)

MET_SF.SetTitle("")
MET_SF.SetMarkerStyle(8)
MET_SF.SetMarkerColor(2)
MET_SF.SetMarkerSize(0.5)
MET_SF.GetYaxis().SetTickSize(0.01)
# MET_SF.SetTitleSize(0)


plotPad.cd()
plotPad.SetFrameLineWidth(1)
# plotPad.SetTickx()
# plotPad.SetTicky()
mcEff.Draw("ap")
dataEff.Draw("same p")
legend = ROOT.TLegend(0.5100287, 0.12, 0.70, 0.22)
legend.SetFillStyle(1001)
legend.SetBorderSize(0)
legend.AddEntry(mcEff, "MC", "ep")
legend.AddEntry(dataEff, "Data", "ep")
legend.Draw("same")
# backgroundStack.SetTitle(variableAxisTitleDictionary[variable])
mcEff.SetTitle("")
mcEff.GetPaintedGraph().GetYaxis().SetRangeUser(-0.05, 1.1)
mcEff.GetPaintedGraph().GetYaxis().SetLimits(-0.05, 1.1)
mcEff.GetPaintedGraph().GetXaxis().SetLimits(80, 500)
mcEff.GetPaintedGraph().GetYaxis().SetTickSize(0.01)
dataEff.GetPaintedGraph().GetYaxis().SetRangeUser(-0.05, 1.1)
dataEff.GetPaintedGraph().GetXaxis().SetLimits(80.0, 500.0)
dataEff.GetPaintedGraph().GetYaxis().SetTickSize(0.01)
# if(args.data):
# DatasetObjects["Data"].HistogramName.GetXaxis().SetErrorX()
mcEff.GetPaintedGraph().GetYaxis().SetTitle("Efficiency")
mcEff.GetPaintedGraph().GetYaxis().SetTitleSize(0.065)
mcEff.GetPaintedGraph().GetYaxis().SetLabelSize(0.05)
# backgroundStack.GetYaxis().SetTitleOffset(1.58)
mcEff.GetPaintedGraph().GetYaxis().SetTitleOffset(0.60)
# backgroundStack.GetYaxis().SetTitleSize(1)
mcEff.GetPaintedGraph().GetXaxis().SetLabelSize(0.0)

cmsLatex = ROOT.TLatex()
cmsLatex.SetTextSize(0.06)
cmsLatex.SetNDC(True)
cmsLatex.SetTextFont(61)
cmsLatex.SetTextAlign(11)
# cmsLatex.DrawLatex(0.1,0.92,"CMS")
cmsLatex.DrawLatex(0.15, 0.87, "CMS")
cmsLatex.SetTextFont(52)
# cmsLatex.DrawLatex(0.1+0.08,0.92,"Preliminary")
cmsLatex.DrawLatex(0.15 + 0.08, 0.87, "Preliminary")
cmsLatex.SetTextAlign(31)
cmsLatex.SetTextFont(42)
lumiText = '19.50 fb^{-1}, 13 TeV'
cmsLatex.DrawLatex(0.85, 0.87, lumiText)

theCanvas.SaveAs("2016APV/Fancy_Eff_and_SF_2016APV.pdf")
theCanvas.SaveAs("2016APV/Fancy_Eff_and_SF_2016APV.png")

# Plot the Numerator and Denominator histograms on the same plot,
# preferably using a log scale
canForhist1 = ROOT.TCanvas("canForhist1", "canForhist1")
canForhist1.cd()
canForhist1.SetLogy()
canForhist1.SetGrid()
# mcNum.SetLogy(1)
# mcDeno.SetLogy(1)
mcNum.SetMaximum(mcDeno.GetMaximum() + 0.4 * mcDeno.GetMaximum())
mcNum.SetMinimum(0.00001)
mcNum.SetStats(0)
mcNum.Draw("HIST")
mcDeno.Draw("HIST SAME")
legend = ROOT.TLegend(0.5100287, 0.12, 0.70, 0.22)
legend.SetFillStyle(1001)
legend.AddEntry(mcNum, "NUM", "ep")
legend.AddEntry(mcDeno, "DEN", "ep")
legend.Draw("same")


canForhist1.SaveAs("2016APV/MonteCarlo_dist.pdf")
canForhist1.SaveAs("2016APV/MonteCarlo_dist.png")
############################################################
canForhist2 = ROOT.TCanvas("canForhist2", "canForhist2")
canForhist2.cd()
canForhist2.SetGrid()
canForhist2.SetLogy(1)
# dataDeno_2016B_ver1_HIPM.SetLogy(1)
dataNum_2016B_ver1_HIPM.SetStats(0)
dataNum_2016B_ver1_HIPM.SetMaximum(
    dataDeno_2016B_ver1_HIPM.GetMaximum() +
    0.4 *
    dataDeno_2016B_ver1_HIPM.GetMaximum())
dataNum_2016B_ver1_HIPM.SetMinimum(0.00001)
dataNum_2016B_ver1_HIPM.Draw("HIST")
dataDeno_2016B_ver1_HIPM.Draw("HIST SAME")
legend = ROOT.TLegend(0.5100287, 0.12, 0.70, 0.22)
legend.SetFillStyle(1001)
legend.AddEntry(dataNum_2016B_ver1_HIPM, "NUM", "ep")
legend.AddEntry(dataDeno_2016B_ver1_HIPM, "DEN", "ep")
legend.Draw("same")

canForhist2.SaveAs("2016APV/Data_dist.pdf")
canForhist2.SaveAs("2016APV/Data_dist.png")
