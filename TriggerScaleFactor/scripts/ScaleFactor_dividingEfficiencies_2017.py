

import ROOT

dataRootFile_2017B = ROOT.TFile("2017/Run2017B_DataHistos.root")
dataDir_2017B = dataRootFile_2017B.GetDirectory("DataHistos")
dataDeno_2017B = dataDir_2017B.Get("h_data_denominator")
dataDeno_2017B.Sumw2()
dataNum_2017B = dataDir_2017B.Get("h_data_numerator")
dataNum_2017B.Sumw2()

dataRootFile_2017C = ROOT.TFile("2017/Run2017C_DataHistos.root")
dataDir_2017C = dataRootFile_2017C.GetDirectory("DataHistos")
dataDeno_2017C = dataDir_2017C.Get("h_data_denominator")
dataDeno_2017C.Sumw2()
dataNum_2017C = dataDir_2017C.Get("h_data_numerator")
dataNum_2017C.Sumw2()

dataRootFile_2017D = ROOT.TFile("2017/Run2017D_DataHistos.root")
dataDir_2017D = dataRootFile_2017D.GetDirectory("DataHistos")
dataDeno_2017D = dataDir_2017D.Get("h_data_denominator")
dataDeno_2017D.Sumw2()
dataNum_2017D = dataDir_2017D.Get("h_data_numerator")
dataNum_2017D.Sumw2()

dataRootFile_2017E = ROOT.TFile("2017/Run2017E_DataHistos.root")
dataDir_2017E = dataRootFile_2017E.GetDirectory("DataHistos")
dataDeno_2017E = dataDir_2017E.Get("h_data_denominator")
dataDeno_2017E.Sumw2()
dataNum_2017E = dataDir_2017E.Get("h_data_numerator")
dataNum_2017E.Sumw2()

dataRootFile_2017F = ROOT.TFile("2017/Run2017F_DataHistos.root")
dataDir_2017F = dataRootFile_2017F.GetDirectory("DataHistos")
dataDeno_2017F = dataDir_2017F.Get("h_data_denominator")
dataDeno_2017F.Sumw2()
dataNum_2017F = dataDir_2017F.Get("h_data_numerator")
dataNum_2017F.Sumw2()


dataNum_2017B.Add(dataNum_2017C)
dataNum_2017B.Add(dataNum_2017D)
dataNum_2017B.Add(dataNum_2017E)
dataNum_2017B.Add(dataNum_2017F)


dataDeno_2017B.Add(dataDeno_2017C)
dataDeno_2017B.Add(dataDeno_2017D)
dataDeno_2017B.Add(dataDeno_2017E)
dataDeno_2017B.Add(dataDeno_2017F)


dataEff = ROOT.TEfficiency(dataNum_2017B, dataDeno_2017B)
dataEff.SetStatisticOption(2)
dataDeno_2017B.Sumw2()
# dataDeno_2017B.Sumw2()
dataEff_hist = dataNum_2017B.Clone()
dataEff_hist.Sumw2()
dataDeno_2017B.Sumw2()
dataEff_hist.Divide(dataDeno_2017B)
dataEff_hist.Print("all")

mcRootFile = ROOT.TFile("2017/2017MC_MCHistos.root")
mcDir = mcRootFile.GetDirectory("MCHistos")

mcDeno = mcDir.Get("h_mc_denominator")
mcNum = mcDir.Get("h_mc_numerator")

# mcNum.GetXaxis().SetRangeUser(60,1000)
# mcDeno.GetXaxis().SetRangeUser(60,1000)

##############################################################
mcNum.SetLineColor(8)
mcDeno.SetLineColor(9)
mcNum.GetXaxis().SetTitle("MET [GeV]")
mcNum.GetYaxis().SetTitle("Entries")
mcDeno.GetXaxis().SetTitle("MET [GeV]")
mcDeno.GetYaxis().SetTitle("Entries")

dataNum_2017B.SetLineColor(8)
dataDeno_2017B.SetLineColor(9)
dataNum_2017B.GetXaxis().SetTitle("MET [GeV]")
dataNum_2017B.GetYaxis().SetTitle("Entries")
dataDeno_2017B.GetXaxis().SetTitle("MET [GeV]")
dataDeno_2017B.GetYaxis().SetTitle("Entries")
#############################################################


mcEff = ROOT.TEfficiency(mcNum, mcDeno)
mcEff.SetStatisticOption(2)
mcEff_hist = mcNum.Clone()
mcEff_hist.Sumw2()
mcDeno.Sumw2()
mcEff_hist.Divide(mcDeno)
mcEff_hist.Print("all")


MET_SF = dataEff_hist.Clone()
MET_SF.Sumw2()
mcEff_hist.Sumw2()
MET_SF.Divide(mcEff_hist)
MET_SF.SetMarkerColor(2)
MET_SF.SetMarkerStyle(8)
MET_SF.Print("all")
# for i in range(1,23):
# MET_SF.SetBinError(i,max(mcEff.GetEfficiencyErrorLow(i),mcEff.GetEfficiencyErrorUp(i),dataEff.GetEfficiencyErrorLow(i),dataEff.GetEfficiencyErrorUp(i)))

# Compute the errors in the overflow bin
# MET_SF.SetBinError(23,max(mcEff.GetEfficiencyErrorLow(23),mcEff.GetEfficiencyErrorUp(23),dataEff.GetEfficiencyErrorLow(23),dataEff.GetEfficiencyErrorUp(23)))


MetTrigSF = ROOT.TFile("2017/2017_MetTriggerSFs.root", "RECREATE")
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
can1.SaveAs("2017/Data_MC_Eff_2017.pdf")
can1.SaveAs("2017/Data_MC_Eff_2017.png")

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
can2.SaveAs("2017/SF_2017.pdf")
can2.SaveAs("2017/SF_2017.png")


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


plotPad.cd()
# plotPad.DrawFrame(100.,-0.5,500.,1.2)
# plotPad.Update()
# plotPad.SetFrameLineWidth(1)
# plotPad.SetTickx()
# plotPad.SetTicky()
# mcEff_hist.Draw("pE")
# dataEff_hist.Draw("same pE")
mcEff.SetTitle("")
mcEff.GetPaintedGraph().GetYaxis().SetRangeUser(-0.05, 1.1)
mcEff.GetPaintedGraph().GetXaxis().SetLimits(80.0, 500.0)
mcEff.GetPaintedGraph().GetYaxis().SetTickSize(0.01)
dataEff.GetPaintedGraph().GetYaxis().SetRangeUser(-0.05, 1.1)
dataEff.GetPaintedGraph().GetXaxis().SetLimits(80.0, 500.0)
dataEff.GetPaintedGraph().GetYaxis().SetTickSize(0.01)
mcEff.GetPaintedGraph().GetYaxis().SetTitle("Efficiency")
mcEff.GetPaintedGraph().GetYaxis().SetTitleSize(0.065)
mcEff.GetPaintedGraph().GetYaxis().SetLabelSize(0.05)
mcEff.GetPaintedGraph().GetYaxis().SetTitleOffset(0.60)
mcEff.GetPaintedGraph().GetXaxis().SetLabelSize(0.0)
mcEff.Draw("ap")
dataEff.Draw("same p")
legend = ROOT.TLegend(0.5100287, 0.12, 0.70, 0.22)
legend.SetFillStyle(1001)
legend.SetBorderSize(0)
legend.AddEntry(mcEff, "MC", "ep")
legend.AddEntry(dataEff, "Data", "ep")
legend.Draw("same")

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
lumiText = '41.48 fb^{-1}, 13 TeV'
cmsLatex.DrawLatex(0.85, 0.87, lumiText)


ratioPad.SetPad("pad2", "ratio", 0, 0, 1, 0.25)
ratioPad.SetFillColor(0)
ratioPad.SetTopMargin(0.02)  # 0.02
ratioPad.SetBottomMargin(0.35)
ratioPad.SetLeftMargin(0.15)
ratioPad.SetRightMargin(0.15)
# ratioPad.SetTickx(0)
# ratioPad.SetTicky(0)
ratioPad.SetFrameLineWidth(1)
ratioPad.SetGridy()


ratioPad.cd()
# ratioPad.DrawFrame(100.,0.4,500.,1.4)
# ratioPad.Update()
MET_SF.GetYaxis().SetNdivisions(5)
MET_SF.GetXaxis().SetTitle("MET[GeV]")
MET_SF.GetYaxis().SetTitle("#epsilon_{Data}/#epsilon_{MC}")
MET_SF.GetYaxis().SetRangeUser(0.5, 1.4)
MET_SF.GetYaxis().SetTitleOffset(0.30)
MET_SF.GetYaxis().SetTickSize(0.01)

MET_SF.GetYaxis().SetTitleSize(0.12)
MET_SF.GetYaxis().SetLabelSize(0.085)
MET_SF.GetXaxis().SetLabelSize(0.1)
MET_SF.GetXaxis().SetTitleSize(0.15)
MET_SF.GetXaxis().SetTitleOffset(0.70)
MET_SF.GetXaxis().SetRangeUser(80.0, 500.0)

MET_SF.SetTitle("")
MET_SF.SetMarkerStyle(8)
MET_SF.SetMarkerColor(2)
MET_SF.SetMarkerSize(0.5)
MET_SF.Draw("")

# MET_SF.SetTitleSize(0)


# backgroundStack.SetTitle(variableAxisTitleDictionary[variable])

# if(args.data):

# mcEff_hist.SetStats(0)
# mcEff_hist.SetTitle("")
# mcEff_hist.GetYaxis().SetRangeUser(0,1.2)
# mcEff_hist.GetXaxis().SetRangeUser(80,1000)

# mcEff_hist.GetYaxis().SetTitle("Efficiency")
# mcEff_hist.GetYaxis().SetTitleSize(0.065)
# mcEff_hist.GetYaxis().SetLabelSize(0.05)
# mcEff_hist.GetYaxis().SetTitleOffset(0.60)
# mcEff_hist.GetXaxis().SetLabelSize(0.0)


theCanvas.SaveAs("2017/Fancy_Eff_and_SF_2017.pdf")
theCanvas.SaveAs("2017/Fancy_Eff_and_SF_2017.png")


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
legend.SetBorderSize(0)
legend.AddEntry(mcNum, "NUM", "ep")
legend.AddEntry(mcDeno, "DEN", "ep")
legend.Draw("same")


canForhist1.SaveAs("2017/MonteCarlo_dist.pdf")
canForhist1.SaveAs("2017/MonteCarlo_dist.png")
############################################################
canForhist2 = ROOT.TCanvas("canForhist2", "canForhist2")
canForhist2.cd()
canForhist2.SetGrid()
canForhist2.SetLogy(1)
# dataDeno_2017B.SetLogy(1)
dataNum_2017B.SetStats(0)
dataNum_2017B.SetMaximum(
    dataDeno_2017B.GetMaximum() +
    0.4 *
    dataDeno_2017B.GetMaximum())
dataNum_2017B.SetMinimum(0.00001)
dataNum_2017B.Draw("HIST")
dataDeno_2017B.Draw("HIST SAME")
legend = ROOT.TLegend(0.5100287, 0.12, 0.70, 0.22)
legend.SetFillStyle(1001)
legend.SetBorderSize(0)
legend.AddEntry(dataNum_2017B, "NUM", "ep")
legend.AddEntry(dataDeno_2017B, "DEN", "ep")
legend.Draw("same")

canForhist2.SaveAs("2017/Data_dist.pdf")
canForhist2.SaveAs("2017/Data_dist.png")
