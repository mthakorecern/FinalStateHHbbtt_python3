import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gStyle.SetPaintTextFormat("4.2f")


import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--wp',help="working points of DeepJetFlavB",choices=["medium","tight","loose"],default="medium")
parser.add_argument('--year','-y',help="Enter the year of data taking",choices=["2016","2016APV","2017","2018"], required=True)
args = parser.parse_args()



#Select settings here.
tagger = 'DeepJetFlavB'   # Choose between CSVv2 and DeepCSV
wp = args.wp      # Choose between loose, medium, and tight
year = args.year       # Choose between 2016, 2017, and 2018


#Remove stats box from histograms
ROOT.gStyle.SetOptStat(0)

# name of root file that contains final btag efficiencies for selected channel

filename = year + "/" + tagger+'_'+year+'_eff_'+wp.upper()+'.root'

outputfilename = year + "/" +wp+"_figs"+ "/" + tagger+'_'+year+'_eff_'+wp+'.root'

for flavor in ['b','c','udsg']:
    c = ROOT.TCanvas('c', 'c', 1400, 800)
    histname_eff = 'eff_%s_%s_%s'%(tagger,flavor,wp)
    file = ROOT.TFile(filename,'READ')
    dir = file.GetDirectory(year)
    hist = dir.Get(histname_eff)
    if flavor == "udsg":
      hist.GetZaxis().SetRangeUser(0, 0.50)
    else:  
      hist.GetZaxis().SetRangeUser(0, 1.0)
    hist.Draw('colz text45 E')

    #histogram text/settings
    hist.SetTitle('')
    title = ROOT.TLatex()
    title.SetTextSize(0.045)
    title.DrawLatexNDC(.12, .91, 'CMS')
    title.SetTextSize(0.03)
    title.DrawLatexNDC(.18, .91, '#bf{#it{Preliminary}}')
    title.SetTextSize(0.03)
    title_x = .76
    title_y = .91
    if year == "2016":
      title.DrawLatexNDC(title_x, title_y, '#bf{16.81 fb^{-1} (13 TeV)}')
    if year == "2016APV":
      title.DrawLatexNDC(title_x, title_y, '#bf{19.52 fb^{-1} (13 TeV)}')
    elif year == "2017":
      title.DrawLatexNDC(title_x, title_y, '#bf{41.5 fb^{-1} (13 TeV)}')
    elif year == "2018":
      title.DrawLatexNDC(title_x, title_y, '#bf{59.7 fb^{-1} (13 TeV)}')

    saveDir = year+'/'
    #if not os.path.exists(saveDir): os.makedirs(saveDir)
    c.SaveAs(outputfilename.replace('.root',histname_eff)+'.pdf')
    c.SaveAs(outputfilename.replace('.root',histname_eff)+'.png')