import os, sys
import ROOT
#from BTagSampleList import *
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
ROOT.gStyle.SetOptStat(0)
parser = argparse.ArgumentParser(description='')
parser.add_argument('--wp',help="working points of DeepJetFlavB",choices=["medium","tight","loose"],default="medium")
parser.add_argument('--year','-y',help="Enter the year of data taking",choices=["2016","2016APV","2017","2018"], required=True)
args = parser.parse_args()





#Select settings here. Filepaths to btag histograms are contained in BTagSampleList. MCsampleList contains path to MC samples (for cross section normalization).
tagger = 'DeepJetFlavB'   # Choose between CSVv2 and DeepCSV
wp = args.wp      # Choose between loose, medium, and tight
year = args.year       # Choose between 2016, 2017, and 2018
# name of root file that contains final btag efficiencies for selected channel
if year == "2016":
  outfilename = year + "/" + tagger+'_'+'2016_eff_'+wp.upper()+'.root'
  lumi = 16.81       
elif year == "2016APV":
  outfilename = year + "/" + tagger+'_'+'2016APV_eff_'+wp.upper()+'.root'
  lumi = 19.52
elif year == "2017":
  outfilename = year + "/" + tagger+'_'+'2017_eff_'+wp.upper()+'.root'
  lumi = 41.48
elif year == 2018:
  outfilename = year + "/" + tagger+'_'+'2018_eff_'+wp.upper()+'.root'
  lumi = 59.83

# Define helper functions for plotting
def makeTitle(tagger,wp,flavor,year):
  flavor = flavor.replace('_',' ')
  if ' b ' in flavor:
    flavor = 'b quark'
  elif ' c ' in flavor:
    flavor = 'c quark'
  else:
    flavor = 'light-flavor'
  string = "%s, %s %s WP for %s"%(flavor,tagger,wp,year)
  return string
  

def ensureTDirectory(file,dirname):
  dir = file.GetDirectory(dirname)
  if not dir:
    dir = file.mkdir(dirname)
    print ">>>   created directory %s in %s" % (dirname,file.GetName())
  dir.cd()
  return dir


# PREPARE numerator and denominator histograms per flavor
hists   = {}

for flavor in ['b','c','udsg']:
    histname = '%s_%s_%s'%(tagger,flavor,wp)
    hists[histname] = None        # numerator
    hists[histname+'_all'] = None # denominator

filename = "%s/btagHists_%s_%s.root"%(year,year,wp.upper())
file = ROOT.TFile(filename,'READ')

for histname in hists:
    hist = file.GetDirectory(year)
    hists[histname] = hist.Get(histname).Clone()
    hists[histname].SetDirectory(0)

file.Close()

# DIVIDE and SAVE histograms
print ">>>   writing to %s..."%(outfilename)
file = ROOT.TFile(outfilename,'RECREATE') 
ensureTDirectory(file,year)
for histname, hist in hists.iteritems():
    if 'all' in histname:
        continue
    histname_all = histname+'_all'
    histname_eff = 'eff_'+histname
    print ">>>     writing %s..."%(histname)
    print ">>>     writing %s..."%(histname_all)
    print ">>>     writing %s..."%(histname_eff)
    hist_all = hists[histname_all]
    hist = hist.Clone(histname_eff)
    drawlist=[histname,histname+'_all']
    for histo in drawlist:
      c = ROOT.TCanvas('c', 'c', 1400, 800)
      hists[histo].Draw('Colz')
      hists[histo].SetTitle('')
      hists[histo].GetZaxis().SetRangeUser(0, 700)
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
      if 'all' in histo:
        c.SaveAs(year+"/"+args.wp+"_figs"+"/Deno_"+tagger+"_"+histname+"_"+year+'_.pdf')
        c.SaveAs(year+"/"+args.wp+"_figs"+"/Deno_"+tagger+"_"+histname+"_"+year+'_.png')
      else:
        c.SaveAs(year+"/"+args.wp+"_figs"+"/Num_"+tagger+"_"+histname+"_"+year+'_.pdf')  
        c.SaveAs(year+"/"+args.wp+"_figs"+"/Num_"+tagger+"_"+histname+"_"+year+'_.png')  

        
      Pt = ROOT.TH1F()
      Eta = ROOT.TH1F()
      Pt =  hists[histo].ProjectionX()
      Pt.SetLineColor(9)
      Pt.GetXaxis().SetTitle("jet Pt (GeV)")
      Pt.GetYaxis().SetTitle("Entries")
      #Pt.GetYaxis().SetRangeUser(0,6000)
      Eta =  hists[histo].ProjectionY()
      Eta.SetLineColor(8)
      Eta.GetXaxis().SetTitle("jet Eta")
      Eta.GetYaxis().SetTitle("Entries")
      #Eta.GetYaxis().SetRangeUser(0,6000)
      c1 = ROOT.TCanvas('c1', 'c1 ', 1400, 800)
      c1.SetLogy()
      c1.SetGrid()
      Pt.Draw("HIST")
      if 'all' in histo:
        c1.SaveAs(year+"/oneDim/"+args.wp+"_figs"+"/PtDeno_"+tagger+"_"+histname+"_"+year+'_.pdf')  
        c1.SaveAs(year+"/oneDim/"+args.wp+"_figs"+"/PtDeno_"+tagger+"_"+histname+"_"+year+'_.png')  
      else:
        c1.SaveAs(year+"/oneDim/"+args.wp+"_figs"+"/PtNum_"+tagger+"_"+histname+"_"+year+'_.pdf')  
        c1.SaveAs(year+"/oneDim/"+args.wp+"_figs"+"/PtNum_"+tagger+"_"+histname+"_"+year+'_.png')  

      c2 = ROOT.TCanvas('c2', 'c2 ', 1400, 800)
      c2.SetGrid()
      Eta.Draw("HIST")
      if 'all' in histo:
        c2.SaveAs(year+"/oneDim/"+args.wp+"_figs"+"/EtaDeno_"+tagger+"_"+histname+"_"+year+'_.pdf')  
        c2.SaveAs(year+"/oneDim/"+args.wp+"_figs"+"/EtaDeno_"+tagger+"_"+histname+"_"+year+'_.png')  
      else:
        c2.SaveAs(year+"/oneDim/"+args.wp+"_figs"+"/EtaNum_"+tagger+"_"+histname+"_"+year+'_.pdf')
        c2.SaveAs(year+"/oneDim/"+args.wp+"_figs"+"/EtaNum_"+tagger+"_"+histname+"_"+year+'_.png')

    hist.SetTitle(makeTitle(tagger,wp,histname_eff,year))
    hist.Divide(hist_all)
    hist.Write(histname_eff,ROOT.TH2F.kOverwrite)
file.Close()
print ">>> "    


