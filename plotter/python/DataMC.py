from numbers import Integral
import ROOT
import argparse  # Importing root and package to take arguments
# from variables import *
import os
from array import array
import glob
from re import search
import ctypes
import math
# if(args.topRW):
from variables_topCR import *
# else:
# from variables import *


class MakeHistograms(object):
    # constructor to initialize the objects
    def __init__(self, RootFilePath, RootFileName, userWeight="1.0"):
        self.RootFileName = ROOT.TFile(RootFilePath + RootFileName + '.root')
        self.HistogramName = None
        self.PassFailHistogramName = ROOT.TH1F(
            "PassFailHist", "PassFailHist", 1, 0, 1)
        self.userWeight = userWeight

    # Cut creating member function

    def CreateCutString(self, standardCutString,
                        otherCuts,
                        weighting):
        # cutString = weighting+'*('+standardCutString+' && '
        if standardCutString is not None:
            cutString = weighting + \
                '*(' + '(' + standardCutString + ')' + ' && '
            if otherCuts is not None:
                for cut in otherCuts:
                    cutString += '(' + cut + ')' + ' && '
        else:
            cutString = weighting + ' && '
        # removing the && at the very end of the final cutstring
        cutString = cutString[:len(cutString) - 3]
        cutString += ')'
        return cutString

    # Histogram Making member function and storing it in an attribute
    def StandardDraw(self, theFile,
                     variable,
                     standardCutString,
                     additionalSelections,
                     histogramName,
                     theWeight='FinalWeighting'):

        theTree = theFile.Get('Events')

        # print ('g'+variable+'>>'+histogramName+'('+variableSettingDictionary[variable]+')',
        #                 self.CreateCutString(standardCutString,
        #                                 additionalSelections,theWeight))
        ##
        # theTree.Draw('g'+variable+'>>'+histogramName+'('+variableSettingDictionary[variable]+')',
        #        self.CreateCutString(standardCutString,
        #                        additionalSelections,theWeight))
        # print ("uhoh No g in it")
        # if (variable=="X_m"):
        # #bin_edges = [750, 900, 1050, 1200, 1350, 1500, 1650, 1800, 1950, 2200, 2450, 5500] # Main Region
        # bin_edges = [750,1050,1350,1650,1950,5500] # validation Region
        # n_bins = len(bin_edges) - 1
        # histogram = ROOT.TH1F(histogramName, "Title", n_bins, array('d', bin_edges))
        # theTree.Draw(variable + ">>" + histogramName, self.CreateCutString(standardCutString, additionalSelections, theWeight))
        #
        # else:
        print(
            (
                variable +
                '>>' +
                histogramName +
                '(' +
                variableSettingDictionary[variable] +
                ')',
                self.CreateCutString(
                    standardCutString,
                    additionalSelections,
                    theWeight)))

        theTree.Draw(
            variable +
            '>>' +
            histogramName +
            '(' +
            variableSettingDictionary[variable] +
            ')',
            self.CreateCutString(
                standardCutString,
                additionalSelections,
                theWeight))

        # print ('g'+variable+'>>'+histogramName+'('+variableSettingDictionary[variable]+')',
        #                     self.CreateCutString(standardCutString,
        #                                     additionalSelections,theWeight))

        # 3
    # so, if the tree has no entries, root doesn't even hand back an empty histogram
    # and therefore this ends up trying to get clone a none type
    # pass the None forward, and we can let the Add handle this
        try:
            theHisto = ROOT.gDirectory.Get(histogramName).Clone()
        except ReferenceError:
            theHisto = None
        # return theHisto
        self.HistogramName = theHisto


def clubHistograms(list, histObjects):
    clubHist = None
    for name in list:
        # print("Name = ",name,"...Raw Entries = ",histObjects[name].HistogramName.GetEntries(),"....Integral = ",histObjects[name].HistogramName.Integral())
        if histObjects[name].HistogramName is not None:
            if clubHist is None:
                clubHist = histObjects[name].HistogramName.Clone()
                continue
            clubHist.Add(histObjects[name].HistogramName)
    return clubHist


def MakeStackErrors(theStack):
    denominatorHistos = theStack.GetHists().At(0).Clone()
    denominatorHistos.Reset()

    for i in range(0, theStack.GetNhists()):
        denominatorHistos.Add(theStack.GetHists().At(i))

    theErrorHisto = denominatorHistos.Clone()
    theErrorHisto.Reset()

    for i in range(0, denominatorHistos.GetNbinsX() + 1):
        theErrorHisto.SetBinContent(i, denominatorHistos.GetBinContent(i))
        theErrorHisto.SetBinError(i, denominatorHistos.GetBinError(i))
    theErrorHisto.SetLineColor(0)
    theErrorHisto.SetLineWidth(0)
    theErrorHisto.SetMarkerStyle(0)
    theErrorHisto.SetFillStyle(3001)
    theErrorHisto.SetFillColor(15)
    return theErrorHisto

    # make the statistical errors on the prediction stack


def MakeStackErrors(theStack):
    denominatorHistos = theStack.GetHists().At(0).Clone()
    denominatorHistos.Reset()

    for i in range(0, theStack.GetNhists()):
        denominatorHistos.Add(theStack.GetHists().At(i))

    theErrorHisto = denominatorHistos.Clone()
    theErrorHisto.Reset()

    for i in range(0, denominatorHistos.GetNbinsX() + 1):
        theErrorHisto.SetBinContent(i, denominatorHistos.GetBinContent(i))
        theErrorHisto.SetBinError(i, denominatorHistos.GetBinError(i))
    theErrorHisto.SetLineColor(0)
    theErrorHisto.SetLineWidth(0)
    theErrorHisto.SetMarkerStyle(0)
    theErrorHisto.SetFillStyle(3008)
    # theErrorHisto.SetFillColorAlpha(ROOT.kViolet+1)
    # theErrorHisto.SetFillStyle(1001)
    # theErrorHisto.SetFillColor(ROOT.TColor.GetColor("#FFE599"))
    theErrorHisto.SetFillColor(ROOT.TColor.GetColor("#545252"))
    # theErrorHisto.SetLineColor(ROOT.kMagenta+2)
    return theErrorHisto


# make the ratio histograms and associated errors
def MakeRatioHistograms(dataHisto, backgroundStack, variable):
    ratioHist = dataHisto.Clone()

    denominatorHistos = dataHisto.Clone()
    denominatorHistos.Reset()
    for i in range(0, backgroundStack.GetNhists()):
        denominatorHistos.Add(backgroundStack.GetHists().At(i))
    ratioHist.Divide(denominatorHistos)
    finalRatioHist = ratioHist.Clone()
    for i in range(1, finalRatioHist.GetNbinsX() + 1):
        try:
            finalRatioHist.SetBinError(
                i,
                (dataHisto.GetBinError(i) /
                 dataHisto.GetBinContent(i)) *
                ratioHist.GetBinContent(i))
        except ZeroDivisionError:
            finalRatioHist.SetBinError(i, 0)

    finalRatioHist.SetMarkerStyle(20)
    finalRatioHist.SetTitle("")
    finalRatioHist.GetYaxis().SetTitle("Data/Prediction")
    # finalRatioHist.GetYaxis().SetTitleSize(0.1)
    # finalRatioHist.GetYaxis().SetTitleSize(0.1)
    finalRatioHist.GetYaxis().SetTitleSize(0.10)
    finalRatioHist.GetYaxis().SetTitleOffset(0.32)
    finalRatioHist.GetYaxis().CenterTitle()
    finalRatioHist.GetYaxis().SetLabelSize(0.1)
    # finalRatioHist.GetYaxis().SetNdivisions(6,0,0)
    finalRatioHist.GetYaxis().SetNdivisions(5)
    # finalRatioHist.GetYaxis().SetRangeUser(1.3*1.05,0.7*0.95) #this doesn't
    # seem to take effect here?
    finalRatioHist.GetXaxis().SetTitleOffset(0.83)
    finalRatioHist.SetMaximum(1.3)
    finalRatioHist.SetMinimum(0.7)
    # finalRatioHist.GetYaxis().SetRangeUser(0.50,1.50)
    finalRatioHist.GetYaxis().SetRangeUser(0.0, 2.1)
    # finalRatioHist.GetYaxis().SetRangeUser(0.50,2.50)

    finalRatioHist.GetXaxis().SetLabelSize(0.15)

    finalRatioHist.GetXaxis().SetTitle(variableAxisTitleDictionary[variable])
    # finalRatioHist.GetXaxis().SetTitleSize(0.14)
    finalRatioHist.GetXaxis().SetTitleSize(0.17)

    MCErrors = ratioHist.Clone()
    MCErrors.Reset()
    for i in range(1, MCErrors.GetNbinsX() + 1):
        MCErrors.SetBinContent(i, 1.0)
        try:
            MCErrors.SetBinError(
                i,
                denominatorHistos.GetBinError(i) /
                denominatorHistos.GetBinContent(i))
        except ZeroDivisionError:
            MCErrors.SetBinError(i, 0)
    MCErrors.SetFillStyle(3008)
    # MCErrors.SetFillStyle(1001)
#    MCErrors.SetFillColor(ROOT.TColor.GetColor("#FFE599"))
    MCErrors.SetFillColor(ROOT.TColor.GetColor("#545252"))
    MCErrors.SetMarkerStyle(0)

    # Draw arrows for bins exceeding the range
    arrows = []
    for i in range(1, finalRatioHist.GetNbinsX() + 1):
        if finalRatioHist.GetBinContent(i) > 2.0:
            x = finalRatioHist.GetXaxis().GetBinCenter(i)
            arrow = ROOT.TArrow(x, 1.95, x, 2.0, 0.02, "|>")
            arrow.SetLineColor(ROOT.kBlack)
            arrow.SetFillColor(ROOT.kBlack)
            arrow.SetLineWidth(1)
            arrows.append(arrow)

    return finalRatioHist, MCErrors, arrows


def main():

    parser = argparse.ArgumentParser(
        description='Generate control plots quick.')
    parser.add_argument(
        '--year',
        nargs='?',
        choices=[
            '2016',
            '2017',
            '2018',
            'test',
            '2016tot',
            '2016APV'],
        help='Use the file\'s fake factor weightings when making plots for these files.',
        required=True)
    parser.add_argument('--batchMode',
                        help='run in batch mode',
                        action='store_true')

    parser.add_argument('--variables',
                        nargs='+',
                        help='Variables to draw the control plots for',
                        default=["FatJet_pt_nom[index_gFatJets[0]]",
                                 "FatJet_eta[index_gFatJets[0]]",
                                 "FatJet_phi[index_gFatJets[0]]",
                                 "FatJet_mass_nom[index_gFatJets[0]]",
                                 "FatJet_msoftdrop_nom[index_gFatJets[0]]",
                                 "FatJet_particleNetLegacy_mass[index_gFatJets[0]]",
                                 "Electron_pt[index_gElectrons[0]]",
                                 "Electron_eta[index_gElectrons[0]]",
                                 "Electron_phi[index_gElectrons[0]]",
                                 "Muon_pt[index_gMuons[0]]",
                                 "Muon_eta[index_gMuons[0]]",
                                 "Muon_phi[index_gMuons[0]]",
                                 "Tau_pt[index_gTaus]",
                                 "Tau_eta[index_gTaus]",
                                 "Tau_phi[index_gTaus]",
                                 "Tau_pt[index_gTaus[0]]",
                                 "Tau_eta[index_gTaus[0]]",
                                 "Tau_phi[index_gTaus[0]]",
                                 "Tau_pt[index_gTaus[1]]",
                                 "Tau_eta[index_gTaus[1]]",
                                 "Tau_phi[index_gTaus[1]]",
                                 "boostedTau_pt[index_gboostedTaus]",
                                 "boostedTau_eta[index_gboostedTaus]",
                                 "boostedTau_phi[index_gboostedTaus]",
                                 "boostedTau_pt[index_gboostedTaus[0]]",
                                 "boostedTau_eta[index_gboostedTaus[0]]",
                                 "boostedTau_phi[index_gboostedTaus[0]]",
                                 "boostedTau_pt[index_gboostedTaus[1]]",
                                 "boostedTau_eta[index_gboostedTaus[1]]",
                                 "boostedTau_phi[index_gboostedTaus[1]]",
                                 "METcorrected_pt",
                                 "METcorrected_phi",
                                 "ngood_Jets",
                                 "ngood_MediumJets",
                                 "ngood_TightJets",
                                 "Jet_pt_nom[index_gJets[0]]",
                                 "Jet_eta[index_gJets[0]]",
                                 "Jet_phi[index_gJets[0]]",
                                 # "HT",
                                 # "MT",
                                 "HTTvis_deltaR",
                                 "HTT_pt",
                                 "HTT_m",
                                 "HTT_phi",
                                 "HTT_eta",
                                 "X_pt",
                                 "X_m",
                                 "X_phi",
                                 "X_eta",
                                 "HTTvis_m",
                                 # "boostedTau_rawDeepTau2018v2p7VSjet[index_gboostedTaus]",
                                 # "Fatjet_pnet_bbvsqcd",
                                 # "Hbb_lep1_deltaR",
                                 # "Hbb_lep2_deltaR",
                                 # "Hbb_met_phi",
                                 # "HTTvis_deltaR_HPS",
                                 # "HTTvis_deltaR_Boosted",
                                 # "Hbb_lep1_deltaR_HPS",
                                 # "Hbb_lep1_deltaR_Boosted",
                                 # "Hbb_lep2_deltaR_HPS",
                                 # "Hbb_lep2_deltaR_Boosted",
                                 # "Hbb_met_phi_HPS",
                                 # "Hbb_met_phi_Boosted"
                                 ])

    # "Electron_pt[1]",
    # "Muon_mvaId",
    # "Muon_pt",
    # "Muon_eta",
    # "Muon_phi",
    # "boostedTau_decayMode",
    # "boostedTau_idAntiEle2018",
    # "Muon_pt[0]",
    # "Muon_pt[1]",
    # "nFatJet",
    # "FatJet_pt",
    # "FatJet_pt[0]",
    # "FatJet_pt[1]",
    # "FatJet_eta",
    # "FatJet_eta[0]",
    # "FatJet_phi",
    # "FatJet_eta[1]",
    # "FatJet_msoftdrop",
    # "DeltaR_LL",
    # "MVis_LL",
    # "FatJet_msoftdrop[0]"]
    # "FatJet_msoftdrop[1]",
    # "FatJet_particleNet_HbbvsQCD",
    # "FatJet_particleNet_HbbvsQCD[0]",
    # "FatJet_particleNet_HbbvsQCD[1]",
    # "FatJet_particleNetMD_Xbb",
    # "FatJet_particleNetMD_Xbb[0]",
    # "FatJet_particleNetMD_Xbb[1]",
    # "FatJet_tau2/FatJet_tau1",
    # "FatJet_tau2[0]/FatJet_tau1[0]",
    # "FatJet_tau2[1]/FatJet_tau1[1]"]
    # )
    # default=['Tau_pt',
    #       'Tau_phi',
    #       'Tau_eta',
    #       'nboostedTau',
    #       'nTau',
    #       'FatJet_pt',
    #       'FatJet_phi',
    #       'FatJet_eta',
    #       'nFatJet',
    #       'Electron_pt',
    #       'Electron_phi',
    #       'Electron_eta',
    #       'nElectron',
    #       'MET_pt',
    #       'MET_phi',
    #       'MET_sumEt',
    #       'Muon_pt',
    #       'Muon_eta',
    #       'Muon_phi',
    #       'nMuon'])

    parser.add_argument('--additionalSelections', '-C2',
                        nargs='+',
                        help='additional region selections',
                        # default=['Tau_idMVAoldDM2017v2 & 4 == 4','nTau==2 || nboostedTau==2','FatJet_btagDeepB > 0.45','nFatJet == 1'])
                        # default=["PV_ndof > 4", "abs(PV_z) < 24","sqrt(PV_x*PV_x+PV_y*PV_y) < 2",
                        # "Flag_goodVertices",
                        # "Flag_globalSuperTightHalo2016Filter",
                        # "Flag_HBHENoiseIsoFilter",
                        # "Flag_HBHENoiseFilter",
                        # "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        # "Flag_BadPFMuonFilter",
                        # "Flag_eeBadScFilter"])
                        # ,"gMVis_LL>0","(gFatJet_particleNetMD_Xbb / (gFatJet_particleNetMD_Xbb + gFatJet_particleNetMD_QCD))>=0.87"
                        # default=["gDeltaR_LL<1.5","gDeltaR_LL>0.05","gMVis_LL>0","fastMTT_RadionLeg_m>800","fastMTT_RadionLeg_m<4750"])
                        # default=["gDeltaR_LL<1.5","gDeltaR_LL>0","fastMTT_RadionLeg_m>=500","fastMTT_RadionLeg_m<=5000"])
                        # default=["gDeltaR_LL<1.5","gDeltaR_LL>0","VisRadion_m>=500","VisRadion_m<=5000"])
                        # default=["RecoGenRadion_mass_pnet>450","RecoGenRadion_mass_pnet<4750"])
                        # default=["RecoGenRadion_mass_pnet>=700","RecoGenRadion_mass_pnet<=4700"])
                        # default=["fastMTT_RadionLeg_m>800","fastMTT_RadionLeg_m<4750"])
                        default=["X_m>750", "X_m<5500"])
    # default=["X_m>800"])
    parser.add_argument(
        '--pause',
        help='pause after drawing each plot to make it easier to view',
        action='store_true')
    parser.add_argument('--standardCutString', '-C1',
                        nargs='?',
                        help='Change the standard cutting definition',
                        default="")
    parser.add_argument(
        '--changeHistogramBounds',
        nargs='?',
        help='Change the standard histogram bounding (affects all histograms)')

    parser.add_argument(
        '--logScale',
        help='make log plots',
        action='store_true')
    # parser.add_argument('--data',help='to include data',action='store_data')
    # parser.add_argument('--massPoint', choices=["1000","2000","3000","3500","4000"], help='to include signal sample in the plot', required=True)
    parser.add_argument(
        '--Channel',
        choices=[
            "tt",
            "et",
            "mt",
            "all",
            "lt"],
        required=True)
    parser.add_argument('--Sub', required=True)

    parser.add_argument('--Path', help='path to the files', required=True)
    parser.add_argument(
        '--Weight',
        help='weight to be added to MC',
        default='FinalWeighting')
    parser.add_argument(
        '--topRW',
        help='Scaling the top MC by a flat factor',
        action='store_true')
    parser.add_argument(
        '--computeTopWt',
        help='Compute the top MC flat factor',
        action='store_true')
    parser.add_argument(
        '--storeshape',
        help='Store the shape in root file',
        action='store_true')
    parser.add_argument(
        '--controlregion',
        '-c',
        help='If it is control region then we donot plot signals',
        action='store_true')

    args = parser.parse_args()

# if(args.controlregion):
# from variables_topCRwider import *
# else:
# from variables_topCR import *

    ROOT.gStyle.SetOptStat(0)
    if args.batchMode:
        ROOT.gROOT.SetBatch(ROOT.kTRUE)
    # change the standard cut definition if that's available

    if args.year == '2016':
        # dataPath = '/data/gparida/Background_Samples/bbtautauAnalysis/2016/ChannelFiles_Camilla/'
        # dataPath = '/data/gparida/Background_Samples/bbtautauAnalysis/2016/ChannelFiles_Camilla_28Jan_2022/'
        dataPath = args.Path
    elif args.year == '2016tot':
        dataPath = args.Path
    elif args.year == '2016APV':
        dataPath = args.Path
    elif args.year == '2017':
        dataPath = args.Path
    elif args.year == '2018':
        dataPath = args.Path
    elif args.year == 'test':
        # dataPath = '/afs/hep.wisc.edu/home/parida/HHbbtt_Analysis_Scripts/Plotting_Scripts/'
        dataPath = "/hdfs/store/user/parida/HHbbtt_Background_Files/Andrew_Script_Skim/"

    # Open all the files that are necessary for plotting .....................
    # for index in range(len(DatasetNameList))

    # fnames = glob.glob(args.Path + "/*.root")
    fnames = list(set(glob.glob(args.Path + "/*.root")) - set(glob.glob(args.Path +
                  "/*Graviton*.root")) - set(glob.glob(args.Path + "/*RadionToHHTo2B2*.root")))
    print(fnames)
    DatasetNameList = []
    DYlow_HistoList = []
    DY_HistoList = []
    ST_HistoList = []
    QCD_HistoList = []
    WJets_HistoList = []
    TT_HistoList = []
    TTHad_HistoList = []
    TTSem_HistoList = []
    TT2L2Nu_HistoList = []
    DiBoson_HistoList = []
    Other_HistoList = []
    SignalNameList = []
    Data_HistoList = []
    # SignalToPlot = []

    for file in fnames:
        nameStrip = file.strip()
        filename = (nameStrip.split('/')[-1]).split('.')[-2]
        if (not search("RadionTohhTohtatahbb", filename)):
            DatasetNameList.append(filename)
        else:
            SignalNameList.append(filename)
            # if (search(args.massPoint,filename)):
            #    SignalToPlot = filename

        if search("TTT", filename):
            TT_HistoList.append(filename)
        if search("TTToHadronic", filename):
            TTHad_HistoList.append(filename)
        if search("TTToSemiLep", filename):
            TTSem_HistoList.append(filename)
        if search("TTTo2L2Nu", filename):
            TT2L2Nu_HistoList.append(filename)
        if ((search("DY", filename) and search("M-10to50", filename))
                or (search("DY", filename) and search("M-4to50", filename))):
            DYlow_HistoList.append(filename)
        if ((search("DY", filename))):
            if ((not search("M-10to50", filename))
                    and ((not search("M-4to50", filename)))):
                DY_HistoList.append(filename)
        if search("WJet", filename):
            WJets_HistoList.append(filename)
        if search("ST_", filename):
            ST_HistoList.append(filename)
        if search("QCD", filename):
            QCD_HistoList.append(filename)
        if (search("WW", filename) or search(
                "WZ", filename) or search("ZZ", filename)):
            DiBoson_HistoList.append(filename)
        if (search("MET", filename) or search("Run", filename)):
            Data_HistoList.append(filename)

    Other_HistoList = QCD_HistoList + DiBoson_HistoList + ST_HistoList
    # Other_HistoList=  DiBoson_HistoList + ST_HistoList
    # print ("Datasets = ",DatasetNameList)
    # print()
    print(("Signal Samples =", SignalNameList))
    print()
    print(("DYlow_HistoList = ", DYlow_HistoList))
    print()
    print(("DY_HistoList =", DY_HistoList))
    print()
    print(("WJets_HistoList =", WJets_HistoList))
    print()
    print(("TT_HistoList =", TT_HistoList))
    print()
    print(("TTHad_HistoList =", TTHad_HistoList))
    print()
    print(("TTSem_HistoList =", TTSem_HistoList))
    print()
    print(("TT2L2Nu_HistoList =", TT2L2Nu_HistoList))
    print()
    print(("DiBoson_HistoList =", DiBoson_HistoList))
    print()
    print(("Other_HistoList =", Other_HistoList))
    print()
    print(("SignalNameList =", SignalNameList))
    print()
    print("SignalToPlot 1TeV, 2TeV, 3TeV, 4Tev")
    print()
    print(("Data = ", Data_HistoList))
    print()
    # print DatasetNameList

    logScaleVarList = ["Fatjet_pnet_bbvsqcd", "X_pt", "X_m",
                       "FatJet_pt_nom[index_gFatJets[0]]", "METcorrected_pt"]

    ##########################################################################
    # If we want to store the shapes then open and recreate a root file
    if (args.storeshape):
        shape_file = ROOT.TFile(
            args.Channel.upper() +
            "Plots/" +
            args.Sub +
            '/' +
            "shaperoot_" +
            args.Channel +
            ".root",
            "RECREATE")
    # For loop to draw histograms
    for variable in args.variables:
        try:
            variableSettingDictionary[variable] is not None
        except KeyError:
            print(("No defined histogram settings for variable: " + variable))
            continue
        try:
            variableAxisTitleDictionary[variable]
        except KeyError:
            print(("No defined title information for variable: " + variable))
            continue

        if args.changeHistogramBounds is not None:
            variableSettingDictionary[variable] = args.changeHistogramBounds

        save_path = MYDIR = os.getcwd() + "/countingData"
        # file_name = "Entries_MTChannel.txt"
        # complete_Name =  os. path. join(save_path, file_name)
        # file = open(complete_Name,"a")
        # file.write("variable "+'\t'+ "Data Counts"+'\n')

        #### Drawing the Histograms#######
        DatasetObjects = {}
        for index in range(len(DatasetNameList)):
            # DatasetObjects[DatasetNameList[index]]=MakeHistograms(dataPath,DatasetNameList[index],str(DatasetNameXSWeightDictionary[DatasetNameList[index]]))
            DatasetObjects[DatasetNameList[index]] = MakeHistograms(
                dataPath, DatasetNameList[index])

        for index in range(len(DatasetNameList)):
            # print DatasetNameList[index]
            # if DatasetNameList[index] == "Data":
            if (search("MET", DatasetNameList[index]) or search(
                    "Run", DatasetNameList[index])):
                DatasetObjects[DatasetNameList[index]].StandardDraw(DatasetObjects[DatasetNameList[index]].RootFileName,
                                                                    variable,
                                                                    # args.standardCutString+"&&(run>=319077)",
                                                                    # "(channel==0)&&(HTTvis_deltaR<1.5)&&(run<319077)",
                                                                    args.standardCutString,
                                                                    args.additionalSelections,
                                                                    DatasetNameList[index], theWeight='1')
                # DatasetNameList[index],theWeight='1'+"*((FatJet_particleNetLegacy_mass[index_gFatJets[0]]
                # > 150) || (FatJet_particleNetLegacy_mass[index_gFatJets[0]] <
                # 100))")

            # Temporary protection for DY-M50 and WJets
            elif (search("DYJetsToLL_M-50", DatasetNameList[index])):
                # print("This is high mass DY : 1.23 weight  -->",DatasetNameList[index])
                DatasetObjects[DatasetNameList[index]].StandardDraw(DatasetObjects[DatasetNameList[index]].RootFileName,
                                                                    variable,
                                                                    args.standardCutString,
                                                                    args.additionalSelections,
                                                                    # DatasetNameList[index],theWeight=args.Weight+"*1.23")
                                                                    DatasetNameList[index], theWeight=args.Weight)
            elif (search("WJets", DatasetNameList[index])):
                # print("This is WJets : 1.21 weight  -->",DatasetNameList[index])
                DatasetObjects[DatasetNameList[index]].StandardDraw(DatasetObjects[DatasetNameList[index]].RootFileName,
                                                                    variable,
                                                                    args.standardCutString,
                                                                    args.additionalSelections,
                                                                    # DatasetNameList[index],theWeight=args.Weight+"*1.21")
                                                                    DatasetNameList[index], theWeight=args.Weight)
            else:
                DatasetObjects[DatasetNameList[index]].StandardDraw(DatasetObjects[DatasetNameList[index]].RootFileName,
                                                                    variable,
                                                                    args.standardCutString,
                                                                    args.additionalSelections,
                                                                    DatasetNameList[index], theWeight=args.Weight)
                # DatasetNameList[index],theWeight=args.Weight+"*((FatJet_particleNetLegacy_mass[index_gFatJets[0]] > 150) || (FatJet_particleNetLegacy_mass[index_gFatJets[0]] < 100))")
                # DatasetObjects[DatasetNameList[index]].userWeight)
                # DatasetObjects[DatasetNameList[index]].FillEvents((DatasetObjects[DatasetNameList[index]].RootFileName),DatasetNameList[index])

        SignalObjects = {}
        for index in range(len(SignalNameList)):
            SignalObjects[SignalNameList[index]] = MakeHistograms(
                dataPath, SignalNameList[index])

        for index in range(len(SignalNameList)):
            print((SignalNameList[index]))
            # weightSig = (1.0/50000)
            SignalObjects[SignalNameList[index]].StandardDraw(SignalObjects[SignalNameList[index]].RootFileName,
                                                              variable,
                                                              args.standardCutString,
                                                              args.additionalSelections,
                                                              # SignalNameList[index],theWeight=str(weightSig)+'*16.81'+'*genWeight*pileupWeighting')
                                                              SignalNameList[index], theWeight=args.Weight)
            # SignalNameList[index],theWeight=args.Weight+"*((FatJet_particleNetLegacy_mass[index_gFatJets[0]]
            # > 150) || (FatJet_particleNetLegacy_mass[index_gFatJets[0]] <
            # 100))")

        ######################## Signal-Histogram#############################
        Signal_Histo = [
            SignalObjects["RadionTohhTohtatahbb_narrow_M-1000"].HistogramName.Clone(),
            SignalObjects["RadionTohhTohtatahbb_narrow_M-2000"].HistogramName.Clone(),
            SignalObjects["RadionTohhTohtatahbb_narrow_M-3000"].HistogramName.Clone(),
            SignalObjects["RadionTohhTohtatahbb_narrow_M-4000"].HistogramName.Clone()]
        #####################################################################

    ########################### CLUB ALL HISTS################################
        totalbkg_Histo = clubHistograms(
            DY_HistoList +
            DYlow_HistoList +
            ST_HistoList +
            TT_HistoList +
            QCD_HistoList +
            WJets_HistoList +
            DiBoson_HistoList,
            DatasetObjects)
        ####################### DYlow-Histograms###############################

        # DYlow_Histo = clubHistograms(DYlow_HistoList,DatasetObjects)

        #######################################################################

        ####################### DY-Histograms##################################

        DY_Histo = clubHistograms(
            DY_HistoList + DYlow_HistoList,
            DatasetObjects)

        #######################################################################
#
        ############################ ST-Histograms#############################
        ST_Histo = clubHistograms(ST_HistoList, DatasetObjects)

        # PF_ST_Histo = DatasetObjects["ST_s-channel_4f"].PassFailHistogramName.Clone()
        #######################################################################
#
        ############################## QCD-Histograms##########################
        print("\n\n\n\n Clubbing QCD")
        QCD_Histo = clubHistograms(QCD_HistoList, DatasetObjects)
        print("\n\n\n\n")

        #######################################################################

        ################################## WJets###############################
        WJets_Histo = clubHistograms(WJets_HistoList, DatasetObjects)

        #######################################################################
#
        #################################### TT-Histograms#####################
        TT_Histo = clubHistograms(TT_HistoList, DatasetObjects)
        TTHad_Histo = clubHistograms(TTHad_HistoList, DatasetObjects)
        TTSem_Histo = clubHistograms(TTSem_HistoList, DatasetObjects)
        TT2L2Nu_Histo = clubHistograms(TT2L2Nu_HistoList, DatasetObjects)
        #######################################################################
#
        ################################ DiBoson-Histograms####################
        DiBoson_Histo = clubHistograms(DiBoson_HistoList, DatasetObjects)

        #######################################################################

        #################################### Combine=ing Backgrounds###########
        Other_Histo = clubHistograms(Other_HistoList, DatasetObjects)
        # new_binning = array('d', [0,(100*6/3),(100*18/3), 1500])
        # Other_Histo = Other_Histo.Rebin(3, '', new_binning ) # for custom
        # binning

        # print ("Number of bins in new histogram = ",Other_Histo.GetNbinsX())
        ################################ Data is represented as points#########

        ################################ Data is represented as points#########

        ######### CLUB Data Histograms################################
        Data_Histo = clubHistograms(Data_HistoList, DatasetObjects)

        # DatasetObjects["Data"].HistogramName.SetMarkerStyle(20)

        Data_Histo.SetMarkerStyle(20)
        Data_Histo.SetMarkerSize(0.7)
        # Data_Histo.Sumw2()
        Data_Histo.SetBinErrorOption(ROOT.TH1.kPoisson)
        # Counting the events contribut

        ########################################### TOP ReWeighting############
        if (args.Channel == "tt" and args.topRW):
            if args.year == "2016":
                # Weights with W,Zpt + HbbLooseTagger and weight applied
                TT_Histo.Scale(1.425416014061129)
                ST_Histo.Scale(1.425416014061129)
                # Weights with W and Z pt reweighting
                # TT_Histo.Scale(1.089220893211982)
                # ST_Histo.Scale(1.089220893211982)
                # Weights without W and Z pt reweighting
                # TT_Histo.Scale(1.10127254991914)
                # ST_Histo.Scale(1.10127254991914)
        elif (args.Channel == "et" and args.topRW):
            if args.year == "2016":
                # Weights with W,Zpt + HbbLooseTagger and weight applied
                TT_Histo.Scale(0.7412743255675137)
                ST_Histo.Scale(0.7412743255675137)
                # Weights with W and Z pt reweighting
                # TT_Histo.Scale(0.9586445902878492)
                # ST_Histo.Scale(0.9586445902878492)
                # Weights without W and Z pt reweighting
                # TT_Histo.Scale(0.9682605849925042)
                # ST_Histo.Scale(0.9682605849925042)
        elif (args.Channel == "mt" and args.topRW):
            if args.year == "2016":
                # Weights with W,Zpt + HbbLooseTagger and weight applied
                TT_Histo.Scale(0.8804393603576028)
                ST_Histo.Scale(0.8804393603576028)
                # Weights with W and Z pt reweighting
                # TT_Histo.Scale( 0.8729334376959164)
                # ST_Histo.Scale( 0.8729334376959164)
                # Weights without W and Z pt reweighting
                # TT_Histo.Scale(0.8827564525729762)
                # ST_Histo.Scale(0.8827564525729762)

#        print ("Number of events in TTBar = ",TT_Histo.Integral(TT_Histo.GetXaxis().FindBin(900),TT_Histo.GetXaxis().FindBin(5000)))
#        print ("Number of events in WJets = ",WJets_Histo.Integral(WJets_Histo.GetXaxis().FindBin(900),WJets_Histo.GetXaxis().FindBin(5000)))
#        print ("Number of events in DY = ",DY_Histo.Integral(DY_Histo.GetXaxis().FindBin(900),DY_Histo.GetXaxis().FindBin(5000)))
#        print ("Number of events in ST = ",ST_Histo.Integral(ST_Histo.GetXaxis().FindBin(900),ST_Histo.GetXaxis().FindBin(5000)))
#        print ("Number of events in QCD = ",QCD_Histo.Integral(QCD_Histo.GetXaxis().FindBin(900),QCD_Histo.GetXaxis().FindBin(5000)))
#        print ("Number of Diboson = ",DiBoson_Histo.Integral(DiBoson_Histo.GetXaxis().FindBin(900),DiBoson_Histo.GetXaxis().FindBin(5000)))
#        print("Total Background = ",TT_Histo.Integral(TT_Histo.GetXaxis().FindBin(900),TT_Histo.GetXaxis().FindBin(5000))+WJets_Histo.Integral(WJets_Histo.GetXaxis().FindBin(900),WJets_Histo.GetXaxis().FindBin(5000))+DY_Histo.Integral(DY_Histo.GetXaxis().FindBin(900),DY_Histo.GetXaxis().FindBin(5000))+Other_Histo.Integral(Other_Histo.GetXaxis().FindBin(900),Other_Histo.GetXaxis().FindBin(5000)))
#        print ("Number of events for the 4 TeV",Signal_Histo[3].Integral(Signal_Histo[3].GetXaxis().FindBin(900),Signal_Histo[3].GetXaxis().FindBin(5000)))
#        print ("Number of events of Data = ",DatasetObjects["Data"].HistogramName.Integral(DatasetObjects["Data"].HistogramName.GetXaxis().FindBin(900),DatasetObjects["Data"].HistogramName.GetXaxis().FindBin(5000)))

#        print ("Number of events in TTBar = ",TT_Histo.Integral())
#        print ("Number of events in WJets = ",WJets_Histo.Integral())
#        print ("Number of events in DY = ",DY_Histo.Integral())
#        print ("Number of events in lowMassDY = ",DYlow_Histo.Integral())
#        print ("Number of events in ST = ",ST_Histo.Integral())
#        print ("Number of events in QCD = ",QCD_Histo.Integral())
#        print ("Number of Diboson = ",DiBoson_Histo.Integral())
#        print ("Total Background = ",TT_Histo.Integral()+WJets_Histo.Integral()+DY_Histo.Integral()+ST_Histo.Integral()+QCD_Histo.Integral()+DiBoson_Histo.Integral()+DYlow_Histo.Integral())
#        print ("Number of events for the 1 TeV",Signal_Histo[0].Integral())
#        print ("Number of events for the 2 TeV",Signal_Histo[1].Integral())
#        print ("Number of events for the 3 TeV",Signal_Histo[2].Integral())
#        print ("Number of events for the 4 TeV",Signal_Histo[3].Integral())
#        print ("Number of events of Data = ",Data_Histo.Integral())

        print(("Number of events in TTBar (incl over and under) = ",
              TT_Histo.Integral(0, TT_Histo.GetNbinsX() + 1)))
        print("TTbar Component Split : >>>>>>>>>>>>>>>>>>")
        print(("Number of events in TTBarHadronic (incl over and under) = ",
              TTHad_Histo.Integral(0, TTHad_Histo.GetNbinsX() + 1)))
        print(("Number of events in TTBarSemileptonic (incl over and under) = ",
              TTSem_Histo.Integral(0, TTSem_Histo.GetNbinsX() + 1)))
        print(("Number of events in TTBar2L2Nu (incl over and under) = ",
              TT2L2Nu_Histo.Integral(0, TT2L2Nu_Histo.GetNbinsX() + 1)))
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print(("Number of events in WJets (incl over and under) = ",
              WJets_Histo.Integral(0, WJets_Histo.GetNbinsX() + 1)))
        print(("Number of events in DY + DY-low (incl over and under) = ",
              DY_Histo.Integral(0, DY_Histo.GetNbinsX() + 1)))
        # print ("Number of events in lowMassDY (incl over and under) = ",DYlow_Histo.Integral(0,DYlow_Histo.GetNbinsX()+1))
        print(("Number of events in ST (incl over and under) = ",
              ST_Histo.Integral(0, ST_Histo.GetNbinsX() + 1)))
        print(("Number of events in QCD (incl over and under) = ",
              QCD_Histo.Integral(0, QCD_Histo.GetNbinsX() + 1)))
        print(("Number of Diboson (incl over and under) = ",
              DiBoson_Histo.Integral(0, DiBoson_Histo.GetNbinsX() + 1)))
        print(("Total Background (incl over and under) = ",
               TT_Histo.Integral(0,
                                 TT_Histo.GetNbinsX() + 1) + WJets_Histo.Integral(0,
                                                                                  WJets_Histo.GetNbinsX() + 1) + DY_Histo.Integral(0,
                                                                                                                                   DY_Histo.GetNbinsX() + 1) + ST_Histo.Integral(0,
                                                                                                                                                                                 ST_Histo.GetNbinsX() + 1) + QCD_Histo.Integral(0,
                                                                                                                                                                                                                                QCD_Histo.GetNbinsX() + 1) + DiBoson_Histo.Integral(0,
                                                                                                                                                                                                                                                                                    DiBoson_Histo.GetNbinsX() + 1)))
        print(("Number of events for the 1 TeV (incl over and under) = ",
              Signal_Histo[0].Integral(0, Signal_Histo[0].GetNbinsX() + 1)))
        print(("Number of events for the 2 TeV (incl over and under) = ",
              Signal_Histo[1].Integral(0, Signal_Histo[1].GetNbinsX() + 1)))
        print(("Number of events for the 3 TeV (incl over and under) = ",
              Signal_Histo[2].Integral(0, Signal_Histo[2].GetNbinsX() + 1)))
        print(("Number of events for the 4 TeV (incl over and under) = ",
              Signal_Histo[3].Integral(0, Signal_Histo[3].GetNbinsX() + 1)))
        print(("Number of events of Data (incl over and under) = ",
              Data_Histo.Integral(0, Data_Histo.GetNbinsX() + 1)))

        print("\n\n\n >>>>>>>>>>>>>>>>>>>Error and Integral<<<<<<<<<<<<<<<<<<<<<<<<<")
        error_bkg = ROOT.Double(0)
        integral_bkg = totalbkg_Histo.IntegralAndError(
            0, totalbkg_Histo.GetNbinsX() + 1, error_bkg, "")
        print(("Background integral = {} +- {}".format(integral_bkg, error_bkg)))
        error_data = ROOT.Double(0)
        integral_data = Data_Histo.IntegralAndError(
            0, Data_Histo.GetNbinsX() + 1, error_data, "")
        print(("Data integral = {} +- {}".format(integral_data, error_data)))
        print(">>>>>>>>>>>>>>>>>>>Error and Integral<<<<<<<<<<<<<<<<<<<<<<<<<\n\n\n")

        if (args.storeshape):
            shape_hist = clubHistograms(
                TT_HistoList +
                ST_HistoList +
                DY_HistoList +
                DYlow_HistoList +
                QCD_HistoList +
                WJets_HistoList +
                DiBoson_HistoList,
                DatasetObjects)
            shape_hist.Sumw2()
            shape_hist.Scale((1.0 / (TT_Histo.Integral(0,
                                                       TT_Histo.GetNbinsX() + 1) + WJets_Histo.Integral(0,
                                                                                                        WJets_Histo.GetNbinsX() + 1) + DY_Histo.Integral(0,
                                                                                                                                                         DY_Histo.GetNbinsX() + 1) + ST_Histo.Integral(0,
                                                                                                                                                                                                       ST_Histo.GetNbinsX() + 1) + QCD_Histo.Integral(0,
                                                                                                                                                                                                                                                      QCD_Histo.GetNbinsX() + 1) + DiBoson_Histo.Integral(0,
                                                                                                                                                                                                                                                                                                          DiBoson_Histo.GetNbinsX() + 1))))
            shape_hist.Sumw2()
            shape_hist.SetName(variable)
            shape_hist.SetTitle(variable + ";;")
            shape_file.cd()
            shape_hist.Write()

        #######################################################################
        if (args.computeTopWt):
            # Now compute the uncertanity
            Data_fe = Data_Histo.Clone()
            Data_fe.Sumw2()
            # Data_fe.Rebin(5)
            MC_fe = TT_Histo.Clone()
            MC_fe.Add(WJets_Histo)
            MC_fe.Add(DY_Histo)
            MC_fe.Add(ST_Histo)
            MC_fe.Add(QCD_Histo)
            MC_fe.Add(DiBoson_Histo)
            # MC_fe.Rebin(5)

            Top_fe = TT_Histo.Clone()
            Top_fe.Add(ST_Histo)
            # Top_fe.Rebin(5)

            dataErr = ctypes.c_double(-999.99)
            TopErr = ctypes.c_double(-999.99)
            MCErr = ctypes.c_double(-999.99)
            Data_Integral = Data_fe.IntegralAndError(0, 5, dataErr)
            Top_Integral = Top_fe.IntegralAndError(0, 5, TopErr)
            MC_Integral = MC_fe.IntegralAndError(0, 5, MCErr)
            print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>REFINED TOP REWEIGHTING FACTOR>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(("Data = ", Data_Integral))
            print(("MC Total = ", MC_Integral))
            print(("Total TT+ST = ", Top_Integral))
            Top_Weight = (
                (Data_Integral -
                 MC_Integral +
                 Top_Integral) /
                (Top_Integral))
            Numerator_fe = (Data_Integral - MC_Integral)
            ratio_fe = (Data_Integral - MC_Integral) / Top_Integral

            print(("Correction factor for channel ",
                  args.Channel, " is = ", Top_Weight))
            print(("dataErr = ", dataErr.value))
            print(("MCErr = ", MCErr.value))
            print(("TopErr = ", TopErr.value))
            Top_Weight_Error = math.sqrt(pow(((ratio_fe * dataErr.value) / (Numerator_fe)), 2) + pow(
                ((ratio_fe * MCErr.value) / (Numerator_fe)), 2) + pow(((ratio_fe * TopErr.value) / (Top_Integral)), 2))
            print(("Propagated Error = ", Top_Weight_Error))

            # a=99.99

            # a= ctypes.c_double(-999.99)
            # a=ROOT.Double()
            # print(TT_Histo.IntegralAndError(0,5,a))
            # print("Error = ",a)

        ################################ Color_Definitions -- Background Fill##
#        color_DiBoson="#CCE4FC"
#        color_TT="#87997C"
#        color_WJets="#65839F"
#        color_QCD="#EEA287"
#        color_ST="#97D870"
#        color_DY="#9e62ee"
#        color_DYlow="#cfb7ef"
#        color_other = "#923814"

        color_TT = "#92cfe0"
        color_ST = "#a5e7fa"
        color_WJets = "#fcd068"
        color_DY = "#d8ed79"
        color_DYlow = "#83e01f"
        color_QCD = "#f29b6f"
        color_DiBoson = "#9d99bd"
        color_other = "#923814"

        # color_jetfake="#f1cde1"

        ################################# Filling Color for Backgrounds########

        # ST_s_channel_4f.SetFillColor(ROOT.TColor.GetColor("#ffcc66"))
        signalXS = 1.0
        # signalXS = 0.000001
        Signal_Histo[0].SetLineColor(ROOT.kRed)
        Signal_Histo[0].Scale(signalXS)
        Signal_Histo[0].SetLineWidth(1)

        Signal_Histo[1].SetLineColor(ROOT.kRed + 2)
        Signal_Histo[1].Scale(signalXS)
        Signal_Histo[1].SetLineWidth(1)

        Signal_Histo[2].SetLineColor(ROOT.kRed + 3)
        Signal_Histo[2].Scale(signalXS)
        Signal_Histo[2].SetLineWidth(1)

        Signal_Histo[3].SetLineColor(ROOT.kRed + 4)
        Signal_Histo[3].Scale(signalXS)
        Signal_Histo[3].SetLineWidth(1)

        DiBoson_Histo.SetFillColor(ROOT.TColor.GetColor(color_DiBoson))
        TT_Histo.SetFillColor(ROOT.TColor.GetColor(color_TT))
        WJets_Histo.SetFillColor(ROOT.TColor.GetColor(color_WJets))
        QCD_Histo.SetFillColor(ROOT.TColor.GetColor(color_QCD))
        ST_Histo.SetFillColor(ROOT.TColor.GetColor(color_ST))
        DY_Histo.SetFillColor(ROOT.TColor.GetColor(color_DY))
        # if (DYlow_Histo!=None):
        # DYlow_Histo.SetFillColor(ROOT.TColor.GetColor(color_DYlow))
        Other_Histo.SetFillColor(ROOT.TColor.GetColor(color_other))

        DiBoson_Histo.SetLineWidth(0)
        TT_Histo.SetLineWidth(0)
        WJets_Histo.SetLineWidth(0)
        QCD_Histo.SetLineWidth(0)
        ST_Histo.SetLineWidth(0)
        DY_Histo.SetLineWidth(0)
        # if (DYlow_Histo!=None):
        # DYlow_Histo.SetLineWidth(0)
        Other_Histo.SetLineWidth(0)

        ######################################## Histograms For Shape Check####
        BackgroundShape = Other_Histo.Clone()
        # BackgroundShape.Add(ST_Histo)
        # BackgroundShape.Add(QCD_Histo)
        BackgroundShape.Add(DY_Histo)
        BackgroundShape.Add(WJets_Histo)
        BackgroundShape.Add(TT_Histo)

        # BackgroundShape.Add(DiBoson_Histo)
        # BackgroundShape.SetOptTitle(0)

        if (BackgroundShape.Integral() != 0):
            ScaleBackground = 1 / BackgroundShape.Integral()
        else:
            ScaleBackground = 1.0

        DataShape = Data_Histo.Clone()

        # file.write(variable + '\t' + str(DataShape.Integral())+'\n')

        if (BackgroundShape.Integral() != 0):
            print(("Tot integral = ", DataShape.Integral()))
            # ScaleData = 1/DataShape.Integral()
            ScaleData = 1
        else:
            ScaleData = 1.0
        DataShape.Scale(ScaleData)

        # Making the Stack of Histogr

        backgroundStack = ROOT.THStack('backgroundStack', 'backgroundstack')

        # This is the original order for the plotting
        backgroundStack.Add(DiBoson_Histo, 'HIST')
        backgroundStack.Add(ST_Histo, 'HIST')
        backgroundStack.Add(TT_Histo, 'HIST')
        backgroundStack.Add(QCD_Histo, 'HIST')
        backgroundStack.Add(WJets_Histo, 'HIST')
        backgroundStack.Add(DY_Histo, 'HIST')

        # For the Combine Inputs for Validation plots
        # backgroundStack.Add(DiBoson_Histo,'HIST')
        # backgroundStack.Add(QCD_Histo,'HIST')
        # backgroundStack.Add(WJets_Histo,'HIST')
        # backgroundStack.Add(DY_Histo,'HIST')
        # backgroundStack.Add(ST_Histo,'HIST')
        # backgroundStack.Add(TT_Histo,'HIST')
#
        backgroundStack_Errors = MakeStackErrors(backgroundStack)
#
        # PassFailStack_Errors = MakeStackErrors(PassFailStack)
        N_DY_Histo = (DY_Histo.Clone())
        N_ST_Histo = (ST_Histo.Clone())
        N_QCD_Histo = (QCD_Histo.Clone())
        N_WJets_Histo = (WJets_Histo.Clone())
        N_TT_Histo = (TT_Histo.Clone())
        N_DiBoson_Histo = (DiBoson_Histo.Clone())
        N_Other_Histo = (Other_Histo.Clone())

        N_DY_Histo.Scale(ScaleBackground)
        # N_ST_Histo.Scale(ScaleBackground)
        # N_QCD_Histo.Scale(ScaleBackground)
        N_WJets_Histo.Scale(ScaleBackground)
        N_TT_Histo.Scale(ScaleBackground)
        # N_DiBoson_Histo.Scale(ScaleBackground)
        N_Other_Histo.Scale(ScaleBackground)

        ShapeStack = ROOT.THStack('ShapeStack', 'ShapeStack')
        ShapeStack.Add(N_DY_Histo, 'HIST')
        # ShapeStack.Add(N_ST_Histo,'HIST')
        # ShapeStack.Add(N_QCD_Histo,'HIST')
        ShapeStack.Add(N_WJets_Histo, 'HIST')
        ShapeStack.Add(N_TT_Histo, 'HIST')
        ShapeStack.Add(N_Other_Histo, 'HIST')
        # ShapeStack.Add(N_DiBoson_Histo,'HIST')

        ShapeStack_Errors = MakeStackErrors(ShapeStack)

        ########################################## Preparing the Canvas########

        theCanvas = ROOT.TCanvas("theCanvas", "theCanvas")
        theCanvas.Divide(1, 2)

        plotPad = ROOT.gPad.GetPrimitive('theCanvas_1')
        ratioPad = ROOT.gPad.GetPrimitive('theCanvas_2')

        plotPad.SetPad("pad1", "plot", 0, 0.25, 1, 1)
        plotPad.SetFillColor(0)
        plotPad.SetBorderMode(0)
        plotPad.SetBorderSize(0)
        plotPad.SetTickx(1)
        plotPad.SetTicky(1)
        # plotPad.SetGridx()
        plotPad.SetLeftMargin(0.15)  # 0.15
        # plotPad.SetLeftMargin(0.07) #0.15
        # plotPad.SetRightMargin(0.15) #0.1
        plotPad.SetRightMargin(0.07)
        plotPad.SetTopMargin(0.14)  # 0.122 if the exponent is not present
        # plotPad.SetBottomMargin(0.025)
        plotPad.SetBottomMargin(0.00)
        plotPad.SetFrameFillStyle(0)
        plotPad.SetFrameLineStyle(0)
        plotPad.SetFrameLineWidth(1)  # 1
        plotPad.SetFrameBorderMode(0)
        plotPad.SetFrameBorderSize(1)
        # plotPad.SetErrorX()
        if ((args.logScale) or (variable in logScaleVarList)):
            plotPad.SetLogy(1)
        # plotPad.SetLogy(1)
        # plotPad.SetOptTitle(0)
        #
        ratioPad.SetPad("pad2", "ratio", 0, 0, 1, 0.25)
        ratioPad.SetFillColor(0)
        ratioPad.SetBorderSize(0)
        # ratioPad.SetTopMargin(0.02)
        ratioPad.SetTopMargin(0.00)
        ratioPad.SetBottomMargin(0.35)
        ratioPad.SetLeftMargin(0.15)
        # ratioPad.SetLeftMargin(0.07)
        # ratioPad.SetRightMargin(0.15)
        ratioPad.SetRightMargin(0.07)
        ratioPad.SetTickx(1)
        ratioPad.SetTicky(1)
        ratioPad.SetFrameLineWidth(1)
        ratioPad.SetGridy()
        # pad2.SetGridx()
#

        ratioHist, ratioError, arrows = MakeRatioHistograms(
            Data_Histo, backgroundStack, variable)
        ratioPad.cd()
        ratioHist.GetXaxis().SetTickLength(0.01)
        ratioHist.GetYaxis().SetTickLength(0.01)
        ratioHist.Draw('ex0')
        ratioHist.Print("all")
        ratioError.Draw('SAME e2')
        for arrow in arrows:
            arrow.Draw()
        # ratioHist.Draw('SAME ex0')
#
        plotPad.cd()
        plotPad.SetFrameLineWidth(1)
        # plotPad.SetTickx()
        # .SetTicky()

#
        # Delete
        maxi = max(
            backgroundStack.GetMaximum(),
            Data_Histo.GetMaximum(),
            Signal_Histo[0].GetMaximum(),
            Signal_Histo[1].GetMaximum(),
            Signal_Histo[2].GetMaximum(),
            Signal_Histo[3].GetMaximum())
        # backgroundStack.SetMaximum(maxi + 0.50*maxi)
        if variable in logScaleVarList:
            backgroundStack.SetMaximum(maxi + 100 * maxi)
        else:
            backgroundStack.SetMaximum(maxi + 0.50 * maxi)
        # SetMinimum(0.00005)
        # backgroundStack.SetMaximum(maxi + 15.0*maxi)
        # backgroundStack.SetMaximum(0.1)
        backgroundStack.SetMinimum(0.1)
        # backgroundStack.SetMinimum(0.00005)
        # backgroundStack.GetYaxis().SetRangeUser(0.001,1000)

        backgroundStack.Draw()
        backgroundStack_Errors.Draw('SAME e2')
        # backgroundStack.SetTitle(variableAxisTitleDictionary[variable])
        backgroundStack.SetTitle("")
        # if(args.data):
        # DatasetObjects["Data"].HistogramName.GetXaxis().SetErrorX()
        # Data_Histo.GetXaxis().SetTickLength(0.1)
        # Data_Histo.GetYaxis().SetTickLength(0.1)
        # Data_Histo.Draw('SAME e1')
        Data_Histo.Draw('SAME ex0')
        if (args.controlregion):
            print("This is a CR plot - NO SIGNALS in the PLOT")
        else:
            Signal_Histo[0].Draw('SAME HIST')
            Signal_Histo[1].Draw('SAME HIST')
            Signal_Histo[2].Draw('SAME HIST')
            Signal_Histo[3].Draw('SAME HIST')
        backgroundStack.GetYaxis().SetTitle("Events")
        backgroundStack.GetYaxis().SetTitleSize(0.065)
        backgroundStack.GetYaxis().SetLabelSize(0.05)
   # backgroundStack.GetYaxis().SetTitleOffset(1.58)
        # backgroundStack.GetYaxis().SetTitleOffset(0.87)
        backgroundStack.GetYaxis().SetTitleOffset(0.70)
        # backgroundStack.GetYaxis().SetTitleSize(1)
        backgroundStack.GetXaxis().SetLabelSize(0.0)
        backgroundStack.GetXaxis().SetTickLength(0.01)
        backgroundStack.GetYaxis().SetTickLength(0.01)

    ############################## Legend############################

        # theLegend = ROOT.TLegend(0.85, 0.45, 1.0, 0.75, "", "brNDC")
        ###
        # theLegend = ROOT.TLegend(0.4, 0.65, 0.85, 0.85, "", "brNDC")
        theLegend = ROOT.TLegend(0.6, 0.6, 0.85, 0.85, "", "brNDC")
        # theLegend.SetNColumns(3)
        theLegend.SetNColumns(2)
        theLegend.SetTextSize(0.08)
        # theLegend.SetTextSize(0.03)
        if (args.Channel == "tt"):
            theLegend.SetHeader("#tau-#tau Channel", "C")
        elif (args.Channel == "et"):
            theLegend.SetHeader("e-#tau Channel", "C")
        elif (args.Channel == "mt"):
            theLegend.SetHeader("#mu-#tau Channel", "C")
        elif (args.Channel == "lt"):
            theLegend.SetHeader("l-#tau Channel", "C")
        elif (args.Channel == "all"):
            theLegend.SetHeader("all Channels", "C")

        else:
            print("Enter a valid channel")

        theLegend.SetTextSize(0.03)
        # theLegend.SetColumnSeparation(0.3)
        theLegend.SetLineWidth(0)
        theLegend.SetLineStyle(1)
        theLegend.SetFillStyle(1001)  # 0
        theLegend.SetFillColor(0)
        theLegend.SetBorderSize(0)
        theLegend.SetTextFont(42)
        # if(args.data):
        theLegend.AddEntry(Data_Histo, 'Observed', 'pe')
        theLegend.AddEntry(DY_Histo, 'Drell-Yan', 'f')
        # if (DYlow_Histo!=None):
        # theLegend.AddEntry(DYlow_Histo,'Drell-Yan-low','f')
        # theLegend.AddEntry(DiBoson_Histo,'DiBoson','f')
        theLegend.AddEntry(WJets_Histo, 'WJets', 'f')
        theLegend.AddEntry(QCD_Histo, 'QCD', 'f')
        # theLegend.AddEntry(ST_Histo,'ST_s_Channel','f')
        theLegend.AddEntry(TT_Histo, 'TTbar', 'f')
        theLegend.AddEntry(ST_Histo, 'STop', 'f')
        # theLegend.AddEntry(Other_Histo,'Others','f')
        # theLegend.AddEntry(QCD_Histo,'QCD','f')
        theLegend.AddEntry(DiBoson_Histo, 'DiBoson', 'f')
        if (args.controlregion):
            print("This is a CR plot - NO SIGNALS in the LEGEND")
        else:
            theLegend.AddEntry(
                Signal_Histo[0],
                '1TeV(' + str(signalXS) + 'pb)',
                'l')
            theLegend.AddEntry(
                Signal_Histo[1],
                '2TeV(' + str(signalXS) + 'pb)',
                'l')
            theLegend.AddEntry(
                Signal_Histo[2],
                '3TeV(' + str(signalXS) + 'pb)',
                'l')
            theLegend.AddEntry(
                Signal_Histo[3],
                '4TeV(' + str(signalXS) + 'pb)',
                'l')

        theLegend.Draw('SAME')

    ##########################################################################

    ##########################################################################
        # also draw the preliminary warnings
        cmsLatex = ROOT.TLatex()
        cmsLatex.SetTextSize(0.06)
        cmsLatex.SetNDC(True)
        cmsLatex.SetTextFont(61)
        cmsLatex.SetTextAlign(11)
        # cmsLatex.DrawLatex(0.1,0.92,"CMS")
        cmsLatex.DrawLatex(0.15, 0.87, "CMS")
        cmsLatex.SetTextFont(52)
        # cmsLatex.DrawLatex(0.1+0.08,0.92,"Preliminary")
        cmsLatex.DrawLatex(0.15 + 0.071, 0.87, "Preliminary")

        cmsLatex.SetTextAlign(31)
        cmsLatex.SetTextFont(42)
        if args.year == '2016':
            lumiText = '16.81 fb^{-1}, 13 TeV'
        elif args.year == '2016APV':
            lumiText = '19.52 fb^{-1}, 13 TeV'
        elif args.year == '2016tot':
            lumiText = '36.33 fb^{-1}, 13 TeV'
        elif args.year == '2017':
            lumiText = '41.48 fb^{-1}, 13 TeV'
        elif args.year == '2018':
            lumiText = '59.83 fb^{-1}, 13 TeV'
        cmsLatex.DrawLatex(0.935, 0.87, lumiText)
        # cmsLatex.DrawLatex(0.85,0.87,lumiText)

    ##########################################################################

    ################################# Shape Check Canvas######################
        # ShapeCanvas = ROOT.TCanvas("ShapeCanvas","ShapeCanvas")
        # ShapeCanvas.Divide(1,2)
        ###
        # Shape_plotPad = ROOT.gPad.GetPrimitive('ShapeCanvas_1')
        # Shape_ratioPad = ROOT.gPad.GetPrimitive('ShapeCanvas_2')
        ###
        # Shape_plotPad.SetPad("Shape_pad1","Shape_plot",0.0,0.20,1.0,1.0,0)
        # Shape_ratioPad.SetPad("Shape_pad2","Shape_ratio",0.0,0.0,1.0,0.25,0)
####
        # Shape_ratioPad.SetTopMargin(0.05)
        # Shape_ratioPad.SetFrameLineWidth(1)
        # Shape_ratioPad.SetBottomMargin(0.27)
        # Shape_plotPad.SetBottomMargin(0.08)
        # Shape_ratioPad.SetGridy()
####
        # Shape_ratioHist, Shape_ratioError = MakeRatioHistograms(DataShape,ShapeStack,variable)
        # Shape_ratioPad.cd()
        # Shape_ratioHist.Draw('ex0')
        # Shape_ratioError.Draw('SAME e2')
        # Shape_ratioHist.Draw('SAME ex0')
####
        # Shape_plotPad.cd()
        # Shape_plotPad.SetFrameLineWidth(1)
        # Shape_plotPad.SetTickx()
        # Shape_plotPad.SetTicky()
####
        # ShapeStack.SetMaximum(max(ShapeStack.GetMaximum(),DataShape.GetMaximum()))
        ###
        # ShapeStack.Draw()
        # ShapeStack_Errors.Draw('SAME e2')
        # ShapeStack.SetTitle(variableAxisTitleDictionary[variable])
        # DataShape.Draw('SAME e1')
        # signalHisto.Draw('SAME HIST')
        # ShapeStack.GetYaxis().SetTitle("Events Normalized")
        # ShapeStack.GetYaxis().SetTitleOffset(1.58)
        # ShapeStack.GetXaxis().SetLabelSize(0.0)
###
        # theLegend2 = ROOT.TLegend(0.85, 0.45, 1.0, 0.75, "", "brNDC")
        # theLegend2.SetHeader("#mu-#tau_{h} Channel","C")
        # theLegend2.SetLineWidth(0)
        # theLegend2.SetLineStyle(1)
        # theLegend2.SetFillStyle(1001) #0
        # theLegend2.SetFillColor(0)
        # theLegend2.SetBorderSize(0)
        # theLegend2.SetTextFont(42)
        # theLegend2.AddEntry(DataShape,'Observed','pe')
        # theLegend2.AddEntry(N_DiBoson_Histo,'DiBoson','f')
        # theLegend2.AddEntry(N_TT_Histo,'TTbar','f')
        # theLegend2.AddEntry(N_WJets_Histo,'WJets','f')
        # theLegend2.AddEntry(N_QCD_Histo,'QCD','f')
        # theLegend2.AddEntry(N_ST_Histo,'ST_s_Channel','f')
        # theLegend2.AddEntry(N_DY_Histo,'Drell-Yan','f')
        # theLegend2.AddEntry(N_Other_Histo,'Others','f')
        # theLegend2.AddEntry(Signal_Histo,'Radion (#times 50)','l')
###
        # theLegend2.Draw('SAME')
        # cmsLatex = ROOT.TLatex()
        # cmsLatex.SetTextSize(0.06)
        # cmsLatex.SetNDC(True)
        # cmsLatex.SetTextFont(61)
        # cmsLatex.SetTextAlign(11)
        # cmsLatex.DrawLatex(0.1,0.92,"CMS")
        # cmsLatex.SetTextFont(52)
        # cmsLatex.DrawLatex(0.1+0.08,0.92,"Preliminary")

    ############################# Saving The Plots############################
        # theCanvas.SaveAs('QuickControlPlots/'+variable+'_'+args.year+'.png')
# if (args.Channel == "tt"):
# theCanvas.SaveAs('TTPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.png')
# theCanvas.SaveAs('TTPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.pdf')
# elif (args.Channel == "et"):
# theCanvas.SaveAs('ETPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.png')
# theCanvas.SaveAs('ETPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.pdf')
# elif (args.Channel == "mt"):
# theCanvas.SaveAs('MTPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.png')
# theCanvas.SaveAs('MTPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.pdf')
# elif (args.Channel == "lt"):
# theCanvas.SaveAs('LTPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.png')
# theCanvas.SaveAs('LTPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.pdf')
# elif (args.Channel == "all"):
# theCanvas.SaveAs('allChannelPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.png')
# theCanvas.SaveAs('allChannelPlots/'+args.Sub+'/'+variableFileNameDictionary[variable]+'_'+args.Channel+'_'+args.year+'.pdf')
# else:
# print ("Enter a valid channel")

        output_dirs = {
            "tt": 'TTPlots/',
            "et": 'ETPlots/',
            "mt": 'MTPlots/',
            "lt": 'LTPlots/',
            "all": 'allChannelPlots/',
        }

        if args.Channel in output_dirs:
            output_dir = output_dirs[args.Channel] + args.Sub + '/'

            # Ensure the output directory exists
            if not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir)
                except OSError as e:
                    print(("Error creating directory {}: {}".format(output_dir, e)))
                    raise

            theCanvas.SaveAs(
                output_dir +
                variableFileNameDictionary[variable] +
                '_' +
                args.Channel +
                '_' +
                args.year +
                '.png')
            theCanvas.SaveAs(
                output_dir +
                variableFileNameDictionary[variable] +
                '_' +
                args.Channel +
                '_' +
                args.year +
                '.pdf')
        else:
            print("Enter a valid channel")

        # theCanvas.SaveAs('TTPlots/CorrWht/'+variableAxisTitleDictionary[variable]+'_'+args.year+'.pdf')
        # theCanvas.SaveAs('QuickControlPlots/'+variable+'_'+args.year+'.root')
        # ShapeCanvas.SaveAs('QuickControlPlots/'+ "Normalized" +variable+'_'+args.year+'.png')
        # ShapeCanvas.SaveAs('TTPlots/CorrWht/'+ "Normalized" +variableAxisTitleDictionary[variable]+'_'+args.year+'.pdf')
        # ShapeCanvas.SaveAs('QuickControlPlots/'+ "Normalized" +variable+'_'+args.year+'.root')
        # PassFailCanvas.SaveAs('QuickControlPlots/'+"Trigg_PF"+'_'+args.year+'.png')
        # PassFailCanvas.SaveAs('QuickControlPlots/'+"Trigg_PF"+'_'+args.year+'.pdf')
        # PassFailCanvas.SaveAs('QuickControlPlots/'+"Trigg_PF"+'_'+args.year+'.root')
    ##########################################################################

        if args.pause:
            eval(input("Press Enter to Continue..."))
        # this causes issues if you don't get rid of the canvas
        # I suspect it to be something to do with modifiying histograms that
        # are alreayd referenced by the canvas in preparing the next one
        del theCanvas
        # del PassFailCanvas
        # del ShapeCanvas

        # Delete the Objects created to avoid memory leaks
        del DatasetObjects
    # file.close()
    if (args.storeshape):
        shape_file.Close()


if __name__ == '__main__':
    main()
