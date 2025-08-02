from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import glob
import argparse
from .BTaggingTool import *
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True


class getBTagHist(Module):
    def __init__(self, tagger, wp, year, channel):
        assert (year in ['2016APV', '2016', '2017', '2018']
                ), "You must choose a year from: 2016, 2017, or 2018."
        assert (tagger in [
                'DeepJetFlavB']), "BTagWeightTool: You must choose a tagger from: CSVv2, DeepCSV!"
        assert (wp in ['loose', 'medium', 'tight']
                ), "BTagWeightTool: You must choose a WP from: loose, medium, tight!"

        self.writeHistFile = True
        self.tagger = tagger
        self.wp = wp
        self.year = year
        self.channel = channel

        threshold = getattr(BTagWPs(tagger, year), wp)
        if 'deepjetflavb' in tagger.lower():
            self.tagged = lambda j: j.btagDeepFlavB > threshold

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.hist_b = createEfficiencyMap(
            '%s_%s_%s' %
            (self.tagger, 'b', self.wp))
        self.hist_c = createEfficiencyMap(
            '%s_%s_%s' %
            (self.tagger, 'c', self.wp))
        self.hist_udsg = createEfficiencyMap(
            '%s_%s_%s' %
            (self.tagger, 'udsg', self.wp))
        self.hist_b_all = createEfficiencyMap(
            '%s_%s_%s_all' %
            (self.tagger, 'b', self.wp))
        self.hist_c_all = createEfficiencyMap(
            '%s_%s_%s_all' %
            (self.tagger, 'c', self.wp))
        self.hist_udsg_all = createEfficiencyMap(
            '%s_%s_%s_all' %
            (self.tagger, 'udsg', self.wp))

        self.addObject(self.hist_b)
        self.addObject(self.hist_c)
        self.addObject(self.hist_udsg)
        self.addObject(self.hist_b_all)
        self.addObject(self.hist_c_all)
        self.addObject(self.hist_udsg_all)

    def analyze(self, event):
        # process event, return True (go to next module) or False (fail, go to
        # next event)
        Jet = Collection(event, 'Jet', 'nJet')

        for index in range(event.ngood_Jets):
            # flavor = flavorToString(jet.hadronFlavour)
            flavor = flavorToString(
                Jet[event.index_gJets[index]].hadronFlavour)
            if flavor == 'b':
                self.hist_b_all.Fill(Jet[event.index_gJets[index]].pt,
                                     Jet[event.index_gJets[index]].eta, event.FinalWeighting)
                if self.tagged(Jet[event.index_gJets[index]]):
                    self.hist_b.Fill(Jet[event.index_gJets[index]].pt,
                                     Jet[event.index_gJets[index]].eta,
                                     event.FinalWeighting)
            elif flavor == 'c':
                self.hist_c_all.Fill(Jet[event.index_gJets[index]].pt,
                                     Jet[event.index_gJets[index]].eta, event.FinalWeighting)
                if self.tagged(Jet[event.index_gJets[index]]):
                    self.hist_c.Fill(Jet[event.index_gJets[index]].pt,
                                     Jet[event.index_gJets[index]].eta,
                                     event.FinalWeighting)
            elif flavor == 'udsg':
                self.hist_udsg_all.Fill(
                    Jet[event.index_gJets[index]].pt, Jet[event.index_gJets[index]].eta, event.FinalWeighting)
                if self.tagged(Jet[event.index_gJets[index]]):
                    self.hist_udsg.Fill(
                        Jet[event.index_gJets[index]].pt, Jet[event.index_gJets[index]].eta, event.FinalWeighting)

        return True


parser = argparse.ArgumentParser(description='')
parser.add_argument(
    '--mcPath',
    '-p',
    help="enter the path to mc files",
    default="")
parser.add_argument(
    '--year',
    '-y',
    help="Enter the year of data taking",
    choices=[
        "2016",
        "2016APV",
        "2017",
        "2018"],
    required=True)
parser.add_argument(
    '--wp',
    help="Enter the wp",
    choices=[
        "loose",
        "medium",
        "tight"],
    required=True)
parser.add_argument(
    '--allwp',
    help='mnake hists for all three wps',
    action='store_true')
args = parser.parse_args()

wplist = ["loose", "medium", "tight"]


fnames_mc = list(set(glob.glob(args.mcPath +
                               "/*.root")) -
                 set(glob.glob(args.mcPath +
                               "/Radion*.root")) -
                 set(glob.glob(args.mcPath +
                               "/*MET*Run*.root")))

print("These are the files that are being considered")
print([file.split('/')[-1] for file in fnames_mc])

if (args.allwp):
    for wp in wplist:
        def getBTagHistLam(): return getBTagHist(
            'DeepJetFlavB', wp, args.year, 'ttbar')
        print(("Creating/Running BtagHist processor for " +
              args.year + " MC : " + wp + " wp"))
        post_MC = PostProcessor(
            args.year +
            "/.",
            fnames_mc,
            cut=None,
            branchsel=None,
            modules=[
                getBTagHistLam()],
            noOut=True,
            histFileName=args.year +
            "/btagHists_" +
            args.year +
            "_" +
            wp.upper() +
            ".root",
            histDirName=args.year)
        post_MC.run()
else:
    print(("Creating/Running BtagHist processor for " +
          args.year + " MC : " + args.wp + " wp"))

    def getBTagHistLam(): return getBTagHist(
        'DeepJetFlavB', args.wp, args.year, 'ttbar')
    post_MC = PostProcessor(
        args.year +
        "/.",
        fnames_mc,
        cut=None,
        branchsel=None,
        modules=[
            getBTagHistLam()],
        noOut=True,
        histFileName=args.year +
        "/btagHists_" +
        args.year +
        "_" +
        args.wp.upper() +
        ".root",
        histDirName=args.year)
    post_MC.run()
