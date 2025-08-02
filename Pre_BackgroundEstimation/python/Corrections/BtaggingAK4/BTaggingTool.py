# https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/BTagCalibration
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
from .helper import modulepath, ensureTFile, warning
from array import array
import ROOT
# ROOT.gROOT.ProcessLine('.L ./BTagCalibrationStandalone.cpp+')
from ROOT import TH2F, BTagCalibration, BTagCalibrationReader
from ROOT.BTagEntry import OP_LOOSE, OP_MEDIUM, OP_TIGHT, OP_RESHAPING
from ROOT.BTagEntry import FLAV_B, FLAV_C, FLAV_UDSG
path = modulepath


class BTagWPs:
    """Contain b tagging working points."""

    def __init__(self, tagger, year=2016):
        assert (year in ["2016", "2016APV", "2017", "2018"]
                ), "You must choose a year from: 2016, 2016APV, 2017, or 2018."
        if year == "2016":
            if 'deepjetflavb' in tagger.lower():
                self.loose = 0.0480
                self.medium = 0.2489
                self.tight = 0.6377

        elif year == "2016APV":
            if 'deepjetflavb' in tagger.lower():
                self.loose = 0.0508
                self.medium = 0.2598
                self.tight = 0.6502

        elif year == "2017":
            if 'deepjetflavb' in tagger.lower():
                self.loose = 0.0532
                self.medium = 0.3040
                self.tight = 0.7476

        elif year == "2018":
            if 'deepjetflavb' in tagger.lower():
                self.loose = 0.0490
                self.medium = 0.2783
                self.tight = 0.7100


class BTagWeightTool:

    def __init__(
            self,
            tagger,
            wp,
            sigmabc='central',
            sigmalight='central',
            year=2016):
        """Load b tag weights from CSV file."""

        assert (year in ["2016", "2016APV", "2017", "2018"]
                ), "You must choose a year from: 2016, 2016APV, 2017, or 2018."
        assert (tagger in [
                'deepjetflavb']), "BTagWeightTool: You must choose a tagger from: deepjetflavb!"
        assert (wp in ['loose', 'medium', 'tight']
                ), "BTagWeightTool: You must choose a WP from: loose, medium, tight!"
        # assert(sigmabc in ['central','up_correlated','down_correlated','up_uncorrelated','down_uncorrelated']), "BTagWeightTool: You must choose a WP for b/c jets from: central, up_correlated, down_correlated, up_uncorrelated, down_uncorrelated!"
        # assert(sigmalight in ['central','up_correlated','down_correlated','up_uncorrelated','down_uncorrelated']), "BTagWeightTool: You must choose a WP for light jets from: central, up_correlated, down_correlated, up_uncorrelated, down_uncorrelated!"
        # assert(channel in ['mutau','eletau','tautau','mumu']), "BTagWeightTool: You must choose a channel from: mutau, eletau, tautau, mumu!"

        # FILE
        if year == "2016":
            if 'deepjetflavb' in tagger.lower():
                # csvnamebc = path+'DeepCSV_2016LegacySF_V1_YearCorrelation-V1.csv'
                csvnamebc = path + '/CSVs/wp_deepJet_2016_tes2.csv'
                csvnamelight = path + 'CSVs/wp_deepJet_2016_test2.csv'
                effname = path + '/' + year + '/DeepJetFlavB_2016_eff_' + wp.upper() + '.root'

        if year == "2016APV":
            if 'deepjetflavb' in tagger.lower():
                # csvnamebc = path+'DeepCSV_2016LegacySF_V1_YearCorrelation-V1.csv'
                csvnamebc = path + '/CSVs/wp_deepJet_2016APV.csv'
                csvnamelight = path + 'CSVs/wp_deepJet_2016APV.csv'
                effname = path + '/' + year + '/DeepJetFlavB_2016APV_eff_' + wp.upper() + \
                    '.root'

        elif year == "2017":
            if 'deepjetflavb' in tagger.lower():
                # csvnamebc = path+'DeepCSV_94XSF_V4_B_F_YearCorrelation-V1.csv'
                csvnamebc = path + '/CSVs/wp_deepJet_2017.csv'
                csvnamelight = path + '/CSVs/wp_deepJet_2017.csv'
                effname = path + '/' + year + '/DeepJetFlavB_2017_eff_' + wp.upper() + '.root'

        elif year == "2018":
            if 'deepjetflavb' in tagger.lower():
                # csvnamebc = path+'DeepCSV_102XSF_V1_YearCorrelation-V1.csv'
                csvnamebc = path + '/CSVs/wp_deepJet_2018.csv'
                csvnamelight = path + '/CSVs/wp_deepJet_2018.csv'
                effname = path + '/' + year + '/DeepJetFlavB_2018_eff_' + wp.upper() + '.root'

        # TAGGING WP
        self.wpname = wp
        self.wp = getattr(BTagWPs(tagger, year), wp)
        if 'deepjetflavb' in tagger.lower():
            def tagged(j): return j.btagDeepFlavB > self.wp

        # CSV READER
        print(("Loading BTagWeightTool for %s (%s WP)..." % (tagger, wp)))
        print((FLAV_B, type(FLAV_B)))
        op = OP_LOOSE if wp == 'loose' else OP_MEDIUM if wp == 'medium' else OP_TIGHT if wp == 'tight' else OP_RESHAPING
        type_udsg = 'incl'
        type_bc = 'mujets'  # 'mujets' for QCD; 'comb' for QCD+TT
        # Load reader for b/c jets
        print(("tagger string : ", tagger, type(tagger),
              " csvnamebc :", csvnamebc, type(csvnamebc)))
        calibbc = BTagCalibration(tagger, csvnamebc)
        readerbc = BTagCalibrationReader(op, sigmabc)
        readerbc.load(calibbc, FLAV_B, type_bc)
        readerbc.load(calibbc, FLAV_C, type_bc)
        # Load reader for light jets
        caliblight = BTagCalibration(tagger, csvnamelight)
        readerlight = BTagCalibrationReader(op, sigmalight)
        readerlight.load(caliblight, FLAV_UDSG, type_udsg)

        # EFFICIENCIES
        effmaps = {}  # b tag efficiencies in MC to compute b tagging weight for an event
        efffile = ensureTFile(effname)
        default = False
        if not efffile:
            warning(
                "File %s with efficiency histograms does not exist! Reverting to default efficiency histogram..." %
                (effname), title="BTagWeightTool")
            default = True
        for flavor in [0, 4, 5]:
            flavor = flavorToString(flavor)
            # histname = "%s_%s_%s"%(tagger,flavor,wp)
            effname = "%s/eff_%s_%s_%s" % (year, tagger, flavor, wp)
            if efffile:
                effmaps[flavor] = efffile.Get(effname)
                if not effmaps[flavor]:
                    warning(
                        "histogram '%s' does not exist in %s! Reverting to default efficiency histogram..." %
                        (effname, efffile.GetName()), title="BTagWeightTool")
                    default = True
                    effmaps[flavor] = createDefaultEfficiencyMap(
                        effname, flavor, wp)
            else:
                effmaps[flavor] = createDefaultEfficiencyMap(
                    effname, flavor, wp)
            effmaps[flavor].SetDirectory(0)
        efffile.Close()

        if default:
            warning(
                "Made use of default efficiency histograms! The b tag weights from this module should be regarded as placeholders only,\n" +
                "and should NOT be used for analyses. B (mis)tag efficiencies in MC are analysis dependent. Please create your own\n" +
                "efficiency histogram with corrections/btag/getBTagEfficiencies.py after running all MC samples with BTagWeightTool.",
                title="BTagWeightTool")

        self.tagged = tagged
        self.calibbc = calibbc
        self.readerbc = readerbc
        self.caliblight = caliblight
        self.readerlight = readerlight
        self.effmaps = effmaps

    def getWeight(self, jetCollection):
        weight = 1.
        for jet in jetCollection:
            weight *= self.getSF(jet)
        return weight

    def getSF(self, jet):
        """Get b tag SF for a single jet."""
        pt = jet.pt
        eta = jet.eta
        flavor = jet.hadronFlavour
        FLAV = flavorToFLAV(flavor)
        if FLAV == FLAV_UDSG:
            SF = self.readerlight.eval(FLAV, abs(eta), pt)
        else:
            SF = self.readerbc.eval(FLAV, abs(eta), pt)
        tagged = self.tagged(jet)
        if tagged:
            weight = SF
        else:
            eff = self.getEfficiency(pt, eta, flavor)
            if eff == 1:
                print(("Warning! BTagWeightTool.getSF: MC efficiency is 1 for pt=%s, eta=%s, flavor=%s, SF=%s" % (
                    pt, eta, flavor, SF)))
                return 1
            else:
                weight = (1 - SF * eff) / (1 - eff)
        return weight

    def getEfficiency(self, pt, eta, flavor):
        """Get b tag efficiency for a single jet in MC."""
        flavor = flavorToString(flavor)
        hist = self.effmaps[flavor]
        xbin = hist.GetXaxis().FindBin(pt)
        ybin = hist.GetYaxis().FindBin(eta)
        if xbin == 0:
            xbin = 1
        elif xbin > hist.GetXaxis().GetNbins():
            xbin -= 1
        if ybin == 0:
            ybin = 1
        elif ybin > hist.GetYaxis().GetNbins():
            ybin -= 1
        eff = hist.GetBinContent(xbin, ybin)
        # if eff==1:
        # print "Warning! BTagWeightTool.getEfficiency: MC efficiency is 1 for
        # pt=%s, eta=%s, flavor=%s, SF=%s"%(pt,eta,flavor,SF)
        return eff


def flavorToFLAV(flavor):
    """Help function to convert an integer flavor ID to a BTagEntry enum value."""
# return FLAV_B if abs(flavor)==5 else FLAV_C if abs(flavor) in [4,15]
# else FLAV_UDSG
    return FLAV_B if abs(flavor) == 5 else FLAV_C if abs(
        flavor) == 4 else FLAV_UDSG


def flavorToString(flavor):
    """Help function to convert an integer flavor ID to a string value."""
    return 'b' if abs(flavor) == 5 else 'c' if abs(flavor) == 4 else 'udsg'


def createEfficiencyMap(histname):
    """Help function to create efficiency maps (TH2D) with uniform binning and layout.
    One method to rule them all."""
#    ptbins  = array('d',[10,20,30,50,70,100,140,200,300,600,1000,1500])
#    etabins = array('d',[-2.5,-1.5,0.0,1.5,2.5])
    ptbins = array('d', [30, 40, 50, 70, 100, 150, 200, 300, 400, 1500])
    etabins = array('d', [-2.5, -1.5, 0.0, 1.5, 2.5])

    bins = (len(ptbins) - 1, ptbins, len(etabins) - 1, etabins)
    hist = TH2F(histname, histname, *bins)
    hist.GetXaxis().SetTitle("jet p_{T} [GeV]")
    hist.GetYaxis().SetTitle("jet #eta")
    hist.SetDirectory(0)
    return hist


def createDefaultEfficiencyMap(histname, flavor, wp='medium'):
    """Create default efficiency histograms. WARNING! Do not use for analysis! Use it as a placeholder,
    until you have made an efficiency map from MC for you analysis."""
    if wp == 'loose':
        eff = 0.75 if flavor == 'b' else 0.11 if flavor == 'c' else 0.01
    elif wp == 'medium':
        eff = 0.85 if flavor == 'b' else 0.42 if flavor == 'c' else 0.10
    else:
        eff = 0.60 if flavor == 'b' else 0.05 if flavor == 'c' else 0.001
    histname = histname.split('/')[-1] + "_default"
    hist = createEfficiencyMap(histname)
    for xbin in range(0, hist.GetXaxis().GetNbins() + 2):
        for ybin in range(0, hist.GetYaxis().GetNbins() + 2):
            hist.SetBinContent(xbin, ybin, eff)
    return hist
