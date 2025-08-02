import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from TauPOG.TauIDSFs.TauIDSFTool import TauESTool


class TauEnergyScaleForHPSandBoosted(Module):
    def __init__(self, year, isData=False, tauID_wp='Loose', ele_wp='VVLoose'):
        print("ACTIVATE:  Tau Energy Scale Producer")
        if (year == '2016'):
            self.year = '2016_postVFP'
        elif (year == '2016APV'):
            self.year = '2016_preVFP'
        else:
            self.year = year
        self.isData = isData
        self.isMC = not self.isData
        self.tauID_wp = tauID_wp
        self.ele_wp = ele_wp
        # Initialize Tau Energy Scale (TES) tool
        self.tesTool = TauESTool(
            'UL' + self.year,
            'DeepTau2018v2p5VSjet',
            wp=tauID_wp,
            wp_vsele=ele_wp)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        # Add new branches for TES-corrected pt, mass and variations
        self.out.branch("Tau_pt_nom", "F", lenVar="nTau")
        self.out.branch("Tau_mass_nom", "F", lenVar="nTau")
        self.out.branch("Tau_pt_tesUp", "F", lenVar="nTau")
        self.out.branch("Tau_pt_tesDown", "F", lenVar="nTau")
        self.out.branch("Tau_mass_tesUp", "F", lenVar="nTau")
        self.out.branch("Tau_mass_tesDown", "F", lenVar="nTau")

        self.out.branch("boostedTau_pt_nom", "F", lenVar="nboostedTau")
        self.out.branch("boostedTau_mass_nom", "F", lenVar="nboostedTau")
        self.out.branch("boostedTau_pt_tesUp", "F", lenVar="nboostedTau")
        self.out.branch("boostedTau_pt_tesDown", "F", lenVar="nboostedTau")
        self.out.branch("boostedTau_mass_tesUp", "F", lenVar="nboostedTau")
        self.out.branch("boostedTau_mass_tesDown", "F", lenVar="nboostedTau")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        # Get the Tau collection
        taus = Collection(event, "Tau")

        # Initialize lists to hold corrected values
        pt_nom = []
        mass_nom = []
        pt_tesUp = []
        pt_tesDown = []
        mass_tesUp = []
        mass_tesDown = []

        # Loop over all taus in the event
        for tau in taus:
            pt = tau.pt
            mass = tau.mass
            dm = tau.decayMode
            # genmatch = tau.genPartFlav  # Gen matching

            # Apply TES only if the tau passes both the vsJet and vsElectron WPs
            # if ((tau.idDeepTau2018v2p5VSjet >= 4 and tau.idDeepTau2018v2p5VSe
            # >= 2) and (self.isMC)):
            if ((self.isMC)):
                # Get TES correction for central value, Up, and Down
                # tes = self.tesTool.getTES(pt, dm, genmatch)
                # tesUp = self.tesTool.getTES(pt, dm, genmatch, unc="Up")
                # tesDown = self.tesTool.getTES(pt, dm, genmatch, unc="Down")

                # Apply the TES to pt and mass
                # pt_nom.append(pt * tes)
                # mass_nom.append(mass * tes)
                # pt_tesUp.append(pt * tesUp)
                # pt_tesDown.append(pt * tesDown)
                # mass_tesUp.append(mass * tesUp)
                # mass_tesDown.append(mass * tesDown)

                # Apply the TES to pt and mass - 3 % up and down
                pt_nom.append(pt * 1.0)
                mass_nom.append(mass * 1.0)
                pt_tesUp.append(pt * 1.03)
                pt_tesDown.append(pt * 0.97)
                mass_tesUp.append(mass * 1.03)
                mass_tesDown.append(mass * 0.97)

            else:  # for data)
                # If the tau doesn't pass the WPs, keep the original pt and
                # mass
                pt_nom.append(pt)
                mass_nom.append(mass)
                pt_tesUp.append(pt)
                pt_tesDown.append(pt)
                mass_tesUp.append(mass)
                mass_tesDown.append(mass)

        # Fill the new branches with the corrected values
        self.out.fillBranch("Tau_pt_nom", pt_nom)
        self.out.fillBranch("Tau_mass_nom", mass_nom)
        self.out.fillBranch("Tau_pt_tesUp", pt_tesUp)
        self.out.fillBranch("Tau_pt_tesDown", pt_tesDown)
        self.out.fillBranch("Tau_mass_tesUp", mass_tesUp)
        self.out.fillBranch("Tau_mass_tesDown", mass_tesDown)

        # Get the Tau collection
        boostedtaus = Collection(event, "boostedTau")

        # Initialize lists to hold corrected values
        pt_nom = []
        mass_nom = []
        pt_tesUp = []
        pt_tesDown = []
        mass_tesUp = []
        mass_tesDown = []

        # Loop over all taus in the event
        for boostedtau in boostedtaus:
            pt = boostedtau.pt
            mass = boostedtau.mass
            dm = boostedtau.decayMode
            # genmatch = boostedtau.genPartFlav  # Gen matching

            # Apply TES only if the tau passes both the vsJet and vsElectron WPs
            # if ((tau.idDeepTau2018v2p5VSjet >= 4 and tau.idDeepTau2018v2p5VSe
            # >= 2) and (self.isMC)):
            if ((self.isMC)):
                # Get TES correction for central value, Up, and Down
                # tes = self.tesTool.getTES(pt, dm, genmatch)
                # tesUp = self.tesTool.getTES(pt, dm, genmatch, unc="Up")
                # tesDown = self.tesTool.getTES(pt, dm, genmatch, unc="Down")

                # Apply the TES to pt and mass
                # pt_nom.append(pt * tes)
                # mass_nom.append(mass * tes)
                # pt_tesUp.append(pt * tesUp)
                # pt_tesDown.append(pt * tesDown)
                # mass_tesUp.append(mass * tesUp)
                # mass_tesDown.append(mass * tesDown)

                # Apply the TES to pt and mass - 3 % up and down
                pt_nom.append(pt * 1.0)
                mass_nom.append(mass * 1.0)
                pt_tesUp.append(pt * 1.03)
                pt_tesDown.append(pt * 0.97)
                mass_tesUp.append(mass * 1.03)
                mass_tesDown.append(mass * 0.97)

            else:  # (for data)
                # If the tau doesn't pass the WPs, keep the original pt and
                # mass
                pt_nom.append(pt)
                mass_nom.append(mass)
                pt_tesUp.append(pt)
                pt_tesDown.append(pt)
                mass_tesUp.append(mass)
                mass_tesDown.append(mass)

        # Fill the new branches with the corrected values
        self.out.fillBranch("boostedTau_pt_nom", pt_nom)
        self.out.fillBranch("boostedTau_mass_nom", mass_nom)
        self.out.fillBranch("boostedTau_pt_tesUp", pt_tesUp)
        self.out.fillBranch("boostedTau_pt_tesDown", pt_tesDown)
        self.out.fillBranch("boostedTau_mass_tesUp", mass_tesUp)
        self.out.fillBranch("boostedTau_mass_tesDown", mass_tesDown)

        return True
