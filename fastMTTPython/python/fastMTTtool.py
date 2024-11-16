#Andrew Loeliger
#Quick tool to make the calling of the fast mtt process a bit easier to handle
#from the python side, so that a long series of incomprehensible numbers don't have to be inserted
#into a less than helpful function descripton

import pluginfastMTT_binding as fastMTT
import re
import ROOT

#There are three things that need to go into the fast mtt algorithm

#The first is met information
# this includes x and y components
# and then the covariance matrix

#the second is the first tau information
#Taus can be decayed into muons, electrons, or hadrons
#This needs to be made explicit
#They all have a pt eta phi m, p4 vector that gets input
#and in the case of the hadronic tau, they have a decay mode

class fastMTTlepton():
    def __init__(self, pt = 0.0, eta = 0.0, phi = 0.0, m = 0.0, leptonType = '', tauDecayMode=0):
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.m = m
        
        if re.match("[eE](lectron)?", leptonType):
            self.leptonType = 3
        elif re.match("[mM](uon|u)?", leptonType):
            self.leptonType = 2
        elif re.match("([hH]ad(ronic)?)?[tT]au", leptonType):
            self.leptonType = 1
        else:
            self.leptonType = 0

        self.tauDecayMode = tauDecayMode
        
    def getPt(self):
        return self.pt
    def getEta(self):
        return self.eta
    def getPhi(self):
        return self.phi
    def getM(self):
        return self.m
    def getLeptonType(self):
        return self.leptonType
    def getTauDecayMode(self):
        return self.tauDecayMode

    def setPt(self, pt):
        self.pt = pt
    def setEta(self, eta):
        self.eta = eta
    def setPhi(self, phi):
        self.phi = phi
    def setM(self, m):
        self.m = m
    def setLeptonType(self, leptonType):
        self.leptonType = leptonType
    def setTauDecayMode(self, tauDecayMode):
        self.tauDecayMode = tauDecayMode

class fastMTTmet():
    def __init__(self, measuredX = 0.0, measuredY = 0.0, xx = 0.0, xy = 0.0, yy = 0.0):
        self.measuredX = measuredX
        self.measuredY = measuredY
        self.xx = xx
        self.xy = xy
        self.yy = yy
        
    def getMeasuredX(self):
        return self.measuredX
    def getMeasuredY(self):
        return self.measuredY
    def getXXTerm(self):
        return self.xx
    def getXYTerm(self):
        return self.xy
    def getYYTerm(self):
        return self.yy

    def setMeasuredX(self):
        return self.measuredX
    def setMeasuredY(self):
        return self.measuredY
    def setXXTerm(self, xx):
        self.xx = xx
    def setXYTerm(self, xy):
        self.xy = xy
    def setYYTerm(self, yy):
        self.yy = yy

class fastMTTtool():
    def __init__(self, firstLepton = fastMTTlepton, secondLepton = fastMTTlepton(), theMET = fastMTTmet()):
        self.firstLepton = firstLepton
        self.secondLepton = secondLepton
        self.theMET = theMET

    def getFirstLepton(self):
        return self.firstLepton
    def getSecondLepton(self):
        return self.secondLepton
    def getTheMET(self):
        return self.theMET
        
    def setFirstLepton(self, firstLepton):
        self.firstLepton = firstLepton
    def setSecondLepton(self, secondLepton):
        self.secondLepton = secondLepton
    def setTheMET(self, theMET):
        self.theMET = theMET

    def getFastMTTmass(self):
        return fastMTT.fastMTTmass(
            self.theMET.getMeasuredX(),
            self.theMET.getMeasuredY(),
            self.theMET.getXXTerm(),
            self.theMET.getXYTerm(),
            self.theMET.getYYTerm(),
            self.firstLepton.getLeptonType(),
            self.firstLepton.getPt(),
            self.firstLepton.getEta(),
            self.firstLepton.getPhi(),
            self.firstLepton.getM(),
            self.firstLepton.getTauDecayMode(),
            self.secondLepton.getLeptonType(),
            self.secondLepton.getPt(),
            self.secondLepton.getEta(),
            self.secondLepton.getPhi(),
            self.secondLepton.getM(),
            self.secondLepton.getTauDecayMode()
        )
    def getFastMTTpt(self):
        return fastMTT.fastMTTpt(
            self.theMET.getMeasuredX(),
            self.theMET.getMeasuredY(),
            self.theMET.getXXTerm(),
            self.theMET.getXYTerm(),
            self.theMET.getYYTerm(),
            self.firstLepton.getLeptonType(),
            self.firstLepton.getPt(),
            self.firstLepton.getEta(),
            self.firstLepton.getPhi(),
            self.firstLepton.getM(),
            self.firstLepton.getTauDecayMode(),
            self.secondLepton.getLeptonType(),
            self.secondLepton.getPt(),
            self.secondLepton.getEta(),
            self.secondLepton.getPhi(),
            self.secondLepton.getM(),
            self.secondLepton.getTauDecayMode()
        )

    def getFastMTTphi(self):
        return fastMTT.fastMTTphi(
            self.theMET.getMeasuredX(),
            self.theMET.getMeasuredY(),
            self.theMET.getXXTerm(),
            self.theMET.getXYTerm(),
            self.theMET.getYYTerm(),
            self.firstLepton.getLeptonType(),
            self.firstLepton.getPt(),
            self.firstLepton.getEta(),
            self.firstLepton.getPhi(),
            self.firstLepton.getM(),
            self.firstLepton.getTauDecayMode(),
            self.secondLepton.getLeptonType(),
            self.secondLepton.getPt(),
            self.secondLepton.getEta(),
            self.secondLepton.getPhi(),
            self.secondLepton.getM(),
            self.secondLepton.getTauDecayMode()
        )
    
    def getFastMTTeta(self):
        return fastMTT.fastMTTeta(
            self.theMET.getMeasuredX(),
            self.theMET.getMeasuredY(),
            self.theMET.getXXTerm(),
            self.theMET.getXYTerm(),
            self.theMET.getYYTerm(),
            self.firstLepton.getLeptonType(),
            self.firstLepton.getPt(),
            self.firstLepton.getEta(),
            self.firstLepton.getPhi(),
            self.firstLepton.getM(),
            self.firstLepton.getTauDecayMode(),
            self.secondLepton.getLeptonType(),
            self.secondLepton.getPt(),
            self.secondLepton.getEta(),
            self.secondLepton.getPhi(),
            self.secondLepton.getM(),
            self.secondLepton.getTauDecayMode()
        )

    def getFastMTTfourvector(self):
        return fastMTT.fastMTTfourvector(
            self.theMET.getMeasuredX(),
            self.theMET.getMeasuredY(),
            self.theMET.getXXTerm(),
            self.theMET.getXYTerm(),
            self.theMET.getYYTerm(),
            self.firstLepton.getLeptonType(),
            self.firstLepton.getPt(),
            self.firstLepton.getEta(),
            self.firstLepton.getPhi(),
            self.firstLepton.getM(),
            self.firstLepton.getTauDecayMode(),
            self.secondLepton.getLeptonType(),
            self.secondLepton.getPt(),
            self.secondLepton.getEta(),
            self.secondLepton.getPhi(),
            self.secondLepton.getM(),
            self.secondLepton.getTauDecayMode()
        )

        