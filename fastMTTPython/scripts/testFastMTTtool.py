#Andrew Loeliger
#quick script for testing the operation of the (hopefully) simplifiying fast mtt tool

from FinalStateHHbbtt.fastMTTPython.fastMTTtool import *
import ROOT
import time

exampleMET = fastMTTmet(
    measuredX = 11.7491,
    measuredY = -51.9172,
    xx = 787.352,
    xy = -178.63,
    yy = 179.545,
)

firstLepton = fastMTTlepton(
    pt = 33.7383,
    eta = 0.9409,
    phi = -0.541458,
    m = 0.511E-3,
    leptonType = 'e'
)

secondLepton = fastMTTlepton(
    pt = 25.7322,
    eta = 0.618228,
    phi = 2.79362,
    m = 0.13957,
    leptonType = 'Tau'
)

theFastMTTtool = fastMTTtool()
theFastMTTtool.setFirstLepton(firstLepton)
theFastMTTtool.setSecondLepton(secondLepton)
theFastMTTtool.setTheMET(exampleMET)

indi_st = time.time()
print theFastMTTtool.getFastMTTmass()
print theFastMTTtool.getFastMTTpt()
print theFastMTTtool.getFastMTTeta()
print theFastMTTtool.getFastMTTphi()
indi_et = time.time()
print ("Time taken for individual retrival = ",(indi_et - indi_st))

full_st = time.time()
print theFastMTTtool.getFastMTTfourvector()
full_et = time.time()

print ("Time taken for Full retrival = ",(full_et - full_st))


