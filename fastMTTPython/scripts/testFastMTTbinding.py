#Andrew Loeliger
#quick script to check that the python version of the fast MTT binding is yielding consistent
#answers to the cpp version of the testing script
import pluginfastMTT_binding as fastMTT

print (fastMTT.fastMTTmass(
    11.7491,
    -51.9172,
    787.352,
    -178.63,
    179.545,
    3,
    33.7383,
    0.9409,
    -0.541458,
    0.511E-3,
    0,
    1,
    25.7322,
    0.618228,
    2.79362,
    0.13957,
    0)
)

print (fastMTT.fastMTTpt(
    11.7491,
    -51.9172,
    787.352,
    -178.63,
    179.545,
    3,
    33.7383,
    0.9409,
    -0.541458,
    0.511E-3,
    0,
    1,
    25.7322,
    0.618228,
    2.79362,
    0.13957,
    0)
)
