def get2016EWKW(pt):
    weight = 1
    if pt < 170:
        weight = 0.94269
    elif 170 <= pt < 200:
        weight = 0.902615
    elif 200 <= pt < 230:
        weight = 0.898827
    elif 230 <= pt < 260:
        weight = 0.959081
    elif 260 <= pt < 290:
        weight = 0.891248
    elif 290 <= pt < 320:
        weight = 0.860188
    elif 320 <= pt < 350:
        weight = 0.884811
    elif 350 <= pt < 390:
        weight = 0.868131
    elif 390 <= pt < 430:
        weight = 0.848655
    elif 430 <= pt < 470:
        weight = 0.806186
    elif 470 <= pt < 510:
        weight = 0.848507
    elif 510 <= pt < 550:
        weight = 0.83763
    elif 550 <= pt < 590:
        weight = 0.792152
    elif 590 <= pt < 640:
        weight = 0.730731
    elif 640 <= pt < 690:
        weight = 0.778061
    elif 690 <= pt < 740:
        weight = 0.771811
    elif 740 <= pt < 790:
        weight = 0.795004
    elif 790 <= pt < 840:
        weight = 0.757859
    elif 840 <= pt < 900:
        weight = 0.709571
    elif 900 <= pt < 960:
        weight = 0.702751
    elif 960 <= pt < 1020:
        weight = 0.657821
    elif 1020 <= pt < 1090:
        weight = 0.762559
    elif 1090 <= pt < 1160:
        weight = 0.845925
    elif pt >= 1160:
        weight = 0.674034

    return weight


def get2016QCDW(pt):
    weight = 1
    if pt < 170:
        weight = 1.43896
    elif 170 <= pt < 200:
        weight = 1.45307
    elif 200 <= pt < 230:
        weight = 1.41551
    elif 230 <= pt < 260:
        weight = 1.42199
    elif 260 <= pt < 290:
        weight = 1.3477
    elif 290 <= pt < 320:
        weight = 1.35302
    elif 320 <= pt < 350:
        weight = 1.34289
    elif 350 <= pt < 390:
        weight = 1.32474
    elif 390 <= pt < 430:
        weight = 1.23267
    elif 430 <= pt < 470:
        weight = 1.22641
    elif 470 <= pt < 510:
        weight = 1.23149
    elif 510 <= pt < 550:
        weight = 1.21593
    elif 550 <= pt < 590:
        weight = 1.16506
    elif 590 <= pt < 640:
        weight = 1.01718
    elif 640 <= pt < 690:
        weight = 1.01575
    elif 690 <= pt < 740:
        weight = 1.05425
    elif 740 <= pt < 790:
        weight = 1.05992
    elif 790 <= pt < 840:
        weight = 1.01503
    elif 840 <= pt < 900:
        weight = 1.01761
    elif 900 <= pt < 960:
        weight = 0.947194
    elif 960 <= pt < 1020:
        weight = 0.932754
    elif 1020 <= pt < 1090:
        weight = 1.00849
    elif 1090 <= pt < 1160:
        weight = 0.94805
    elif pt >= 1160:
        weight = 0.86956

    return weight


def get2016RenUpW(pt):
    weight = 1
    if pt < 170:
        weight = 1.33939
    elif 170 <= pt < 200:
        weight = 1.34649
    elif 200 <= pt < 230:
        weight = 1.3115
    elif 230 <= pt < 260:
        weight = 1.30782
    elif 260 <= pt < 290:
        weight = 1.24619
    elif 290 <= pt < 320:
        weight = 1.24531
    elif 320 <= pt < 350:
        weight = 1.23472
    elif 350 <= pt < 390:
        weight = 1.20709
    elif 390 <= pt < 430:
        weight = 1.13368
    elif 430 <= pt < 470:
        weight = 1.12227
    elif 470 <= pt < 510:
        weight = 1.12612
    elif 510 <= pt < 550:
        weight = 1.1118
    elif 550 <= pt < 590:
        weight = 1.06549
    elif 590 <= pt < 640:
        weight = 0.931838
    elif 640 <= pt < 690:
        weight = 0.929282
    elif 690 <= pt < 740:
        weight = 0.959553
    elif 740 <= pt < 790:
        weight = 0.955823
    elif 790 <= pt < 840:
        weight = 0.920614
    elif 840 <= pt < 900:
        weight = 0.917243
    elif 900 <= pt < 960:
        weight = 0.855649
    elif 960 <= pt < 1020:
        weight = 0.84587
    elif 1020 <= pt < 1090:
        weight = 0.906862
    elif 1090 <= pt < 1160:
        weight = 0.858763
    elif pt >= 1160:
        weight = 0.794909

    return weight


def get2016RenDownW(pt):
    weight = 1
    if pt < 170:
        weight = 1.52924
    elif 170 <= pt < 200:
        weight = 1.55189
    elif 200 <= pt < 230:
        weight = 1.5104
    elif 230 <= pt < 260:
        weight = 1.53402
    elif 260 <= pt < 290:
        weight = 1.44328
    elif 290 <= pt < 320:
        weight = 1.45859
    elif 320 <= pt < 350:
        weight = 1.4486
    elif 350 <= pt < 390:
        weight = 1.44818
    elif 390 <= pt < 430:
        weight = 1.33058
    elif 430 <= pt < 470:
        weight = 1.33239
    elif 470 <= pt < 510:
        weight = 1.33976
    elif 510 <= pt < 550:
        weight = 1.32336
    elif 550 <= pt < 590:
        weight = 1.26854
    elif 590 <= pt < 640:
        weight = 1.10489
    elif 640 <= pt < 690:
        weight = 1.1055
    elif 690 <= pt < 740:
        weight = 1.15559
    elif 740 <= pt < 790:
        weight = 1.17682
    elif 790 <= pt < 840:
        weight = 1.118
    elif 840 <= pt < 900:
        weight = 1.13097
    elif 900 <= pt < 960:
        weight = 1.04988
    elif 960 <= pt < 1020:
        weight = 1.02796
    elif 1020 <= pt < 1090:
        weight = 1.12438
    elif 1090 <= pt < 1160:
        weight = 1.04704
    elif pt >= 1160:
        weight = 0.94791

    return weight


def get2016FacUpW(pt):
    weight = 1
    if pt < 170:
        weight = 1.44488
    elif 170 <= pt < 200:
        weight = 1.45847
    elif 200 <= pt < 230:
        weight = 1.41004
    elif 230 <= pt < 260:
        weight = 1.39864
    elif 260 <= pt < 290:
        weight = 1.3391
    elif 290 <= pt < 320:
        weight = 1.33663
    elif 320 <= pt < 350:
        weight = 1.32073
    elif 350 <= pt < 390:
        weight = 1.3076
    elif 390 <= pt < 430:
        weight = 1.20904
    elif 430 <= pt < 470:
        weight = 1.20066
    elif 470 <= pt < 510:
        weight = 1.20462
    elif 510 <= pt < 550:
        weight = 1.18286
    elif 550 <= pt < 590:
        weight = 1.12586
    elif 590 <= pt < 640:
        weight = 0.990615
    elif 640 <= pt < 690:
        weight = 0.984473
    elif 690 <= pt < 740:
        weight = 1.0171
    elif 740 <= pt < 790:
        weight = 1.01706
    elif 790 <= pt < 840:
        weight = 0.986107
    elif 840 <= pt < 900:
        weight = 0.972452
    elif 900 <= pt < 960:
        weight = 0.910183
    elif 960 <= pt < 1020:
        weight = 0.885284
    elif 1020 <= pt < 1090:
        weight = 0.950662
    elif 1090 <= pt < 1160:
        weight = 0.89605
    elif pt >= 1160:
        weight = 0.823212

    return weight


def get2016FacDownW(pt):
    weight = 1
    if pt < 170:
        weight = 1.43694
    elif 170 <= pt < 200:
        weight = 1.45181
    elif 200 <= pt < 230:
        weight = 1.42649
    elif 230 <= pt < 260:
        weight = 1.4513
    elif 260 <= pt < 290:
        weight = 1.3624
    elif 290 <= pt < 320:
        weight = 1.37472
    elif 320 <= pt < 350:
        weight = 1.3718
    elif 350 <= pt < 390:
        weight = 1.34637
    elif 390 <= pt < 430:
        weight = 1.26276
    elif 430 <= pt < 470:
        weight = 1.25762
    elif 470 <= pt < 510:
        weight = 1.26389
    elif 510 <= pt < 550:
        weight = 1.25504
    elif 550 <= pt < 590:
        weight = 1.2116
    elif 590 <= pt < 640:
        weight = 1.04846
    elif 640 <= pt < 690:
        weight = 1.052
    elif 690 <= pt < 740:
        weight = 1.09788
    elif 740 <= pt < 790:
        weight = 1.10967
    elif 790 <= pt < 840:
        weight = 1.046
    elif 840 <= pt < 900:
        weight = 1.06961
    elif 900 <= pt < 960:
        weight = 0.989001
    elif 960 <= pt < 1020:
        weight = 0.987722
    elif 1020 <= pt < 1090:
        weight = 1.07552
    elif 1090 <= pt < 1160:
        weight = 1.00802
    elif pt >= 1160:
        weight = 0.922172

    return weight

##########################################################################


def get2016EWKZ(pt):
    weight = 1
    if pt < 170:
        weight = 0.970592
    elif 170 <= pt < 200:
        weight = 0.964424
    elif 200 <= pt < 230:
        weight = 0.956695
    elif 230 <= pt < 260:
        weight = 0.948747
    elif 260 <= pt < 290:
        weight = 0.941761
    elif 290 <= pt < 320:
        weight = 0.934246
    elif 320 <= pt < 350:
        weight = 0.927089
    elif 350 <= pt < 390:
        weight = 0.919181
    elif 390 <= pt < 430:
        weight = 0.909926
    elif 430 <= pt < 470:
        weight = 0.900911
    elif 470 <= pt < 510:
        weight = 0.892561
    elif 510 <= pt < 550:
        weight = 0.884353
    elif 550 <= pt < 590:
        weight = 0.8761
    elif 590 <= pt < 640:
        weight = 0.867687
    elif 640 <= pt < 690:
        weight = 0.858047
    elif 690 <= pt < 740:
        weight = 0.849014
    elif 740 <= pt < 790:
        weight = 0.840317
    elif 790 <= pt < 840:
        weight = 0.832017
    elif 840 <= pt < 900:
        weight = 0.823545
    elif 900 <= pt < 960:
        weight = 0.814596
    elif 960 <= pt < 1020:
        weight = 0.806229
    elif 1020 <= pt < 1090:
        weight = 0.798038
    elif 1090 <= pt < 1160:
        weight = 0.789694
    elif pt >= 1160:
        weight = 0.781163

    return weight


def get2016QCDZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.47528
    elif 170 <= pt < 200:
        weight = 1.5428
    elif 200 <= pt < 230:
        weight = 1.49376
    elif 230 <= pt < 260:
        weight = 1.39119
    elif 260 <= pt < 290:
        weight = 1.40538
    elif 290 <= pt < 320:
        weight = 1.44661
    elif 320 <= pt < 350:
        weight = 1.38176
    elif 350 <= pt < 390:
        weight = 1.37381
    elif 390 <= pt < 430:
        weight = 1.29145
    elif 430 <= pt < 470:
        weight = 1.33452
    elif 470 <= pt < 510:
        weight = 1.25765
    elif 510 <= pt < 550:
        weight = 1.24265
    elif 550 <= pt < 590:
        weight = 1.24331
    elif 590 <= pt < 640:
        weight = 1.16187
    elif 640 <= pt < 690:
        weight = 1.07349
    elif 690 <= pt < 740:
        weight = 1.10748
    elif 740 <= pt < 790:
        weight = 1.06617
    elif 790 <= pt < 840:
        weight = 1.05616
    elif 840 <= pt < 900:
        weight = 1.1149
    elif 900 <= pt < 960:
        weight = 1.03164
    elif 960 <= pt < 1020:
        weight = 1.06872
    elif 1020 <= pt < 1090:
        weight = 0.981645
    elif 1090 <= pt < 1160:
        weight = 0.81729
    elif pt >= 1160:
        weight = 0.924246

    return weight


def get2016RenUpZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.38864
    elif 170 <= pt < 200:
        weight = 1.43444
    elif 200 <= pt < 230:
        weight = 1.39055
    elif 230 <= pt < 260:
        weight = 1.2934
    elif 260 <= pt < 290:
        weight = 1.30618
    elif 290 <= pt < 320:
        weight = 1.33453
    elif 320 <= pt < 350:
        weight = 1.27239
    elif 350 <= pt < 390:
        weight = 1.25981
    elif 390 <= pt < 430:
        weight = 1.18738
    elif 430 <= pt < 470:
        weight = 1.21686
    elif 470 <= pt < 510:
        weight = 1.14351
    elif 510 <= pt < 550:
        weight = 1.13305
    elif 550 <= pt < 590:
        weight = 1.12346
    elif 590 <= pt < 640:
        weight = 1.05153
    elif 640 <= pt < 690:
        weight = 0.98056
    elif 690 <= pt < 740:
        weight = 1.00651
    elif 740 <= pt < 790:
        weight = 0.96642
    elif 790 <= pt < 840:
        weight = 0.956161
    elif 840 <= pt < 900:
        weight = 0.998079
    elif 900 <= pt < 960:
        weight = 0.927421
    elif 960 <= pt < 1020:
        weight = 0.954094
    elif 1020 <= pt < 1090:
        weight = 0.883954
    elif 1090 <= pt < 1160:
        weight = 0.729348
    elif pt >= 1160:
        weight = 0.833531

    return weight


def get2016RenDownZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.55029
    elif 170 <= pt < 200:
        weight = 1.64816
    elif 200 <= pt < 230:
        weight = 1.59095
    elif 230 <= pt < 260:
        weight = 1.48394
    elif 260 <= pt < 290:
        weight = 1.49876
    elif 290 <= pt < 320:
        weight = 1.55923
    elif 320 <= pt < 350:
        weight = 1.49283
    elif 350 <= pt < 390:
        weight = 1.49279
    elif 390 <= pt < 430:
        weight = 1.39829
    elif 430 <= pt < 470:
        weight = 1.46136
    elif 470 <= pt < 510:
        weight = 1.38252
    elif 510 <= pt < 550:
        weight = 1.36081
    elif 550 <= pt < 590:
        weight = 1.37848
    elif 590 <= pt < 640:
        weight = 1.28575
    elif 640 <= pt < 690:
        weight = 1.17335
    elif 690 <= pt < 740:
        weight = 1.21841
    elif 740 <= pt < 790:
        weight = 1.17707
    elif 790 <= pt < 840:
        weight = 1.16846
    elif 840 <= pt < 900:
        weight = 1.25225
    elif 900 <= pt < 960:
        weight = 1.15246
    elif 960 <= pt < 1020:
        weight = 1.20461
    elif 1020 <= pt < 1090:
        weight = 1.09389
    elif 1090 <= pt < 1160:
        weight = 0.921666
    elif pt >= 1160:
        weight = 1.02897

    return weight


def get2016FacUpZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.48202
    elif 170 <= pt < 200:
        weight = 1.54551
    elif 200 <= pt < 230:
        weight = 1.49889
    elif 230 <= pt < 260:
        weight = 1.37781
    elif 260 <= pt < 290:
        weight = 1.39565
    elif 290 <= pt < 320:
        weight = 1.42918
    elif 320 <= pt < 350:
        weight = 1.35817
    elif 350 <= pt < 390:
        weight = 1.34966
    elif 390 <= pt < 430:
        weight = 1.26377
    elif 430 <= pt < 470:
        weight = 1.30245
    elif 470 <= pt < 510:
        weight = 1.22426
    elif 510 <= pt < 550:
        weight = 1.20589
    elif 550 <= pt < 590:
        weight = 1.20726
    elif 590 <= pt < 640:
        weight = 1.12468
    elif 640 <= pt < 690:
        weight = 1.04086
    elif 690 <= pt < 740:
        weight = 1.06898
    elif 740 <= pt < 790:
        weight = 1.02577
    elif 790 <= pt < 840:
        weight = 1.01591
    elif 840 <= pt < 900:
        weight = 1.06486
    elif 900 <= pt < 960:
        weight = 0.986133
    elif 960 <= pt < 1020:
        weight = 1.01224
    elif 1020 <= pt < 1090:
        weight = 0.937427
    elif 1090 <= pt < 1160:
        weight = 0.77774
    elif pt >= 1160:
        weight = 0.87435

    return weight


def get2016FacDownZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.47257
    elif 170 <= pt < 200:
        weight = 1.54471
    elif 200 <= pt < 230:
        weight = 1.49357
    elif 230 <= pt < 260:
        weight = 1.41154
    elif 260 <= pt < 290:
        weight = 1.42113
    elif 290 <= pt < 320:
        weight = 1.46933
    elif 320 <= pt < 350:
        weight = 1.41145
    elif 350 <= pt < 390:
        weight = 1.40379
    elif 390 <= pt < 430:
        weight = 1.32597
    elif 430 <= pt < 470:
        weight = 1.37278
    elif 470 <= pt < 510:
        weight = 1.2973
    elif 510 <= pt < 550:
        weight = 1.28653
    elif 550 <= pt < 590:
        weight = 1.28624
    elif 590 <= pt < 640:
        weight = 1.20542
    elif 640 <= pt < 690:
        weight = 1.11142
    elif 690 <= pt < 740:
        weight = 1.1518
    elif 740 <= pt < 790:
        weight = 1.11242
    elif 790 <= pt < 840:
        weight = 1.1022
    elif 840 <= pt < 900:
        weight = 1.17247
    elif 900 <= pt < 960:
        weight = 1.08313
    elif 960 <= pt < 1020:
        weight = 1.13457
    elif 1020 <= pt < 1090:
        weight = 1.03169
    elif 1090 <= pt < 1160:
        weight = 0.862166
    elif pt >= 1160:
        weight = 0.981254

    return weight
