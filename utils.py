import collections
import numpy as np
import awkward as ak

lumi = {
    2016 : 35.5,
    2017 : 41.5,
    2018 : 59.2,
}

# Cross-sections in pb
xsecs = {
    #"QCD_Pt-15to20_MuEnrichedPt5"    : 2799000.0,
    #"QCD_Pt-20to30_MuEnrichedPt5"    : 2526000.0,
    #"QCD_Pt-30to50_MuEnrichedPt5"    : 1362000.0,
    #"QCD_Pt-50to80_MuEnrichedPt5"    : 376600.0,
    #"QCD_Pt-80to120_MuEnrichedPt5"   : 88930.0,
    "QCD_Pt-120to170_MuEnrichedPt5"  : 21230.0,
    "QCD_Pt-170to300_MuEnrichedPt5"  : 7055.0,
    "QCD_Pt-300to470_MuEnrichedPt5"  : 619.3,
    "QCD_Pt-470to600_MuEnrichedPt5"  : 59.24,
    "QCD_Pt-600to800_MuEnrichedPt5"  : 18.21,
    "QCD_Pt-800to1000_MuEnrichedPt5" : 3.275,
    "QCD_Pt-1000toInf_MuEnrichedPt5" : 1.078,
    "QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8"  : 7055.0,
    "QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8"  : 619.3,
    "QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8"  : 59.24,
    "QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8"  : 18.21,
    "QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8" : 3.275,
    "QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8" : 1.078,

    "GluGluHToBB_M-125_13TeV"  : 27.8,
    "GluGluHToCC_M-125_13TeV"  : 27.8,
}

histogram_settings = {
    'variables' : {
        'fatjet_pt'    : {'binning' : {'n_or_arr' : 125, 'lo' : 0,      'hi' : 1000},  'xlim' : {'xmin' : 0, 'xmax' : 1000}},
        'fatjet_eta'   : {'binning' : {'n_or_arr' : 30,  'lo' : -3,     'hi' : 3},     'xlim' : {'xmin' : -3, 'xmax' : 3}},
        'fatjet_phi'   : {'binning' : {'n_or_arr' : 30,  'lo' : -np.pi, 'hi' : np.pi}, 'xlim' : {'xmin' : -np.pi, 'xmax' : np.pi}},
        #'fatjet_phi'   : {'binning' : {'n_or_arr' : 30,  'lo' : -3,     'hi' : 3},     'xlim' : {'xmin' : -3, 'xmax' : 3}},
        'fatjet_mass'  : {'binning' : {'n_or_arr' : 150, 'lo' : 0,      'hi' : 300},   'xlim' : {'xmin' : 0, 'xmax' : 300}},
        'nsv1'         : {'binning' : {'n_or_arr' : 30,  'lo' : 0,      'hi' : 30},    'xlim' : {'xmin' : 0, 'xmax' : 10}},
        'nsv2'         : {'binning' : {'n_or_arr' : 30,  'lo' : 0,      'hi' : 30},    'xlim' : {'xmin' : 0, 'xmax' : 10}},
        'nmusj1'       : {'binning' : {'n_or_arr' : 30,  'lo' : 0,      'hi' : 30},    'xlim' : {'xmin' : 0, 'xmax' : 10}},
        'nmusj2'       : {'binning' : {'n_or_arr' : 30,  'lo' : 0,      'hi' : 30},    'xlim' : {'xmin' : 0, 'xmax' : 10}},
    }
}

def rescale(accumulator, xsecs=xsecs, lumi=lumi, data="BTagMu"):
#def rescale(accumulator, xsecs=xsecs, data="BTagMu"):
    """Scale by lumi"""
    #lumi = 1000*lumi    # Convert lumi to pb^-1
    from coffea import hist
    scale = {}
    sumxsecs = ak.sum(xsecs.values())
    #N_data = accumulator['nbtagmu_event_level'][data]
    #print("N_data =", N_data)
    for dataset, N_mc in collections.OrderedDict(sorted(accumulator['sumw'].items())).items():
        if dataset in xsecs:
            print(" ", dataset, "\t", N_mc, "events\t", xsecs[dataset], "pb")
            #scale[dataset] = (xsecs[dataset]/sumxsecs)*(N_data/N_mc)
            scale[dataset] = (xsecs[dataset]*lumi)/N_mc
        else:
            print(" ", "X ", dataset)
            scale[dataset] = 0#lumi / N_mc
    print(scale)

    datasets_mc = [item for item in list(xsecs.keys()) if not 'GluGlu' in item]
    for h in accumulator.values():
        if isinstance(h, hist.Hist):
            h.scale(scale,       axis="dataset")
            N_data = ak.sum(h[data].values().values())
            N_mc = ak.sum(h[datasets_mc].sum('dataset', 'flavor').values().values())
            #scaletodata = dict(zip(scale.keys(), len(scale)*[1./N_data]))
            scaletodata = dict(zip(scale.keys(), len(scale)*[N_data/N_mc]))
            h.scale(scaletodata, axis="dataset")
    return accumulator

def get_nsv(sj, sv, R=0.4):

    sv_dr = sj.delta_r(sv)
    nsv = ak.count(sv_dr[sv_dr < R], axis=1)

    return nsv

"""
def xSecReader(fname):
   # Probably unsafe
    with open(fname) as file:
        lines = file.readlines()
    lines = [l.strip("\n") for l in lines if not l.startswith("#")]
    lines = [l.split("#")[0] for l in lines if len(l) > 5]

    _dict = {}
    for line in lines:
        key = line.split()[0]
        valuex = line.split()[1:]
        if len(valuex) > 1:
            value = "".join(valuex)
        else:
            value = valuex[0]
        _dict[key] = float(eval(value))
    return _dict

xsecs = {
    "/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"    : 2799000.0,
    "/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"    : 2526000.0,
    "/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"    : 1362000.0,
    "/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"    : 376600.0,
    "/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"   : 88930.0,
    "/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  : 21230.0,
    "/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  : 7055.0,
    "/QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  : 619.3,
    "/QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  : 59.24,
    "/QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  : 18.21,
    "/QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" : 3.275,
    "/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" : 1.078,

    "/GluGluHToCC_M-125_13TeV_powheg_MINLO_NNLOPS_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"  : 27.8,
}

"""
