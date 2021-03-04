import collections
import awkward1 as ak

lumi = {
    2016 : 35.5,
    2017 : 41.5,
    2018 : 59.2,
}

# Cross-sections in pb
xsecs = {
    "QCD_Pt-15to20_MuEnrichedPt5"    : 2799000.0,
    "QCD_Pt-20to30_MuEnrichedPt5"    : 2526000.0,
    "QCD_Pt-30to50_MuEnrichedPt5"    : 1362000.0,
    "QCD_Pt-50to80_MuEnrichedPt5"    : 376600.0,
    "QCD_Pt-80to120_MuEnrichedPt5"   : 88930.0,
    "QCD_Pt-120to170_MuEnrichedPt5"  : 21230.0,
    "QCD_Pt-170to300_MuEnrichedPt5"  : 7055.0,
    "QCD_Pt-300to470_MuEnrichedPt5"  : 619.3,    
    "QCD_Pt-470to600_MuEnrichedPt5"  : 59.24,
    "QCD_Pt-600to800_MuEnrichedPt5"  : 18.21,
    "QCD_Pt-800to1000_MuEnrichedPt5" : 3.275,
    "QCD_Pt-1000toInf_MuEnrichedPt5" : 1.078,

    "GluGluHToBB_M-125_13TeV"  : 27.8,
    "GluGluHToCC_M-125_13TeV"  : 27.8,
}

def rescale(accumulator, xsec, lumi):
    """Scale by lumi"""
    lumi = 1000*lumi    # Convert lumi to pb^-1
    from coffea import hist
    scale = {}
    print("Scaling:")
    for dataset, dataset_sumw in collections.OrderedDict(sorted(accumulator['sumw'].items())).items():
        #dataset_key = dataset.lstrip("/").split("/")[0]
        dataset_key = dataset
        if dataset_key in xsecs:
            print(" ", dataset_key, "\t", xsecs[dataset_key], "pb")
            scale[dataset] = lumi*xsecs[dataset_key]/dataset_sumw
        else:
            print(" ", "X ", dataset_key)
            scale[dataset] = 0#lumi / dataset_sumw

    for h in accumulator.values():
        if isinstance(h, hist.Hist):
            h.scale(scale, axis="dataset")
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