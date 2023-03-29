from copy import deepcopy
from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.fatjet_base import fatjetBaseProcessor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
from pocket_coffea.lib.categorization import StandardSelection, CartesianSelection
from pocket_coffea.lib.columns_manager import ColOut
from config.fatjet_base.custom.cuts import twojets_presel, mutag_sel, get_ptmsd, get_ptmsdtau, get_nObj_minmsd, get_flavor
from config.fatjet_base.custom.functions import get_HLTsel
from config.fatjet_base.custom.parameters.parameters import AK8Taggers

common_cats = {
    "inclusive" : [passthrough],
    "pt450msd40" : [get_ptmsd(450., 40.)],
    "pt450msd40_mutag" : [get_ptmsd(450., 40.), mutag_sel],
    "pt450msd40tau20" : [get_ptmsdtau(450., 40., 0.20)],
    "pt450msd40tau20_mutag" : [get_ptmsdtau(450., 40., 0.20), mutag_sel],
    "pt450msd40tau25" : [get_ptmsdtau(450., 40., 0.25)],
    "pt450msd40tau25_mutag" : [get_ptmsdtau(450., 40., 0.25), mutag_sel],
    "pt450msd40tau30" : [get_ptmsdtau(450., 40., 0.30)],
    "pt450msd40tau30_mutag" : [get_ptmsdtau(450., 40., 0.30), mutag_sel],
    "pt450msd40tau35" : [get_ptmsdtau(450., 40., 0.35)],
    "pt450msd40tau35_mutag" : [get_ptmsdtau(450., 40., 0.35), mutag_sel],
    "pt450msd40tau40" : [get_ptmsdtau(450., 40., 0.40)],
    "pt450msd40tau40_mutag" : [get_ptmsdtau(450., 40., 0.40), mutag_sel],
    "pt450msd40tau45" : [get_ptmsdtau(450., 40., 0.45)],
    "pt450msd40tau45_mutag" : [get_ptmsdtau(450., 40., 0.45), mutag_sel],
}

samples = ["QCD_MuEnriched",
           "VJets",
           "GluGluHToBB",
           "GluGluHToCC",
           ]
subsamples = {}
for s in filter(lambda x: 'DATA' not in x, samples):
    subsamples[s] = {f"{s}_{f}" : [get_flavor(f)] for f in ['l', 'c', 'b', 'cc', 'bb']}

cfg =  {
    "dataset" : {
        "jsons": ["datasets/skim/datasets_definition_skim.json"],
        "filter" : {
            "samples": samples,
            "samples_exclude" : [],
            "year": ['2018']
        },
        "subsamples": subsamples
    },
    

    # Input and output files
    "workflow" : fatjetBaseProcessor,
    "output"   : "output/pocket_coffea/ggH_proxy/ggH_proxy_2018UL_tau21",
    "workflow_options" : {},

    "run_options" : {
        "executor"       : "dask/slurm",
        "workers"        : 1,
        "scaleout"       : 200,
        "queue"          : "standard",
        "walltime"       : "12:00:00",
        "mem_per_worker" : "12GB", # GB
        "exclusive"      : False,
        "chunk"          : 100000,
        "retries"        : 50,
        "treereduction"  : 10,
        "max"            : None,
        "skipbadfiles"   : None,
        "voms"           : None,
        "limit"          : None,
        "adapt"          : False,
        "env"            : "conda",
    },

    # Cuts and plots settings
    "finalstate" : "mutag",
    "skim": [get_nObj_min(1, 200., "FatJet"),
             get_nObj_minmsd(1, 30., "FatJet"),
             get_nObj_min(2, 3., "Muon"),
             get_HLTsel("mutag")],
    "save_skimmed_files" : None,
    "preselections" : [twojets_presel],
    "categories": common_cats,

    "weights": {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup"#, "sf_L1prefiring"
                          ],
            "bycategory" : {
            }
        },
        "bysample": {
        }
    },

    "variations": {
        "weights": {
            "common": {
                "inclusive": [ ],#, "sf_L1prefiring" ],
                "bycategory" : {
                }
            },
        "bysample": {
        }    
        },
        "shape": {
            "common":{
                "inclusive": [ ]
            }
        }
    },

   "variables":
    {
        **muon_hists(coll="MuonGood", pos=0),
        **muon_hists(coll="MuonGood", pos=1),
        #**jet_hists(coll="JetGood"),
        #**jet_hists(coll="JetGood", pos=0),
        **fatjet_hists(coll="FatJetGood"),
        **fatjet_hists(coll="FatJetGood", pos=0),
        **fatjet_hists(coll="FatJetGood", pos=1),
        "FatJetGood_hadronFlavour": HistConf(
            [ Axis(coll="FatJetGood", field="hadronFlavour", label="hadronFlavour", bins=10, start=0, stop=10) ]
        ),
        "FatJetGood_nBHadrons": HistConf(
            [ Axis(coll="FatJetGood", field="nBHadrons", label="nBHadrons", bins=10, start=0, stop=10) ]
        ),
        "FatJetGood_nCHadrons": HistConf(
            [ Axis(coll="FatJetGood", field="nCHadrons", label="nCHadrons", bins=10, start=0, stop=10) ]
        ),
        **sv_hists(coll="events"),
        **sv_hists(coll="events", pos=0),
        **sv_hists(coll="events", pos=1),
        #**count_hist(name="nElectronGood", coll="ElectronGood",bins=10, start=0, stop=10),
        **count_hist(name="nMuonGood", coll="MuonGood",bins=10, start=0, stop=10),
        #**count_hist(name="nJets", coll="JetGood",bins=10, start=0, stop=10),
        **count_hist(name="nFatJets", coll="FatJetGood",bins=10, start=0, stop=10),
        **count_hist(name="nSV", coll="SV",bins=10, start=0, stop=10),
        "nmusj_fatjet1": HistConf(
            [ Axis(coll="events", field="nmusj_fatjet1", label=r"$N_{\mu, J1}$", bins=10, start=0, stop=10) ]
        ),
        "nmusj_fatjet2": HistConf(
            [ Axis(coll="events", field="nmusj_fatjet2", label=r"$N_{\mu, J2}$", bins=10, start=0, stop=10) ]
        ),
        "FatJetGood_tau21_particleNetMD_Xcc_QCD": HistConf(
            [ Axis(coll="FatJetGood", field="particleNetMD_Xcc_QCD", label="particleNetMD Xbb/(Xbb + QCD)", bins=20, start=0, stop=1),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", bins=20, start=0, stop=1) ]
        ),
        "FatJetGood_tau21_particleNetMD_Xbb_QCD": HistConf(
            [ Axis(coll="FatJetGood", field="particleNetMD_Xbb_QCD", label="particleNetMD Xbb/(Xbb + QCD)", bins=20, start=0, stop=1),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", bins=20, start=0, stop=1) ]
        ),
        "FatJetGood_tau21_btagDDBvLV2": HistConf(
            [ Axis(coll="FatJetGood", field="btagDDBvLV2", label="btagDDBvLV2", bins=20, start=0, stop=1),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", bins=20, start=0, stop=1) ]
        ),
        "FatJetGood_tau21_btagDDCvLV2": HistConf(
            [ Axis(coll="FatJetGood", field="btagDDCvLV2", label="btagDDCvLV2", bins=20, start=0, stop=1),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", bins=20, start=0, stop=1) ]
        ),
        "FatJetGood_tau21_deepTagMD_ZHbbvsQCD": HistConf(
            [ Axis(coll="FatJetGood", field="deepTagMD_ZHbbvsQCD", label="deepTagMD_ZHbbvsQCD", bins=20, start=0, stop=1),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", bins=20, start=0, stop=1) ]
        ),
        "FatJetGood_tau21_deepTagMD_ZHccvsQCD": HistConf(
            [ Axis(coll="FatJetGood", field="deepTagMD_ZHccvsQCD", label="deepTagMD_ZHccvsQCD", bins=20, start=0, stop=1),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", bins=20, start=0, stop=1) ]
        ),
        "FatJetGood_tau21_btagHbb": HistConf(
            [ Axis(coll="FatJetGood", field="btagHbb", label="btagHbb", bins=20, start=0, stop=1),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", bins=20, start=0, stop=1) ]
        ),
    },

    "columns" : {
        "common" : {
            "inclusive": [
                ColOut("FatJetGood", ["pt", "eta", "phi", "msoftdrop",
                                    "tau21", "hadronFlavour", "nBHadrons", "nCHadrons",
                                    "particleNetMD_Xbb_QCD", "particleNetMD_Xcc_QCD",
                                    "btagDDBvLV2", "btagDDCvLV2"],
                ),
            ]
        }
    }

}

for tagger in AK8Taggers:
    setting = deepcopy(default_axis_settings[f"fatjet_{tagger}"])
    setting["coll"] = "FatJetGood"
    for nbins in [10, 20, 40]:
        setting["bins"] = nbins
        cfg["variables"][f"FatJetGood_{tagger}_{nbins}bins"] = HistConf(axes=[Axis(**setting)])
