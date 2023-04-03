from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.pt_reweighting import ptReweightProcessor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
from pocket_coffea.lib.weights_manager import WeightCustom
from config.fatjet_base.custom.cuts import twojets_presel, mutag_sel, get_ptmsd, get_ptmsdtau, get_nObj_minmsd, get_flavor
from config.fatjet_base.custom.functions import get_HLTsel
import numpy as np

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
    "pt450msd40tau50" : [get_ptmsdtau(450., 40., 0.50)],
    "pt450msd40tau50_mutag" : [get_ptmsdtau(450., 40., 0.50), mutag_sel],
    "pt450msd40tau55" : [get_ptmsdtau(450., 40., 0.55)],
    "pt450msd40tau55_mutag" : [get_ptmsdtau(450., 40., 0.55), mutag_sel],
    "pt450msd40tau60" : [get_ptmsdtau(450., 40., 0.60)],
    "pt450msd40tau60_mutag" : [get_ptmsdtau(450., 40., 0.60), mutag_sel],
}

samples = ["QCD_MuEnriched",
           "VJets",
           "DATA"
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
    "workflow" : ptReweightProcessor,
    "output"   : "output/pocket_coffea/pt_reweighting/pt_eta_reweighting_2018UL_twojets",
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
                "inclusive": [ ],
                "bycategory" : {
                }
            },
        "bysample": {
        }    
        },
        
    },

   "variables":
    {
        **fatjet_hists(coll="FatJetGood"),
        **sv_hists(coll="events"),
        **count_hist(name="nFatJets", coll="FatJetGood",bins=10, start=0, stop=10),
        "FatJetGood_pt_eta": HistConf(
            [ Axis(name="FatJetGood_pt", coll="FatJetGood", field="pt", label=r"Leading FatJet $p_{T}$ [GeV]", bins=150, start=0, stop=1500),
              Axis(name="FatJetGood_eta", coll="FatJetGood", field="eta", label=r"Leading FatJet $\eta$", bins=40, start=-4, stop=4) ]
        ),
        "FatJetGood_pt_eta_bineta0p40": HistConf(
            [ Axis(name="FatJetGood_pt", coll="FatJetGood", field="pt", label=r"Leading FatJet $p_{T}$ [GeV]", bins=150, start=0, stop=1500),
              Axis(name="FatJetGood_eta", coll="FatJetGood", field="eta", label=r"Leading FatJet $\eta$", bins=20, start=-4, stop=4) ]
        ),
        "FatJetGood_pt_eta_binpt20": HistConf(
            [ Axis(name="FatJetGood_pt", coll="FatJetGood", field="pt", label=r"Leading FatJet $p_{T}$ [GeV]", bins=75, start=0, stop=1500),
              Axis(name="FatJetGood_eta", coll="FatJetGood", field="eta", label=r"Leading FatJet $\eta$", bins=40, start=-4, stop=4) ]
        ),
        "FatJetGood_pt_eta_binpt20_bineta0p40": HistConf(
            [ Axis(name="FatJetGood_pt", coll="FatJetGood", field="pt", label=r"Leading FatJet $p_{T}$ [GeV]", bins=75, start=0, stop=1500),
              Axis(name="FatJetGood_eta", coll="FatJetGood", field="eta", label=r"Leading FatJet $\eta$", bins=20, start=-4, stop=4) ]
        ),
    },

    "columns" : {}

}
