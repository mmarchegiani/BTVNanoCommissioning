from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.pt_reweighting import ptReweightProcessor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
from config.fatjet_base.custom.cuts import twojets_presel, mutag_fatjet_sel, mutag_subjet_sel, get_ptmsd, get_ptmsdtau, get_nObj_minmsd, get_flavor
from config.fatjet_base.custom.functions import get_HLTsel

common_cats = {
    "inclusive" : [passthrough],
    "pt450msd40_mutag_fatjet_nmu-1" : [get_ptmsd(450., msd), mutag_fatjet_sel(nmu=1)],
    "pt450msd40_mutag_fatjet_nmu-2" : [get_ptmsd(450., msd), mutag_fatjet_sel(nmu=2)],
    "pt450msd40_mutag_subjet" : [get_ptmsd(450., msd), mutag_subjet_sel(unique_matching=False)],
    "pt450msd40_mutag_subjet_unique" : [get_ptmsd(450., msd), mutag_subjet_sel(unique_matching=True)],
}

samples = ["QCD_MuEnriched",
           "VJets",
           "SingleTop_ttbar",
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
    "output"   : "output/pocket_coffea/pt_reweighting/pt_eta_tau21_reweighting_2018",
    "workflow_options" : {"histograms_to_reweigh" : []},

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
             get_nObj_min(1, 3., "Muon"),
             get_HLTsel("mutag")],
    "save_skimmed_files" : None,
    "preselections" : [twojets_presel(pt=250, msd=40)],
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
            [ Axis(name="FatJetGood_pt", coll="FatJetGood", field="pt", label=r"FatJet $p_{T}$ [GeV]", bins=150, start=0, stop=1500),
              Axis(name="FatJetGood_eta", coll="FatJetGood", field="eta", label=r"FatJet $\eta$", bins=40, start=-4, stop=4) ]
        ),
        "FatJetGood_pt_eta_bineta0p40": HistConf(
            [ Axis(name="FatJetGood_pt", coll="FatJetGood", field="pt", label=r"FatJet $p_{T}$ [GeV]", bins=150, start=0, stop=1500),
              Axis(name="FatJetGood_eta", coll="FatJetGood", field="eta", label=r"FatJet $\eta$", bins=20, start=-4, stop=4) ]
        ),
        "FatJetGood_pt_eta_binpt20": HistConf(
            [ Axis(name="FatJetGood_pt", coll="FatJetGood", field="pt", label=r"FatJet $p_{T}$ [GeV]", bins=75, start=0, stop=1500),
              Axis(name="FatJetGood_eta", coll="FatJetGood", field="eta", label=r"FatJet $\eta$", bins=40, start=-4, stop=4) ]
        ),
        "FatJetGood_pt_eta_binpt20_bineta0p40": HistConf(
            [ Axis(name="FatJetGood_pt", coll="FatJetGood", field="pt", label=r"FatJet $p_{T}$ [GeV]", bins=75, start=0, stop=1500),
              Axis(name="FatJetGood_eta", coll="FatJetGood", field="eta", label=r"FatJet $\eta$", bins=20, start=-4, stop=4) ]
        ),
        "FatJetGood_pt_eta_tau21": HistConf(
            [ Axis(name="FatJetGood_pt", coll="FatJetGood", field="pt", label=r"FatJet $p_{T}$ [GeV]", bins=150, start=0, stop=1500),
              Axis(name="FatJetGood_eta", coll="FatJetGood", field="eta", label=r"FatJet $\eta$", bins=40, start=-4, stop=4),
              Axis(name="FatJetGood_tau21", coll="FatJetGood", field="tau21", type="variable", label=r"FatJet $\tau_{21}$",
                   bins=[0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1]) ]
        ),
    },

    "columns" : {}

}
