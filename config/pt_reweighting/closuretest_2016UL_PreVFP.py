from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.mutag_processor import mutagAnalysisProcessor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
from config.fatjet_base.custom.cuts import mutag_fatjet_sel, mutag_subjet_sel, get_nObj_minmsd, get_flavor
from config.fatjet_base.custom.functions import get_HLTsel

samples = ["QCD_MuEnriched",
           "VJets",
           "SingleTop_ttbar",
           "GluGluHToBB",
           "GluGluHToCC",
           "DATA"
           ]

pt_min = 350.0
msd = 40.0

common_cats = {
    "inclusive" : [passthrough],
    "mutag_fatjet_nmu-1" : [mutag_fatjet_sel(nmu=1)],
    "mutag_fatjet_nmu-2" : [mutag_fatjet_sel(nmu=2)],
    "mutag_subjet_nmu-2" : [mutag_subjet_sel(unique_matching=False)],
    "mutag_subjet_unique_nmu-2" : [mutag_subjet_sel(unique_matching=True)]
}

subsamples = {}
for s in filter(lambda x: 'DATA' not in x, samples):
    subsamples[s] = {f"{s}_{f}" : [get_flavor(f)] for f in ['l', 'c', 'b', 'cc', 'bb']}

cfg =  {
    "dataset" : {
        "jsons": ["datasets/MC_QCD_MuEnriched_RunIISummer20UL_local.json",
                  "datasets/MC_VJets_RunIISummer20UL.json",
                  "datasets/MC_top_RunIISummer20UL_local.json",
                  "datasets/MC_GluGluH_RunIISummer20UL_local.json",
                  "datasets/DATA_BTagMu_RunIISummer20UL_local.json"
                  ],
        "filter" : {
            "samples": samples,
            "samples_exclude" : [],
            "year": ['2016_PreVFP']
        },
        "subsamples" : subsamples
    },

    # Input and output files
    "workflow" : mutagAnalysisProcessor,
    "output"   : "output/pocket_coffea/closure_test/closuretest_2016_PreVFP_bintau05",
    "workflow_options" : {
        "histograms_to_reweigh" : {
            "by_pos" : {
                'all' : [],
                '1' : [],
                '2' : [],
            }
        },
    },

    "run_options" : {
        "executor"       : "dask/slurm",
        "workers"        : 1,
        "scaleout"       : 200,
        "queue"          : "standard",
        "walltime"       : "8:00:00",
        "mem_per_worker" : "6GB", # GB
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
    "preselections" : [get_nObj_min(1, 350., "FatJetGood")],
    "categories": common_cats,

    "weights": {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup", "sf_L1prefiring",
                          "sf_ptetatau21_reweighting"
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
                "inclusive": ["pileup", "sf_L1prefiring",
                              "sf_ptetatau21_reweighting"
                              ],
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
        "FatJetGood_pt" : HistConf([Axis(name=f"FatJetGood_pt", coll="FatJetGood", field="pt",
                                         label=r"FatJet $p_{T}$ [GeV]", bins=list(range(350, 1010, 10)))]),
        "FatJetGood_msoftdrop" : HistConf([Axis(name=f"FatJetGood_msoftdrop", coll="FatJetGood", field="msoftdrop",
                                         label=r"FatJet $m_{SD}$ [GeV]", bins=list(range(40, 410, 10)))]),
    },

    "columns" : {}

}

axis_pos = Axis(name=f"FatJetGood_pos", coll="FatJetGood", field="pos", type="int", label=r"FatJet position", bins=2, start=0, stop=2)

hists_2d = {}
for histname, hist_cfg in cfg["variables"].items():
    hists_2d[f"{histname}_pos"] = HistConf([axis_pos] + hist_cfg.axes)

for histname, hist_cfg in hists_2d.items():
    cfg["variables"][histname] = hists_2d[histname]

# We reweigh histograms differently depending whether:
# - The histogram contains the whole jet collection ("all")
# - The histogram contains the leading jet collection ("1")
# - The histogram contains the subleading jet collection ("2")
cfg["workflow_options"]["histograms_to_reweigh"]["by_pos"]["all"] = [name for name in cfg["variables"].keys() if name.startswith("FatJetGood_") and not name.endswith(("_1", "_2"))]
