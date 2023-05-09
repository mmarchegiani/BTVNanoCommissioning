from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.templates import templatesProcessor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
from config.fatjet_base.custom.cuts import twojets_presel, mutag_fatjet_sel, get_ptmsd, get_nObj_minmsd, get_flavor
from config.fatjet_base.custom.functions import get_HLTsel

samples = ["QCD_MuEnriched",
           "VJets",
           "SingleTop_ttbar",
           "DATA"
           ]

pt_min = 350.0
msd = 40.0

common_cats = {
    "inclusive" : [passthrough],
    f"pt{int(pt_min)}msd{int(msd)}_mutag_fatjet_nmu-1" : [get_ptmsd(pt_min, msd), mutag_fatjet_sel(nmu=1)]
}

subsamples = {}
for s in filter(lambda x: 'DATA' not in x, samples):
    subsamples[s] = {f"{s}_{f}" : [get_flavor(f)] for f in ['l', 'c', 'b', 'cc', 'bb']}

cfg =  {
    "dataset" : {
        "jsons": ["datasets/MC_QCD_MuEnriched_RunIISummer20UL_local.json",
                  "datasets/MC_VJets_RunIISummer20UL.json",
                  "datasets/MC_top_RunIISummer20UL_local.json",
                  "datasets/DATA_BTagMu_RunIISummer20UL_local.json"],
        "filter" : {
            "samples": samples,
            "samples_exclude" : [],
            "year": ['2016_PostVFP']
        },
        "subsamples" : subsamples
    },

    # Input and output files
    "workflow" : templatesProcessor,
    "output"   : "output/pocket_coffea/closure_test/closuretest_2016_PostVFP",
    "workflow_options" : {"histograms_to_reweigh" : [],
                          "reweighting_scheme"    : "ptetatau21",
                          "reweighting_map"       : "/work/mmarcheg/BTVNanoCommissioning/config/fatjet_base/custom/parameters/pt_reweighting/pt_eta_tau21_reweighting_2016_PostVFP_twojets_pt350/FatJetGoodNMuon1_pt_eta_tau21_2016_PostVFP_reweighting.json"},

    "run_options" : {
        "executor"       : "dask/slurm",
        "workers"        : 1,
        "scaleout"       : 200,
        "queue"          : "standard",
        "walltime"       : "6:00:00",
        "mem_per_worker" : "4GB", # GB
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
    "categories": {
        "inclusive" : [passthrough]
    },

    "weights": {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup", "sf_L1prefiring"
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
        "FatJetGood_pt" : HistConf([Axis(name=f"FatJetGood_pt", coll="FatJetGood", field="pt",
                                         label=r"FatJet $p_{T}$ [GeV]", bins=list(range(350, 1010, 10)))]),
        "FatJetGood_msoftdrop" : HistConf([Axis(name=f"FatJetGood_msoftdrop", coll="FatJetGood", field="msoftdrop",
                                         label=r"FatJet $m_{SD}$ [GeV]", bins=list(range(40, 410, 10)))]),
        #**sv_hists(coll="events"),
        **count_hist(name="nFatJetGood", coll="FatJetGood",bins=10, start=0, stop=10),
        "nMuonGoodMatchedToFatJetGood": HistConf(
            [ Axis(coll="FatJetGood", field="nMuonGoodMatchedToFatJetGood", label="Number of muons inside AK8 jet", bins=4, start=0, stop=4) ]
        ),
        "nMuonGoodMatchedToSubJet": HistConf(
            [ Axis(coll="FatJetGood", field="nMuonGoodMatchedToSubJet", label="Number of muons matched to subjet", bins=3, start=0, stop=3) ]
        ),
        "nMuonGoodMatchedUniquelyToSubJet": HistConf(
            [ Axis(coll="FatJetGood", field="nMuonGoodMatchedUniquelyToSubJet", label="Number of muons matched uniquely to subjet", bins=3, start=0, stop=3) ]
        ),
    },

    "columns" : {}

}

cfg["workflow_options"]["histograms_to_reweigh"] = [name for name in cfg["variables"].keys() if name.startswith("FatJetGood_")]
