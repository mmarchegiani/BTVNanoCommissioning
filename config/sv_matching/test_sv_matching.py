from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.mutag_oneMuAK8_processor import mutagAnalysisOneMuonInAK8Processor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
from config.fatjet_base.custom.cuts import mutag_fatjet_sel, mutag_subjet_sel, get_nObj_minmsd, get_flavor
from config.fatjet_base.custom.functions import get_HLTsel

pt_min = 450.

samples = ["QCD_MuEnriched",
           "QCD_HT",
           "VJets",
           "SingleTop_ttbar",
           "GluGluHToBB",
           "GluGluHToCC",
           "DATA"
           ]

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
                  "datasets/MC_QCD_HT_RunIISummer20UL.json",
                  "datasets/MC_VJets_RunIISummer20UL.json",
                  "datasets/MC_top_RunIISummer20UL_local.json",
                  "datasets/DATA_BTagMu_RunIISummer20UL_local.json"
                  ],
        "filter" : {
            "samples": samples,
            "samples_exclude" : [],
            "year": ['2018']
        },
        "subsamples" : subsamples
    },

    # Input and output files
    "workflow" : mutagAnalysisOneMuonInAK8Processor,
    "output"   : "output/test/sv_matching",
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
        "walltime"       : "12:00:00",
        "mem_per_worker" : "12GB", # GB
        "exclusive"      : False,
        "chunk"          : 10000,
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
    "preselections" : [get_nObj_min(1, pt_min, "FatJetGood")],
    "categories": common_cats,

    "weights": {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup",
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
                "inclusive": ["pileup",
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
                                         label=r"FatJet $p_{T}$ [GeV]", bins=list(range(int(pt_min), 1010, 10)))]),
        "FatJetGood_msoftdrop" : HistConf([Axis(name=f"FatJetGood_msoftdrop", coll="FatJetGood", field="msoftdrop",
                                         label=r"FatJet $m_{SD}$ [GeV]", bins=list(range(40, 410, 10)))]),
        "FatJetGood_logsumcorrmass": HistConf(
            [ Axis(coll="FatJetGood", field="logsumcorrSVmass", label=r"log($\sum({m^{corr}_{SV}})$)", bins=42, start=-2.4, stop=6) ]
        ),
        "nSVMatchedToFatJetGood": HistConf(
            [ Axis(coll="FatJetGood", field="nSVMatchedToFatJetGood", label="Number of SV inside AK8 jet", bins=10, start=0, stop=10) ]
        ),
        "FatJetGood_logsumcorrSVmass_tau21": HistConf(
            [ Axis(coll="FatJetGood", field="logsumcorrSVmass", label=r"log($\sum({m^{corr}_{SV}})$)", bins=42, start=-2.4, stop=6),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", type="variable", bins=[0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1]) ]
        ),
    },

    "columns" : {}

}

# We reweigh histograms differently depending whether:
# - The histogram contains the whole jet collection ("all")
#cfg["workflow_options"]["histograms_to_reweigh"]["by_pos"]["all"] = [name for name in cfg["variables"].keys() if name.startswith("FatJetGood_") and not name.endswith(("_1", "_2"))]
