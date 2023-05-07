from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.pt_reweighting import ptReweightProcessor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
from config.fatjet_base.custom.cuts import get_nObj_minmsd
from config.fatjet_base.custom.functions import get_HLTsel

samples = ["QCD_MuEnriched",
           "VJets",
           "SingleTop_ttbar",
           "DATA"
           ]

cfg =  {
    "dataset" : {
        "jsons": ["datasets/MC_QCD_MuEnriched_RunIISummer20UL_local.json",
                  "datasets/MC_VJets_RunIISummer20UL.json",
                  "datasets/MC_top_RunIISummer20UL_local.json",
                  "datasets/DATA_BTagMu_RunIISummer20UL_local.json"],
        "filter" : {
            "samples": samples,
            "samples_exclude" : [],
            "year": ['2016_PreVFP']
        }
    },

    # Input and output files
    "workflow" : ptReweightProcessor,
    "output"   : "output/pocket_coffea/pt_reweighting/pt_eta_tau21_reweighting_2016_PreVFP_twojets_pt350_leadsublead",
    "workflow_options" : {},

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
    "preselections" : [passthrough],
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
        **count_hist(name="nFatJetGood", coll="FatJetGood",bins=10, start=0, stop=10),
        **count_hist(coll="FatJetGoodNMuon1",bins=10, start=0, stop=10),
        **count_hist(coll="FatJetGoodNMuon2",bins=10, start=0, stop=10),
        **count_hist(coll="FatJetGoodNMuonSJ1",bins=10, start=0, stop=10),
        **count_hist(coll="FatJetGoodNMuonSJUnique1",bins=10, start=0, stop=10),
    },

    "columns" : {}

}

collections = ["FatJetGoodNMuon1", "FatJetGoodNMuon2", "FatJetGoodNMuonSJ1", "FatJetGoodNMuonSJUnique1"]

for coll in collections:
    cfg["variables"].update(**fatjet_hists(coll=coll))
    cfg["variables"][f"{coll}_pt"] = HistConf([Axis(name=f"{coll}_pt", coll=coll, field="pt",
                                                    label=r"FatJet $p_{T}$ [GeV]", bins=list(range(350, 1010, 10)))])
    cfg["variables"][f"{coll}_msoftdrop"] = HistConf([Axis(name=f"{coll}_msoftdrop", coll=coll, field="msoftdrop",
                                                           label=r"FatJet $m_{SD}$ [GeV]", bins=list(range(40, 410, 10)))])
    cfg["variables"][f"{coll}_pt_eta"] = HistConf(
        [ Axis(name=f"{coll}_pos", coll=coll, field="pos", type="int", label=r"FatJet position", bins=2, start=0, stop=2),
          Axis(name=f"{coll}_pt", coll=coll, field="pt", type="variable", label=r"FatJet $p_{T}$ [GeV]",
               bins=[350., 400., 450., 500., 550., 600., 700., 800., 900., 2500.]),
          Axis(name=f"{coll}_eta", coll=coll, field="eta", type="variable", label=r"FatJet $\eta$",
               bins=[-5, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 5]) ]
    )
    cfg["variables"][f"{coll}_pt_eta_tau21"] = HistConf(
        [ Axis(name=f"{coll}_pos", coll=coll, field="pos", type="int", label=r"FatJet position", bins=2, start=0, stop=2),
          Axis(name=f"{coll}_pt", coll=coll, field="pt", type="variable", label=r"FatJet $p_{T}$ [GeV]",
               bins=[350., 400., 450., 500., 550., 600., 700., 800., 900., 2500.]),
          Axis(name=f"{coll}_eta", coll=coll, field="eta", type="variable", label=r"FatJet $\eta$",
               bins=[-5, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 5]),
          Axis(name=f"{coll}_tau21", coll=coll, field="tau21", type="variable", label=r"FatJet $\tau_{21}$",
               bins=[0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1]) ]
    )
    cfg["variables"][f"{coll}_pt_eta_tau21_bintau05"] = HistConf(
        [ Axis(name=f"{coll}_pos", coll=coll, field="pos", type="int", label=r"FatJet position", bins=2, start=0, stop=2),
          Axis(name=f"{coll}_pt", coll=coll, field="pt", type="variable", label=r"FatJet $p_{T}$ [GeV]",
               bins=[350., 400., 450., 500., 550., 600., 700., 800., 900., 2500.]),
          Axis(name=f"{coll}_eta", coll=coll, field="eta", type="variable", label=r"FatJet $\eta$",
               bins=[-5, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 5]),
          Axis(name=f"{coll}_tau21", coll=coll, field="tau21", type="variable", label=r"FatJet $\tau_{21}$",
               bins=[0, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 1]) ]
    )
