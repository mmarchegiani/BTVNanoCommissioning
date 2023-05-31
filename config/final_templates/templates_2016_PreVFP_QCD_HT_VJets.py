from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.mutag_oneMuAK8_processor import mutagAnalysisOneMuonInAK8Processor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
from pocket_coffea.lib.categorization import CartesianSelection, MultiCut
from config.fatjet_base.custom.cuts import get_nObj_minmsd, get_flavor, get_ptbin
from config.fatjet_base.custom.functions import get_HLTsel, get_inclusive_wp
from parameters import PtBinning, AK8TaggerWP, AK8Taggers

PtBinning = PtBinning['UL']['2016_PreVFP']
wps = AK8TaggerWP['UL']['2016_PreVFP']
pt_min = 350.
msd = 40.

common_cats = {
    "inclusive" : [passthrough],
}

cuts_pt = []
cuts_names_pt = []
for pt_low, pt_high in PtBinning.values():
    cuts_pt.append(get_ptbin(pt_low, pt_high))
    cuts_names_pt.append(f'Pt-{pt_low}to{pt_high}')
cuts_tagger = []
cuts_names_tagger = []
for tagger in AK8Taggers:
    for wp in ["L", "M", "H"]:
        for region in ["pass", "fail"]:
            cuts_tagger.append(get_inclusive_wp(tagger, wps[tagger][wp], region))
            cuts_names_tagger.append(f"msd{int(msd)}{tagger}{region}{wp}wp")

multicuts = [
    MultiCut(name="tagger",
             cuts=cuts_tagger,
             cuts_names=cuts_names_tagger),
    MultiCut(name="pt",
             cuts=cuts_pt,
             cuts_names=cuts_names_pt),
]

samples = ["QCD_HT",
           "VJets",
           ]

subsamples = {}
for s in filter(lambda x: 'DATA' not in x, samples):
    subsamples[s] = {f"{s}_{f}" : [get_flavor(f)] for f in ['l', 'c', 'b', 'cc', 'bb']}

cfg =  {
    "dataset" : {
        "jsons": ["datasets/MC_QCD_HT_RunIISummer20UL.json",
                  "datasets/MC_VJets_RunIISummer20UL.json",
                  ],
        "filter" : {
            "samples": samples,
            "samples_exclude" : [],
            "year": ['2016_PreVFP']
        },
        "subsamples" : subsamples
    },

    # Input and output files
    "workflow" : mutagAnalysisOneMuonInAK8Processor,
    "output"   : "output/pocket_coffea/final_templates/templates_2016_PreVFP_QCD_HT_VJets_withJES",
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
        "scaleout"       : 250,
        "queue"          : "standard",
        "walltime"       : "12:00:00",
        "mem_per_worker" : "12GB", # GB
        "exclusive"      : False,
        "chunk"          : 400000,
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
    "categories": CartesianSelection(multicuts=multicuts, common_cats=common_cats),

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
        "shape": {
            "common":{
                "inclusive": [ "JES_Total", "JER" ]
            }
        }
    },

   "variables":
    {
        **fatjet_hists(coll="FatJetGood"),
        "FatJetGood_pt" : HistConf([Axis(name=f"FatJetGood_pt", coll="FatJetGood", field="pt",
                                         label=r"FatJet $p_{T}$ [GeV]", bins=list(range(int(pt_min), 1010, 10)))]),
        "FatJetGood_msoftdrop" : HistConf([Axis(name=f"FatJetGood_msoftdrop", coll="FatJetGood", field="msoftdrop",
                                         label=r"FatJet $m_{SD}$ [GeV]", bins=list(range(int(msd), 410, 10)))]),
        "FatJetGood_logsumcorrSVmass": HistConf(
            [ Axis(coll="FatJetGood", field="logsumcorrSVmass", label=r"log($\sum({m^{corr}_{SV}})$)", bins=42, start=-2.4, stop=6) ]
        ),
        "nSVMatchedToFatJetGood": HistConf(
            [ Axis(coll="FatJetGood", field="nSVMatchedToFatJetGood", label="Number of SV inside AK8 jet", bins=10, start=0, stop=10) ]
        ),
        "FatJetGood_logsumcorrSVmass_tau21": HistConf(
            [ Axis(coll="FatJetGood", field="logsumcorrSVmass", label=r"log($\sum({m^{corr}_{SV}})$)", bins=42, start=-2.4, stop=6),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", type="variable", bins=[0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1]) ]
        ),
        "FatJetGood_pt_tau21": HistConf(
            [ Axis(coll="FatJetGood", field="pt", label=r"FatJet $p_{T}$ [GeV]", bins=list(range(int(pt_min), 1010, 10))),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", type="variable", bins=[0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1]) ]
        ),
        "FatJetGood_msoftdrop_tau21": HistConf(
            [ Axis(coll="FatJetGood", field="msoftdrop", label=r"FatJet $m_{SD}$ [GeV]", bins=list(range(int(msd), 410, 10))),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", type="variable", bins=[0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1]) ]
        ),
        "FatJetGood_msoftdrop_pt": HistConf(
            [ Axis(coll="FatJetGood", field="msoftdrop", label=r"FatJet $m_{SD}$ [GeV]", bins=list(range(int(msd), 410, 10))),
              Axis(coll="FatJetGood", field="pt", label=r"FatJet $p_{T}$ [GeV]", bins=list(range(int(pt_min), 1010, 10))) ]
        ),
    },

    "columns" : {}

}

# We reweigh histograms differently depending whether:
# - The histogram contains the whole jet collection ("all")
cfg["workflow_options"]["histograms_to_reweigh"]["by_pos"]["all"] = [name for name in cfg["variables"].keys() if name.startswith("FatJetGood_") and not name.endswith(("_1", "_2"))]
