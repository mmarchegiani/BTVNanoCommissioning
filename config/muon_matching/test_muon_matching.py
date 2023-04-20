from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.fatjet_base import fatjetBaseProcessor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
from config.fatjet_base.custom.cuts import twojets_presel, mutag_fatjet_sel, mutag_subjet_sel, get_ptmsd, get_nObj_minmsd, get_flavor
from config.fatjet_base.custom.functions import get_HLTsel

msd = 40.

common_cats = {
    "inclusive" : [passthrough],
    "pt350msd40_mutag_fatjet_nmu-1" : [get_ptmsd(350., msd), mutag_fatjet_sel(nmu=1)],
    "pt450msd40_mutag_fatjet_nmu-1" : [get_ptmsd(450., msd), mutag_fatjet_sel(nmu=1)],
    "pt350msd40_mutag_fatjet_nmu-2" : [get_ptmsd(350., msd), mutag_fatjet_sel(nmu=2)],
    "pt450msd40_mutag_fatjet_nmu-2" : [get_ptmsd(450., msd), mutag_fatjet_sel(nmu=2)],
    "pt350msd40_mutag_subjet" : [get_ptmsd(350., msd), mutag_subjet_sel(unique_matching=False)],
    "pt450msd40_mutag_subjet" : [get_ptmsd(450., msd), mutag_subjet_sel(unique_matching=False)],
    "pt350msd40_mutag_subjet_unique" : [get_ptmsd(350., msd), mutag_subjet_sel(unique_matching=True)],
    "pt450msd40_mutag_subjet_unique" : [get_ptmsd(450., msd), mutag_subjet_sel(unique_matching=True)],
}

samples = ["QCD_MuEnriched"]
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
    "output"   : "output/test/muon_matching",
    "workflow_options" : {
        "histograms_to_reweigh" : [],
        "reweighting_scheme"    : "ptetatau21"},

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
             get_nObj_min(2, 3., "Muon"),
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
                "inclusive": [ "pileup" ],#, "sf_L1prefiring" ],
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
        #**muon_hists(coll="MuonGoodMatchedToFatJetGood", pos=0),
        #**muon_hists(coll="MuonGoodMatchedToFatJetGood", pos=1),
        #**muon_hists(coll="MuonGoodMatchedToLeadingSubJet"),
        #**muon_hists(coll="MuonGoodMatchedToSubleadingSubJet"),
        #**muon_hists(coll="MuonGoodMatchedUniquelyToLeadingSubJet"),
        #**muon_hists(coll="MuonGoodMatchedUniquelyToSubleadingSubJet"),
        **fatjet_hists(coll="FatJetGood"),
        **sv_hists(coll="events"),
        "nMuonGoodMatchedToFatJetGood": HistConf(
            [ Axis(coll="FatJetGood", field="nMuonGoodMatchedToFatJetGood", label=r"$N_{\mu, AK8}$", bins=4, start=0, stop=4) ]
        ),
        "nMuonGoodMatchedToLeadingSubJet": HistConf(
            [ Axis(coll="FatJetGood", field="nMuonGoodMatchedToLeadingSubJet", label=r"$N_{\mu, sj1}$", bins=4, start=0, stop=4) ]
        ),
        "nMuonGoodMatchedToSubleadingSubJet": HistConf(
            [ Axis(coll="FatJetGood", field="nMuonGoodMatchedToSubleadingSubJet", label=r"$N_{\mu, sj2}$", bins=4, start=0, stop=4) ]
        ),
        "nMuonGoodMatchedUniquelyToLeadingSubJet": HistConf(
            [ Axis(coll="FatJetGood", field="nMuonGoodMatchedUniquelyToLeadingSubJet", label=r"$N_{\mu, sj1}$", bins=4, start=0, stop=4) ]
        ),
        "nMuonGoodMatchedUniquelyToSubleadingSubJet": HistConf(
            [ Axis(coll="FatJetGood", field="nMuonGoodMatchedUniquelyToSubleadingSubJet", label=r"$N_{\mu, sj2}$", bins=4, start=0, stop=4) ]
        ),
        #"dimuon_pt": HistConf(
        #    [ Axis(coll="dimuon", field="pt", label=r"Di-muon $p_{T}$ [GeV]", bins=50, start=0, stop=200) ]
        #),
        **count_hist(name="nMuonGood", coll="MuonGood",bins=5, start=0, stop=5),
        **count_hist(name="nFatJets", coll="FatJetGood",bins=5, start=0, stop=5),
        **count_hist(name="nSV", coll="SV",bins=10, start=0, stop=10),
        "FatJetGood_tau21_logsumcorrmass": HistConf(
            [ Axis(coll="FatJetGood", field="logsumcorrmass", label=r"log($\sum({m^{corr}_{SV}})$)", bins=43, start=-2.6, stop=6),
              Axis(coll="FatJetGood", field="tau21", label=r"$\tau_{21}$", bins=20, start=0, stop=1) ]
        ),
    },

    "columns" : {}

}

# Here we update the dictionary of the histograms to reweigh such that
# only the histograms with the same shape as the FatJetGood collection are reweighed
cfg["workflow_options"]["histograms_to_reweigh"] = list(cfg["variables"].keys())
