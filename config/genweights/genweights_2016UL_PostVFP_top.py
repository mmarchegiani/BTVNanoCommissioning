from pocket_coffea.parameters.cuts.preselection_cuts import passthrough
from config.genweights.workflow import genWeightsProcessor
from pocket_coffea.lib.cut_functions import get_nObj_min
from config.fatjet_base.custom.cuts import get_nObj_minmsd
from config.fatjet_base.custom.functions import get_HLTsel

samples = ["SingleTop_ttbar"]

cfg =  {
    "dataset" : {
        "jsons": ["datasets/MC_top_RunIISummer20UL_local.json"],
        "filter" : {
            "samples": samples,
            "samples_exclude" : [],
            "year": ['2016_PostVFP']
        }
    },

    # Input and output files
    "workflow" : genWeightsProcessor,
    "output"   : "output/pocket_coffea/genweights/genweights_2016UL_PostVFP_top",
    "workflow_options" : {},

    "run_options" : {
        "executor"       : "dask/slurm",
        "workers"        : 1,
        "scaleout"       : 125,
        "queue"          : "standard",
        "walltime"       : "8:00:00",
        "mem_per_worker" : "4GB", # GB
        "exclusive"      : False,
        "chunk"          : 1000000,
        "retries"        : 50,
        "treereduction"  : 20,
        "max"            : None,
        "skipbadfiles"   : None,
        "voms"           : None,
        "limit"          : None,
        "adapt"          : False,
        "env"            : "conda"
    },


    # Cuts and plots settings
    "finalstate" : "mutag",
    "skim": [get_nObj_min(1, 200., "FatJet"),
             get_nObj_minmsd(1, 30., "FatJet"),
             get_nObj_min(1, 3., "Muon"),
             get_HLTsel("mutag")],

    "save_skimmed_files": None,
    
    "preselections" : [passthrough],
    "categories": {
        "inclusive" : [passthrough]
    },

    "weights": {
        "common": {
            "inclusive": [ ],
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
        "shape": {
            "common":{
                "inclusive": [ ]
            }
        }
    },

    "variables": {},

    "columns" : {}

}
