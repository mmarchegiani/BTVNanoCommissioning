from pocket_coffea.parameters.cuts.preselection_cuts import *
from workflows.genweights import genWeightsProcessor
from pocket_coffea.lib.cut_functions import get_nObj_min
from pocket_coffea.parameters.histograms import *
import numpy as np

samples = [
          "QCD_MuEnriched",
          "DATA"
           ]

cfg =  {
    "dataset" : {
        "jsons": ["datasets/MC_QCD_MuEnriched_RunIISummer20UL_redirector.json",
                  "datasets/DATA_BTagMu_RunIISummer20UL_redirector.json"],
        "filter" : {
            "samples": samples,
            "samples_exclude" : [],
            "year": ['2017']
        },
    },
    

    # Input and output files
    "workflow" : genWeightsProcessor,
    "output"   : "output/pocket_coffea/genweights/genweights_2017UL_NanoAODv9",
    "workflow_options" : {},

    "run_options" : {
        "executor"       : "dask/slurm",
        "workers"        : 1,
        "scaleout"       : 125,
        "queue"          : "short",
        "walltime"       : "1:00:00",
        "mem_per_worker" : "2GB", # GB
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
    "skim": [],

    "save_skimmed_files": None,
    
    "preselections" : [],
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
                "inclusive": [ "pileup", "sf_L1prefiring" ],
                "bycategory" : {
                }
            },
        "bysample": {
        }    
        },
        "shape": {
            "common":{
                "inclusive": [ "JES_Total" ]
            }
        }
    },

    "variables": {},

    "columns" : {}

}
