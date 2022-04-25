cfg =  {
    # Dataset parameters
    "dataset"  : "datasets/DAS/datasets_btag2016ULPreVFP.txt",
    "json"     : "datasets/RunIISummer20UL16-PreVFP.json",
    "storage_prefix" : "/pnfs/psi.ch/cms/trivcat/store/user/mmarcheg/BTVNanoCommissioning",
    "campaign" : "UL",
    "year"     : "2016",

    # PU files https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/LUMI_puWeights_Run2_UL/
    'puFile'  : '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2016preVFP_UL/puWeights.json.gz',
    'puJSON'   : 'Collisions16_UltraLegacy_goldenJSON',

    # JEC
    "JECfolder": "correction_files/tmp",

    # Input and output files
    "workflow" : "fatjet_tagger",
    "input"    : "datasets/RunIISummer20UL16-PreVFP_local.json",
    "output"   : "histograms/RunIISummer20UL16-PreVFP.coffea",
    "plots"    : "plots/test",

    # Executor parameters
    "run_options" : {
        "executor"       : "parsl/slurm",
        "workers"        : 12,
        "scaleout"       : 50,
        "partition"      : "standard",
        "walltime"       : "6:00:00",
        "mem_per_worker" : 2, # GB
        "exclusive"      : False,
        "chunk"          : 50000,
        "max"            : None,
        "skipbadfiles"   : None,
        "voms"           : None,
        "limit"          : None,
    },

    # Processor parameters
    "checkOverlap" : False,
    "hist2d"       : False,
    "mupt"         : 5,
}
