cfg =  {
    # Dataset parameters
    "dataset"  : "datasets/DAS/RunIISummer20UL17.txt",
    "json"     : "datasets/RunIISummer20UL17.json",
    "storage_prefix" : "/pnfs/psi.ch/cms/trivcat/store/user/mmarcheg/BTVNanoCommissioning",
    "campaign" : "UL",
    "year"     : "2017",

    # PU files https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/LUMI_puWeights_Run2_UL/
    'puFile'  : '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2017_UL/puWeights.json.gz',
    'puJSON'   : 'Collisions17_UltraLegacy_goldenJSON',
    'nTrueFile' : '',

    # JEC
    "JECfolder": "correction_files/tmp",

    # Input and output files
    "workflow" : "fatjet_tagger",
    "input"    : "datasets/RunIISummer20UL17.json",
    "output"   : "histograms/RunIISummer20UL17.coffea",
    "plots"    : "plots/test",

    # Executor parameters
    "run_options" : {
        "executor"     : "futures",
        "workers"      : 12,
        "scaleout"     : 10,
        "chunk"        : 50000,
        "max"          : None,
        "skipbadfiles" : None,
        "voms"         : None,
        "limit"        : 1,
    },

    # Processor parameters
    "checkOverlap" : False,
    "hist2d"       : False,
    "mupt"         : 5,
}
