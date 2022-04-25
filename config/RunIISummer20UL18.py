cfg =  {
    # Dataset parameters
    "dataset"  : "datasets/DAS/RunIISummer20UL18.txt",
    "json"     : "datasets/RunIISummer20UL18.json",
    "storage_prefix" : "/pnfs/psi.ch/cms/trivcat/store/user/mmarcheg/BTVNanoCommissioning",
    "campaign" : "UL",
    "year"     : "2018",

    # PU files https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/LUMI_puWeights_Run2_UL/
    'puFile'  : '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2018_UL/puWeights.json.gz',
    'puJSON'   : 'Collisions18_UltraLegacy_goldenJSON',
    'nTrueFile' : '',   # for backcompatibility with EOY

    # JEC
    "JECfolder": "correction_files/tmp",

    # Input and output files
    "workflow" : "fatjet_tagger",
    "input"    : "datasets/RunIISummer20UL18_local.json",
    "output"   : "histograms/RunIISummer20UL18_FIXED_SV.coffea",
    "plots"    : "plots/RunIISummer20UL18_FIXED_SV",

    # Executor parameters
    "run_options" : {
        "executor"       : "parsl/slurm",
        "workers"        : 12,
        "scaleout"       : 10,
        "partition"      : "standard",
        "walltime"       : "12:00:00",
        "mem_per_worker" : None, # GB
        "exclusive"      : True,
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
