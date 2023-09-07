import os
import argparse
import json
import numpy as np
from coffea.util import load

from pocket_coffea.utils.dataset import Dataset

from config.fatjet_base.custom.parameters.genweights.genweights import genweights_files

parser = argparse.ArgumentParser(description='Append genweights to the datasets_definition_skim.json file where the skimmed datasets are defined.')
parser.add_argument(
    '--cfg',
    default=os.getcwd() + "/datasets/skim/datasets_definition_skim.json",
    help='Config file with datasets parameters',
    required=False,
)
parser.add_argument(
    '-o',
    '--overwrite',
    action='store_true',
    help="Overwrite existing file definition json",
    default=False,
)
parser.add_argument(
    '-gf',
    '--genweights-file',
    help="Coffea file with sum_genweights",
    required=True
)
parser.add_argument(
    '-y',
    '--year',
    help="Data-taking year of dataset to update",
    required=True
)
args = parser.parse_args()
config = json.load(open(args.cfg))
genweights_dict = load(args.genweights_file)['sum_genweights']

#breakpoint()

for dataset, dict_dataset in config.items():
    sample = dict_dataset["sample"]
    if sample != "QCD_MuEnriched": continue
    for file in dict_dataset["files"]:
        year = file["metadata"]["year"]
        isMC = file["metadata"]["isMC"]
        if year != args.year: continue
        if isMC:
            try:
                part = file["metadata"]["part"]
            except:
                breakpoint()
            key = f"{sample}_{part}"
            if key in genweights_dict:
                file["metadata"]["sum_genweights"] = np.float64(genweights_dict[key])
            else:
                raise Exception(f"The metadata dict of {dataset} has no key named `sum_genweights`.")

if args.overwrite:
    with open(args.cfg, 'w') as fp:
        json.dump(config, fp, indent=4)
        fp.close()
