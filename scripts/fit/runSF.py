import os
import sys
import argparse
import json
import pandas as pd

sys.path.append('/work/mmarcheg/BTVNanoCommissioning')

from utils.Fit import Fit
from parameters import categories, AK8Taggers_bb, AK8Taggers_cc

parser = argparse.ArgumentParser(description='Save histograms in pickle format for combine fit')
parser.add_argument('-i', '--input', default=None, help='Input templates', required=True)
parser.add_argument('-o', '--output', default="/work/mmarcheg/BTVNanoCommissioning/output/fit", help='Output folder', required=True)
parser.add_argument('--only', type=str, default=None, help='Filter categories by key', required=False)
parser.add_argument('--scheme', type=str, choices=['3f', '5f'],  help='3-flavor scheme', required=True)
parser.add_argument('--binwidth', type=float, default=0.2, choices=[0.1, 0.2, 0.4],  help='Specify the binwidth of the logsumcorrmass distribution', required=False)
parser.add_argument('--year', type=str, choices=["2016_PreVFP", "2016_PostVFP", "2017", "2018"], help='Specify the data-taking year', required=True)
parser.add_argument('-m', '--mode', type=str, default="FitDiagnostics", choices=["FitDiagnostics", "MultiDimFit", "all"], help='Specify combine mode', required=False)
parser.add_argument('--no-jobs', action="store_true", required=False)
args = parser.parse_args()

if args.only:
    args.only = args.only.split('*')
    categories = list(filter(lambda cat : all(f in cat for f in args.only), categories))

if os.path.exists(args.output):
    sys.exit("The output folder {} is already existing. Please choose a different folder name.".format(args.output))

failed_fits = []
for var in ['FatJetGood_logsumcorrSVmass']:
    fit = Fit(args.input, args.output, categories, var, args.year, xlim=(-2.4, 6.0), binwidth=args.binwidth, scheme=args.scheme)
    if args.mode == "all":
        fit.run_fits("FitDiagnostics")
        fit.run_fits("MultiDimFit")
    else:
        fit.run_fits(args.mode, job=not args.no_jobs)
        fit.save_results(args.mode)
    first_bb = True
    first_cc = True
    for model_name, folder in fit.fitdirs.items():
        file_results = os.path.join(folder, "fitResults.csv")
        if not os.path.exists(file_results):
            continue
        df = pd.read_csv(file_results)
        kwargs = {'mode' : 'a', 'header' : False}
        if any(tagger in model_name for tagger in AK8Taggers_bb):
            file_results_all = os.path.join(args.output, "fitResults_bb.csv")
            if first_bb:
                kwargs = {'mode' : 'w', 'header' : True}
                first_bb = False
        elif any(tagger in model_name for tagger in AK8Taggers_cc):
            file_results_all = os.path.join(args.output, "fitResults_cc.csv")
            if first_cc:
                kwargs = {'mode' : 'w', 'header' : True}
                first_cc = False
        df.to_csv(file_results_all, **kwargs)
