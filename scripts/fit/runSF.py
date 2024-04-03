import os
import sys
import re
from time import sleep
import argparse
import json
import pandas as pd

sys.path.append('/work/mmarcheg/BTV/BTVNanoCommissioning')

from utils.Fit import Fit
from parameters_fit import fit_parameters
from parameters import categories, AK8Taggers, AK8Taggers_bb, AK8Taggers_cc

def find_matching_strings(pattern, strings):
    matching_strings = []
    regex = re.compile(pattern)
    for string in strings:
        if regex.match(string):
            matching_strings.append(string)
    return matching_strings

parser = argparse.ArgumentParser(description='Save histograms in pickle format for combine fit')
parser.add_argument('-i', '--input', default=None, help='Input templates', required=True)
parser.add_argument('-o', '--output', default="/work/mmarcheg/BTVNanoCommissioning/output/fit", help='Output folder', required=True)
parser.add_argument('--only', type=str, default=None, nargs='+', help='Filter categories by list of regular expressions', required=False)
parser.add_argument('--scheme', type=str, choices=['3f', '5f'],  help='3-flavor scheme', required=True)
parser.add_argument('--npoi', type=int, choices=[1,3],  help='Number of POIs to be used in the fit. Available choices: 1 POI (only dominant flavour will be floating) or 3 POIs (all 3 POIs will be floating)', required=True)
parser.add_argument('--binwidth', type=float, default=0.2, choices=[0.1, 0.2, 0.4],  help='Specify the binwidth of the logsumcorrmass distribution', required=False)
parser.add_argument('--var', type=str, default="FatJetGood_logsumcorrSVmass", choices=["FatJetGood_logsumcorrSVmass", "FatJetGood_logsv1mass"], help='Specify the variable to be used for the fit', required=False)
parser.add_argument('--year', type=str, choices=["2016_PreVFP", "2016_PostVFP", "2017", "2018"], help='Specify the data-taking year', required=True)
parser.add_argument('--frac', type=float, default=1.2, help='Relative variation of the frac_* parameters', required=False)
parser.add_argument('--setParameterRanges', type=str, default=None, help='Set the range of the POIs with a string of the form "POI1:min,max:POI2:min,max:POI3:min,max"', required=False)
parser.add_argument('-m', '--mode', type=str, default="FitDiagnostics", choices=["FitDiagnostics", "MultiDimFit", "all"], help='Specify combine mode', required=False)
parser.add_argument('--dim', type=int, default=1, required=False)
parser.add_argument('--passonly', action="store_true", required=False)
parser.add_argument('--no-jobs', action="store_true", required=False)
parser.add_argument('--xlim', type=float, nargs=2, default=(-2.4, 6.0), required=False)

args = parser.parse_args()

filtered_categories = []
if args.only:
    for key in args.only:
        filtered_categories += find_matching_strings(key, categories)
    categories = filtered_categories

if os.path.exists(args.output):
    sys.exit("The output folder {} is already existing. Please choose a different folder name.".format(args.output))

failed_fits = []
if args.dim == 1:
    var = args.var
elif args.dim == 2:
    var = args.var + '_tau21'
else:
    raise NotImplementedError

parameter_ranges_default = {
    "l"    : {"value" : 1, "lo" : 0.5, "hi" : 2},
    "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
    "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
}
parameters = {}
if not args.setParameterRanges:
    parameter_ranges = parameter_ranges_default
else:
    # Get the dictionary of parameters from the string
    parameter_ranges = {}
    for par in args.setParameterRanges.split(':'):
        poi, _range = par.split('=')
        lo, hi = (float(x) for x in _range.split(','))
        parameter_ranges[poi] = {'value' : 1.0, 'lo' : lo, 'hi' : hi}

if args.passonly:
    regions = ['pass']
else:
    regions = ['pass', 'fail']

if args.scheme == "3f":
    pois = ['l', 'b_bb', 'c_cc']
elif args.scheme == "5f":
    pois = ['l', 'c', 'b', 'cc', 'bb']
wps = ['L', 'M', 'H']
pts = ['450to500', '500to600', '600toInf']
# Set the POIs range in the dictionary
for year, pars in fit_parameters.items():
    parameters[year] = {}
    for tagger in AK8Taggers:
        for wp in wps:
            for pt in pts:
                cat = "msd40{}{}wp_Pt-{}".format(tagger, wp, pt)
                parameters[year][cat] = {}
                for poi in pois:
                    if poi in parameter_ranges:
                        parameters[year][cat][poi] = parameter_ranges[poi]
                    else:
                        parameters[year][cat][poi] = parameter_ranges_default[poi]

"""
print("Parameters:")
print(json.dumps(parameters, indent=4))
sys.exit()
"""

fit = Fit(args.input, args.output, categories, var, args.year, xlim=args.xlim, binwidth=args.binwidth, regions=regions, scheme=args.scheme, npoi=args.npoi, frac_effect=args.frac, parameters=parameters)
if args.mode == "all":
    fit.run_fits("FitDiagnostics")
    fit.run_fits("MultiDimFit")
elif args.mode == "FitDiagnostics":
    fit.run_fits(args.mode, job=not args.no_jobs)
    fit.save_results(args.mode)
elif args.mode == "MultiDimFit":
    fit.run_fits(args.mode, job=not args.no_jobs)

if args.mode == "FitDiagnostics":
    first_bb = True
    first_cc = True
    folders_corrupted = []

    # Wait until all fits are done
    while not fit.all_done():
        continue

    for model_name, folder in fit.fitdirs.items():
        file_results = os.path.join(folder, "fitResults.csv")
        if not os.path.exists(file_results):
            folders_corrupted.append(folder)
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

    print("Corrupted folders:")
    for folder in folders_corrupted:
        print(folder)
