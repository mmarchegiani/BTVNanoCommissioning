import os
import sys
import re
from copy import deepcopy
import argparse
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
parser.add_argument('--xlim', type=float, nargs=2, default=(-1.2, 5.2), required=False)
parser.add_argument('--adapt_light', action="store_true", required=False)
parser.add_argument('--step_light', type=float, default=0.1, required=False)
parser.add_argument('--retries', type=int, default=8, required=False)
parser.add_argument('--freeze_light_extra', action="store_true", required=False)
parser.add_argument('--threshold_light', type=float, default=0.05, required=False)
parser.add_argument('--freeze_frac_l', action="store_true", required=False)
parser.add_argument('--freeze_bkg', action="store_true", required=False)
parser.add_argument('--threshold_bkg', type=float, default=1.0, required=False)
parser.add_argument('--threshold_chi2', type=float, default=1.5, required=False)

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
parameters[args.year] = {}
for tagger in AK8Taggers:
    for wp in wps:
        for pt in pts:
            cat = "msd40{}{}wp_Pt-{}".format(tagger, wp, pt)
            parameters[args.year][cat] = {}
            for poi in pois:
                if poi in parameter_ranges:
                    parameters[args.year][cat][poi] = parameter_ranges[poi]
                else:
                    parameters[args.year][cat][poi] = parameter_ranges_default[poi]

"""
if args.year == "2016_PreVFP":
    parameters[args.year]["msd40btagDDBvLV2Hwp_Pt-450to500"]["c_cc"] = {"value" : 1, "lo" : 1, "hi" : 1}
    parameters[args.year]["msd40particleNetMD_Xbb_QCDLwp_Pt-600toInf"]["c_cc"] = {"value" : 1, "lo" : 1, "hi" : 1}
"""
if args.year == "2016_PostVFP":
    if ("fit_tau21p25_3poi_range_bb05_cc02_l05" in args.output) & ("reweighted" not in args.output):
        parameters[args.year]["msd40deepTagMD_ZHccvsQCDLwp_Pt-600toInf"]["b_bb"] = {"value" : 1, "lo" : 1, "hi" : 1}
    if ("fit_tau21p30_3poi_range_bb05_cc02_l05" in args.output) & ("reweighted" not in args.output):
        parameters[args.year]["msd40deepTagMD_ZHccvsQCDLwp_Pt-600toInf"]["b_bb"] = {"value" : 1, "lo" : 1, "hi" : 1}
    if ("fit_tau21p30_3poi_range_bb02_cc05_l05" in args.output) & ("reweighted" not in args.output):
        parameters[args.year]["msd40btagHbbLwp_Pt-450to500"]["c_cc"] = {"value" : 1, "lo" : 1, "hi" : 1}
    if ("fit_tau21p35_3poi_range_bb05_cc02_l05" in args.output) & ("reweighted" not in args.output):
        parameters[args.year]["msd40deepTagMD_ZHccvsQCDLwp_Pt-600toInf"]["b_bb"] = {"value" : 1, "lo" : 1, "hi" : 1}
    if ("fit_tau21p40_3poi_range_bb05_cc02_l05" in args.output) & ("reweighted" not in args.output):
        parameters[args.year]["msd40deepTagMD_ZHccvsQCDLwp_Pt-600toInf"]["b_bb"] = {"value" : 1, "lo" : 1, "hi" : 1}
if args.year == "2018":
    if ("fit_tau21p35_3poi_range_bb02_cc05_l05" in args.output) & ("reweighted" not in args.output):
        parameters[args.year]["msd40particleNetMD_Xbb_QCDLwp_Pt-450to500"]["c_cc"] = {"value" : 1, "lo" : 1, "hi" : 1}
        parameters[args.year]["msd40btagHbbHwp_Pt-450to500"]["c_cc"] = {"value" : 1, "lo" : 1, "hi" : 1}

fit = Fit(args.input,
          args.output,
          categories,
          var,
          args.year,
          xlim=args.xlim,
          binwidth=args.binwidth,
          regions=regions,
          scheme=args.scheme,
          npoi=args.npoi,
          frac_effect=args.frac,
          parameters=parameters,
          freeze_light=True,
          threshold_light=args.threshold_light,
          freeze_frac_l=args.freeze_frac_l,
          freeze_bkg=args.freeze_bkg,
          threshold_bkg=args.threshold_bkg,
          threshold_chi2=args.threshold_chi2
          )
if args.mode == "all":
    fit.run_fits("FitDiagnostics")
    fit.run_fits("MultiDimFit")
elif args.mode == "FitDiagnostics":
    fit.run_fits(args.mode, job=not args.no_jobs)
elif args.mode == "MultiDimFit":
    fit.run_fits(args.mode, job=not args.no_jobs)

if args.mode == "FitDiagnostics":
    first_bb = True
    first_cc = True
    folders_corrupted = []

    # Wait until all fits are done
    while not fit.all_done():
        continue

    idx_iteration = 0
    if args.adapt_light:
        categories_at_boundary = [cat for cat in categories if fit.is_at_boundary(cat.replace("pass", "").replace("fail", ""))]
        for idx_iteration in range(args.retries):
            l_lo = 0 + args.step_light*idx_iteration
            l_hi = 2 - args.step_light*idx_iteration
            print("Adapting light POI range: l={},{}".format(l_lo, l_hi))
            parameters_alternative = deepcopy(parameters)
            for year in parameters_alternative:
                for cat in categories_at_boundary:
                    model_name = cat.replace("pass", "").replace("fail", "")
                    parameters_alternative[year][model_name]['l']['lo'] = l_lo
                    parameters_alternative[year][model_name]['l']['hi'] = l_hi

            fit_alternative = Fit(args.input,
                                  os.path.join(args.output, "iteration_{}".format(idx_iteration)),
                                  categories_at_boundary,
                                  var,
                                  args.year,
                                  xlim=args.xlim,
                                  binwidth=args.binwidth,
                                  regions=regions,
                                  scheme=args.scheme,
                                  npoi=args.npoi,
                                  frac_effect=args.frac,
                                  parameters=parameters_alternative,
                                  freeze_frac_l=args.freeze_frac_l,
                                  freeze_bkg=args.freeze_bkg,
                                  threshold_bkg=args.threshold_bkg,
                                  threshold_chi2=args.threshold_chi2
                                  )
            fit_alternative.run_fits(args.mode, job=not args.no_jobs)
            if fit_alternative.any_at_boundary():
                categories_at_boundary = [cat for cat in categories_at_boundary if fit_alternative.is_at_boundary(cat.replace("pass", "").replace("fail", ""))]
            else:
                break
    if args.freeze_light_extra:
        categories_at_boundary = [cat for cat in categories if fit.is_at_boundary(cat.replace("pass", "").replace("fail", ""))]
        categories_failed = [cat for cat in categories if fit.is_failed(cat.replace("pass", "").replace("fail", ""))]
        if len(categories_at_boundary + categories_failed) > 0:
            fit_alternative = Fit(args.input,
                                os.path.join(args.output, "freeze_light"),
                                categories_at_boundary + categories_failed,
                                var,
                                args.year,
                                xlim=args.xlim,
                                binwidth=args.binwidth,
                                regions=regions,
                                scheme=args.scheme,
                                npoi=args.npoi,
                                frac_effect=args.frac,
                                parameters=parameters,
                                freeze_light=True,
                                freeze_frac_l=args.freeze_frac_l,
                                freeze_bkg=args.freeze_bkg,
                                threshold_bkg=args.threshold_bkg,
                                threshold_chi2=args.threshold_chi2
                                )
            fit_alternative.run_fits(args.mode, job=not args.no_jobs)

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
        print("Results saved in file {}".format(file_results_all))

    print("Corrupted folders:")
    for folder in folders_corrupted:
        print(folder)

    if args.adapt_light:
        first_bb = True
        first_cc = True
        folders_corrupted = []
        for model_name, folder in fit.fitdirs.items():
            file_results = os.path.join(folder, "fitResults.csv")

            # Get results from the fitResults.csv file in the latest iteration folder
            indices = range(args.retries)
            indices.reverse()
            for idx_iteration in indices:
                folder_alternative = folder.replace("fitdir", "iteration_{}/fitdir".format(idx_iteration))
                if os.path.exists(folder_alternative):
                    print("Reading results in alternative folder {}".format(folder_alternative))
                    file_results = os.path.join(folder_alternative, "fitResults.csv")
                    break

            if not os.path.exists(file_results):
                folders_corrupted.append(folder)
                continue
            print("Reading results in file {}".format(file_results))
            df = pd.read_csv(file_results)
            kwargs = {'mode' : 'a', 'header' : False}
            if any(tagger in model_name for tagger in AK8Taggers_bb):
                file_results_alternative = os.path.join(args.output, "fitResults_alternative_bb.csv")
                if first_bb:
                    kwargs = {'mode' : 'w', 'header' : True}
                    first_bb = False
            elif any(tagger in model_name for tagger in AK8Taggers_cc):
                file_results_alternative = os.path.join(args.output, "fitResults_alternative_cc.csv")
                if first_cc:
                    kwargs = {'mode' : 'w', 'header' : True}
                    first_cc = False
            df.to_csv(file_results_alternative, **kwargs)
            print("Results saved in file {}".format(file_results_alternative))

        print("Corrupted folders:")
        for folder in folders_corrupted:
            print(folder)

    if args.freeze_light_extra:
        first_bb = True
        first_cc = True
        folders_corrupted = []
        for model_name, folder in fit.fitdirs.items():
            file_results = os.path.join(folder, "fitResults.csv")

            # Get results from the fitResults.csv file in the latest iteration folder
            indices = range(args.retries)
            indices.reverse()
            folder_alternative = folder.replace("fitdir", "freeze_light/fitdir")
            if os.path.exists(folder_alternative):
                print("Reading results in alternative folder {}".format(folder_alternative))
                file_results = os.path.join(folder_alternative, "fitResults.csv")

            if not os.path.exists(file_results):
                folders_corrupted.append(folder)
                continue
            print("Reading results in file {}".format(file_results))
            df = pd.read_csv(file_results)
            columns_to_drop = ["corr_b_bb_l", "corr_c_cc_l"]
            columns_filtered = list(filter(lambda x : x not in columns_to_drop, df.columns))
            df = df[columns_filtered]
            kwargs = {'mode' : 'a', 'header' : False}
            if any(tagger in model_name for tagger in AK8Taggers_bb):
                file_results_alternative = os.path.join(args.output, "fitResults_freeze_light_bb.csv")
                if first_bb:
                    kwargs = {'mode' : 'w', 'header' : True}
                    first_bb = False
            elif any(tagger in model_name for tagger in AK8Taggers_cc):
                file_results_alternative = os.path.join(args.output, "fitResults_freeze_light_cc.csv")
                if first_cc:
                    kwargs = {'mode' : 'w', 'header' : True}
                    first_cc = False
            df.to_csv(file_results_alternative, **kwargs)
            print("Results saved in file {}".format(file_results_alternative))

        print("Corrupted folders:")
        for folder in folders_corrupted:
            print(folder)
