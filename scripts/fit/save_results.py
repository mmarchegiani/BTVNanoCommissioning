import os
import sys
import argparse
import pandas as pd
import uproot
import ROOT

AK8Taggers_bb = ['btagDDBvLV2', 'particleNetMD_Xbb_QCD', 'deepTagMD_ZHbbvsQCD', 'btagHbb']
AK8Taggers_cc = ['btagDDCvLV2', 'particleNetMD_Xcc_QCD', 'deepTagMD_ZHccvsQCD']

sys.path.append('/work/mmarcheg/BTVNanoCommissioning')

def get_signal_name(model_name):
    for tagger in AK8Taggers_bb:
        if model_name.startswith("msd40{}".format(tagger)):
            return 'b+bb', tagger
    for tagger in AK8Taggers_cc:
        if model_name.startswith("msd40{}".format(tagger)):
            return 'c+cc', tagger

parser = argparse.ArgumentParser(description='Save histograms in pickle format for combine fit')
parser.add_argument('-i', '--input', default=None, help='Fit folder', required=True)
args = parser.parse_args()

if not os.path.exists(args.input):
    sys.exit("The input folder {} does not exist.".format(args.input))

parent_folder = os.path.join(args.input, "fitdir")
ls = os.listdir(parent_folder)
folder_results = os.path.join(args.input, "results")
if not os.path.exists(folder_results):
    os.makedirs(folder_results)

failed_fits = []
for model_name in ls:
    if model_name == "msd40deepTagMD_ZHbbvsQCDLwp_Pt-450to500_testPOISfrozen": continue
    fitdir = os.path.join(parent_folder, model_name)
    for file in os.listdir(fitdir):
        if file.startswith("higgsCombine"):
            try:
                combineFile = uproot.open(os.path.join(fitdir, file))
            except:
                raise Exception("Error in opening file {}".format(os.path.join(fitdir, file)))
            break
    if not 'limit;1' in combineFile.keys():
        failed_fits.append(model_name)
        continue
    combineTree = combineFile['limit']
    combineBranches = combineTree.arrays()
    results = combineBranches['limit']
    if len(results) < 4:
        failed_fits.append(model_name)
        continue
    combineCont, low, high, temp = results
    combineErrUp = high - combineCont
    combineErrDown = combineCont - low

    d = {}

    POI, tagger = get_signal_name(model_name)
    wp = model_name.split("wp")[0][-1]
    wpt = model_name.split("Pt-")[-1]
    columns = ['year', 'selection', 'wp', 'pt',
        POI, '{}ErrUp'.format(POI), '{}ErrDown'.format(POI),
        'SF({})'.format(POI)]
    columns_for_latex = ['year', 'pt', 'SF({})'.format(POI)]
    d = {'selection' : [model_name], 'tagger' : [tagger],
        'wp' : [wp], 'pt' : [wpt],
        POI : [combineCont], '{}ErrUp'.format(POI) : [combineErrUp], '{}ErrDown'.format(POI) : [combineErrDown],
        'SF({})'.format(POI) : ['{}$^{{+{}}}_{{-{}}}$'.format(combineCont, combineErrUp, combineErrDown)]}

    fitResults = ROOT.TFile.Open(os.path.join(fitdir, "fitDiagnostics_{}.root".format(model_name)))
    fit_s = fitResults.Get('fit_s')

    for flavor in ['b+bb', 'c+cc', 'l']:
        if (flavor == POI): continue
        par_result = fit_s.floatParsFinal().find(flavor)
        columns.append(flavor)
        columns.append('{}Err'.format(flavor))
        columns.append('SF({})'.format(flavor))
        columns_for_latex.append('SF({})'.format(flavor))
        if par_result == None:
            d.update({flavor : -999, '{}Err'.format(flavor) : -999, 'SF({})'.format(flavor) : r'{}$\pm${}'.format(-999, -999)})
            continue
        parVal = par_result.getVal()
        parErr = par_result.getAsymErrorHi()
        #columns.append(flavor)
        #columns.append('{}Err'.format(flavor))
        #columns.append('SF({})'.format(flavor))
        #columns_for_latex.append('SF({})'.format(flavor))
        d.update({flavor : parVal, '{}Err'.format(flavor) : parErr, 'SF({})'.format(flavor) : r'{}$\pm${}'.format(parVal, parErr)})
    df = pd.DataFrame(data=d)
    filename = os.path.join(folder_results, "fitResults_{}.csv".format(model_name))
    print("Saving fit output in {}".format(filename))
    df.to_csv(filename, columns=columns, mode='w', header=True)

with open(os.path.join(folder_results, "FAILED_FITS.txt"), 'w') as f:
    for model_name in failed_fits:
        f.write("{}\n".format(model_name))
print("Failed fits:")
print(failed_fits)
