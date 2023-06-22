import sys, os
import argparse
import numpy as np
from coffea.util import load
import hist
import pickle

def dense_axes(h):
    '''Returns the list of dense axes of a histogram.'''
    dense_axes = []
    if type(h) == dict:
        h = h[list(h.keys())[0]]
    for ax in h.axes:
        if not type(ax) in [hist.axis.StrCategory, hist.axis.IntCategory]:
            dense_axes.append(ax)
    return dense_axes

def get_axis_items(h, axis_name):
    axis = h.axes[axis_name]
    return list(axis.value(range(axis.size)))

parser = argparse.ArgumentParser(description='Save histograms in pickle format for combine fit')
parser.add_argument('-i', '--input', type=str, help='Input file with histograms', required=True)
parser.add_argument('-o', '--output', type=str, default=None, help='Output folder to save templates', required=False)
parser.add_argument('-l', '--label', type=str, default=None, help='Extra label for output file', required=False)

args = parser.parse_args()

if os.path.isfile( args.input ): accumulator = load(args.input)
else: sys.exit(f"Input file '{args.input}' does not exist")

if args.output == None:
    args.output = os.path.dirname(args.input)
    if args.output == '':
        args.output = os.getcwd()

flavors = {'l', 'c', 'b', 'cc', 'bb'}
output_templates = {}
fit_variables = [ 'FatJetGood_logsumcorrSVmass' ]
sf_label = 'sf_ptetatau21_reweighting'

for histname in fit_variables:
    h = accumulator['variables'][histname]
    h_5f = { f : sum({k : val for k, val in h.items() if k.endswith(f"_{f}")}.values()) for f in flavors }
    #h_5f.update({'DATA' : h['DATA']})
    h_3f = {'l' : h_5f['l'], 'b+bb' : h_5f['b'] + h_5f['bb'], 'c+cc' : h_5f['c'] + h_5f['cc']}
    #h_3f = {'l' : h_5f['l'], 'b+bb' : h_5f['b'] + h_5f['bb'], 'c+cc' : h_5f['c'] + h_5f['cc'], 'DATA' : h_5f['DATA']}
    for scheme, h in zip(['3f', '5f'], [h_3f, h_5f]):
        print(f"Histogram: {histname}\tScheme: {scheme}")
        output_templates[scheme] = {}
        samples = h.keys()
        samples_data = list(filter(lambda d : 'DATA' in d, samples))
        flavors   = list(filter(lambda d : 'DATA' not in d, samples))

        h_mc = h[flavors[0]]

        dense_axis = dense_axes(h_mc)[0]
        years      = get_axis_items(h_mc, 'year')
        categories = get_axis_items(h_mc, 'cat')
        variations = get_axis_items(h_mc, 'variation')
        variations_reweighting = [var for var in variations if sf_label in var]

        for year in years:
            for cat in categories:
                slicing = {'cat' : cat, 'year' : year}
                slicing_nominal = {'cat' : cat, 'year' : year, 'variation' : 'nominal'}
                # Store data shapes
                for f in samples_data:
                    sumw, sumw2 = h[f][slicing].values(), h[f][slicing].variances()
                    output_templates[scheme][f"{histname}_{year}_{cat}_{f}"] = [sumw, sumw2]
                # Extract the reweighting uncertainty
                slicing_reweighting_up = {'cat' : cat, 'year' : year, 'variation' : f'{sf_label}Up'}
                slicing_reweighting_down = {'cat' : cat, 'year' : year, 'variation' : f'{sf_label}Down'}
                for f in flavors:
                    # Store the nominal shape and its variance
                    sumw_nominal, sumw2_nominal = h[f][slicing_nominal].values(), h[f][slicing_nominal].variances()
                    sumw_up = h[f][slicing_reweighting_up].values()
                    sumw_down = h[f][slicing_reweighting_down].values()
                    variance_up = (sumw_up - sumw_nominal)**2
                    variance_down = (sumw_down - sumw_nominal)**2
                    # Compute the variance associated with the 3D reweighting as the maximum bin by bin between the up and down variances
                    variance_reweighting = np.max((variance_up, variance_down), axis=0)
                    # Save the nominal shape with the corrected variance now accounting for the 3D reweighting uncertainty
                    output_templates[scheme][f"{histname}_{year}_{cat}_QCD_{f}_nominal"] = [sumw_nominal, sumw2_nominal]
                    #output_templates[scheme][f"{histname}_{year}_{cat}_QCD_{f}_nominal"] = [sumw_nominal, sumw2_nominal + variance_reweighting]
                # Save the remaining varied shapes in the output dictionary
                for var in variations:
                    if (var == "nominal") | (var in variations_reweighting):
                        continue
                    for f in flavors:
                        slicing.update({'variation' : var})
                        sumw, sumw2 = h[f][slicing].values(), h[f][slicing].variances()
                        output_templates[scheme][f"{histname}_{year}_{cat}_QCD_{f}_{var}"] = [sumw, sumw2]

#### Saving into pickle
if not os.path.exists(args.output):
    os.makedirs(args.output)

for scheme, templates in output_templates.items():
    filename = args.input.replace('.coffea', f'_templates_{scheme}.pkl')

    if args.label:
        filename = filename.replace('.pkl', f'_{args.label}.pkl')
    filepath = os.path.join(args.output, filename)
    print(f"Saving templates file to {filepath}")
    outfile = open( filepath, 'wb' )
    pickle.dump( templates, outfile, protocol=2 )
    outfile.close()
