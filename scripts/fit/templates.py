import sys, os
from copy import deepcopy
import argparse
import numpy as np
from coffea.util import load
import hist
import pickle

def get_dense_axes(h):
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

def get_sumw(h, flavor, slicing, sumw2=True, is2D=False, index_cut=None):
    h_s = deepcopy(h[flavor][slicing])
    if is2D & (index_cut != None):
        assert type(index_cut) in [int, np.int64], "The index to slice the histogram must be an integer."
        h_s = h_s[:,:index_cut:sum]
    if sumw2:
        return h_s.values(), h_s.variances()
    else:
        return h_s.values()

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
samples_qcd_muenriched = ('QCD_MuEnriched')
samples_vjets_top = ('VJets', 'SingleTop_ttbar')
samples_madgraph = ('QCD_HT')
fit_variables = [ 'FatJetGood_logsumcorrSVmass', 'FatJetGood_logsumcorrSVmass_tau21' ]
sf_label = 'sf_ptetatau21_reweighting'
variations_psweight = ['psWeight_isrUp', 'psWeight_isrDown', 'psWeight_fsrUp', 'psWeight_fsrDown']
cuts_tau21_2d = [0.2, 0.25, 0.3, 0.35, 0.4]
output_templates = {k : {'3f' : {}, '5f' : {}} for k in ['inclusive'] + cuts_tau21_2d}

for histname in fit_variables:
    h = accumulator['variables'][histname]

    # Define the templates split by dataset: this is needed because we cannot sum QCD_MuEnriched and the other samples since they do not have the same number of variations
    h_data = { 'DATA' : h['DATA'] }
    h_5f_qcd_muenriched = { f : sum({k : val for k, val in h.items() if k.startswith(samples_qcd_muenriched) and k.endswith(f"_{f}")}.values()) for f in flavors }
    h_5f_vjets_top = { f : sum({k : val for k, val in h.items() if k.startswith(samples_vjets_top) and k.endswith(f"_{f}")}.values()) for f in flavors }
    h_5f_madgraph = { f : sum({k : val for k, val in h.items() if k.startswith(samples_madgraph) and k.endswith(f"_{f}")}.values()) for f in flavors }
    h_3f_qcd_muenriched = {'l' : h_5f_qcd_muenriched['l'], 'b+bb' : h_5f_qcd_muenriched['b'] + h_5f_qcd_muenriched['bb'], 'c+cc' : h_5f_qcd_muenriched['c'] + h_5f_qcd_muenriched['cc']}
    h_3f_vjets_top = {'l' : h_5f_vjets_top['l'], 'b+bb' : h_5f_vjets_top['b'] + h_5f_vjets_top['bb'], 'c+cc' : h_5f_vjets_top['c'] + h_5f_vjets_top['cc']}
    h_3f_madgraph = {'l' : h_5f_madgraph['l'], 'b+bb' : h_5f_madgraph['b'] + h_5f_madgraph['bb'], 'c+cc' : h_5f_madgraph['c'] + h_5f_madgraph['cc']}
    for scheme, h_qcd_muenriched, h_vjets_top, h_madgraph in zip( ['3f', '5f'],
                                                                  [h_3f_qcd_muenriched, h_5f_qcd_muenriched],
                                                                  [h_3f_vjets_top, h_5f_vjets_top],
                                                                  [h_3f_madgraph, h_5f_madgraph] ):
        print(f"Histogram: {histname}\tScheme: {scheme}")
        samples = h_qcd_muenriched.keys()
        samples_data = h_data.keys()
        flavors   = list(filter(lambda d : 'DATA' not in d, samples))

        h_l = h_qcd_muenriched[flavors[0]]

        dense_axes = get_dense_axes(h_l)
        ndim = len(dense_axes)
        is1D = ndim == 1
        is2D = ndim == 2
        if is2D:
            if not dense_axes[-1].name == "FatJetGood.tau21":
                raise Exception("The last axis of the 2D histogram must be 'FatJetGood.tau21'.")
            axis_tau21 = dense_axes[-1]
            tau21 = axis_tau21.centers
            hi = int(np.where(tau21 < 0.4)[0][-1])
        years      = get_axis_items(h_l, 'year')
        categories = get_axis_items(h_l, 'cat')
        variations = get_axis_items(h_l, 'variation')
        variations_reweighting = [var for var in variations if sf_label in var]

        for year in years:
            for cat in categories:
                slicing_baseline = {'cat' : cat, 'year' : year}
                if is1D:
                    cuts_tau21 = ['inclusive']
                elif is2D:
                    cuts_tau21 = deepcopy(cuts_tau21_2d)
                for tau21_cut in cuts_tau21:
                    print(f"Saving templates for tau21 cut at {tau21_cut}")
                    if is1D:
                        kwargs = {}
                    elif is2D:
                        index_cut = np.where(axis_tau21.edges == tau21_cut)[0][0]
                        kwargs = {'is2D' : True, 'index_cut' : index_cut}
                    slicing_nominal = deepcopy(slicing_baseline)
                    slicing_nominal['variation'] = 'nominal'
                    # Store data shapes
                    for f in samples_data:
                        sumw, sumw2 = get_sumw(h_data, f, slicing_baseline, **kwargs)
                        output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_{f}"] = [sumw, sumw2]
                    # Extract the reweighting uncertainty
                    slicing_reweighting_up = deepcopy(slicing_baseline)
                    slicing_reweighting_down = deepcopy(slicing_baseline)
                    slicing_reweighting_up.update({'variation' : f'{sf_label}Up'})
                    slicing_reweighting_down.update({'variation' : f'{sf_label}Down'})
                    for f in flavors:
                        # Store the nominal shape and its variance
                        sumw_nominal_qcd_muenriched, sumw2_nominal_qcd_muenriched = get_sumw(h_qcd_muenriched, f, slicing_nominal, **kwargs)
                        sumw_nominal_vjets_top, sumw2_nominal_vjets_top = get_sumw(h_vjets_top, f, slicing_nominal, **kwargs)
                        sumw_up_qcd_muenriched = get_sumw(h_qcd_muenriched, f, slicing_reweighting_up, sumw2=False, **kwargs)
                        sumw_down_qcd_muenriched = get_sumw(h_qcd_muenriched, f, slicing_reweighting_down, sumw2=False, **kwargs)
                        variance_up = (sumw_up_qcd_muenriched - sumw_nominal_qcd_muenriched)**2
                        variance_down = (sumw_down_qcd_muenriched - sumw_nominal_qcd_muenriched)**2
                        # Compute the variance associated with the 3D reweighting as the maximum bin by bin between the up and down variances
                        variance_reweighting = np.max((variance_up, variance_down), axis=0)
                        # Save the nominal shape with the corrected variance now accounting for the 3D reweighting uncertainty
                        sumw = sumw_nominal_qcd_muenriched + sumw_nominal_vjets_top
                        sumw2 = sumw2_nominal_qcd_muenriched + sumw2_nominal_vjets_top + variance_reweighting
                        output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_MC_{f}_nominal"] = [sumw, sumw2]
                    # Save the remaining varied shapes in the output dictionary
                    slicing_var = deepcopy(slicing_baseline)
                    for var in variations:
                        if (var == "nominal") | (var in variations_reweighting):
                            continue
                        for f in flavors:
                            slicing_var.update({'variation' : var})
                            sumw_qcd_muenriched, sumw2_qcd_muenriched = get_sumw(h_qcd_muenriched, f, slicing_var, **kwargs)
                            sumw_vjets_top, sumw2_vjets_top = get_sumw(h_vjets_top, f, slicing_var, **kwargs)
                            sumw = sumw_qcd_muenriched + sumw_vjets_top
                            sumw2 = sumw2_qcd_muenriched + sumw2_vjets_top
                            output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_MC_{f}_{var}"] = [sumw, sumw2]
                    # Save the ISR/FSR varied shapes taken from the Madgraph samples and the QCD flavor composition variation
                    for var in variations_psweight:
                        for f in flavors:
                            # Store the nominal shape and its variance for both the PYTHIA and Madgraph samples
                            sumw_nominal_qcd_muenriched, sumw2_nominal_qcd_muenriched = get_sumw(h_qcd_muenriched, f, slicing_nominal, **kwargs)
                            sumw_nominal_madgraph, sumw2_nominal_madgraph = get_sumw(h_madgraph, f, slicing_nominal, **kwargs)
                            slicing_var.update({'variation' : var})
                            sumw_madgraph, sumw2_madgraph = get_sumw(h_madgraph, f, slicing_var, **kwargs)
                            sumw_vjets_top, sumw2_vjets_top = get_sumw(h_vjets_top, f, slicing_var, **kwargs)
                            ratio_madgraph = sumw_madgraph / sumw_nominal_madgraph
                            # Save the ISR/FSR varied shapes
                            sumw = ratio_madgraph * sumw_nominal_qcd_muenriched + sumw_vjets_top
                            sumw2 = ratio_madgraph**2 * sumw2_nominal_qcd_muenriched + sumw2_vjets_top
                            output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_MC_{f}_{var}"] = [sumw, sumw2]
                    for f in flavors:
                        # Store the nominal shape and its variance for both the PYTHIA and Madgraph samples
                        sumw_nominal_qcd_muenriched, sumw2_nominal_qcd_muenriched = get_sumw(h_qcd_muenriched, f, slicing_nominal, **kwargs)
                        sumw_nominal_madgraph, sumw2_nominal_madgraph = get_sumw(h_madgraph, f, slicing_nominal, **kwargs)
                        # Here we save the QCD flavor composition variation: the Up variation is defined as the sumw of the nominal Madgraph shape, while
                        # the Down variation is defined as the inverse of the sumw of the nominal Magraph shape.
                        # In this way, a 20% up variation will correspond to a 16.7% down variation since 1/1.2 = 0.833
                        sumw_nominal_vjets_top, sumw2_nominal_vjets_top = get_sumw(h_vjets_top, f, slicing_nominal, **kwargs)
                        sumw = sumw_nominal_madgraph + sumw_nominal_vjets_top
                        sumw2 = sumw2_nominal_madgraph + sumw2_nominal_vjets_top
                        ratio_madgraph_pythia = sumw_nominal_madgraph / sumw_nominal_qcd_muenriched
                        output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_MC_{f}_QCDFlvComposUp"] = [sumw, sumw2]
                        sumw = sumw_nominal_qcd_muenriched / ratio_madgraph_pythia
                        #sumw2 = 0.0
                        output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_MC_{f}_QCDFlvComposDown"] = [sumw, sumw2]

#### Saving into pickle
if not os.path.exists(args.output):
    os.makedirs(args.output)

label_tau21 = {tau21 : f'tau21p{int(100*tau21)}' for tau21 in cuts_tau21_2d}
label_tau21.update({'inclusive' : 'inclusive'})

for tau21_cut, templates_dict in output_templates.items():
    for scheme, templates in templates_dict.items():
        filename = args.input.replace('.coffea', f'_templates_{scheme}_{label_tau21[tau21_cut]}.pkl')

        if args.label:
            filename = filename.replace('.pkl', f'_{args.label}.pkl')
        filepath = os.path.join(args.output, filename)
        print(f"Saving templates file to {filepath}")
        outfile = open( filepath, 'wb' )
        pickle.dump( templates, outfile, protocol=2 )
        outfile.close()
