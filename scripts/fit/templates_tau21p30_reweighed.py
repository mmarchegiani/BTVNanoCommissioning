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

def get_model_name(cat):
    if "pass" in cat:
        region = "pass"
    elif "fail" in cat:
        region = "fail"
    return ''.join(cat.split(region))

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
fit_variables = [ 'FatJetGood_logsumcorrSVmass_tau21' ]
sf_label = 'sf_ptetatau21_reweighting'
variations_psweight = ['psWeight_isrUp', 'psWeight_isrDown', 'psWeight_fsrUp', 'psWeight_fsrDown']
# Only tau21 < 0.3 cut is considered in this script
cuts_tau21_2d = [0.3]
schemes = ['3f', '5f']
output_templates = {k : {'3f' : {}, '5f' : {}} for k in ['inclusive'] + cuts_tau21_2d}
sumw_passfail_data = {k : {'3f' : {}, '5f' : {}} for k in ['inclusive'] + cuts_tau21_2d}
sumw_passfail_nominal = {k : {'3f' : {}, '5f' : {}} for k in ['inclusive'] + cuts_tau21_2d}
map_reweighting = {k : {'3f' : {}, '5f' : {}} for k in ['inclusive'] + cuts_tau21_2d}

# Dictionary to save the numerator and denominator of the renormalization factor for the ISR, FSR, JES, JER and pileup shape variations
shapes_to_renormalize = ["psWeight_isrUp", "psWeight_isrDown", "psWeight_fsrUp", "psWeight_fsrDown", "JES_TotalUp", "JES_TotalDown", "JERUp", "JERDown", "pileupUp", "pileupDown"]
sumw_passfail_shape = {k : {'3f' : {}, '5f' : {}} for k in ['inclusive'] + cuts_tau21_2d}
map_renormalization = {k : {'3f' : {}, '5f' : {}} for k in ['inclusive'] + cuts_tau21_2d}

"""
Structure of the dictionary:
den_renormalization = {
    "inclusive" : {
        "3f" : {
            "msd40btagDDCvLV2Hwp_Pt-450to500" : {
                "JES_TotalUp" : [ THIS IS THE LIST WHERE WE STORE THE NUMERATOR OF THE REWEIGHTING FACTOR ]
            }
        }
    },
    0.2 : {...}
}

Same is for the denominator dict, and extract the reweighting factor from their ratio.
"""

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
    for scheme, h_qcd_muenriched, h_vjets_top, h_madgraph in zip( schemes,
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
        
        # Define the list of the tau21 cuts: 'inclusive' only for 1D histograms, while a list of cuts for the 2D histograms
        if is1D:
            cuts_tau21 = ['inclusive']
        elif is2D:
            cuts_tau21 = deepcopy(cuts_tau21_2d)

        # Initialize the dictionaries to compute the sumw in the pass+fail region, for each scheme and tau21_cut
        # N.B.: The list `models` is emptied at each loop over the schemes, and is filled only at the first 'pass' or 'fail' category
        # encountered in the loop for the corresponding model.
        models = []
        for year in years:
            for cat in categories:
                if not any(['pass' in cat, 'fail' in cat]):
                    continue
                else:
                    suffix = get_model_name(cat)
                if not suffix in models:
                    models.append(suffix)
                    for tau21_cut in cuts_tau21:
                        sumw_passfail_data[tau21_cut][scheme][f"{histname}_{year}_{suffix}"] = 0.0
                        sumw_passfail_nominal[tau21_cut][scheme][f"{histname}_{year}_{suffix}"] = 0.0
                        sumw_passfail_shape[tau21_cut][scheme][f"{histname}_{year}_{suffix}"] = {}
                        map_renormalization[tau21_cut][scheme][f"{histname}_{year}_{suffix}"] = {}
                        
                        for shape_name in shapes_to_renormalize:
                            sumw_passfail_shape[tau21_cut][scheme][f"{histname}_{year}_{suffix}"][shape_name] = 0.0

        for year in years:
            # Initialize the reweighting dictionaries to zero
            for cat in categories:
                if not any(['pass' in cat, 'fail' in cat]):
                    continue
                else:
                    model_name = get_model_name(cat)
                slicing_baseline = {'cat' : cat, 'year' : year}
                
                for tau21_cut in cuts_tau21:
                    print(f"Scheme: {scheme}\t Category: {cat}\t tau21 < {tau21_cut}")
                    if is1D:
                        kwargs = {}
                    elif is2D:
                        index_cut = np.where(axis_tau21.edges == tau21_cut)[0][0]
                        kwargs = {'is2D' : True, 'index_cut' : index_cut}
                    slicing_nominal = deepcopy(slicing_baseline)
                    slicing_nominal['variation'] = 'nominal'
                    # Store data shapes
                    for f in samples_data:
                        sumw_data, sumw2_data = get_sumw(h_data, f, slicing_baseline, **kwargs)
                        output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_{f}"] = [sumw_data, sumw2_data]
                        print(cat, model_name, tau21_cut)
                        sumw_passfail_data[tau21_cut][scheme][f"{histname}_{year}_{model_name}"] += sumw_data
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
                        sumw_passfail_nominal[tau21_cut][scheme][f"{histname}_{year}_{model_name}"] += sumw

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
                            if var in shapes_to_renormalize:
                                sumw_passfail_shape[tau21_cut][scheme][f"{histname}_{year}_{model_name}"][var] += sumw
                    # Save the ISR/FSR varied shapes taken from the Madgraph samples and the QCD flavor composition variation
                    for var in variations_psweight:
                        for f in flavors:
                            # Store the nominal shape and its variance for both the PYTHIA and Madgraph samples
                            sumw_nominal_qcd_muenriched, sumw2_nominal_qcd_muenriched = get_sumw(h_qcd_muenriched, f, slicing_nominal, **kwargs)
                            sumw_nominal_madgraph, sumw2_nominal_madgraph = get_sumw(h_madgraph, f, slicing_nominal, **kwargs)
                            slicing_var.update({'variation' : var})
                            sumw_madgraph, sumw2_madgraph = get_sumw(h_madgraph, f, slicing_var, **kwargs)
                            sumw_vjets_top, sumw2_vjets_top = get_sumw(h_vjets_top, f, slicing_var, **kwargs)
                            ratio_madgraph = np.nan_to_num(sumw_madgraph / sumw_nominal_madgraph, nan=1.0)
                            # Save the ISR/FSR varied shapes
                            sumw = ratio_madgraph * sumw_nominal_qcd_muenriched + sumw_vjets_top
                            sumw2 = ratio_madgraph**2 * sumw2_nominal_qcd_muenriched + sumw2_vjets_top
                            output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_MC_{f}_{var}"] = [sumw, sumw2]
                            if var in shapes_to_renormalize:
                                sumw_passfail_shape[tau21_cut][scheme][f"{histname}_{year}_{model_name}"][var] += sumw
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
                        ratio_madgraph_pythia = np.nan_to_num(sumw_nominal_madgraph / sumw_nominal_qcd_muenriched, nan=1.0)
                        output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_MC_{f}_QCDFlvComposUp"] = [sumw, sumw2]
                        sumw = sumw_nominal_qcd_muenriched / ratio_madgraph_pythia
                        output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_MC_{f}_QCDFlvComposDown"] = [sumw, sumw2]
for tau21_cut in cuts_tau21_2d:
    for scheme in schemes:
        # Reweighting on fit variable in pass+fail
        for template_name, num in sumw_passfail_data[tau21_cut][scheme].items():
            den = sumw_passfail_nominal[tau21_cut][scheme][template_name]
            map_reweighting[tau21_cut][scheme][template_name] = np.nan_to_num(num / den, nan=1.0)

        # Renormalization of shape variations such that the normalization matches the nominal normalization
        for template_name, var_dict in sumw_passfail_shape[tau21_cut][scheme].items():
            num = sumw_passfail_nominal[tau21_cut][scheme][template_name]
            for var, den in var_dict.items():
                map_renormalization[tau21_cut][scheme][template_name][var] = sum(np.nan_to_num(num, nan=0)) / sum(np.nan_to_num(den, nan=0))

if not os.path.exists(args.output):
    os.makedirs(args.output)

label_tau21 = {tau21 : f'tau21p{int(100*tau21)}' for tau21 in cuts_tau21_2d}
label_tau21.update({'inclusive' : 'inclusive'})

"""
#### Saving bare templates into pickle file
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
"""

"""
# Save templates with MC-to-data reweighting in the pass+fail region
for histname in ["FatJetGood_logsumcorrSVmass"]:
    for scheme in ['3f', '5f']:
        for year in years:
            for tau21_cut, templates_dict in output_templates.items():
                filename = args.input.replace('.coffea', f'_templates_{scheme}_{label_tau21[tau21_cut]}_reweighed.pkl')
                if tau21_cut == "inclusive":
                    templates_reweighed = {}
                    for cat in categories:
                        model_name = get_model_name(cat)
                        reweighting_factor = map_reweighting[tau21_cut][scheme][f"{histname}_{year}_{model_name}"]
                        for template_name, (sumw, sumw2) in output_templates[tau21_cut][scheme].items():
                            if 'DATA' in template_name:
                                templates_reweighed[template_name] = [ sumw, sumw2 ]
                            else:
                                templates_reweighed[template_name] = [ reweighting_factor * sumw, reweighting_factor**2 * sumw2 ]
                filepath = os.path.join(args.output, filename)
                print(f"Saving templates file to {filepath}")
                outfile = open( filepath, 'wb' )
                pickle.dump( templates_reweighed, outfile, protocol=2 )
                outfile.close()
"""

# Save templates with MC-to-data reweighting in the pass+fail region and
# with renormalized varied shapes to extract the overall normalization effect from the shape variation
# for tau21 < 0.3
for histname in ["FatJetGood_logsumcorrSVmass_tau21"]:
    for scheme in ['3f', '5f']:
        if scheme == '3f':
            flavors = {'l', 'c+cc', 'b+bb'}
        elif scheme == '5f':
            flavors = {'l', 'c', 'b', 'cc', 'bb'}
        for year in years:
            for tau21_cut, templates_dict in output_templates.items():
                if not tau21_cut == 0.3:
                    continue
                templates_fit_variable_reweighed = {}
                for cat in categories:
                    if not any(['pass' in cat, 'fail' in cat]):
                        continue
                    model_name = get_model_name(cat)
                    reweighting_factor = map_reweighting[tau21_cut][scheme][f"{histname}_{year}_{model_name}"]
                    templates_fit_variable_reweighed[f"{histname}_{year}_{cat}_DATA"] = output_templates[tau21_cut][scheme][f"{histname}_{year}_{cat}_DATA"]
                    for f in flavors:
                        for var in variations + variations_psweight + ["QCDFlvComposUp", "QCDFlvComposDown"]:
                            if var in variations_reweighting:
                                continue
                            template_name = f"{histname}_{year}_{cat}_MC_{f}_{var}"
                            sumw, sumw2 = output_templates[tau21_cut][scheme][template_name]
                            if var in shapes_to_renormalize:
                                renormalizing_factor = map_renormalization[tau21_cut][scheme][f"{histname}_{year}_{model_name}"][var]
                            else:
                                renormalizing_factor = 1.0
                            templates_fit_variable_reweighed[template_name] = [ reweighting_factor * renormalizing_factor * sumw, reweighting_factor**2 * renormalizing_factor**2 * sumw2 ]
                filename = args.input.replace('.coffea', f'_templates_{scheme}_{label_tau21[tau21_cut]}_fit_variable_reweighed.pkl')
                filepath = os.path.join(args.output, filename)
                print(f"Saving templates file to {filepath}")
                outfile = open( filepath, 'wb' )
                pickle.dump( templates_fit_variable_reweighed, outfile, protocol=2 )
                outfile.close()
