import os
import argparse
from collections import defaultdict

import numpy as np
import hist
from coffea.util import save, load
from coffea.processor import accumulate

import correctionlib, rich
import correctionlib.convert


pt_eta_2d_maps = [
    'FatJetGoodNMuon1_pt_eta',
    'FatJetGoodNMuon2_pt_eta',
    'FatJetGoodNMuonSJ1_pt_eta',
    'FatJetGoodNMuonSJUnique1_pt_eta',
]
pt_eta_tau21_3d_maps = [
    'FatJetGoodNMuon1_pt_eta_tau21', 'FatJetGoodNMuon1_pt_eta_tau21_bintau05',
    'FatJetGoodNMuon2_pt_eta_tau21', 'FatJetGoodNMuon2_pt_eta_tau21_bintau05',
    'FatJetGoodNMuonSJ1_pt_eta_tau21', 'FatJetGoodNMuonSJ1_pt_eta_tau21_bintau05',
    'FatJetGoodNMuonSJUnique1_pt_eta_tau21', 'FatJetGoodNMuonSJUnique1_pt_eta_tau21_bintau05',
]

def dense_axes(h):
    '''Returns the list of dense axes of a histogram.'''
    dense_axes = []
    if type(h) == dict:
        h = h[list(h.keys())[0]]
    for ax in h.axes:
        if not type(ax) in [hist.axis.StrCategory, hist.axis.IntCategory]:
            dense_axes.append(ax)
    return dense_axes

def stack_sum(stack):
    '''Returns the sum histogram of a stack (`hist.stack.Stack`) of histograms.'''
    if len(stack) == 1:
        return stack[0]
    else:
        htot = stack[0]
        for h in stack[1:]:
            htot = htot + h
        return htot

def get_axis_items(h, axis_name):
    axis = h.axes[axis_name]
    return list(axis.value(range(axis.size)))

def get_data_mc_ratio(h_data, h_qcd, h_diff):
    if type(h_data) == hist.Stack:
        h_data = stack_sum(h_data)
    if type(h_qcd) == hist.Stack:
        h_qcd = stack_sum(h_qcd)
    if type(h_diff) == hist.Stack:
        h_diff = stack_sum(h_diff)
    data = h_data.values()
    qcd = h_qcd.values()
    diff = h_diff.values()
    num = data - diff
    den = qcd
    ratio = num / den
    sumw2_num = h_data.variances() + h_diff.variances()
    sumw2_den = h_qcd.variances()
    # Statistical uncertainty on the reweighting SF taking into account
    # the uncertainty on data, QCD, top and WJets
    unc = np.sqrt( sumw2_num/den**2 + (num**2/den**4)*sumw2_den )
    # Statistical uncertainty on the reweighting SF taking into account
    # the uncertainty on data and QCD only
    unc_no_diff = np.sqrt( data/den**2 + (num**2/den**4)*sumw2_den )
    unc[np.isnan(unc)] = np.inf
    unc_no_diff[np.isnan(unc_no_diff)] = np.inf

    return ratio, unc, unc_no_diff

def overwrite_check(outfile):
    path = outfile
    version = 1
    while os.path.exists(path):
        tag = str(version).rjust(2, '0')
        path = outfile.replace('.json', f'_v{tag}.json')
        version += 1
    if path != outfile:
        print(f"The output will be saved to {path}")
    return path

def pt_reweighting(accumulator, year, histname, args):
    h = accumulator['variables'][histname]
    samples = list(filter(lambda d: 'QCD_MuEnriched' not in d, h.keys())) # Exclude the QCD_MuEnriched datasets
    samples_data = list(filter(lambda d: 'DATA' in d, samples))
    samples_mc = list(filter(lambda d: 'DATA' not in d, samples))
    samples_qcd = list(filter(lambda d: 'QCD_HT' in d, samples_mc))
    samples_vjets_top = list(filter(lambda d: (('VJets' in d) | ('SingleTop_ttbar' in d)), samples_mc))

    h_qcd = h[samples_qcd[0]]

    axes = dense_axes(h_qcd)
    categories = get_axis_items(h_qcd, 'cat')

    ratio_dict = defaultdict(float)

    for cat in categories:
        ratio_dict[cat] = {}
        for var_shape in ["nominal", "JES_TotalUp", "JES_TotalDown", "JERUp", "JERDown"]:
            slicing_mc = {'year': year, 'cat': cat, 'variation': var_shape}
            dict_qcd = {d: h[d][slicing_mc] for d in samples_qcd}
            dict_vjets_top = {d: h[d][slicing_mc] for d in samples_vjets_top}
            stack_qcd = hist.Stack.from_dict(dict_qcd)
            stack_vjets_top = hist.Stack.from_dict(dict_vjets_top)

            if 'era' in h[samples_data[0]].axes.name:
                slicing_data = {'year': year, 'cat': cat, 'era': sum}
            else:
                slicing_data = {'year': year, 'cat': cat}
            dict_data = {d: h[d][slicing_data] for d in samples_data}
            stack_data = hist.Stack.from_dict(dict_data)
            if len(stack_data) > 1:
                raise NotImplementedError
            ratio, unc, unc_no_diff = get_data_mc_ratio(stack_data, stack_qcd, stack_vjets_top)
            mod_ratio  = np.nan_to_num(ratio)
            mod_unc = np.nan_to_num(unc)
            mod_unc_no_diff = np.nan_to_num(unc_no_diff)

            ratio_dict[cat][var_shape] = {}
            ratio_dict[cat][var_shape].update({ "nominal" : mod_ratio })
            ratio_dict[cat][var_shape].update({ "statUp" : mod_ratio + mod_unc })
            ratio_dict[cat][var_shape].update({ "statDown" : mod_ratio - mod_unc })

    categories = list(ratio_dict.keys())
    shape_variations = list(ratio_dict[categories[0]].keys())
    variations = list(ratio_dict[categories[0]][shape_variations[0]].keys())
    axis_category = hist.axis.StrCategory(categories, name="cat")
    axis_shape_variation = hist.axis.StrCategory(shape_variations, name="shape_variation")
    axis_variation = hist.axis.StrCategory(variations, name="variation")
    # Stack nominal, statUp, statDown maps for each category
    stack_map = np.stack([[list(ratio_dict[cat][var_shape].values()) for var_shape in shape_variations] for cat in categories])
    sfhist = hist.Hist(axis_category, axis_shape_variation, axis_variation, *axes, data=stack_map)
    sfhist.label = "out"
    sfhist.name = f"{histname}_corr_{year}"
    description = "Reweighting SF matching the leading fatjet pT and eta MC distribution to data."
    clibcorr = correctionlib.convert.from_histogram(sfhist, flow="clamp")
    clibcorr.description = description
    cset = correctionlib.schemav2.CorrectionSet(
        schema_version=2,
        description="MC to data reweighting SF",
        corrections=[clibcorr],
    )
    rich.print(cset)

    outfile_reweighting = os.path.join(args.output, f'{histname}_{year}_reweighting.json')
    outfile_reweighting = overwrite_check(outfile_reweighting)
    print(f"Saving pt reweighting factors in {outfile_reweighting}")
    with open(outfile_reweighting, "w") as fout:
        fout.write(cset.json(exclude_unset=True))
    fout.close()
    print(f"Loading correction from {outfile_reweighting}...")
    cset = correctionlib.CorrectionSet.from_file(outfile_reweighting)
    print("(cat): inclusive,", "(var): nominal")
    pt_corr = cset[sfhist.name]
    pos  = np.array([0, 1, 0, 1, 0], dtype=int)
    pt  = np.array([50, 100, 400, 500, 1000], dtype=float)
    eta = np.array([-2, -1, 0, 1, 2], dtype=float)
    print("pos =", pos)
    print("pt =", pt)
    print("eta =", eta)
    if histname in pt_eta_2d_maps:
        args = (pos, pt, eta)

    elif histname in pt_eta_tau21_3d_maps:
        tau21 = np.array([0.1, 0.35, 0.65, 0.8, 0.9], dtype=float)
        print("tau21 =", tau21)
        args = (pos, pt, eta, tau21)
    for var_shape in shape_variations:
        categorical_args = ['inclusive', var_shape, 'nominal']
        print(categorical_args)
        print(pt_corr.evaluate(*categorical_args, *args))
        print()

parser = argparse.ArgumentParser(description='Accumulate coffea outputs')
# Inputs
parser.add_argument('-i','--inputfiles', required=True, type=str, nargs="+",
                    help='List of coffea input files')
parser.add_argument("-o", "--output", required=True, type=str,
                    help="Output folder")
parser.add_argument('-y', '--year', required=True, type=str,
                    help='Data-taking year of dataset to update')
args = parser.parse_args()

files = list(set([f for f in args.inputfiles]))

accumulator = accumulate([ load(f) for f in files])

for histname in pt_eta_2d_maps + pt_eta_tau21_3d_maps:
    pt_reweighting(accumulator=accumulator, year=args.year, histname=histname, args=args)

save(accumulator, os.path.join(args.output, "output_accumulated_with_QCD_HT.coffea"))
