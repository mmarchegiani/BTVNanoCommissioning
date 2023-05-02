import os
from collections import defaultdict

import awkward as ak
import numpy as np
import hist
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.util import save, load

import correctionlib, rich
import correctionlib.convert

from workflows.fatjet_base import fatjetBaseProcessor
from pocket_coffea.utils.configurator import Configurator
from config.fatjet_base.custom.cuts import get_ptmsd, mutag_fatjet_sel, mutag_subjet_sel

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

def get_data_mc_ratio(h1, h2, h_diff):
    if type(h1) == hist.Stack:
        h1 = stack_sum(h1)
    if type(h2) == hist.Stack:
        h2 = stack_sum(h2)
    if type(h_diff) == hist.Stack:
        h_diff = stack_sum(h_diff)
    num = h1.values()
    den = h2.values()
    diff = h_diff.values()
    ratio = (num - diff) / den
    unc = np.sqrt(num - diff) / den
    unc[np.isnan(unc)] = np.inf

    return ratio, unc

def overwrite_check(outfile):
    path = outfile
    version = 1
    while os.path.exists(path):
        tag = str(version).rjust(2, '0')
        path = outfile.replace('.coffea', f'_v{tag}.coffea')
        version += 1
    if path != outfile:
        print(f"The output will be saved to {path}")
    return path

class ptReweightProcessor(fatjetBaseProcessor):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.pt_eta_2d_maps = [
            'FatJetGoodNMuon1_pt_eta',# 'FatJetGoodNMuon1_pt_eta_bineta0p40', 'FatJetGoodNMuon1_pt_eta_binpt20', 'FatJetGoodNMuon1_pt_eta_binpt20_bineta0p40',
            'FatJetGoodNMuon2_pt_eta',# 'FatJetGoodNMuon2_pt_eta_bineta0p40', 'FatJetGoodNMuon2_pt_eta_binpt20', 'FatJetGoodNMuon2_pt_eta_binpt20_bineta0p40',
            'FatJetGoodNMuonSJ1_pt_eta',# 'FatJetGoodNMuonSJ1_pt_eta_bineta0p40', 'FatJetGoodNMuonSJ1_pt_eta_binpt20', 'FatJetGoodNMuonSJ1_pt_eta_binpt20_bineta0p40',
            'FatJetGoodNMuonSJUnique1_pt_eta',# 'FatJetGoodNMuonSJUnique1_pt_eta_bineta0p40', 'FatJetGoodNMuonSJUnique1_pt_eta_binpt20', 'FatJetGoodNMuonSJUnique1_pt_eta_binpt20_bineta0p40',
        ]
        self.pt_eta_tau21_3d_maps = [
            'FatJetGoodNMuon1_pt_eta_tau21',# 'FatJetGoodNMuon1_pt_eta_tau21_bineta0p40', 'FatJetGoodNMuon1_pt_eta_tau21_binpt20', 'FatJetGoodNMuon1_pt_eta_tau21_binpt20_bineta0p40',
            'FatJetGoodNMuon2_pt_eta_tau21',# 'FatJetGoodNMuon2_pt_eta_tau21_bineta0p40', 'FatJetGoodNMuon2_pt_eta_tau21_binpt20', 'FatJetGoodNMuon2_pt_eta_tau21_binpt20_bineta0p40',
            'FatJetGoodNMuonSJ1_pt_eta_tau21',# 'FatJetGoodNMuonSJ1_pt_eta_tau21_bineta0p40', 'FatJetGoodNMuonSJ1_pt_eta_tau21_binpt20', 'FatJetGoodNMuonSJ1_pt_eta_tau21_binpt20_bineta0p40',
            'FatJetGoodNMuonSJUnique1_pt_eta_tau21',# 'FatJetGoodNMuonSJUnique1_pt_eta_tau21_bineta0p40', 'FatJetGoodNMuonSJUnique1_pt_eta_tau21_binpt20', 'FatJetGoodNMuonSJUnique1_pt_eta_tau21_binpt20_bineta0p40'
        ]
        for histname in self.pt_eta_2d_maps + self.pt_eta_tau21_3d_maps:
            if not histname in self.cfg.variables.keys():
                raise Exception(f"'{histname}' is not present in the histogram keys.")

    def apply_object_preselection(self, variation):
        super().apply_object_preselection(variation)
        pt_min = 450.
        msd = 40.
        cuts_dict = {
            "FatJetGoodNMuon1" : [get_ptmsd(pt_min, msd), mutag_fatjet_sel(nmu=1)],
            "FatJetGoodNMuon2" : [get_ptmsd(pt_min, msd), mutag_fatjet_sel(nmu=2)],
            "FatJetGoodNMuonSJ1" : [get_ptmsd(pt_min, msd), mutag_subjet_sel(unique_matching=False)],
            "FatJetGoodNMuonSJUnique1" : [get_ptmsd(pt_min, msd), mutag_subjet_sel(unique_matching=True)],
        }
        self._ak8jet_collections = list(cuts_dict.keys())
        for coll, cuts in cuts_dict.items():
            masks = [cut.get_mask(self.events) for cut in cuts]
            mask_total = ak.ones_like(masks[0], dtype=bool)
            for mask in masks:
                mask_total = mask_total & mask
            self.events[coll] = self.events.FatJetGood[mask_total]
            self.events[coll] = ak.with_field(self.events[coll], ak.local_index(self.events[coll].pt), "pos")

    def count_objects(self, variation):
        super().count_objects(variation)
        for coll in self._ak8jet_collections:
            self.events[f"n{coll}"] = ak.num(self.events[coll])

    def pt_reweighting(self, accumulator, year, histname):
        h = accumulator['variables'][histname]
        samples = h.keys()
        samples_data = list(filter(lambda d: 'DATA' in d, samples))
        samples_mc = list(filter(lambda d: 'DATA' not in d, samples))
        samples_qcd = list(filter(lambda d: 'QCD' in d, samples_mc))
        samples_vjets_top = list(filter(lambda d: (('VJets' in d) | ('SingleTop_ttbar' in d)), samples_mc))

        h_qcd = h[samples_qcd[0]]

        axes = dense_axes(h_qcd)
        categories = get_axis_items(h_qcd, 'cat')

        ratio_dict = defaultdict(float)

        for cat in categories:
            slicing_mc_nominal = {'year': year, 'cat': cat, 'variation': 'nominal'}
            dict_qcd_nominal = {d: h[d][slicing_mc_nominal] for d in samples_qcd}
            dict_vjets_top_nominal = {d: h[d][slicing_mc_nominal] for d in samples_vjets_top}
            stack_qcd_nominal = hist.Stack.from_dict(dict_qcd_nominal)
            stack_vjets_top_nominal = hist.Stack.from_dict(dict_vjets_top_nominal)

            if 'era' in h[samples_data[0]].axes.name:
                slicing_data = {'year': year, 'cat': cat, 'era': sum}
            else:
                slicing_data = {'year': year, 'cat': cat}
            dict_data = {d: h[d][slicing_data] for d in samples_data}
            stack_data = hist.Stack.from_dict(dict_data)
            if len(stack_data) > 1:
                raise NotImplementedError
            ratio, unc = get_data_mc_ratio(stack_data, stack_qcd_nominal, stack_vjets_top_nominal)
            mod_ratio  = np.nan_to_num(ratio)
            if histname in self.pt_eta_2d_maps:
                mod_ratio[mod_ratio == 0.0] = 1

            ratio_dict.update({ cat : mod_ratio })

        axis_category = hist.axis.StrCategory(list(ratio_dict.keys()), name="cat")
        sfhist = hist.Hist(axis_category, *axes, data=np.stack(list(ratio_dict.values())))
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

        outfile_reweighting = os.path.join(self.cfg.output, f'{histname}_{year}_reweighting.json')
        outfile_reweighting = overwrite_check(outfile_reweighting)
        print(f"Saving pt reweighting factors in {outfile_reweighting}")
        with open(outfile_reweighting, "w") as fout:
            fout.write(cset.json(exclude_unset=True))
        fout.close()
        print(f"Loading correction from {outfile_reweighting}...")
        cset = correctionlib.CorrectionSet.from_file(outfile_reweighting)
        print("inclusive:")
        pt_corr = cset[sfhist.name]
        pos  = np.array([0, 1, 0, 1, 0], dtype=int)
        pt  = np.array([50, 100, 400, 500, 1000], dtype=float)
        eta = np.array([-2, -1, 0, 1, 2], dtype=float)
        print("pos =", pos)
        print("pt =", pt)
        print("eta =", eta)
        if histname in self.pt_eta_2d_maps:
            args = (pos, pt, eta)

        elif histname in self.pt_eta_tau21_3d_maps:
            tau21 = np.array([0.1, 0.35, 0.65, 0.8, 0.9], dtype=float)
            print("tau21 =", tau21)
            args = (pos, pt, eta, tau21)
        print(pt_corr.evaluate('inclusive', *args))

    def postprocess(self, accumulator):
        '''
        Rescale MC histograms by the total sum of the genweights, read from the
        output computed before skimming.
        '''

        years = self.cfg.dataset["filter"]["year"]
        if len(years) > 1:
            raise Exception("Only one data-taking year can be processed at a time.")
        else:
            year = years[0]

        for histname in self.pt_eta_2d_maps + self.pt_eta_tau21_3d_maps:
            self.pt_reweighting(accumulator=accumulator, year=year, histname=histname)

        return accumulator
