import os
import awkward as ak
import correctionlib

from workflows.fatjet_base import fatjetBaseProcessor
from pocket_coffea.utils.configurator import Configurator
from lib.sv import *
from config.fatjet_base.custom.cuts import get_ptmsd, mutag_fatjet_sel

class templatesProcessor(fatjetBaseProcessor):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        if not "histograms_to_reweigh" in self.cfg.workflow_options.keys():
            raise Exception("The entry of the config file 'workflow_options' does not contain a key 'histograms_to_reweigh'. Please specify it in the config file.")
        if not "reweighting_scheme" in self.cfg.workflow_options.keys():
            raise Exception("The entry of the config file 'workflow_options' does not contain a key 'reweighting_scheme'. Please specify it in the config file.")
        if not "reweighting_map" in self.cfg.workflow_options.keys():
            raise Exception("The entry of the config file 'workflow_options' does not contain a key 'reweighting_map'. Please specify it in the config file.")
        self.histograms_to_reweigh = self.cfg.workflow_options["histograms_to_reweigh"]
        self.reweighting_scheme = self.cfg.workflow_options["reweighting_scheme"]
        self.reweighting_map = self.cfg.workflow_options["reweighting_map"]
        possible_schemes = ["pteta", "ptetatau21"]
        if not self.reweighting_scheme in possible_schemes:
            raise Exception(f"The reweighting scheme '{self.reweighting_scheme}' does not exist. Please specify a reweighting scheme among {possible_schemes}.")
        if self.reweighting_map == None:
            raise Exception(f"The reweighting map file is not specified. Please specify the key `reweighting_map` in the config file.")
        else:
            if not os.path.exists(self.reweighting_map):
                raise Exception(f"The reweighting map file '{self.reweighting_map}' does not exist. Please check the path of the reweighting map file.")

    def apply_object_preselection(self, variation):
        super().apply_object_preselection(variation)

        pt_min = 350.
        msd = 40.
        cuts = [get_ptmsd(pt_min, msd), mutag_fatjet_sel(nmu=1)]
        masks = [cut.get_mask(self.events) for cut in cuts]
        mask_total = ak.ones_like(masks[0], dtype=bool)
        for mask in masks:
            mask_total = mask_total & mask

        # Apply (pt, msd) cuts and muon-tagging of the AK8 jet requiring at least one muon inside the AK8 jet
        self.events["FatJetGood"] = self.events.FatJetGood[mask_total]
        
        # Restrict analysis to leading and subleading jets only
        self.events["FatJetGood"] = self.events.FatJetGood[ak.local_index(self.events.FatJetGood, axis=1) < 2]

    def pteta_reweighting(self):
        '''Reweighting scale factor based on the leading fatjet pT, eta'''
        cat = 'inclusive'
        cset = correctionlib.CorrectionSet.from_file(self.reweighting_map)
        keys = list(cset.keys())
        assert len(keys) == 1, f"The correction has {len(keys)} keys. The choice of the key is ambiguous."
        key = list(cset.keys())[0]
        pteta_corr = cset[key]

        nfatjet  = ak.num(self.events.FatJetGood.pt)
        pos = ak.flatten(ak.local_index(self.events.FatJetGood.pt))
        pt  = ak.flatten(self.events.FatJetGood.pt)
        eta = ak.flatten(self.events.FatJetGood.eta)

        weight = pteta_corr.evaluate(cat, pos, pt, eta)
        weight = ak.unflatten(weight, nfatjet)

        self.weight_2d = weight

    def ptetatau21_reweighting(self):
        '''Reweighting scale factor based on the leading fatjet pT, eta, tau21'''
        cat = 'inclusive'
        cset = correctionlib.CorrectionSet.from_file(self.reweighting_map)
        keys = list(cset.keys())
        assert len(keys) == 1, f"The correction has {len(keys)} keys. The choice of the key is ambiguous."
        key = list(cset.keys())[0]
        ptetatau21_corr = cset[key]

        nfatjet  = ak.num(self.events.FatJetGood.pt)
        pos = ak.flatten(ak.local_index(self.events.FatJetGood.pt))
        pt = ak.flatten(self.events.FatJetGood.pt)
        eta = ak.flatten(self.events.FatJetGood.eta)
        tau21 = ak.flatten(self.events.FatJetGood.tau21)

        weight = ptetatau21_corr.evaluate(cat, pos, pt, eta, tau21)
        weight = ak.unflatten(weight, nfatjet)

        self.weight_3d = weight

    def process_extra_after_presel(self, variation):
        if self._sample == "QCD_MuEnriched":
            if self.reweighting_scheme == "pteta":
                self.pteta_reweighting()
                self.custom_histogram_weights = {k : self.weight_2d for k in self.histograms_to_reweigh}
            elif self.reweighting_scheme == "ptetatau21":
                self.ptetatau21_reweighting()
                self.custom_histogram_weights = {k : self.weight_3d for k in self.histograms_to_reweigh}

    def fill_column_accumulators(self, variation):
        pass
