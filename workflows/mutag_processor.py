from collections import defaultdict

import awkward as ak

import correctionlib

from workflows.fatjet_base import fatjetBaseProcessor
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.categorization import StandardSelection
from pocket_coffea.parameters.jet_scale_factors import ptetatau21_reweighting
from lib.sv import *
from config.fatjet_base.custom.cuts import get_ptmsd

class mutagAnalysisProcessor(fatjetBaseProcessor):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        if not "histograms_to_reweigh" in self.cfg.workflow_options.keys():
            raise Exception("The entry of the config file 'workflow_options' does not contain a key 'histograms_to_reweigh'. Please specify it in the config file.")
        self.histograms_to_reweigh = self.cfg.workflow_options["histograms_to_reweigh"]
        self.weight_3d = defaultdict(dict)

    def apply_object_preselection(self, variation, pt_min=350., msd=40.):
        super().apply_object_preselection(variation)

        mask_name = f"pt{int(pt_min)}msd{int(msd)}"
        cuts_fatjet = {mask_name : [get_ptmsd(pt_min, msd)]}
        selection_fatjet = StandardSelection(cuts_fatjet)
        selection_fatjet.prepare(
            events=self.events,
            year=self._year,
            sample=self._sample,
            isMC=self._isMC,
        )
        mask_fatjet = selection_fatjet.get_mask(mask_name)

        # Apply (pt, msd) cuts
        self.events["FatJetGood"] = self.events.FatJetGood[mask_fatjet]

        # Restrict analysis to leading and subleading jets only
        self.events["FatJetGood"] = self.events.FatJetGood[ak.local_index(self.events.FatJetGood, axis=1) < 2]

        # Label leading and subleading AK8 jets BEFORE muon tagging selection
        # Leading: pos=0, Subleading: pos=1
        self.events["FatJetGood"] = ak.with_field(self.events["FatJetGood"], ak.local_index(self.events["FatJetGood"], axis=1), "pos")

    def ptetatau21_reweighting(self, variation):
        '''Correction of jets observable by a 3D reweighting based on (pT, eta, tau21).
        The function stores the nominal, up and down weights in self.weight_3d,
        where the up/down variations are computed considering the statistical uncertainty on data and MC.'''
        cset = correctionlib.CorrectionSet.from_file(ptetatau21_reweighting[self._sample][self._year])
        key = list(cset.keys())[0]
        corr = cset[key]

        cat = "inclusive"
        nfatjet  = ak.num(self.events.FatJetGood.pt)
        pos = ak.flatten(self.events.FatJetGood.pos)
        pt = ak.flatten(self.events.FatJetGood.pt)
        eta = ak.flatten(self.events.FatJetGood.eta)
        tau21 = ak.flatten(self.events.FatJetGood.tau21)

        weight_dict = {"all" : {}, "1" : {}, "2" : {}}
        for var in ["nominal", "statUp", "statDown"]:
            w = corr.evaluate(cat, variation, var, pos, pt, eta, tau21)
            weight_dict["all"][var] = ak.unflatten(w, nfatjet)
            # Here we build the flattened custom weights for the leading and subleading jet collections.
            # In order for the length of the weights array to match the number of the per-event mask,
            # we set the weight to be 1 for the events that does not contain a jet with pos=0(1)
            weight_dict["1"][var] = ak.fill_none(ak.firsts(weight_dict["all"]["nominal"][self.events.FatJetGood.pos == 0]), 1)
            weight_dict["2"][var] = ak.fill_none(ak.firsts(weight_dict["all"]["nominal"][self.events.FatJetGood.pos == 1]), 1)

        # Here we store the dictionary for the custom weights, with the following content:
        # - pos = "all": the jagged array of weights of the full jet collection
        # - pos = "1": the flat array of weights of the leading jet collection
        # - pos = "2": the flat array of weights of the subleading jet collection
        # For each of the 3 cases, the nominal, up and down variations are stored.
        for pos, weight in weight_dict.items():
            self.weight_3d[pos] = {
                "nominal" : weight["nominal"],
                "sf_ptetatau21_reweightingUp" : weight["statUp"],
                "sf_ptetatau21_reweightingDown" : weight["statDown"]
            }

    def process_extra_after_presel(self, variation):
        if self._sample in ["QCD_MuEnriched", "QCD_HT"]:
            self.ptetatau21_reweighting(variation)
            for pos, hists in self.histograms_to_reweigh["by_pos"].items():
                for histname in hists:
                    self.custom_histogram_weights[histname] = self.weight_3d[pos]

    def fill_column_accumulators(self, variation):
        pass
