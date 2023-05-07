import os
import awkward as ak
from collections import defaultdict
import correctionlib

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.leptons import lepton_selection

from config.fatjet_base.custom.leptons import lepton_selection_noniso
from config.fatjet_base.custom.jets import jet_selection
from lib.sv import *
from lib.muon_matching import muons_matched_to_fatjet, muon_matched_to_subjet
from pocket_coffea.lib.hist_manager import Axis

class fatjetBaseProcessor(BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        # Define dictionary to save fatjet JER seeds
        self.output_format.update({"seed_fatjet_chunk": defaultdict(str)})

        # Additional axis for the year
        self.custom_axes.append(
            Axis(
                coll="metadata",
                field="year",
                name="year",
                bins=set(sorted(self.cfg.years)),
                type="strcat",
                growth=False,
                label="Year",
            )
        )

    def apply_object_preselection(self, variation):
        '''
        The ttHbb processor cleans
          - Electrons
          - Muons
          - Jets -> JetGood
          - BJet -> BJetGood

        '''
        # Include the supercluster pseudorapidity variable
        electron_etaSC = self.events.Electron.eta + self.events.Electron.deltaEtaSC
        self.events["Electron"] = ak.with_field(
            self.events.Electron, electron_etaSC, "etaSC"
        )
        # Build masks for selection of muons, electrons, jets, fatjets

        ################################################
        # Dedicated Muon selection for mutag final state
        self.events["MuonGood"] = lepton_selection_noniso(
            self.events, "Muon", self.cfg.finalstate
        )
        ################################################
        self.events["ElectronGood"] = lepton_selection(
            self.events, "Electron", self.cfg.finalstate
        )
        leptons = ak.with_name(
            ak.concatenate((self.events.MuonGood, self.events.ElectronGood), axis=1),
            name='PtEtaPhiMCandidate',
        )
        self.events["LeptonGood"] = leptons[ak.argsort(leptons.pt, ascending=False)]

        # Apply JEC + JER
        #self.apply_JERC()
        self.events["JetGood"], self.jetGoodMask = jet_selection(
            self.events, "Jet", self.cfg.finalstate
        )

        self.events["FatJetGood"], self.fatjetGoodMask = jet_selection(
            self.events, "FatJet", self.cfg.finalstate
        )

        # Select here events with at least one FatJetGood
        self.events = self.events[ak.num(self.events.FatJetGood) >= 1]

        # Uniquely match muons to leading and subleading subjets
        # The shape of these masks is the same as the self.events.FatJetGood collection:
        # the first object are the muons matched to the leading subjet,
        # the second object are the muons that are matched to the subleading subjet, while
        # the third object are the dimuon pairs constructed from the matched muons

        self.events["MuonGoodMatchedToFatJetGood"] = muons_matched_to_fatjet(self.events)
        muon_matched_to_leading_subjet = muon_matched_to_subjet(self.events, pos=0, unique=False)
        muon_matched_to_subleading_subjet = muon_matched_to_subjet(self.events, pos=1, unique=False)
        self.events["MuonGoodMatchedToSubJet"] = ak.concatenate((muon_matched_to_leading_subjet, muon_matched_to_subleading_subjet), axis=2)
        muon_matched_uniquely_to_leading_subjet = muon_matched_to_subjet(self.events, pos=0, unique=True)
        muon_matched_uniquely_to_subleading_subjet = muon_matched_to_subjet(self.events, pos=1, unique=True)
        self.events["MuonGoodMatchedUniquelyToSubJet"] = ak.concatenate((muon_matched_uniquely_to_leading_subjet, muon_matched_uniquely_to_subleading_subjet), axis=2)

        fatjet_fields = {
            "tau21" : self.events.FatJetGood.tau2 / self.events.FatJetGood.tau1,
            #"nSubJet" : ak.count(events.FatJetGood.subjets.pt, axis=2),
            "nMuonGoodMatchedToFatJetGood" : ak.count(self.events["MuonGoodMatchedToFatJetGood"].pt, axis=2),
            "nMuonGoodMatchedToSubJet" : ak.count(self.events["MuonGoodMatchedToSubJet"].pt, axis=2),
            "nMuonGoodMatchedUniquelyToSubJet" : ak.count(self.events["MuonGoodMatchedUniquelyToSubJet"].pt, axis=2)
        }
        for field, value in fatjet_fields.items():
            self.events["FatJetGood"] = ak.with_field(self.events.FatJetGood, value, field)

    def count_objects(self, variation):
        self.events["nMuonGood"] = ak.num(self.events.MuonGood)
        self.events["nElectronGood"] = ak.num(self.events.ElectronGood)
        self.events["nLeptonGood"] = (
            self.events["nMuonGood"] + self.events["nElectronGood"]
        )
        self.events["nJetGood"] = ak.num(self.events.JetGood)
        self.events["nFatJetGood"] = ak.num(self.events.FatJetGood)
        self.events["nSV"] = ak.num(self.events.SV)

    def define_common_variables_after_presel(self, variation):
        
        fatjet_sv_fields = {
            "nsv1"    : get_nsv(self.events.FatJetGood, self.events.SV, pos=0),
            "nsv2"    : get_nsv(self.events.FatJetGood, self.events.SV, pos=1),
        }

        for field, value in fatjet_sv_fields.items():
            self.events = ak.with_field(self.events, value, field)

        Xbb = self.events.FatJetGood.particleNetMD_Xbb
        Xcc = self.events.FatJetGood.particleNetMD_Xcc
        QCD = self.events.FatJetGood.particleNetMD_QCD
        fatjet_fields = {
            "subjet1" : self.events.FatJetGood.subjets[:, :, 0],
            "subjet2" : self.events.FatJetGood.subjets[:, :, 1],
            "particleNetMD_Xbb_QCD" : Xbb / (Xbb + QCD),
            "particleNetMD_Xcc_QCD" : Xcc / (Xcc + QCD),
        }

        for field, value in fatjet_fields.items():
            self.events["FatJetGood"] = ak.with_field(self.events.FatJetGood, value, field)

        self.events.FatJetGood.subjet1 = ak.with_field(self.events.FatJetGood.subjet1, self.events.FatJetGood.subjet1.nearest(self.events.MuonGood, threshold=0.4)[:, 0], "MuonLeading")
        self.events.FatJetGood.subjet2 = ak.with_field(self.events.FatJetGood.subjet2, self.events.FatJetGood.subjet2.nearest(self.events.MuonGood, threshold=0.4)[:, 0], "MuonLeading")

        # Define SV observables

        #self.events.SV = self.events.SV[get_sv_in_jet(self.events.FatJetGood[:,0], self.events.SV)]
        sv_in_jet1 = self.events.SV[get_sv_in_jet(self.events.FatJetGood, self.events.SV, pos=0)]
        sv_in_jet2 = self.events.SV[get_sv_in_jet(self.events.FatJetGood, self.events.SV, pos=1)]
        i1_maxPt     = ak.argsort(sv_in_jet1.pt, ascending=False)
        i2_maxPt     = ak.argsort(sv_in_jet2.pt, ascending=False)
        sv_in_jet1_sorted = sv_in_jet1[i1_maxPt]
        sv_in_jet2_sorted = sv_in_jet2[i2_maxPt]

        sv_fields = {}
        position = {'jet1' : 0, 'jet2' : 1}
        for jet_key, sv_in_jet_sorted in zip(['jet1', 'jet2'], [sv_in_jet1_sorted, sv_in_jet2_sorted]):
            summass, logsummass = get_summass(sv_in_jet_sorted)
            #projmass, logprojmass = get_projmass(self.events.FatJetGood, sv_in_jet_sorted, pos=position[jet_key])
            sv1mass, logsv1mass = get_sv1mass(sv_in_jet_sorted)
            sumcorrmass, logsumcorrmass = get_sumcorrmass(sv_in_jet_sorted)
            sv_fields[jet_key] = {
                "summass" : summass,
                "logsummass" : logsummass,
                #"projmass" : projmass,
                #"logprojmass" : logprojmass,
                "sv1mass" : sv1mass,
                "logsv1mass" : logsv1mass,
                "sumcorrmass" : sumcorrmass,
                "logsumcorrmass" : logsumcorrmass,
            }

        for field in sv_fields['jet1'].keys():
            value1_unflattened = ak.unflatten( sv_fields['jet1'][field], counts=1 )
            value2_unflattened = ak.unflatten( sv_fields['jet2'][field], counts=1 )
            value_concat = ak.concatenate((value1_unflattened, value2_unflattened), axis=1)
            # Pad to FatJetGood dims
            padded = ak.concatenate((value_concat, 
                                     ak.zeros_like(self.events['FatJetGood'].pt)[:, 2:]), axis=1)
            # Clip to FatJetGood size
            padded = padded[ak.local_index(self.events['FatJetGood'].pt)]
            self.events = ak.with_field(self.events, value_concat, field)
            self.events['FatJetGood'] = ak.with_field(self.events['FatJetGood'], padded, field)

    def fill_column_accumulators(self, variation):
        pass
