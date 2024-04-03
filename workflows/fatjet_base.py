import awkward as ak
from collections import defaultdict

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.leptons import lepton_selection

from config.fatjet_base.custom.leptons import lepton_selection_noniso
from config.fatjet_base.custom.jets import jet_selection
from lib.sv import get_corrmass, sv_matched_to_fatjet, get_sumcorrmass, get_sv1mass
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
        # WARNING: Here we are applying a per-event selection asking for at least one AK8 jet in the event
        # In this way, we can compute the number of MuonGood matched to the FatJetGood and its subjets
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
        
        # Correct SV mass for mis-aligment between the SV momentum and the PV-SV direction.
        # This takes into account particles that are not reconstructed in the SV reconstruction.
        # The corrected mass is assigned to the SV.mass branch
        self.events.SV = ak.with_field(self.events.SV, get_corrmass(self.events.SV), "mass")
        # Match SV to AK8 jets
        self.events["SVMatchedToFatJetGood"] = sv_matched_to_fatjet(self.events)

        # Compute final particleNetMD scores
        Xbb = self.events.FatJetGood.particleNetMD_Xbb
        Xcc = self.events.FatJetGood.particleNetMD_Xcc
        QCD = self.events.FatJetGood.particleNetMD_QCD
        sumcorrSVmass, logsumcorrSVmass = get_sumcorrmass(self.events.SVMatchedToFatJetGood)
        # Order SV by dxySig and compute the leading SV mass and its log
        #index_max_pt = ak.argsort(self.events.SVMatchedToFatJetGood.pt, ascending=False)
        index_max_dxySig = ak.argsort(self.events.SVMatchedToFatJetGood.dxySig, ascending=False)
        sv1mass, logsv1mass = get_sv1mass(self.events.SVMatchedToFatJetGood[index_max_dxySig])
        fatjet_fields = {
            "nSVMatchedToFatJetGood": ak.count(self.events["SVMatchedToFatJetGood"].pt, axis=2),
            "particleNetMD_Xbb_QCD" : Xbb / (Xbb + QCD),
            "particleNetMD_Xcc_QCD" : Xcc / (Xcc + QCD),
            "sumcorrSVmass" : sumcorrSVmass,
            "logsumcorrSVmass" : logsumcorrSVmass,
            "sv1mass" : sv1mass,
            "logsv1mass" : logsv1mass,
        }

        for field, value in fatjet_fields.items():
            self.events["FatJetGood"] = ak.with_field(self.events.FatJetGood, value, field)

    def fill_column_accumulators(self, variation):
        pass
