import os
import tarfile

import numpy as np
import awkward as ak
import uproot
import coffea
from coffea import hist, processor, lookup_tools
from coffea.util import load
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory

from lib.luminosity import rescale
from lib.sv import get_nsv, get_sv_in_jet
from parameters import JECversions, jecTarFiles, FinalMask, PtBinning, AK8Taggers, AK8TaggerWP

class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self, cfg):
        self.cfg = cfg
        self._year = self.cfg['year']
        self._campaign = self.cfg['campaign']
        self._mask_fatjets = {
            'basic'       : {
                'pt_cut' : 350.,
                'eta_cut': 2.4,
                'jetId_cut': 2,
                'mass_cut' : 20.,
                'tau21_cut' : 1.1,
                'DDCvB_cut' : -1
                    },
            'msd100tau06'       : {
                'pt_cut' : 350.,
                'eta_cut': 2.4,
                'jetId_cut': 2,
                'mass_cut' : 100.,
                'tau21_cut' : 0.6,
                'DDCvB_cut' : -1
                    },
            'msd100tau06ggHcc'       : {
                'pt_cut' : 350.,
                'eta_cut': 2.4,
                'jetId_cut': 2,
                'mass_cut' : 100.,
                'tau21_cut' : 0.6,
                'DDCvB_cut' : 0.03
                    },
        }
        self._final_mask = FinalMask
        #self._final_mask = ['msd100tau06', 'msd100tau03']
        self._AK8TaggerWP = AK8TaggerWP[self._campaign][self._year]
        self._PtBinning = PtBinning[self._campaign][self._year]
        self.mupt = cfg['mupt']
        self.corrJECfolder = cfg['JECfolder']
        self.hist2d = self.cfg['hist2d']
        self.checkOverlap = self.cfg['checkOverlap']
        if self.checkOverlap:
            self.eventTags = {'run' : None, 'lumi' : None, 'event' : None}

        ############


        ##############
        # Trigger level
        self.triggers = [
            "HLT_BTagMu_AK8Jet300_Mu5",
            "HLT_BTagMu_AK4Jet300_Mu5",
        ]

        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        flavor_axis  = hist.Cat("flavor",   "Flavor")

        # Events
        #nel_axis     = hist.Bin("nel",   r"N electrons",     [0,1,2,3,4,5,6,7,8,9,10])
        #nmu_axis     = hist.Bin("nmu",   r"N muons",         [0,1,2,3,4,5,6,7,8,9,10])
        #njet_axis    = hist.Bin("njet",  r"N jets",          [0,1,2,3,4,5,6,7,8,9,10])
        #nbjet_axis   = hist.Bin("nbjet", r"N b-jets",        [0,1,2,3,4,5,6,7,8,9,10])
        nfatjet_axis = hist.Bin("nfatjet",  r"N fatjets",    [0,1,2,3,4,5,6,7,8,9,10])
        nmusj1_axis  = hist.Bin("nmusj1",  r"$N_{mu}$(sj1)", 30, 0, 30)
        nmusj2_axis  = hist.Bin("nmusj2",  r"$N_{mu}$(sj2)", 30, 0, 30)
        nsv1_axis    = hist.Bin("nsv1",  r"$N_{SV}$(sj1)",   30, 0, 30)
        nsv2_axis    = hist.Bin("nsv2",  r"$N_{SV}$(sj2)",   30, 0, 30)

        # Muon
        leadmuon_pt_axis   = hist.Bin("pt",   r"lead. Muon $p_{T}$ [GeV]", 200, 0, 200)
        #leadmuon_eta_axis  = hist.Bin("eta",  r"lead. Muon $\eta$", 60, -3, 3)
        #leadmuon_phi_axis  = hist.Bin("phi",  r"lead. Muon $\phi$", 60, -np.pi, np.pi)
        subleadmuon_pt_axis   = hist.Bin("pt",   r"sublead. Muon $p_{T}$ [GeV]", 200, 0, 200)
        #subleadmuon_eta_axis  = hist.Bin("eta",  r"sublead. Muon $\eta$", 60, -3, 3)
        #subleadmuon_phi_axis  = hist.Bin("phi",  r"sublead. Muon $\phi$", 60, -np.pi, np.pi)
        leadmuonsj1_pt_axis   = hist.Bin("pt",   r"lead. Muon (sj1) $p_{T}$ [GeV]", 200, 0, 200)
        leadmuonsj2_pt_axis   = hist.Bin("pt",   r"lead. Muon (sj2) $p_{T}$ [GeV]", 200, 0, 200)

        # Jet
        #jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 100, 20, 400)
        #jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 60, -3, 3)
        #jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 60, -3, 3)
        #jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 100, 0, 50)
        #ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)

        # FatJet
        fatjet_tau21_axis = hist.Bin("tau21", r"lead. FatJet $\tau_{21}$", 50, 0, 1)
        fatjet_n2b1_axis  = hist.Bin("n2b1", r"lead. FatJet $N_{2}^{(\beta=1)}$", 50, 0, 0.5)
        fatjet_pt_axis    = hist.Bin("pt",   r"lead. FatJet $p_{T}$ [GeV]", 600, 0, 3000)
        fatjet_eta_axis   = hist.Bin("eta",  r"lead. FatJet $\eta$", 60, -3, 3)
        fatjet_phi_axis   = hist.Bin("phi",  r"lead. FatJet $\phi$", 60, -np.pi, np.pi)
        fatjet_mass_axis  = hist.Bin("mass", r"lead. FatJet $m_{SD}$ [GeV]", 1000, 0, 1000)
        fatjet_jetproba_axis = hist.Bin("Proba", r"lead. FatJet JP", 50, 0, 2.5)
        #fatjet_vertexmass_axis  = hist.Bin("vertexmass", r"lead. FatJet tau1 vertex $m_{SD}$ [GeV]", 1000, 0, 1000)

        # SV
        sv_sv1mass_axis             = hist.Bin("sv1mass", r"lead. FatJet $m_{SV,1}$ [GeV]", 1000, 0, 1000)
        sv_logsv1mass_axis          = hist.Bin("logsv1mass", r"lead. FatJet log($m_{SV,1}$/GeV)", 80, -4, 4)
        sv_sv1mass_maxdxySig_axis    = hist.Bin("sv1mass_maxdxySig", r"lead. FatJet $m_{SV,1~for~max(\sigma_{d_{xy}})}$ [GeV]", 1000, 0, 1000)
        sv_logsv1mass_maxdxySig_axis = hist.Bin("logsv1mass_maxdxySig", r"lead. FatJet log($m_{SV,1~for~max(\sigma_{d_{xy}})}$/GeV)", 80, -4, 4)
        sv_logsv1massratio_axis     = hist.Bin("logsv1massratio", r"log($m_{SV_1~for~max(\sigma_{d_{xy}})}$/GeV) / log($m_{SV,1~for~max(p_T)}$/GeV)", 200, -100, 100)
        sv_logsv1massres_axis     = hist.Bin("logsv1massres", r"(log($m_{SV_1~for~max(\sigma_{d_{xy}})}$/GeV) - log($m_{SV,1~for~max(p_T)}$/GeV)) / log($m_{SV,1~for~max(p_T)}$/GeV))", 100, -1, 1)

        # Define similar axes dynamically
        disc_list = ["btagCMVA", "btagCSVV2", 'btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC',]
        disc_list_fj = AK8Taggers
        btag_axes = []
        btag_axes_fj = []
        for d in disc_list:
            btag_axes.append(hist.Bin(d, d, 40, 0, 1))
        for d in disc_list_fj:
            btag_axes_fj.append(hist.Bin(d, d, 40, 0, 1))

        # Define histograms from axes
        _hist_muon_dict = {
                'leadmuon_pt'  : hist.Hist("Events", dataset_axis, flavor_axis, leadmuon_pt_axis),
                #'leadmuon_eta' : hist.Hist("Events", dataset_axis, flavor_axis, leadmuon_eta_axis),
                #'leadmuon_phi' : hist.Hist("Events", dataset_axis, flavor_axis, leadmuon_phi_axis),
                'subleadmuon_pt'  : hist.Hist("Events", dataset_axis, flavor_axis, subleadmuon_pt_axis),
                #'subleadmuon_eta' : hist.Hist("Events", dataset_axis, flavor_axis, subleadmuon_eta_axis),
                #'subleadmuon_phi' : hist.Hist("Events", dataset_axis, flavor_axis, subleadmuon_phi_axis),
                'leadmuonsj1_pt'  : hist.Hist("Events", dataset_axis, flavor_axis, leadmuonsj1_pt_axis),
                'leadmuonsj2_pt'  : hist.Hist("Events", dataset_axis, flavor_axis, leadmuonsj2_pt_axis),
            }

        #_hist_jet_dict = {
        #        'jet_pt'  : hist.Hist("Events", dataset_axis, jet_pt_axis),
        #        'jet_eta' : hist.Hist("Events", dataset_axis, jet_eta_axis),
        #        'jet_phi' : hist.Hist("Events", dataset_axis, jet_phi_axis),
        #        'jet_mass': hist.Hist("Events", dataset_axis, jet_mass_axis),
        #    }

        _hist_fatjet_dict = {
                'fatjet_tau21' : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_tau21_axis),
                'fatjet_n2b1'  : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_n2b1_axis),
                'fatjet_pt'  : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_pt_axis),
                'fatjet_eta' : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_eta_axis),
                'fatjet_phi' : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_phi_axis),
                'fatjet_mass': hist.Hist("Events", dataset_axis, flavor_axis, fatjet_mass_axis),
                'fatjet_nmusj1' : hist.Hist("Events", dataset_axis, flavor_axis, nmusj1_axis),
                'fatjet_nmusj2' : hist.Hist("Events", dataset_axis, flavor_axis, nmusj2_axis),
                'fatjet_nsv1'   : hist.Hist("Events", dataset_axis, flavor_axis, nsv1_axis),
                'fatjet_nsv2'   : hist.Hist("Events", dataset_axis, flavor_axis, nsv2_axis),
                'fatjet_jetproba' : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_jetproba_axis),
                #'fatjet_DDX_tau1_vertexMass' : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_vertexmass_axis),
            }
        _hist_sv_dict = {
                'sv_sv1mass'              : hist.Hist("Events", dataset_axis, flavor_axis, sv_sv1mass_axis),
                'sv_logsv1mass'           : hist.Hist("Events", dataset_axis, flavor_axis, sv_logsv1mass_axis),
                'sv_sv1mass_maxdxySig'    : hist.Hist("Events", dataset_axis, flavor_axis, sv_sv1mass_maxdxySig_axis),
                'sv_logsv1mass_maxdxySig' : hist.Hist("Events", dataset_axis, flavor_axis, sv_logsv1mass_maxdxySig_axis),
                'sv_logsv1massratio'      : hist.Hist("Events", dataset_axis, flavor_axis, sv_logsv1massratio_axis),
                'sv_logsv1massres'        : hist.Hist("Events", dataset_axis, flavor_axis, sv_logsv1massres_axis),
            }

        for (i, disc) in enumerate(disc_list_fj):
            _hist_fatjet_dict['fatjet_' + disc] = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i])

        # Define 2D histograms
        if self.hist2d:
            _hist2d_dict = {
                'hist2d_sv_logsv1mass_maxdxySig_vs_maxPt' : hist.Hist("Events", dataset_axis, flavor_axis, sv_logsv1mass_maxdxySig_axis, sv_logsv1mass_axis),
                'hist2d_sv_logsv1massratio_vs_maxPt'      : hist.Hist("Events", dataset_axis, flavor_axis, sv_logsv1massratio_axis, sv_logsv1mass_axis),
                'hist2d_sv_logsv1massres_vs_maxPt'        : hist.Hist("Events", dataset_axis, flavor_axis, sv_logsv1massres_axis, sv_logsv1mass_axis),
            }
            for (i, disc) in enumerate(disc_list_fj):
                _hist2d_dict['hist2d_fatjet_pt_vs_' + disc]    = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], fatjet_pt_axis)
                _hist2d_dict['hist2d_fatjet_mass_vs_' + disc]  = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], fatjet_mass_axis)
                _hist2d_dict['hist2d_fatjet_tau21_vs_' + disc] = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], fatjet_tau21_axis)
                _hist2d_dict['hist2d_fatjet_n2b1_vs_' + disc]  = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], fatjet_n2b1_axis)
                _hist2d_dict['hist2d_fatjet_nsv1_vs_' + disc]         = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], nsv1_axis)
                _hist2d_dict['hist2d_fatjet_nsv2_vs_' + disc]         = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], nsv2_axis)
                _hist2d_dict['hist2d_fatjet_nmusj1_vs_' + disc]       = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], nmusj1_axis)
                _hist2d_dict['hist2d_fatjet_nmusj2_vs_' + disc]       = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], nmusj2_axis)

        _hist_event_dict = {
                #'njet'   : hist.Hist("Events", dataset_axis, njet_axis),
                #'nbjet'  : hist.Hist("Events", dataset_axis, nbjet_axis),
                #'nel'    : hist.Hist("Events", dataset_axis, nel_axis),
                #'nmu'    : hist.Hist("Events", dataset_axis, nmu_axis),
                'nfatjet': hist.Hist("Events", dataset_axis, flavor_axis, nfatjet_axis),
            }
        _sumw_dict = {'sumw': processor.defaultdict_accumulator(float),
                      'nbtagmu': processor.defaultdict_accumulator(float),
            }

        self.muon_hists = list(_hist_muon_dict.keys())
        #self.jet_hists = list(_hist_jet_dict.keys())
        self.fatjet_hists = list(_hist_fatjet_dict.keys())
        self.sv_hists = list(_hist_sv_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())

        #_hist_dict = {**_hist_jet_dict, **_hist_fatjet_dict, **_hist2d_dict, **_hist_event_dict, **_sumw_dict}
        if self.hist2d:
            self._hist_dict = {**_hist_muon_dict, **_hist_fatjet_dict, **_hist_sv_dict, **_hist2d_dict, **_hist_event_dict}
        else:
            self._hist_dict = {**_hist_muon_dict, **_hist_fatjet_dict, **_hist_sv_dict, **_hist_event_dict}
        self.append_mask()
        self._hist_dict.update({**_sumw_dict})
        self._accumulator = processor.dict_accumulator(self._hist_dict)
        if self.checkOverlap:
            for var in self.eventTags.keys():
                self._accumulator.add(processor.dict_accumulator({var : processor.column_accumulator(np.array([]))}))

    @property
    def accumulator(self):
        return self._accumulator

    # Function to load year-dependent parameters
    def load_metadata(self):
        self._dataset = self.events.metadata["dataset"]
        self._sample = self.events.metadata["sample"]
        self._year = self.events.metadata["year"]
        self._campaign = self.events.metadata["campaign"]
        if (self._campaign == 'UL') & (self._year == '2016'):
            self._VFP = self.events.metadata["VFP"]

    def load_era_specific_parameters(self):
        # JEC files
        # Correction files in https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
        self.jesInputFilePath = os.getcwd()+"/correction_files/tmp"
        if not os.path.exists(self.jesInputFilePath):
            os.makedirs(self.jesInputFilePath)
        if (self._campaign == 'UL') & (self._year == '2016'):
            files = jecTarFiles[self._campaign][f"{self._year}_{self._VFP}"]
        else:
            files = jecTarFiles[self._campaign][self._year]
        for itar in files:
            jecFile = os.getcwd()+itar
            jesArchive = tarfile.open( jecFile, "r:gz")
            jesArchive.extractall(self.jesInputFilePath)

        # PU files
        self.puFile    = self.cfg['puFile']
        self.nTrueFile = self.cfg['nTrueFile']

    def append_mask(self):
        masks = list(self._mask_fatjets.keys())
        d = {}
        for histname in self._hist_dict.keys():
            h = self._hist_dict[histname]
            d[f'{histname}_{masks[0]}'] = h
            for maskname in masks[1:]:
                d[f'{histname}_{maskname}'] = h.copy()
                if maskname in self._final_mask:
                    for tagger in self._AK8TaggerWP.keys():
                        for wp in self._AK8TaggerWP[tagger].keys():
                            d[f'{histname}_{maskname}{tagger}pass{wp}wp'] = h.copy()
                            d[f'{histname}_{maskname}{tagger}fail{wp}wp'] = h.copy()
                            for wpt in self._PtBinning.keys():
                                pt_low, pt_high = self._PtBinning[wpt]
                                pt_low, pt_high = (str(pt_low), str(pt_high))
                                d[f'{histname}_{maskname}{tagger}pass{wp}wpPt-{pt_low}to{pt_high}'] = h.copy()
                                d[f'{histname}_{maskname}{tagger}fail{wp}wpPt-{pt_low}to{pt_high}'] = h.copy()
        self._hist_dict = d.copy()

        #for attr, hists in zip(["fatjet_hists", "sv_hists", "event_hists"], [self.fatjet_hists, self.sv_hists, self.event_hists]):
        for attr in ["muon_hists", "fatjet_hists", "sv_hists", "event_hists"]:
            attr_updated = []
            for histname in getattr(self, attr):
                for maskname in masks:
                    attr_updated.append(f'{histname}_{maskname}')
                    if maskname in self._final_mask:
                        for tagger in self._AK8TaggerWP.keys():
                            for wp in self._AK8TaggerWP[tagger].keys():
                                attr_updated.append(f'{histname}_{maskname}{tagger}pass{wp}wp')
                                attr_updated.append(f'{histname}_{maskname}{tagger}fail{wp}wp')
                                for wpt in self._PtBinning.keys():
                                    pt_low, pt_high = self._PtBinning[wpt]
                                    pt_low, pt_high = (str(pt_low), str(pt_high))
                                    attr_updated.append(f'{histname}_{maskname}{tagger}pass{wp}wpPt-{pt_low}to{pt_high}')
                                    attr_updated.append(f'{histname}_{maskname}{tagger}fail{wp}wpPt-{pt_low}to{pt_high}')
            setattr(self, attr, attr_updated)

        return self._hist_dict

    def puReweight(self, puFile, nTrueFile, dataset ):
        '''Based on https://github.com/andrzejnovak/coffeandbacon/blob/master/analysis/compile_corrections.py#L166-L192'''

        nTrueIntLoad = load(nTrueFile)
        print([y for x,y in nTrueIntLoad[dataset].sum('dataset').values().items()])
        nTrueInt = [y for x,y in nTrueIntLoad[dataset].sum('dataset').values().items()][0]  ## not sure is the best way

        with uproot.open(puFile) as file_pu:
            norm = lambda x: x / x.sum()
            data = norm(file_pu['pileup'].counts())
            mc_pu = norm(nTrueInt)
            mask = mc_pu > 0.
            corr = data.copy()
            corr[mask] /= mc_pu[mask]
            pileup_corr = lookup_tools.dense_lookup.dense_lookup(corr, file_pu["pileup"].axis().edges())
        return pileup_corr

    def applyJEC( self, jets, fixedGridRhoFastjetAll, events_cache, typeJet, isData, JECversion ):
        '''Based on https://coffeateam.github.io/coffea/notebooks/applying_corrections.html#Applying-energy-scale-transformations-to-Jets'''

        ext = lookup_tools.extractor()
        JECtypes = [ 'L1FastJet', 'L2Relative', 'L2Residual', 'L3Absolute', 'L2L3Residual' ]
        jec_stack_names = [ JECversion+'_'+k+'_'+typeJet for k in JECtypes ]
        JECtypesfiles = [ '* * '+self.corrJECfolder+'/'+k+'.txt' for k in jec_stack_names ]
        ext.add_weight_sets( JECtypesfiles )
        ext.finalize()
        evaluator = ext.make_evaluator()

        #print("available evaluator keys:")
        #for key in evaluator.keys():
        #    print("\t", key)

        jec_inputs = {name: evaluator[name] for name in jec_stack_names}
        corrector = FactorizedJetCorrector( **jec_inputs )
        for i in jec_inputs: print(i,'\n',evaluator[i])

        #print(dir(evaluator))
        #print()
        jec_stack = JECStack(jec_inputs)
        name_map = jec_stack.blank_name_map
        name_map['JetPt'] = 'pt'
        name_map['JetMass'] = 'mass'
        name_map['JetEta'] = 'eta'
        name_map['JetA'] = 'area'

        jets['pt_raw'] = (1 - jets['rawFactor']) * jets['pt']
        jets['mass_raw'] = (1 - jets['rawFactor']) * jets['mass']
        jets['rho'] = ak.broadcast_arrays(fixedGridRhoFastjetAll, jets.pt)[0]
        name_map['ptRaw'] = 'pt_raw'
        name_map['massRaw'] = 'mass_raw'
        name_map['Rho'] = 'rho'
        if not isData:
            jets['pt_gen'] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
            name_map['ptGenJet'] = 'pt_gen'


        jet_factory = CorrectedJetsFactory(name_map, jec_stack)
        corrected_jets = jet_factory.build(jets, lazy_cache=events_cache)
        #print()
        #print('starting columns:',ak.fields(jets))
        #print()

        #print('untransformed pt ratios',jets.pt/jets.pt_raw)
        #print('untransformed mass ratios',jets.mass/jets.mass_raw)

        #print('transformed pt ratios',corrected_jets.pt/corrected_jets.pt_raw)
        #print('transformed mass ratios',corrected_jets.mass/corrected_jets.mass_raw)

        #print()
        #print('transformed columns:', ak.fields(corrected_jets))
        return corrected_jets

    def process(self, events):
        output = self.accumulator.identity()
        if len(events) == 0: return output

        self.events = events
        self.load_metadata()
        self.load_era_specific_parameters()

        isRealData = 'genWeight' not in events.fields
        if not isRealData:
            output['sumw'][self._sample] += sum(events.genWeight)
            if (self._campaign == 'UL') & (self._year == '2016'):
                JECversion = JECversions[self._campaign][f"{self._year}_{self._VFP}"]['MC']
            else:
                JECversion = JECversions[self._campaign][self._year]['MC']
        else:
            output['nbtagmu'][self._sample] += ak.count(events.event)
            if (self._campaign == 'UL') & (self._year == '2016'):
                JECversion = JECversions[self._campaign][f"{self._year}_{self._VFP}"]['Data'][self._sample.split('BTagMu')[1]]
            else:
                JECversion = JECversions[self._campaign][self._year]['Data'][self._sample.split('BTagMu')[1]]

        ############
        # Basic Cleaning
        events = events[ events.PV.npvsGood>0 ]
        METFilters = [ 'goodVertices','globalSuperTightHalo2016Filter', 'HBHENoiseFilter', 'HBHENoiseIsoFilter', 'EcalDeadCellTriggerPrimitiveFilter', 'BadPFMuonFilter' ]
        if self._campaign.startswith('UL'):
            if self._year in ['2016']:
                METFilters = METFilters + [ 'BadPFMuonDzFilter', 'eeBadScFilter', 'ecalBadCalibFilter']
            elif self._year in ['2017', '2018']:
                METFilters = METFilters + [ 'eeBadScFilter', 'ecalBadCalibFilter']
        if self._campaign.startswith('EOY') and isRealData: METFilters.append('eeBadScFilter')
        for imet in METFilters: events = events[ getattr( events.Flag, imet )==1 ]

        ############
        # Some corrections
        weights = processor.Weights(len(events))
        corrections = {}
        if not isRealData:
            weights.add( 'genWeight', events.genWeight)
            weights.add( 'pileup_weight', self.puReweight( self.puFile, self.nTrueFile, self._dataset )( events.Pileup.nPU )  )

        events.FatJet = self.applyJEC( events.FatJet, events.fixedGridRhoFastjetAll, events.caches[0], 'AK8PFPuppi', isRealData, JECversion )

        #cuts = {}
        #for tagger in self._AK8TaggerWP.keys():
        #    cuts[tagger] = processor.PackedSelection()
        cuts = processor.PackedSelection()

        ############
        # Trigger selection
        if self._year == '2016':
            if 'BTagMu_AK4Jet300_Mu5' not in events.HLT.fields:
                self.triggers = [trigger.replace('AK4', '') for trigger in self.triggers]
            if 'BTagMu_AK8Jet300_Mu5' not in events.HLT.fields:
                self.triggers = [trigger.replace('AK8', '') for trigger in self.triggers]
            #print("self.triggers", self.triggers)
            #print("events.HLT.fields", [item for item in events.HLT.fields if 'BTagMu' in item])
        elif self._year == '2018':
            for (i, trigger) in enumerate(self.triggers):
                if trigger.strip("HLT_") not in events.HLT.fields:
                    self.triggers[i] = trigger + "_noalgo"

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in self.triggers]
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t
        cuts.add('trigger', ak.to_numpy(req_trig))

        ############
        # Basic cuts
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[(events.Muon.pt > self.mupt) & (abs(events.Muon.eta < 2.4)) & (events.Muon.tightId != 1) & (events.Muon.pfRelIso04_all > 0.15)]
        events.Muon = ak.pad_none(events.Muon, 2, axis=1)

        ## Jet cuts  (not used)
        events.Jet = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.5)]
        #req_jets = (ak.count(events.Jet.pt, axis=1) >= 2)

        ## FatJet cuts
        events.FatJet = events.FatJet[(events.FatJet.pt > self._mask_fatjets['basic']['pt_cut']) & (abs(events.FatJet.eta) <= self._mask_fatjets['basic']['eta_cut']) & (events.FatJet.jetId > self._mask_fatjets['basic']['jetId_cut'])  & (ak.count(events.FatJet.subjets.pt, axis=2) >= 2) ]  ## subjet sel to crosscheck
        #print(events['FatJetSVs'])

        ## Event level variables
        eventVariables = {}
        eventVariables['nfatjet'] = ak.num(events.FatJet)

        ## Leading/subleading muon variables
        leadmuon = events.Muon[:, 0]
        subleadmuon = events.Muon[:, 1]

        ## Leading jet variables
        leadfatjet = ak.firsts(events.FatJet)
        leadfatjet['tau21'] = leadfatjet.tau2 / leadfatjet.tau1
        subjet1 = ak.pad_none(leadfatjet.subjets, 2)[:, 0]
        subjet2 = ak.pad_none(leadfatjet.subjets, 2)[:, 1]
        leadfatjet['nsv1'] = get_nsv( subjet1, events.SV )
        leadfatjet['nsv2'] = get_nsv( subjet2, events.SV )
        leadfatjet['nmusj1'] = ak.num(subjet1.delta_r(events.Muon) < 0.4)
        leadfatjet['nmusj2'] = ak.num(subjet2.delta_r(events.Muon) < 0.4)
        leadmuonsj1 = ak.pad_none(events.Muon[subjet1.delta_r(events.Muon) < 0.4], 1)[:,0]
        leadmuonsj2 = ak.pad_none(events.Muon[subjet2.delta_r(events.Muon) < 0.4], 1)[:,0]
        muonCollection = {'leadmuon' : leadmuon, 'subleadmuon' : subleadmuon,
                          'leadmuonsj1' : leadmuonsj1, 'leadmuonsj2' : leadmuonsj2,}

        events.SV = events.SV[get_sv_in_jet(leadfatjet, events.SV)]
        i_maxPt     = ak.argsort(events.SV.pt, ascending=False)
        i_maxdxySig = ak.argsort(events.SV.dxySig, ascending=False)

        try: events.SV[i_maxPt]
        except: return output

        leadsv = ak.firsts(events.SV[i_maxPt])
        leadsv_dxySig = ak.firsts(events.SV[i_maxdxySig])
        leadsv['sv1mass'] = leadsv.mass
        leadsv['logsv1mass'] = np.log(leadsv.mass)
        leadsv['sv1mass_maxdxySig'] = leadsv_dxySig.mass
        leadsv['logsv1mass_maxdxySig'] = np.log(leadsv_dxySig.mass)
        leadsv['logsv1massratio'] = leadsv['logsv1mass_maxdxySig'] / leadsv['logsv1mass']
        leadsv['logsv1massres'] = (leadsv['logsv1mass_maxdxySig'] - leadsv['logsv1mass']) / leadsv['logsv1mass']

        fatjet_mutag = (leadfatjet.nmusj1 >= 1) & (leadfatjet.nmusj2 >= 1)
        cuts.add( 'fatjet_mutag', ak.to_numpy(fatjet_mutag) )

        for tagger in self._AK8TaggerWP.keys():
            for wp, (cut_low, cut_high) in self._AK8TaggerWP[tagger].items():
                tag_pass = (leadfatjet[tagger] > cut_low) & (leadfatjet[tagger] <= cut_high)
                tag_fail = ~tag_pass & (leadfatjet[tagger] >= 0) & (leadfatjet[tagger] <= 1)
                cuts.add( f'{tagger}pass{wp}wp', ak.to_numpy(tag_pass) )
                cuts.add( f'{tagger}fail{wp}wp', ak.to_numpy(tag_fail) )
        for wpt, (pt_low, pt_high) in self._PtBinning.items():
            if pt_high == 'Inf':
                tag_pt = (leadfatjet.pt >= pt_low)
            else:
                tag_pt = (leadfatjet.pt >= pt_low) & (leadfatjet.pt < pt_high)
            cuts.add( f'Pt-{pt_low}to{pt_high}', ak.to_numpy(tag_pt) )

        flavors = {}
        if not isRealData:
            flavors['b'] = (leadfatjet.hadronFlavour == 5)
            flavors['c'] = (leadfatjet.hadronFlavour == 4)
            flavors['l'] = (leadfatjet.hadronFlavour < 4)
            flavors['bb'] = abs(leadfatjet.hadronFlavour == 5) & (leadfatjet.nBHadrons >= 2) #& (leadfatjet.nCHadrons == 0)
            flavors['cc'] = abs(leadfatjet.hadronFlavour == 4) & (leadfatjet.nBHadrons == 0) & (leadfatjet.nCHadrons >= 2)
            #flavors['ll'] = abs(leadfatjet.hadronFlavour < 4) & (leadfatjet.nBHadrons == 0) & (leadfatjet.nCHadrons == 0)
            flavors['b'] = flavors['b'] & ~flavors['bb']
            flavors['c'] = flavors['c'] & ~flavors['cc']
            flavors['l'] = flavors['l'] & ~flavors['bb'] & ~flavors['cc'] & ~flavors['b'] & ~flavors['c']
            #flavors['others'] = ~flavors['l'] & ~flavors['bb'] & ~flavors['cc'] & ~flavors['b'] & ~flavors['c']
        else:
            flavors['Data'] = np.ones(len(events), dtype='bool')

        for selname, cut in self._mask_fatjets.items():

            sel = (leadfatjet.pt > cut['pt_cut']) & \
                    (leadfatjet.msoftdrop > cut['mass_cut']) & \
                    (abs(leadfatjet.eta) < cut['eta_cut']) & \
                    (leadfatjet.jetId >= cut['jetId_cut']) & \
                    (leadfatjet.tau21 < cut['tau21_cut']) & \
                    (leadfatjet.btagDDCvBV2 > cut['DDCvB_cut'])

            cuts.add( selname, ak.to_numpy( sel ) )

        selection = {}
        selection['basic'] = { 'trigger', 'basic' }
        selection['msd100tau06'] = { 'trigger', 'fatjet_mutag', 'msd100tau06' }
        selection['msd100tau06ggHcc'] = { 'trigger', 'fatjet_mutag', 'msd100tau06ggHcc' }
        for mask_f in self._final_mask:
            for tagger in self._AK8TaggerWP.keys():
                for wp, cut in self._AK8TaggerWP[tagger].items():
                    selection[f'{mask_f}{tagger}pass{wp}wp'] = selection[mask_f].copy()
                    selection[f'{mask_f}{tagger}pass{wp}wp'].add(f'{tagger}pass{wp}wp')
                    selection[f'{mask_f}{tagger}fail{wp}wp'] = selection[mask_f].copy()
                    selection[f'{mask_f}{tagger}fail{wp}wp'].add(f'{tagger}fail{wp}wp')
                    for wpt, (pt_low, pt_high) in self._PtBinning.items():
                        #selection[f'{mask_f}{tagger}pass{wp}wpPt-{pt_low}to{pt_high}'] = selection[mask_f].copy()
                        selection[f'{mask_f}{tagger}pass{wp}wpPt-{pt_low}to{pt_high}'] = selection[f'{mask_f}{tagger}pass{wp}wp'].copy()
                        selection[f'{mask_f}{tagger}pass{wp}wpPt-{pt_low}to{pt_high}'].add(f'Pt-{pt_low}to{pt_high}')
                        selection[f'{mask_f}{tagger}fail{wp}wpPt-{pt_low}to{pt_high}'] = selection[f'{mask_f}{tagger}fail{wp}wp'].copy()
                        selection[f'{mask_f}{tagger}fail{wp}wpPt-{pt_low}to{pt_high}'].add(f'Pt-{pt_low}to{pt_high}')


        for histname, h in output.items():
            sel = [mask + histname.split(mask)[-1] for mask in self._mask_fatjets.keys() if mask in histname]
            if histname in self.muon_hists:
                #print('array printout:', np.array(list(muonCollection.keys())))
                #print('list printout:', np.array(list(muonCollection.keys()))[[k in histname for k in muonCollection.keys()]])
                muonKey = np.array(list(muonCollection.keys()))[[k in histname for k in muonCollection.keys()]][0]
                muon = muonCollection[muonKey]
                for flav, mask in flavors.items():
                    weight = weights.weight() * cuts.all(*selection[sel[0]]) * ak.to_numpy(mask)
                    fields = {k: ak.fill_none(muon[k], -9999) for k in h.fields if k in dir(muon)}
                    h.fill(dataset=self._sample, flavor=flav, **fields, weight=weight)
            if ((histname in self.fatjet_hists) | ('hist2d_fatjet' in histname)):
                for flav, mask in flavors.items():
                    weight = weights.weight() * cuts.all(*selection[sel[0]]) * ak.to_numpy(mask)
                    fields = {k: ak.fill_none(leadfatjet[k], -9999) for k in h.fields if k in dir(leadfatjet)}
                    h.fill(dataset=self._sample, flavor=flav, **fields, weight=weight)
            if histname in self.event_hists:
                for flav, mask in flavors.items():
                    weight = weights.weight() * cuts.all(*selection[sel[0]]) * ak.to_numpy(mask)
                    fields = {k: ak.fill_none(eventVariables[k], -9999) for k in h.fields if k in eventVariables.keys() }
                    h.fill(dataset=self._sample, flavor=flav, **fields, weight=weight)
            if ((histname in self.sv_hists) | ('hist2d_sv' in histname)):
                for flav, mask in flavors.items():
                    weight = weights.weight() * cuts.all(*selection[sel[0]]) * ak.to_numpy(mask)
                    fields = {k: ak.fill_none(leadsv[k], -9999) for k in h.fields if k in dir(leadsv) }
                    h.fill(dataset=self._sample, flavor=flav, **fields, weight=weight)

        #if isRealData & (self.checkOverlap is not None):
        if self.checkOverlap:
            mask = self._final_mask[0]
            self.eventTags['run'] = events.run[cuts.all(*selection[mask])]
            self.eventTags['lumi'] = events.luminosityBlock[cuts.all(*selection[mask])]
            self.eventTags['event'] = events.event[cuts.all(*selection[mask])]
            for var in self.eventTags.keys():
                output[var] = output[var] + processor.column_accumulator(ak.to_numpy(self.eventTags[var]))

        return output

    def postprocess(self, accumulator):

        if self.checkOverlap:
            mask = self._final_mask[0]
            self.checkOverlap = self.checkOverlap.replace('.txt', f'_{mask}.txt')
            run = accumulator['run'].value
            lumi = accumulator['lumi'].value
            event = accumulator['event'].value
            #print(run)
            with open(self.checkOverlap, 'w') as file:
                for (iev,r) in enumerate(run):
                    if r==1:
                        continue
                    else:
                        file.write(f'{int(run[iev])}:{int(lumi[iev])}:{int(event[iev])}\n')
            file.close()
            print(f"Saving run:lumi:event to {self.checkOverlap}")

        return accumulator
