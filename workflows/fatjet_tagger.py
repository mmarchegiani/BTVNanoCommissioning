import coffea
from coffea import hist, processor, lookup_tools
from coffea.util import load
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory
import os
import numpy as np
import awkward as ak
import uproot
from utils import rescale, get_nsv, lumi, xsecs, JECversions


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self, year=2017, JECfolder='correction_files'):
        self.year = year
        self._mask_fatjets = {
          'basic'       : None,
          'pt350msd50'  : None,
          'msd100tau06' : None,
        }
        self.year = year
        self.corrJECfolder = JECfolder
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

        # Jet
        #jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 100, 20, 400)
        #jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 60, -3, 3)
        #jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 60, -3, 3)
        #jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 100, 0, 50)
        #ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)

        # FatJet
        #fatjet_tau1_axis  = hist.Bin("tau1",  r"lead. FatJet $\tau_{1}$", 50, 0, 1)
        #fatjet_tau2_axis  = hist.Bin("tau2",  r"lead. FatJet $\tau_{2}$", 50, 0, 1)
        fatjet_tau21_axis = hist.Bin("tau21", r"lead. FatJet $\tau_{21}$", 50, 0, 1)
        fatjet_n2b1_axis  = hist.Bin("n2b1", r"lead. FatJet $N_{2}^{(\beta=1)}$", 50, 0, 0.5)
        fatjet_pt_axis    = hist.Bin("pt",   r"lead. FatJet $p_{T}$ [GeV]", 600, 0, 3000)
        fatjet_eta_axis   = hist.Bin("eta",  r"lead. FatJet $\eta$", 60, -3, 3)
        fatjet_phi_axis   = hist.Bin("phi",  r"lead. FatJet $\phi$", 60, -np.pi, np.pi)
        fatjet_mass_axis  = hist.Bin("mass", r"lead. FatJet $m_{SD}$ [GeV]", 1000, 0, 1000)
        #lfjpt_axis     = hist.Bin("lfjpt", r"Leading fatjet $p_{T}$ [GeV]", 250, 0, 1000)

        # Define similar axes dynamically
        disc_list = ["btagCMVA", "btagCSVV2", 'btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC',]
        disc_list_fj = ['btagDDBvLV2', 'btagDDCvLV2', 'btagDDCvBV2',]
        btag_axes = []
        btag_axes_fj = []
        for d in disc_list:
            btag_axes.append(hist.Bin(d, d, 40, 0, 1))
        for d in disc_list_fj:
            btag_axes_fj.append(hist.Bin(d, d, 40, 0, 1))

        # Define histograms from axes
        #_hist_jet_dict = {
        #        'jet_pt'  : hist.Hist("Events", dataset_axis, jet_pt_axis),
        #        'jet_eta' : hist.Hist("Events", dataset_axis, jet_eta_axis),
        #        'jet_phi' : hist.Hist("Events", dataset_axis, jet_phi_axis),
        #        'jet_mass': hist.Hist("Events", dataset_axis, jet_mass_axis),
        #    }

        _hist_fatjet_dict = {
                #'fatjet_tau1'  : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_tau1_axis),
                #'fatjet_tau2'  : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_tau2_axis),
                'fatjet_tau21' : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_tau21_axis),
                'fatjet_n2b1'  : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_n2b1_axis),
                'fatjet_pt'  : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_pt_axis),
                'fatjet_eta' : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_eta_axis),
                'fatjet_phi' : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_phi_axis),
                'fatjet_mass': hist.Hist("Events", dataset_axis, flavor_axis, fatjet_mass_axis),
            }

        for (i, disc) in enumerate(disc_list_fj):
            _hist_fatjet_dict['fatjet_' + disc] = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i])

        # Define 2D histograms
        _hist2d_dict = {}
        for (i, disc) in enumerate(disc_list_fj):
            _hist2d_dict['hist2d_fatjet_pt_vs_' + disc]    = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], fatjet_pt_axis)
            _hist2d_dict['hist2d_fatjet_mass_vs_' + disc]  = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], fatjet_mass_axis)
            _hist2d_dict['hist2d_fatjet_tau21_vs_' + disc] = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], fatjet_tau21_axis)
            _hist2d_dict['hist2d_fatjet_n2b1_vs_' + disc]  = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], fatjet_n2b1_axis)
            _hist2d_dict['hist2d_nsv1_vs_' + disc]         = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], nsv1_axis)
            _hist2d_dict['hist2d_nsv2_vs_' + disc]         = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], nsv2_axis)
            _hist2d_dict['hist2d_nmusj1_vs_' + disc]       = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], nmusj1_axis)
            _hist2d_dict['hist2d_nmusj2_vs_' + disc]       = hist.Hist("Events", dataset_axis, flavor_axis, btag_axes_fj[i], nmusj2_axis)

        # Generate some histograms dynamically
        #for disc, axis in zip(disc_list, btag_axes):
        #    _hist_jet_dict[disc] = hist.Hist("Events", dataset_axis, axis)
        for disc, axis in zip(disc_list_fj, btag_axes_fj):
            _hist_fatjet_dict[disc] = hist.Hist("Events", dataset_axis, flavor_axis, axis)

        _hist_event_dict = {
                #'njet'   : hist.Hist("Events", dataset_axis, njet_axis),
                #'nbjet'  : hist.Hist("Events", dataset_axis, nbjet_axis),
                #'nel'    : hist.Hist("Events", dataset_axis, nel_axis),
                #'nmu'    : hist.Hist("Events", dataset_axis, nmu_axis),
                'nfatjet': hist.Hist("Events", dataset_axis, flavor_axis, nfatjet_axis),
                'nmusj1' : hist.Hist("Events", dataset_axis, flavor_axis, nmusj1_axis),
                'nmusj2' : hist.Hist("Events", dataset_axis, flavor_axis, nmusj2_axis),
                'nsv1'   : hist.Hist("Events", dataset_axis, flavor_axis, nsv1_axis),
                'nsv2'   : hist.Hist("Events", dataset_axis, flavor_axis, nsv2_axis),
            }
        _sumw_dict = {'sumw': processor.defaultdict_accumulator(float),
                      'nbtagmu': processor.defaultdict_accumulator(float),
                      'nbtagmu_event_level': processor.defaultdict_accumulator(float),
            }

        #self.jet_hists = list(_hist_jet_dict.keys())
        self.fatjet_hists = list(_hist_fatjet_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())

        ############
        # PU files
        if self.year == 2016:
            self.puFile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/PrelLum15And1613TeV/PileupHistogram-goldenJSON-13tev-2016-69200ub-99bins.root'
            self.nTrueFile = os.getcwd()+'/correction_files/nTrueInt_datasets_btag2017_2016.coffea'
        if self.year == 2017:
            self.puFile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root'
            self.nTrueFile = os.getcwd()+'/correction_files/nTrueInt_datasets_btag2017_2017.coffea'
        if self.year == 2018:
            self.puFile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root'
            self.nTrueFile = os.getcwd()+'/correction_files/nTrueInt_datasets_btag2017_2018.coffea'

        #_hist_dict = {**_hist_jet_dict, **_hist_fatjet_dict, **_hist2d_dict, **_hist_event_dict, **_sumw_dict}
        self._hist_dict = {**_hist_fatjet_dict, **_hist2d_dict, **_hist_event_dict}
        self.append_mask()
        self._hist_dict.update({**_sumw_dict})
        self._accumulator = processor.dict_accumulator(self._hist_dict)

    def append_mask(self):
        masks = list(self._mask_fatjets.keys())
        d = {}
        for histname in self._hist_dict.keys():
            h = self._hist_dict[histname]
            d[f'{histname}_{masks[0]}'] = h
            for maskname in masks[1:]:
                d[f'{histname}_{maskname}'] = h.copy()
        self._hist_dict = d.copy()

        l = []
        for histname in self.fatjet_hists:
            for maskname in masks:
                l.append(f'{histname}_{maskname}')
        self.fatjet_hists = l
        l = []
        for histname in self.event_hists:
            for maskname in masks:
                l.append(f'{histname}_{maskname}')
        self.event_hists = l

        return self._hist_dict

    def puReweight(self, puFile, nTrueFile, dataset ):

        nTrueIntLoad = load(nTrueFile)
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

        ext = lookup_tools.extractor()
        JECtypes = [ 'L1FastJet', 'L2Relative', 'L2Residual', 'L3Absolute', 'L2L3Residual' ]
        jec_stack_names = [ JECversion+'_'+k+'_'+typeJet for k in JECtypes ]
        JECtypesfiles = [ '* * '+self.corrJECfolder+'/'+k+'.txt' for k in jec_stack_names ]
        ext.add_weight_sets( JECtypesfiles )
        #        [
        #    "* * "+self.corrJECfolder+"/Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi.txt",
        #    "* * "+self.corrJECfolder+"/Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFPuppi.txt",
        #])
        ext.finalize()
        evaluator = ext.make_evaluator()

        print("available evaluator keys:")
        for key in evaluator.keys():
            print("\t", key)

        #print()
        #print("Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi:")
        #print(evaluator['Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi'])

        #jec_stack_names = ["Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi",
        #           "Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFPuppi"]

        jec_inputs = {name: evaluator[name] for name in jec_stack_names}
        jec_stack = JECStack(jec_inputs)

        print(dir(evaluator))
        print()
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

        #corrector = FactorizedJetCorrector(
        #    Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi=evaluator['Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi'],
        #    Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFPuppi=evaluator['Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFPuppi'],
        #)

        jet_factory = CorrectedJetsFactory(name_map, jec_stack)
        corrected_jets = jet_factory.build(jets, lazy_cache=events_cache)
        print()
        print('starting columns:',ak.fields(jets))
        print()

        print('untransformed pt ratios',jets.pt/jets.pt_raw)
        print('untransformed mass ratios',jets.mass/jets.mass_raw)

        print('transformed pt ratios',corrected_jets.pt/corrected_jets.pt_raw)
        print('transformed mass ratios',corrected_jets.mass/corrected_jets.mass_raw)

        print()
        print('transformed columns:', ak.fields(corrected_jets))
        return corrected_jets


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        dataset = events.metadata['dataset']

        isRealData = 'genWeight' not in events.fields
        if not isRealData:
            output['sumw'][dataset] += sum(events.genWeight)
            JECversion = JECversions[str(self.year)]['MC']
        else:
            output['nbtagmu'][dataset] += ak.count(events.event)
            JECversion = JECversions[str(self.year)]['Data'][dataset.split('BTagMu')[1]]

        weights = processor.Weights(len(events))
        corrections = {}
        if not isRealData:
            weights.add( 'genWeight', events.genWeight)
            weights.add( 'pileup_weight', self.puReweight( self.puFile, self.nTrueFile, dataset )( events.Pileup.nPU )  )

        events.FatJet = self.applyJEC( events.FatJet, events.fixedGridRhoFastjetAll, events.caches[0], 'AK8PFPuppi', isRealData, JECversion )

        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        ##############
        # Trigger level
        triggers = [
        #"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        #"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_BTagMu_AK8Jet300_Mu5",
        "HLT_BTagMu_AK4Jet300_Mu5",
        ]

        if self.year == 2016:
            jetId_cut = 3
            #if 'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ' not in events.HLT:
            #    triggers = [trigger.replace('IsoVL_DZ', 'IsoVL') for trigger in triggers]
            if 'BTagMu_AK4Jet300_Mu5' not in events.HLT.fields:
                triggers = [trigger.replace('AK4', '') for trigger in triggers]
        elif self.year == 2017:
            jetId_cut = 2
        elif self.year == 2018:
            jetId_cut = 2
            for (i, trigger) in enumerate(triggers):
                if trigger.strip("HLT_") not in events.HLT.fields:
                    triggers[i] = trigger + "_noalgo"

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        #req_trig = np.ones(len(events), dtype='bool')
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        baseline_jet    = {var : flatten(events.Jet[var]) for var in ['pt', 'eta', 'phi', 'mass']}
        baseline_fatjet = {var : flatten(events.FatJet[var]) for var in ['pt', 'eta', 'phi', 'msoftdrop']}

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        #events.Muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta < 2.4))] # & (events.Muon.tightId > .5)
        events.Muon = events.Muon[(events.Muon.pt > 5) & (abs(events.Muon.eta < 2.4)) & (events.Muon.tightId != 1) & (events.Muon.pfRelIso04_all > 0.15)]
        events.Muon = ak.pad_none(events.Muon, 2, axis=1)
        #req_muons =(ak.count(events.Muon.pt, axis=1) >= 2)

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        #events.Electron = events.Electron[(events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4)]
        events.Electron = events.Electron[(events.Electron.pt > 10) & (abs(events.Electron.eta) < 2.4)]
        events.Electron = ak.pad_none(events.Electron, 1, axis=1)
        #req_ele = (ak.count(events.Electron.pt, axis=1) == 1)

        ## Jet cuts
        events.Jet = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.5)]
        #req_jets = (ak.count(events.Jet.pt, axis=1) >= 2)

        for selname in self._mask_fatjets.keys():
            ## FatJet cuts
            if selname == 'basic':
                pt_cut    = 250
                mass_cut  = 20
                tau21_cut = 1.1
            elif selname == 'pt350msd50':
                pt_cut    = 350
                mass_cut  = 50
                tau21_cut = 1.1
            elif selname == 'msd100tau06':
                pt_cut    = 350
                mass_cut  = 100
                tau21_cut = 0.6
            sfatjets = events.FatJet[(events.FatJet.pt > pt_cut) &
                                     (events.FatJet.msoftdrop > mass_cut) &
                                     (abs(events.FatJet.eta) < 2.4) &
                                     (events.FatJet.jetId >= jetId_cut) &
                                     ((events.FatJet.tau2/events.FatJet.tau1) < tau21_cut)]
            sfatjets = ak.firsts(sfatjets)
            tmp = ak.is_none(sfatjets)
            #sfatjets = sfatjets[~tmp]
            #sfatjets['tau21'] = sfatjets.tau2/sfatjets.tau1
            req_fatjets = ~tmp  #(ak.count(sfatjets.pt, axis=0) >= 1)
            req_subjets = ak.any(ak.count(sfatjets.subjets.pt, axis=1) >= 2)


            #req_opposite_charge = events.Electron[:, 0].charge * events.Muon[:, 0].charge == -1

            event_level = req_trig & req_fatjets & req_subjets #& req_muons
            #event_level = req_trig & req_fatjets

            # Selected
            selev = events[event_level]
            nEvents = ak.count(selev.event)

            #########

            # Per electron
            el_eta   = (abs(selev.Electron.eta) <= 2.4)
            el_pt    = selev.Electron.pt > 10
            el_level = el_eta & el_pt

            # Per muon
            mu_eta   = (abs(selev.Muon.eta) <= 2.4)
            mu_pt    = selev.Muon.pt > 10
            mu_not_iso = (selev.Muon.pfRelIso04_all > 0.15)
            mu_not_tight = (selev.Muon.tightId != 1)
            mu_level = mu_eta & mu_pt & mu_not_iso & mu_not_tight

            # Per jet
            jet_eta    = (abs(selev.Jet.eta) <= 2.4)
            jet_pt     = selev.Jet.pt > 25
            jet_pu     = selev.Jet.puId > 6
            jet_level  = jet_pu & jet_eta & jet_pt

            # Per fatjet
            fatjet_pt      = selev.FatJet.pt > pt_cut
            fatjet_mass    = selev.FatJet.msoftdrop > mass_cut
            fatjet_tau21   = (selev.FatJet.tau2/selev.FatJet.tau1) < tau21_cut
            fatjet_subjets = (ak.count(selev.FatJet.subjets.pt, axis=2) >= 2)

            #fatjet_level = fatjet_pt & fatjet_mass & fatjet_tau21
            fatjet_level = fatjet_pt & fatjet_mass & fatjet_tau21 & fatjet_subjets
            self._mask_fatjets[selname] = fatjet_level
            selection = self._mask_fatjets[selname]
            # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
            bjet_disc  = selev.Jet.btagDeepB > 0.7264 # L=0.0494, M=0.2770, T=0.7264
            bjet_level = jet_level & bjet_disc

            sel      = selev.Electron[el_level]
            smu      = selev.Muon[mu_level]
            sjets    = selev.Jet[jet_level]
            sbjets   = selev.Jet[bjet_level]
            nfatjet  = ak.num(selev.FatJet[selection])

            #sfatjets = ak.pad_none(selev.FatJet[selection], 1)[:,0]
            sfatjets = ak.firsts(selev.FatJet[selection])
            tmp = ak.is_none(sfatjets)
            sfatjets = sfatjets[~tmp]
            sfatjets['tau21'] = sfatjets.tau2/sfatjets.tau1
            subjet1  = ak.firsts(sfatjets.subjets) #ak.pad_none(sfatjets.subjets, 2)[:, 0]
            subjet2  = ak.firsts(sfatjets.subjets, 1) #ak.pad_none(sfatjets.subjets, 2)[:, 1]
            SV       = selev.SV[~tmp]
            nsv1     = get_nsv(subjet1, SV)
            nsv2     = get_nsv(subjet2, SV)
            nmusj1   = ak.num(subjet1.delta_r(smu[~tmp]) < 0.4)
            nmusj2   = ak.num(subjet2.delta_r(smu[~tmp]) < 0.4)

            fatjet_mutag = (nmusj1 >= 1) & (nmusj2 >= 1)
            selection = fatjet_mutag

            #fatjet_tau21 = sfatjets.tau21 < 0.75
            #fatjet_nsv1 = nsv1 > 0
            #fatjet_nsv12 = (nsv1 > 0) & (nsv2 > 0)
            #fatjet_mutag = fatjet_mutag & fatjet_nsv12
            sfatjets = sfatjets[fatjet_mutag]
            sfatjets['tau21'] = sfatjets.tau2/sfatjets.tau1
            subjet1  = subjet1[fatjet_mutag]
            subjet2  = subjet2[fatjet_mutag]
            nsv1     = nsv1[fatjet_mutag]
            nsv2     = nsv2[fatjet_mutag]
            nmusj1   = nmusj1[fatjet_mutag]
            nmusj2   = nmusj2[fatjet_mutag]

            fatjet_mutag_level = fatjet_mutag  #(ak.any(selection[~tmp], axis=1) & fatjet_mutag)

            # Flavor matching
            if not isRealData:
                _b = (sfatjets.hadronFlavour == 5)
                _c = (sfatjets.hadronFlavour == 4)
                _l = (sfatjets.hadronFlavour < 4)
                _bb = abs(sfatjets.hadronFlavour == 5) & (sfatjets.nBHadrons >= 2) #& (sfatjets.nCHadrons == 0)
                _cc = abs(sfatjets.hadronFlavour == 4) & (sfatjets.nBHadrons == 0) & (sfatjets.nCHadrons >= 2)
                #_ll = abs(sfatjets.hadronFlavour < 4) & (sfatjets.nBHadrons == 0) & (sfatjets.nCHadrons == 0)
                _b = _b & ~_bb
                _c = _c & ~_cc
                _l = _l & ~_bb & ~_cc & ~_b & ~_c
                #_others = ~_l & ~_bb & ~_cc & ~_b & ~_c
                flavor = _bb*5 + _cc*4 + _b*3 + _c*2 + _l*1
            else:
                output['nbtagmu_event_level'][dataset] += ak.count_nonzero(event_level)

            # Fill histograms dynamically
            for histname, h in output.items():
                if not selname in histname: continue
                #if histname in self.jet_hists:
                #    fields = {k: flatten(sjets[k]) for k in h.fields if k in dir(sjets)}
                #    fields.update({k: flatten(sjets[k]) for k in h.fields if k.split('jet_')[-1] in ['pt', 'eta', 'phi', 'mass']})
                #    if isRealData:
                #        h.fill(dataset=dataset, **fields)
                #    else:
                #        h.fill(dataset=dataset, **fields, weight=weights.weight()[event_level][fatjet_mutag_level])
                if ((histname in self.fatjet_hists) | ('hist2d_fatjet' in histname)):
                    fields = {k: flatten(sfatjets[k]) for k in h.fields if k in dir(sfatjets)}
                    #fields.update({k: flatten(sfatjets[k]) for k in h.fields if k.split('fatjet_')[-1] in ['pt', 'eta', 'phi', 'msoftdrop']})
                    #h.fill(dataset=dataset, flavor="inclusive", **fields, weight=weights.weight()[event_level][fatjet_mutag_level][mask])
                    if isRealData:
                        h.fill(dataset=dataset, flavor="Data", **fields)
                    else:
                        #for flav, mask in zip(['light', 'c', 'b', 'cc', 'bb', 'others'], [_l, _c, _b, _cc, _bb, _others]):
                        for flav, mask in zip(['light', 'c', 'b', 'cc', 'bb'], [_l, _c, _b, _cc, _bb]):
                            mask= ak.is_none(mask, False)
                            sfatjets_flavor = sfatjets[mask]
                            fields = {k: flatten(sfatjets_flavor[k]) for k in h.fields if k in dir(sfatjets_flavor)}
                            #fields.update({k: flatten(sfatjets_flavor[k]) for k in h.fields if k.split('fatjet_')[-1] in ['pt', 'eta', 'phi', 'msoftdrop']})
                            if len(weights.weight()[event_level][~tmp][fatjet_mutag_level][mask])>0: print('jdhfsdhsklh')
                            h.fill(dataset=dataset, flavor=flav, **fields, weight=weights.weight()[event_level][~tmp][fatjet_mutag_level][mask])

                elif (((histname in self.event_hists) | ('hist2d_nsv' in histname) | ('hist2d_nmusj' in histname)) & (not histname in ['njet', 'nbjet', 'nel', 'nmu'])):
                    fields = {k: flatten(sfatjets[k]) for k in h.fields if k in dir(sfatjets)}
                    for varname, values in zip(['nfatjet', 'nsv1', 'nsv2', 'nmusj1', 'nmusj2'], [nfatjet, nsv1, nsv2, nmusj1, nmusj2]):
                        if varname in histname:
                            fields.update({varname: flatten(values)})
                    #h.fill(dataset=dataset, flavor="inclusive", **fields, weight=weights.weight()[event_level][fatjet_mutag_level][mask])
                    if isRealData:
                        h.fill(dataset=dataset, flavor="Data", **fields)
                    else:
                        for flav, mask in zip(['light', 'c', 'b', 'cc', 'bb'], [_l, _c, _b, _cc, _bb]):
                            sfatjets_flavor = sfatjets[mask]
                            fields = {k: flatten(sfatjets_flavor[k]) for k in h.fields if k in dir(sfatjets_flavor)}
                            for varname, values in zip(['nfatjet', 'nsv1', 'nsv2', 'nmusj1', 'nmusj2'], [nfatjet, nsv1, nsv2, nmusj1, nmusj2]):
                                if varname in histname:
                                    fields.update({varname: flatten(values[mask])})
                            h.fill(dataset=dataset, flavor=flav, **fields) #, weight=weights.weight()[event_level][fatjet_mutag_level][mask])
                else: continue

        #output['njet'].fill(dataset=dataset,    njet=ak.num(sjets))
        #output['nbjet'].fill(dataset=dataset,   nbjet=ak.num(sbjets))
        #output['nel'].fill(dataset=dataset,     nel=ak.num(sel))
        #output['nmu'].fill(dataset=dataset,     nmu=ak.num(smu))
        #output['nfatjet'].fill(dataset=dataset, flavor="inclusive", nfatjet=num(sfatjets))
        #output['nmusj1'].fill(dataset=dataset,  flavor="inclusive", nmusj1=flatten(nmusj1))
        #output['nmusj2'].fill(dataset=dataset,  flavor="inclusive", nmusj2=flatten(nmusj2))
        #output['nsv1'].fill(dataset=dataset,    flavor="inclusive", nsv1=ak.fill_none(nsv1, -1))
        #output['nsv2'].fill(dataset=dataset,    flavor="inclusive", nsv2=ak.fill_none(nsv2, -1))

        return output

    def postprocess(self, accumulator):

        #isSplit = (len(accumulator['sumw'].keys()) <= 1)
        #if not isSplit:
            #accumulator = rescale(accumulator, xsecs, lumi[self.year])
            #accumulator = rescale(accumulator, xsecs, lumi[2017])

        return accumulator
