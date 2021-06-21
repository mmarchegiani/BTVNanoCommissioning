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
    def __init__(self, year=2017, JECfolder='correction_files', nTrueFile=''):
        self.year = year
        self._mask_fatjets = {
            'basic'       : {
                'pt_cut' : 250.,
                'eta_cut': 2.4,
                'jetId_cut': 3 if self.year==2016 else 2,
                'mass_cut' : 20.,
                'tau21_cut' : 1.1
                    },
            'pt350msd50'       : {
                'pt_cut' : 350.,
                'eta_cut': 2.4,
                'jetId_cut': 3 if self.year==2016 else 2,
                'mass_cut' : 50.,
                'tau21_cut' : 1.1
                    },
            'msd100tau06'       : {
                'pt_cut' : 350.,
                'eta_cut': 2.4,
                'jetId_cut': 3 if self.year==2016 else 2,
                'mass_cut' : 100.,
                'tau21_cut' : 0.6
                    },
        }
        self.year = year
        self.corrJECfolder = JECfolder

        ############
        # PU files
        if self.year == 2016:
            self.puFile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/PrelLum15And1613TeV/PileupHistogram-goldenJSON-13tev-2016-69200ub-99bins.root'
            self.nTrueFile = os.getcwd()+'/correction_files/nTrueInt_datasets_btag2017_2016.coffea'
        if self.year == 2017:
            self.puFile = os.getcwd()+'/correction_files/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root'
            #self.puFile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root'
            self.nTrueFile = os.getcwd()+'/correction_files/nTrueInt_datasets_btag2017_2017.coffea'
        if self.year == 2018:
            #self.puFile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root'
            self.puFile = os.getcwd()+'/correction_files/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root'
            self.nTrueFile = os.getcwd()+'/correction_files/nTrueInt_datasets_btag2017_2018.coffea'
        if nTrueFile: self.nTrueFile = nTrueFile

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

        #self.jet_hists = list(_hist_jet_dict.keys())
        self.fatjet_hists = list(_hist_fatjet_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())

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
        '''Based on https://github.com/andrzejnovak/coffeandbacon/blob/master/analysis/compile_corrections.py#L166-L192'''

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
        '''Based on https://coffeateam.github.io/coffea/notebooks/applying_corrections.html#Applying-energy-scale-transformations-to-Jets'''

        ext = lookup_tools.extractor()
        JECtypes = [ 'L1FastJet', 'L2Relative', 'L2Residual', 'L3Absolute', 'L2L3Residual' ]
        jec_stack_names = [ JECversion+'_'+k+'_'+typeJet for k in JECtypes ]
        JECtypesfiles = [ '* * '+self.corrJECfolder+'/'+k+'.txt' for k in jec_stack_names ]
        ext.add_weight_sets( JECtypesfiles )
        ext.finalize()
        evaluator = ext.make_evaluator()

        print("available evaluator keys:")
        for key in evaluator.keys():
            print("\t", key)

        jec_inputs = {name: evaluator[name] for name in jec_stack_names}
        corrector = FactorizedJetCorrector( **jec_inputs )
        for i in jec_inputs: print(i,'\n',evaluator[i])

        print(dir(evaluator))
        print()
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

        ############
        # Some corrections
        weights = processor.Weights(len(events))
        corrections = {}
        if not isRealData:
            weights.add( 'genWeight', events.genWeight)
            weights.add( 'pileup_weight', self.puReweight( self.puFile, self.nTrueFile, dataset )( events.Pileup.nPU )  )

        events.FatJet = self.applyJEC( events.FatJet, events.fixedGridRhoFastjetAll, events.caches[0], 'AK8PFPuppi', isRealData, JECversion )

        cuts = processor.PackedSelection()

        ############
        # Trigger selection
        if self.year == 2016:
            if 'BTagMu_AK4Jet300_Mu5' not in events.HLT.fields:
                self.triggers = [trigger.replace('AK4', '') for trigger in self.triggers]
        elif self.year == 2018:
            for (i, trigger) in enumerate(self.triggers):
                if trigger.strip("HLT_") not in events.HLT.fields:
                    self.triggers[i] = trigger + "_noalgo"

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in self.triggers]
        if isRealData:
            req_trig = np.zeros(len(events), dtype='bool')
            for t in trig_arrs:
                req_trig = req_trig | t
        else:
            req_trig = np.ones(len(events), dtype='bool')
        cuts.add('trigger', ak.to_numpy(req_trig))

        ############
        # Basic cuts
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[(events.Muon.pt > 5) & (abs(events.Muon.eta < 2.4)) & (events.Muon.tightId != 1) & (events.Muon.pfRelIso04_all > 0.15)]
        events.Muon = ak.pad_none(events.Muon, 2, axis=1)

        ## Jet cuts  (not used)
        events.Jet = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.5)]
        #req_jets = (ak.count(events.Jet.pt, axis=1) >= 2)

        ## FatJet cuts
        events.FatJet = events.FatJet[(events.FatJet.pt > self._mask_fatjets['basic']['pt_cut']) & (abs(events.FatJet.eta) <= self._mask_fatjets['basic']['eta_cut']) & (events.FatJet.jetId > self._mask_fatjets['basic']['jetId_cut'])  & (ak.count(events.FatJet.subjets.pt, axis=2) >= 2) ]  ## subjet sel to crosscheck

        ## Event level variables
        eventVariables = {}
        eventVariables['nfatjet'] = ak.num(events.FatJet)

        ## Leading jet variables
        leadfatjet = ak.firsts(events.FatJet)
        leadfatjet['tau21'] = leadfatjet.tau2/leadfatjet.tau1
        subjet1 = ak.pad_none(leadfatjet.subjets, 2)[:, 0]
        subjet2 = ak.pad_none(leadfatjet.subjets, 2)[:, 1]
        leadfatjet['nsv1'] = get_nsv( subjet1, events.SV )
        leadfatjet['nsv2'] = get_nsv( subjet2, events.SV )
        leadfatjet['nmusj1'] = ak.num(subjet1.delta_r(events.Muon) < 0.4)
        leadfatjet['nmusj2'] = ak.num(subjet2.delta_r(events.Muon) < 0.4)

        fatjet_mutag = (leadfatjet.nmusj1 >= 1) & (leadfatjet.nmusj2 >= 1)
        cuts.add( 'fatjet_mutag', ak.to_numpy(fatjet_mutag) )

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
                    (leadfatjet.tau21 < cut['tau21_cut'])

            cuts.add( selname, ak.to_numpy( sel ) )

        selection = {}
        selection['basic'] = { 'trigger', 'basic' }
        selection['pt350msd50'] = { 'trigger', 'fatjet_mutag', 'pt350msd50' }
        selection['msd100tau06'] = { 'trigger', 'fatjet_mutag', 'msd100tau06' }

        for histname, h in output.items():
            sel = [ r for r in selection.keys() if r in histname.split('_') ]
            if ((histname in self.fatjet_hists) | ('hist2d_fatjet' in histname)):
                for flav, mask in flavors.items():
                    weight = weights.weight() * cuts.all(*selection[sel[0]]) * ak.to_numpy(mask)
                    fields = {k: ak.fill_none(leadfatjet[k], -9999) for k in h.fields if k in dir(leadfatjet)}
                    h.fill(dataset=dataset, flavor=flav, **fields, weight=weight)
            if histname in self.event_hists:
                for flav, mask in flavors.items():
                    weight = weights.weight() * cuts.all(*selection[sel[0]]) * ak.to_numpy(mask)
                    fields = {k: ak.fill_none(eventVariables[k], -9999) for k in h.fields if k in eventVariables.keys() }
                    h.fill(dataset=dataset, flavor=flav, **fields, weight=weight)

        return output

    def postprocess(self, accumulator):
        return accumulator
