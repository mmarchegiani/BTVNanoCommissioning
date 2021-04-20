import coffea
from coffea import hist, processor
import numpy as np
#import awkward1 as ak
import awkward as ak
from utils import rescale, get_nsv, lumi, xsecs


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self):
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
        fatjet_tau1_axis  = hist.Bin("tau1",  r"lead. FatJet $\tau_{1}$", 50, 0, 1)
        fatjet_tau2_axis  = hist.Bin("tau2",  r"lead. FatJet $\tau_{2}$", 50, 0, 1)
        fatjet_tau21_axis = hist.Bin("tau21", r"lead. FatJet $\tau_{21}$", 50, 0, 1)
        fatjet_n2b1_axis  = hist.Bin("n2b1", r"lead. FatJet $N_{2}^{(\beta=1)}$", 50, 0, 0.5)
        fatjet_pt_axis    = hist.Bin("pt",   r"lead. FatJet $p_{T}$ [GeV]", 250, 0, 1000)
        fatjet_eta_axis   = hist.Bin("eta",  r"lead. FatJet $\eta$", 60, -3, 3)
        fatjet_phi_axis   = hist.Bin("phi",  r"lead. FatJet $\phi$", 60, -np.pi, np.pi)
        fatjet_mass_axis  = hist.Bin("mass", r"lead. FatJet $m_{SD}$ [GeV]", 300, 0, 300)
        #lfjpt_axis     = hist.Bin("lfjpt", r"Leading fatjet $p_{T}$ [GeV]", 250, 0, 1000)

        # Define similar axes dynamically
        disc_list = ["btagCMVA", "btagCSVV2", 'btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC',]
        disc_list_fj = ['btagDDBvLV2', 'btagDDCvLV2', 'btagDDCvBV2',]
        btag_axes = []
        btag_axes_fj = []
        for d in disc_list:
            btag_axes.append(hist.Bin(d, d, 50, 0, 1))
        for d in disc_list_fj:
            btag_axes_fj.append(hist.Bin(d, d, 50, 0, 1))

        # Define histograms from axes
        #_hist_jet_dict = {
        #        'jet_pt'  : hist.Hist("Events", dataset_axis, jet_pt_axis),
        #        'jet_eta' : hist.Hist("Events", dataset_axis, jet_eta_axis),
        #        'jet_phi' : hist.Hist("Events", dataset_axis, jet_phi_axis),
        #        'jet_mass': hist.Hist("Events", dataset_axis, jet_mass_axis),
        #    }
        
        _hist_fatjet_dict = {
                'fatjet_tau1'  : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_tau1_axis),
                'fatjet_tau2'  : hist.Hist("Events", dataset_axis, flavor_axis, fatjet_tau2_axis),
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

        #_hist_dict = {**_hist_jet_dict, **_hist_fatjet_dict, **_hist2d_dict, **_hist_event_dict, **_sumw_dict}
        _hist_dict = {**_hist_fatjet_dict, **_hist2d_dict, **_hist_event_dict, **_sumw_dict}
        self._accumulator = processor.dict_accumulator(_hist_dict)


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        
        dataset = events.metadata['dataset']

        isRealData = 'genWeight' not in events.fields
        if not isRealData:
            output['sumw'][dataset] += sum(events.genWeight)
        else:
            output['nbtagmu'][dataset] += ak.count(events.event)

        ##############
        # Trigger level
        triggers = [
        #"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_BTagMu_AK8Jet300_Mu5",
        "HLT_BTagMu_AK4Jet300_Mu5",
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        #req_trig = np.ones(len(events), dtype='bool')
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        baseline_jet    = {var : ak.flatten(events.Jet[var], axis=None) for var in ['pt', 'eta', 'phi', 'mass']}
        baseline_fatjet = {var : ak.flatten(events.FatJet[var], axis=None) for var in ['pt', 'eta', 'phi', 'msoftdrop']}

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

        ## FatJet cuts
        events.FatJet = events.FatJet[(events.FatJet.pt > 250) & (events.FatJet.mass > 20)]
        #events.FatJet['tau21'] = events.FatJet.tau2/events.FatJet.tau1
        req_fatjets = (ak.count(events.FatJet.pt, axis=1) >= 1)
        req_subjets = ak.any(ak.count(events.FatJet.subjets.pt, axis=2) >= 2, axis=1)

        #req_opposite_charge = events.Electron[:, 0].charge * events.Muon[:, 0].charge == -1

        event_level = req_trig & req_fatjets & req_subjets #& req_muons 
        #event_level = req_trig & req_fatjets 

        # Selected
        selev = events[event_level]
        sweight = None
        if not isRealData:
            sweight = selev.genWeight

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
        fatjet_pt    = selev.FatJet.pt > 250
        fatjet_mass  = selev.FatJet.mass > 20
        #fatjet_tau21 = (selev.FatJet.tau2/selev.FatJet.tau1) < 0.5
        fatjet_subjets = (ak.count(selev.FatJet.subjets.pt, axis=2) >= 2)

        #fatjet_level = fatjet_pt & fatjet_mass & fatjet_tau21
        fatjet_level = fatjet_pt & fatjet_mass & fatjet_subjets
        #fatjet_level = ak.ones_like(fatjet_pt)

        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc  = selev.Jet.btagDeepB > 0.7264 # L=0.0494, M=0.2770, T=0.7264
        bjet_level = jet_level & bjet_disc

        sel      = selev.Electron[el_level]
        smu      = selev.Muon[mu_level]
        sjets    = selev.Jet[jet_level]
        sbjets   = selev.Jet[bjet_level]
        sfatjets = ak.pad_none(selev.FatJet[fatjet_level], 1)[:,0]
        sfatjets['tau21'] = sfatjets.tau2/sfatjets.tau1
        subjet1  = ak.pad_none(sfatjets.subjets, 2)[:, 0]
        subjet2  = ak.pad_none(sfatjets.subjets, 2)[:, 1]
        SV       = selev.SV
        nsv1     = get_nsv(subjet1, SV)
        nsv2     = get_nsv(subjet2, SV)
        nmusj1   = ak.num(subjet1.delta_r(smu) < 0.4)
        nmusj2   = ak.num(subjet2.delta_r(smu) < 0.4)

        fatjet_mutag = (nmusj1 >= 1) & (nmusj2 >= 1)
        #fatjet_tau21 = sfatjets.tau21 < 0.75
        #fatjet_nsv1 = nsv1 > 0
        #fatjet_nsv12 = (nsv1 > 0) & (nsv2 > 0)
        #fatjet_mutag = fatjet_mutag & fatjet_nsv12
        sfatjets = sfatjets[fatjet_mutag]
        sfatjets['tau21'] = sfatjets.tau2/sfatjets.tau1
        subjet1  = subjet1[fatjet_mutag]
        subjet2  = subjet2[fatjet_mutag]
        #print('nmusj1:', len(nmusj1))
        #print('fatjet_mutag:', len(fatjet_mutag))
        nmusj1   = nmusj1[fatjet_mutag]
        nmusj2   = nmusj2[fatjet_mutag]
        nsv1     = nsv1[fatjet_mutag]
        nsv2     = nsv2[fatjet_mutag]

        sweight_jets = None
        sweight_fatjets = None
        if not isRealData:
            #sweight_fatjets = sweight[ak.any(fatjet_level, axis=1) & fatjet_mutag]
            sweight_fatjets = sweight[ak.any(fatjet_level, axis=1)]
            sweight_fatjets = ak.flatten(sweight_fatjets, axis=None)

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

        """
        if not isRealData:
            sbbfatjets = ak.pad_none(selev.FatJet[_bb], 1)[:,0]
            sbbfatjets['tau21'] = sbbfatjets.tau2/sbbfatjets.tau1
            sccfatjets = ak.pad_none(selev.FatJet[_cc], 1)[:,0]
            sccfatjets['tau21'] = sccfatjets.tau2/sccfatjets.tau1
        """

        # output['pt'].fill(dataset=dataset, pt=selev.Jet.pt.flatten())
        # Fill histograms dynamically
        for histname, h in output.items():
            #if histname in self.jet_hists:
            #    fields = {k: ak.flatten(sjets[k], axis=None) for k in h.fields if k in dir(sjets)}
            #    fields.update({k: ak.flatten(sjets[k], axis=None) for k in h.fields if k.split('jet_')[-1] in ['pt', 'eta', 'phi', 'mass']})
            #    if isRealData:
            #        h.fill(dataset=dataset, **fields)
            #    else:
            #        h.fill(dataset=dataset, **fields, weight=sweight_jets)
            if ((histname in self.fatjet_hists) | ('hist2d_fatjet' in histname)):
                fields = {k: ak.flatten(sfatjets[k], axis=None) for k in h.fields if k in dir(sfatjets)}
                #fields.update({k: ak.flatten(sfatjets[k], axis=None) for k in h.fields if k.split('fatjet_')[-1] in ['pt', 'eta', 'phi', 'msoftdrop']})
                #h.fill(dataset=dataset, flavor="inclusive", **fields, weight=sweight_fatjets)
                if isRealData:
                    h.fill(dataset=dataset, flavor="Data", **fields)
                else:
                    #for flav, mask in zip(['light', 'c', 'b', 'cc', 'bb', 'others'], [_l, _c, _b, _cc, _bb, _others]):
                    for flav, mask in zip(['light', 'c', 'b', 'cc', 'bb'], [_l, _c, _b, _cc, _bb]):
                        sfatjets_flavor = sfatjets[mask]
                        sweight_fatjets_flavor = ak.flatten(sweight_fatjets[mask], axis=None)
                        fields = {k: ak.flatten(sfatjets_flavor[k], axis=None) for k in h.fields if k in dir(sfatjets_flavor)}
                        #fields.update({k: ak.flatten(sfatjets_flavor[k], axis=None) for k in h.fields if k.split('fatjet_')[-1] in ['pt', 'eta', 'phi', 'msoftdrop']})
                        h.fill(dataset=dataset, flavor=flav, **fields, weight=sweight_fatjets_flavor)

            elif (((histname in self.event_hists) | ('hist2d_nsv' in histname) | ('hist2d_nmusj' in histname)) & (not histname in ['njet', 'nbjet', 'nel', 'nmu', 'nfatjet'])):
                fields = {k: ak.flatten(sfatjets[k], axis=None) for k in h.fields if k in dir(sfatjets)}
                for varname, values in zip(['nsv1', 'nsv2', 'nmusj1', 'nmusj2'], [nsv1, nsv2, nmusj1, nmusj2]):
                    if varname in histname:
                        fields.update({varname: ak.flatten(values, axis=None)})
                #h.fill(dataset=dataset, flavor="inclusive", **fields, weight=sweight_fatjets)
                if isRealData:
                    h.fill(dataset=dataset, flavor="Data", **fields)
                else:
                    for flav, mask in zip(['light', 'c', 'b', 'cc', 'bb'], [_l, _c, _b, _cc, _bb]):
                        sfatjets_flavor = sfatjets[mask]
                        sweight_fatjets_flavor = ak.flatten(sweight_fatjets[mask], axis=None)                        
                        fields = {k: ak.flatten(sfatjets_flavor[k], axis=None) for k in h.fields if k in dir(sfatjets_flavor)}
                        for varname, values in zip(['nsv1', 'nsv2', 'nmusj1', 'nmusj2'], [nsv1, nsv2, nmusj1, nmusj2]):
                            if varname in histname:
                                fields.update({varname: ak.flatten(values[mask], axis=None)})
                        h.fill(dataset=dataset, flavor=flav, **fields, weight=sweight_fatjets_flavor)
            else: continue

        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)

        #output['njet'].fill(dataset=dataset,    njet=ak.num(sjets))
        #output['nbjet'].fill(dataset=dataset,   nbjet=ak.num(sbjets))
        #output['nel'].fill(dataset=dataset,     nel=ak.num(sel))
        #output['nmu'].fill(dataset=dataset,     nmu=ak.num(smu))
        output['nfatjet'].fill(dataset=dataset, flavor="inclusive", nfatjet=num(sfatjets))
        #output['nmusj1'].fill(dataset=dataset,  flavor="inclusive", nmusj1=flatten(nmusj1))
        #output['nmusj2'].fill(dataset=dataset,  flavor="inclusive", nmusj2=flatten(nmusj2))
        #output['nsv1'].fill(dataset=dataset,    flavor="inclusive", nsv1=ak.fill_none(nsv1, -1))
        #output['nsv2'].fill(dataset=dataset,    flavor="inclusive", nsv2=ak.fill_none(nsv2, -1))

        return output

    def postprocess(self, accumulator):

        #isSplit = (len(accumulator['sumw'].keys()) <= 1)
        #if not isSplit:
            #accumulator = rescale(accumulator, xsecs, lumi[args.year])
            #accumulator = rescale(accumulator, xsecs, lumi[2017])

        return accumulator
