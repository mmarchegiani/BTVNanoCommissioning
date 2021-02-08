import coffea
from coffea import hist, processor
import numpy as np
import awkward1 as ak


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self):
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        cutflow_axis   = hist.Cat("cut",   "Cut")

        # Events
        nel_axis   = hist.Bin("nel",   r"N electrons", [0,1,2,3,4,5,6,7,8,9,10])
        nmu_axis   = hist.Bin("nmu",   r"N muons",     [0,1,2,3,4,5,6,7,8,9,10])
        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3,4,5,6,7,8,9,10])
        nbjet_axis = hist.Bin("nbjet", r"N b-jets",    [0,1,2,3,4,5,6,7,8,9,10])

        # Electron
        el_pt_axis   = hist.Bin("pt",    r"Electron $p_{T}$ [GeV]", 100, 20, 400)
        el_eta_axis  = hist.Bin("eta",   r"$\eta$", 60, -3, 3)
        el_phi_axis  = hist.Bin("phi",   r"$\phi$", 60, -3, 3)
        lelpt_axis   = hist.Bin("lelpt", r"Leading electron $p_{T}$ [GeV]", 100, 20, 200)

        # Muons
        mu_pt_axis   = hist.Bin("pt",    r"Muon $p_{T}$ [GeV]", 100, 20, 400)
        mu_eta_axis  = hist.Bin("eta",   r"$\eta$", 60, -3, 3)
        mu_phi_axis  = hist.Bin("phi",   r"$\phi$", 60, -3, 3)
        lmupt_axis   = hist.Bin("lmupt", r"Leading muon $p_{T}$ [GeV]", 100, 20, 200)

        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 100, 20, 400)
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 60, -3, 3)
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 60, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 100, 0, 50)
        ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)

        # FatJet
        fatjet_tau1_axis  = hist.Bin("tau1",   r"FatJet $\tau_{1}$", 50, 0, 1)
        fatjet_tau2_axis  = hist.Bin("tau2",   r"FatJet $\tau_{2}$", 50, 0, 1)
        fatjet_tau21_axis = hist.Bin("tau21", r"FatJet $\tau_{21}$", 50, 0, 1)
        fatjet_pt_axis    = hist.Bin("pt",   r"FatJet $p_{T}$ [GeV]", 250, 0, 1000)
        fatjet_eta_axis   = hist.Bin("eta",  r"FatJet $\eta$", 60, -3, 3)
        fatjet_phi_axis   = hist.Bin("phi",  r"FatJet $\phi$", 60, -3, 3)
        fatjet_mass_axis  = hist.Bin("mass", r"FatJet $m_{SD}$ [GeV]", 300, 0, 300)
        #lfjpt_axis     = hist.Bin("lfjpt", r"Leading fatjet $p_{T}$ [GeV]", 250, 0, 1000)

        # cc FatJet
        ccfatjet_tau1_axis  = hist.Bin("tau1",   r"cc FatJet $\tau_{1}$", 50, 0, 1)
        ccfatjet_tau2_axis  = hist.Bin("tau2",   r"cc FatJet $\tau_{2}$", 50, 0, 1)
        ccfatjet_tau21_axis = hist.Bin("tau21", r"cc FatJet $\tau_{21}$", 50, 0, 1)
        ccfatjet_pt_axis    = hist.Bin("pt",   r"cc FatJet $p_{T}$ [GeV]", 250, 0, 1000)
        ccfatjet_eta_axis   = hist.Bin("eta",  r"cc FatJet $\eta$", 60, -3, 3)
        ccfatjet_phi_axis   = hist.Bin("phi",  r"cc FatJet $\phi$", 60, -3, 3)
        ccfatjet_mass_axis  = hist.Bin("mass", r"cc FatJet $m_{SD}$ [GeV]", 300, 0, 300)

        # Define similar axes dynamically
        disc_list = ["btagCMVA", "btagCSVV2", 'btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC',]
        disc_list_fj = ['btagDDBvLV2', 'btagDDCvLV2', 'btagDDCvBV2',]
        btag_axes = []
        btag_axes_fj = []
        btag_axes_ccfj = []
        for d in disc_list:
            btag_axes.append(hist.Bin(d, d, 50, 0, 1))
        for d in disc_list_fj:
            btag_axes_fj.append(hist.Bin(d, d, 50, 0, 1))
        for d in disc_list_fj:
            btag_axes_ccfj.append(hist.Bin(d, d, 50, 0, 1))

        # Define histograms from axes
        _hist_jet_dict = {
                'jet_pt'  : hist.Hist("Counts", dataset_axis, jet_pt_axis),
                'jet_eta' : hist.Hist("Counts", dataset_axis, jet_eta_axis),
                'jet_phi' : hist.Hist("Counts", dataset_axis, jet_phi_axis),
                'jet_mass': hist.Hist("Counts", dataset_axis, jet_mass_axis),
            }
        
        _hist_fatjet_dict = {
                'fatjet_tau1'  : hist.Hist("Counts", dataset_axis, fatjet_tau1_axis),
                'fatjet_tau2'  : hist.Hist("Counts", dataset_axis, fatjet_tau2_axis),
                'fatjet_tau21' : hist.Hist("Counts", dataset_axis, fatjet_tau21_axis),
                'fatjet_pt'  : hist.Hist("Counts", dataset_axis, fatjet_pt_axis),
                'fatjet_eta' : hist.Hist("Counts", dataset_axis, fatjet_eta_axis),
                'fatjet_phi' : hist.Hist("Counts", dataset_axis, fatjet_phi_axis),
                'fatjet_mass': hist.Hist("Counts", dataset_axis, fatjet_mass_axis),
            }

        _hist_ccfatjet_dict = {
                'ccfatjet_tau1'  : hist.Hist("Counts", dataset_axis, ccfatjet_tau1_axis),
                'ccfatjet_tau2'  : hist.Hist("Counts", dataset_axis, ccfatjet_tau2_axis),
                'ccfatjet_tau21' : hist.Hist("Counts", dataset_axis, ccfatjet_tau21_axis),
                'ccfatjet_pt'  : hist.Hist("Counts", dataset_axis, ccfatjet_pt_axis),
                'ccfatjet_eta' : hist.Hist("Counts", dataset_axis, ccfatjet_eta_axis),
                'ccfatjet_phi' : hist.Hist("Counts", dataset_axis, ccfatjet_phi_axis),
                'ccfatjet_mass': hist.Hist("Counts", dataset_axis, ccfatjet_mass_axis),
            }
        

        # Generate some histograms dynamically
        for disc, axis in zip(disc_list, btag_axes):
            _hist_jet_dict[disc] = hist.Hist("Counts", dataset_axis, axis)
        for disc, axis in zip(disc_list_fj, btag_axes_fj):
            _hist_fatjet_dict[disc] = hist.Hist("Counts", dataset_axis, axis)
        for disc, axis in zip(disc_list_fj, btag_axes_ccfj):
            _hist_ccfatjet_dict[disc] = hist.Hist("Counts", dataset_axis, axis)

        _hist_event_dict = {
                'njet'  : hist.Hist("Counts", dataset_axis, njet_axis),
                'nbjet' : hist.Hist("Counts", dataset_axis, nbjet_axis),
                'nel'   : hist.Hist("Counts", dataset_axis, nel_axis),
                'nmu'   : hist.Hist("Counts", dataset_axis, nmu_axis),
                #'lelpt' : hist.Hist("Counts", dataset_axis, lelpt_axis),
                #'lmupt' : hist.Hist("Counts", dataset_axis, lmupt_axis),
                #'ljpt'  : hist.Hist("Counts", dataset_axis, ljpt_axis),
                #'lfjpt' : hist.Hist("Counts", dataset_axis, lfjpt_axis),
            }

        self.jet_hists = list(_hist_jet_dict.keys())
        self.fatjet_hists = list(_hist_fatjet_dict.keys())
        self.ccfatjet_hists = list(_hist_ccfatjet_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())

        _hist_dict = {**_hist_jet_dict, **_hist_fatjet_dict, **_hist_ccfatjet_dict, **_hist_event_dict}
        self._accumulator = processor.dict_accumulator(_hist_dict)


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        dataset = events.metadata['dataset']

        ##############
        # Trigger level
        triggers = [
        #"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        baseline_jet    = {var : ak.flatten(events.Jet[var], axis=None) for var in ['pt', 'eta', 'phi', 'mass']}
        baseline_fatjet = {var : ak.flatten(events.FatJet[var], axis=None) for var in ['pt', 'eta', 'phi', 'msoftdrop']}

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta < 2.4))] # & (events.Muon.tightId > .5)
        events.Muon = ak.pad_none(events.Muon, 1, axis=1)
        #req_muon =(ak.count(events.Muon.pt, axis=1) == 1)

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[(events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4)]
        events.Electron = ak.pad_none(events.Electron, 1, axis=1)
        #req_ele = (ak.count(events.Electron.pt, axis=1) == 1)

        ## Jet cuts
        events.Jet = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.5)]
        #req_jets = (ak.count(events.Jet.pt, axis=1) >= 2)

        ## FatJet cuts
        events.FatJet = events.FatJet[(events.FatJet.pt > 250) & (events.FatJet.mass > 20)]
        #events.FatJet['tau21'] = events.FatJet.tau2/events.FatJet.tau1
        req_fatjets = (ak.count(events.FatJet.pt, axis=1) >= 1)

        #req_opposite_charge = events.Electron[:, 0].charge * events.Muon[:, 0].charge == -1

        #event_level = req_trig & req_muon & req_ele & req_opposite_charge & req_jets & req_fatjets
        event_level = req_trig & req_fatjets

        # Selected
        selev = events[event_level]

        #########

        # Per electron
        el_eta   = (abs(selev.Electron.eta) <= 2.4)
        el_pt    = selev.Electron.pt > 30
        el_level = el_eta & el_pt

        # Per muon
        mu_eta   = (abs(selev.Muon.eta) <= 2.4)
        mu_pt    = selev.Muon.pt > 30
        mu_level = mu_eta & mu_pt

        # Per jet
        jet_eta    = (abs(selev.Jet.eta) <= 2.4)
        jet_pt     = selev.Jet.pt > 25
        jet_pu     = selev.Jet.puId > 6
        jet_level  = jet_pu & jet_eta & jet_pt

        # Per fatjet
        fatjet_pt    = selev.FatJet.pt > 250
        fatjet_mass  = selev.FatJet.mass > 20
        #fatjet_tau21 = (selev.FatJet.tau2/selev.FatJet.tau1) < 0.5
        #fatjet_level = fatjet_pt & fatjet_mass & fatjet_tau21
        fatjet_level = fatjet_pt & fatjet_mass

        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc  = selev.Jet.btagDeepB > 0.7264 # L=0.0494, M=0.2770, T=0.7264
        bjet_level = jet_level & bjet_disc

        # Selecting cc FatJet
        hadronFlavour = (abs(selev.FatJet.hadronFlavour == 4))
        nBHadrons = (selev.FatJet.nBHadrons == 0)
        nCHadrons = (selev.FatJet.nCHadrons >= 2)

        ccfatjet_level = fatjet_level & hadronFlavour & nBHadrons & nCHadrons

        sel    = selev.Electron[el_level]
        smu    = selev.Muon[mu_level]
        sjets  = selev.Jet[jet_level]
        sbjets = selev.Jet[bjet_level]
        sfatjets = selev.FatJet[fatjet_level]
        sccfatjets = selev.FatJet[ccfatjet_level]
        sfatjets['tau21'] = sfatjets.tau2/sfatjets.tau1
        sccfatjets['tau21'] = sccfatjets.tau2/sccfatjets.tau1

        # output['pt'].fill(dataset=dataset, pt=selev.Jet.pt.flatten())
        # Fill histograms dynamically
        for histname, h in output.items():
            if histname in self.jet_hists:
                fields = {k: ak.flatten(sjets[k], axis=None) for k in h.fields if k in dir(sjets)}
                fields.update({k: ak.flatten(sjets[k], axis=None) for k in h.fields if k.split('jet_')[-1] in ['pt', 'eta', 'phi', 'mass']})
                h.fill(dataset=dataset, **fields)
            elif histname in self.fatjet_hists:
                fields = {k: ak.flatten(sfatjets[k], axis=None) for k in h.fields if k in dir(sfatjets)}
                fields.update({k: ak.flatten(sfatjets[k], axis=None) for k in h.fields if k.split('fatjet_')[-1] in ['pt', 'eta', 'phi', 'msoftdrop']})
                h.fill(dataset=dataset, **fields)
            elif histname in self.ccfatjet_hists:
                fields = {k: ak.flatten(sccfatjets[k], axis=None) for k in h.fields if k in dir(sccfatjets)}
                fields.update({k: ak.flatten(sccfatjets[k], axis=None) for k in h.fields if k.split('ccfatjet_')[-1] in ['pt', 'eta', 'phi', 'msoftdrop']})
                h.fill(dataset=dataset, **fields)
            else: continue

        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)

        output['njet'].fill(dataset=dataset,  njet=num(sjets))
        output['nbjet'].fill(dataset=dataset, nbjet=num(sbjets))
        output['nel'].fill(dataset=dataset,   nel=num(sel))
        output['nmu'].fill(dataset=dataset,   nmu=num(smu))

        #output['lelpt'].fill(dataset=dataset, lelpt=flatten(selev.Electron[:, 0].pt))
        #output['lmupt'].fill(dataset=dataset, lmupt=flatten(selev.Muon[:, 0].pt))
        #output['ljpt'].fill(dataset=dataset,  ljpt=flatten(selev.Jet[:, 0].pt))
        #output['lfjpt'].fill(dataset=dataset,  ljpt=flatten(selev.FatJet[:, 0].pt))

        return output

    def postprocess(self, accumulator):
        return accumulator
