import awkward as ak

def run_deltar_matching(obj1, obj2, radius=0.4): # NxM , NxG arrays
    '''
    Doing this you can keep the assignment on the obj2 collection unique, 
    but you are not checking the uniqueness of the matching to the first collection. 
    '''
    _, obj2 = ak.unzip(ak.cartesian([obj1, obj2], nested=True)) # Obj2 is now NxMxG
    obj2['dR'] = obj1.delta_r(obj2)  # Calculating delta R
    t_index = ak.argmin(obj2.dR, axis=-2) # Finding the smallest dR (NxG array)
    s_index = ak.local_index(obj1.eta, axis=-1) #  NxM array
    _, t_index = ak.unzip(ak.cartesian([s_index, t_index], nested=True)) 
    obj2 = obj2[s_index == t_index] # Pairwise comparison to keep smallest delta R
    # Cutting on delta R
    obj2 = obj2[obj2.dR < radius] #Additional cut on delta R, now a NxMxG' array 
    return obj2

def muons_matched_to_fatjet(events):
    '''This function returns the collection of muons matched to the fatjets.
    The output array has the same shape as the events.FatJetGood collection.
    '''
    return run_deltar_matching(events.FatJetGood, events.MuonGood, radius=0.8)

def muons_matched_to_subjet(events, pos, unique=True):
    '''This function returns the collection of muons matched to the subjet in the position `pos` contained in the AK8 jet.
    If pos=0, the muons matched to the leading subjet are returned.
    If pos=1, the muons matched to the leading subjet are returned.
    The output array has the same shape as the events.FatJetGood collection.
    '''
    R = 0.4
    sj = events.FatJetGood.subjets[:,:,pos]
    if unique:
        fatjet = ak.pad_none(events.FatJetGood, 2)
        sj1 = events.FatJetGood.subjets[:,:,0]
        sj2 = events.FatJetGood.subjets[:,:,1]
        dr_sj1_sj2 = sj1.delta_r(sj2)
        dr_non_overlapping_cone = 0.5 * dr_sj1_sj2
        radius = ak.where(dr_non_overlapping_cone > R, R, dr_non_overlapping_cone)
    else:
        radius = R
    return run_deltar_matching(sj, events.MuonGood, radius=radius)

def _muon_match(jet, muon, pos, R=0.4, unique=True):
    '''This function returns the mask of muons matched to the one of the subjet.
    The index of the subjet is specified with the argument `pos` (leading: 0, sub-leading: 1).
    The matching is unique since muons that are contained in both subjets are excluded,
    by varying the deltaR threshold as min(R, 0.5*deltaR(subjet_1, subjet_2))'''

    if pos == 0:
        jet = jet[:,0]
    elif pos == 1:
        jet = ak.pad_none(jet,2)[:,1]
    else:
        raise Exception("Only the leading and subleading jets can be considered.")

    sj1 = jet.subjets[:, 0]
    sj2 = jet.subjets[:, 1]
    dr_sj1_sj2 = sj1.delta_r(sj2)
    if unique:
        dr_non_overlapping_cone = 0.5 * dr_sj1_sj2
        dr_threshold = ak.where(dr_non_overlapping_cone > R, R, dr_non_overlapping_cone)
    else:
        dr_threshold = R
    dr_mu_sj1 = sj1.delta_r(muon)
    dr_mu_sj2 = sj2.delta_r(muon)

    mask_matched_muon_sj1 = ak.fill_none(dr_mu_sj1 < dr_threshold, [], axis=0)
    mask_matched_muon_sj2 = ak.fill_none(dr_mu_sj2 < dr_threshold, [], axis=0)

    return mask_matched_muon_sj1, mask_matched_muon_sj2

def _muons_matched_to_subjets_of_leadfatjet(jet, muon, R=0.4, unique=True):
    return _muon_match(jet, muon, pos=0, R=R, unique=unique)

def _muons_matched_to_subjets_of_subleadfatjet(jet, muon, R=0.4, unique=True):
    return _muon_match(jet, muon, pos=1, R=R, unique=unique)

def muon_matched_selection(jet, muon, R=0.4, unique=True):
    mask_muon_leadfatjet_sj1, mask_muon_leadfatjet_sj2 = _muons_matched_to_subjets_of_leadfatjet(jet, muon, R=R, unique=unique)
    mask_muon_subleadfatjet_sj1, mask_muon_subleadfatjet_sj2 = _muons_matched_to_subjets_of_subleadfatjet(jet, muon, R=R, unique=unique)
    assert ~ak.any(mask_muon_sj1_1 & mask_muon_sj2_1), "The muon matching of subjets inside the leading fatjet is not unique. Please review the implementation of the matching."
    assert ~ak.any(mask_muon_sj1_2 & mask_muon_sj2_2), "The muon matching of subjets inside the subleading fatjet is not unique. Please review the implementation of the matching."

    muons_in_sj1_1 = muon[mask_muon_sj1_1]
    muons_in_sj2_1 = muon[mask_muon_sj2_1]
    muons_in_sj1_2 = muon[mask_muon_sj1_2]
    muons_in_sj2_2 = muon[mask_muon_sj2_2]
    dimuon_1 = ak.firsts(muons_in_sj1_1) + ak.firsts(muons_in_sj2_1)
    dimuon_2 = ak.firsts(muons_in_sj1_2) + ak.firsts(muons_in_sj2_2)

    muons_in_sj1_1_unflatten = ak.unflatten(muons_in_sj1_1, counts=1)
    muons_in_sj2_1_unflatten = ak.unflatten(muons_in_sj2_1, counts=1)
    muons_in_sj1_2_unflatten = ak.unflatten(muons_in_sj1_2, counts=1)
    muons_in_sj2_2_unflatten = ak.unflatten(muons_in_sj2_2, counts=1)
    dimuon_1_unflatten = ak.unflatten(dimuon_1, counts=1)
    dimuon_2_unflatten = ak.unflatten(dimuon_2, counts=1)
    # Get the dimuon collection. One dimuon pair per jet
    dimuon = ak.concatenate((dimuon_1_unflatten, dimuon_2_unflatten), axis=1)
    # Broadcast dimuon collection to have the same shape as jet
    dimuon = dimuon[ak.local_index(jet, axis=1)]

    # Get the muon collections: muons in the leading and subleading subjets
    muons_in_sj1 = ak.concatenate((muons_in_sj1_1_unflatten, muons_in_sj1_2_unflatten), axis=1)
    muons_in_sj2 = ak.concatenate((muons_in_sj2_1_unflatten, muons_in_sj2_2_unflatten), axis=1)

    breakpoint()

    return muons_in_sj1, muons_in_sj2, dimuon

def get_nmu_in_subjet(jet, muon, pos, R=0.4):

    # Compute the number of muons inside the subjet of the leading and subleading FatJet (nmusj1, nmusj2).
    # The concatenated array is returned.
    if pos == 0:
        jet = jet[:,0]
    elif pos == 1:
        jet = ak.pad_none(jet,2)[:,1]
    else:
        raise Exception("Only the leading and subleading jets can be considered.")

    sj1 = jet.subjets[:, 0]
    sj2 = jet.subjets[:, 1]
    dr12 = sj1.delta_r(sj2)
    mu_dr1 = sj1.delta_r(muon)
    mu_dr2 = sj2.delta_r(muon)
    nmusj1_flat = ak.fill_none(ak.count(mu_dr1[mu_dr1 < R], axis=1), 0, axis=0)
    nmusj2_flat = ak.fill_none(ak.count(mu_dr2[mu_dr2 < R], axis=1), 0, axis=0)
    nmusj1 = ak.unflatten( nmusj1_flat, counts=1 )
    nmusj2 = ak.unflatten( nmusj2_flat, counts=1 )

    assert not ak.any(ak.is_none(nmusj1, axis=1)), nmusj1[ak.any(ak.is_none(nmusj1, axis=1), axis=1)]
    assert not ak.any(ak.is_none(nmusj2, axis=1)), nmusj2[ak.any(ak.is_none(nmusj2, axis=1), axis=1)]

    return ak.concatenate((nmusj1, nmusj2), axis=1)
