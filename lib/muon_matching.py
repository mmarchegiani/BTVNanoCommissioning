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

def muon_matched_to_subjet(events, pos, unique=True):
    '''This function returns the collection of muons matched to the subjet in the position `pos` contained in the AK8 jet.
    If pos=0, the muons matched to the leading subjet are returned.
    If pos=1, the muons matched to the leading subjet are returned.
    The output array has the same shape as the events.FatJetGood collection.
    '''
    R = 0.4
    sj = events.FatJetGood.subjets[:,:,pos]
    if unique:
        sj1 = events.FatJetGood.subjets[:,:,0]
        sj2 = events.FatJetGood.subjets[:,:,1]
        dr_sj1_sj2 = sj1.delta_r(sj2)
        dr_non_overlapping_cone = 0.5 * dr_sj1_sj2
        radius = ak.where(dr_non_overlapping_cone > R, R, dr_non_overlapping_cone)
    else:
        radius = R

    # This collection of muons will contain all the muons contained within the dR cone
    muons_matched = run_deltar_matching(sj, events.MuonGood, radius=radius)

    # Of all the muons contained in the dR cone, we only consider the leading muon to be matched to the subjet
    # N.B.: the slicing syntax `[:,:,None]` is needed in order for the output array to have a 3 dimensions
    return ak.firsts(muons_matched, axis=2)[:,:,None]
