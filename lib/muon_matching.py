import awkward as ak
from .deltar_matching import run_deltar_matching

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
