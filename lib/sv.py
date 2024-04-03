import sys

import collections
import numpy as np
import awkward as ak

from .deltar_matching import run_deltar_matching

def project(a, b):
    return a.dot(b)/b.dot(b) * b

def sv_matched_to_fatjet(events):
    '''This function returns the collection of SV matched to the fatjets.
    The output array has the same shape as the events.FatJetGood collection.
    '''
    return run_deltar_matching(events.FatJetGood, events.SV, radius=0.8)

# N.B.: In the following the logarithm of the mass-like variables is set to -5 as default value,
# when the corresponding mass value is 0. This way, the log(mass) histograms will be filled with
# -5 when the mass is 0. In any case, for the final fit only the range above -2.5 is considered.

def get_projmass(jet, sv, pos, log=True):

    if pos == 0:
        jet = jet[:,0]
    elif pos == 1:
        jet = ak.pad_none(jet,2)[:,1]
    else:
        raise Exception("Only the leading and subleading jets can be considered.")
    
    nsv = ak.count(sv.pt, axis=1)
    projmass = project(sv.p4.sum(), jet).mass

    if log:
        logprojmass = ak.where(nsv < 1, -5, np.log(projmass))
        return projmass, logprojmass
    else:
        return projmass

"""
def get_sv1mass(sv, log=True):

    nsv = ak.count(sv.pt, axis=1)
    sv1mass = ak.firsts(sv).mass

    if log:
        logsv1mass = ak.where(nsv < 1, -5, np.log(sv1mass))
        return sv1mass, logsv1mass
    else:
        return sv1mass
"""

def get_corrmass(sv):

    corrmass = np.sqrt(sv.p4.mass**2 + sv.p4.p**2 * np.sin(sv.pAngle)**2) + sv.p4.p * np.sin(sv.pAngle)

    return corrmass

def get_sv1mass(sv, log=True):

    sv1mass = ak.firsts(sv.mass, axis=-1)[:,0]
    sv1mass = ak.where(ak.is_none(sv1mass), 0, sv1mass)

    if log:
        logsv1mass = np.where(sv1mass == 0, -5, np.log(sv1mass))
        return sv1mass, logsv1mass
    else:
        return sv1mass


#def get_sumcorrmass(sv, log=True):
def get_sumcorrmass(sv, log=True):

    sumcorrmass = sv.p4.sum().mass

    if log:
        logsumcorrmass = np.log(sumcorrmass)
        return sumcorrmass, logsumcorrmass
    else:
        return sumcorrmass
