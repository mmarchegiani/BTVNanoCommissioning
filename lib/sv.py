import sys

import collections
import numpy as np
import awkward as ak

def project(a, b):
    return a.dot(b)/b.dot(b) * b

def get_nsv(jet, sv, pos, R=0.4):

    # Compute the number of SV inside the subjet of the leading and subleading FatJet (nsv1, nsv2)
    # and return the concatenated array
    if pos == 0:
        jet = jet[:,0]
    elif pos == 1:
        jet = ak.pad_none(jet,2)[:,1]
    else:
        raise Exception("Only the leading and subleading jets can be considered.")

    sj1 = jet.subjets[:, 0]
    sj2 = jet.subjets[:, 1]
    sv_dr1 = sj1.delta_r(sv.p4)
    sv_dr2 = sj2.delta_r(sv.p4)
    nsv1 = ak.unflatten( ak.count(sv_dr1[sv_dr1 < R], axis=1), counts=1 )
    nsv2 = ak.unflatten( ak.count(sv_dr2[sv_dr2 < R], axis=1), counts=1 )

    return ak.concatenate((nsv1, nsv2), axis=1)

def get_sv_in_jet(jet, sv, pos, R=0.8):

    if pos == 0:
        jet = jet[:,0]
    elif pos == 1:
        jet = ak.pad_none(jet,2)[:,1]
    else:
        raise Exception("Only the leading and subleading jets can be considered.")

    sv_dr = jet.delta_r(sv.p4)
    sv_in_jet = sv_dr < R
    empty = ak.Array([ak.count(sv_in_jet)*[]])
    sv_in_jet = ak.where(ak.is_none(sv_in_jet), empty, sv_in_jet)

    return sv_in_jet

# N.B.: In the following the logarithm of the mass-like variables is set to -5 as default value,
# when the corresponding mass value is 0. This way, the log(mass) histograms will be filled with
# -5 when the mass is 0. In any case, for the final fit only the range above -2.5 is considered.

def get_summass(sv, log=True):

    nsv = ak.count(sv.pt, axis=1)
    summass = sv.p4.sum().mass

    if log:
        logsummass = ak.where(nsv < 1, -5, np.log(summass))
        return summass, logsummass
    else:
        return summass

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

def get_sv1mass(sv, log=True):

    nsv = ak.count(sv.pt, axis=1)
    sv1mass = ak.firsts(sv).mass

    if log:
        logsv1mass = ak.where(nsv < 1, -5, np.log(sv1mass))
        return sv1mass, logsv1mass
    else:
        return sv1mass

def get_sumcorrmass(sv, log=True):
    
    nsv = ak.count(sv.pt, axis=1)
    corrmass = np.sqrt(sv.p4.mass**2 + sv.p4.pt**2 * np.sin(sv.pAngle)**2) + sv.p4.pt * np.sin(sv.pAngle)
    sv['mass'] = corrmass
    sumcorrmass = sv.p4.sum().mass

    if log:
        logsumcorrmass = ak.where(nsv < 1, -5, np.log(sumcorrmass))
        return sumcorrmass, logsumcorrmass
    else:
        return sumcorrmass
