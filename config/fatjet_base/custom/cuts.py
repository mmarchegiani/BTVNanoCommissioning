# Per-event cuts applied to each event
from pocket_coffea.lib.cut_definition import Cut
from config.fatjet_base.custom.functions import twojets_ptmsd, mutag_fatjet, mutag_subjet, ptbin, ptbin_mutag, msoftdrop, ptmsd, ptmsdtau, min_nObj_minmsd, flavor_mask

def twojets_presel(pt, msd, name=None):
    if name == None:
        name = f"twojets_pt{pt}msd{msd}"
    return Cut(
    name=name,
    params={
        "pt" : pt,
        "msd" : msd,
    },
    function=twojets_ptmsd
)

def mutag_fatjet_sel(nmu, name=None):
    if name == None:
        name = f"mutag_fatjet_nmu-{nmu}"
    return Cut(
        name=name,
        params={
            "nmu" : nmu,
        },
        collection="FatJetGood",
        function=mutag_fatjet
    )

def mutag_subjet_sel(unique_matching, name=None):
    if name == None:
        name = "mutag_subjet"
        if unique_matching:
            name += "_unique"
    return Cut(
        name=name,
        params={
            "nsubjet" : 2,
            "nmusj" : 1,
            "unique_matching": unique_matching
        },
        collection="FatJetGood",
        function=mutag_subjet
    )

def get_ptbin(pt_low, pt_high, name=None):
    if name == None:
        name = f"Pt-{pt_low}to{pt_high}"
    return Cut(
        name=name,
        params= {"pt_low" : pt_low, "pt_high" : pt_high},
        function=ptbin,
        collection="FatJetGood"
    )

def get_ptbin_mutag(pt_low, pt_high, name=None):
    if name == None:
        name = f"Pt-{pt_low}to{pt_high}"
    return Cut(
        name=name,
        params= {"pt_low" : pt_low,
                 "pt_high" : pt_high,
                 "nsubjet" : 2,
                 "nmusj" : 1,
                 "dimuon_pt_ratio": 0.6
        },
        function=ptbin_mutag,
        collection="FatJetGood"
    )

def get_msd(msd, name=None):
    if name == None:
        name = f"msd{msd}"
    return Cut(
        name=name,
        params= {"msd" : msd},
        function=msoftdrop,
        collection="FatJetGood"
    )

def get_ptmsd(pt, msd, name=None):
    if name == None:
        name = f"pt{pt}msd{msd}"
    return Cut(
        name=name,
        params= {"pt" : pt, "msd" : msd},
        function=ptmsd,
        collection="FatJetGood"
    )

def get_ptmsdtau(pt, msd, tau21, name=None):
    if name == None:
        name = f"msd{msd}tau{tau21}"
    return Cut(
        name=name,
        params= {"pt" : pt, "msd" : msd, "tau21" : tau21},
        function=ptmsdtau
    )

def get_nObj_minmsd(N, minmsd=None, coll="JetGood", name=None):
    '''
    Factory function which creates a cut for minimum number of objects.
    Optionally a minimum msd is requested.
    :param N: request >= N objects
    :param coll: collection to use
    :param minmsd: minimum msd
    :param name: name for the cut, by defaul it is built as n{coll}_min{N}_msd{minmsd}
    :returns: a Cut object
    '''
    if name == None:
        if minmsd:
            name = f"n{coll}_min{N}_msd{minmsd}"
        else:
            name = f"n{coll}_min{N}"
    if minmsd:
        return Cut(
            name=name,
            params={"N": N, "coll": coll, "minmsd": minmsd},
            function=min_nObj_minmsd,
        )
    else:
        raise NotImplementedError
        #return Cut(name=name, params={"N": N, "coll": coll}, function=min_nObj)

def get_flavor(flavor):
    return Cut(
        name=flavor,
        params={"flavor": flavor},
        function=flavor_mask,
        collection="FatJetGood"
    )
