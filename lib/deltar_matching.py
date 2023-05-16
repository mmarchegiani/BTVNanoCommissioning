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
