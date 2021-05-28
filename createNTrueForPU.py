#####################################
##### This script sum the genWeights from nanoAOD
##### To run: python computeGenWeights.py -d NAMEOFDATASET
#####################################

#!/usr/bin/env python
import argparse, os, shutil, sys
import numpy as np
import uproot
import coffea.hist as hist
from coffea.util import save
import concurrent.futures
import json
#from root_numpy import root2array, tree2array



#####################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--samples', '--json', dest='samplejson', default='dummy_samples.json', help='JSON file containing dataset and file locations (default: %(default)s)')
    parser.add_argument('--year', type=int, choices=[2016, 2017, 2018], help='Year of data/MC samples', required=True)
    parser.add_argument('--only', type=str, default=None, help='Only process specific dataset or file')
    parser.add_argument('--dataset', type=str, default=None, help='Dataset in the JSON file to process')

    try: args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    # load dataset
    with open(args.samplejson) as f:
        sample_dict = json.load(f)

    if args.dataset != parser.get_default('dataset'):
        sample_dict = {args.dataset : sample_dict[args.dataset]}

    # For debugging
    if args.only is not None:
        if args.only in sample_dict.keys():  # is dataset
            sample_dict = dict([(args.only, sample_dict[args.only])])
        if "*" in args.only: # wildcard for datasets
            _new_dict = {}
            print("Will only proces the following datasets:")
            for k, v in sample_dict.items():
                if k.lstrip("/").startswith(args.only.rstrip("*")):
                    print("    ", k)
                    _new_dict[k] = v
            sample_dict = _new_dict
        else:  # is file
            for key in sample_dict.keys():
                if args.only in sample_dict[key]:
                    sample_dict = dict([(key, [args.only])])

    if len(sample_dict)==0: print ('No sample found. \n Have a nice day :)')

    puHisto = hist.Hist( 'nTrueInt', hist.Cat('dataset', 'dataset name'), hist.Bin('x', 'nTrueInt', 99, 0, 99) )

    executor = concurrent.futures.ThreadPoolExecutor()
    for isample, ifiles in sample_dict.items():
        print(f'processing {isample}')
        if isample.startswith('BTag'): continue
        if not ifiles:
            print('|--------> Sample '+isample+' does not have files to process')
            continue
        filesToProcess = { k:'Events' for k in ifiles[:10] }

        for array in uproot.iterate([filesToProcess], ['Pileup_nTrueInt']):
            puHisto.fill(dataset=isample, x=array['Pileup_nTrueInt'])
    save( puHisto, os.getcwd()+'/correction_files/nTrueInt_'+args.samplejson.split('/')[-1].split('.json')[0]+'_'+str(args.year)+'.coffea' )


