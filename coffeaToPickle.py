import argparse
import numpy as np
from coffea.util import load
from coffea.hist import plot
import coffea.hist as hist
import itertools
import os
import pickle
from utils import histogram_settings, lumi, rescale, xsecs

parser = argparse.ArgumentParser(description='Plot histograms from coffea file')
parser.add_argument('-i', '--input', type=str, help='Input histogram filename', required=True)
parser.add_argument('-o', '--output', type=str, default='', help='Output file')
parser.add_argument('--outputDir', type=str, default=None, help='Output directory')
parser.add_argument('--year', type=int, choices=[2016, 2017, 2018], help='Year of data/MC samples', required=True)
parser.add_argument('--only', action='store', default='', help='Plot only one histogram')
parser.add_argument('--test', action='store_true', default=False, help='Test with lower stats.')
parser.add_argument('--data', type=str, default='BTagMu', help='Data sample name')
parser.add_argument('--selection', type=str, default='all', help='Plot only plots with this selection. ("all" to plot all the selections in file)')

args = parser.parse_args()

if os.path.isfile( args.input ): accumulator = load(args.input)
else:
    files_list = [ifile for ifile in os.listdir(args.input) if ifile != args.output]
    accumulator = load(args.input + files_list[0])
    histograms = accumulator.keys()
    for ifile in files_list[1:]:
        output = load(args.input + ifile)
        for histname in histograms:
            accumulator[histname].add(output[histname])

scaleXS = {}
for isam in accumulator[next(iter(accumulator))].identifiers('dataset'):
    isam = str(isam)
    scaleXS[isam] = 1 if isam.startswith('BTag') else xsecs[isam]/accumulator['sumw'][isam]

outputDict = {}
for ivar in [ 'fatjet_jetproba' ]:
    for isel in [ 'msd100tau06' ]:
        for DDX in [ 'DDB', 'DDC' ]:
            for wp in [ 'M' ]:
                for passfail in ['pass', 'fail']:

                    histname=f'{ivar}_{isel}{DDX}{passfail}{wp}wp'
                    h = accumulator[histname]
                    h.scale( scaleXS, axis='dataset' )
                    h = h.rebin(h.fields[-1], hist.Bin(h.fields[-1], h.axis(h.fields[-1]).label, **histogram_settings['variables']['_'.join(histname.split('_')[:-1])]['binning']))

                    ##### grouping data and QCD histos
                    datasets = [str(s) for s in h.axis('dataset').identifiers() if str(s) != 'dataset']
                    mapping = {
                        r'QCD ($\mu$ enriched)' : [dataset for dataset in datasets if 'QCD_Pt' in dataset],
                        r'BTagMu': [ idata for idata in datasets if args.data in idata ],
                    }
                    datasets = mapping.keys()
                    datasets_data  = [dataset for dataset in datasets if args.data in dataset]
                    datasets_QCD = [dataset for dataset in datasets if ((args.data not in dataset) & ('GluGlu' not in dataset))]

                    h = h.group("dataset", hist.Cat("dataset", "Dataset"), mapping)

                    #### rescaling QCD to data
                    dataSum = np.sum( h[args.data].sum('flavor').values()[('BTagMu',)] )
                    QCDSum = np.sum( h[datasets_QCD].sum('dataset', 'flavor').values()[()] )
                    QCD_rescaled = h[datasets_QCD].sum('dataset')
                    QCD_rescaled.scale( dataSum/QCDSum )

                    #### storing into dict
                    for iflav in QCD_rescaled.values():
                        tmpValue, sumw2 = QCD_rescaled[iflav].sum('flavor').values(sumw2=True)[()]
                        #outputDict[ histname+'_QCD_'+iflav[0] ] = np.array( value )
                        outputDict[ histname+'_QCD_'+iflav[0] ] = [ tmpValue, sumw2  ]
                    tmpValue, sumw2 = h[args.data].sum('flavor').values(sumw2=True)[('BTagMu',)]
                    outputDict[ histname+'_BtagMu' ] = [ tmpValue, sumw2 ]

#### Saving into pickle
output_dir = args.outputDir if args.outputDir else os.getcwd()+"/histograms/" + args.output + "/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
outputFileName = output_dir + ( args.output if args.output else args.input.split('/')[-1].replace('coffea7', 'pkl')  )
outputFile = open( outputFileName, 'wb'  )
pickle.dump( outputDict, outputFile, protocol=2 )
outputFile.close()

