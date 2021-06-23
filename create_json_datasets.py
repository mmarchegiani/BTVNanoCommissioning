import os
import argparse
import json
import subprocess
import uproot

parser = argparse.ArgumentParser(description='Check for broken files')
parser.add_argument('-i', '--input', default='datasets.txt', help='Input file with datasets lists')
parser.add_argument('-o', '--output', default='', help='Output file with file lists')
parser.add_argument('--outputDir', default='', help='Output directory')
parser.add_argument('-d', '--dataset', default='', help='Dataset to check')
parser.add_argument('-l', '--local', action='store_true', default=False, help='Create json file for local files')
parser.add_argument('--checkFiles', action='store_true', default=False, help='Check if files contain keys')

args = parser.parse_args()

if args.local: prefix='/pnfs/psi.ch/cms/trivcat/store/user/mmarcheg/'
else: prefix = 'root://xrootd-cms.infn.it/'

outDir = args.outputDir if args.outputDir else os.getcwd() + '/datasets/'
if not os.path.exists(outDir): os.makedirs(outDir)
outName = outDir+ ( args.output if args.output else 'datasets_'+args.input.split('/')[-1].replace('.txt', '.json' ) )
if args.local: outName = outName.replace('.json', '_local.json')
if args.dataset: outName = outName.replace('.json', '_'+args.dataset+'.json')

outputjson = {}
inputDatasets = open( args.input, 'r' )
for isam in inputDatasets.read().splitlines():
    if args.dataset:
        if not (args.dataset in isam): continue
    print(f'Adding {isam}')
    cmd = 'dasgoclient --query="'+isam+' instance=prod/phys03 file"'
    allFiles = subprocess.check_output( cmd, shell=True  ).decode("utf-8")
    name = isam.split('/')[1]
    if name.startswith('BTagMu'): name = name+isam.split('/')[2].split('-')[1].split('Run20')[1][2]
    if args.checkFiles:
        outputjson[ name ] = []
        for ifile in allFiles.split('\n'):
            if ifile:
                f = uproot.open(prefix+ifile)['Events'].keys()
                if len(f) > 10: outputjson[ name ].append( prefix+ifile )
    else:
        outputjson[ name ] = [ prefix+ifile for ifile in allFiles.split('\n') if ifile ]

print(f'Saving file {outName}')
outFile = open(outName, 'w')
json.dump(outputjson, outFile, indent=4)
outFile.close()
