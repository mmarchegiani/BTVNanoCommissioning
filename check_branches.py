import argparse
import uproot
import json
import sys

parser = argparse.ArgumentParser(description='Check for broken files')
parser.add_argument('-i', '--input', type=str, default='datasets_local.json', help='Input file with file lists')
parser.add_argument('-d', '--dataset', type=str, default='BTagMu', help='Dataset to check')

args = parser.parse_args()

if not args.input.endswith(".json"):
	sys.exit("Only json files allowed as input.")
files = []
with open(args.input) as f:
	files = json.load(f)
files_to_check = files[args.dataset]

keys_broken = []
for file in files_to_check:
	f = uproot.open(file)
	keys = f['Events'].keys()
	nkeys = len(keys)
	if nkeys > 10:
		print('"'+file+'",')
		keys_broken.append(keys)

#for item in keys_broken:
#	print(item)
