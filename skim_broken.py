import argparse
import uproot
import json
import sys
import os

parser = argparse.ArgumentParser(description='Save json files without broken files')
parser.add_argument('-i', '--input', type=str, default='datasets_local.json', help='Input file with file list')
parser.add_argument('-b', '--broken', type=str, default='datasets_broken.json', help='Input file with broken files list')
parser.add_argument('-o', '--output', type=str, default='datasets_fixed.json', help='Output file without broken files')

args = parser.parse_args()

if ((not args.input.endswith(".json")) | (not args.broken.endswith(".json")) | (not args.output.endswith(".json"))):
	sys.exit("Only json files allowed as input/output.")
if (os.path.exists(args.output)):
	sys.exit(f"Output file {args.output} is already existing.")
files = None
brokenfiles = None
with open(args.input) as f:
	files = json.load(f)
with open(args.broken) as g:
	brokenfiles = json.load(g)
f.close()
g.close()

for key in brokenfiles.keys():
	filestoskim = []
	for f in files[key]:
		if f in brokenfiles[key]:
			filestoskim.append(f)
			#print(f"Removing {f}")
	files[key] = [f for f in files[key] if not f in filestoskim]

with open(args.output, 'w') as outfile:
    json.dump(files, outfile, sort_keys=True, indent=4)
