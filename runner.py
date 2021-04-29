import os
import sys
import json
import argparse
import time
import gc

import numpy as np

#import uproot4 as uproot
import uproot
from coffea import hist
from coffea.nanoevents import NanoEventsFactory
from coffea.util import load, save
from coffea import processor
from utils import rescale, lumi, xsecs

def validate(file):
	try:
		fin = uproot.open(file)
		return fin['Events'].num_entries
	except:
		print("Corrupted file: {}".format(file))
		return 


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Run analysis on baconbits files using processor coffea files')
	# Inputs
	parser.add_argument( '--wf', '--workflow', dest='workflow', choices=['ttcom', 'fattag'], help='Which processor to run', required=True)
	parser.add_argument('-o', '--output', default=r'hists.coffea', help='Output histogram filename (default: %(default)s)')
	parser.add_argument('--samples', '--json', dest='samplejson', default='dummy_samples.json', help='JSON file containing dataset and file locations (default: %(default)s)')
	parser.add_argument('--year', type=int, choices=[2016, 2017, 2018], help='Year of data/MC samples', required=True)
	parser.add_argument('--trigger', type=str, choices=['HLT_BTagMu_AK8Jet300_Mu5', 'HLT_BTagMu_AK4Jet300_Mu5'], default=None, help='Single trigger to use', required=False)

	# Scale out
	parser.add_argument('--executor', choices=['iterative', 'futures', 'parsl/condor', 'parsl/slurm', 'dask/condor', 'dask/slurm'], default='futures', help='The type of executor to use (default: %(default)s)')
	parser.add_argument('-j', '--workers', type=int, default=12, help='Number of workers (cores/threads) to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
	parser.add_argument('-s', '--scaleout', type=int, default=6, help='Number of nodes to scale out to if using slurm/condor. Total number of concurrent threads is ``workers x scaleout`` (default: %(default)s)')
	parser.add_argument('--voms', default=None, type=str, help='Path to voms proxy, accsessible to worker nodes. By default a copy will be made to $HOME.')
	parser.add_argument('--splitdataset', action='store_true', help='Process each dataset separately.')

	# Debugging
	parser.add_argument('--validate', action='store_true', help='Do not process, just check all files are accessible')
	parser.add_argument('--skipbadfiles', action='store_true', help='Skip bad files.')
	parser.add_argument('--only', type=str, default=None, help='Only process specific dataset or file')
	parser.add_argument('--limit', type=int, default=None, metavar='N', help='Limit to the first N files of each dataset in sample JSON')
	parser.add_argument('--chunk', type=int, default=500000, metavar='N', help='Number of events per process chunk')
	parser.add_argument('--max', type=int, default=None, metavar='N', help='Max number of chunks to run in total')
	parser.add_argument('--offset', type=int, default=None, metavar='N', help='Offset in JSON reading')
	parser.add_argument('--dataset', type=str, default=None, help='Dataset in the JSON file to process')

	args = parser.parse_args()
	if args.output == parser.get_default('output'):
		args.output = f'hists_{args.workflow}_{(args.samplejson).rstrip(".json")}.coffea'

	# load dataset
	with open(args.samplejson) as f:
		sample_dict = json.load(f)
	
	if args.dataset != parser.get_default('dataset'):
		sample_dict = {args.dataset : sample_dict[args.dataset]}

	if args.offset != parser.get_default('offset'):
		for key in sample_dict.keys():
			sample_dict[key] = sample_dict[key][args.offset:args.offset+args.limit]
	else:
		for key in sample_dict.keys():
			sample_dict[key] = sample_dict[key][:args.limit]

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

	hist_dir = os.getcwd() + "/histograms/"
	if not os.path.exists(hist_dir):
		os.makedirs(hist_dir)

	# Scan if files can be opened
	if args.validate:
		start = time.time()
		from p_tqdm import p_map
		all_invalid = []
		for sample in sample_dict.keys():
			_rmap = p_map(validate, sample_dict[sample], num_cpus=args.workers,
					  desc=f'Validating {sample[:20]}...')
			_results = list(_rmap)
			counts = np.sum([r for r in _results if np.isreal(r)])
			all_invalid += [r for r in _results if type(r) == str]
			print("Events:", np.sum(counts))
		print("Bad files:")
		for fi in all_invalid:
			print(f"  {fi}")
		end = time.time()
		print("TIME:", time.strftime("%H:%M:%S", time.gmtime(end-start)))
		if input("Remove bad files? (y/n)") == "y": 
			print("Removing:")
			for fi in all_invalid:
				print(f"Removing: {fi}")            
				os.system(f'rm {fi}')
		sys.exit(0)
	
	# load workflow
	if args.workflow == "ttcom":
		from workflows.ttbar_validation import NanoProcessor
		processor_instance = NanoProcessor()
	elif args.workflow == "fattag":
		from workflows.fatjet_tagger import NanoProcessor
		processor_instance = NanoProcessor(year=args.year, trigger=args.trigger)
	else:
		raise NotImplemented

	if args.executor not in ['futures', 'iterative']:
		# dask/parsl needs to export x509 to read over xrootd
		if args.voms is not None:
			_x509_path = args.voms
		else:
			_x509_localpath = [l for l in os.popen('voms-proxy-info').read().split("\n") if l.startswith('path')][0].split(":")[-1].strip()
			_x509_path = os.environ['HOME'] + f'/.{_x509_localpath.split("/")[-1]}'
			os.system(f'cp {_x509_localpath} {_x509_path}')

		env_extra = [
			'export XRD_RUNFORKHANDLER=1',
			f'export X509_USER_PROXY={_x509_path}',
			f'export X509_CERT_DIR={os.environ["X509_CERT_DIR"]}',
			'ulimit -u 32768',
		]

		wrk_init = [
			f'export X509_USER_PROXY={_x509_path}',
			f'export X509_CERT_DIR={os.environ["X509_CERT_DIR"]}',
			'source /etc/profile.d/conda.sh',
			'export PATH=$CONDA_PREFIX/bin:$PATH',
			'conda activate btv',
			'cd /afs/cern.ch/work/m/mmarcheg/BTVNanoCommissioning/',
		]

		condor_cfg = '''
		getenv      =  True
		+JobFlavour =  "nextweek"
		'''
		#process_worker_pool = os.environ['CONDA_PREFIX'] + "/bin/process_worker_pool.py"

	#########
	# Execute
	output_split = []
	if args.executor in ['futures', 'iterative']:
		if args.executor == 'iterative':
			_exec = processor.iterative_executor
		else:
			_exec = processor.futures_executor
		if not args.splitdataset:
			output = processor.run_uproot_job(sample_dict,
										treename='Events',
										processor_instance=processor_instance,
										executor=_exec,
										executor_args={
											'skipbadfiles':args.skipbadfiles,
											'schema': processor.NanoAODSchema, 
											'workers': args.workers},
										chunksize=args.chunk, maxchunks=args.max
										)
		else:
			hist_dir = hist_dir + args.output.split(".coffea")[0] + "/"
			if not os.path.exists(hist_dir):
				os.makedirs(hist_dir)
			for dataset in sample_dict.keys():
				output = processor.run_uproot_job({dataset : sample_dict[dataset]},
											treename='Events',
											processor_instance=processor_instance,
											executor=_exec,
											executor_args={
												'skipbadfiles':args.skipbadfiles,
												'schema': processor.NanoAODSchema, 
												'workers': args.workers},
											chunksize=args.chunk, maxchunks=args.max
											)
				filepath = hist_dir + args.output.replace(".coffea", "_" + dataset + ".coffea")
				save(output, filepath)
				print(f"Saving output to {filepath}")
				del output
				#output_split.append(output)
	#elif args.executor == 'parsl/slurm':
	elif 'parsl' in args.executor:
		import parsl
		from parsl.providers import LocalProvider, CondorProvider, SlurmProvider
		from parsl.channels import LocalChannel
		from parsl.config import Config
		from parsl.executors import HighThroughputExecutor
		from parsl.launchers import SrunLauncher, SingleNodeLauncher
		from parsl.addresses import address_by_hostname

		if 'slurm' in args.executor:
			slurm_htex = Config(
				executors=[
					HighThroughputExecutor(
						label="coffea_parsl_slurm",
						address=address_by_hostname(),
						prefetch_capacity=0,
						provider=SlurmProvider(
							channel=LocalChannel(script_dir='logs_parsl'),
							launcher=SrunLauncher(),
							#launcher=SingleNodeLauncher(),
							max_blocks=(args.scaleout)+10,
							init_blocks=args.scaleout,
							partition='wn',
							worker_init="\n".join(env_extra) + "\nexport PYTHONPATH=$PYTHONPATH:$PWD",
							walltime='02:00:00'
						),
					)
				],
				retries=20,
			)
			dfk = parsl.load(slurm_htex)

			if not args.splitdataset:
				output = processor.run_uproot_job(sample_dict,
											treename='Events',
											processor_instance=processor_instance,
											executor=processor.parsl_executor,
											executor_args={
												'skipbadfiles':True,
												'schema': processor.NanoAODSchema, 
												'config': None,
											},
											chunksize=args.chunk, maxchunks=args.max
											)
			else:
				hist_dir = hist_dir + args.output.split(".coffea")[0] + "/"
				if not os.path.exists(hist_dir):
					os.makedirs(hist_dir)
				for dataset in sample_dict.keys():
					print("Processing " + dataset)
					output = processor.run_uproot_job({dataset : sample_dict[dataset]},
												treename='Events',
												processor_instance=processor_instance,
												executor=processor.parsl_executor,
												executor_args={
													'skipbadfiles':True,
													'schema': processor.NanoAODSchema, 
													'config': None,
												},
												chunksize=args.chunk, maxchunks=args.max
												)
					filepath = hist_dir + args.output.replace(".coffea", "_" + dataset + ".coffea")
					save(output, filepath)
					print(f"Saving output to {filepath}")
					del output
					#output_split.append(output)
		elif 'condor' in args.executor:
			#xfer_files = [process_worker_pool, _x509_path]
			#print(xfer_files)

			condor_htex = Config(
				executors=[
					HighThroughputExecutor(
						label="coffea_parsl_slurm",
						#address=address_by_hostname(),
						worker_ports=(8786,8785),
						prefetch_capacity=0,
						provider=CondorProvider(
							channel=LocalChannel(script_dir='logs_parsl'),
							launcher=SingleNodeLauncher(),
							max_blocks=(args.scaleout)+10,
							init_blocks=args.scaleout,
							worker_init="\n".join(wrk_init),
							#transfer_input_files=xfer_files,
							scheduler_options=condor_cfg,
							walltime='00:30:00'
						),
					)
				],
				#retries=20,
			)
			dfk = parsl.load(condor_htex)

			if not args.splitdataset:
				output = processor.run_uproot_job(sample_dict,
											treename='Events',
											processor_instance=processor_instance,
											executor=processor.parsl_executor,
											executor_args={
												'skipbadfiles':True,
												'schema': processor.NanoAODSchema, 
												'config': None,
											},
											chunksize=args.chunk, maxchunks=args.max
											)
			else:
				raise NotImplementedError
		
	elif 'dask' in args.executor:
		from dask_jobqueue import SLURMCluster, HTCondorCluster
		from distributed import Client
		from dask.distributed import performance_report

		if 'slurm' in args.executor:
			cluster = SLURMCluster(
				queue='all',
				cores=args.workers,
				processes=args.workers,
				memory="200 GB",
				retries=10,
				walltime='00:30:00',
				env_extra=env_extra,
			)
		elif 'condor' in args.executor:
			cluster = HTCondorCluster(
				 cores=args.workers, 
				 memory='2GB', 
				 disk='2GB', 
				 env_extra=env_extra,
			)
		cluster.scale(jobs=args.scaleout)

		client = Client(cluster)
		with performance_report(filename="dask-report.html"):
			output = processor.run_uproot_job(sample_dict,
										treename='Events',
										processor_instance=processor_instance,
										executor=processor.dask_executor,
										executor_args={
											'client': client,
											'skipbadfiles':args.skipbadfiles,
											'schema': processor.NanoAODSchema, 
										},
										chunksize=args.chunk, maxchunks=args.max
							)

	if not args.splitdataset:
		if args.offset == parser.get_default("offset"):
			if len(sample_dict.keys()) > 1:
				output = rescale(output, xsecs, lumi[args.year])
				#output = rescale(output, xsecs, lumi[args.year], "JetHT")
			save(output, hist_dir + args.output)
			print(output)
			print(f"Saving output to {hist_dir + args.output}")
		else:
			# In this case the MC is not rescaled yet
			print("No MC rescaling applied")			
			hist_dir = hist_dir + args.output.split(".coffea")[0] + "/"
			if not os.path.exists(hist_dir):
				os.makedirs(hist_dir)
			args.output = args.output.replace(".coffea", "_0" + str(args.offset) + ".coffea")
			save(output, hist_dir + args.output)
			print(output)
			print(f"Saving output to {hist_dir + args.output}")
		
	else:
		files_list = [file for file in os.listdir(hist_dir) if file != args.output]
		#accumulator = output_split[0]
		accumulator = load(hist_dir + files_list[0])
		histograms = accumulator.keys()
		for histname in histograms:
			for file in files_list[1:]:
				output = load(hist_dir + file)
				accumulator[histname].add(output[histname])
				del output

		if not os.path.exists(hist_dir):
			os.makedirs(hist_dir)
		if len(sample_dict.keys()) > 1:
			accumulator = rescale(accumulator, xsecs, lumi[args.year])
		save(accumulator, hist_dir + args.output)
		print(accumulator)
		print(f"Saving output to {hist_dir + args.output}")

