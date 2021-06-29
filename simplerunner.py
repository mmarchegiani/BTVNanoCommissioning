import os
import sys
import json
import argparse
import time
import tarfile
import tempfile

import numpy as np

import uproot
from coffea.util import load, save
from coffea import processor
from coffea.nanoevents import PFNanoAODSchema
from distributed import Client
from dask_jobqueue import HTCondorCluster
import socket


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
    parser.add_argument('--outputDir', type=str, default=None, help='Output directory')
    parser.add_argument('--nTrueFile', type=str, default='', help='To specify nTrue file. To use the default leave it empty')

    # Scale out
    parser.add_argument('--executor', choices=['iterative', 'futures', 'parsl/condor', 'parsl/slurm', 'dask/condor', 'dask/slurm'], default='futures', help='The type of executor to use (default: %(default)s)')
    parser.add_argument('-j', '--workers', type=int, default=12, help='Number of workers (cores/threads) to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    parser.add_argument('-s', '--scaleout', type=int, default=6, help='Number of nodes to scale out to if using slurm/condor. Total number of concurrent threads is ``workers x scaleout`` (default: %(default)s)')

    # Debugging
    parser.add_argument('--validate', action='store_true', help='Do not process, just check all files are accessible')
    parser.add_argument('--skipbadfiles', action='store_true', help='Skip bad files.')
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

    hist_dir = os.getcwd() + "/histograms/" if args.outputDir is None else args.outputDir
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

    ##### Untar JECs
    ##### Correction files in https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
    jesInputFilePath = tempfile.mkdtemp()
    if args.year==2017:
        jecTarFiles = [
            '/correction_files/JEC/Fall17_17Nov2017B_V32_DATA.tar.gz',
            '/correction_files/JEC/Fall17_17Nov2017C_V32_DATA.tar.gz',
            '/correction_files/JEC/Fall17_17Nov2017DE_V32_DATA.tar.gz',
            '/correction_files/JEC/Fall17_17Nov2017F_V32_DATA.tar.gz',
            '/correction_files/JEC/Fall17_17Nov2017_V32_MC.tar.gz',
            ]
    if args.year==2018:
        jecTarFiles = [
            '/correction_files/JEC/Autumn18_RunA_V19_DATA.tar.gz',
            '/correction_files/JEC/Autumn18_RunB_V19_DATA.tar.gz',
            '/correction_files/JEC/Autumn18_RunC_V19_DATA.tar.gz',
            '/correction_files/JEC/Autumn18_RunD_V19_DATA.tar.gz',
            '/correction_files/JEC/Autumn18_V19_MC.tar.gz',
            ]
    for itar in jecTarFiles:
        jecFile = os.getcwd()+itar
        jesArchive = tarfile.open( jecFile, "r:gz")
        jesArchive.extractall(jesInputFilePath)

    # load workflow
    if args.workflow == "ttcom":
        from workflows.ttbar_validation import NanoProcessor
        processor_instance = NanoProcessor()
    elif args.workflow == "fattag":
        from workflows.fatjet_tagger import NanoProcessor
        processor_instance = NanoProcessor(year=args.year, JECfolder=jesInputFilePath, nTrueFile=args.nTrueFile)
    else:
        raise NotImplemented


    #########
    # Execute
    _exec_args = {
                    'skipbadfiles':args.skipbadfiles,
                    'schema': PFNanoAODSchema,
                    'workers': args.workers}
    if args.executor.startswith('iterative'): _exec = processor.iterative_executor
    elif args.executor.startswith('futures'): _exec = processor.futures_executor
    elif args.executor.startswith('dask/condor'):
        _exec = processor.dask_executor
        n_port = 8786
        cluster = HTCondorCluster(
                    cores=1,
                    memory='2000MB',
                    disk='1000MB',
                    death_timeout = '60',
                    nanny = False,
                    scheduler_options={
                        'port': n_port,
                        'host': socket.gethostname()
                        },
                    job_extra={
                        'log': 'dask_job_output.log',
                        'output': 'dask_job_output.out',
                        'error': 'dask_job_output.err',
                        'should_transfer_files': 'Yes',
                        'when_to_transfer_output': 'ON_EXIT',
                        '+JobFlavour': '"tomorrow"',
                        },
                    extra = ['--worker-port {}'.format(n_port)]
                )
        cluster.scale(jobs=args.scaleout)
        print(cluster.job_script())
        client = Client(cluster)
        _exec_args = {
                    'schema': PFNanoAODSchema,
                    'align_clusters' : True,
                    'client': client
                    }

    output = processor.run_uproot_job(sample_dict,
                                                treename='Events',
                                                processor_instance=processor_instance,
                                                executor=_exec,
                                                executor_args=_exec_args,
                                                chunksize=args.chunk, maxchunks=args.max
                                                )
    filepath = hist_dir + args.output
    save(output, filepath)
    print(f"Saving output to {filepath}")
