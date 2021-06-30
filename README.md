# BTVNanoCommissioning
Repository for Commissioning studies in the BTV POG based on (custom) nanoAOD samples


## Structure
Example worfkflow for ttbar is included. 

Each workflow can be a separete "processor" file, creating the mapping from NanoAOD to
the histograms we need. Workflow processors can be passed to the `runner.py` script 
along with the fileset these should run over. Multiple executors can be chosen 
(for now iterative - one by one, uproot/futures - multiprocessing and dask-slurm). 

To run the example, run:
```
python runner.py --workflow BTVNanoCommissioning
```

Example plots can be found in ` make_some_plots.ipynb` though we might want to make
that more automatic in the end.

## Requirements
### Coffea installation with Miniconda
For installing Miniconda, see also https://hackmd.io/GkiNxag0TUmHnnCiqdND1Q#Local-or-remote
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# Run and follow instructions on screen
bash Miniconda3-latest-Linux-x86_64.sh
```
NOTE: always make sure that conda, python, and pip point to local Miniconda installation (`which conda` etc.).

You can either use the default environment`base` or create a new one:
```
# create new environment from `conda_env.yml` file:
conda env create -f conda_env.yml
# activate environment `btv`:
conda activate btv
# install additional modules required on T3:
pip install bokeh
conda install dask
```
## How to run
### Execution on local machine with Futures Executor
Run the example locally:
```
python runner.py --workflow fattag --executor futures --samples datasets_local_fixed.json --output hists_fattag.coffea --workers 16
```
An alternative method is implemented, submitting dedicated jobs for each dataset, with the option `--splitdataset`:
```
python runner.py --workflow fattag --executor futures --samples datasets_local_fixed.json --output hists_fattag.coffea --workers 16
```
### Execution on Slurm provider with Parsl
Run the example on a Slurm cluster:
```
python runner.py --workflow fattag --executor parsl/slurm --samples datasets_local_fixed.json --output hists_fattag.coffea -s 10
```
### Execution on Condor provider with Parsl
Run the example on a Condor cluster:
```
python runner.py --workflow fattag --executor parsl/condor --samples datasets_local_fixed.json --output hists_fattag.coffea -s 10
```
## Apply corrections 

### Pileup reweighting
To apply pileup reweighting in your MC samples, one need to create first the pileup profile of the MC sample. This can be done with the script `createNTrueForPU.py`. This will take some of the files from your dataset and create a coffea file with that profile. To run it:
```
python createNTrueForPU.py --samples datasets/datasets_btag2017.json --year 2017 
```
The output will be stored in the `correction_files` folder with the name like: `nTrueInt_datasets_btag2017_2017.coffea`. In addition, you need to properly set the `self.puFile` and `self.nTrueFile` in your workflow file.

### Jet energy corrections
You need to download the tar files needed for the JECs from [this twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECdataMC), in the `correction_files/JEC` folder and properly set these names in the `jecTarFiles` of `runner.py` and in `JECversions` of the `utils.py` script.


## Scale Factors

### Installation

There are known conflicts when using a `coffea` conda environment and [combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/). Try to keep two different environments. One to create the histograms in ``coffea`` and another to run `combine`. 

To install combine, follow the instructions from their [documentation](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#for-end-users-that-dont-need-to-commit-or-do-any-development). If you dont have `root` by default in your machine, you might need to install the version inside CMSSW. Otherwise you can use the standalone version.

In addition, you need to download the [rhalphalib](https://github.com/nsmith-/rhalphalib) package
```
pip install --user https://github.com/nsmith-/rhalphalib/archive/master.zip
```

### Create datacards for combine

First we need to translate the coffea histograms into numpy arrays. For this in the coffea environment one can run the script: [coffeaToPickle.py](coffeaToPickle.py), which creates a pickle file with the inputs needed. As example:
```
python coffeaToPickle.py -i histograms/hists_fattag_pileupJEC_2017_WPcuts_v01.coffea7 --year 2017
```
For the rest, migrate to the environment with combine. The script [scaleFactorComputation.py](scaleFactorComputation.py) can create the datacard needed in combine and also run the `FitDiagnostics` step in combine. As example:
```
python scaleFactorComputation.py --year 2017 --tpf histograms/hists_fattag_pileupJEC_2017_WPcuts_v01.pkl --selection msd100tau06DDB
```
