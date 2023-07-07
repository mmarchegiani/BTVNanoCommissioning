import os
import sys
import json

import numpy as np
import pandas as pd
import rhalphalib as rl
import uproot
import ROOT

from parameters_fit import fit_parameters
from parameters import AK8Taggers, AK8Taggers_bb, AK8Taggers_cc

regions = ['pass', 'fail']
unc_lumi = {
    '2016_PreVFP'  : 1.012,
    '2016_PostVFP' : 1.012,
    '2017'         : 1.023,
    '2018'         : 1.025
}

epsilon = 1e-3

def merge_bins(values, rebin_size):
    values_rebinned = []
    for i in np.arange(0, len(values), rebin_size):
        values_sum = 0
        for j in range(rebin_size):
            values_sum += values[i+j]
        values_rebinned.append(values_sum)

    return np.array(values_rebinned)

def iserror(func, *args, **kw):
    try:
        func(*args, **kw)
        return False
    except Exception:
        return True

class Fit():
    def __init__(self, input, output, categories, var, year, xlim, binwidth, scheme='3f'):
        self.input = input
        self.output = output
        self.scheme = scheme
        self.year = year
        self.var = var
        self.lo = xlim[0]
        self.hi = xlim[1]
        self.binwidth = binwidth
        self.define_flavors()
        self.templates = np.load(os.path.abspath(self.input), allow_pickle=True)
        self.categories = categories
        """
        if self.scheme == '3f':
            self.fitResults = {
                'b+bb' : os.path.join(os.getcwd(), self.output, 'fitResults_bb.csv'),
                'c+cc' : os.path.join(os.getcwd(), self.output, 'fitResults_cc.csv'),
            }
        elif self.scheme == '5f':
            self.fitResults = {
                'bb' : os.path.join(os.getcwd(), self.output, 'fitResults_bb.csv'),
                'cc' : os.path.join(os.getcwd(), self.output, 'fitResults_cc.csv'),
            }
        """
        #self.fitResults = os.path.join(os.getcwd(), self.output)
        self.signal_name = {}
        self.define_bins()
        self.define_observable()
        self.initialize_models_dict()
        self.load_parameters()
        self.define_independent_parameters()
        self.define_nuisance_parameters()
        self.define_frozen_parameters()
        self.build_fit_models()
        self.build_combine_script()
        self.build_job_submission_script()
        #self.run_fits()
        #self.save_results()

    def define_flavors(self):
        if self.scheme == '3f':
            self.flavors = ['b+bb', 'c+cc', 'l']
        elif self.scheme == '5f':
            self.flavors = ['bb', 'cc', 'b', 'c', 'l']
        else:
            raise Exception("Flavor scheme '{}' is not allowed.".format(self.scheme))

    def define_bins(self):
        # Here we assume a 1D histogram therefore we take the first axis
        #num, start, stop = ( self.cfg['variables'][self.var]['axes'][0][key] for key in ['bins', 'start', 'stop'] )
        nbins = 12/self.binwidth
        self.bins = np.arange(self.lo, self.hi + self.binwidth, self.binwidth)
        #if (self.lo not in self.bins) | (self.hi not in self.bins):
        #    raise Exception("The low and high bounds of the interval have to be present in the new bin edges.")

    def define_observable(self):
        self.observable = rl.Observable(self.var, self.bins)

    def initialize_models_dict(self):
        self.categories_region = {}
        self.categories_region['pass'] = list( filter( lambda c : 'pass' in c, self.categories ) )
        if len(self.categories_region['pass']) == 0:
            sys.exit("No 'pass' region found in categories' list")
        self.categories_region['fail'] = [c.replace('pass', 'fail') for c in self.categories_region['pass']]
        self.models = {}
        for region in ['pass']:
            categories = self.categories_region[region]

            for cat in categories:
                model_name = cat.replace(region, '')
                self.models[model_name] = rl.Model("sfModel")

    def load_parameters(self):
        self.parameters = fit_parameters
        for model_name in self.models.keys():
            if not model_name in self.parameters.keys():
                self.parameters[model_name] = {}
                for flavor in self.flavors:
                    self.parameters[model_name][flavor] = {'value' : 1., 'lo' : 0.5, 'hi' : 2.}

    def define_independent_parameters(self):
        self.indep_pars = {}
        for model_name in self.models.keys():
            self.indep_pars[model_name] = {}
            for flavor in self.flavors:
                self.indep_pars[model_name][flavor] = rl.IndependentParameter(flavor, **self.parameters[model_name][flavor])
                #self.indep_pars[model_name][flavor] = rl.IndependentParameter(flavor.replace('+', '_'), **self.parameters[model_name][flavor])

    def define_nuisance_parameters(self):
        self.nuisance_shapes = {}
        for model_name in self.models.keys():
            self.nuisance_shapes[model_name] = {}
            self.nuisance_shapes[model_name]["pileup"] = rl.NuisanceParameter("pileup", "shape")
            self.nuisance_shapes[model_name]["JES_Total"] = rl.NuisanceParameter("JES_Total", "shape")
            self.nuisance_shapes[model_name]["JER"] = rl.NuisanceParameter("JER", "shape")
            self.nuisance_shapes[model_name]["psWeight_isr"] = rl.NuisanceParameter("psWeight_isr", "shape")
            self.nuisance_shapes[model_name]["psWeight_fsr"] = rl.NuisanceParameter("psWeight_fsr", "shape")
            self.nuisance_shapes[model_name]["QCDFlvCompos"] = rl.NuisanceParameter("QCDFlvCompos", "shape")
            self.nuisance_shapes[model_name]["frac_bb"] = rl.NuisanceParameter("frac_bb", "shape")
            self.nuisance_shapes[model_name]["frac_cc"] = rl.NuisanceParameter("frac_cc", "shape")
            self.nuisance_shapes[model_name]["frac_l"] = rl.NuisanceParameter("frac_l", "shape")
            self.frac_effect = 0.2
            if self.year != '2018':
                self.nuisance_shapes[model_name]["sf_L1prefiring"] = rl.NuisanceParameter("sf_L1prefiring", "shape")

        self.nuisance_lnN = {}
        self.nuisance_lnN_effect = {}
        for model_name in self.models.keys():
            self.nuisance_lnN[model_name] = {}
            self.nuisance_lnN_effect[model_name] = {}
            self.nuisance_lnN[model_name]["lumi"] = rl.NuisanceParameter("lumi", "lnN")
            self.nuisance_lnN_effect[model_name]["lumi"] = unc_lumi[self.year]

    def define_frozen_parameters(self):
        self.freeze = {}
        for model_name in self.models.keys():
            self.freeze[model_name] = []
            for name, par in self.parameters[model_name].items():
                if par['lo'] == par['hi']:
                    self.freeze[model_name].append(name)
            #print(model_name, ": Freeze parameters :", self.freeze)

    def define_signal_name(self, model_name):
        self.tagger = [tagger for tagger in AK8Taggers if tagger in model_name][0]
        if self.tagger in AK8Taggers_bb:
            self.signal_name[model_name] = ["b+bb" if self.scheme == '3f' else "bb"][0]
        elif self.tagger in AK8Taggers_cc:
            self.signal_name[model_name] = ["c+cc" if self.scheme == '3f' else "cc"][0]
        else:
            raise Exception("There is no known tagger to calibrate in the given category")

    def get_templ(self, histname, sumw2=True):
        h_vals = self.templates[histname][0]
        bins = self.observable.binning

        if not sumw2:
            return (h_vals, bins, self.observable.name)
        else:
            h_sumw2 = self.templates[histname][1]
            return (h_vals, bins, self.observable.name, h_sumw2)

    def build_fit_models(self):
        self.fitdirs = {}

        # Define fit templates for MC and observation
        for region in regions:
            for cat in self.categories_region[region]:
                model_name = cat.replace(region, '')
                self.define_signal_name(model_name)
                channel = rl.Channel("sf{}".format(region))
                for flavor in self.flavors:
                    histname = "{}_{}_{}_MC_{}_{}".format(self.var, self.year, cat, flavor, "nominal")
                    template = self.get_templ(histname)

                    is_signal = True if flavor == self.signal_name[model_name] else False
                    sType = rl.Sample.SIGNAL if is_signal else rl.Sample.BACKGROUND
                    sample = rl.TemplateSample("{}_{}".format(channel.name, flavor), sType, template)
                    sample.autoMCStats(uncertainty_type='poisson')
                    channel.addSample(sample)

                template_obs = self.get_templ("{}_{}_{}_DATA".format(self.var, self.year, cat), sumw2=False)
                channel.setObservation(template_obs)

                self.models[model_name].addChannel(channel)

        # Define effect of independent parameters on the templates
        for model_name, model in self.models.items():
            for flavor, SF in self.indep_pars[model_name].items():
                sample_pass = self.models[model_name]['sfpass'][flavor]
                sample_fail = self.models[model_name]['sffail'][flavor]
                pass_fail   = sample_pass.getExpectation(nominal=True).sum() / sample_fail.getExpectation(nominal=True).sum()
                sample_pass.setParamEffect(SF, 1.0 * SF)
                sample_fail.setParamEffect(SF, (1 - SF) * pass_fail + 1)

        # Define effect of nuisance parameters on the templates
        # N.B.: the nuisance parameters are fully correlated between flavors and in pass/fail regions
        for region in regions:
            for cat in self.categories_region[region]:
                model_name = cat.replace(region, '')
                for nuisance_name, nuisance in self.nuisance_shapes[model_name].items():
                    if nuisance_name == "frac_bb":
                        sample = self.models[model_name]["sf{}".format(region)][["b+bb" if self.scheme == '3f' else "bb"][0]]
                        sample.setParamEffect(nuisance, 1+self.frac_effect, 1-self.frac_effect)
                        continue
                    elif nuisance_name == "frac_cc":
                        sample = self.models[model_name]["sf{}".format(region)][["c+cc" if self.scheme == '3f' else "cc"][0]]
                        sample.setParamEffect(nuisance, 1+self.frac_effect, 1-self.frac_effect)
                        continue
                    elif nuisance_name == "frac_l":
                        sample = self.models[model_name]["sf{}".format(region)]["l"]
                        sample.setParamEffect(nuisance, 1+self.frac_effect, 1-self.frac_effect)
                        continue
                    else:
                        for flavor in self.flavors:
                            sample = self.models[model_name]["sf{}".format(region)][flavor]
                            histname_nominal = "{}_{}_{}_MC_{}_{}".format(self.var, self.year, cat, flavor, "nominal")
                            histname_up   = "{}_{}_{}_MC_{}_{}".format(self.var, self.year, cat, flavor, nuisance.name+"Up")
                            histname_down = "{}_{}_{}_MC_{}_{}".format(self.var, self.year, cat, flavor, nuisance.name+"Down")
                            h_nominal, _bins, _obs_name = self.get_templ(histname_nominal, sumw2=False)
                            h_up, _bins, _obs_name = self.get_templ(histname_up, sumw2=False)
                            h_down, _bins, _obs_name = self.get_templ(histname_down, sumw2=False)
                            r_up = h_up / h_nominal
                            r_down = h_down / h_nominal
                            effect_up = np.where(~(np.isnan(r_up) | np.isinf(r_up)), r_up, 1)
                            effect_down = np.where(~(np.isnan(r_down) | np.isinf(r_down)), r_down, 1)
                            sample.setParamEffect(nuisance, effect_up, effect_down)

                for nuisance_name, nuisance in self.nuisance_lnN[model_name].items():
                    for flavor in self.flavors:
                        sample = self.models[model_name]["sf{}".format(region)][flavor]
                        sample.setParamEffect(nuisance, self.nuisance_lnN_effect[model_name][nuisance_name])
        # Save the combine fit models to output folder
        for model_name, model in self.models.items():
            fitdir = os.path.abspath(os.path.join(self.output, 'fitdir', model_name))
            if not os.path.exists(fitdir):
                os.makedirs(fitdir)
            model.renderCombine(fitdir)
            self.fitdirs[model_name] = fitdir

    def build_combine_script(self):
        for model_name, model in self.models.items():
            script_FitDiagnostics = os.path.join(self.fitdirs[model_name], 'build.sh')
            script_MultiDimFit = os.path.join(self.fitdirs[model_name], 'build_MultiDimFit.sh')
            os.system("cp {} {}".format(script_FitDiagnostics, script_MultiDimFit))
            with open(script_FitDiagnostics, 'a') as file:
                extra_args = ""
                #extra_args = "--robustHesse=1 "
                extra_args += "--stepSize=0.001 --X-rtd=MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=9999999 --cminFallbackAlgo Minuit2,Migrad,0:0.2 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND"
                #POIs = self.signal_name[model_name].replace('+', '_')
                POIs = ','.join([self.signal_name[model_name]] + [f for f in self.flavors if not f == self.signal_name[model_name]])
                combineCommand = '\ncombine -M FitDiagnostics -d model_combined.root --saveWorkspace --name _{} --cminDefaultMinimizerStrategy 0 --robustFit=1 --saveShapes --saveWithUncertainties --saveOverallShapes --redefineSignalPOIs={} --setParameters r=1 --freezeParameters r --rMin 1 --rMax 1 {}'.format(model_name, POIs, extra_args)
                combineCommand_MultiDimFit = '\ncombine -M MultiDimFit -d model_combined.root --saveWorkspace --name _{} --cminDefaultMinimizerStrategy 0 --robustFit=1 --redefineSignalPOIs={} --setParameters r=1 --freezeParameters r --rMin 1 --rMax 1 {}'.format(model_name, POIs, extra_args)
                setParameters = 'r=1'
                freezeParameters = 'r'
                for par in self.freeze[model_name]:
                    setParameters += ',{}=1'.format(par)
                    freezeParameters += ',{}'.format(par)
                combineCommand = combineCommand.replace('--setParameters r=1', '--setParameters {}'.format(setParameters))
                combineCommand = combineCommand.replace('--freezeParameters r', '--freezeParameters {}'.format(freezeParameters))
                combineCommand_MultiDimFit = combineCommand_MultiDimFit.replace('--setParameters r=1', '--setParameters {}'.format(setParameters))
                combineCommand_MultiDimFit = combineCommand_MultiDimFit.replace('--freezeParameters r', '--freezeParameters {}'.format(freezeParameters))
                file.write(combineCommand)
            with open(script_MultiDimFit, 'a') as file:
                file.write(combineCommand_MultiDimFit)

    def build_job_submission_script(self):
        for model_name, model in self.models.items():
            script_job = os.path.join(self.fitdirs[model_name], 'job.sub')
            firstlines = ['#!/bin/bash\n', '#\n', '#SBATCH -p short\n', '#SBATCH --account=t3\n',
                          '#SBATCH --job-name=fit_mutag\n', '#SBATCH --mem=3000M\n', '#SBATCH --time 00:30:00\n', '#SBATCH -o %x-%j.out\n', '#SBATCH -e %x-%j.err\n', '\n',
                          'echo HOME: $HOME\n', 'echo USER: $USER\n', 'echo SLURM_JOB_ID: $SLURM_JOB_ID\n', 'echo HOSTNAME: $HOSTNAME\n', '\n',
                          'mkdir -p /scratch/$USER/${SLURM_JOB_ID}\n', 'export TMPDIR=/scratch/$USER/${SLURM_JOB_ID}\n', '\n']
            lastlines = ['rm  -rf /scratch/$USER/${SLURM_JOB_ID}\n', '\n', 'date\n']
            with open(script_job, 'a') as file:
                file.writelines(firstlines)
                file.writelines(["cd /work/mmarcheg/CMSSW_10_2_13\n", "cmsenv\n", "cd {}\n".format(self.fitdirs[model_name]), "bash build.sh\n"])
                file.writelines(lastlines)
    def run_fits(self, mode, job=True):
        if mode == "FitDiagnostics":
            command = 'bash build.sh'
            if job == True:
                command = 'sbatch job.sub'
        elif mode == "MultiDimFit":
            command = 'bash build_{}.sh'.format(mode)
            if job == True:
                raise NotImplementedError

        parent_dir = os.getcwd()
        if self.scheme == '3f':
            self.first = {'b+bb' : True, 'c+cc' : True}
        elif self.scheme == '5f':
            self.first = {'bb' : True, 'cc' : True}
        for model_name, fitdir in self.fitdirs.items():
            if fitdir:
                os.chdir(fitdir)
            else:
                sys.exit("The fit directory is not well defined or does not exist")
            print("Running fit of model '{}'".format(model_name))
            print("parameters:", self.parameters[model_name])
            os.system(command)
        os.chdir(parent_dir)

    def save_scans(self, model_name, fitdir):
        combineFile = "higgsCombine_{}.MultiDimFit.mH120.root".format(model_name)
        combineCommand = "combine {} -M MultiDimFit -n .scan.total --algo grid\
                          --snapshotName MultiDimFit --setParameterRanges r=0,2".format(combineFile)
        print("Scanning all nuisances...")
        os.system(combineCommand)
        nuisances = list(self.nuisance_shapes[model_name].keys()) + list(self.nuisance_lnN[model_name].keys())
        nuisances_to_freeze = []
        for nuisance_name in nuisances:
            nuisances_to_freeze.append(nuisance_name)
            combineCommand = "combine {} -M MultiDimFit -n .freeze.{} --algo grid\
                              --snapshotName MultiDimFit --setParameterRanges r=0,2\
                              --freezeParameters {}".format(combineFile, nuisance_name, ','.join(nuisances_to_freeze))
            print("Freezing {}...".format(','.join(nuisances_to_freeze)))
            os.system(combineCommand)
        files = ["higgsCombine.freeze.{}.MultiDimFit.mH120.root".format(nuisance_name) for nuisance_name in nuisances]
        others = ' '.join(["{}:{}:{}".format(files[i], nuisances[i], i+2) for i in range(len(files))])
        scriptName = '/work/mmarcheg/BTVNanoCommissioning/scripts/fit/plot1DScanWithOutput.py'
        combineCommand = 'python {} higgsCombine.scan.total.MultiDimFit.mH120.root --main-label "Total Uncert."\
                          --others {} --output breakdown --y-max 10 --y-cut 40 --breakdown "{},stat" --POI {}'.format(scriptName, others, ','.join(nuisances), self.signal_name[model_name].replace('+', '_'))
        print(combineCommand)
        os.system(combineCommand)        
        #nuisances_to_freeze = nuisances
        #combineCommand = "combine {} -M MultiDimFit -n scan.freeze_all --algo grid --snapshotName MultiDimFit --setParameterRanges r=0,2 --freezeParameters {}".format(combineFile, ','.join(nuisances_to_freeze))
        #print(combineCommand)
        #print("Freezing all nuisances...".format(nuisance_name))

    def save_results(self, mode):
        for model_name, fitdir in self.fitdirs.items():
            wp = model_name.split('wp')[0][-1]
            wpt = model_name.split('Pt-')[1]
            filename = os.path.join(fitdir, "higgsCombine_{}.{}.mH120.root".format(model_name, mode))

            # Wait until the file has been fully written by the job before reading it
            while iserror(uproot.open, filename):
                pass
            combineFile = uproot.open(filename)

            combineTree = combineFile['limit']
            combineBranches = combineTree.arrays()
            results = combineBranches['limit']

            if len(results) < 4:
                print("FIT FAILED : ", model_name)
                return

            combineCont, low, high, temp = results
            combineErrUp = high - combineCont
            combineErrDown = combineCont - low
            d = {}

            POI = self.signal_name[model_name]
            columns = ['year', 'selection', 'wp', 'pt',
                POI, '{}ErrUp'.format(POI), '{}ErrDown'.format(POI),
                'SF({})'.format(POI)]
            columns_for_latex = ['year', 'pt', 'SF({})'.format(POI)]
            d = {'year' : [self.year], 'selection' : [model_name], 'tagger' : [self.tagger],
                'wp' : [wp], 'pt' : [wpt],
                POI : [combineCont], '{}ErrUp'.format(POI) : [combineErrUp], '{}ErrDown'.format(POI) : [combineErrDown],
                'SF({})'.format(POI) : ['{}$^{{+{}}}_{{-{}}}$'.format(combineCont, combineErrUp, combineErrDown)]}

            value, lo, hi = (self.parameters[model_name][POI]['value'], self.parameters[model_name][POI]['lo'], self.parameters[model_name][POI]['hi'])
            f = open(os.path.join(fitdir, "fitResults_{}Pt.txt".format(model_name)), 'w')
            lineIntro = 'Best fit '
            firstline = '{}{}: {}  -{}/+{}  (68%  CL)  range = [{}, {}]\n'.format(lineIntro, POI, combineCont, combineErrDown, combineErrUp, lo, hi)
            f.write(firstline)
            fitResults = ROOT.TFile.Open(os.path.join(fitdir, "fitDiagnostics_{}.root".format(model_name)))
            fit_s = fitResults.Get('fit_s')

            for flavor in self.flavors:
                if (flavor == POI): continue
                par_result = fit_s.floatParsFinal().find(flavor)
                columns.append(flavor)
                columns.append('{}Err'.format(flavor))
                columns.append('SF({})'.format(flavor))
                columns_for_latex.append('SF({})'.format(flavor))
                if par_result == None:
                    d.update({flavor : -999, '{}Err'.format(flavor) : -999, 'SF({})'.format(flavor) : r'{}$\pm${}'.format(-999, -999)})
                    continue
                parVal = par_result.getVal()
                parErr = par_result.getAsymErrorHi()
                value, lo, hi = (self.parameters[model_name][flavor]['value'], self.parameters[model_name][flavor]['lo'], self.parameters[model_name][flavor]['hi'])
                gapSpace = ''.join( (len(lineIntro) + len(POI) - len(flavor) )*[' '])
                lineResult = '{}{}: {}  -+{}'.format(gapSpace, flavor, parVal, parErr)
                gapSpace2 = ''.join( (firstline.find('(') - len(lineResult) )*[' '])
                line = lineResult + gapSpace2 + '(68%  CL)  range = [{}, {}]\n'.format(lo, hi)
                f.write(line)
                #columns.append(flavor)
                #columns.append('{}Err'.format(flavor))
                #columns.append('SF({})'.format(flavor))
                #columns_for_latex.append('SF({})'.format(flavor))
                d.update({flavor : parVal, '{}Err'.format(flavor) : parErr, 'SF({})'.format(flavor) : r'{}$\pm${}'.format(parVal, parErr)})
            f.close()
            df = pd.DataFrame(data=d)
            filename = os.path.join(self.fitdirs[model_name], "fitResults.csv")
            print("Saving fit output in {}".format(filename))
            df.to_csv(filename, columns=columns, mode='w', header=True)
