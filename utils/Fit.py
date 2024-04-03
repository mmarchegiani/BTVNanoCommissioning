import os
import sys
from time import sleep

import numpy as np
import pandas as pd
import rhalphalib as rl
import uproot
import ROOT

from parameters import AK8Taggers, AK8Taggers_bb, AK8Taggers_cc

unc_lumi = {
    '2016_PreVFP'  : 1.012,
    '2016_PostVFP' : 1.012,
    '2017'         : 1.023,
    '2018'         : 1.025
}

epsilon = 1e-3

def iserror(func, *args, **kw):
    try:
        func(*args, **kw)
        return False
    except Exception:
        return True
    
def get_correlation_matrix(fit_s, scheme='3f'):
    parameters = fit_s.floatParsFinal()
    parameter_names = [parameters.at(i).GetName() for i in range(parameters.getSize())]
    if scheme == '3f':
        pois = ['b_bb', 'c_cc', 'l']
    elif scheme == '5f':
        pois = ['bb', 'cc', 'b', 'c', 'l']
    else:
        raise NotImplementedError
    rows = []
    for poi_x in pois:
        column = []
        i = parameter_names.index(poi_x)
        for poi_y in l:
            j = parameter_names.index(poi_y)
            column.append(fit_s.correlation(parameters.at(i), parameters.at(j)))
        rows.append(column)
    C = np.array(rows)
    return C

def get_correlation(fit_s, par1, par2):
    parameters = fit_s.floatParsFinal()
    parameter_names = [parameters.at(i).GetName() for i in range(parameters.getSize())]
    i = parameter_names.index(par1)
    j = parameter_names.index(par2)
    return fit_s.correlation(parameters.at(i), parameters.at(j))

class Fit():
    def __init__(self, input, output, categories, var, year, xlim, binwidth, regions=['pass', 'fail'], scheme='3f', npoi=3, frac_effect=1.2, parameters={}):
        self.input = input
        self.output = output
        self.scheme = scheme
        assert type(npoi) == int, "The number of POIs has to be an integer number!"
        assert npoi in [1, 3], "The number of POIs has to be 1 or 3!"
        self.npoi = npoi
        self.year = year
        self.var = var
        self.lo = xlim[0]
        self.hi = xlim[1]
        self.binwidth = binwidth
        self.define_flavors()
        self.templates = np.load(os.path.abspath(self.input), allow_pickle=True)
        self.categories = categories
        self.regions = regions
        if (len(self.regions) == 1) & (self.regions[0] == 'pass'):
            self.passonly = True
        elif (len(self.regions) == 2):
            self.passonly = False
        else:
            raise NotImplementedError
        self.frac_effect = frac_effect
        self.parameters = parameters[self.year]
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
        self.define_signal_name()
        self.load_parameters()
        self.define_independent_parameters()
        self.define_nuisance_parameters()
        self.define_frozen_parameters()
        self.build_fit_models()
        self.build_combine_script()
        self.build_job_submission_script()
        #self.run_fits()
        #self.save_results()

    def status(self, model_name):
        file_status = os.path.join(self.fitdirs[model_name], 'STATUS')
        i = 0
        while not os.path.exists(file_status):
            i+=1
            continue
        with open(file_status, 'r') as f:
            status = f.readline().strip('\n')
        return status

    def is_done(self, model_name):
        return self.status(model_name) in ["DONE", "FAILED"]

    def all_done(self):
        return all(self.is_done(model_name) for model_name in self.models.keys())

    def define_flavors(self):
        if self.scheme == '3f':
            self.flavors = ['b_bb', 'c_cc', 'l']
        elif self.scheme == '5f':
            self.flavors = ['bb', 'cc', 'b', 'c', 'l']
        else:
            raise Exception("Flavor scheme '{}' is not allowed.".format(self.scheme))

    def define_bins(self):
        # Here we define the bin edges according to the low and high limits of the observable and the bin width
        self.bins = np.round(np.arange(self.lo, self.hi + self.binwidth, self.binwidth), 3)
        bins_old = np.round(np.arange(-2.4, 6.0 + self.binwidth, self.binwidth), 3)
        print(bins_old)
        print(self.bins)
        i_lo = np.where(bins_old == self.lo)[0][0]
        i_hi = np.where(bins_old == self.hi)[0][0]
        # Here we adjust the h_sumw and h_sumw2 according to the bin edges
        for histname, h in self.templates.items():
            h_sumw = h[0]
            h_sumw2 = h[1]
            # Here we slice the histograms according to the bin edges
            h_sumw = h_sumw[i_lo:i_hi]
            h_sumw2 = h_sumw2[i_lo:i_hi]
            self.templates[histname] = [h_sumw, h_sumw2]

    def define_observable(self):
        self.observable = rl.Observable(self.var, self.bins)

    def initialize_models_dict(self):
        self.categories_region = {}
        self.categories_region['pass'] = list( filter( lambda c : 'pass' in c, self.categories ) )
        if len(self.categories_region['pass']) == 0:
            sys.exit("No 'pass' region found in categories' list")
        if not self.passonly:
            self.categories_region['fail'] = [c.replace('pass', 'fail') for c in self.categories_region['pass']]
        self.models = {}
        for region in ['pass']:
            categories = self.categories_region[region]

            for cat in categories:
                model_name = cat.replace(region, '')
                self.models[model_name] = rl.Model("sfModel")

    def load_parameters(self):
        for model_name in self.models.keys():
            if not model_name in self.parameters.keys():
                print(self.parameters.keys())
                raise Exception("The model '{}' is not defined in the parameters dictionary".format(model_name))
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
            #self.nuisance_shapes[model_name]["reweigh"] = rl.NuisanceParameter("reweigh", "shape")
            if self.year != '2018':
                self.nuisance_shapes[model_name]["sf_L1prefiring"] = rl.NuisanceParameter("sf_L1prefiring", "shape")

        self.nuisance_lnN = {}
        self.nuisance_lnN_effect = {}
        for model_name in self.models.keys():
            self.nuisance_lnN[model_name] = {}
            self.nuisance_lnN_effect[model_name] = {}
            self.nuisance_lnN[model_name]["lumi"] = rl.NuisanceParameter("lumi", "lnN")
            self.nuisance_lnN[model_name]["frac_bb"] = rl.NuisanceParameter("frac_bb", "lnN")
            self.nuisance_lnN[model_name]["frac_cc"] = rl.NuisanceParameter("frac_cc", "lnN")
            self.nuisance_lnN[model_name]["frac_l"] = rl.NuisanceParameter("frac_l", "lnN")
            self.nuisance_lnN_effect[model_name]["lumi"] = unc_lumi[self.year]
            self.nuisance_lnN_effect[model_name]["frac_bb"] = self.frac_effect
            self.nuisance_lnN_effect[model_name]["frac_cc"] = self.frac_effect
            self.nuisance_lnN_effect[model_name]["frac_l"] = self.frac_effect

    def define_frozen_parameters(self):
        self.freeze = {}
        if self.npoi == 3:
            for model_name in self.models.keys():
                self.freeze[model_name] = []
                for name, par in self.parameters[model_name].items():
                    if par['lo'] == par['hi']:
                        self.freeze[model_name].append(name)
                #print(model_name, ": Freeze parameters :", self.freeze)
        elif self.npoi == 1:
            for model_name in self.models.keys():
                self.freeze[model_name] = []
                for name, par in self.parameters[model_name].items():
                    if not name == self.signal_name[model_name]:
                        self.freeze[model_name].append(name)

    def define_signal_name(self):
        for model_name in self.models.keys():
            self.tagger = [tagger for tagger in AK8Taggers if tagger in model_name][0]
            if self.tagger in AK8Taggers_bb:
                self.signal_name[model_name] = ["b_bb" if self.scheme == '3f' else "bb"][0]
            elif self.tagger in AK8Taggers_cc:
                self.signal_name[model_name] = ["c_cc" if self.scheme == '3f' else "cc"][0]
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
        for region in self.regions:
            for cat in self.categories_region[region]:
                model_name = cat.replace(region, '')
                channel = rl.Channel("sf{}".format(region))
                for flavor in self.flavors:
                    flavor_key = flavor.replace('_', '+')
                    histname = "{}_{}_{}_MC_{}_{}".format(self.var, self.year, cat, flavor_key, "nominal")
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
                sample_pass.setParamEffect(SF, 1.0 * SF)
                if not self.passonly:
                    sample_fail = self.models[model_name]['sffail'][flavor]
                    pass_fail   = sample_pass.getExpectation(nominal=True).sum() / sample_fail.getExpectation(nominal=True).sum()
                    sample_fail.setParamEffect(SF, (1 - SF) * pass_fail + 1)

        # Define effect of nuisance parameters on the templates
        # N.B.: the nuisance parameters are fully correlated between flavors and in pass/fail regions
        for region in self.regions:
            for cat in self.categories_region[region]:
                model_name = cat.replace(region, '')
                for nuisance_name, nuisance in self.nuisance_shapes[model_name].items():
                    for flavor in self.flavors:
                        sample = self.models[model_name]["sf{}".format(region)][flavor]
                        flavor_key = flavor.replace('_', '+')
                        histname_nominal = "{}_{}_{}_MC_{}_{}".format(self.var, self.year, cat, flavor_key, "nominal")
                        histname_up   = "{}_{}_{}_MC_{}_{}".format(self.var, self.year, cat, flavor_key, nuisance.name+"Up")
                        histname_down = "{}_{}_{}_MC_{}_{}".format(self.var, self.year, cat, flavor_key, nuisance.name+"Down")
                        h_nominal, _bins, _obs_name = self.get_templ(histname_nominal, sumw2=False)
                        h_up, _bins, _obs_name = self.get_templ(histname_up, sumw2=False)
                        h_down, _bins, _obs_name = self.get_templ(histname_down, sumw2=False)
                        r_up = h_up / h_nominal
                        r_down = h_down / h_nominal
                        effect_up = np.where(~(np.isnan(r_up) | np.isinf(r_up)), r_up, 1)
                        effect_down = np.where(~(np.isnan(r_down) | np.isinf(r_down)), r_down, 1)
                        effect_down = np.where(effect_down < 0.0, 0.0, effect_down)
                        sample.setParamEffect(nuisance, effect_up, effect_down)

                for nuisance_name, nuisance in self.nuisance_lnN[model_name].items():
                    for flavor in self.flavors:
                        sample = self.models[model_name]["sf{}".format(region)][flavor]
                        if ( ( nuisance_name == "lumi" ) | 
                             ( (nuisance_name == "frac_bb") & (flavor == "b_bb") ) |
                             ( (nuisance_name == "frac_cc") & (flavor == "c_cc") ) |
                             ( (nuisance_name == "frac_l") & (flavor == "l") ) ):
                            sample.setParamEffect(nuisance, self.nuisance_lnN_effect[model_name][nuisance_name])
                        else:
                            continue

        # Save the combine fit models to output folder
        for model_name, model in self.models.items():
            fitdir = os.path.abspath(os.path.join(self.output, 'fitdir', model_name))
            if not os.path.exists(fitdir):
                os.makedirs(fitdir)
            model.renderCombine(fitdir)
            # Create STATUS flag file
            with open(os.path.join(fitdir, 'STATUS'), 'w') as f:
                f.write("PENDING\n")
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
                POIs = ','.join([self.signal_name[model_name]] + [f.replace('+', '_') for f in self.flavors if not f.replace('+', '_') == self.signal_name[model_name]])
                combineCommand = '\ncombine -M FitDiagnostics -d model_combined.root --saveWorkspace --name _{} --cminDefaultMinimizerStrategy 0 --robustFit=1 --saveShapes --saveWithUncertainties --saveOverallShapes --redefineSignalPOIs={} --setParameters r=1 --freezeParameters r --rMin 1 --rMax 1 {}'.format(model_name, POIs, extra_args)
                # N.B.: In the MultiDimFit we set the POI to be only the signal SF, so that we can profile it and extract the breakdown of uncertainties
                combineCommand_MultiDimFit = '\ncombine -M MultiDimFit -d model_combined.root --saveWorkspace --name _{} --algo=singles --cminDefaultMinimizerStrategy 0 --robustFit=1 --redefineSignalPOIs={} --setParameters r=1 --freezeParameters r --rMin 1 --rMax 1 {}'.format(model_name, POIs, extra_args)
                setParameters = 'r=1'
                freezeParameters = 'r'
                for par in self.freeze[model_name]:
                    setParameters += ',{}=1'.format(par)
                    freezeParameters += ',{}'.format(par)
                combineCommand = combineCommand.replace('--setParameters r=1', '--setParameters {}'.format(setParameters))
                combineCommand = combineCommand.replace('--freezeParameters r', '--freezeParameters {}'.format(freezeParameters))
                combineCommand_MultiDimFit = combineCommand_MultiDimFit.replace('--setParameters r=1', '--setParameters {}'.format(setParameters))
                combineCommand_MultiDimFit = combineCommand_MultiDimFit.replace('--freezeParameters r', '--freezeParameters {}'.format(freezeParameters))
                # Save a STATUS flag file
                file.write(combineCommand)
                # Save a STATUS flag file
                file.write('\necho "DONE" > STATUS\n')
            # Write the commands for the breakdown of uncertainties in the bash script
            parameter_ranges = ":".join("{}=0.5,2".format(poi) for poi in self.flavors)
            combineFile = "higgsCombine_{}.MultiDimFit.mH120.root".format(model_name)
            breakdown_lines = [
                "combine {} -M MultiDimFit -n .scan.total --algo grid --snapshotName MultiDimFit --redefineSignalPOIs={} --setParameters r=1 --freezeParameters r --setParameterRanges {}\n".format(combineFile, self.signal_name[model_name], parameter_ranges),
            ]
            nuisances = list(self.nuisance_shapes[model_name].keys()) + list(self.nuisance_lnN[model_name].keys())
            nuisances_to_freeze = []
            for nuisance_name in nuisances:
                nuisances_to_freeze.append(nuisance_name)
                combineCommand = "combine {} -M MultiDimFit -n .freeze.{} --algo grid --snapshotName MultiDimFit --redefineSignalPOIs={} --setParameters r=1 --freezeParameters r,{} --setParameterRanges {}\n".format(combineFile, nuisance_name, self.signal_name[model_name], ','.join(nuisances_to_freeze), parameter_ranges)
                breakdown_lines.append(combineCommand + '\n')
            files = ["higgsCombine.freeze.{}.MultiDimFit.mH120.root".format(nuisance_name) for nuisance_name in nuisances]
            others = ' '.join(["{}:{}:{}".format(files[i], nuisances[i], i+2) for i in range(len(files))])
            scriptName = '/work/mmarcheg/BTVNanoCommissioning/scripts/fit/plot1DScanWithOutput.py'
            combineCommand = 'python {} higgsCombine.scan.total.MultiDimFit.mH120.root --main-label "Total Uncert." --others {} --output breakdown --y-max 10 --y-cut 40 --breakdown "{},stat" --POI {}'.format(scriptName, others, ','.join(nuisances), self.signal_name[model_name].replace('+', '_'))
            breakdown_lines.append(combineCommand + '\n')
            
            with open(script_MultiDimFit, 'a') as file:
                file.write(combineCommand_MultiDimFit + '\n')
                file.writelines(breakdown_lines)

    def build_job_submission_script(self):
        for model_name, model in self.models.items():
            script_job = os.path.join(self.fitdirs[model_name], 'job.sub')
            script_job_multidimfit = os.path.join(self.fitdirs[model_name], 'job_MultiDimFit.sub')
            firstlines = ['#!/bin/bash\n', '#\n', '#SBATCH -p short\n', '#SBATCH --account=t3\n',
                          '#SBATCH --job-name=fit_mutag\n', '#SBATCH --mem=3000M\n', '#SBATCH --time 00:30:00\n', '#SBATCH -o %x-%j.out\n', '#SBATCH -e %x-%j.err\n', '\n',
                          'echo HOME: $HOME\n', 'echo USER: $USER\n', 'echo SLURM_JOB_ID: $SLURM_JOB_ID\n', 'echo HOSTNAME: $HOSTNAME\n', '\n',
                          'mkdir -p /scratch/$USER/${SLURM_JOB_ID}\n', 'export TMPDIR=/scratch/$USER/${SLURM_JOB_ID}\n', '\n']
            lastlines = ['rm  -rf /scratch/$USER/${SLURM_JOB_ID}\n', '\n', 'date\n']
            with open(script_job, 'a') as file:
                file.writelines(firstlines)
                file.writelines(["cd /work/mmarcheg/CMSSW_10_2_13\n", "cmsenv\n", "cd {}\n".format(self.fitdirs[model_name]), "bash build.sh\n"])
                file.writelines(lastlines)
            
            with open(script_job_multidimfit, 'a') as file:
                file.writelines(firstlines)
                file.writelines(["cd /work/mmarcheg/CMSSW_10_2_13\n", "cmsenv\n", "cd {}\n".format(self.fitdirs[model_name]), "bash build_MultiDimFit.sh\n"])
                file.writelines(lastlines)
    def run_fits(self, mode, job=True):
        if mode == "FitDiagnostics":
            command = 'bash build.sh'
            if job == True:
                command = 'sbatch job.sub'
        elif mode == "MultiDimFit":
            command = 'bash build_MultiDimFit.sh'
            if job == True:
                command = 'sbatch job_MultiDimFit.sub'

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
                          --others {} --output breakdown --y-max 10 --y-cut 40 --breakdown "{},stat" --POI {}'.format(scriptName, others, ','.join(nuisances), self.signal_name[model_name])
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
                print("{}: Waiting 1 second before reading the output file...".format(model_name))
                sleep(1)
            combineFile = uproot.open(filename)

            combineTree = combineFile['limit']
            combineBranches = combineTree.arrays()
            results = combineBranches['limit']

            if len(results) < 4:
                lines = ["FAILED_FIT\n", "{}\n".format(len(results)), "{}\n".format(results)]
                with open(os.path.join(self.fitdirs[model_name], "FAILED_FIT"), 'w') as f:
                    f.writelines(lines)
                # UPDATE STATUS
                with open(os.path.join(fitdir, 'STATUS'), 'w') as f:
                    f.write("FAILED\n")
                print("FIT FAILED : ", model_name)
                continue

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
                if (flavor == POI):
                    for flavor_y in self.flavors:
                        corr = get_correlation(fit_s, flavor, flavor_y)
                        columns = columns + ['corr_{}_{}'.format(flavor, flavor_y)]
                        d.update({'corr_{}_{}'.format(flavor, flavor_y) : corr})
                    continue
                par_result = fit_s.floatParsFinal().find(flavor)
                columns = columns + [flavor, '{}Err'.format(flavor), 'SF({})'.format(flavor)]
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

            for nuisance_name in ['frac_bb', 'frac_cc', 'frac_l']:
                par_result = fit_s.floatParsFinal().find(nuisance_name)
                columns = columns + [nuisance_name, '{}Err'.format(nuisance_name)]
                columns_for_latex.append(nuisance_name)
                if par_result == None:
                    d.update({nuisance_name : -999, '{}Err'.format(nuisance_name) : -999})
                    continue
                parVal = par_result.getVal()
                parErr = par_result.getAsymErrorHi()
                d.update({nuisance_name : parVal, '{}Err'.format(nuisance_name) : parErr})

            f.close()
            df = pd.DataFrame(data=d)
            filename = os.path.join(self.fitdirs[model_name], "fitResults.csv")
            try:
                print("Saving fit output in {}".format(filename))
                df.to_csv(filename, columns=columns, mode='w', header=True)
            except Exception as e:
                print("Error saving fit output in {}".format(filename))
                print(e)
