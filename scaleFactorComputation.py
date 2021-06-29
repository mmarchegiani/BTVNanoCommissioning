from __future__ import print_function, division
import sys
import os
import rhalphalib as rl
import numpy as np
import scipy.stats
import pickle
import uproot

def exec_me(command, dryRun=False, folder=False):
    print(command)
    if not dryRun:
        if folder: os.chdir(folder)
        os.system(command)


def get_templ(f, sample, obs, syst=None, sumw2=True):
    hist_name = sample
    if syst is not None:
        hist_name += "_" + syst
    h_vals = f[hist_name][0]
    h_sumw2 = f[hist_name][1]
    if not sumw2:
        return (h_vals, obs.binning, obs.name)
    else:
        return (h_vals, obs.binning, obs.name, h_sumw2)


def test_sfmodel(tmpdir, inputFile, sel, wp, fittype='single', scale=1, smear=0.1):
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')
    jecs = rl.NuisanceParameter('CMS_jecs', 'lnN')
    pu = rl.NuisanceParameter('CMS_pu', 'lnN')

    # Indeps
    indep_bb = rl.IndependentParameter('bb', 1., -20, 20)
    indep_cc = rl.IndependentParameter('cc', 1., -20, 20)
    indep_b = rl.IndependentParameter('b', 1., -20, 20)
    indep_c = rl.IndependentParameter('c', 1., -20, 20)
    indep_l = rl.IndependentParameter('l', 1., -20, 20)
    indep_bbf = rl.IndependentParameter('bbf', 1., -20, 20)
    indep_ccf = rl.IndependentParameter('ccf', 1., -20, 20)
    indep_bf = rl.IndependentParameter('bf', 1., -20, 20)
    indep_cf = rl.IndependentParameter('cf', 1., -20, 20)
    indep_lf = rl.IndependentParameter('lf', 1., -20, 20)

    jetprobabins = np.linspace(0, 2.55, 51)
    jetproba = rl.Observable('jetproba', jetprobabins)
    model = rl.Model("sfModel")

    regions = ['pass', 'fail']
    fout = np.load(inputFile, allow_pickle=True)
    sample_names = ['bb', 'cc', 'b', 'c', 'l']
    for region in regions:
        ch = rl.Channel("sf{}".format(region))
        for sName in sample_names:
            template = get_templ(fout, 'fatjet_jetproba_{}{}{}wp_QCD_{}'.format(sel, region, wp, sName), jetproba)

            isSignal = True if sName == ('bb' if sel.endswith('DDB') else 'cc') else False
            sType = rl.Sample.SIGNAL if isSignal else rl.Sample.BACKGROUND
            sample = rl.TemplateSample("{}_{}".format(ch.name, sName), sType, template)
            sample.setParamEffect(lumi, 1.023)
            sample.setParamEffect(jecs, 1.02)
            sample.setParamEffect(pu, 1.05)
            #sample.autoMCStats()
            ch.addSample(sample)

        data_obs = get_templ(fout, 'fatjet_jetproba_{}{}{}wp_BtagMu'.format(sel, region, wp), jetproba)[:-1]
        ch.setObservation(data_obs)

        model.addChannel(ch)

    # for sample, SF in zip(sample_names, [indep_bb, indep_cc, indep_o]):
    #     pass_sample = model['sfpass'][sample]
    #     fail_sample = model['sffail'][sample]
    #     pass_fail = pass_sample.getExpectation(nominal=True).sum() / fail_sample.getExpectation(nominal=True).sum()
    #     pass_sample.setParamEffect(SF, 1.0 * SF)
    #     fail_sample.setParamEffect(SF, (1 - SF) * pass_fail + 1)
    for sample, SF in zip(sample_names, [indep_bb, indep_cc, indep_b, indep_c, indep_l]):
        pass_sample = model['sfpass'][sample]
        pass_sample.setParamEffect(SF, 1.0 * SF)
    for sample, SF in zip(sample_names, [indep_bbf, indep_ccf, indep_bf, indep_cf, indep_lf]):
        fail_sample = model['sffail'][sample]
        fail_sample.setParamEffect(SF, 1.0 * SF)

    model.renderCombine(tmpdir)
    with open(tmpdir+'/build.sh', 'a') as ifile:
        ifile.write('\ncombine -M FitDiagnostics --expectSignal 1 -d sfModel_combined.root --cminDefaultMinimizerStrategy 0 --robustFit=1 --saveShapes  --rMin 0.5 --rMax 1.5')

    exec_me( 'bash build.sh', folder=tmpdir )


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--outputDir', type=str, default=None, help='Output directory')
    parser.add_argument('--year', type=str, choices=['2016', '2017', '2018'], help='Year of data/MC samples')
    parser.add_argument('--selection', type=str, default='msd100tau06DDB', help='Selection to compute SF.')
    parser.add_argument('--wp', type=str, default='M', help='Working point')
    parser.add_argument("--fit", type=str, choices={'single', 'double'}, default='double',
                        help="Fit type")  ##not used
    parser.add_argument("--scale", type=float, default='1',
                        help="Datacard magnitude for scale shift.")  ##not used yet
    parser.add_argument("--smear", type=float, default='0.1',
                        help="Datacard magnitude for smear shift.")  ##not used yet

    parser.add_argument("--tp", "--template-pass", dest='tp', type=str,
                        default='histograms/hists_fattag_pileupJEC_2017_WPcuts_v01.pkl',
                        help="Pass(Pass/Pass) templates")  ##not used

    parser.add_argument("--tpf", "--template-passfail", dest='tpf', type=str,
                        default='histograms/hists_fattag_pileupJEC_2017_WPcuts_v01.pkl',
                        help="Pass/Fail templates, only for `fit=double`")

    parser.add_argument("--tf", "--template-fail", dest='tf', type=str,
                        default='histograms/hists_fattag_pileupJEC_2017_WPcuts_v01.pkl',
                        help="Fail templates")  ##not used

    args = parser.parse_args()
    print("Running with options:")
    print("    ", args)

    output_dir = args.outputDir if args.outputDir else os.getcwd()+"/fitdir/"+args.year+'/'+args.selection+'/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    test_sfmodel(output_dir, args.tpf, args.selection, args.wp)
