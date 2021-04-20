import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import mplhep as hep
from coffea.util import load
from coffea.hist import plot
import coffea.hist as hist
import itertools
import os
from utils import histogram_settings

parser = argparse.ArgumentParser(description='Plot histograms from coffea file')
parser.add_argument('-i', '--input', type=str, help='Input histogram filename', required=True)
parser.add_argument('-o', '--output', type=str, help='Output directory', required=True)
parser.add_argument('-s', '--scale', type=str, default='linear', help='Plot y-axis scale', required=False)
parser.add_argument('-d', '--dense', action='store_true', help='Normalized plots')
parser.add_argument('--data', type=str, default='BTagMu', help='Data sample name')

args = parser.parse_args()

accumulator = load(args.input)

data_err_opts = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}

qcd_opts = {
    'facecolor': 'yellow',
    'edgecolor': 'black',
    'alpha': 1.0
}

signal_opts = {
    'facecolor': 'None',
    'edgecolor': ['green', 'red'],
    'linestyle': ['--', '-'],
    'linewidth': 2,
    'alpha': 0.7
}

flavor_opts = {
    'facecolor': ['cyan', 'magenta', 'red', 'green', 'blue'],
    'edgecolor': 'black',
    'alpha': 1.0
}

plt.style.use([hep.style.ROOT, {'font.size': 16}])
plot_dir = "plots/" + args.output + "/"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

for histname in accumulator:
    hist1d = not 'hist2d' in histname
    if histname in ["sumw", "nbtagmu", "nbtagmu_event_level", "nfatjet"]: continue
    #if (not 'fatjet' in histname) & (not histname in ['nmusj1', 'nmusj2', 'nsv1', 'nsv2']): continue
    if (not 'fatjet' in histname) & (not 'nmusj' in histname) & (not 'nsv' in histname): continue
    print("Plotting", histname)
    #fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    
    if any([histname.startswith('cutflow')]): break
    h = accumulator[histname]
    if histname in histogram_settings['variables'].keys():
        varname = h.fields[-1]
        varlabel = h.axis(varname).label
        h = h.rebin(varname, hist.Bin(varname, varlabel, **histogram_settings['variables'][histname]['binning']))
    datasets = [str(s) for s in h.axis('dataset').identifiers() if str(s) is not 'dataset']
    mapping = {
        r'QCD ($\mu$ enriched)' : [dataset for dataset in datasets if 'QCD_Pt' in dataset],
        args.data : [args.data],
    }
    for dataset in datasets:
        if 'QCD' in dataset: continue
        if args.data in dataset: continue
        mapping[dataset] = [dataset]
    datasets = mapping.keys()
    #datasets_mc  = [dataset for dataset in datasets if args.data not in dataset]
    datasets_QCD = [dataset for dataset in datasets if ((args.data not in dataset) & ('GluGlu' not in dataset))]
    datasets_ggH = [dataset for dataset in datasets if 'GluGlu' in dataset]

    h = h.group("dataset", hist.Cat("dataset", "Dataset"), mapping)
    flavors = ['bb', 'cc', 'b', 'c', 'light']
    #flavors = ['bb', 'cc', 'b', 'c', 'light', 'others']
    if hist1d:
        fig, (ax, rax) = plt.subplots(2, 1, figsize=(12,12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=.07)
        plot.plot1d(h[datasets_QCD].sum('flavor'), ax=ax, legend_opts={'loc':1}, density=args.dense, fill_opts=qcd_opts, stack=True)
        ggH_rescaled = h[datasets_ggH].sum('flavor')
        scale_ggH = 1000
        ggH_rescaled.scale(scale_ggH)
        ax_ggH = plot.plot1d(ggH_rescaled, ax=ax, legend_opts={'loc':1}, density=args.dense, fill_opts=signal_opts, stack=False, clear=False)
        plot.plot1d(h[args.data].sum('flavor'), ax=ax, legend_opts={'loc':1}, density=args.dense, error_opts=data_err_opts, clear=False)
        plot.plotratio(num=h[args.data].sum('dataset', 'flavor'), denom=h[datasets_QCD].sum('dataset', 'flavor'), ax=rax,
                       error_opts=data_err_opts, denom_fill_opts={}, guide_opts={}, unc='num')
        handles, labels = ax.get_legend_handles_labels()
        for (i, label) in enumerate(labels):
            if "GluGlu" in label:
                if "BB" in label:
                    labels[i] = r"ggH$\rightarrow$bb $\times$" + str(scale_ggH)
                if "CC" in label:
                    labels[i] = r"ggH$\rightarrow$cc $\times$" + str(scale_ggH)
        ax.legend(handles, labels)
        ax.set_yscale(args.scale)
        rax.set_ylabel('Data/MC')
        #.rax.set_yscale(args.scale)
        rax.set_ylim(0.5,1.5)
        if histname in histogram_settings['variables'].keys():
            ax.set_xlim(**histogram_settings['variables'][histname]['xlim'])
            rax.set_xlim(**histogram_settings['variables'][histname]['xlim'])
        at = AnchoredText(r"$\geq$1 AK8 jets"+"\n"+
                          r"$p_T > 250 GeV$"+"\n"+
                          r"$m_{SD} > 20 GeV$"+"\n"+
                          r"$\geq$2 $\mu$-tagged AK4 subjets"+"\n",
                              loc=2, frameon=False)
        ax.add_artist(at)
        if histname.startswith("btag"):
            ax.semilogy()
        if (not args.dense) & (args.scale == "log"):
            ax.set_ylim(0.1, 10**7)
        #hep.mpl_magic(ax)
        filepath = plot_dir + histname + ".png"
        if args.scale != parser.get_default('scale'):
            #rax.set_ylim(0.1,10)
            filepath = filepath.replace(".png", "_" + args.scale + ".png")
        print("Saving", filepath)
        plt.savefig(filepath, dpi=300, format="png")
        plt.close(fig)

        #if 'fatjet' in histname:
        fig, (ax, rax) = plt.subplots(2, 1, figsize=(12,12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=.07)
        #iterables = list(itertools.product(datasets_QCD,flavors))   # Shortcut to avoid conisdering the flavor 'inclusive'
        #print(h[('QCD_MuEnrichedPt5', ['b', 'c', 'cc'])].values())
        #print(h[[(dataset, flavor) for dataset in datasets_QCD for flavor in flavors]].values())
        #plot.plot1d(h[[(dataset, flavor) for dataset in datasets_QCD for flavor in flavors]].sum('dataset'), ax=ax, legend_opts={'loc':1}, density=args.dense, stack=True)
        plot.plot1d(h[(datasets_QCD, flavors)].sum('dataset'), ax=ax, legend_opts={'loc':1}, density=args.dense, fill_opts=flavor_opts, order=flavors, stack=True)
        plot.plot1d(h[args.data].sum('dataset'), ax=ax, legend_opts={'loc':1}, density=args.dense, error_opts=data_err_opts, clear=False)
        plot.plotratio(num=h[args.data].sum('dataset', 'flavor'), denom=h[(datasets_QCD, flavors)].sum('dataset', 'flavor'), ax=rax,
                       error_opts=data_err_opts, denom_fill_opts={}, guide_opts={}, unc='num')
        ax.set_yscale(args.scale)
        #print(list(h[args.data].sum('dataset', 'flavor').values().values()))
        #ax.set_ylim(0.001,1.05*max(np.array(h[args.data].sum('dataset', 'flavor').values().values())))
        rax.set_ylabel('Data/MC')
        #rax.set_yscale(args.scale)
        rax.set_ylim(0.5,1.5)
        if histname in histogram_settings['variables'].keys():
            ax.set_xlim(**histogram_settings['variables'][histname]['xlim'])
            rax.set_xlim(**histogram_settings['variables'][histname]['xlim'])
        at = AnchoredText(r"$\geq$1 AK8 jets"+"\n"+
                          r"$p_T > 250 GeV$"+"\n"+
                          r"$m_{SD} > 20 GeV$"+"\n"+
                          r"$\geq$2 $\mu$-tagged AK4 subjets"+"\n",
                              loc=2, frameon=False)
        ax.add_artist(at)
        if histname.startswith("btag"):
            ax.semilogy()
        if (not args.dense) & (args.scale == "log"):
            ax.set_ylim(0.1, 10**7)
        #hep.mpl_magic(ax)
        filepath = plot_dir + histname + ".png"
        filepath = filepath.replace(".png", "_flavormatch.png")
        if args.scale != parser.get_default('scale'):
            #rax.set_ylim(0.1,10)
            filepath = filepath.replace(".png", "_" + args.scale + ".png")
        print("Saving", filepath)
        plt.savefig(filepath, dpi=300, format="png")
        plt.close(fig)
    else:
        for dataset in datasets:
            if 'QCD' in dataset:
                histo_QCD = h[dataset].sum('dataset')
                histo_QCD_bb = h[dataset].sum('dataset')['bb']
                histo_QCD_cc = h[dataset].sum('dataset')['cc']
            if 'GluGluHToBB' in dataset:
                histo_BB    = h[dataset].sum('dataset')
                histo_BB_bb = h[dataset].sum('dataset')['bb']
            if 'GluGluHToCC' in dataset:
                histo_CC    = h[dataset].sum('dataset')
                histo_CC_cc = h[dataset].sum('dataset')['cc']
        xaxis = [axis.name for axis in histo_QCD.axes() if 'btag' in axis.name][0]

        for histo_GluGlu, dataset_GluGlu in zip([histo_BB, histo_CC], ['GluGluHToBB', 'GluGluHToCC']):
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24,9))
            plot.plot2d(histo_QCD.sum('flavor'), xaxis=xaxis, ax=ax1)
            plot.plot2d(histo_GluGlu.sum('flavor'), xaxis=xaxis, ax=ax2)
            ax1.set_title('QCD')
            ax2.set_title(dataset_GluGlu)
            filepath = plot_dir + histname + "_" + dataset_GluGlu + ".png"
            print("Saving", filepath)
            plt.savefig(filepath, dpi=300, format="png")
            plt.close(fig)

        for histo_GluGlu_xx, dataset_GluGlu_xx in zip([histo_BB_bb, histo_CC_cc], ['GluGluHToBB (bb)', 'GluGluHToCC (cc)']):
            if dataset_GluGlu_xx == 'GluGluHToBB (bb)':
                histo_QCD_xx = histo_QCD_bb
                qcd_label = 'QCD (bb)'
            elif dataset_GluGlu_xx == 'GluGluHToCC (cc)':
                histo_QCD_xx = histo_QCD_cc
                qcd_label = 'QCD (cc)'
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24,9))
            plot.plot2d(histo_QCD_xx.sum('flavor'), xaxis=xaxis, ax=ax1)
            plot.plot2d(histo_GluGlu_xx.sum('flavor'), xaxis=xaxis, ax=ax2)
            ax1.set_title(qcd_label)
            ax2.set_title(dataset_GluGlu_xx)
            filepath = plot_dir + histname + "_" + '_'.join(dataset_GluGlu_xx.strip(')').split(' (')) + ".png"
            print("Saving", filepath)
            plt.savefig(filepath, dpi=300, format="png")
            plt.close(fig)

        for flavor in flavors:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24,9))
            plot.plot2d(h.sum('dataset')['light'].sum('flavor'), xaxis=xaxis, ax=ax1)
            plot.plot2d(h.sum('dataset')[flavor].sum('flavor'), xaxis=xaxis, ax=ax2)
            ax1.set_title('light')
            ax2.set_title(flavor)
            filepath = plot_dir + histname + "_" + flavor + ".png"
            print("Saving", filepath)
            plt.savefig(filepath, dpi=300, format="png")
            plt.close(fig)
