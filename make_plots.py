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
    datasets = [str(s) for s in h.axis('dataset').identifiers() if str(s) is not 'dataset']
    mapping = {
        'QCD_MuEnrichedPt5' : [dataset for dataset in datasets if 'QCD_Pt' in dataset],
        args.data : [args.data],
    }
    for dataset in datasets:
        if 'QCD' in dataset: continue
        if args.data in dataset: continue
        mapping[dataset] = [dataset]
    datasets = mapping.keys()
    datasets_mc = [dataset for dataset in datasets if args.data not in dataset]

    h = h.group("dataset", hist.Cat("dataset", "Dataset"), mapping)
    #flavors = ['bb', 'cc', 'b', 'c', 'light']
    flavors = ['bb', 'cc', 'b', 'c', 'light', 'others']
    if hist1d:
        fig, (ax, rax) = plt.subplots(2, 1, figsize=(12,12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=.07)
        plot.plot1d(h[datasets_mc].sum('flavor'), ax=ax, legend_opts={'loc':1}, density=args.dense, stack=True)
        plot.plot1d(h[args.data].sum('flavor'), ax=ax, legend_opts={'loc':1}, density=args.dense, error_opts=data_err_opts, clear=False)
        plot.plotratio(num=h[args.data].sum('dataset', 'flavor'), denom=h[datasets_mc].sum('dataset', 'flavor'), ax=rax,
                       error_opts=data_err_opts, denom_fill_opts={}, guide_opts={}, unc='num')
        ax.set_yscale(args.scale)
        rax.set_ylabel('data/MC')
        rax.set_yscale(args.scale)
        rax.set_ylim(0.1,3)
        at = AnchoredText(r"$\geq$1 AK8 jets"+"\n"+
                          r"$p_T > 250 GeV$"+"\n"+
                          r"$m_{SD} > 20 GeV$"+"\n"+
                          r"$\geq$2 $\mu$-tagged AK4 subjets"+"\n",
                              loc=2, frameon=False)
        ax.add_artist(at)
        if histname.startswith("btag"):
            ax.semilogy()
        if (not args.dense) & (args.scale == "log"):
            ax.set_ylim(0.00001, None)
        #hep.mpl_magic(ax)
        filepath = plot_dir + histname + ".png"
        if args.scale != parser.get_default('scale'):
            rax.set_ylim(0.1,10)
            filepath = filepath.replace(".png", "_" + args.scale + ".png")
        print("Saving", filepath)
        plt.savefig(filepath, dpi=300, format="png")
        plt.close(fig)

        #if 'fatjet' in histname:
        fig, (ax, rax) = plt.subplots(2, 1, figsize=(12,12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=.07)
        #iterables = list(itertools.product(datasets_mc,flavors))   # Shortcut to avoid conisdering the flavor 'inclusive'
        #print(h[('QCD_MuEnrichedPt5', ['b', 'c', 'cc'])].values())
        #print(h[[(dataset, flavor) for dataset in datasets_mc for flavor in flavors]].values())
        #plot.plot1d(h[[(dataset, flavor) for dataset in datasets_mc for flavor in flavors]].sum('dataset'), ax=ax, legend_opts={'loc':1}, density=args.dense, stack=True)
        plot.plot1d(h[(datasets_mc, flavors)].sum('dataset'), ax=ax, legend_opts={'loc':1}, density=args.dense, order=flavors, stack=True)
        plot.plot1d(h[args.data].sum('dataset'), ax=ax, legend_opts={'loc':1}, density=args.dense, error_opts=data_err_opts, clear=False)
        plot.plotratio(num=h[args.data].sum('dataset', 'flavor'), denom=h[(datasets_mc, flavors)].sum('dataset', 'flavor'), ax=rax,
                       error_opts=data_err_opts, denom_fill_opts={}, guide_opts={}, unc='num')
        ax.set_yscale(args.scale)
        #print(list(h[args.data].sum('dataset', 'flavor').values().values()))
        #ax.set_ylim(0.001,1.05*max(np.array(h[args.data].sum('dataset', 'flavor').values().values())))
        rax.set_ylabel('data/MC')
        rax.set_yscale(args.scale)
        rax.set_ylim(0.1,3)
        at = AnchoredText(r"$\geq$1 AK8 jets"+"\n"+
                          r"$p_T > 250 GeV$"+"\n"+
                          r"$m_{SD} > 20 GeV$"+"\n"+
                          r"$\geq$2 $\mu$-tagged AK4 subjets"+"\n",
                              loc=2, frameon=False)
        ax.add_artist(at)
        if histname.startswith("btag"):
            ax.semilogy()
        if (not args.dense) & (args.scale == "log"):
            ax.set_ylim(0.0000001, 1000000)
        #hep.mpl_magic(ax)
        filepath = plot_dir + histname + ".png"
        filepath = filepath.replace(".png", "_flavormatch.png")
        if args.scale != parser.get_default('scale'):
            rax.set_ylim(0.1,10)
            filepath = filepath.replace(".png", "_" + args.scale + ".png")
        print("Saving", filepath)
        plt.savefig(filepath, dpi=300, format="png")
        plt.close(fig)
    else:
        for dataset in datasets:
            if 'QCD' in dataset:
                histo_QCD = h[dataset].sum('dataset')
            if 'GluGluHToBB' in dataset:
                histo_BB = h[dataset].sum('dataset')
            if 'GluGluHToCC' in dataset:
                histo_CC = h[dataset].sum('dataset')
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
