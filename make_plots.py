import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import mplhep as hep
from coffea.util import load
from coffea.hist import plot
import coffea.hist as hist
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
    if histname in ["sumw", "nbtagmu", "nbtagmu_event_level"]: continue
    if 'ccfatjet' in histname: continue
    print("Plotting", histname)
    #fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(12,12), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)
    if any([histname.startswith('cutflow')]): break
    h = accumulator[histname]
    datasets = [str(s) for s in h.axis('dataset').identifiers() if str(s) is not 'dataset']
    mapping = {
        'QCD_MuEnrichedPt5' : [dataset for dataset in datasets if 'QCD' in dataset],
        args.data : [args.data],
    }
    for dataset in datasets:
        if 'QCD' in dataset: continue
        if args.data in dataset: continue
        mapping[dataset] = [dataset]
    datasets = mapping.keys()
    h = h.group("dataset", hist.Cat("dataset", "Dataset"), mapping)
    plot.plot1d(h[[dataset for dataset in datasets if args.data not in dataset]], ax=ax, legend_opts={'loc':1}, density=args.dense, stack=True)
    plot.plot1d(h[args.data], ax=ax, legend_opts={'loc':1}, density=args.dense, error_opts=data_err_opts, clear=False)
    plot.plotratio(num=h[args.data].sum('dataset'), denom=h[[dataset for dataset in datasets if args.data not in dataset]].sum('dataset'), ax=rax,
                   error_opts=data_err_opts, denom_fill_opts={}, guide_opts={} )#, unc='num')
    ax.set_yscale(args.scale)
    rax.set_ylabel('data/MC')
    rax.set_yscale(args.scale)
    rax.set_ylim(0.000001,2)
    at = AnchoredText(r"$\geq$1 AK8 jets"+"\n"+
                          r"$p_T > 250 GeV$"+"\n"+
                          r"$m_{SD} > 20 GeV$",
                          loc=2, frameon=False)
    ax.add_artist(at)
    if histname.startswith("btag"):
        ax.semilogy()
    if not args.dense:
        ax.set_ylim(0.001, None)
    else:
        plt.set_ylabel('a.u.')
        ax.set_ylim(0.000001, None)
    #hep.mpl_magic(ax)
    filepath = plot_dir + histname + ".png"
    if args.scale != parser.get_default('scale'):
    	filepath.replace(".png", "_" + args.scale + ".png")
    plt.savefig(plot_dir + histname + ".png", dpi=300, format="png")
