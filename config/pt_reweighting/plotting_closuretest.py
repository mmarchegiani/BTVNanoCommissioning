hatch_density = 4

style_cfg = {
    "fontsize" : 22,
    "fontsize_legend_ratio" : 12,
    "opts_figure" : {
        "datamc" : {
            'figsize' : (12,9),
        },
        "datamc_ratio" : {
            'figsize' : (12,12),
            'gridspec_kw' : {"height_ratios": (3, 1)},
            'sharex' : True,
        },
        "partial" : {
            'figsize' : (12,15),
            'gridspec_kw' : {"height_ratios": (3, 1)},
            'sharex' : True,
        },
    },
    "opts_mc" : {
        'histtype': 'fill',
        'stack' : True
    },
    "opts_data" : {
        'linestyle': 'solid',
        'linewidth': 0,
        'marker': '.',
        'markersize': 5.0,
        'color': 'black',
        'elinewidth': 1,
        'label': 'Data',
    },
    "opts_unc" : {
        "total" : {
            "step": "post",
            "color": (0, 0, 0, 0.4),
            "facecolor": (0, 0, 0, 0.0),
            "linewidth": 0,
            "hatch": '/' * hatch_density,
            "zorder": 2,
        },
        'Up': {
            'linestyle': 'dashed',
            'linewidth': 1,
            'marker': '.',
            'markersize': 1.0,
            #'color': 'red',
            'elinewidth': 1,
        },
        'Down': {
            'linestyle': 'dotted',
            'linewidth': 1,
            'marker': '.',
            'markersize': 1.0,
            #'color': 'red',
            'elinewidth': 1,
        },
    },
    "opts_syst" : {
        'nominal': {
            'linestyle': 'solid',
            'linewidth': 1,
            'color': 'black',
        },
        'up': {
            'linestyle': 'dashed',
            'linewidth': 1,
            'color': 'red',
        },
        'down': {
            'linestyle': 'dotted',
            'linewidth': 1,
            'color': 'blue',
        },
    },
    "samples_map": {},
    "samples_exclude": [],
    "labels": {
        "l": "light",
        "c": "c",
        "b": "b",
        "cc": "cc",
        "bb": "bb",
    },
    "colors" : {
        'l': (0.2, 0.60, 1.0),  # blue
        'c': (0.0, 0.8, 0.8), # green/blue
        'b' : (1.0, 0.71, 0.24), # orange
        'cc': (0.4, 0.8, 0.0),  # green
        'bb' : (1.0, 0.4, 0.4), #red
    }
}

#samples = ["QCD_MuEnriched", "SingleTop_ttbar", "VJets", "GluGluHToBB", "GluGluHToCC"]
samples = ["QCD_MuEnriched", "SingleTop_ttbar", "VJets"]
samples_exclude = ["GluGluHToBB", "GluGluHToCC"]
flavors = ["l", "c", "b", "cc", "bb"]

for f in flavors:
    style_cfg["samples_map"][f] = []
    for s in samples:
        style_cfg["samples_map"][f].append(f"{s}_{f}")

for f in flavors:
    for s in samples_exclude:
        style_cfg["samples_exclude"].append(f"{s}_{f}")
