from parameters import PtBinning, AK8TaggerWP

"""
parameters_2016_PostVFP =  {}
for tagger in ["deepTagMD_ZHccvsQCD"]:
    for wp in ['L', 'M', 'H']:
        for wpt in ['450to500', '500to600', '600toInf']:
            model_name = "msd40{}{}wp_Pt-{}".format(tagger, wp, wpt)
            parameters_2016_PostVFP[model_name] = {
                "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
                "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
                "c_cc" : {"value" : 1, "lo" : 0, "hi" : 2}
            },
"""

fit_parameters ={
    "2016_PreVFP" : {},
    "2016_PostVFP" : {},
    "2017" : {
        "msd40deepTagMD_ZHccvsQCDLwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDLwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDLwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDMwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDMwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDMwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDHwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDHwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDHwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        # tau21p30_all
        "msd40btagDDCvLV2Lwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : -2, "hi" : 2}
        },
        "msd40btagDDCvLV2Mwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 1.5}
        },
        # tau21p25_all
        "msd40btagDDCvLV2Mwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 1.5}
        },
        # tau21p30_all
        "msd40btagDDCvLV2Hwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : -10, "hi" : 10}
        },
        # tau21p25_all
        "msd40btagDDCvLV2Hwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : -10, "hi" : 10}
        },
        # tau21p35_all
        "msd40particleNetMD_Xcc_QCDHwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : -20, "hi" : 20}
        },
        # tau21p20_all
        "msd40deepTagMD_ZHccvsQCDLwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : -2, "hi" : 2}
        },
        # tau21p20_all
        "msd40deepTagMD_ZHccvsQCDHwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0, "hi" : 2}
        },
        # tau21p30_all
        "msd40btagDDBvLV2Lwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : -10, "hi" : 10},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
        # tau21p25_all
        "msd40btagDDBvLV2Lwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : -20, "hi" : 20},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
        # tau21p25_all
        "msd40btagDDBvLV2Mwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
        # tau21p35_all
        "msd40particleNetMD_Xbb_QCDMwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0, "hi" : 10},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
        # tau21p30_all
        "msd40deepTagMD_ZHbbvsQCDLwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : -10, "hi" : 10},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
        # tau21p35_all
        "msd40deepTagMD_ZHbbvsQCDHwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : -2, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
    },
    "2018" : {
        # cc fits
        "msd40btagDDCvLV2Mwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40btagDDCvLV2Mwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40btagDDCvLV2Mwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40btagDDCvLV2Hwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40btagDDCvLV2Hwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40btagDDCvLV2Hwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40particleNetMD_Xcc_QCDHwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDLwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDLwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDLwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDMwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDMwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDMwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDHwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDHwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        "msd40deepTagMD_ZHccvsQCDHwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.2, "hi" : 2}
        },
        # bb fits
        "msd40btagDDBvLV2Hwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.2, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
        "msd40particleNetMD_Xbb_QCDHwp_Pt-450to500" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.2, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
        "msd40deepTagMD_ZHbbvsQCDMwp_Pt-600toInf" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.2, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
        "msd40btagHbbMwp_Pt-500to600" : {
            "l"  : {"value" : 1, "lo" : 0.5, "hi" : 2},
            "b_bb" : {"value" : 1, "lo" : 0.2, "hi" : 2},
            "c_cc" : {"value" : 1, "lo" : 0.5, "hi" : 2}
        },
        
    }
}
