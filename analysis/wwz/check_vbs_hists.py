import argparse
import pickle
import json
import gzip
import os
import shutil
import numpy as np
import math
import matplotlib.pyplot as plt
import copy
import hist

from topcoffea.scripts.make_html import make_html
from topcoffea.modules import utils
import topcoffea.modules.MakeLatexTable as mlt

import ewkcoffea.modules.yield_tools as yt
import ewkcoffea.modules.sample_groupings as sg
import ewkcoffea.modules.plotting_tools as plt_tools

import get_wwz_yields as gwy

HTML_PC = "/home/users/kmohrman/ref_scripts/html_stuff/index.php"

GRP_DICT_FULL = {

    "Signal" : [
        "VBSWWH_OS_VBSCuts",
        "VBSWWH_SS_VBSCuts",
        "VBSWZH_VBSCuts",
        "VBSZZH_VBSCuts",
    ],

    "QCD" : [
        "QCD_HT1000to1500",
        "QCD_HT100to200",
        "QCD_HT1500to2000",
        "QCD_HT2000toInf",
        "QCD_HT200to300",
        "QCD_HT300to500",
        "QCD_HT500to700",
        "QCD_HT50to100",
        "QCD_HT700to1000",
    ],

    "ttbar" : [
        "TTToHadronic",
        "TTToSemiLeptonic",
    ],

    "single-t" : [
        "ST_t-channel_antitop_4f",
        "ST_t-channel_top_4f",
        "ST_tW_antitop_5f",
        "ST_tW_top_5f",
    ],

    "ttX" : [
        "ttHTobb_M125",
        "ttHToNonbb_M125",

        "TTWJetsToQQ",
        "TTWW",
        "TTWZ",

    ],

    "Vjets" : [
        "ZJetsToQQ_HT-200to400",
        "ZJetsToQQ_HT-400to600",
        "ZJetsToQQ_HT-600to800",
        "ZJetsToQQ_HT-800toInf",

        "WJetsToQQ_HT-200to400",
        "WJetsToQQ_HT-400to600",
        "WJetsToQQ_HT-600to800",
        "WJetsToQQ_HT-800toInf",

        "EWKWminus2Jets_WToQQ_dipoleRecoilOn",
        "EWKWplus2Jets_WToQQ_dipoleRecoilOn",
        "EWKZ2Jets_ZToLL_M-50",
        "EWKZ2Jets_ZToNuNu_M-50",
        "EWKZ2Jets_ZToQQ_dipoleRecoilOn",
    ],

    "VV" : [
        "WWTo1L1Nu2Q",
        "WWTo4Q",
        "WZJJ_EWK_InclusivePolarization",
        "WZTo1L1Nu2Q",
        "WZTo2Q2L",

        "ZZTo2Nu2Q",
        "ZZTo2Q2L",
        "ZZTo4Q",
    ],

    "VH" : [
        "ZH_HToBB_ZToQQ_M-125",
        "WminusH_HToBB_WToLNu_M-125",
        "WplusH_HToBB_WToLNu_M-125",
    ],

    "VVV" : [
        "WWW_4F",
        "WWZ_4F",
        "WZZ",
        "ZZZ",
        "VHToNonbb_M125",
    ],
}


GRP_DICT_INDIV = {
    "WWH_OS" : ["UL18_VBSWWH_OS"],
    "WWH_SS" : ["UL18_VBSWWH_SS"],
    "WZH" : ["UL18_VBSWZH"],
    "ZZH" : ["UL18_VBSZZH"],
    "ttbar" : ["UL18_TTToSemiLeptonic"],
}

GRP_DICT_SUM = {
    "sig" : [
        "UL18_VBSWWH_OS",
        "UL18_VBSWWH_SS",
        "UL18_VBSWZH",
        "UL18_VBSZZH",
    ],
    "ttbar" : ["UL18_TTToSemiLeptonic"],
}

#cat_lst = ["all_events", "filters", "exactly1lep"]
CAT_LST = [
    "all_events",
    "filters",
    "exactly1lep",
    "exactly1lep_exactly1fj",
    #"exactly1lep_exactly1fj550",
    #"exactly1lep_exactly1fj550_2j",

    #"exactly1lep_exactly1fj_2j",

    #"exactly1lep_exactly1fj1100",
    #"exactly1lep_exactly1fj800_0j1j",
    #"exactly1lep_exactly1fj700_0jcent1jcent",
    #"exactly1lep_exactly1fj700_0j",

    #"exactly1lep_exactly1fj_STmet900",
    "exactly1lep_exactly1fj_STmet1100",
    #"exactly1lep_exactly1fj_ST600",

    #"exactly1lep_exactly1fj_STmetFjpt1000",
    #"exactly1lep_exactly1fj_STmetFjpt1500",
    #"exactly1lep_exactly1fj_STmetFjpt1700",
]


# Append the years to sample names dict
def append_years(sample_dict_base,year_lst):
    out_dict = {}
    for proc_group in sample_dict_base.keys():
        out_dict[proc_group] = []
        for proc_base_name in sample_dict_base[proc_group]:
            for year_str in year_lst:
                out_dict[proc_group].append(f"{year_str}_{proc_base_name}")
    return out_dict


def main():

    get_ylds = 0
    make_plots = 1
    check_rwgt = 0

    # Regular histos
    #pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/tmp_vvh_plots_8.pkl.gz"
    #pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/tmp.pkl.gz"
    #pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/jun12_ref.pkl.gz"
    pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/tmp_full.pkl.gz"

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(pkl_file_path))



    #### Print the yields ####
    if get_ylds:

        #var_name = "njets"
        var_name = "njets_counts"

        tot_raw = sum(sum(histo_dict["njets_counts"][{"systematic":"nominal", "category":"all_events"}].values(flow=True)))
        print("Tot raw events:",tot_raw)
        #exit()

        for cat in CAT_LST:
            cat_yld = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat}].values(flow=True))
            cat_yld_raw = sum(histo_dict["njets_counts"][{"systematic":"nominal", "category":cat}].values(flow=True))

            yld_sig_WWH_OS = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_VBSWWH_OS"}].values(flow=True))
            yld_sig_WWH_SS = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_VBSWWH_SS"}].values(flow=True))
            yld_sig_WZH    = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_VBSWZH"}].values(flow=True))
            yld_sig_ZZH    = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_VBSZZH"}].values(flow=True))
            yld_bkg_ttbar  = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_TTToSemiLeptonic"}].values(flow=True))

            sumw2_sig_WWH_OS = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_VBSWWH_OS"}].variances(flow=True))
            sumw2_sig_WWH_SS = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_VBSWWH_SS"}].variances(flow=True))
            sumw2_sig_WZH    = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_VBSWZH"}].variances(flow=True))
            sumw2_sig_ZZH    = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_VBSZZH"}].variances(flow=True))
            sumw2_bkg_ttbar  = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "process":"UL18_TTToSemiLeptonic"}].variances(flow=True))

            yld_sig = yld_sig_WWH_OS + yld_sig_WWH_SS + yld_sig_WZH + yld_sig_ZZH
            yld_bkg = yld_bkg_ttbar
            #metric = yld_sig/(yld_bkg**0.5)

            sumw2_sig = sumw2_sig_WWH_OS + sumw2_sig_WWH_SS + sumw2_sig_WZH + sumw2_sig_ZZH
            sumw2_bkg = sumw2_bkg_ttbar

            #print(yld_sig_WWH_OS)
            #print(yld_sig_WWH_SS)
            #print(yld_sig_WZH)
            #print(yld_sig_ZZH)
            #print(metric)

            perr_sig = 100*(((sumw2_sig)**0.5)/yld_sig)
            perr_bkg = 100*(((sumw2_bkg)**0.5)/yld_bkg)
            print("\n",cat)
            print(f"{yld_sig} +- {np.round(perr_sig,2)}%")
            print(f"{yld_bkg} +- {np.round(perr_bkg,2)}%")
            #print(f"{cat}: {yld_sig}")
            #print(f"{cat}: {yld_bkg}")



    #### Make the plots ####
    if make_plots:

        #cat_lst = ["exactly1lep_exactly1fj", "exactly1lep_exactly1fj550", "exactly1lep_exactly1fj550_2j", "exactly1lep_exactly1fj_2j"]

        grouping_dict = append_years(GRP_DICT_FULL,["UL16APV","UL16","UL17","UL18"])
        #grouping_dict = GRP_DICT_SUM

        for cat in CAT_LST:
            print("\nCat:",cat)
            for var in histo_dict.keys():
                print("\nVar:",var)
                if var != "njets": continue # TMP

                histo = copy.deepcopy(histo_dict[var][{"systematic":"nominal", "category":cat}])

                # Clean up a bit (rebin, regroup, and handle overflow)
                if var not in ["njets","nleps","nbtagsl","nbtagsm","njets_counts","nleps_counts","nfatjets","njets_forward","njets_tot"]:
                    histo = plt_tools.rebin(histo,6)
                histo = plt_tools.group(histo,"process","process_grp",grouping_dict)
                histo = plt_tools.merge_overflow(histo)

                # Get one hist of just sig and one of just bkg
                grp_names_bkg_lst = list(grouping_dict.keys()) # All names, still need to drop signal
                grp_names_bkg_lst.remove("Signal")
                histo_sig = histo[{"process_grp":["Signal"]}]
                histo_bkg = plt_tools.group(histo,"process_grp","process_grp",{"Background":grp_names_bkg_lst})
                sig_yld = sum(sum(histo_sig.values(flow=True)))
                bkg_yld = sum(sum(histo_bkg.values(flow=True)))
                #scale_dict = {"ttbar" : 0.0001}
                #histo = plt_tools.scale(histo,"process_grp",scale_dict)

                # Make the figure
                title = f"{cat}__{var}"
                fig,ext_tup = gwy.make_vvh_fig(
                    histo_mc = histo,
                    histo_mc_sig = histo_sig,
                    histo_mc_bkg = histo_bkg,
                    title=title
                )

                # Save
                save_dir_path = "plots"
                save_dir_path_cat = os.path.join(save_dir_path,cat)
                if not os.path.exists(save_dir_path_cat): os.mkdir(save_dir_path_cat)
                fig.savefig(os.path.join(save_dir_path_cat,title+".png"),bbox_extra_artists=ext_tup,bbox_inches='tight')
                shutil.copyfile(HTML_PC, os.path.join(save_dir_path_cat,"index.php"))


    ##### Sanity check of the different reweight points (for a hist that has extra axis to store that) ####
    if check_rwgt:

        #pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/check_wgt_genw.pkl.gz"
        #pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/check_wgt_sm.pkl.gz"
        #pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/check_wgt_rwgtscan.pkl.gz"

        var_name = "njets"
        #var_name = "njets_counts"
        cat = "exactly1lep_exactly1fj_STmet1100"

        #cat_yld = sum(sum(histo_dict[var_name][{"systematic":"nominal", "category":cat}].values(flow=True)))
        #cat_err = (sum(sum(histo_dict[var_name][{"systematic":"nominal", "category":cat}].variances(flow=True))))**0.5
        #print(cat_yld, cat_err)
        #exit()

        wgts = []
        for i in range(120):
            idx_name = f"idx{i}"
            cat_yld = sum(sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "rwgtidx":idx_name}].values(flow=True)))
            cat_err = (sum(sum(histo_dict[var_name][{"systematic":"nominal", "category":cat, "rwgtidx":idx_name}].variances(flow=True))))**0.5
            wgts.append(cat_yld)
            print(i,cat_yld)

        print(min(wgts))




main()







