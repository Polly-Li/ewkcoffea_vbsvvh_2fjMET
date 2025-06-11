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

import get_wwz_yields as plt_tools

HTML_PC = "/home/users/kmohrman/ref_scripts/html_stuff/index.php"

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
    "exactly1lep_exactly1fj550",
    "exactly1lep_exactly1fj550_2j",

    "exactly1lep_exactly1fj_2j",

    "exactly1lep_exactly1fj1100",
    "exactly1lep_exactly1fj800_0j1j",
    "exactly1lep_exactly1fj700_0jcent1jcent",
    "exactly1lep_exactly1fj700_0j",

    "exactly1lep_exactly1fj_STmet900",
    "exactly1lep_exactly1fj_STmet1100",
    "exactly1lep_exactly1fj_ST600",

    "exactly1lep_exactly1fj_STmetFjpt1000",
    "exactly1lep_exactly1fj_STmetFjpt1500",
    "exactly1lep_exactly1fj_STmetFjpt1700",
]

def main():

    #pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/tmp_all_njets_njetscounts.pkl.gz"
    #pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/tmp_vvh_plots_0.pkl.gz"
    pkl_file_path = "/home/users/kmohrman/vbs_vvh/ewkcoffea_for_vbs_vvh/ewkcoffea/analysis/vbs_vvh/histos/tmp_vvh_plots_8.pkl.gz"

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(pkl_file_path))

    get_ylds = 1
    make_plots = 0

    if get_ylds:

        var_name = "njets"
        #var_name = "njets_counts"

        for cat in CAT_LST:
            cat_yld = sum(histo_dict[var_name][{"systematic":"nominal", "category":cat}].values(flow=True))
            #cat_yld_raw = sum(histo_dict["njets_counts"][{"systematic":"nominal", "category":cat}].values(flow=True))

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

            yld_bkg = yld_bkg_ttbar
            yld_sig = yld_sig_WWH_OS + yld_sig_WWH_SS + yld_sig_WZH + yld_sig_ZZH
            metric = yld_sig/(yld_bkg**0.5)

            sumw2_bkg = sumw2_bkg_ttbar
            sumw2_sig = sumw2_sig_WWH_OS + sumw2_sig_WWH_SS + sumw2_sig_WZH + sumw2_sig_ZZH

            print("\n",cat)
            #print(yld_sig_WWH_OS)
            #print(yld_sig_WWH_SS)
            #print(yld_sig_WZH)
            #print(yld_sig_ZZH)
            #print(metric)
            #print(yld_sig)
            #print(yld_bkg)
            perr_sig = 100*(((sumw2_sig)**0.5)/yld_sig)
            perr_bkg = 100*(((sumw2_bkg)**0.5)/yld_bkg)
            print(f"{yld_sig} +- {np.round(perr_sig,2)}%")
            print(f"{yld_bkg} +- {np.round(perr_bkg,2)}%")


    if make_plots:

        #cat_lst = ["exactly1lep_exactly1fj", "exactly1lep_exactly1fj550", "exactly1lep_exactly1fj550_2j", "exactly1lep_exactly1fj_2j"]

        for cat in CAT_LST:
            print("\nCat:",cat)
            for var in histo_dict.keys():
                print("\nVar:",var)

                histo = copy.deepcopy(histo_dict[var][{"systematic":"nominal", "category":cat}])

                if var not in ["njets","nleps","nbtagsl","nbtagsm","njets_counts","nleps_counts","nfatjets","njets_forward","njets_tot"]:
                    histo = plt_tools.rebin(histo,6)

                # Regroup
                histo = plt_tools.group(histo,"process","process_grp",GRP_DICT_INDIV)
                #histo = plt_tools.group(histo,"process","process_grp",GRP_DICT_SUM)

                # Handle overflow
                histo = plt_tools.merge_overflow(histo)

                # Scale down bkg
                scale_dict = {"ttbar" : 0.0001}
                for i, name in enumerate(histo.axes["process_grp"]):
                    # Scale the hist, see https://github.com/CoffeaTeam/coffea/discussions/705
                    histo.view(flow=True)[i] *= scale_dict.get(name,1) # Scale by 1 if the process is not ttZ or ZZ

                title = f"{cat}__{var}"
                fig,ext_tup = plt_tools.make_cr_fig(histo,title=title)

                # Save
                save_dir_path = "plots"
                save_dir_path_cat = os.path.join(save_dir_path,cat)
                if not os.path.exists(save_dir_path_cat): os.mkdir(save_dir_path_cat)
                #fig.savefig(os.path.join(save_dir_path_cat,title+".pdf"))
                #fig.savefig(os.path.join(save_dir_path_cat,title+".png"))
                fig.savefig(os.path.join(save_dir_path_cat,title+".png"),bbox_extra_artists=ext_tup,bbox_inches='tight')
                shutil.copyfile(HTML_PC, os.path.join(save_dir_path_cat,"index.php"))


main()







