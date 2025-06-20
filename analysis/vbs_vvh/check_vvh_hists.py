import argparse
import pickle
import gzip
import json
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import copy

import ewkcoffea.modules.plotting_tools as plt_tools


HTML_PC = "/home/users/kmohrman/ref_scripts/html_stuff/index.php"

CLR_LST = ['#d55e00', '#e69f00', '#f0e442', '#009e73', '#0072b2', '#56b4e9', '#cc79a7', '#6e3600', '#a17500'] #, '#a39b2f', '#00664f', '#005d87', '#999999', '#8c5d77']

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
        "QCD_HT50to100", # Has a spike
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


CAT_LST = [
    # The ones we store in the ref file
    #"all_events",
    #"filters",
    #"exactly1lep",
    #"exactly1lep_exactly1fj",
    #"exactly1lep_exactly1fj_STmet1000",
    #"exactly1lep_exactly1fj_STmet1000_msd170",
    #"exactly1lep_exactly1fj_STmet1000_msd170_NjCentralLessThan3",

    "exactly1lep_exactly1fj",
    "exactly1lep_exactly1fj_STmet600",

]


########################
### Helper functions ###

# Append the years to sample names dict
def append_years(sample_dict_base,year_lst):
    out_dict = {}
    for proc_group in sample_dict_base.keys():
        out_dict[proc_group] = []
        for proc_base_name in sample_dict_base[proc_group]:
            for year_str in year_lst:
                out_dict[proc_group].append(f"{year_str}_{proc_base_name}")
    return out_dict


# Get sig and bkg yield in all categories
def get_yields_per_cat(histo_dict,var_name):
    out_dict = {}

    # Get the initial grouping dict
    grouping_dict = append_years(GRP_DICT_FULL,["UL16APV","UL16","UL17","UL18"])

    # Get list of all of the backgrounds together
    bkg_lst = []
    for grp in grouping_dict:
        if grp != "Signal":
            bkg_lst = bkg_lst + grouping_dict[grp]

    # Make the dictionary to get yields for, it includes what's in grouping_dict, plus the backgrounds grouped as one
    groups_to_get_yields_for_dict = copy.deepcopy(grouping_dict)
    groups_to_get_yields_for_dict["Background"] = bkg_lst

    # Loop over cats and fill dict of sig and bkg
    for cat in CAT_LST:
        out_dict[cat] = {}
        histo_base = histo_dict[var_name][{"systematic":"nominal", "category":cat}]

        # Get values per proc
        for group_name,group_lst in groups_to_get_yields_for_dict.items():
            histo = plt_tools.group(histo_base,"process","process",{group_name:group_lst})
            yld = sum(sum(histo.values(flow=True)))
            var = sum(sum(histo.variances(flow=True)))
            out_dict[cat][group_name] = [yld,(var)**0.5]

        # Get the metric
        sig = out_dict[cat]["Signal"][0]
        bkg = out_dict[cat]["Background"][0]
        metric = sig/(bkg)**0.5
        out_dict[cat]["metric"] = [metric,None] # Don't bother propagating error

    return out_dict


# Make the figures for the vvh sudy
def make_vvh_fig(histo_mc,histo_mc_sig,histo_mc_bkg,title="test",axisrangex=None):

    # Create the figure
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(
        nrows=4,
        ncols=1,
        figsize=(7,10),
        gridspec_kw={"height_ratios": (3, 1, 1, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.07)

    # Plot the stack plot
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        color=CLR_LST,
        ax=ax1,
        zorder=100,
    )

    # Get the errs on MC and plot them by hand on the stack plot
    histo_mc_sum = histo_mc[{"process_grp":sum}]
    mc_arr = histo_mc_sum.values()
    mc_err_arr = np.sqrt(histo_mc_sum.variances())
    err_p = np.append(mc_arr + mc_err_arr, 0)
    err_m = np.append(mc_arr - mc_err_arr, 0)
    bin_edges_arr = histo_mc_sum.axes[0].edges
    bin_centers_arr = histo_mc_sum.axes[0].centers
    ax1.fill_between(bin_edges_arr,err_m,err_p, step='post', facecolor='none', edgecolor='gray', alpha=0.5, linewidth=0.0, label='MC stat', hatch='/////')


    ## Draw the normalized shapes ##

    # Get normalized hists of sig and bkg
    yld_sig = sum(sum(histo_mc_sig.values(flow=True)))
    yld_bkg = sum(sum(histo_mc_bkg.values(flow=True)))
    metric = yld_sig/(yld_bkg**0.5)
    histo_mc_sig_scale_to_bkg = plt_tools.scale(copy.deepcopy(histo_mc_sig), "process_grp", {"Signal":yld_bkg/yld_sig})
    histo_mc_sig_norm         = plt_tools.scale(copy.deepcopy(histo_mc_sig), "process_grp", {"Signal":1.0/yld_sig})
    histo_mc_bkg_norm         = plt_tools.scale(copy.deepcopy(histo_mc_bkg), "process_grp", {"Background":1.0/yld_bkg})

    histo_mc_sig_scale_to_bkg.plot1d(color=["red"], ax=ax1, zorder=100)
    histo_mc_sig_norm.plot1d(color="red",  ax=ax2, zorder=100)
    histo_mc_bkg_norm.plot1d(color="gray", ax=ax2, zorder=100)


    ## Draw the significance ##

    # Get the sig and bkg arrays (Not including flow bins here, overflow should already be handled, and if we have underflow, why?)
    yld_sig_arr = sum(histo_mc_sig.values())
    yld_bkg_arr = sum(histo_mc_bkg.values())

    # Get the cumulative signifiance, starting from left
    yld_sig_arr_cum = np.cumsum(yld_sig_arr)
    yld_bkg_arr_cum = np.cumsum(yld_bkg_arr)
    metric_cum = yld_sig_arr_cum/np.sqrt(yld_bkg_arr_cum)
    metric_cum = np.nan_to_num(metric_cum,nan=0,posinf=0) # Set the nan (from sig and bkg both being 0) to 0

    # Get the cumulative signifiance, starting from right
    yld_sig_arr_cum_ud = np.cumsum(np.flipud(yld_sig_arr))
    yld_bkg_arr_cum_ud = np.cumsum(np.flipud(yld_bkg_arr))
    metric_cum_ud = np.flipud(yld_sig_arr_cum_ud/np.sqrt(yld_bkg_arr_cum_ud))
    metric_cum_ud = np.nan_to_num(metric_cum_ud,nan=0,posinf=0) # Set the nan (from sig and bkg both being 0) to 0
    yld_sig_arr_cum_ud = np.flipud(yld_sig_arr_cum_ud) # Flip back so the order is as expected for later use
    yld_bkg_arr_cum_ud = np.flipud(yld_bkg_arr_cum_ud) # Flip back so the order is as expected for later use

    # Draw it on the third plot
    ax3.scatter(bin_centers_arr,metric_cum,   facecolor='none',edgecolor='black',marker=">",label="Cum. from left", zorder=100)
    ax3.scatter(bin_centers_arr,metric_cum_ud,facecolor='none',edgecolor='black',marker="<",label="Cum. from right", zorder=100)

    # Write the max values on the plot
    max_metric_from_left_idx  = np.argmax(metric_cum)
    max_metric_from_right_idx = np.argmax(metric_cum_ud)
    left_max_y  = metric_cum[max_metric_from_left_idx]
    right_max_y = metric_cum_ud[max_metric_from_right_idx]
    left_max_x  = bin_centers_arr[max_metric_from_left_idx]
    right_max_x = bin_centers_arr[max_metric_from_right_idx]
    left_s_at_max  = yld_sig_arr_cum[max_metric_from_left_idx]
    right_s_at_max = yld_sig_arr_cum_ud[max_metric_from_right_idx]
    left_b_at_max  = yld_bkg_arr_cum[max_metric_from_left_idx]
    right_b_at_max = yld_bkg_arr_cum_ud[max_metric_from_right_idx]
    plt.text(0.15,0.35, f"Max from left:  {np.round(left_max_y,3)} (at x={np.round(left_max_x,2)}, sig: {np.round(left_s_at_max,2)}, bkg: {np.round(left_b_at_max,1)})", fontsize=9, transform=fig.transFigure)
    plt.text(0.15,0.33, f"Max from right: {np.round(right_max_y,3)} (at x={np.round(right_max_x,2)} , sig: {np.round(right_s_at_max,2)}, bkg: {np.round(right_b_at_max,1)})", fontsize=9, transform=fig.transFigure)


    ## Draw on the fraction of signal retained ##
    yld_sig_arr_cum_frac    = np.cumsum(yld_sig_arr)/yld_sig
    yld_sig_arr_cum_frac_ud = np.flipud(np.cumsum(np.flipud(yld_sig_arr)))/yld_sig
    ax4.scatter(bin_centers_arr,yld_sig_arr_cum_frac,   facecolor='none',edgecolor='black',marker=">",label="Cum. from left", zorder=100)
    ax4.scatter(bin_centers_arr,yld_sig_arr_cum_frac_ud,facecolor='none',edgecolor='black',marker="<",label="Cum. from right", zorder=100)


    ## Legend, scale the axis, set labels, etc ##

    extr = ax1.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize="12", frameon=False)
    extr = ax2.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize="10", frameon=False)
    extr = ax3.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize="10", frameon=False)
    extr = ax4.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize="10", frameon=False)
    plt.text(0.15,0.85, f"Sig. yield: {np.round(yld_sig,2)}", fontsize = 11, transform=fig.transFigure)
    plt.text(0.15,0.83, f"Bkg. yield: {np.round(yld_bkg,2)}", fontsize = 11, transform=fig.transFigure)
    plt.text(0.15,0.81, f"Metric: {np.round(metric,3)}", fontsize = 11, transform=fig.transFigure)
    plt.text(0.15,0.79, f"[Note: sig. overlay scaled {np.round(yld_bkg/yld_sig,1)}x]", fontsize = 12, transform=fig.transFigure)

    extt = ax1.set_title(title)
    ax1.set_xlabel(None)
    ax2.set_xlabel(None)
    extb = ax3.set_xlabel(None)
    # Plot a dummy hist on ax4 to get the label to show up
    histo_mc.plot1d(alpha=0, ax=ax4)

    extl = ax2.set_ylabel('Shapes')
    ax3.set_ylabel('Significance')
    ax4.set_ylabel('Signal kept (%)')
    ax1.tick_params(axis='y', labelsize=16)
    ax2.tick_params(axis='x', labelsize=16)
    ax3.axhline(0.0,linestyle="-",color="k",linewidth=0.5)
    ax4.axhline(0.0,linestyle="-",color="k",linewidth=0.5)
    #ax1.grid() # Note: grid does not respect z order :(
    #ax2.grid()
    #ax3.grid()
    #ax4.grid()

    shapes_ymax = max( max(sum(histo_mc_sig_norm.values(flow=True))) , max(sum(histo_mc_bkg_norm.values(flow=True))) )
    significance_max = max(max(metric_cum),max(metric_cum_ud))
    significance_min = 0-0.1*significance_max
    ax1.autoscale(axis='y')
    ax2.set_ylim(0.0,1.5*shapes_ymax)
    ax3.set_ylim(significance_min,2.5*significance_max)
    ax4.set_ylim(-0.1,1.2)
    #ax1.set_yscale('log')

    if axisrangex is not None:
        ax1.set_xlim(axisrangex[0],axisrangex[1])
        ax2.set_xlim(axisrangex[0],axisrangex[1])


    return (fig,(extt,extr,extb,extl))



##############################################################
### Wrapper functions for each of the main functionalities ###


### Sanity check of the different reweight points (for a hist that has extra axis to store that) ###
# Old
def check_rwgt(histo_dict):

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



### Dumps the yields and counts for a couple categories into a json ###
# The output of this is used for the CI check
def dump_json_simple(histo_dict,out_name="vvh_yields_simple"):
    out_dict = {}
    hist_to_use = "njets"
    cats_to_check = ["all_events","exactly1lep_exactly1fj","exactly1lep_exactly1fj_STmet1000","exactly1lep_exactly1fj_STmet1000_msd170"]
    for proc_name in histo_dict[hist_to_use].axes["process"]:
        out_dict[proc_name] = {}
        for cat_name in cats_to_check:
            yld = sum(sum(histo_dict[hist_to_use][{"systematic":"nominal", "category":cat_name}].values(flow=True)))
            out_dict[proc_name][cat_name] = [yld,None]

    # Dump counts dict to json
    output_name = f"{out_name}.json"
    with open(output_name,"w") as out_file: json.dump(out_dict, out_file, indent=4)
    print(f"\nSaved json file: {output_name}\n")



### Get the sig and bkg yields and print or dump to json ###
def print_yields(histo_dict,roundat=None,print_counts=False,dump_to_json=True,quiet=False,out_name="yields"):

    # Get ahold of the yields
    yld_dict    = get_yields_per_cat(histo_dict,"njets")
    counts_dict = get_yields_per_cat(histo_dict,"njets_counts")

    group_lst_order = ['Signal', 'Background', 'ttbar', 'VV', 'Vjets', 'QCD', 'single-t', 'ttX', 'VH', 'VVV']

    # Print to screen
    if not quiet:

        ### Print readably ###
        print(f"\n--- Yields ---")
        for cat in yld_dict:
            print(f"\n{cat}")
            for group_name in group_lst_order:
                if group_name == "metric": continue
                yld, err = yld_dict[cat][group_name]
                perr = 100*(err/yld)
                print(f"    {group_name}:  {np.round(yld,roundat)} +- {np.round(perr,2)}%")
            print(f"    -> Metric: {np.round(yld_dict[cat]['metric'][0],3)}")


        ### Print csv, build op as an out string ###

        # Append the header
        out_str = ""
        header = "cat name"
        for proc_name in group_lst_order:
            header = header + f", {proc_name}_val , pm , {proc_name}_pm"
        header = header + ", metric"
        out_str = out_str + header

        # Appead a line for each category, with yields and metric
        for cat in yld_dict:
            line_str = cat
            for group_name in yld_dict[cat]:
                if group_name == "metric": continue
                yld, err = yld_dict[cat][group_name]
                perr = 100*(err/yld)
                line_str = line_str + f" , {np.round(yld,roundat)} , ± , {np.round(perr,2)}%"
            # And also append the metric
            metric = yld_dict[cat]["metric"][0]
            line_str = line_str + f" , {np.round(metric,3)}"
            # Append the string for this line to the out string
            out_str = out_str + f"\n{line_str}"

        # Print the out string to the screen
        print(f"\n\n--- Yields CSV formatted ---\n")
        print(out_str)


    # Dump directly to json
    if dump_to_json:
        out_dict = {"yields":yld_dict, "counts":counts_dict}
        output_name = f"{out_name}.json"
        with open(output_name,"w") as out_file: json.dump(out_dict, out_file, indent=4)
        if not quiet:
            print(f"\n\n--- Yields json formatted ---")
            print(f"\nSaved json file: {output_name}\n")




### Make the plots ###
def make_plots(histo_dict):

    #cat_lst = ["exactly1lep_exactly1fj", "exactly1lep_exactly1fj550", "exactly1lep_exactly1fj550_2j", "exactly1lep_exactly1fj_2j"]

    grouping_dict = append_years(GRP_DICT_FULL,["UL16APV","UL16","UL17","UL18"])
    #grouping_dict = append_years(GRP_DICT_FULL,["UL18"])

    cat_lst = CAT_LST
    var_lst = histo_dict.keys()
    #cat_lst = ["exactly1lep_exactly1fj_STmet1000"]
    #var_lst = ["scalarptsum_lepmet"]

    for cat in cat_lst:
        print("\nCat:",cat)
        for var in var_lst:
            print("\nVar:",var)
            #if var not in ["njets","njets_counts","scalarptsum_lepmet"]: continue # TMP

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

            # Make the figure
            title = f"{cat}__{var}"
            fig,ext_tup = make_vvh_fig(
                histo_mc = histo,
                histo_mc_sig = histo_sig,
                histo_mc_bkg = histo_bkg,
                title=title
            )

            # Save
            save_dir_path = "plots"
            if not os.path.exists("./plots"): os.mkdir("./plots")
            save_dir_path_cat = os.path.join(save_dir_path,cat)
            if not os.path.exists(save_dir_path_cat): os.mkdir(save_dir_path_cat)
            fig.savefig(os.path.join(save_dir_path_cat,title+".png"),bbox_extra_artists=ext_tup,bbox_inches='tight')
            shutil.copyfile(HTML_PC, os.path.join(save_dir_path_cat,"index.php"))



##################################### Main #####################################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument('-y', "--get-yields", action='store_true', help = "Get yields from the pkl file")
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    parser.add_argument('-j', "--dump-json", action='store_true', help = "Dump some yield numbers into a json file")
    parser.add_argument('-o', "--output-name", default='vvh', help = "What to name the outputs")
    args = parser.parse_args()

    # Get the dictionary of histograms from the input pkl file
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))

    # Print total raw events
    #tot_raw = sum(sum(histo_dict["njets_counts"][{"systematic":"nominal", "category":"all_events"}].values(flow=True)))
    #print("Tot raw events:",tot_raw)

    # Which main functionalities to run
    if args.dump_json:
        dump_json_simple(histo_dict,args.output_name)
    if args.get_yields:
        print_yields(histo_dict,out_name=args.output_name+"_yields_sig_bkg",roundat=4,print_counts=False,dump_to_json=True)
    if args.make_plots:
        make_plots(histo_dict)


main()

