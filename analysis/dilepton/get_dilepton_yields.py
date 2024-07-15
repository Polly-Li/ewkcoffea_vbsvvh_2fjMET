import argparse
import pickle
import gzip
import os
import numpy as np
import matplotlib.pyplot as plt
import copy
import hist

from topcoffea.scripts.make_html import make_html

import ewkcoffea.modules.yield_tools as yt
import ewkcoffea.modules.sample_groupings_dilepton as sg

# This script opens a pkl file of histograms produced by wwz processor
# Reads the histograms and dumps out the yields for each group of processes
# Example usage: python get_yld_check.py histos/tmp_histo.pkl.gz -y

CLR_LST = ["red","blue","#F09B9B","#00D091"]

# Names of the cut-based and BDT SRs
SR_SF_mumu = ["dil_presel_sf_mumu"]
SR_SF_ee = ["dil_presel_sf_ee"]


TMP_VAR_LST = [
    "j0pt",
    "njets",
    "nbtagsl",
    "nleps",
    "met",
    "l0pt",
]


SOVERROOTB = "$S/\sqrt{B}$"
SOVERROOTSPLUSB = "$S/\sqrt{S+B}$"


# Get the list of categories on the sparese axis
def get_axis_cats(histo,axis_name):
    process_list = [x for x in histo.axes[axis_name]]
    return process_list


# Merges the last bin (overflow) into the second to last bin, zeros the content of the last bin, returns a new hist
# Note assumes just one axis!
def merge_overflow(hin):
    hout = copy.deepcopy(hin)
    for cat_idx,arr in enumerate(hout.values(flow=True)):
        hout.values(flow=True)[cat_idx][-2] += hout.values(flow=True)[cat_idx][-1]
        hout.values(flow=True)[cat_idx][-1] = 0
        hout.variances(flow=True)[cat_idx][-2] += hout.variances(flow=True)[cat_idx][-1]
        hout.variances(flow=True)[cat_idx][-1] = 0
    return hout


# Rebin according to https://github.com/CoffeaTeam/coffea/discussions/705
def rebin(histo,factor):
    return histo[..., ::hist.rebin(factor)]


# Regroup categories (e.g. processes)
def group(h, oldname, newname, grouping):

    # Build up a grouping dict that drops any proc that is not in our h
    grouping_slim = {}
    proc_lst = get_axis_cats(h,oldname)
    for grouping_name in grouping.keys():
        for proc in grouping[grouping_name]:
            if proc in proc_lst:
                if grouping_name not in grouping_slim:
                    grouping_slim[grouping_name] = []
                grouping_slim[grouping_name].append(proc)
            #else:
            #    print(f"WARNING: process {proc} not in this hist")

    # From Nick: https://github.com/CoffeaTeam/coffea/discussions/705#discussioncomment-4604211
    hnew = hist.Hist(
        hist.axis.StrCategory(grouping_slim, name=newname),
        *(ax for ax in h.axes if ax.name != oldname),
        storage=h.storage_type(),
    )
    for i, indices in enumerate(grouping_slim.values()):
        hnew.view(flow=True)[i] = h[{oldname: indices}][{oldname: sum}].view(flow=True)

    return hnew


# Takes a mc hist and data hist and plots both
def make_cr_fig(histo_mc,histo_data=None,title="test",unit_norm_bool=False):

    # Create the figure
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(7,7),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.07)

    # Plot the mc
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        color=CLR_LST,
        ax=ax,
    )
    # Plot the data
    if histo_data is not None:
        histo_data.plot1d(
            stack=False,
            histtype="errorbar",
            color="k",
            ax=ax,
            w2=histo_data.variances(),
            w2method="sqrt",
        )
    # Plot a dummy hist on rax to get the label to show up
    histo_mc.plot1d(alpha=0, ax=rax)

    ### Get the errs on MC and plot them by hand ###
    histo_mc_sum = histo_mc[{"process_grp":sum}]
    mc_arr = histo_mc_sum.values()
    mc_err_arr = np.sqrt(histo_mc_sum.variances())
    err_p = np.append(mc_arr + mc_err_arr, 0)
    err_m = np.append(mc_arr - mc_err_arr, 0)
    bin_edges_arr = histo_mc_sum.axes[0].edges
    bin_centers_arr = histo_mc_sum.axes[0].centers
    ax.fill_between(bin_edges_arr,err_m,err_p, step='post', facecolor='none', edgecolor='gray', alpha=0.5, linewidth=0.0, label='MC stat', hatch='/////')

    ### Get the errs on data and ratios and plot them by hand ###
    if histo_data is not None:
        histo_data_sum = histo_data[{"process_grp":sum}]

        data_arr = histo_data_sum.values()
        data_err_arr = np.sqrt(histo_data_sum.variances())

        err_ratio_p = np.append(1+mc_err_arr/mc_arr,1)
        err_ratio_m = np.append(1-mc_err_arr/mc_arr,1)

        data_ratio_err_p = (data_arr + data_err_arr)/mc_arr
        data_ratio_err_m = (data_arr - data_err_arr)/mc_arr

        rax.fill_between(bin_edges_arr,err_ratio_m,err_ratio_p,step='post', facecolor='none',edgecolor='gray', label='MC stat', linewidth=0.0, hatch='/////',alpha=0.5)
        rax.scatter(bin_centers_arr,data_arr/mc_arr,facecolor='black',edgecolor='black',marker="o")
        rax.vlines(bin_centers_arr,data_ratio_err_p,data_ratio_err_m,color='k')

    # Scale the y axis and labels
    ax.legend(fontsize="12")
    ax.set_title(title)
    ax.autoscale(axis='y')
    ax.set_xlabel(None)
    rax.set_ylabel('Ratio')
    rax.set_ylim(0.0,2.0)
    rax.axhline(1.0,linestyle="-",color="k",linewidth=1)
    ax.tick_params(axis='y', labelsize=16)
    rax.tick_params(axis='x', labelsize=16)
    #ax.set_yscale('log')

    return fig


# Main function for making CR plots
def make_plots(histo_dict,grouping_mc,grouping_data,save_dir_path):

    histo = histo_dict["njets"]
    yld_dict_mc   = yt.get_yields(histo,grouping_mc,blind=False)
    yld_dict_data = yt.get_yields(histo,grouping_data,blind=False)

    for var_name in histo_dict.keys():
        # Skip over some variables if you want to
        #if "counts" in var_name: continue

        # Get the relevant histogram from the input dict
        #print(f"\n{var_name}")
        histo_orig = histo_dict[var_name]

        # Just plot nominal syst for now
        histo_orig = histo_orig[{"systematic":"nominal"}]

        # Loop over categories and make plots for each
        for cat_name in histo_orig.axes["category"]:
            # Skip some of the cats if you want to
            #if "bdt" in cat_name: continue
            #if cat_name not in ["sr_4l_sf_incl", "sr_4l_of_incl"]: continue
            #print(cat_name)

            # Make a copy so changes to binning do not propagate to next loop
            histo = copy.deepcopy(histo_orig)

            # Rebin if continous variable
            if var_name not in ["nleps","njets","nbtagsl","nbtagsm"]:
                histo = rebin(histo,6)

            histo_cat = histo[{"category":cat_name}]

            # Group the mc and data samples
            histo_grouped_mc = group(histo_cat,"process","process_grp",grouping_mc)
            histo_grouped_data = group(histo_cat,"process","process_grp",grouping_data)

            ######
            # Print stuff if you want to
            if (cat_name == "mumu_2l_sf" or cat_name == "ee_2l_sf") and var_name == "nleps":
                print(f"\n{cat_name} {var_name}:")
                print("Yields")
                data = sum(sum(histo_grouped_data.values(flow=True)))
                data_error = sum(sum(histo_grouped_data.variances(flow=True))) ** 0.5
                mc = sum(sum(histo_grouped_mc.values(flow=True)))
                mc_error = sum(sum(histo_grouped_mc.variances(flow=True))) ** 0.5
                data_over_mc_ratio = data / mc
                data_over_mc_ratio_error = data_over_mc_ratio * ((data_error / data)**2 + (mc_error / mc)**2)**0.5
                print("mc:",mc)
                print("mc Error:",mc_error)
                print("data:",data)
                print("data Error:",data_error)
                print("data/mc:", data_over_mc_ratio)
                print("data/mc Error:", data_over_mc_ratio_error)
            continue
            #####

            # Merge overflow into last bin (so it shows up in the plot)
            histo_grouped_data = merge_overflow(histo_grouped_data)
            histo_grouped_mc = merge_overflow(histo_grouped_mc)

            # Make figure
            title = f"{cat_name}_{var_name}"
            print("Making: ",title)
            fig = make_cr_fig(histo_grouped_mc,histo_grouped_data,title=title)

            # Save
            save_dir_path_cat = os.path.join(save_dir_path,cat_name)
            if not os.path.exists(save_dir_path_cat): os.mkdir(save_dir_path_cat)
            fig.savefig(os.path.join(save_dir_path_cat,title+".pdf"))
            fig.savefig(os.path.join(save_dir_path_cat,title+".png"))

            make_html(os.path.join(os.getcwd(),save_dir_path_cat))

################### Main ###################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default="plots", help = "The path the output files should be saved to")
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    parser.add_argument('-u', "--ul-year", default='run2', help = "Which year to process", choices=["all","run2","run3","yr22","yr23","UL16APV","UL16","UL17","UL18","2022","2022EE","2023","2023BPix"])
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))
    sample_dict_mc = sg.create_mc_sample_dict(args.ul_year)
    sample_dict_data = sg.create_data_sample_dict(args.ul_year)
    out_path = args.output_path

    # Make plots
    if args.make_plots:
        make_plots(histo_dict,sample_dict_mc,sample_dict_data,save_dir_path=out_path)

if __name__ == "__main__":
    main()

