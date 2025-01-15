import argparse
import pickle
import gzip
import os
import numpy as np
import matplotlib.pyplot as plt
#import copy
#import hist
#import json

from topcoffea.scripts.make_html import make_html
import ewkcoffea.modules.yield_tools as yt
import ewkcoffea.modules.sample_groupings as sg
import get_wwz_yields as gy

STYLE_DICT = {


    # Input vars in OF SR
    "input_var_of" : {
        "cat_of_interest" : "sr_4l_bdt_of_trn",
        "rebin" : {"run2": 18, "run3" : 30},
        "var_dict" : {
            "mll_wl0_wl1" : {},
            "mllll" : {},
            "absdphi_4l_met" : {},
            "absdphi_wleps_met" : {},
            "absdphi_wl0_met" : {},
            "absdphi_wl1_met" : {},
            "dr_wl0_wl1" : {},
            "dr_zl0_zl1" : {},
            "dr_wleps_zleps" : {},
            "met" : {},
            "mt2" : {},
            "ptl4" : {},
            "scalarptsum_jet" : {},
            "scalarptsum_lepmet" : {},
            "z_lep0_pt" : {},
            "z_lep1_pt" : {},
            "w_lep0_pt" : {},
            "w_lep1_pt" : {},
            "njets" : {"rebin":{"run2":None, "run3":None}},
            "cos_helicity_x" : {},
            "mt_wl0_met" : {},
            "mt_wl1_met" : {},
            "mt_wleps_met" : {},
            "mt_4l_met" : {},
            "dr_wl0_j_min" : {},
            "dr_wl1_j_min" : {},
        },
    },

    # Input vars in SF SR
    "input_var_sf" : {
        "cat_of_interest" : "sr_4l_bdt_sf_trn",
        "rebin" : {"run2": 18, "run3" : 30},
        "var_dict" : {
            "mll_wl0_wl1" : {},
            "mllll" : {},
            "absdphi_4l_met" : {},
            "absdphi_zleps_met" : {},
            "absdphi_wleps_met" : {},
            "absdphi_wl0_met" : {},
            "absdphi_wl1_met" : {},
            "dr_wl0_wl1" : {},
            "dr_zl0_zl1" : {},
            "dr_wleps_zleps" : {},
            "met" : {},
            "mt2" : {},
            "ptl4" : {},
            "scalarptsum_jet" : {},
            "scalarptsum_lepmet" : {},
            "z_lep0_pt" : {},
            "z_lep1_pt" : {},
            "w_lep0_pt" : {},
            "w_lep1_pt" : {},
            "njets" : {"rebin": {"run2":None, "run3":None}},
            "cos_helicity_x" : {},
            "mt_wl0_met" : {},
            "mt_wl1_met" : {},
            "mt_wleps_met" : {},
            "mt_4l_met" : {},
            "dr_wl0_j_min" : {},
            "dr_wl1_j_min" : {},
        },
    },

    # BDT scores in OF SR
    "score_of" : {
        "cat_of_interest" : "sr_4l_bdt_of_trn",
        "rebin" : {"run2": 18, "run3" : 30},
        "var_dict" : {
            "bdt_of_wwz" : {},
            "bdt_of_zh" : {},
            "bdt_of_bkg" : {},
            "bdt_of_wwz_m_zh" : {},
        },
    },

    # BDT scores in OF SR
    "score_sf" : {
        "cat_of_interest" : "sr_4l_bdt_sf_trn",
        "rebin" : {"run2": 18, "run3" : 30},
        "var_dict" : {
            "bdt_sf_wwz" : {},
            "bdt_sf_zh" : {},
            "bdt_sf_bkg" : {},
            "bdt_sf_wwz_m_zh" : {},
        },
    },
}

# Takes a mc hist and data hist and plots both
def make_public_fig(histo_mc,histo_data=None,title="test",unit_norm_bool=False,axisrangex=None):

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
        color=gy.CLR_LST,
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

    # Scale the axis and set labels
    if axisrangex is not None:
        ax.set_xlim(axisrangex[0],axisrangex[1])
        rax.set_xlim(axisrangex[0],axisrangex[1])
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
def make_plots(histo_dict,grouping_mc,grouping_data,save_dir_path,apply_nsf_to_cr=False,year="run2"):

    # Set up the output dir if it does not exist
    if not os.path.exists(save_dir_path):
        os.mkdir(save_dir_path)

    # Loop over the groups of plots to make
    for group_name in STYLE_DICT.keys():

        # Make a sub dir for this group of plots
        save_dir_path_group = os.path.join(save_dir_path,group_name)
        if not os.path.exists(save_dir_path_group):
            os.mkdir(save_dir_path_group)

        # Which SR or CR is relevant
        cat_name = STYLE_DICT[group_name]["cat_of_interest"]

        # Loop over the plots in this group
        for var_name in STYLE_DICT[group_name]["var_dict"].keys():
            print(f"\nVar name: {var_name}")

            histo_cat = histo_dict[var_name][{"systematic":"nominal", "category":cat_name}]

            # Rebin and set x range
            rangex = None
            rangex = [50,150]
            if var_name not in ["njets"]:
                rebin_factor = STYLE_DICT[group_name]["rebin"][year]
                if "rebin" in STYLE_DICT[group_name]["var_dict"][var_name]:
                    rebin_factor = STYLE_DICT[group_name]["var_dict"][var_name]["rebin"][year]
                histo_cat = gy.rebin(histo_cat,rebin_factor)

            # Group the mc and data samples
            histo_grouped_mc = gy.group(histo_cat,"process","process_grp",grouping_mc)
            histo_grouped_data = gy.group(histo_cat,"process","process_grp",grouping_data)

            # Merge overflow into last bin (so it shows up in the plot)
            histo_grouped_data = gy.merge_overflow(histo_grouped_data)
            histo_grouped_mc = gy.merge_overflow(histo_grouped_mc)

            # Make figure
            title = f"{group_name}_{var_name}"
            print("Making: ",title)
            #fig = make_public_fig(histo_grouped_mc,histo_grouped_data,axisrangex=rangex,title=title)

            # Save figure
            #fig.savefig(os.path.join(save_dir_path_group,title+".pdf"))
            #fig.savefig(os.path.join(save_dir_path_group,title+".png"))

        # Make html for this sub dir
        make_html(os.path.join(os.getcwd(),save_dir_path_group))



################### Main ###################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default="plots_public", help = "The path the output files should be saved to")
    parser.add_argument('-y', "--get-yields", action='store_true', help = "Get yields from the pkl file")
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    parser.add_argument('-u', "--ul-year", default='run2', help = "Which year to process", choices=["all","run2","run3","y22","y23","UL16APV","UL16","UL17","UL18","2022","2022EE","2023","2023BPix"])
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))

    sample_dict_mc = sg.create_mc_sample_dict(args.ul_year)
    sample_dict_data = sg.create_data_sample_dict(args.ul_year)
    out_path = args.output_path

    make_plots(histo_dict,sample_dict_mc,sample_dict_data,save_dir_path=out_path,apply_nsf_to_cr=False,year=args.ul_year)

if __name__ == "__main__":
    main()
