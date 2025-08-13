import argparse
import pickle
import gzip
import json
import os, sys
import shutil
import numpy as np
import matplotlib.pyplot as plt
import copy
import warnings
warnings.filterwarnings("ignore", message="List indexing selection is experimental.*")
warnings.filterwarnings("ignore", message="All sumw are zero!*")

import csv
from collections import defaultdict

import utils.plotting_tools as plt_tools
from variables.paths import objsel_cf,get_cutflow,cutflow_yamls_dir,default_cutflow_yaml,default_output_dir
from utils.file_handling import save_array_to_csv

default_CAT_LST = objsel_cf
CAT_LST = default_CAT_LST
rebin = False

clr_map = {
    'EWK': '#d55e00',
    'Other': '#e377c2',
    'QCD': '#8c564b',
    'ST': '#9467bd',
    'WJets': '#d62728',
    'ZJets': '#2ca02c',
    'ttbar': '#ff7f0e',
    'ttx': '#1f77b4',
}

#for sample name csv
sample_name_csv = "/data/userdata/pyli/projects/VVHjj/coffea/vvhjj_coffea/analysis/vbsvvh/variables/sample_names.csv"
def parse_sample_csv(csv_path):
    GRP_DICT_FULL = defaultdict(list)
    bkg_color_map = {}

    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["sample_year"] != "2016preVFP": #temp fix: change to get year automatically in the future
                continue
            cat = row["sample_category"].strip().lower()
            sample_type = row["sample_type"].strip()
            sample_name = row["sample_name"].strip()
            color = row["plotting_colour"].strip()

            if cat == "sig":
                GRP_DICT_FULL["Signal"].append(sample_name)
            elif cat == "bkg":
                GRP_DICT_FULL.setdefault(sample_type, []).append(sample_name)
                if sample_type not in bkg_color_map:
                    bkg_color_map[sample_type] = color

    # Get background color list in the same order as background keys
    bkg_CLR_LST = [bkg_color_map[grp] for grp in GRP_DICT_FULL if grp != "Signal"]

    return dict(GRP_DICT_FULL), bkg_CLR_LST

GRP_DICT_FULL, CLR_LST = parse_sample_csv(sample_name_csv)

# Get sig and bkg yield in all categories
def get_yields_per_cat(histo_dict,var_name):
    out_dict = {}

    # Get the initial grouping dict
    grouping_dict = GRP_DICT_FULL

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
        histo_base = plt_tools.group(histo_base,"year", "temp_year", {"all_year":['2016preVFP','2016postVFP','2017','2018']})
        histo_base = histo_base[{"temp_year":"all_year"}]

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
        punzi = sig/(1+(bkg)**0.5)
        punzi_p1 = sig/(0.11+(bkg)**0.5)
        out_dict[cat]["metric"] = [metric,None] # Don't bother propagating error
        out_dict[cat]["punzi"] = [punzi,None] # Don't bother propagating error
        out_dict[cat]["punzi_p1"] = [punzi_p1,None] # Don't bother propagating error

    return out_dict

def get_ABCD(histo_dict,var_name):
    return 0

# Make the figures for the vvh sudy
def make_vvh_fig(histo_mc,histo_mc_sig,histo_mc_bkg,title="test",axisrangex=None, csv_file = None, var = None,sig_coupling = "nominal"):
    varinfo = {} #this is to store and return varinfo we gathered during make fig
        #[name of quantity, max_quantity, s, b]
        
    
    # Create the figure
    fig, (ax1, ax2, punzi_plot, ax3, ax4) = plt.subplots(
        nrows=5,
        ncols=1,
        figsize=(7,11),
        gridspec_kw={"height_ratios": (3, 1, 1, 1, 1)},
        sharex=True
    )
    #histo_mc = histo_mc_bkg #debug: removed this to try to avoid including sig into filled histos
    fig.subplots_adjust(hspace=.07)
    bkg_samples = list(histo_mc.axes[0])
    colors = [clr_map[sample] for sample in bkg_samples]

    # Plot the stack plot
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        color=colors, #to fix this per sample instead of first n elements)
        ax=ax1,
        zorder=100,
    )
    with np.errstate(divide='ignore', invalid='ignore'):
        # Get the errs on MC and plot them by hand on the stack plot
        histo_mc_sum = histo_mc[{"process_grp":sum}]
        mc_arr = histo_mc_sum.values()
        mc_err_arr = np.sqrt(histo_mc_sum.variances())
        err_p = np.append(mc_arr + mc_err_arr, 0)
        err_m = np.append(mc_arr - mc_err_arr, 0)
        bin_edges_arr = histo_mc_sum.axes[0].edges
        bin_centers_arr = histo_mc_sum.axes[0].centers
        ax1.fill_between(bin_edges_arr,err_m,err_p, step='post', facecolor='none', edgecolor='gray', alpha=0.5, linewidth=0.0, label='MC stat', hatch='/////')

        #calculate yield for later use
        yld_sig = sum(sum(histo_mc_sig.values(flow=True)))
        yld_bkg = sum(sum(histo_mc_bkg.values(flow=True)))

        #plot sig scale to bkg in fig 1
        histo_mc_sig_scale_to_bkg = plt_tools.scale(copy.deepcopy(histo_mc_sig), "process_grp", {"Signal":yld_bkg/yld_sig})
        histo_mc_sig_scale_to_bkg.plot1d(color=["red"], ax=ax1, zorder=100)

        ##2nd plot: normalized s, b
        ## Draw the normalized shapes ##
        # Get normalized hists of sig and bkg
        metric = yld_sig/(yld_bkg**0.5)
        histo_mc_sig_norm         = plt_tools.scale(copy.deepcopy(histo_mc_sig), "process_grp", {"Signal":1.0/yld_sig})
        histo_mc_bkg_norm         = plt_tools.scale(copy.deepcopy(histo_mc_bkg), "process_grp", {"Background":1.0/yld_bkg})
        histo_mc_sig_norm.plot1d(color="red",  ax=ax2, zorder=100)
        histo_mc_bkg_norm.plot1d(color="gray", ax=ax2, zorder=100)

        ##3rd plot: s/sqrt(b + b_error)
        ## Draw the significance ##
        # Get the sig and bkg arrays (Not including flow bins here, overflow should already be handled, and if we have underflow, why?)
        mc_arr = np.where(mc_arr >= 0, mc_arr, 0)
        yld_sig_arr = sum(histo_mc_sig.values())
        #yld_bkg_arr = sum(histo_mc_bkg.values())
        yld_bkg_arr = sum(histo_mc_bkg.values()) + mc_err_arr #changed to use s/sqrt(b + error)
        
        # Get the cumulative signifiance, starting from left
        yld_sig_arr_cum = np.cumsum(yld_sig_arr)
        yld_bkg_arr_cum = np.cumsum(abs(yld_bkg_arr))
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
        plt.text(0.15,0.31, f"Max from left:  {np.round(left_max_y,3)} (at x={np.round(left_max_x,2)}, sig: {np.round(left_s_at_max,2)}, bkg: {np.round(left_b_at_max,1)})", fontsize=9, transform=fig.transFigure)
        plt.text(0.15,0.29, f"Max from right: {np.round(right_max_y,3)} (at x={np.round(right_max_x,2)} , sig: {np.round(right_s_at_max,2)}, bkg: {np.round(right_b_at_max,1)})", fontsize=9, transform=fig.transFigure)

        max_ssqrtb_info = ['s_sqrt_b',right_max_y,right_max_x,right_s_at_max,right_b_at_max]
        varinfo['s_sqrt_b'] = max_ssqrtb_info

        ##4th plot: signal cumulative 
        ## Draw on the fraction of signal retained ##
        yld_sig_arr_cum_frac    = np.cumsum(yld_sig_arr)/yld_sig
        yld_sig_arr_cum_frac_ud = np.flipud(np.cumsum(np.flipud(yld_sig_arr)))/yld_sig
        ax4.scatter(bin_centers_arr,yld_sig_arr_cum_frac,   facecolor='none',edgecolor='black',marker=">",label="Cum. from left", zorder=100)
        ax4.scatter(bin_centers_arr,yld_sig_arr_cum_frac_ud,facecolor='none',edgecolor='black',marker="<",label="Cum. from right", zorder=100)

        ##punzi plot: to get s/(a/2+sqrt(b))for a 
        # Get the cumulative signifiance, starting from left
        a1 = 0
        a2 = 2
        metric_cum_1 = yld_sig_arr_cum/(a1 + np.sqrt(yld_bkg_arr_cum + yld_sig_arr_cum))
        metric_cum_1 = np.nan_to_num(metric_cum_1,nan=0,posinf=0) # Set the nan (from sig and bkg both being 0) to 0
        
        metric_cum_2 = yld_sig_arr_cum/(a2 + np.sqrt(yld_bkg_arr_cum))
        metric_cum_2 = np.nan_to_num(metric_cum_2,nan=0,posinf=0) # Set the nan (from sig and bkg both being 0) to 0

        s_plus_b = yld_sig_arr_cum + yld_bkg_arr_cum
        metric_cum_3 = yld_sig_arr_cum * (s_plus_b -0.333)/np.power(s_plus_b,1.5)
        metric_cum_4 = yld_sig_arr_cum * (yld_bkg_arr_cum -0.333)/np.power(yld_bkg_arr_cum,1.5)
        
        yld_sig_arr_cum_ud = np.flipud(yld_sig_arr_cum_ud) # Flip back so the order is as expected for later use
        yld_bkg_arr_cum_ud = np.flipud(yld_bkg_arr_cum_ud) # Flip back so the order is as expected for later use
        metric_cum_ud_1 = np.flipud(yld_sig_arr_cum_ud/(a1 + np.sqrt(yld_bkg_arr_cum_ud +yld_sig_arr_cum_ud)))
        metric_cum_ud_1 = np.nan_to_num(metric_cum_ud_1,nan=0,posinf=0) # Set the nan (from sig and bkg both being 0) to 0
        metric_cum_ud_2 = np.flipud(yld_sig_arr_cum_ud/(a2 + np.sqrt(yld_bkg_arr_cum_ud)))
        metric_cum_ud_2 = np.nan_to_num(metric_cum_ud_2,nan=0,posinf=0) # Set the nan (from sig and bkg both being 0) to 0
        yld_sig_arr_cum_ud = np.flipud(yld_sig_arr_cum_ud) # Flip back so the order is as expected for later use
        yld_bkg_arr_cum_ud = np.flipud(yld_bkg_arr_cum_ud) # Flip back so the order is as expected for later use


        # Draw it on the plot
        punzi_plot.scatter(bin_centers_arr,metric_cum_1,   facecolor='none',edgecolor='red',marker=">",label="s/sqrt (s+b) Cum. from left", zorder=100)
        punzi_plot.scatter(bin_centers_arr,metric_cum_ud_1,facecolor='none',edgecolor='red',marker="<",label="s/sqrt (s+b) Cum. from right", zorder=100)
        punzi_plot.scatter(bin_centers_arr,metric_cum_2,   facecolor='none',edgecolor='green',marker=">",label="Punzi Cum. from left", zorder=100)
        punzi_plot.scatter(bin_centers_arr,metric_cum_ud_2,facecolor='none',edgecolor='green',marker="<",label="Punzi Cum. from right", zorder=100)

        # Write the max values on the plot
        max_metric_from_left_idx_1  = np.argmax(metric_cum_1)
        max_metric_from_right_idx_1 = np.argmax(metric_cum_ud_1)
        left_max_y_1  = metric_cum_1[max_metric_from_left_idx_1]
        right_max_y_1 = metric_cum_ud_1[max_metric_from_right_idx_1]
        left_max_x_1  = bin_centers_arr[max_metric_from_left_idx_1]
        right_max_x_1 = bin_centers_arr[max_metric_from_right_idx_1]
        left_s_at_max_1  = yld_sig_arr_cum[max_metric_from_left_idx_1]
        right_s_at_max_1 = yld_sig_arr_cum_ud[max_metric_from_right_idx_1]
        left_b_at_max_1  = yld_bkg_arr_cum[max_metric_from_left_idx_1]
        right_b_at_max_1 = yld_bkg_arr_cum_ud[max_metric_from_right_idx_1]
        plt.text(0.15,0.43, f"For a = {a1}: Max from left:  {np.round(left_max_y_1,3)} (at x={np.round(left_max_x_1,2)}, sig: {np.round(left_s_at_max_1,2)}, bkg: {np.round(left_b_at_max_1,1)})", fontsize=7, transform=fig.transFigure)
        plt.text(0.15,0.42, f"For a = {a1}: Max from right: {np.round(right_max_y_1,3)} (at x={np.round(right_max_x_1,2)} , sig: {np.round(right_s_at_max_1,2)}, bkg: {np.round(right_b_at_max_1,1)})", fontsize=7, transform=fig.transFigure)

        max_s_sqrt_splusb = ['s_sqrt_splusb',right_max_y_1,right_max_x_1,right_s_at_max,right_b_at_max]
        varinfo['s_sqrt_splusb'] = max_s_sqrt_splusb

        # for a = 2
        max_metric_from_left_idx_2  = np.argmax(metric_cum_2)
        max_metric_from_right_idx_2 = np.argmax(metric_cum_ud_2)
        left_max_y_2  = metric_cum_2[max_metric_from_left_idx_2]
        right_max_y_2 = metric_cum_ud_2[max_metric_from_right_idx_2]
        left_max_x_2  = bin_centers_arr[max_metric_from_left_idx_2]
        right_max_x_2 = bin_centers_arr[max_metric_from_right_idx_2]
        left_s_at_max_2  = yld_sig_arr_cum[max_metric_from_left_idx_2]
        right_s_at_max_2 = yld_sig_arr_cum_ud[max_metric_from_right_idx_2]
        left_b_at_max_2  = yld_bkg_arr_cum[max_metric_from_left_idx_2]
        right_b_at_max_2 = yld_bkg_arr_cum_ud[max_metric_from_right_idx_2]
        plt.text(0.15,0.41, f"For a = {a2}: Max from left:  {np.round(left_max_y_2,3)} (at x={np.round(left_max_x_2,2)}, sig: {np.round(left_s_at_max_2,2)}, bkg: {np.round(left_b_at_max_2,1)})", fontsize=7, transform=fig.transFigure)
        plt.text(0.15,0.40, f"For a = {a2}: Max from right: {np.round(right_max_y_2,3)} (at x={np.round(right_max_x_2,2)} , sig: {np.round(right_s_at_max_2,2)}, bkg: {np.round(right_b_at_max_2,1)})", fontsize=7, transform=fig.transFigure)

        max_punzi_2 = ['punzi_2',right_max_y_2,right_max_x_2,right_s_at_max_2,right_b_at_max_2]
        varinfo['punzi_2'] = max_punzi_2


        yld_sig_arr_cum_ud = np.flipud(yld_sig_arr_cum_ud) # Flip back so the order is as expected for later use
        yld_bkg_arr_cum_ud = np.flipud(yld_bkg_arr_cum_ud) # Flip back so the order is as expected for later use
        s_plus_b_ud = yld_sig_arr_cum + yld_bkg_arr_cum
        metric_cum_ud_3 = np.flipud(yld_sig_arr_cum_ud * (s_plus_b_ud -0.333)/np.power(s_plus_b_ud,1.5))
        metric_cum_ud_3 = np.nan_to_num(metric_cum_ud_3,nan=0,posinf=0) # Set the nan (from sig and bkg both being 0) to 0
        metric_cum_ud_4 = np.flipud(yld_sig_arr_cum_ud * (yld_bkg_arr_cum_ud -0.333)/np.power(yld_bkg_arr_cum_ud,1.5))
        metric_cum_ud_4 = np.nan_to_num(metric_cum_ud_4,nan=0,posinf=0) # Set the nan (from sig and bkg both being 0) to 0
        yld_sig_arr_cum_ud = np.flipud(yld_sig_arr_cum_ud) # Flip back so the order is as expected for later use
        yld_bkg_arr_cum_ud = np.flipud(yld_bkg_arr_cum_ud) # Flip back so the order is as expected for later use

        #punzi_plot.scatter(bin_centers_arr,metric_cum_3,   facecolor='none',edgecolor='green',marker=">",label="custom1 Cum. from left", zorder=100)
        punzi_plot.scatter(bin_centers_arr,metric_cum_ud_3,facecolor='none',edgecolor='grey',marker="<",label="custom1 Cum. from right", zorder=100)
        #punzi_plot.scatter(bin_centers_arr,metric_cum_4,   facecolor='none',edgecolor='green',marker=">",label="custom1 Cum. from left", zorder=100)
        punzi_plot.scatter(bin_centers_arr,metric_cum_ud_4,facecolor='none',edgecolor='royalblue',marker="<",label="custom2 Cum. from right", zorder=100)



        # custom
        max_metric_from_left_idx_3  = np.argmax(metric_cum_3)
        max_metric_from_right_idx_3 = np.argmax(metric_cum_ud_3)
        left_max_y_3  = metric_cum_3[max_metric_from_left_idx_3]
        right_max_y_3 = metric_cum_ud_3[max_metric_from_right_idx_3]
        left_max_x_3  = bin_centers_arr[max_metric_from_left_idx_3]
        right_max_x_3 = bin_centers_arr[max_metric_from_right_idx_3]
        left_s_at_max_3  = yld_sig_arr_cum[max_metric_from_left_idx_3]
        right_s_at_max_3 = yld_sig_arr_cum_ud[max_metric_from_right_idx_3]
        left_b_at_max_3  = yld_bkg_arr_cum[max_metric_from_left_idx_3]
        right_b_at_max_3 = yld_bkg_arr_cum_ud[max_metric_from_right_idx_3]
        #plt.text(0.15,0.39, f"For custom1: Max from left:  {np.round(left_max_y_3,3)} (at x={np.round(left_max_x_3,2)}, sig: {np.round(left_s_at_max_3,2)}, bkg: {np.round(left_b_at_max_3,1)})", fontsize=7, transform=fig.transFigure)
        plt.text(0.15,0.38, f"For custom1: Max from right: {np.round(right_max_y_3,3)} (at x={np.round(right_max_x_3,2)} , sig: {np.round(right_s_at_max_3,2)}, bkg: {np.round(right_b_at_max_3,1)})", fontsize=7, transform=fig.transFigure)

        max_custom1 = ['custom1',right_max_y_3,right_max_x_3,right_s_at_max_3,right_b_at_max_3]
        varinfo['custom1'] = max_custom1

        max_metric_from_left_idx_4  = np.argmax(metric_cum_4)
        max_metric_from_right_idx_4 = np.argmax(metric_cum_ud_4)
        left_max_y_4  = metric_cum_4[max_metric_from_left_idx_4]
        right_max_y_4 = metric_cum_ud_4[max_metric_from_right_idx_4]
        left_max_x_4  = bin_centers_arr[max_metric_from_left_idx_4]
        right_max_x_4 = bin_centers_arr[max_metric_from_right_idx_4]
        left_s_at_max_4  = yld_sig_arr_cum[max_metric_from_left_idx_4]
        right_s_at_max_4 = yld_sig_arr_cum_ud[max_metric_from_right_idx_4]
        left_b_at_max_4  = yld_bkg_arr_cum[max_metric_from_left_idx_4]
        right_b_at_max_4 = yld_bkg_arr_cum_ud[max_metric_from_right_idx_4]
        #plt.text(0.15,0.37, f"For custom2: Max from left:  {np.round(left_max_y_4,3)} (at x={np.round(left_max_x_4,2)}, sig: {np.round(left_s_at_max_4,2)}, bkg: {np.round(left_b_at_max_4,1)})", fontsize=7, transform=fig.transFigure)
        plt.text(0.15,0.36, f"For custom2: Max from right: {np.round(right_max_y_4,3)} (at x={np.round(right_max_x_4,2)} , sig: {np.round(right_s_at_max_4,2)}, bkg: {np.round(right_b_at_max_4,1)})", fontsize=7, transform=fig.transFigure)

        max_custom2 = ['custom2',right_max_y_4,right_max_x_4,right_s_at_max_4,right_b_at_max_4]
        varinfo['custom2'] = max_custom2

        first_half = right_max_y_1 / 2
        first_half_idx = np.where(metric_cum_ud_1[:max_metric_from_right_idx_1+1] <= first_half)[0]
        if len(first_half_idx) > 0:
            first_half_bin = bin_centers_arr[first_half_idx[-1]]
        else:
            first_half_bin = bin_centers_arr[0]  # fallback if it never drops below half

        # From the right, same logic
        second_half = left_max_y_1 / 2
        second_half_idx = np.where(metric_cum_1[max_metric_from_left_idx_1:] <= second_half)[-1]
        if len(second_half_idx) > 0:
            second_half_bin = bin_centers_arr[max_metric_from_left_idx_1 + second_half_idx[0]]
        else:
            second_half_bin = bin_centers_arr[-1]
        # Add to plot
        #plt.text(0.15, 0.35, f"a={a2}| Half-max range: [{np.round(first_half_bin, 2)}, {np.round(second_half_bin, 2)}]", fontsize=7, transform=fig.transFigure)

        ## Legend, scale the axis, set labels, etc ##
        extr = ax1.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize="12", frameon=False)
        extr = ax2.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize="10", frameon=False)
        extr = ax3.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize="10", frameon=False)
        extr = ax4.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize="10", frameon=False)
        plt.text(0.15,0.86, f"Signal coupling: {sig_coupling}", fontsize = 11, transform=fig.transFigure)
        plt.text(0.15,0.84, f"Sig. yield: {np.round(yld_sig,2)}", fontsize = 11, transform=fig.transFigure)
        plt.text(0.15,0.82, f"Bkg. yield: {np.round(yld_bkg,2)}", fontsize = 11, transform=fig.transFigure)
        plt.text(0.15,0.80, f"Metric: {np.round(metric,3)}", fontsize = 11, transform=fig.transFigure)
        plt.text(0.15,0.78, f"[Note: sig. overlay scaled {np.round(yld_bkg/yld_sig,1)}x]", fontsize = 12, transform=fig.transFigure)

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

        if csv_file is not None:
            row = [var, np.round(right_max_y, 3), np.round(right_max_x, 2), np.round(left_max_y, 3), np.round(left_max_x, 2), 
                    np.round(right_max_y_1, 3), np.round(right_max_x_1, 2), np.round(left_max_y_1, 3), np.round(left_max_x_1, 2), round(right_max_y_2, 3), np.round(right_max_x_2, 2), np.round(left_max_y_2, 3), np.round(left_max_x_2, 2), ]
            with open(csv_file, mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(row)

    return fig,(extt,extr,extb,extl),varinfo


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
    hist_to_use = "nGoodAK4"
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
    yld_dict    = get_yields_per_cat(histo_dict,"nGoodAK4")
    #counts_dict = get_yields_per_cat(histo_dict,"nAK4_counts")

    #group_lst_order = ['Signal', 'Background', 'ttbar', 'VV', 'Vjets', 'QCD', 'single-t', 'ttX', 'VH', 'VVV']
    group_lst_order = list(GRP_DICT_FULL.keys()) + ['Background']
    print("debug",group_lst_order)

    # Print to screen
    if not quiet:

        ### Print readably ###
        print("\n--- Yields ---")
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
                if group_name == "punzi": continue
                if group_name == "punzi_p1": continue
                yld, err = yld_dict[cat][group_name]
                perr = 100*(err/yld)
                line_str = line_str + f" , {np.round(yld,roundat)} , ± , {np.round(perr,2)}%"
            # And also append the metric
            metric = yld_dict[cat]["metric"][0]
            line_str = line_str + f" , {np.round(metric,3)}"
            # Append the string for this line to the out string
            out_str = out_str + f"\n{line_str}"

        # Print the out string to the screen
        print("\n\n--- Yields CSV formatted ---\n")
        print(out_str)


    # Dump directly to json
    if dump_to_json:
        #out_dict = {"yields":yld_dict, "counts":counts_dict}
        out_dict = {"yields": yld_dict}
        output_name = f"{out_name}.json"  # e.g., 'results/my_output.json'

        # Ensure the parent directory exists
        output_dir = os.path.dirname(output_name)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Now write the file
        with open(output_name, "w") as out_file:
            json.dump(out_dict, out_file, indent=4)
        if not quiet:
            print("\n\n--- Yields json formatted ---")
            print(f"\nSaved json file: {output_name}\n")




### Make the plots ###
def make_plots(histo_dict,output_name,sig_coupling):

    grouping_dict = GRP_DICT_FULL

    cat_lst = CAT_LST
    var_lst = histo_dict.keys()

    # Save dir
    save_dir_path = f"{output_name}/"
    if not os.path.exists(f"{output_name}/"): os.makedirs(f"{output_name}/", exist_ok=True)
    
    csv_file = f'{save_dir_path}/var_summary.csv'
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['var', 'right_max', 'right_max_bin', 'left_max', 'left_max_bin', 'right_max_p1', 'right_max_bin_p1', 'left_max_p1', 'left_max_bin_p1','right_max_1p0', 'right_max_bin_1p0', 'left_max_1p0', 'left_max_bin_1p0'])

    var_info = {}
    for cat in cat_lst:
        #print("\nCat:",cat)
        for var in var_lst:
            #if var not in ["njets","njets_counts","scalarptsum_lepmet"]: continue # TMP

            #histo = copy.deepcopy(histo_dict[var][{"systematic":"nominal", "category":cat}]) 
            histo = copy.deepcopy(histo_dict[var][{"category":cat}]) #keep systematic before group as bkg and sig may use different ones

            #can't sum directly over year
            histo = plt_tools.group(histo,"year", "temp_year", {"all_year":['2016preVFP','2016postVFP','2017','2018']})
            histo = histo[{"temp_year":"all_year"}]
            #histo = histo.sum("year")

            # Clean up a bit (rebin, regroup, and handle overflow)

            # if var not in ["njets","nleps","nbtagsl","nbtagsm","njets_counts","nleps_counts","nfatjets","njets_forward","njets_tot"]:
            #    histo = plt_tools.rebin(histo,6)
            
            histo = plt_tools.group(histo,"process","process_grp",grouping_dict)


            # Get one hist of just sig and one of just bkg
            grp_names_bkg_lst = list(grouping_dict.keys()) # All names, still need to drop signal
            
            #print(grp_names_bkg_lst)

            grp_names_bkg_lst.remove("Signal")

            histo_sig = histo[{"process_grp":["Signal"], "systematic":sig_coupling}] #specify coupling for signal
            #histo = histo[{"systematic":"nominal"}] #poor naming, nominal for bkg is really just weight
            histo_bkg = plt_tools.group(histo[{"systematic":"nominal"}],"process_grp","process_grp",{"Background":grp_names_bkg_lst})
            #print("grp_names_bkg_lst",grp_names_bkg_lst)
            #print("histo key", histo)
            
            available_groups = list(histo.axes["process_grp"])


            # Filter only those group names that exist in both the list and the histo
            valid_bkg_groups = [grp for grp in grp_names_bkg_lst if grp in available_groups]
            # Now safely subset


            histo = plt_tools.merge_overflow(histo)
            hist_bkg_test = histo[{"systematic":"nominal","process_grp": valid_bkg_groups}]
            #same binning for plotting
            if rebin:
                histo_bkg,rebin_factor = plt_tools.rebin_to_nonnegative(histo_bkg,"process_grp")
                if rebin_factor >1:
                    hist_bkg_test = plt_tools.rebin(hist_bkg_test,rebin_factor)
                    histo_sig = plt_tools.rebin(histo_sig,rebin_factor)

            # Make the figure
            title = f"{cat}__{var}"
            fig,ext_tup,var_info[cat,var] = make_vvh_fig(
                #histo_mc = histo,
                histo_mc = hist_bkg_test, #debug: added this to input hist list only in bkg
                histo_mc_sig = histo_sig,
                histo_mc_bkg = histo_bkg,
                title=title,
                csv_file = csv_file if cat == list(cat_lst)[-1] else None,
                var = var,
                sig_coupling = sig_coupling, #passing coupling to write on plot
            )
            save_dir_path_cat = os.path.join(save_dir_path,cat)
            if not os.path.exists(save_dir_path_cat): os.mkdir(save_dir_path_cat)
            fig.savefig(os.path.join(save_dir_path_cat,title+".png"),bbox_extra_artists=ext_tup,bbox_inches='tight')

    sqrtb=[['var','max','cut','s','b']]
    splusb=[['var','max','cut','s','b']]
    punzi2=[['var','max','cut','s','b']]
    custom1=[['var','max','cut','s','b']]
    custom2=[['var','max','cut','s','b']]
    for cat in [list(cat_lst)[-1]]:
        #print("\nCat:",cat)
        for var in var_lst:
            sqrtb.append([var]+var_info[cat,var]['s_sqrt_b'][1:])
            splusb.append([var]+var_info[cat,var]['s_sqrt_splusb'][1:])
            punzi2.append([var]+var_info[cat,var]['punzi_2'][1:])
            custom1.append([var]+var_info[cat,var]['custom1'][1:])
            custom2.append([var]+var_info[cat,var]['custom2'][1:])
    save_array_to_csv(sqrtb,os.path.join(save_dir_path_cat,"sqrtb.csv"))
    save_array_to_csv(splusb,os.path.join(save_dir_path_cat,"splusb.csv"))
    save_array_to_csv(punzi2,os.path.join(save_dir_path_cat,"punzi2.csv"))
    save_array_to_csv(custom1,os.path.join(save_dir_path_cat,"custom1.csv"))
    save_array_to_csv(custom2,os.path.join(save_dir_path_cat,"custom2.csv"))



def read_histo_dict(file_path='./', prefix_lst=None, cat_lst=None):
    """
    Reads multiple category histo_dict files and merges into a single dict
    with the same structure as the original combined histo_dict.
    
    Expected filenames: <file_path>/<prefix><cat><suffix>
    """
    suffix='.pkl.gz' #probably always have this file ext
    
    if prefix_lst is None:
        raise ValueError("prefix_list must be provided to know which files to load")
    if cat_lst is None:
        raise ValueError("cat_lst must be provided to know which files to load")

    merged_histo_dict = {}

    for cat in cat_lst:
        if cat == 'all_events':continue
        prefix = prefix_lst[cat]
        file_name = os.path.join(file_path, f"{prefix}_lastcut{suffix}")
        print(f"debug {cat}: {file_name}")
        if not os.path.isfile(file_name):
            print(f"Warning: File not found for category '{cat}': {file_name}")
            continue

        # Load individual category histo dict
        cat_histo_dict = pickle.load(gzip.open(file_name, "rb"))
        print(f"cat:{cat},file {file_name}")
        # Merge into the main histo_dict
        for var, histo in cat_histo_dict.items():
            if var not in merged_histo_dict:
                print(f"var {var}")
                print(histo)
                merged_histo_dict[var] = histo.copy()
            else:
                merged_histo_dict[var] = histo.copy()
                #merged_histo_dict[var][{"category": cat}] = histo.copy()

    return merged_histo_dict

##################################### Main #####################################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--pkl_file_path","-f",default=None, help = "The path to the pkl file or dir if separated")
    parser.add_argument('-y', "--get-yields", action='store_true', help = "Get yields from the pkl file")
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    parser.add_argument('-j', "--dump-json", action='store_true', help = "Dump some yield numbers into a json file")
    parser.add_argument('-o', "--output-name", default=None, help = "What to name the output dir")
    parser.add_argument('-c', "--cutflow", default=None, help = "cutflow to use")
    parser.add_argument('-s', "--sig_coupling", default="nominal", help = "sig_coupling (systematic axis)")
    parser.add_argument('-x','--separate',action='store_true',default=False,help="if pkl files are separated per cut")
    parser.add_argument('--project',default=None,help='add if needed')
    parser.add_argument('--cfyaml',default=None,help='add if need to specify a file')
    parser.add_argument('--histo',default=None,help='none, minus (n_minus_1), one_only, etc.')


    args = parser.parse_args()
    
    
    try:
        cutflow_dict = get_cutflow(args.cfyaml)
        print(f"read cutflow from {args.cfyaml}")
    except:
        try:
            cutflow_dict = get_cutflow(f'{cutflow_yamls_dir}/{args.project}.yaml')
            print(f"read cutflow from project: {cutflow_yamls_dir}/{args.project}.yaml")
        except:
            try:
                cutflow_dict = get_cutflow(default_cutflow_yaml)
                print(f"read from all.yaml")
            except:
                print("cannot find any cutflow to use")
                sys.exit(1)

    if args.cutflow is not None:
        global CAT_LST
        print(cutflow_dict)
        CAT_LST = cutflow_dict[args.cutflow].keys()

    sig_coupling = args.sig_coupling

    #possible pkl file name
    if args.histo is None:
        try_file_name = f'{args.cutflow}.pkl.gz'
        plots_name = 'plots'
    elif args.histo == 'minus' or args.histo == 'm' or args.histo == 'n_minus_1':
        try_file_name = f'{args.cutflow}_m.pkl.gz'
        plots_name = 'n_minus_1'
    elif args.histo == 'last' or args.histo == 'l' or args.histo == 'last_cut':
        try_file_name = f'{args.cutflow}_lastcut.pkl.gz'
        plots_name = 'lastcut'


    # Get the dictionary of histograms from the input pkl file
    if args.separate: 
        print("skip") #skip until I make something up
        #this assume there is a sequence of cutflow in the same project.yaml with +1 cut after the previous cutflow
        '''project_dir = f'/home/users/pyli/vvhjj_coffea/analysis/vbsvvh/outputs/{args.project}/histos/'
        
        file_name_lst = {}
        cat_temp = list(CAT_LST)
        for cat in cat_temp[1:]:  # skip the first since it has no prefix before it
            prefix = cat_temp[:cat_temp.index(cat) + 1]  # e.g. up to current cat
            for target, steps in cutflow_dict.items():
                if list(steps.keys()) == prefix:
                    file_name_lst[cat] = target
                    break  # stop at first match (remove if multiple matches possible)

        histo_dict = read_histo_dict(project_dir,file_name_lst,CAT_LST)
        try: 
            histo_dict = read_histo_dict(project_dir,file_name_lst,CAT_LST)
        except:
            print(f'try using pkl file path {args.pkl_file_path}')
            try:
                histo_dict = read_histo_dict(args.pkl_file_path,file_name_lst,CAT_LST)
            except:
                print(f"cannot get histo: project{project_dir}, pkl_path{args.pkl_file_path}, cutflow: {args.cutflow}")
                print("cut list:",CAT_LST)
                print("file name used: ",file_name_lst)
                sys.exit(1)'''
    else:
        try:    
            histo_dict = pickle.load(gzip.open(args.pkl_file_path))
            print(f"load file {args.pkl_file_path}")
        except:
            try:
                file_path = f'{default_output_dir}/{args.project}/histos/{try_file_name}'
                histo_dict = pickle.load(gzip.open(file_path))
                print(f"load file {file_path}")
            except:
                print(f"cannot get file by {args.pkl_file_path} nor {default_output_dir}/{args.project}/histos/{args.cutflow}.pkl.gz")
                sys.exit(1)

    #print(histo_dict)
    #return coffea.processor.accumulator.dict_accumulator (accumulated hists from processor))

    # Print total raw events
    #tot_raw = sum(sum(histo_dict["njets_counts"][{"systematic":"nominal", "category":"all_events"}].values(flow=True)))
    #print("Tot raw events:",tot_raw)
    if args.output_name is None:
        output_dir = f'{default_output_dir}/{args.project}/{args.cutflow}/{plots_name}/'
    else:
        output_dir = args.output_name
    

    # Which main functionalities to run
    print(f"Setting: cutflow-{list(CAT_LST)} weight-{sig_coupling}")
    if args.dump_json:
        dump_json_simple(histo_dict,output_dir)
    if args.get_yields:
        print_yields(histo_dict,out_name=output_dir+"/summary",roundat=4,print_counts=False,dump_to_json=True)
    if args.make_plots:
        make_plots(histo_dict,output_dir,sig_coupling)
        print(f"Saved all plots in {output_dir}")

main()