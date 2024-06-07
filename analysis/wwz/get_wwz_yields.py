import argparse
import pickle
import json
import gzip
import os
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

import yld_dicts_for_comp as yd

# This script opens a pkl file of histograms produced by wwz processor
# Reads the histograms and dumps out the yields for each group of processes
# Example usage: python get_yld_check.py histos/tmp_histo.pkl.gz -y

# Colors in VVV observation
#ZZ    = (240, 155, 205)  #F09B9B
#ttZ   = (0, 208, 145) #00D091
#WZ    = (163, 155, 47) #A39B2F
#tWZ   = (205, 240, 155) #CDF09B
#Other = (205, 205, 205) #CDCDCD
CLR_LST = ["red","blue","#F09B9B","#00D091","#CDF09B","#A39B2F","#CDCDCD"]
#CLR_LST = ["#F09B9B","#00D091","#CDF09B"]

# Names of the cut-based and BDT SRs
SR_SF_CB = ["sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C"]
SR_OF_CB = ["sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4"]
SR_SF_BDT = ["sr_4l_bdt_sf_1", "sr_4l_bdt_sf_2", "sr_4l_bdt_sf_3", "sr_4l_bdt_sf_4", "sr_4l_bdt_sf_5", "sr_4l_bdt_sf_6", "sr_4l_bdt_sf_7", "sr_4l_bdt_sf_8"]
SR_OF_BDT = ["sr_4l_bdt_of_1", "sr_4l_bdt_of_2", "sr_4l_bdt_of_3", "sr_4l_bdt_of_4", "sr_4l_bdt_of_5", "sr_4l_bdt_of_6", "sr_4l_bdt_of_7", "sr_4l_bdt_of_8"]

BDT_INPUT_LST = [
    "mll_wl0_wl1",
    "mllll",
    "absdphi_4l_met",
    "absdphi_zleps_met",
    "absdphi_wleps_met",
    "absdphi_wl0_met",
    "absdphi_wl1_met",
    "dr_wl0_wl1",
    "dr_zl0_zl1",
    "dr_wleps_zleps",
    "met",
    "mt2",
    "ptl4",
    "scalarptsum_jet",
    "scalarptsum_lepmet",
    "z_lep0_pt",
    "z_lep1_pt",
    "w_lep0_pt",
    "w_lep1_pt",
    "njets",
    "cos_helicity_x",
    "mt_wl0_met",
    "mt_wl1_met",
    "mt_wleps_met",
    "mt_4l_met",
    "dr_wl0_j_min",
    "dr_wl1_j_min",
]

BDT_SCORE_LST = [
    "bdt_of_wwz",
    "bdt_sf_wwz",
    "bdt_of_zh",
    "bdt_sf_zh",
    "bdt_of_bkg",
    "bdt_sf_bkg",
    "bdt_of_wwz_m_zh",
    "bdt_sf_wwz_m_zh",
    #"bdt_of_bin",
    #"bdt_sf_bin",
]

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







################### Getting and printing yields ###################

# Get the yields in the SR
def get_yields(histos_dict,sample_dict,raw_counts=False,quiet=True,blind=True,systematic_name="nominal"):

    yld_dict = {}

    # Look at the yields in one histo (e.g. njets)
    if raw_counts: dense_axis = "njets_counts"
    else: dense_axis = "njets"
    for proc_name in sample_dict.keys():
        yld_dict[proc_name] = {}
        for cat_name in histos_dict[dense_axis].axes["category"]:
            #if "bdt" in cat_name: continue # TMP!!!
            if blind and (("data" in proc_name) and (not cat_name.startswith("cr_"))):
                # If this is data and we're not in a CR category, put placeholder numbers for now
                yld_dict[proc_name][cat_name] = [-999,-999]
            else:
                val = sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name],"systematic":systematic_name}].values(flow=True)))
                var = sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name],"systematic":systematic_name}].variances(flow=True)))
                yld_dict[proc_name][cat_name] = [val,var]

    # Print to screen
    if not quiet:
        for proc in yld_dict.keys():
            print(f"\n{proc}:")
            for cat in yld_dict[proc].keys():
                val = yld_dict[proc][cat]
                print(f"\t{cat}: {val}")

    return yld_dict


# Gets the process sums for S and B and gets metrics e.g. S/sqrt(B) and puts it into the dict
# Hard coded for the summed values (e.g. looking for "ZH" not "GluGluZH","qqToZHToZTo2L")
def put_proc_row_sums(yld_dict,sr_cat_lst, sig_lst=sg.SIG_LST,bkg_lst=sg.BKG_LST):
    # Build up the empty dicts for sig and bkg that we can later fill and then put into the yld dict
    # Will look something like this: {"sr_4l_sf_A":[0,0], "sr_4l_sf_B":[0,0], "sr_4l_sf_C":[0,0], "sr_4l_of_1":[0,0], "sr_4l_of_2":[0,0], "sr_4l_of_3":[0,0], "sr_4l_of_4":[0,0]}
    sig_sum = {}
    bkg_sum = {}
    for sr_cat_name in sr_cat_lst:
        sig_sum[sr_cat_name] = [0,0]
        bkg_sum[sr_cat_name] = [0,0]

    # Finding sums
    for proc in yld_dict.keys():
        print(proc)
        for cat in yld_dict[proc].keys():
            if cat not in sig_sum: continue
            val,var = yld_dict[proc][cat]
            print("   ",cat,val)
            if proc in sig_lst:
                sig_sum[cat][0] += val
                sig_sum[cat][1] += var
            if proc in bkg_lst:
                bkg_sum[cat][0] += val
                bkg_sum[cat][1] += var

    # Finding metrics, and putting sums and metrics into the yld dict
    if SOVERROOTB not in yld_dict:      yld_dict[SOVERROOTB] = {}
    if SOVERROOTSPLUSB not in yld_dict: yld_dict[SOVERROOTSPLUSB] = {}
    if "Sig" not in yld_dict: yld_dict["Sig"] = {}
    if "Bkg" not in yld_dict: yld_dict["Bkg"] = {}
    if "Zmetric" not in yld_dict: yld_dict["Zmetric"] = {}
    for cat in sig_sum.keys():
        s = sig_sum[cat][0]
        b = bkg_sum[cat][0]
        s_var = sig_sum[cat][1]
        b_var = bkg_sum[cat][1]
        yld_dict[SOVERROOTB][cat]      = [s/math.sqrt(b) , None]
        yld_dict[SOVERROOTSPLUSB][cat] = [s/math.sqrt(s+b) , None]
        yld_dict["Zmetric"][cat] = [math.sqrt(2 * ((s + b) * math.log(1 + s / b) - s)), None] # Eq 18 https://cds.cern.ch/record/2203244/files/1087459_109-114.pdf
        yld_dict["Sig"][cat] = [s, s_var]
        yld_dict["Bkg"][cat] = [b, b_var]


# Gets the sums of categoreis (assumed to be columns in the input dict) and puts them into the dict
# Special handling for rows that are metrics (e.g. s/sqrt(b)), sums these in quadrature
def put_cat_col_sums(yld_dict,sr_sf_lst,sr_of_lst,metrics_names_lst=["Zmetric",SOVERROOTB,SOVERROOTSPLUSB],tag=""):

    # The full SR list should be the sf and of together
    sr_lst = sr_sf_lst + sr_of_lst

    # Loop over rows (processes) and sum columns together, fill the result into new_dict
    new_dict = {}
    for proc in yld_dict:
        sr_sf_val = 0
        sr_sf_var = 0
        sr_of_val = 0
        sr_of_var = 0
        sr_val = 0
        sr_var = 0
        for cat in yld_dict[proc]:
            val = yld_dict[proc][cat][0]
            var = yld_dict[proc][cat][1]
            if cat in sr_sf_lst:
                if proc in metrics_names_lst:
                    sr_sf_val += val*val
                else:
                    sr_sf_val += val
                    sr_sf_var += var
            if cat in sr_of_lst:
                if proc in metrics_names_lst:
                    sr_of_val += val*val
                else:
                    sr_of_val += val
                    sr_of_var += var
            if cat in sr_lst:
                if proc in metrics_names_lst:
                    sr_val += val*val
                else:
                    sr_val += val
                    sr_var += var

        # Fill our new_dict with what we've computed
        new_dict[proc] = {}
        if proc in metrics_names_lst:
            new_dict[proc]["sr_sf_all"] = [np.sqrt(sr_sf_val),None]
            new_dict[proc]["sr_of_all"] = [np.sqrt(sr_of_val),None]
            new_dict[proc]["sr_all"]    = [np.sqrt(sr_val),None]
        else:
            new_dict[proc]["sr_sf_all"] = [sr_sf_val, sr_sf_var]
            new_dict[proc]["sr_of_all"] = [sr_of_val, sr_of_var]
            new_dict[proc]["sr_all"]    = [sr_val, sr_var]

    # Put the columns into the yld_dict
    for proc in new_dict:
        yld_dict[proc][f"sr_sf_all{tag}"] = new_dict[proc]["sr_sf_all"]
        yld_dict[proc][f"sr_of_all{tag}"] = new_dict[proc]["sr_of_all"]
        yld_dict[proc][f"sr_all{tag}"] = new_dict[proc]["sr_all"]


# Print yields
def print_yields(ul_year,yld_dict_in,cats_to_print,procs_to_print,ref_dict=yd.EWK_REF_NOSF,print_fom=True,hlines=[]):

    # Get err from var
    def get_err_from_var(in_dict):
        out_dict = {}
        for proc in in_dict:
            out_dict[proc] = {}
            for cat in in_dict[proc]:
                if in_dict[proc][cat][1] is None: var = None
                else: var = np.sqrt(in_dict[proc][cat][1])
                out_dict[proc][cat] = [in_dict[proc][cat][0],var]
        return out_dict

    yld_dict = get_err_from_var(yld_dict_in)

    # Print the yields directly
    mlt.print_latex_yield_table(
        yld_dict,
        tag="All yields",
        key_order=procs_to_print,
        subkey_order=cats_to_print,
        print_begin_info=True,
        print_end_info=True,
        print_errs=True,
        column_variable="subkeys",
        size="tiny",
        hz_line_lst=[6],
    )
    #exit()

    ### Compare with other yields, print comparison ###

    tag1 = "New"
    tag2 = "Ref"

    yld_dict_comp = get_err_from_var(ref_dict)

    yld_dict_1 = copy.deepcopy(yld_dict)
    yld_dict_2 = copy.deepcopy(yld_dict_comp)

    pdiff_dict = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="percent_diff",inpercent=True)
    diff_dict  = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="absolute_diff")

    mlt.print_begin()
    mlt.print_latex_yield_table(yld_dict_1,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag1,hz_line_lst=hlines,print_errs=True,size="tiny",column_variable="keys")
    mlt.print_latex_yield_table(yld_dict_2,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag2,hz_line_lst=hlines,print_errs=True,size="tiny",column_variable="keys")
    mlt.print_latex_yield_table(pdiff_dict,key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Percent diff between {tag1} and {tag2}",hz_line_lst=hlines,size="tiny",column_variable="keys")
    mlt.print_latex_yield_table(diff_dict, key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Diff between {tag1} and {tag2}",hz_line_lst=hlines,size="tiny",column_variable="keys")
    mlt.print_end()


# Dump the counts dict to a latex table
def print_counts(counts_dict):

    cats_to_print = ["all_events", "4l_presel", "sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]

    # Print the yields directly
    mlt.print_latex_yield_table(
        counts_dict,
        tag="Raw MC counts (ewkcoffea)",
        key_order=counts_dict.keys(),
        subkey_order=cats_to_print,
        print_begin_info=True,
        print_end_info=True,
        column_variable="subkeys",
    )


# This should maybe be in a different script
################### Hist manipulation and plotting ###################


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


# Plots a hist
def make_single_fig(histo_mc,ax_to_overlay="process",err_p=None,err_m=None,ylims=None,unit_norm_bool=False,title=None,fig_size=(12,7)):
    #print("\nPlotting values:",histo.values())
    fig, ax = plt.subplots(1, 1, figsize=fig_size)

    # Plot the mc
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        color=CLR_LST,
        yerr=True,
        overlay=ax_to_overlay
    )

    # Set title and y lims
    if title is not None: plt.title(title)
    if ylims is None:
        ax.autoscale(axis='y')
    else:
        ax.set_ylim(ylims)

    # Draw errors
    if (err_p is not None) and (err_m is not None):
        bin_edges_arr = histo_mc.axes[1].edges
        bin_centers_arr = histo_mc.axes[1].centers
        ax.fill_between(bin_edges_arr,err_m,err_p, step='post', facecolor='none', edgecolor='gray', alpha=0.5, linewidth=0.0, label='MC stat', hatch='/////')

    plt.legend()
    return fig


# Takes a mc hist and data hist and plots both
def make_syst_fig(histo_mc,mc_up_arr,mc_do_arr,syst,histo_data=None,title="test",unit_norm_bool=False):

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

    # Plot the syst
    histo_mc_sum = histo_mc[{"process_grp":sum}]
    bin_edges_arr = histo_mc_sum.axes[0].edges
    bin_centers_arr = histo_mc_sum.axes[0].centers
    ax.stairs(mc_up_arr, bin_edges_arr, color="cyan", linestyle="--", label=f'{syst} up')
    ax.stairs(mc_do_arr, bin_edges_arr, color="magenta", linestyle="--", label=f'{syst} down')

    # Plot the syst on ratio plots
    mc_arr = histo_mc_sum.values()
    rax.scatter(bin_centers_arr,mc_up_arr/mc_arr,facecolor='cyan',edgecolor='cyan',marker="o")
    rax.scatter(bin_centers_arr,mc_do_arr/mc_arr,facecolor='magenta',edgecolor='magenta',marker="o")

    # Scale the y axis and labels
    ax.legend(fontsize="12")
    ax.set_title(title)
    ax.autoscale(axis='y')
    ax.set_xlabel(None)
    rax.set_ylabel('Ratio')
    rax.set_ylim(0.9,1.1)
    rax.axhline(1.0,linestyle="-",color="k",linewidth=1)
    ax.tick_params(axis='y', labelsize=16)
    rax.tick_params(axis='x', labelsize=16)

    return fig


# IN PROGRESS
# Main function for checking individual systematics
def make_syst_plots(histo_dict,grouping_mc,grouping_data,save_dir_path,year):

    for var_name in histo_dict.keys():
        #print(f"\n{var_name}")
        if var_name not in TMP_VAR_LST: continue
        histo = histo_dict[var_name]

        cat_lst = [
            "sr_4l_sf_A",
            "sr_4l_sf_B",
            "sr_4l_sf_C",
            "sr_4l_of_1",
            "sr_4l_of_2",
            "sr_4l_of_3",
            "sr_4l_of_4",
            #"cr_4l_sf",
            #"cr_4l_btag_sf_offZ_met80",
            #"cr_4l_btag_of",
            #"sr_4l_of_incl",
            #"sr_4l_sf_incl",
        ]

        # Rebin if continous variable
        if var_name not in ["njets","nbtagsl","nleps"]:
            histo = rebin(histo,6)

        # Get the list of systematic base names (i.e. without the up and down tags)
        # Assumes each syst has a "systnameUp" and a "systnameDown" category on the systematic axis
        syst_var_lst = []
        all_syst_var_lst = histo.axes["systematic"]
        for syst_var_name in all_syst_var_lst:
            if syst_var_name.endswith("Up"):
                syst_name_base = syst_var_name.replace("Up","")
                if syst_name_base not in syst_var_lst:
                    syst_var_lst.append(syst_name_base)

        for cat in cat_lst:
            print("\n",cat)
            if "cr_4l_of" not in cat and var_name == "j0pt": continue
            histo_cat = histo[{"category":cat}]
            histo_grouped_mc = group(histo_cat,"process","process_grp",grouping_mc)
            histo_grouped_data = group(histo_cat,"process","process_grp",grouping_data)

            mc_nom   = merge_overflow(histo_grouped_mc[{"systematic":"nominal"}])
            data_nom = merge_overflow(histo_grouped_data[{"systematic":"nominal"}])

            for syst in syst_var_lst:
                #if "btag" not in syst: continue
                #if "uncorrelated" not in syst: continue
                #if "lepSF" not in syst: continue
                #if "PreFiring" not in syst: continue
                #if "PU" not in syst: continue
                #if "ISR" not in syst and "FSR" not in syst: continue
                if "renorm" not in syst and "fact" not in syst: continue

                # Skip the variations that don't apply (TODO: why are these in the hist to begin with??)
                if year == "UL16APV": blacklist_years = ["2016","2017","2018","2022","2022EE"]
                if year == "UL16": blacklist_years = ["2016APV","2017","2018","2022","2022EE"]
                if year == "UL17": blacklist_years = ["2016APV","2016","2018","2022","2022EE"]
                if year == "UL18": blacklist_years = ["2016APV","2016","2017","2022","2022EE"]
                if year == "2022": blacklist_years = ["2016APV","2016","2017","2018","2022EE"]
                if year == "2022EE": blacklist_years = ["2016APV","2016","2017","2018","2022"]
                if year == "all": blacklist_years = []
                skip = False
                for y in blacklist_years:
                    if syst.endswith(y):
                        skip = True
                if skip: continue

                if year == "all": yeartag = "FullR2"
                else: yeartag = year

                mc_up     = merge_overflow(histo_grouped_mc[{"systematic":f"{syst}Up"}])
                mc_down   = merge_overflow(histo_grouped_mc[{"systematic":f"{syst}Down"}])
                data_up   = merge_overflow(histo_grouped_data[{"systematic":f"{syst}Up"}])
                data_down = merge_overflow(histo_grouped_data[{"systematic":f"{syst}Down"}])

                mc_up_arr = mc_up[{"process_grp":sum}].values()
                mc_down_arr = mc_down[{"process_grp":sum}].values()

                # Print individual syst numbers
                #if var_name != "nleps": continue
                #n = sum(sum(mc_nom.values()))
                #u = sum(mc_up_arr)
                #d = sum(mc_down_arr)
                #print("\n",syst)
                #print("nom",n)
                #print("up",u)
                #print("do",d)
                #r_up = abs((n-u)/n)
                #r_do = abs((n-d)/n)
                #r = (r_up+r_do)/2
                #print("err",np.round(100*abs(n-u)/n,1),"%")
                #print("err up",np.round(100*r_up,1),"%")
                #print("err do",np.round(100*r_do,1),"%")
                #print("err do",np.round(100*r,1),"%")
                #continue

                fig = make_syst_fig(mc_nom,mc_up_arr,mc_down_arr,syst,title=f"{var_name}_{yeartag}_{cat}_{syst}")

                out_path_for_this_cat = os.path.join(save_dir_path,os.path.join(yeartag,cat))
                if not os.path.exists(out_path_for_this_cat): os.makedirs(out_path_for_this_cat)
                fig.savefig(f"{out_path_for_this_cat}/{var_name}_{yeartag}_{cat}_{syst}.png")

            make_html(os.path.join(os.getcwd(),out_path_for_this_cat))


# A function for making a summary plot of SR yields
def make_sr_comb_plot(histo_dict,grouping_mc,grouping_data,year,ana_type="cb"):

    # Set variables based on cut based or bdt
    if ana_type == "cb":
        sr_lst  = sg.CAT_LST_CB
        hist_label = "Cut-based SRs"
        y_max = 9 # 9 is good for R2, 4 is good for R3
        fig_size = (12,7)
    elif ana_type == "bdt":
        sr_lst  = sg.CAT_LST_BDT
        hist_label = "BDT-based SRs"
        y_max = 20
        fig_size = (24,7)
    else:
        raise Exception("Unknown analysis type.")

    proc_lst  = sg.PROC_LST

    # Declare the hist we'll be filling
    histo_comb = hist.Hist(
        hist.axis.StrCategory(proc_lst, name="process", label="process"),
        hist.axis.StrCategory(sr_lst,   name="cat",     label=hist_label),
    )

    # Get the yield dict
    histo = histo_dict["njets"]
    sample_names_dict_mc   = sg.create_mc_sample_dict(year)
    sample_names_dict_data = sg.create_data_sample_dict(year)
    yld_dict_mc   = yt.get_yields(histo,sample_names_dict_mc)
    yld_dict_data = yt.get_yields(histo,sample_names_dict_data)

    # Apply the data-driven normalization for ZZ, ttZ
    yld_dict_mc, _, _ = yt.do_tf(yld_dict_mc,yld_dict_data,None,sg.BKG_TF_MAP,quiet=False)

    # Get the values and fill the combined hist
    histo = histo_dict["nleps"][{"systematic":"nominal"}]
    err_lst_p = []
    err_lst_m = []
    for cat_name in sr_lst:
        val_sum = 0
        var_sum = 0
        # Loop over processes, get yields to fill hist with (along with summing errs)
        for proc_name in proc_lst:
            val = yld_dict_mc[cat_name]["nominal"][proc_name][0]
            histo_comb[{"process": proc_name, "cat": cat_name}] = val
            # Sum the variances so we can have stat uncertainty (also need sum of nom)
            var_sum += yld_dict_mc[cat_name]["nominal"][proc_name][1]
            val_sum += val
        err_lst_p.append(val_sum + np.sqrt(var_sum))
        err_lst_m.append(val_sum - np.sqrt(var_sum))

    # Append a 0 err for overflow bin
    err_lst_p.append(0)
    err_lst_m.append(0)

    # Make plot
    fig = make_single_fig(histo_comb,err_p=err_lst_p,err_m=err_lst_m,ylims=[0,y_max],fig_size=fig_size)
    fig.savefig("sr_comb_plot.png")


# Main function for making CR plots
def make_plots(histo_dict,grouping_mc,grouping_data,save_dir_path,apply_nsf_to_cr=False):

    # Get the yield dict to get NSF dict (if we want to scale hists by NSFs)
    histo = histo_dict["njets"]
    yld_dict_mc   = yt.get_yields(histo,grouping_mc)
    yld_dict_data = yt.get_yields(histo,grouping_data)
    nsf_dict = yt.get_nsf_dict(yld_dict_mc,yld_dict_data,sg.BKG_TF_MAP)


    for var_name in histo_dict.keys():
        # Skip over some variables if you want to
        if "counts" in var_name: continue
        if var_name == "nbtagsm": continue
        #if var_name not in ["bdt_of_bin","bdt_sf_bin"]: continue
        #if var_name not in BDT_INPUT_LST: continue
        #if var_name not in BDT_SCORE_LST: continue

        # Get the relevant histogram from the input dict
        #print(f"\n{var_name}")
        histo_orig = histo_dict[var_name]

        # Just plot nominal syst for now
        histo_orig = histo_orig[{"systematic":"nominal"}]

        # Loop over categories and make plots for each
        for cat_name in histo_orig.axes["category"]:
            # Skip some of the cats if you want to
            #if "bdt" in cat_name: continue
            #if cat_name not in ["sr_4l_sf_incl", "sr_4l_of_incl", "cr_4l_btag_of", "cr_4l_btag_sf_offZ_met80", "cr_4l_sf", "sr_4l_bdt_sf_trn", "sr_4l_bdt_of_trn"]: continue # TMP
            if cat_name not in ["cr_4l_btag_of", "cr_4l_btag_sf_offZ_met80", "cr_4l_sf"]: continue
            #print(cat_name)

            # Make a copy so changes to binning do not propagate to next loop
            histo = copy.deepcopy(histo_orig)

            # Rebin if continous variable
            if var_name not in ["njets","nbtagsl","nleps","bdt_of_bin","bdt_sf_bin"]:
                if cat_name in ["cr_4l_btag_sf_offZ_met80","cr_4l_btag_of"]:
                    histo = rebin(histo,15)
                else:
                    histo = rebin(histo,6)

            histo_cat = histo[{"category":cat_name}]

            # Group the mc and data samples
            histo_grouped_mc = group(histo_cat,"process","process_grp",grouping_mc)
            histo_grouped_data = group(histo_cat,"process","process_grp",grouping_data)

            # Apply the NSF (the NSF dict is set up for SRs, not CRs but we can just grab the ones we need)
            if apply_nsf_to_cr:
                nsf_zz = 1
                nsf_ttz = 1
                if cat_name in ["cr_4l_btag_of"]:
                    nsfs = {"ZZ" : nsf_dict["sr_4l_of_1"]["ZZ"]["nsf"][0], "ttZ" : nsf_dict["sr_4l_of_1"]["ttZ"]["nsf"][0]}
                elif cat_name in ["cr_4l_btag_sf_offZ_met80", "cr_4l_sf"]:
                    nsfs = {"ZZ" : nsf_dict["sr_4l_sf_A"]["ZZ"]["nsf"][0], "ttZ" : nsf_dict["sr_4l_sf_A"]["ttZ"]["nsf"][0]}
                for i, name in enumerate(histo_grouped_mc.axes["process_grp"]):
                    # Scale the hist, see https://github.com/CoffeaTeam/coffea/discussions/705
                    histo_grouped_mc.view(flow=True)[i] *= nsfs.get(name,1) # Scale by 1 if the process is not ttZ or ZZ

            ######
            # Print stuff if you want to
            #if (cat_name == "cr_4l_sf" or cat_name == "cr_4l_btag_of" or cat_name=="cr_4l_btag_sf_offZ_met80") and var_name == "nleps":
                #print(f"\n{cat_name} {var_name}:")
                #print("Yields")
                #data = sum(sum(histo_grouped_data.values(flow=True)))
                #data_error = sum(sum(histo_grouped_data.variances(flow=True))) ** 0.5
                #mc = sum(sum(histo_grouped_mc.values(flow=True)))
                #mc_error = sum(sum(histo_grouped_mc.variances(flow=True))) ** 0.5
                #data_over_mc_ratio = data / mc
                #data_over_mc_ratio_error = data_over_mc_ratio * ((data_error / data)**2 + (mc_error / mc)**2)**0.5
                #print("mc:",mc)
                #print("mc Error:",mc_error)
                #print("data:",data)
                #print("data Error:",data_error)
                #print("data/mc:", data_over_mc_ratio)
                #print("data/mc Error:", data_over_mc_ratio_error)
            #continue
            #####

            # Merge overflow into last bin (so it shows up in the plot)
            histo_grouped_data = merge_overflow(histo_grouped_data)
            histo_grouped_mc = merge_overflow(histo_grouped_mc)

            # Make figure
            title = f"{cat_name}_{var_name}"
            print("Making: ",title)
            if "cr" in title:
                fig = make_cr_fig(histo_grouped_mc,histo_grouped_data,title=title)
            else:
                fig = make_cr_fig(histo_grouped_mc,title=title)

            # Save
            save_dir_path_cat = os.path.join(save_dir_path,cat_name)
            if not os.path.exists(save_dir_path_cat): os.mkdir(save_dir_path_cat)
            fig.savefig(os.path.join(save_dir_path_cat,title+".pdf"))
            fig.savefig(os.path.join(save_dir_path_cat,title+".png"))

            make_html(os.path.join(os.getcwd(),save_dir_path_cat))



###### Transfer factors for background ######

# TODO Get rid of this it's old
# Function for getting a dict with NSF and TF etc
def get_background_dict(yld_dict_mc,yld_dict_data,bkg_proc,cr_name,sr_name):

    # Get the sum of all other contributions
    bkg_all_but_bkg_of_interest = [0,0]
    for proc in yld_dict_mc.keys():
        if proc != bkg_proc:
            bkg_all_but_bkg_of_interest[0] += yld_dict_mc[proc][cr_name][0]
            bkg_all_but_bkg_of_interest[1] += yld_dict_mc[proc][cr_name][1]

    n_cr = yld_dict_data["data"][cr_name][0] - bkg_all_but_bkg_of_interest[0]
    m_cr = yld_dict_mc[bkg_proc][cr_name][0]
    m_sr = yld_dict_mc[bkg_proc][sr_name][0]

    n_cr_err = np.sqrt(yld_dict_data["data"][cr_name][1] + bkg_all_but_bkg_of_interest[1])
    m_cr_err = np.sqrt(yld_dict_mc[bkg_proc][cr_name][1])
    m_sr_err = np.sqrt(yld_dict_mc[bkg_proc][sr_name][1])

    out_dict = {
        "n_sr_est" : [n_cr*(m_sr/m_cr) , (n_cr*(m_sr/m_cr))*np.sqrt((n_cr_err/n_cr)**2 + (m_sr_err/m_sr)**2 + (m_cr_err/m_cr)**2)],
        "m_sr"     : [m_sr , m_sr_err],
        "n_cr"     : [n_cr , n_cr_err],
        "m_cr"     : [m_cr , m_cr_err],
        "tf"       : [m_sr/m_cr , (m_sr/m_cr)*np.sqrt((m_sr_err/m_sr)**2+(m_cr_err/m_cr)**2)],
        "nsf"      : [n_cr/m_cr , (n_cr/m_cr)*np.sqrt((n_cr_err/n_cr)**2+(m_cr_err/m_cr)**2)],
    }
    return out_dict

# CAUTION This function should probably not be used, just use the version in the datacard maker
# TODO move to same function as used in datacard maker
# Wrapper around the background estimation of TFs and yields
def do_background_estimation(yld_dict_mc,yld_dict_data,ul_year):

    # Map between short name and name to display in table
    kname_dict = {
        "n_sr_est" : "$N_{SR \\rm \\; est} = TF \\cdot N_{CR}$",
        "m_sr"     : "$MC_{SR}$",
        "n_cr"     : "$N_{CR}$",
        "m_cr"     : "$MC_{CR}$",
        "tf"       : "TF",
        "nsf"      : "NSF",
    }


    print_dict = {}

    # Do the ttZ and ZZ estimation for cut-based SRs
    print_dict["ttZ SR_OF"] = get_background_dict(yld_dict_mc,yld_dict_data,"ttZ","cr_4l_btag_of","sr_of_all")
    print_dict["ZZ SR_OF"]  = get_background_dict(yld_dict_mc,yld_dict_data,"ZZ","cr_4l_sf","sr_of_all")
    print_dict["ttZ SR_SF"] = get_background_dict(yld_dict_mc,yld_dict_data,"ttZ","cr_4l_btag_sf_offZ_met80","sr_sf_all")
    print_dict["ZZ SR_SF"]  = get_background_dict(yld_dict_mc,yld_dict_data,"ZZ","cr_4l_sf","sr_sf_all")

    # Do the ttZ and ZZ estimation for BDT SRs
    if (ul_year == "run2") or (("UL" in ul_year) and ("2022" not in ul_year)):
        for bdt_sr in SR_OF_BDT:
            print_dict[f"ttZ {bdt_sr}"] = get_background_dict(yld_dict_mc,yld_dict_data,"ttZ","cr_4l_btag_of",bdt_sr)
            print_dict[f"ZZ {bdt_sr}"] = get_background_dict(yld_dict_mc,yld_dict_data,"ZZ","cr_4l_sf",bdt_sr)
        for bdt_sr in SR_SF_BDT:
            print_dict[f"ttZ {bdt_sr}"] = get_background_dict(yld_dict_mc,yld_dict_data,"ttZ","cr_4l_btag_sf_offZ_met80",bdt_sr)
            print_dict[f"ZZ {bdt_sr}"] = get_background_dict(yld_dict_mc,yld_dict_data,"ZZ", "cr_4l_sf",bdt_sr)

    for k in print_dict:
        print(k,"\n",print_dict[k])

    # Print the dicts so we can look at the values
    # First replace key names with more descriptive names for printing
    for tf_dict_name in print_dict.keys():
        for kname in kname_dict.keys():
            print_dict[tf_dict_name][kname_dict[kname]] = print_dict[tf_dict_name].pop(kname)
    mlt.print_latex_yield_table(
        print_dict,
        tag="NSFs and TFs for ttZ and ZZ SR estimations",
        print_begin_info=True,
        print_end_info=True,
        roundat=3,
        print_errs=True,
        size="footnotesize",
        hz_line_lst=[1,3,3,19]
    )



################### Main ###################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default="plots", help = "The path the output files should be saved to")
    parser.add_argument('-y', "--get-yields", action='store_true', help = "Get yields from the pkl file")
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    parser.add_argument('-b', "--get-backgrounds", action='store_true', help = "Get background estimations")
    parser.add_argument('-u', "--ul-year", default='run2', help = "Which year to process", choices=["all","run2","run3","UL16APV","UL16","UL17","UL18","2022","2022EE"])
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))

    sample_dict_mc = sg.create_mc_sample_dict(args.ul_year)
    sample_dict_data = sg.create_data_sample_dict(args.ul_year)
    out_path = args.output_path


    # Wrapper around the code for getting the raw counts and dump to latex table
    #counts_dict = get_yields(histo_dict,sample_dict_mc,raw_counts=True)
    #print_counts(counts_dict)
    #exit()


    # Wrapper around the code for getting the TFs and background estimation factors
    # TODO move to same function as used in datacard maker
    if args.get_backgrounds:
        yld_dict_data = get_yields(histo_dict,sample_dict_data,quiet=True,blind=True)
        yld_dict_mc   = get_yields(histo_dict,sample_dict_mc,quiet=True)
        put_cat_col_sums(yld_dict_mc, sr_sf_lst=SR_SF_CB, sr_of_lst=SR_OF_CB)
        do_background_estimation(yld_dict_mc,yld_dict_data,args.ul_year)


    # Wrapper around the code for getting the yields for sr and bkg samples
    if args.get_yields:

        # Get the yield dict and put the extra columns and rows into it
        yld_dict = get_yields(histo_dict,sample_dict_mc)
        put_proc_row_sums(yld_dict, SR_SF_CB+SR_OF_CB)
        put_cat_col_sums(yld_dict, sr_sf_lst=SR_SF_CB, sr_of_lst=SR_OF_CB, tag="_cutbased")
        if (args.ul_year == "run2") or (("UL" in args.ul_year) and ("2022" not in args.ul_year)):
            put_proc_row_sums(yld_dict, SR_SF_BDT+SR_OF_BDT)
            put_cat_col_sums(yld_dict, sr_sf_lst=SR_SF_BDT, sr_of_lst=SR_OF_BDT, tag="_bdt")
        #print(yld_dict)
        #exit()

        # Dump latex table for cut based
        hlines = [2,3,7,8]
        sr_cats_to_print = SR_SF_CB + ["sr_sf_all_cutbased"] + SR_OF_CB + ["sr_of_all_cutbased","sr_all_cutbased"]
        #sr_cats_to_print = ["sr_sf_all_cutbased" , "sr_of_all_cutbased" , "sr_all_cutbased" , "sr_4l_sf_presel" , "sr_4l_sf_trn" , "sr_4l_of_presel"] # Preselection SR categories
        procs_to_print = ["WWZ","ZH","Sig","ZZ","ttZ","tWZ","WZ","other","Bkg",SOVERROOTB,SOVERROOTSPLUSB,"Zmetric"]
        print_yields(args.ul_year,yld_dict,sr_cats_to_print,procs_to_print,hlines=hlines,ref_dict=yd.EWK_REF) # Or e.g. for 2022 comp use yd.EWK_REF_2022

        # Dump latex table for BDT
        #hlines = [6,7,15,16]
        #sr_cats_to_print = SR_SF_BDT + ["sr_sf_all_bdt"] + SR_OF_BDT + ["sr_of_all_bdt","sr_all_bdt"]
        #procs_to_print = ["WWZ","ZH","Sig","ZZ","ttZ","tWZ","WZ","other","Bkg",SOVERROOTB,SOVERROOTSPLUSB,"Zmetric"]
        #print_yields(args.ul_year,yld_dict,sr_cats_to_print,procs_to_print,ref_dict=yd.EWK_REF,hlines=hlines)

        # Compare BDT yields against Keegan ref yields
        #keegan_ref = utils.put_none_errs(copy.deepcopy(yd.KEEGAN_BDT_YLDS))
        #sr_cats_to_print = SR_SF_BDT + SR_OF_BDT
        #procs_to_print = ["WWZ","ZH","ZZ","ttZ","tWZ","WZ","other"]
        #print_yields(args.ul_year,yld_dict,sr_cats_to_print,procs_to_print,ref_dict=keegan_ref,print_fom=False)


        # Dump yield dict to json
        json_name = "process_yields.json" # Could be an argument
        json_name = os.path.join(out_path,json_name)
        with open(json_name,"w") as out_file: json.dump(yld_dict, out_file, indent=4)
        print(f"\nSaved json file: {json_name}\n")


    # Make plots
    if args.make_plots:
        make_plots(histo_dict,sample_dict_mc,sample_dict_data,save_dir_path=out_path,apply_nsf_to_cr=False)
        #make_syst_plots(histo_dict,sample_dict_mc,sample_dict_data,out_path,args.ul_year) # Check on individual systematics
        #make_sr_comb_plot(histo_dict,sample_dict_mc,sample_dict_data,args.ul_year,ana_type="cb") # Make plot of all SR yields in one plot




if __name__ == "__main__":
    main()

