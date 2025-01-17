import argparse
import pickle
import gzip
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mp

from topcoffea.scripts.make_html import make_html
import ewkcoffea.modules.yield_tools as yt
import ewkcoffea.modules.sample_groupings as sg
import get_wwz_yields as gy


LABEL_MAP = {
    "mll_wl0_wl1"        : "$\mathrm{m_{\ell\ell}(\ell^W_0, \; \ell^W_1) \;[GeV]}$",
    "mllll"              : "$\mathrm{m_{\ell\ell\ell\ell} \;[GeV]}$",
    "absdphi_4l_met"     : "$|\Delta \phi \mathrm{(4\ell, \; p_T^{miss})}|$",
    "absdphi_zleps_met"  : "$|\Delta \phi \mathrm{(\ell^Z_0 + \ell^Z_1, \; p_T^{miss})}|$",
    "absdphi_wleps_met"  : "$|\Delta \phi \mathrm{(\ell^W_0 + \ell^W_1, \; p_T^{miss})}|$",
    "absdphi_wl0_met"    : "$|\Delta \phi \mathrm{(\ell^W_0, \; p_T^{miss})}|$",
    "absdphi_wl1_met"    : "$|\Delta \phi \mathrm{(\ell^W_1, \; p_T^{miss})}|$",
    "dr_wl0_wl1"         : "$\Delta \mathrm{R (\ell^W_0, \; \ell^W_1)}$",
    "dr_zl0_zl1"         : "$\Delta \mathrm{R (\ell^Z_0, \; \ell^Z_1)}$",
    "dr_wleps_zleps"     : "$\Delta \mathrm{R (\ell^Z_0 + \ell^Z_1, \; \ell^W_0 + \ell^W_1)}$",
    "met"                : "$\mathrm{p_T^{miss} \;[GeV]}$",
    "mt2"                : "$\mathrm{m_{T2} \;[GeV]}$",
    "ptl4"               : "$\mathrm{p_T^{4\ell} \;[GeV]}$",
    "scalarptsum_jet"    : "$\Sigma \mathrm{p_T^{j} \;[GeV]}$",
    "scalarptsum_lepmet" : "$\Sigma \mathrm{p_T^\ell + p_T^{miss} \;[GeV]}$",
    "z_lep0_pt"          : "$\mathrm{p_T({\ell^Z_0}) \;[GeV]}$",
    "z_lep1_pt"          : "$\mathrm{p_T({\ell^Z_1}) \;[GeV]}$",
    "w_lep0_pt"          : "$\mathrm{p_T({\ell^W_0}) \;[GeV]}$",
    "w_lep1_pt"          : "$\mathrm{p_T({\ell^W_1}) \;[GeV]}$",
    "njets"              : "$\mathrm{n_{j}}$",
    "cos_helicity_x"     : "$cos \; helicity \; X$",
    "mt_wl0_met"         : "$\mathrm{m_T(\ell^W_0, p_T^{miss}) \;[GeV]}$",
    "mt_wl1_met"         : "$\mathrm{m_T(\ell^W_1, p_T^{miss}) \;[GeV]}$",
    "mt_wleps_met"       : "$\mathrm{m_T(\ell^W_0 + \ell^W_1, p_T^{miss}) \;[GeV]}$",
    "mt_4l_met"          : "$\mathrm{m_T(4\ell, p_T^{miss}) \;[GeV]}$",
    "dr_wl0_j_min"       : "$\Delta \mathrm{R ( \ell^W_0, j)^{min}}$",
    "dr_wl1_j_min"       : "$\Delta \mathrm{R ( \ell^W_1, j)^{min}}$",

    "bdt_of_wwz"         : "$\mathrm{BDT \; Score} \; _{\mathrm{WWZ}}$",
    "bdt_of_zh"          : "$\mathrm{BDT \; Score} \; _{\mathrm{ZH}}$",
    "bdt_of_bkg"         : "$\mathrm{BDT \; Score} \; _{\mathrm{Background}}$",
    "bdt_of_wwz_m_zh"    : "$\mathrm{BDT \; Score} \; _{\mathrm{WWZ}} \; - \; \mathrm{BDT \; score} \; _{\mathrm{ZH}}$",

    "bdt_sf_wwz"         : "$\mathrm{BDT \; Score} \; _{\mathrm{WWZ}}$",
    "bdt_sf_zh"          : "$\mathrm{BDT \; Score} \; _{\mathrm{ZH}}$",
    "bdt_sf_bkg"         : "$\mathrm{BDT \; Score} \; _{\mathrm{Background}}$",
    "bdt_sf_wwz_m_zh"    : "$\mathrm{BDT \; Score} \; _{\mathrm{WWZ}} \; - \; \mathrm{BDT \; score} \; _{\mathrm{ZH}}$",

    "nbtagsl"            : "$\mathrm{n_{b}}$",
}

STYLE_DICT = {

    # Input vars in OF SR
    "input_vars_of" : {
        "cats_of_interest" : ["sr_4l_bdt_of_trn", "cr_4l_btag_of"],
        "rebin" : {"run2": 18, "run3" : 30},
        "var_dict" : {
            "mll_wl0_wl1" : {
            },
            "mllll" : {
            },
            "absdphi_4l_met" : {
            },
            "absdphi_wleps_met" : {
            },
            "absdphi_wl0_met" : {
            },
            "absdphi_wl1_met" : {
            },
            "dr_wl0_wl1" : {
            },
            "dr_zl0_zl1" : {
            },
            "dr_wleps_zleps" : {
            },
            "met" : {
            },
            "mt2" : {
            },
            "ptl4" : {
            },
            "scalarptsum_jet" : {
            },
            "scalarptsum_lepmet" : {
            },
            "z_lep0_pt" : {
            },
            "z_lep1_pt" : {
            },
            "w_lep0_pt" : {
            },
            "w_lep1_pt" : {
            },
            "njets" : {
                "rebin":{"run2":None, "run3":None}
            },
            "cos_helicity_x" : {
            },
            "mt_wl0_met" : {
            },
            "mt_wl1_met" : {
            },
            "mt_wleps_met" : {
            },
            "mt_4l_met" : {
            },
            "dr_wl0_j_min" : {
            },
            "dr_wl1_j_min" : {
            },
        },
    },

    # Input vars in SF SR
    "input_vars_sf" : {
        "cats_of_interest" : ["sr_4l_bdt_sf_trn", "cr_4l_sf", "cr_4l_btag_sf_offZ_met80"],
        "rebin" : {"run2": 18, "run3" : 30},
        "var_dict" : {
            "mll_wl0_wl1" : {
            },
            "mllll" : {
            },
            "absdphi_4l_met" : {
            },
            "absdphi_zleps_met" : {
            },
            "absdphi_wleps_met" : {
            },
            "absdphi_wl0_met" : {
            },
            "absdphi_wl1_met" : {
            },
            "dr_wl0_wl1" : {
            },
            "dr_zl0_zl1" : {
            },
            "dr_wleps_zleps" : {
            },
            "met" : {
            },
            "mt2" : {
            },
            "ptl4" : {
            },
            "scalarptsum_jet" : {
            },
            "scalarptsum_lepmet" : {
            },
            "z_lep0_pt" : {
            },
            "z_lep1_pt" : {
            },
            "w_lep0_pt" : {
            },
            "w_lep1_pt" : {
            },
            "njets" : {
                "rebin": {"run2":None, "run3":None}
            },
            "cos_helicity_x" : {
            },
            "mt_wl0_met" : {
            },
            "mt_wl1_met" : {
            },
            "mt_wleps_met" : {
            },
            "mt_4l_met" : {
            },
            "dr_wl0_j_min" : {
            },
            "dr_wl1_j_min" : {
            },
        },
    },

    # BDT scores in OF SR
    "scores_of" : {
        "cats_of_interest" : ["sr_4l_bdt_of_trn", "cr_4l_btag_of"],
        "rebin" : {"run2": 18, "run3" : 30},
        "var_dict" : {
            "bdt_of_wwz" : {
            },
            "bdt_of_zh" : {
            },
            "bdt_of_bkg" : {
            },
            "bdt_of_wwz_m_zh" : {
            },
        },
    },

    # BDT scores in OF SR
    "scores_sf" : {
        "cats_of_interest" : ["sr_4l_bdt_sf_trn", "cr_4l_sf", "cr_4l_btag_sf_offZ_met80"],
        "rebin" : {"run2": 18, "run3" : 30},
        "var_dict" : {
            "bdt_sf_wwz" : {
            },
            "bdt_sf_zh" : {
            },
            "bdt_sf_bkg" : {
            },
            "bdt_sf_wwz_m_zh" : {
            },
        },
    },

    # Variables we cut on shown in pseudo preselection regions
    "kinematic_cut_vars_of" : {
        "cats_of_interest" : ["sr_4l_bdt_of_presel_nobreq"],
        "rebin" : {"run2": 6, "run3" : 12},
        "var_dict" : {
            "nbtagsl" : {
            },
        },
    },
    "kinematic_cut_vars_sf" : {
        "cats_of_interest" : ["sr_4l_bdt_sf_presel"],
        "rebin" : {"run2": 12, "run3" : 18},
        "var_dict" : {
            "mt2" : {
                "logscale" : True,
                "rangex" : [0,100],
            },
        },
    },
}


# Takes a mc hist and data hist and plots both
def make_public_fig(histo_mc,histo_data=None,title="test",unit_norm_bool=False,axisrangex=None,xlabel=None,year="run2",logscale=False):

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
    ax.fill_between(bin_edges_arr,err_m,err_p, step='post', facecolor='none', edgecolor='gray', alpha=0.5, linewidth=0.0, label='Stat', hatch='/////')

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
        # TODO: NOTE This does not touch overlflow, simply cuts the axis range that is displayed
        ax.set_xlim(axisrangex[0],axisrangex[1])
        rax.set_xlim(axisrangex[0],axisrangex[1])

    # CMS text
    plt.text(0,1.02,"CMS",fontsize=23,weight="bold",transform=ax.transAxes)
    plt.text(0.15,1.02,"$\it{Supplementary}$",fontsize=19,transform=ax.transAxes)
    if year == "run2":
        extt = plt.text(0.59,1.02,"138 $\mathrm{fb^{{-}1}}$ (13 TeV)",fontsize=18,transform=ax.transAxes)
    elif year == "run3":
        extt = plt.text(0.57,1.02,"62 $\mathrm{fb^{{-}1}}$ (13.6 TeV)",fontsize=18,transform=ax.transAxes)

    # Internal legend
    extr = ax.legend(loc="upper left",bbox_to_anchor=(1,1),fontsize="16", frameon=False)

    # Set style things on main plot
    ax.autoscale(axis='y')
    ax.set_xlabel(None)
    ax.tick_params(axis='y', labelsize=16)
    extl = ax.set_ylabel('Events',fontsize=22,loc="top")

    # Set style things on ratio plot
    if xlabel is not None:
        extb = rax.set_xlabel(xlabel,fontsize=18,loc="right")
    rax.set_ylabel('Data/Pred.',fontsize=18)
    rax.set_ylim(0.0,2.0)
    rax.axhline(1.0,linestyle="-",color="k",linewidth=1)
    rax.tick_params(axis='x', labelsize=16)
    #rax.xaxis.set_label_coords(0.82, -0.40)
    #rax.yaxis.set_label_coords(-0.09, 0.5)

    # No more than 5 main ticks for the x axis
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    rax.xaxis.set_major_locator(plt.MaxNLocator(5))

    if logscale:
        ax.set_yscale('log')

    return (fig,(extt,extr,extb,extl))


# Make a standalone legend
def make_legend():
    lfig, lax = plt.subplots(figsize=(1.5,2.5))
    leg_elements = [
        mp.Patch(facecolor=gy.CLR_LST[0], lw=9, label="WWZ"),
        mp.Patch(facecolor=gy.CLR_LST[1], lw=9, label="ZH"),
        mp.Patch(facecolor=gy.CLR_LST[2], lw=9, label="ZZ"),
        mp.Patch(facecolor=gy.CLR_LST[3], lw=9, label="ttZ"),
        mp.Patch(facecolor=gy.CLR_LST[4], lw=9, label="tWZ"),
        mp.Patch(facecolor=gy.CLR_LST[5], lw=9, label="WZ"),
        mp.Patch(facecolor=gy.CLR_LST[6], lw=9, label="Other"),
        mp.Patch(facecolor="none",edgecolor="gray",hatch="/////", alpha=0.5, label="Stat"),
        Line2D([0],[0],color="black", marker="o", label="Data"),
    ]
    lax.axis('off')
    lax.legend(handles=leg_elements,fontsize="12",loc="center",frameon=False)
    return lfig


# Main function for making CR plots
def make_plots(histo_dict,grouping_mc,grouping_data,out_dir_map,save_dir_path,year="run2",max_nplots=9999):

    # Set up the output dir if it does not exist
    save_dir_path_year = os.path.join(save_dir_path,year)
    if not os.path.exists(save_dir_path): os.mkdir(save_dir_path)
    if not os.path.exists(save_dir_path_year): os.mkdir(save_dir_path_year)

    # Make a standalone legend
    lfig = make_legend()
    lfig.savefig(os.path.join(save_dir_path_year,"legend.pdf"))

    ### Loop over the groups of plots ###
    for group_name in STYLE_DICT.keys():

        ### Loop over the relevant regions ###
        for cat_name in STYLE_DICT[group_name]["cats_of_interest"]:

            # Make a subdir for this set
            set_tag = f"{group_name}_{cat_name}"
            if set_tag not in out_dir_map: continue
            else:
                save_dir_path_year_set = os.path.join(save_dir_path_year,out_dir_map[set_tag])
                if not os.path.exists(save_dir_path_year_set): os.mkdir(save_dir_path_year_set)

            ### Loop over the variabls and make plots ###
            n_plotted = 0
            for var_name in STYLE_DICT[group_name]["var_dict"].keys():
                print(f"\nVar name: {var_name}")

                histo_cat = histo_dict[var_name][{"systematic":"nominal", "category":cat_name}]

                # Rebin and set x range
                rangex = None
                if "rangex" in STYLE_DICT[group_name]["var_dict"][var_name]:
                    rangex = STYLE_DICT[group_name]["var_dict"][var_name]["rangex"]
                if var_name not in ["njets", "nbtagsl"]:
                    rebin_factor = STYLE_DICT[group_name]["rebin"][year]
                    if "rebin" in STYLE_DICT[group_name]["var_dict"][var_name]:
                        rebin_factor = STYLE_DICT[group_name]["var_dict"][var_name]["rebin"][year]
                    histo_cat = gy.rebin(histo_cat,rebin_factor)

                logscale=False
                if "logscale" in STYLE_DICT[group_name]["var_dict"][var_name]:
                    logscale = STYLE_DICT[group_name]["var_dict"][var_name]["logscale"]

                # Group the mc and data samples
                histo_grouped_mc = gy.group(histo_cat,"process","process_grp",grouping_mc)
                histo_grouped_data = gy.group(histo_cat,"process","process_grp",grouping_data)

                # Merge overflow into last bin (so it shows up in the plot)
                histo_grouped_data = gy.merge_overflow(histo_grouped_data)
                histo_grouped_mc = gy.merge_overflow(histo_grouped_mc)

                # Make and save figure
                title = f"{group_name}_{var_name}"
                print("Making: ",title)
                fig, ext_tup = make_public_fig(histo_grouped_mc,histo_grouped_data,title=title,xlabel=LABEL_MAP[var_name],year=year,logscale=logscale)
                fig.savefig(os.path.join(save_dir_path_year_set,title+".pdf"),bbox_extra_artists=ext_tup,bbox_inches='tight')
                fig.savefig(os.path.join(save_dir_path_year_set,title+".png"),bbox_extra_artists=ext_tup,bbox_inches='tight')

                # Break if we've hit the max we'd like to plot (mainly for testing small changes)
                n_plotted = n_plotted+1
                if n_plotted >= max_nplots: break

            # Make html for this sub dir
            make_html(os.path.join(os.getcwd(),save_dir_path_year_set),width=445, height=355)



################### Main ###################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default="plots_public", help = "The path the output files should be saved to")
    parser.add_argument('-u', "--ul-year", default='run2', help = "Which year to process", choices=["run2","run3"])
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))

    sample_dict_mc = sg.create_mc_sample_dict(args.ul_year)
    sample_dict_data = sg.create_data_sample_dict(args.ul_year)
    out_path = args.output_path

    # This defines what the output dirs will be called
    # To skip plotting one of these groups, just comment it
    max_n_to_plot = 1 # Plot only a limited number per set (for testing small subsets)
    dir_name_map = {

        # Input variables (in SRs and CRs)
        "input_vars_of_sr_4l_bdt_of_trn"         : "input_vars_SR_OF",
        "input_vars_sf_sr_4l_bdt_sf_trn"         : "input_vars_SR_SF",
        "input_vars_sf_cr_4l_sf"                 : "input_vars_CR_ZZ",
        "input_vars_of_cr_4l_btag_of"            : "input_vars_CR_ttZ_OF",
        "input_vars_sf_cr_4l_btag_sf_offZ_met80" : "input_vars_CR_ttZ_SF",

        # BDT score distributions (in SRs and CRs)
        "scores_of_sr_4l_bdt_of_trn"             : "bdt_scores_SR_OF",
        "scores_sf_sr_4l_bdt_sf_trn"             : "bdt_scores_SR_SF",
        "scores_sf_cr_4l_sf"                     : "bdt_scores_CR_ZZ",
        "scores_of_cr_4l_btag_of"                : "bdt_scores_CR_ttZ_OF",
        "scores_sf_cr_4l_btag_sf_offZ_met80"     : "bdt_scores_CR_ttZ_SF",

        # Plots to motivate kinematic cuts
        "kinematic_cut_vars_of_sr_4l_bdt_of_presel_nobreq" : "kinematic_cut_vars_OF_nbtag",
        "kinematic_cut_vars_sf_sr_4l_bdt_sf_presel"        : "kinematic_cut_vars_SF_mt2",
    }

    # Make the plots
    make_plots(histo_dict,sample_dict_mc,sample_dict_data,dir_name_map,save_dir_path=out_path,year=args.ul_year,max_nplots=max_n_to_plot)


if __name__ == "__main__":
    main()

