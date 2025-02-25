import numpy as np
import os
import argparse

DC_NAMES_DICT = {
    "run2" : [
        "wwz4l_card_sr_4l_bdt_of_1_run2.txt",
        "wwz4l_card_sr_4l_bdt_of_2_run2.txt",
        "wwz4l_card_sr_4l_bdt_of_3_run2.txt",
        "wwz4l_card_sr_4l_bdt_of_4_run2.txt",
        "wwz4l_card_sr_4l_bdt_of_5_run2.txt",
        "wwz4l_card_sr_4l_bdt_of_6_run2.txt",
        "wwz4l_card_sr_4l_bdt_of_7_run2.txt",
        "wwz4l_card_sr_4l_bdt_of_8_run2.txt",
        "wwz4l_card_sr_4l_bdt_sf_1_run2.txt",
        "wwz4l_card_sr_4l_bdt_sf_2_run2.txt",
        "wwz4l_card_sr_4l_bdt_sf_3_run2.txt",
        "wwz4l_card_sr_4l_bdt_sf_4_run2.txt",
        "wwz4l_card_sr_4l_bdt_sf_5_run2.txt",
        "wwz4l_card_sr_4l_bdt_sf_6_run2.txt",
        "wwz4l_card_sr_4l_bdt_sf_7_run2.txt",
        "wwz4l_card_sr_4l_bdt_sf_8_run2.txt",
    ],
    "run3" : [
        "wwz4l_card_sr_4l_bdt_of_coarse_1_run3.txt",
        "wwz4l_card_sr_4l_bdt_of_coarse_2_run3.txt",
        "wwz4l_card_sr_4l_bdt_of_coarse_3_run3.txt",
        "wwz4l_card_sr_4l_bdt_of_coarse_4_run3.txt",
        "wwz4l_card_sr_4l_bdt_sf_coarse_1_run3.txt",
        "wwz4l_card_sr_4l_bdt_sf_coarse_2_run3.txt",
        "wwz4l_card_sr_4l_bdt_sf_coarse_3_run3.txt",
        "wwz4l_card_sr_4l_bdt_sf_coarse_4_run3.txt",
    ]
}

SYST_GRP = {

    "run2" : {
        "btag" : [
            "CMS_btag_fixedWP_comb_bc_correlated",
            "CMS_btag_fixedWP_comb_bc_uncorrelated_2016postVFP",
            "CMS_btag_fixedWP_comb_bc_uncorrelated_2016preVFP",
            "CMS_btag_fixedWP_comb_bc_uncorrelated_2017",
            "CMS_btag_fixedWP_comb_bc_uncorrelated_2018",
            "CMS_btag_fixedWP_incl_light_correlated",
            "CMS_btag_fixedWP_incl_light_uncorrelated_2016postVFP",
            "CMS_btag_fixedWP_incl_light_uncorrelated_2016preVFP",
            "CMS_btag_fixedWP_incl_light_uncorrelated_2017",
            "CMS_btag_fixedWP_incl_light_uncorrelated_2018",
        ],
        "ele" : ["CMS_eff_e_13TeV",],
        "muo" : ["CMS_eff_m_13TeV",],
        "prefire" : ["CMS_l1_ecal_prefiring",],
        "pu" : ["CMS_pileup",],
        "jerc" : [
            "CMS_res_j_2016postVFP",
            "CMS_res_j_2016preVFP",
            "CMS_res_j_2017",
            "CMS_res_j_2018",
            "CMS_scale_j_2016postVFP",
            "CMS_scale_j_2016preVFP",
            "CMS_scale_j_2017",
            "CMS_scale_j_2018",
        ],
        "met" : [
            "CMS_scale_met_unclustered_energy_2016postVFP",
            "CMS_scale_met_unclustered_energy_2016preVFP",
            "CMS_scale_met_unclustered_energy_2017",
            "CMS_scale_met_unclustered_energy_2018",
        ],
        "renormfact" : [
            "QCDscale_fac_WWZ",
            "QCDscale_fac_WZ",
            "QCDscale_fac_ZH",
            "QCDscale_fac_ZZ",
            "QCDscale_fac_other",
            "QCDscale_fac_tWZ",
            "QCDscale_fac_ttZ",
            "QCDscale_ren_WWZ",
            "QCDscale_ren_WZ",
            "QCDscale_ren_ZH",
            "QCDscale_ren_ZZ",
            "QCDscale_ren_other",
            "QCDscale_ren_tWZ",
            "QCDscale_ren_ttZ",
        ],
        "ps" : [
            "ps_fsr_WWZ", "ps_fsr_ZH", "ps_fsr_ZZ", "ps_fsr_ttZ", "ps_fsr_tWZ", "ps_fsr_WZ", "ps_fsr_other",
            "ps_isr_WWZ", "ps_isr_ZH", "ps_isr_ZZ", "ps_isr_ttZ", "ps_isr_tWZ", "ps_isr_WZ", "ps_isr_other",
        ],

    },

    "run3" : {

        "btag" : [
            "CMS_btag_fixedWP_comb_bc_correlated",
            "CMS_btag_fixedWP_comb_bc_uncorrelated_2022",
            "CMS_btag_fixedWP_comb_bc_uncorrelated_2022EE",
            "CMS_btag_fixedWP_comb_bc_uncorrelated_2023",
            "CMS_btag_fixedWP_comb_bc_uncorrelated_2023BPix",
            "CMS_btag_fixedWP_incl_light_correlated",
        ],
        "ele" : ["CMS_eff_e_13p6TeV",],
        "muo" : ["CMS_eff_m_13p6TeV",],
        "pu" : ["CMS_pileup",],
        "jerc" : [
            "CMS_res_j_2022",
            "CMS_res_j_2022EE",
            "CMS_res_j_2023",
            "CMS_res_j_2023BPix",
            "CMS_scale_j_2022",
            "CMS_scale_j_2022EE",
            "CMS_scale_j_2023",
            "CMS_scale_j_2023BPix",
        ],
        "met" : [
            "CMS_scale_met_unclustered_energy_2022",
            "CMS_scale_met_unclustered_energy_2022EE",
            "CMS_scale_met_unclustered_energy_2023",
            "CMS_scale_met_unclustered_energy_2023BPix",
        ],

        "renormfact" : [
            "QCDscale_fac_WWZ",
            "QCDscale_fac_WZ",
            "QCDscale_fac_ZH",
            "QCDscale_fac_ZZ",
            "QCDscale_fac_other",
            "QCDscale_fac_tWZ",
            "QCDscale_fac_ttZ",
            "QCDscale_ren_WWZ",
            "QCDscale_ren_WZ",
            "QCDscale_ren_ZH",
            "QCDscale_ren_ZZ",
            "QCDscale_ren_other",
            "QCDscale_ren_tWZ",
            "QCDscale_ren_ttZ",
        ],
        "ps" : [
            "ps_fsr_WWZ", "ps_fsr_ZH", "ps_fsr_ZZ", "ps_fsr_ttZ", "ps_fsr_tWZ", "ps_fsr_WZ", "ps_fsr_other",
            "ps_isr_WWZ", "ps_isr_ZH", "ps_isr_ZZ", "ps_isr_ttZ", "ps_isr_tWZ", "ps_isr_WZ", "ps_isr_other",
        ],
    }

}


# Read a text file
def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


# Manipulate syst row in datacard to get list of symmetrized variations, more specifically:
#     - Turn ["syst","lnN","a/b","c/d"] into [e,f]
#     - Where e is symmetrized avg of a and b, similar for f
def reformat_syst_lin(in_line):
    out_line = []
    syst_name = in_line[0]
    vals = in_line[2:]
    for kappa in vals:

        # We have just one kappa val
        if len(kappa.split("/")) == 1:
            # The kappa val is NA
            if kappa== "-":
                out_line.append(1.0) # TODO how to handle later
            else:
                if not float(kappa) >= 1: raise Exception # Don't think this should happen
                out_line.append(float(kappa))

        # Otherwise we have an up down kappa pair
        else:
            do = float(kappa.split("/")[0])
            up = float(kappa.split("/")[1])
            av = 1 + (abs(1-do) + abs(1-up))/2
            out_line.append(av)

    return out_line


# From the list of lines in the datacard, get the proc list and rate list and dict of syst lists
def get_proc_rate_syst(dc_lines):
    # Loop over lines in card to get procs and rates and systs
    syst_dict = {}
    rate_lst = []
    for line in dc_lines:
        line_split = line.split()

        # Proc line
        if (line_split[0] == "process") and "-1" not in line_split:
            proc_lst = line_split[1:]

        # Rate line
        if (line_split[0] == "rate"):
            for rate in line_split[1:]:
                rate_lst.append(float(rate))

        # Get the syst lines
        if "lnN" in line_split:
            syst_dict[line_split[0]] = reformat_syst_lin(line_split)

    # Sanity check of lenghts
    for k,v in syst_dict.items():
        if len(v) != len(proc_lst): raise Exception("Wrong len proc")
        if len(v) != len(rate_lst): raise Exception("Wrong len rate")

    return [proc_lst,rate_lst,syst_dict]



# Get average sizes of systematics
#     - Inputs are for a single card
#     - Inputs are process list, rate list, variation list for each syst
#     - Also input the info about which systematics to average together
#     - Averages the sizes of the syst effects for all systs in the group
def get_avg_sizes(proc_lst,rate_lst,syst_dict,syst_group_tag,syst_names_to_group,verbose=True):

    # Get nominal yield (sum over processes)
    yld_nom = 0
    for i,proc in enumerate(proc_lst):
        rate = rate_lst[i]
        yld_nom = yld_nom + rate

    tot_up_vars = []
    for syst_name in syst_names_to_group:
        syst_row = syst_dict[syst_name]
        yld_up = 0
        for i,proc in enumerate(proc_lst):
            yld_up  = yld_up  + rate_lst[i]*syst_row[i]
        tot_up_vars.append(yld_up)

    # Get average and print stuff

    ### For quad sum comparison ###
    '''
    print("syst_group_tag", syst_group_tag)
    x = 0
    for up_var in tot_up_vars:
        x = x + (up_var/yld_nom - 1)**2
    x = 100*(x**0.5)
    #avg = x
    '''
    ###


    percent = 100*(np.array(tot_up_vars)/yld_nom - 1)
    avg = sum(percent)/len(percent)
    if verbose:
        print(f"  {syst_group_tag}: {round(avg,2)} %")
        #print(f"  {syst_group_tag}: {avg} %")

    return avg



#########################################
############# Main function #############

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("base_path", help = "The path to the cards")
    parser.add_argument('-u', "--ul-year", default='run2', help = "Which year to process", choices=["run2","run3"])
    args = parser.parse_args()
    run_name = args.ul_year
    base_path = args.base_path
    #base_path = "/home/users/kmohrman/public_html/WWZ/smp24015_ref/ref_datacards/analysis_test_02_v00_oct06/bdt_withSys"

    syst_grp_dict = SYST_GRP[run_name]
    dc_names_lst = DC_NAMES_DICT[run_name]


    # Loop over datacards in the R2 or R3 set, to get the info per card
    syst_summary_dict = {}
    for dc_name in dc_names_lst:

        # Get lines from this datacard, set up entry in summary dict for it
        lines = read_file(os.path.join(base_path,dc_name))
        dc_tag = dc_name[:-4] # Drop the .txt from dc name
        syst_summary_dict[dc_tag] = {}

        # Get the proc list, the rate list, and the list for each syst
        proc_lst, rate_lst, syst_lst_dict = get_proc_rate_syst(lines)

        # Get the average syst variation
        print("\nCheck avg sizes for:",dc_tag)
        for systs_tag,syst_names_to_group in syst_grp_dict.items():
            syst_summary_dict[dc_tag][systs_tag] = get_avg_sizes(proc_lst,rate_lst,syst_lst_dict,systs_tag,syst_names_to_group)


    # Print meta averages for each group across all datacards
    print(syst_summary_dict)
    for grp in syst_grp_dict:
        tmp_lst = []
        for card_name in syst_summary_dict.keys():
            tmp_lst.append(syst_summary_dict[card_name][grp])
        avg = sum(np.array(tmp_lst))/len(tmp_lst)
        print(f"\n{grp}")
        print(avg)


main()

