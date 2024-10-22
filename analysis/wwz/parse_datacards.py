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
        "pu" : ['PU'],
        "prefire" : ['PreFiring'],
        "scale" : ['renorm', 'fact'],
        "ps": ['ISR', 'FSR',],
        "btag" : [
            'btagSFbc_correlated',
            'btagSFlight_correlated',
            'btagSFbc_uncorrelated_2018',
            'btagSFlight_uncorrelated_2018',
            'btagSFbc_uncorrelated_2016APV',
            'btagSFlight_uncorrelated_2016APV',
            'btagSFbc_uncorrelated_2017',
            'btagSFlight_uncorrelated_2017',
            'btagSFbc_uncorrelated_2016',
            'btagSFlight_uncorrelated_2016',
        ],
        "ele" : ['lepSF_elec_run2'],
        "mu" : ['lepSF_muon_run2'],
        "jerc" : [
            'JEC_2018',
            'JER_2018',
            'JEC_2017',
            'JER_2017',
            'JEC_2016',
            'JER_2016',
            'JEC_2016APV',
            'JER_2016APV',
        ],
        #['lumi'],
        #['theory_norm_other_other'],
        #['fake_WZ_run2']
    },

    "run3" : {
        "pu" : ['PU'],
        "scale" : ['renorm', 'fact'],
        "ps": ['ISR', 'FSR',],
        "btag" : [
            'btagSFbc_correlated',
            'btagSFlight_correlated',
            'btagSFbc_uncorrelated_2022',
            'btagSFbc_uncorrelated_2022EE',
            'btagSFbc_uncorrelated_2023',
            'btagSFbc_uncorrelated_2023BPix',
        ],
        "ele" : ['lepSF_elec_run3'],
        "mu" : ['lepSF_muon_run3'],
        "jerc" : [
            'JEC_2022',
            'JER_2022',
            'JEC_2022EE',
            'JER_2022EE',
            'JEC_2023',
            'JER_2023',
            'JEC_2023BPix',
            'JER_2023BPix',
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
    percent = 100*(np.array(tot_up_vars)/yld_nom - 1)
    avg = sum(percent)/len(percent)
    if verbose:
        print(f"  {syst_group_tag}: {round(avg,2)} %")

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

