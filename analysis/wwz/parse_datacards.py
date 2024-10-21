import numpy as np
import os
import argparse

# Check average size of systs

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


# Read the file
def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content

# Turn ["syst","lnN","a/b","c/d"] into [e,f] where e is symmetrized avg of a and b, similar for f
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

######################
# Messy
# Get average sizes
def get_sizes(proc_lst,rate_lst,tag,systs_to_group,syst_dict):

    # Get nominal
    yld_nom = 0
    for i,proc in enumerate(proc_lst):
        rate = rate_lst[i]
        yld_nom = yld_nom + rate

    tot_up_vars = []
    for syst_name in systs_to_group:
        syst_row = syst_dict[syst_name]
        yld_up = 0
        for i,proc in enumerate(proc_lst):
            yld_up  = yld_up  + rate_lst[i]*syst_row[i]
        tot_up_vars.append(yld_up)

    # Sum the relative errors in quad
    quad_sum_sq = 0
    for up_var in tot_up_vars:
        rel_err = up_var/yld_nom
        quad_sum_sq = quad_sum_sq + (1-rel_err)**2
    quad_sum = quad_sum_sq**0.5

    ## Average
    ##avg = (sum(np.array(tot_up_vars)))/len(tot_up_vars)
    ##print("tot_up_vars:",tot_up_vars)
    ##print("nom:",yld_nom)
    ##print("tot_up_vars/nom:",np.array(tot_up_vars)/yld_nom)
    ##print("quad",quad_sum,"->",quad_sum*100,"%")

    percent = 100*(np.array(tot_up_vars)/yld_nom - 1)
    avg = sum(percent)/len(percent)
    #print(systs_to_group,"tot_up_vars/nom:",np.array(tot_up_vars)/yld_nom,percent,"%")
    #print("\n",systs_to_group,"\n",round(avg,2),"%")
    print(f"  {tag}: {round(avg,2)} %")

    return avg


def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("base_path", help = "The path to the cards")
    parser.add_argument('-u', "--ul-year", default='run2', help = "Which year to process", choices=["run2","run3"])
    args = parser.parse_args()

    run_name = args.ul_year
    base_path = args.base_path
    syst_grp_dict = SYST_GRP[run_name]

    # Get the lines fsrom the in files
    dc_names_dict = {
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


    out_dict = {}
    for dc_name in dc_names_dict[run_name]:
        lines = read_file(os.path.join(base_path,dc_name))

        dc_tag = dc_name.split("/")[-1][:-4]
        out_dict[dc_tag] = {}

        # Loop over lines in card to get procs and rates and systs
        syst_dict = {}
        rate_lst = []
        for line in lines:
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


        # Get the syst vars
        print("\nCheck avg sizes for:",dc_tag)
        for systs_tag,systs_to_group in syst_grp_dict.items():
            out_dict[dc_tag][systs_tag] = get_sizes(proc_lst,rate_lst,systs_tag,systs_to_group,syst_dict)

    # After
    print(out_dict)

    for grp in syst_grp_dict:
        print("\n",grp)
        tmp_lst = []
        for card_name in out_dict.keys():
            tmp_lst.append(out_dict[card_name][grp])
        avg = sum(np.array(tmp_lst))/len(tmp_lst)
        #print(tmp_lst,avg)
        print(avg)


main()

