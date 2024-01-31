import os
import json
import pickle
import gzip
import argparse

from topcoffea.modules import utils
import get_wwz_yields as gy # Note the fact that we're using functions from here means they probably belongs in ewkcoffea/ewkcoffea/modules


# Global variables
PRECISION = 6   # Decimal point precision in the text datacard output
PROC_LST = ["WWZ","ZH","ZZ","ttZ","tWZ","other"]
SIG_LST = ["WWZ","ZH"]
CAT_LST = ["sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]

# Systs that are not correlated across years
SYSTS_SPECIAL = {
    "btagSFlight_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
    "btagSFbc_uncorrelated_2016APV"    : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
    "btagSFlight_uncorrelated_2016"    : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
    "btagSFbc_uncorrelated_2016"       : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
    "btagSFlight_uncorrelated_2017"    : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
    "btagSFbc_uncorrelated_2017"       : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
    "btagSFlight_uncorrelated_2018"    : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
    "btagSFbc_uncorrelated_2018"       : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
}

# Takes two yield dictionaries, returns a dict of their ratio (for a given category)
# Drops mc stats uncertainty
def get_uncty_dict(cat,nom_dict,up_dict,do_dict):
    ratio_dict = {}
    for proc in nom_dict.keys():
        ratio_dict[proc] = {}
        n = nom_dict[proc][cat][0]
        u = up_dict[proc][cat][0]
        d = do_dict[proc][cat][0]
        ratio_dict[proc] = [u/n, d/n]
    return ratio_dict


# Sum the predicted yields over categorires to get asimov data number
def get_fake_data_for_ch(yld_dict,ch):
    data_obs = 0
    for proc in yld_dict.keys():
        data_obs += yld_dict[proc][ch][0]
    return data_obs


# Make the datacard for a given channel
def make_ch_card(ch,ch_ylds,ch_kappas=None,out_dir="."):

    # Building blocks we'll need to build the card formatting
    bin_str = f"bin_{ch}"
    syst_width = 0
    col_width = max(PRECISION*2+5,len(bin_str))
    line_break = "##----------------------------------\n"
    left_width = len(line_break) + 2
    left_width = max(syst_width+len("shape")+1,left_width)

    # The output name, location
    outf_card_name = f"test_card_{ch}.txt"
    print(f"Generating text file: {out_dir}/{outf_card_name}")
    outf_card_name = os.path.join(out_dir,outf_card_name)

    # Create the card for this channel
    with open(outf_card_name,"w") as f:

        # Shapes rows, not sure of the purpose of the shapes lines when we have no shape templates
        #f.write("shapes *        * {ch} FAKE\n")
        f.write("shapes *        * FAKE\n")
        f.write(line_break)
        f.write(f"bin         {bin_str}\n")
        f.write(f"observation {ch_ylds['data_obs']}\n")
        f.write(line_break)
        f.write(line_break)

        # Note: This list is what controls the columns in the text datacard, if a process appears
        proc_order = PROC_LST

        # Bin row
        row = [f"{'bin':<{left_width}}"] # Python string formatting is pretty great!
        for p in proc_order:
            row.append(f"{bin_str:>{col_width}}")
        row = " ".join(row) + "\n"
        f.write(row)

        # 1st process row
        row = [f"{'process':<{left_width}}"]
        for p in proc_order:
            row.append(f"{p:>{col_width}}")
        row = " ".join(row) + "\n"
        f.write(row)

        # 2nd process row
        row = [f"{'process':<{left_width}}"]
        bkgd_count =  1
        sgnl_count = -1
        for p in proc_order:
            if any([x in p for x in SIG_LST]): # Check for if the process is signal or not
                row.append(f"{sgnl_count:>{col_width}}")
                sgnl_count += -1
            else:
                row.append(f"{bkgd_count:>{col_width}}")
                bkgd_count += 1
        row = " ".join(row) + "\n"
        f.write(row)

        # Rate row
        row = [f"{'rate':<{left_width}}"]
        for p in proc_order:
            r = ch_ylds[p][0]
            if r < 0:
                raise Exception(f"\nERROR: Process {p} has negative total rate: {r:.3f}.\n")
            row.append(f"{r:>{col_width}.{PRECISION}f}")
        row = " ".join(row) + "\n"
        f.write(row)
        f.write(line_break)

        # Systematics rows
        for syst_name in ch_kappas:
            row = [f"{syst_name} lnN"]
            for p in proc_order:
                kappa_str = f"{ch_kappas[syst_name][p][0]}/{ch_kappas[syst_name][p][1]}"
                row.append(kappa_str)
            row = " ".join(row) + "\n"
            f.write(row)
        f.write(line_break)



def main():

    # Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file_name",help="Either json file of yields or pickle file with scikit hists")
    parser.add_argument("--out-dir","-d",default="./testcards",help="Output directory to write root and text datacard files to")
    parser.add_argument("--do-nuisance",action="store_true",help="Include nuisance parameters")
    parser.add_argument("--unblind",action="store_true",help="If set, use real data, otherwise use asimov data")

    args = parser.parse_args()
    in_file = args.in_file_name
    out_dir = args.out_dir
    do_nuis = args.do_nuisance
    unblind = args.unblind

    # Check args
    if do_nuis:
        raise Exception("Nuisance params not implimented yet.")
    if out_dir != "." and not os.path.exists(out_dir):
        print(f"Making dir \"{out_dir}\"")
        os.makedirs(out_dir)

    # Set up the dict of samples for full R2 and for each year individually
    sample_dict_mc = gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"all")
    sample_dict_mc_byyear = {
        "UL16APV" : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL16APV"),
        "UL16"    : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL16"),
        "UL17"    : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL17"),
        "UL18"    : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL18"),
    }

    # Get the yields dict from the input file
    #     - We can load a scikit hist (produced by wwz4l.py processor) and get the yields from that, dumpt to a dict
    #     - We can also load a json that contains the yields directly
    if in_file.endswith(".pkl.gz"):
        f = pickle.load(gzip.open(in_file))
        yld_dict = gy.get_yields(f,sample_dict_mc)
    elif in_file.endswith(".json"):
        with open(in_file) as f:
            yld_dict = json.load(f)
    else:
        raise Exception(f"ERROR: This script can only take hists or jsons, not files of type \"{in_file.split('.')[-1]}\".")

    ################
    # Systematics

    histo = f["njets"]

    # Get the list of systematic base names (i.e. without the up and down tags)
    # Assumes each syst has a "systnameUp" and a "systnameDown" category on the systematic axis
    syst_var_lst = []
    all_syst_var_lst = histo.axes["systematic"]
    for syst_var_name in all_syst_var_lst:
        if syst_var_name.endswith("Up"):
            syst_name_base = syst_var_name.replace("Up","")
            if syst_name_base not in syst_var_lst:
                syst_var_lst.append(syst_name_base)

    # Get just the systs that are correlated across years, for now
    sys_yr_correlated_lst = []
    for sys in syst_var_lst:
        if ("2016" in sys) or ("2016APV" in sys) or ("2017" in sys) or ("2018" in sys):
            pass
        else:
            sys_yr_correlated_lst.append(sys)

    # Find the syst variations for each
    kappa_dict = {}
    for cat in CAT_LST:
        kappa_dict[cat] = {}
        for sys in sys_yr_correlated_lst:
            yld_up = gy.get_yields(f,sample_dict_mc,systematic_name=f"{sys}Up")
            yld_do = gy.get_yields(f,sample_dict_mc,systematic_name=f"{sys}Down")
            kappa_dict[cat][sys] = get_uncty_dict(cat,yld_dict,yld_do,yld_up)

        for sys in SYSTS_SPECIAL:
            yr_rel = SYSTS_SPECIAL[sys]["yr_rel"]
            yld_yr_up = gy.get_yields(f,sample_dict_mc_byyear[yr_rel],systematic_name=f"{sys}Up")
            yld_yr_do = gy.get_yields(f,sample_dict_mc_byyear[yr_rel],systematic_name=f"{sys}Down")

            # Get nominal yld for all years other than the relevant one
            yld_yr_nom = gy.get_yields(f,sample_dict_mc_byyear[yr_rel])
            yld_all_nom = gy.get_yields(f,sample_dict_mc)
            yld_nom_all_but_yr  = utils.get_diff_between_nested_dicts(yld_all_nom,yld_yr_nom,difftype="absolute_diff") # This is: x = sum(nom yld for all years except relevant year)

            # Get the yield with just the relevant year's up/down variation varied (with nominal yld for all other years)
            yld_up = utils.get_diff_between_nested_dicts(yld_nom_all_but_yr, yld_yr_up, difftype="sum") # This is: x + (up   yld for relevant year)
            yld_do = utils.get_diff_between_nested_dicts(yld_nom_all_but_yr, yld_yr_do, difftype="sum") # This is: x + (down yld for relevant year)

            kappa_dict[cat][sys] = get_uncty_dict(cat,yld_dict,yld_do,yld_up)

    ################

    # Make the cards for each channel
    print(f"Making cards for {CAT_LST}. \nPutting in {out_dir}.")
    for ch in CAT_LST:

        # Get nominal yields for this category
        ch_ylds = {}
        for proc_name in yld_dict.keys():
            ch_ylds[proc_name] = yld_dict[proc_name][ch]

        # Overwrite data if necessary
        if not unblind:
            ch_ylds["data_obs"] = get_fake_data_for_ch(yld_dict,ch) # Make asimov data for now
        else:
            raise Exception("We are not unblinded yet.")

        # The kappas for this channel
        ch_kappas = kappa_dict[ch]

        # Make the card for this chan
        make_ch_card(ch,ch_ylds,ch_kappas,out_dir)


    print("Finished!")


if __name__ == "__main__":
    main()
