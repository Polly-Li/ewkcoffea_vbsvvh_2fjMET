import numpy as np
import os
import json
import pickle
import gzip
import argparse

from topcoffea.modules import utils
from ewkcoffea.modules.paths import ewkcoffea_path

import get_wwz_yields as gy # Note the fact that we're using functions from here means they probably belongs in ewkcoffea/ewkcoffea/modules

# TMP!!! To match the old script's order
TMP_SYS_ORDER = [
    "btagSFlight_correlated",
    "btagSFbc_correlated",
    "lepSF_elec",
    "lepSF_muon",
    "PreFiring",
    "PU",
    "renorm",
    "fact",
    "ISR",
    "FSR",
    "btagSFlight_uncorrelated_2016APV",
    "btagSFbc_uncorrelated_2016APV",
    "btagSFlight_uncorrelated_2016",
    "btagSFbc_uncorrelated_2016",
    "btagSFlight_uncorrelated_2017",
    "btagSFbc_uncorrelated_2017",
    "btagSFlight_uncorrelated_2018",
    "btagSFbc_uncorrelated_2018",
]


# Global variables
PRECISION = 6   # Decimal point precision in the text datacard output
PROC_LST = ["WWZ","ZH","ZZ","ttZ","tWZ","other"]
SIG_LST = ["WWZ","ZH"]
CAT_LST_CB = ["sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]

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

BKG_TF_MAP = {

    "ZZ" : {
        "sr_4l_sf_A" : "cr_4l_sf",
        "sr_4l_sf_B" : "cr_4l_sf",
        "sr_4l_sf_C" : "cr_4l_sf",
        "sr_4l_of_1" : "cr_4l_sf",
        "sr_4l_of_2" : "cr_4l_sf",
        "sr_4l_of_3" : "cr_4l_sf",
        "sr_4l_of_4" : "cr_4l_sf",

    },
    "ttZ" : {
        "sr_4l_sf_A" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_sf_B" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_sf_C" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_of_1" : "cr_4l_btag_of",
        "sr_4l_of_2" : "cr_4l_btag_of",
        "sr_4l_of_3" : "cr_4l_btag_of",
        "sr_4l_of_4" : "cr_4l_btag_of",
    }
}


# Make the datacard for a given channel
def make_ch_card(ch,proc_order,ch_ylds,ch_kappas=None,out_dir="."):

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
        f.write("shapes *        * FAKE\n")
        f.write(line_break)
        f.write(f"bin         {bin_str}\n")
        f.write(f"observation {ch_ylds['data_obs']}\n")
        f.write(line_break)
        f.write(line_break)

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
            r = ch_ylds[p]
            row.append(r)
        row = " ".join(row) + "\n"
        f.write(row)
        f.write(line_break)

        # Systematics rows
        if ch_kappas is not None:
            #for syst_name in ch_kappas:
            ### TMP so we can match old card order for the diffs ###
            if set(TMP_SYS_ORDER) != set(ch_kappas.keys()):
                raise Exception("THIS IS HERE")
            for syst_name in TMP_SYS_ORDER:
            ###
                row = [f"{syst_name} lnN"]
                for p in proc_order:
                    kappa_str = ch_kappas[syst_name][p]
                    row.append(kappa_str)
                row = " ".join(row) + "\n"
                f.write(row)
            f.write(line_break)


# Get the yields (nested in the order: year, category, syst, proc)
def get_yields(histo,sample_dict,blind=True,systematic_name=None):

    yld_dict = {}

    if systematic_name is None: syst_lst = histo.axes["systematic"]
    else: syst_lst = [systematic_name]

    # Look at the yields in the histo
    for cat_name in histo.axes["category"]:
        yld_dict[cat_name] = {}
        for syst_name in syst_lst:
            #if syst_name not in ["nominal", "PUUp", "PUDown", "btagSFbc_uncorrelated_2016APVUp", "btagSFbc_uncorrelated_2016APVDown"]: continue # TMP !!!
            yld_dict[cat_name][syst_name ] = {}
            for proc_name in sample_dict.keys():
                if blind and (("data" in proc_name) and (not cat_name.startswith("cr_"))):
                    # If this is data and we're not in a CR category, put placeholder numbers for now
                    yld_dict[cat_name][proc_name] = [-999,-999]
                else:
                    val = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":syst_name }].values(flow=True)))
                    var = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":syst_name }].variances(flow=True)))
                    yld_dict[cat_name][syst_name ][proc_name] = [val,var]

    return yld_dict


# Modifies the "FR2" subdict of the yield dict to have proper up and down variations for the per year systematics
# Because of how we fill in the processor, the yields for per year systs come _only_ from that year
# So this function adds the nominal yields from the other three years to the up/down variation for the relevant year
# Note the in_dict is modifed in place (we do not return a copy of the dict)
def handle_per_year_systs_for_fr2(in_dict, systs_special=SYSTS_SPECIAL):
    for cat in in_dict["FR2"].keys():
        for sys in systs_special:
            # Find up/down variation for the year relevant to that syst
            yrrel = systs_special[sys]["yr_rel"] # The relevant year for this special syst
            yld_yrrel_up = in_dict[yrrel][cat][f"{sys}Up"]
            yld_yrrel_do = in_dict[yrrel][cat][f"{sys}Down"]
            # Get nominal yld for all years other than the relevant one
            yld_yrrel_nom = in_dict[yrrel][cat]["nominal"]
            yld_yrall_nom = in_dict["FR2"][cat]["nominal"]
            yld_allbutyrrel_nom = utils.get_diff_between_dicts(yld_yrall_nom,yld_yrrel_nom,difftype="absolute_diff") # This is: x = sum(nom yld for all years except relevant year)
            # Get the yield with just the relevant year's up/down variation varied (with nominal yld for all other years)
            yld_up = utils.get_diff_between_dicts(yld_allbutyrrel_nom, yld_yrrel_up, difftype="sum") # This is: x + (up   yld for relevant year)
            yld_do = utils.get_diff_between_dicts(yld_allbutyrrel_nom, yld_yrrel_do, difftype="sum") # This is: x + (down yld for relevant year)
            in_dict["FR2"][cat][f"{sys}Up"] = yld_up
            in_dict["FR2"][cat][f"{sys}Down"] = yld_do



#####################################################
# Get kappa dict
def get_kappa_dict(in_dict_mc,in_dict_data,bkg_tf_map):

    # Get the list of systematic base names (i.e. without the up and down tags)
    #     - Assumes each syst has a "systnameUp" and a "systnameDown"
    #     - This will drop nominal (since there is no "nominalUp" to tag on)
    def get_syst_base_name_lst(in_lst):
        out_lst = []
        for syst in in_lst:
            if syst.endswith("Up"):
                syst_name_base = syst.replace("Up","")
                if syst_name_base not in out_lst:
                    out_lst.append(syst_name_base)
        return out_lst


    # Takes two pairs [val1,var1] and [val2,var2], returns the product with error propagated
    # Note the var is like sumw2 i.e. alreayd squared (in both input and output)
    def valvar_op(valvar_1, valvar_2, op):
        #valvar_out = [None,None]

        val1 = valvar_1[0]
        var1 = valvar_1[1]
        val2 = valvar_2[0]
        var2 = valvar_2[1]

        if op == "prod":
            val = val1*val2
        elif op == "div":
            val = val1/val2
        else:
            raise Exception("Unknown operatation")

        var = (val**2) * ( (np.sqrt(var1)/val1)**2 + (np.sqrt(var2)/val2)**2 )

        return [val,var]


    kappa_dict = {}
    for cat in in_dict_mc.keys():
        kappa_dict[cat] = {}
        for sys in get_syst_base_name_lst(list(in_dict_mc[cat].keys())):
            kappa_dict[cat][sys] = {}
            for proc in in_dict_mc[cat]["nominal"]:
                kappa_dict[cat][sys][proc] = {}
                valvar_up = in_dict_mc[cat][f"{sys}Up"][proc]
                valvar_do = in_dict_mc[cat][f"{sys}Down"][proc]
                valvar_nom = in_dict_mc[cat]["nominal"][proc]
                valvar_kappa_up = valvar_op(valvar_up,valvar_nom,"div")
                valvar_kappa_do = valvar_op(valvar_do,valvar_nom,"div")

                kappa_dict[cat][sys][proc]["Up"] = valvar_kappa_up
                kappa_dict[cat][sys][proc]["Down"] = valvar_kappa_do

    return kappa_dict


## VERY IN PROGRESS
#def do_tf():
#    # Do the TF calculation
#    if proc in bkg_tf_map:
#        if cat not in bkg_tf_map[proc]:
#            # Skip the TF calculation for categories that we have not defined it
#            #print(f"Warning: the \"{cat}\" category does not have a TF CR defined. Skipping.")
#            continue
#        cr_name = bkg_tf_map[proc][cat]
#        valvar_cr_mc   = in_dict_mc[cr_name]["nominal"][proc]
#        valvar_cr_data = in_dict_data[cr_name]["nominal"]["data"]
#
#        print("PROC!!!",proc)
#        print("mc",valvar_cr_mc)
#        print("da",valvar_cr_data)
#####################################################


# Get just the numbers we want for rate row for datacard
# Also sum all MC rates together into asimov number
# Assumes in_dict has nested keys: cat,syst,proc
def get_rate_for_dc(in_dict):
    out_dict = {}
    for cat in in_dict:
        out_dict[cat] = {}
        asimov_data = 0
        for proc in in_dict[cat]["nominal"]:
            rate = in_dict[cat]["nominal"][proc][0]
            if rate < 0:
                print(f"\nWarning: Process \"{proc}\" in \"{cat}\" has negative total rate: {rate}.\n")
            out_dict[cat][proc] = str(rate)
            asimov_data += rate
        out_dict[cat]["data_obs"] = str(asimov_data)
    return out_dict


# Assumes in_dict has nested keys: cat,syst_base,proc,"Up"(or "Down)
# Takes just the val (dropps uncty)
def get_kappa_for_dc(in_dict):
    out_dict = {}
    for cat in in_dict:
        out_dict[cat] = {}
        for systname_base in in_dict[cat]:
            out_dict[cat][systname_base] = {}
            for proc in in_dict[cat][systname_base]:
                u = in_dict[cat][systname_base][proc]['Up'][0]
                d = in_dict[cat][systname_base][proc]['Down'][0]
                out_dict[cat][systname_base][proc] = f"{d}/{u}"
    return(out_dict)


#####################################################


def main():

    # Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file_name",help="Either json file of yields or pickle file with scikit hists")
    parser.add_argument("--out-dir","-d",default="./testcards",help="Output directory to write root and text datacard files to")
    parser.add_argument("-s","--do-nuisance",action="store_true",help="Include nuisance parameters")
    parser.add_argument("--unblind",action="store_true",help="If set, use real data, otherwise use asimov data")

    args = parser.parse_args()
    in_file = args.in_file_name
    out_dir = args.out_dir
    do_nuis = args.do_nuisance
    unblind = args.unblind

    # Check args
    if out_dir != "." and not os.path.exists(out_dir):
        print(f"Making dir \"{out_dir}\"")
        os.makedirs(out_dir)

    # Get the histo
    f = pickle.load(gzip.open(in_file))
    histo = f["njets"] # Let's use njets

    # Get the dictionary defining the mc sample grouping
    sample_names_dict_data = {"FR2" : gy.create_data_sample_dict("all")}
    sample_names_dict_mc = {
        "UL16APV" : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL16APV"),
        "UL16"    : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL16"),
        "UL17"    : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL17"),
        "UL18"    : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL18"),
        "FR2"     : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"all"),
    }

    # Get yield dictionary (nested in the order: year, category, syst, proc)
    yld_dict_mc = {}
    for year in sample_names_dict_mc:
        yld_dict_mc[year] = get_yields(histo,sample_names_dict_mc[year])
    handle_per_year_systs_for_fr2(yld_dict_mc)

    # We're only looking at Full R2 for now
    yld_dict_mc = yld_dict_mc["FR2"]
    yld_dict_data = get_yields(histo,sample_names_dict_data["FR2"])

    # Get the ratios to nominal
    kappa_dict = get_kappa_dict(yld_dict_mc,yld_dict_data,BKG_TF_MAP)

    # Get just the info we want to put in the card in str form
    yld_rate_for_dc = get_rate_for_dc(yld_dict_mc)
    kappa_for_dc = get_kappa_for_dc(kappa_dict)

    #print(kappa_for_dc)
    #exit()


    # Make the cards for each channel
    print(f"Making cards for {CAT_LST_CB}. \nPutting in {out_dir}.")
    for ch in CAT_LST_CB:

        # Make the card for this chan
        make_ch_card(
            ch,
            PROC_LST,
            yld_rate_for_dc[ch],
            kappa_for_dc[ch],
            out_dir,
        )


    print("Finished!")


if __name__ == "__main__":
    main()
