import os
import json
import pickle
import gzip
import argparse

from topcoffea.modules import utils
from ewkcoffea.modules.paths import ewkcoffea_path

import get_wwz_yields as gy # Note the fact that we're using functions from here means they probably belongs in ewkcoffea/ewkcoffea/modules


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


# Get the yields
def get_yields(histo,sample_dict,blind=True,systematic_name=None):

    yld_dict = {}

    if systematic_name is None: syst_lst = histo.axes["systematic"]
    else: syst_lst = [systematic_name]

    # Look at the yields in the histo
    for cat_name in histo.axes["category"]:
        yld_dict[cat_name] = {}
        for syst_name in syst_lst:
            yld_dict[cat_name][syst_name ] = {}
            for proc_name in sample_dict.keys():
                if blind and (("data" in proc_name) and (not cat_name.startswith("cr_"))):
                    # If this is data and we're not in a CR category, put placeholder numbers for now
                    yld_dict[proc_name][cat_name] = [-999,-999]
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
    sample_names_dict_data = gy.create_data_sample_dict("all")
    sample_names_dict_mc = {
        "UL16APV" : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL16APV"),
        "UL16"    : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL16"),
        "UL17"    : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL17"),
        "UL18"    : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"UL18"),
        "FR2"     : gy.create_mc_sample_dict(gy.SAMPLE_DICT_BASE,"all"),
    }

    # Get yield dictionary (nested in the order: year, category, syst, proc)
    yld_dict = {}
    for year in sample_names_dict_mc:
        yld_dict[year] = get_yields(histo,sample_names_dict_mc[year])
    handle_per_year_systs_for_fr2(yld_dict)



    #exit()
    ### Print stuff
    print("")
    print(yld_dict)
    for yr in yld_dict:
        print(yr)
        if yr != "FR2": continue
        for cat in yld_dict[yr]:
            print("\t",cat)
            for sys in yld_dict[yr][cat]:
                #print("\t\t",sys,yld_dict[yr][cat][sys].keys())
                #if cat == "sr_4l_of_4" and sys == "btagSFbc_uncorrelated_2018Up":
                print("\t\t",sys,yld_dict[yr][cat][sys])

main()
