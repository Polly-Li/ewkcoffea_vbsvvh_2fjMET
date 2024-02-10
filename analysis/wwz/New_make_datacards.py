import numpy as np
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
    #"btagSFlight_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
    "btagSFbc_uncorrelated_2016APV"    : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
    #"btagSFlight_uncorrelated_2016"    : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
    #"btagSFbc_uncorrelated_2016"       : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
    #"btagSFlight_uncorrelated_2017"    : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
    #"btagSFbc_uncorrelated_2017"       : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
    #"btagSFlight_uncorrelated_2018"    : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
    #"btagSFbc_uncorrelated_2018"       : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
}


# Get the yields (nested in the order: year, category, syst, proc)
def get_yields(histo,sample_dict,blind=True,systematic_name=None):

    yld_dict = {}

    if systematic_name is None: syst_lst = histo.axes["systematic"]
    else: syst_lst = [systematic_name]

    # Look at the yields in the histo
    for cat_name in histo.axes["category"]:
        yld_dict[cat_name] = {}
        for syst_name in syst_lst:
            #if syst_name not in ["nominal", "btagSFlight_correlatedUp", "btagSFlight_correlatedDown", "PUUp", "PUDown", "btagSFbc_uncorrelated_2016APVUp", "btagSFbc_uncorrelated_2016APVDown"]: continue # TMP !!!
            if syst_name not in ["nominal", "PUUp", "PUDown", "btagSFbc_uncorrelated_2016APVUp", "btagSFbc_uncorrelated_2016APVDown"]: continue # TMP !!!
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


# Takes a full yield dict and returns the yields for just the procs we care about for datacard
# Specify the year of interest, and category list, (and nominal syst is returned by default)
def get_rate_dict(in_dict,yr_key,cat_lst):
    out_dict = {}
    for cat in cat_lst:
        out_dict[cat] = in_dict[yr_key][cat]["nominal"]
    return out_dict


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

    # We're only looking at Full R2 for now
    yld_dict = yld_dict["FR2"]


    # Get just the yields we care about for now

    def do_bkg_tf(in_dict):
        print(in_dict.keys())

    def get_kappa_dict(in_dict):

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


        for cat in in_dict.keys():
            print("\n",cat)
            for sys in get_syst_base_name_lst(list(in_dict[cat].keys())):
                print("\t",sys)
                for proc in in_dict[cat]["nominal"]:
                    valvar_up = in_dict[cat][f"{sys}Up"][proc]
                    valvar_do = in_dict[cat][f"{sys}Down"][proc]
                    valvar_nom = in_dict[cat]["nominal"][proc]
                    valvar_kappa_up = valvar_op(valvar_up,valvar_nom,"div")
                    valvar_kappa_do = valvar_op(valvar_do,valvar_nom,"div")
                    print("\t\t",proc,valvar_nom)
                    print("\t\t",proc,valvar_up)
                    print("\t\t",proc,valvar_do)
                    print("\t\t",proc,valvar_kappa_up)
                    print("\t\t",proc,valvar_kappa_do)


    tf_info_dict = {

        'sr_4l_sf_A' : {
            "ZZ" : {
                "sr" : "sr_sf_all",
                "cr" : "cr_4l_sf",
            },
            "ttZ" : {
                "sr" : "sr_sf_all",
                "cr" : "cr_4l_btag_sf_offZ_met80",
            },
        },
    }

    #for cat in tf_info_dict:
        #zz_yld = n_cr * mc_sr / mc_cr

    #'sr_4l_sf_B' :
    #'sr_4l_sf_C' :

    #'sr_4l_of_1' :
    #'sr_4l_of_2' :
    #'sr_4l_of_3' :
    #'sr_4l_of_4' :

    #'cr_4l_btag_of'
    #'cr_4l_btag_sf_offZ_met80'
    #'cr_4l_sf'


    #x = do_bkg_tf(yld_dict)
    x = get_kappa_dict(yld_dict)




    #rate_dict = get_rate_dict(yld_dict,yr_key="FR2",cat_lst=CAT_LST_CB)
    #kappa_dict = get_kappa_dict(yld_didct
    #print("Before",yld_dict)
    #print("After",rate_dict)


    exit()
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
