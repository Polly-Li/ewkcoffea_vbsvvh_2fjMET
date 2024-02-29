import numpy as np
import os
import copy
import json
import pickle
import gzip
import argparse

from topcoffea.modules import utils
from ewkcoffea.modules.paths import ewkcoffea_path

import sample_groupings as sg # Note the fact that we're using functions from here means they probably belongs in ewkcoffea/ewkcoffea/modules

SMALL = 0.000001

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
    "lumi",
    "qcd_scale_ttZ",
]


# Global variables
PRECISION = 6   # Decimal point precision in the text datacard output
PROC_LST = ["WWZ","ZH","ZZ","ttZ","tWZ","WZ","other"]
SIG_LST = ["WWZ","ZH"]
CAT_LST_CB = ["sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]

# Systs that are not correlated across years
SYSTS_SPECIAL = {
    "btagSFlight_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
    "btagSFbc_uncorrelated_2016APV"    : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]}, # !!!
    "btagSFlight_uncorrelated_2016"    : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
    "btagSFbc_uncorrelated_2016"       : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
    "btagSFlight_uncorrelated_2017"    : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
    "btagSFbc_uncorrelated_2017"       : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
    "btagSFlight_uncorrelated_2018"    : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
    "btagSFbc_uncorrelated_2018"       : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
}


########### Writing the datacard ###########

# Make the datacard for a given channel
def make_ch_card(ch,proc_order,ch_ylds,ch_kappas=None,ch_gmn=None,out_dir="."):

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
            ### TMP so we can match old card order for the diffs !!! ###
            #if set(TMP_SYS_ORDER) != set(ch_kappas.keys()):
            #    raise Exception("THIS IS BAD HERE")
            #    pass
            #for syst_name in TMP_SYS_ORDER:
            ###
            for syst_name in ch_kappas:
                row = [f"{syst_name} lnN"]
                for p in proc_order:
                    kappa_str = ch_kappas[syst_name][p]
                    row.append(kappa_str)
                row = " ".join(row) + "\n"
                f.write(row)

        # Systematics rows for gmN
        if ch_gmn is not None:
            for name in ch_gmn:
                row = [name]
                for p in proc_order:
                    alpha_str = ch_gmn[name][p]
                    row.append(alpha_str)
                row = " ".join(row) + "\n"
                f.write(row)
            f.write(line_break)
        else:
            f.write(line_break)


########### General ###########

# Takes two pairs [val1,var1] and [val2,var2], returns the product or sum with error propagated
# Note the var is like sumw2 i.e. alreayd squared (in both input and output)
def valvar_op(valvar_1, valvar_2, op):
    val1 = valvar_1[0]
    var1 = valvar_1[1]
    val2 = valvar_2[0]
    var2 = valvar_2[1]

    if op == "prod":
        val = val1*val2
        var = (val**2) * ( (np.sqrt(var1)/val1)**2 + (np.sqrt(var2)/val2)**2 )
    elif op == "div":
        val = val1/val2
        var = (val**2) * ( (np.sqrt(var1)/val1)**2 + (np.sqrt(var2)/val2)**2 )
    elif op == "sum":
        val = val1 + val2
        var = var1 + var2
    elif op == "diff":
        val = val1 - val2
        var = var1 + var2
    else:
        raise Exception("Unknown operatation")

    return [val,var]


########### Getting and manipulating yields ###########

# Get the yields (nested in the order: year,cat,syst,proc)
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
                    yld_dict[cat_name][syst_name][proc_name] = [-999,-999]
                else:
                    val = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":syst_name }].values(flow=True)))
                    var = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":syst_name }].variances(flow=True)))
                    yld_dict[cat_name][syst_name][proc_name] = [val,var]

    return yld_dict


# Modify the yields dict to properly calculate the per-year systs
#   - Needs all years in input dict, expects input dict with keys: year,cat,syst,proc
#   - Modifies the "FR2" subdict of the yield dict to have proper up and down variations for the per year systematics
#   - Because of how we fill in the processor, the yields for per year systs come _only_ from that year
#   - So this function adds the nominal yields from the other three years to the up/down variation for the relevant year
#   - Note the in_dict is modifed in place (we do not return a copy of the dict)
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


# Get rid of negative values in the yld dict
#   - Replace the value with SMALL
#   - And add |value| to the stat error to be more conservative
#   - Shift the up/down variations to be centered around SMALL (does not touch stat uncertainty on up/down)
def handle_negatives(in_dict):
    out_dict = copy.deepcopy(in_dict)
    for cat in in_dict:
        for proc in in_dict[cat]["nominal"]:
            val = in_dict[cat]["nominal"][proc][0]
            var = in_dict[cat]["nominal"][proc][1]
            if val < 0:
                print(f"WARNING: Process \"{proc}\" in cat \"{cat}\" is negative ({val}), replacing with {SMALL} and shifting up/down systematic variations accordingly.")
                out_dict[cat]["nominal"][proc][0] = SMALL
                out_dict[cat]["nominal"][proc][1] = (abs(val) + np.sqrt(var))**2
                for syst in out_dict[cat]:
                    if syst == "nominal": continue # Already handled this one
                    syst_var_orig = out_dict[cat][syst][proc]
                    out_dict[cat][syst][proc][0] = (syst_var_orig - val) + SMALL # Center around SMALL

    return out_dict



########### Regarding uncertainties ###########

# Get the rate systs from the intput json, dump into dict (with nested keys: syst, proc)
# Outputs strings, ready to be dumped into datacard
def get_rate_systs(proc_lst):
    syst_json = ewkcoffea_path("params/rate_systs.json")
    with open(syst_json) as f_systs: rate_systs_dict = json.load(f_systs)

    # Build up the dictionary
    out_dict = {}

    # Make row for rate uncertainties that impact all processes (for now just lumi)
    for uncty_name in rate_systs_dict["rate_uncertainties_all_proc"]:
        out_dict[uncty_name] = {}
        for proc in proc_lst:
            out_dict[uncty_name][proc] = str(rate_systs_dict["rate_uncertainties_all_proc"][uncty_name])

    # Make rows for rate uncertainties that impact a subset of the processes
    for uncty_name in rate_systs_dict["rate_uncertainties_some_proc"]:
        for proc_of_interest in rate_systs_dict["rate_uncertainties_some_proc"][uncty_name]["procs"]:
            out_dict[f"{uncty_name}_{proc_of_interest}"] = {}
            for proc_itr in proc_lst:
                if proc_itr == proc_of_interest:
                    uncty = str(rate_systs_dict["rate_uncertainties_some_proc"][uncty_name]["val"])
                else:
                    uncty = "-"
                out_dict[f"{uncty_name}_{proc_of_interest}"][proc_itr] = uncty

    return out_dict


# Get kappa dict (e.g. up/nom ratios) from the dict of all histograms
def get_kappa_dict(in_dict_mc,in_dict_data):

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

    # Get the kappa dict
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

                # Handle negative cases
                if (valvar_kappa_up[0]<0) and (valvar_kappa_do[0]<0): raise Exception("Both kappas negative, should not be possible.")
                if valvar_kappa_up[0] < 0:
                    print(f"WARNING: Up var for {sys} for {proc} for {cat} is negative, setting to {SMALL}.")
                    valvar_kappa_up[0] = SMALL
                if valvar_kappa_do[0] < 0:
                    valvar_kappa_do[0] = SMALL
                    print(f"WARNING: Down var for {sys} for {proc} for {cat} is negative, setting to {SMALL}.")

                kappa_dict[cat][sys][proc]["Up"] = valvar_kappa_up
                kappa_dict[cat][sys][proc]["Down"] = valvar_kappa_do


    return kappa_dict


# Calculate the background estimation from relevant CRs
def do_tf(yld_mc,yld_data,kappas,tf_map,quiet=True):

    yld_mc_out = copy.deepcopy(yld_mc)
    kappas_out = copy.deepcopy(kappas)
    gmn_alpha_out = {}

    # Loop over cat and do NSF calculation for each relevant proc in each cat
    for cat in yld_mc:
        # Just nonimal for now
        for proc_of_interest in yld_mc[cat]["nominal"]:

            # Skip procs that do not get TF calculations
            if proc_of_interest not in tf_map: continue
            elif cat not in tf_map[proc_of_interest]:
                #print(f"Warning, cat \"{cat}\" not defined for this proc.")
                continue

            # Otherwise we go ahead and do the background estimation stuff
            else:

                # Get the nominal mc and data yields in the CR
                cr_name = tf_map[proc_of_interest][cat]
                valvar_cr_mc   = yld_mc[cr_name]["nominal"][proc_of_interest]
                valvar_cr_data = yld_data[cr_name]["nominal"]["data"]

                # Sum up all contributions in the CR besides bkg of interest
                valvar_bkg_all_but_bkg_of_interest = [0,0]
                for p in yld_mc[cr_name]["nominal"]:
                    if p != proc_of_interest:
                        valvar_bkg_all_but_bkg_of_interest[0] += yld_mc[cr_name]["nominal"][p][0]
                        valvar_bkg_all_but_bkg_of_interest[1] += yld_mc[cr_name]["nominal"][p][1]

                # Subtract those extra background contributions from the data
                valvar_cr_data_corrected = valvar_op(valvar_cr_data, valvar_bkg_all_but_bkg_of_interest, "diff")

                # Calculate the NSF = N_CR_with_other_bkg_subtracted / MC_CR, use this to scale the yld
                valvar_nsf = valvar_op(valvar_cr_data_corrected, valvar_cr_mc, "div")
                valvar_bkg = valvar_op(valvar_nsf, yld_mc[cat]["nominal"][proc_of_interest], "prod")

                # Dump NSF value
                if not quiet:
                    relerr = 100*(np.sqrt(valvar_nsf[1])/valvar_nsf[0])
                    print(f"NSF for {cat} {proc_of_interest}: {np.round(valvar_nsf[0],2)} +- {np.round(relerr,2)}%")
                    #print(f"NSF for {cat} {proc_of_interest}: {np.round(valvar_nsf[0],2)} +- {np.round(np.sqrt(valvar_nsf[1]),2)}")

                # Put the old yield times the NSF into the out dict
                yld_mc_out[cat]["nominal"][proc_of_interest] = valvar_bkg

                # Get the gmN numbers
                valvar_alpha = valvar_op(valvar_bkg,valvar_cr_data, "div")
                if cat not in gmn_alpha_out: gmn_alpha_out[cat] = {}
                if cr_name not in gmn_alpha_out[cat]: gmn_alpha_out[cat][cr_name] = {"N": valvar_cr_data[0], "proc_alpha": {}}
                gmn_alpha_out[cat][cr_name]["proc_alpha"][proc_of_interest] = valvar_alpha[0]

                ### Handle the kappas ###

                if kappas is not None:

                    # Loop over syst and replace the kappas with the e.g. MC_SR_up/MC_CR_up
                    for syst_base_name in kappas[cat]:

                        sr_up = kappas[cat][syst_base_name][proc_of_interest]["Up"]
                        sr_do = kappas[cat][syst_base_name][proc_of_interest]["Down"]
                        cr_up = kappas[cr_name][syst_base_name][proc_of_interest]["Up"]
                        cr_do = kappas[cr_name][syst_base_name][proc_of_interest]["Down"]
                        new_kappa_up = valvar_op(sr_up,cr_up,"div")
                        new_kappa_do = valvar_op(sr_do,cr_do,"div")

                        kappas_out[cat][syst_base_name][proc_of_interest]["Up"] = new_kappa_up
                        kappas_out[cat][syst_base_name][proc_of_interest]["Down"] = new_kappa_do

    return [yld_mc_out, kappas_out, gmn_alpha_out]




# Get the MC stats from the yld dict and put it into kappa dict
# Can optionally skip processes (e.g. skip data driven ZZ and ttZ since do these different via gmN)
def add_stats_kappas(yld_mc, kappas, skip_procs=[]):

    kappas_out = copy.deepcopy(kappas)

    # Loop over cat and do NSF calculation for each relevant proc in each cat
    for cat in yld_mc:

        # Make a row for proc_of_interest's stats
        for proc_of_interest in yld_mc[cat]["nominal"]:
            if proc_of_interest in skip_procs: continue
            kappas_out[cat][f"stats_{cat}_{proc_of_interest}"] = {}
            # Now fill the columns for proc_of_interest's row
            for proc_itr in yld_mc[cat]["nominal"]:
                # Most columns are not relevant (will just be a "-" in the datacard)
                kappas_out[cat][f"stats_{cat}_{proc_of_interest}"][proc_itr] = {"Up": [None,None], "Down": [None,None]}
                # But for the column that goes with this proc_of_interest, fill the actual numbers
                if proc_itr == proc_of_interest:
                    valvar = yld_mc[cat]["nominal"][proc_of_interest]
                    up = (valvar[0] + np.sqrt(valvar[1]))/valvar[0]
                    do = (valvar[0] - np.sqrt(valvar[1]))/valvar[0]
                    # Clip at 0
                    if do < 0:
                        print(f"WARNING: For cat \"{cat}\" and proc \"{proc_of_interest}\", the uncertainty {np.sqrt(valvar[1])} is larger than the value {valvar[0]}. Clipping down variation to {SMALL}.")
                        do = SMALL
                    kappas_out[cat][f"stats_{cat}_{proc_of_interest}"][proc_itr]["Up"]   = [up, None] # Rel err up, do not include error on the error (just leave as None)
                    kappas_out[cat][f"stats_{cat}_{proc_of_interest}"][proc_itr]["Down"] = [do, None] # Rel err down, do not include error on the error (just leave as None)

    return kappas_out



########### Put stuff into form to pass to the function to write out cards ###########

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
            #out_dict[cat][proc] = "{:.6f}".format(np.round(rate,6)) ### TMP!!!
            asimov_data += rate
        out_dict[cat]["data_obs"] = str(asimov_data)
        #out_dict[cat]["data_obs"] = str(np.round(asimov_data,6)) ### TMP!!!
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
                if (u is None) and (d is None):
                    out_dict[cat][systname_base][proc] = "-"
                else:
                    out_dict[cat][systname_base][proc] = f"{d}/{u}"
    return out_dict


# Prepares the gmN line for the datacard
# The proc_lst shoudl be all processes
# Assumes this is for a single category, so dict looks like {cr_name: {"N": 60, "proc_alpha": {"proc": value}}}
def get_gmn_for_dc(in_dict,proc_lst):
    out_dict = {}
    for cr_name in in_dict:
        N = in_dict[cr_name]['N']
        if int(N)!= N: raise Exception(f"ERROR: Why is the number of events in your CR ({N}) not an int?")
        row_name = f"stats_{cr_name} gmN {int(N)}"
        out_dict[row_name] = {}
        for p_itr in proc_lst:
            if p_itr in in_dict[cr_name]["proc_alpha"]:
                out = str(in_dict[cr_name]["proc_alpha"][p_itr])
            else:
                out = "-"
            out_dict[row_name][p_itr] = out
    return out_dict


#####################################
########### Main function ###########

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
    sample_names_dict_data = {"FR2" : sg.create_data_sample_dict("all")}
    sample_names_dict_mc = {
        "UL16APV" : sg.create_mc_sample_dict(sg.SAMPLE_DICT_BASE,"UL16APV"),
        "UL16"    : sg.create_mc_sample_dict(sg.SAMPLE_DICT_BASE,"UL16"),
        "UL17"    : sg.create_mc_sample_dict(sg.SAMPLE_DICT_BASE,"UL17"),
        "UL18"    : sg.create_mc_sample_dict(sg.SAMPLE_DICT_BASE,"UL18"),
        "FR2"     : sg.create_mc_sample_dict(sg.SAMPLE_DICT_BASE,"all"),
    }

    # Get yield dictionary (nested in the order: year,cat,syst,proc)
    yld_dict_mc_allyears = {}
    for year in sample_names_dict_mc:
        yld_dict_mc_allyears[year] = get_yields(histo,sample_names_dict_mc[year])
    handle_per_year_systs_for_fr2(yld_dict_mc_allyears)

    # We're only looking at Full R2 for now
    yld_dict_mc = yld_dict_mc_allyears["FR2"]
    yld_dict_data = get_yields(histo,sample_names_dict_data["FR2"])

    # Print info about a bin
    printinfo = 1
    if printinfo:
        s = "renorm"
        p = "WWZ"
        c = "sr_4l_of_2"
        cr = "cr_4l_btag_of"
        print(p,c,s,cr)
        print("mc sr n",yld_dict_mc[c]["nominal"][p])
        print("mc sr u",yld_dict_mc[c][f"{s}Up"][p])
        print("mc cr n",yld_dict_mc[cr]["nominal"][p])
        print("mc cr u",yld_dict_mc[cr][f"{s}Up"][p])
        print("da cr n",yld_dict_data[cr]["nominal"]["data"])
        print("da cr u",yld_dict_data[cr][f"{s}Up"]["data"])
        #exit()


    # Get the syst ratios to nominal (i.e. kappas)
    kappa_dict = get_kappa_dict(yld_dict_mc,yld_dict_data)

    # Do the TF calculation
    yld_dict_mc, kappa_dict, gmn_dict = do_tf(yld_dict_mc,yld_dict_data,kappa_dict,sg.BKG_TF_MAP)

    # Get rid of negative yields (and recenter syst variations around SMALL), should happen before computing kappas
    yld_dict_mc = handle_negatives(yld_dict_mc)

    # Add mc stats to kappa dict (important to be after TF calculation so that data driven bkg stats include stats of CRs)
    kappa_dict = add_stats_kappas(yld_dict_mc,kappa_dict,skip_procs=["ZZ","ttZ"])




    # Get just the info we want to put in the card in str form
    # Maybe move this to inside the card making loop
    rate_for_dc = get_rate_for_dc(yld_dict_mc)
    kappa_for_dc = get_kappa_for_dc(kappa_dict)


    # Make the cards for each channel
    print(f"Making cards for {CAT_LST_CB}. \nPutting in {out_dir}.")
    for ch in CAT_LST_CB:

        # The rates for this channel
        rate_for_dc_ch = rate_for_dc[ch]

        # Get the kappa and gamma dict for this channel if we are doing systs
        kappa_for_dc_ch = None
        gmn_for_dc_ch   = None
        if do_nuis:
            kappa_for_dc_ch = kappa_for_dc[ch]
            kappa_for_dc_ch.update(get_rate_systs(PROC_LST)) # Append in the ones from rate json
            gmn_for_dc_ch = get_gmn_for_dc(gmn_dict[ch],proc_lst=list(sg.SAMPLE_DICT_BASE.keys()))


        # Make the card for this chan
        make_ch_card(
            ch,
            PROC_LST,
            rate_for_dc_ch,
            kappa_for_dc_ch,
            gmn_for_dc_ch,
            out_dir,
        )


    print("Finished!")


if __name__ == "__main__":
    main()
