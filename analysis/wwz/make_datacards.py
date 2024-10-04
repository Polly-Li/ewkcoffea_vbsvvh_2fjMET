import numpy as np
import os
import copy
import json
import pickle
import gzip
import argparse

from topcoffea.modules import utils
from ewkcoffea.modules.paths import ewkcoffea_path
import ewkcoffea.modules.sample_groupings as sg
import ewkcoffea.modules.yield_tools as yt

SMALL = 0.000001

# Global variables
PRECISION = 6   # Decimal point precision in the text datacard output

#JERC Syst List
jerc_list = [
"AbsoluteMPFBias_correlated","AbsoluteScale_correlated","FlavorQCD_correlated","Fragmentation_correlated","PileUpDataMC_correlated",
"PileUpPtBB_correlated","PileUpPtEC1_correlated","PileUpPtEC2_correlated","PileUpPtHF_correlated","PileUpPtRef_correlated",
"RelativeFSR_correlated","RelativeJERHF_correlated","RelativePtBB_correlated","RelativePtHF_correlated","RelativeBal_correlated",
"SinglePionECAL_correlated","SinglePionHCAL_correlated"
]

# What each recognized year grouping consists of
ALL_YEARS_LST = ["UL16","UL16APV","UL17","UL18", "2022","2022EE", "2023","2023BPix"]

# Systs that are not correlated across years

SYSTS_SPECIAL = {

    "run2" : {
        "btagSFlight_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "btagSFbc_uncorrelated_2016APV"    : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "btagSFlight_uncorrelated_2016"    : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "btagSFbc_uncorrelated_2016"       : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "btagSFlight_uncorrelated_2017"    : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "btagSFbc_uncorrelated_2017"       : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "btagSFlight_uncorrelated_2018"    : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "btagSFbc_uncorrelated_2018"       : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "AbsoluteStat_uncorrelated_2016APV"   : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "RelativeJEREC1_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "RelativeJEREC2_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "RelativePtEC1_uncorrelated_2016APV"  : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "RelativePtEC2_uncorrelated_2016APV"  : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "TimePtEta_uncorrelated_2016APV"      : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "RelativeSample_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "RelativeStatEC_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "RelativeStatFSR_uncorrelated_2016APV": {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "RelativeStatHF_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "JER_2016APV"                         : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "AbsoluteStat_uncorrelated_2016"   : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "RelativeJEREC1_uncorrelated_2016" : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "RelativeJEREC2_uncorrelated_2016" : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "RelativePtEC1_uncorrelated_2016"  : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "RelativePtEC2_uncorrelated_2016"  : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "TimePtEta_uncorrelated_2016"      : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "RelativeSample_uncorrelated_2016" : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "RelativeStatEC_uncorrelated_2016" : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "RelativeStatFSR_uncorrelated_2016": {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "RelativeStatHF_uncorrelated_2016" : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "JER_2016"                         : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "AbsoluteStat_uncorrelated_2017"   : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "RelativeJEREC1_uncorrelated_2017" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "RelativeJEREC2_uncorrelated_2017" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "RelativePtEC1_uncorrelated_2017"  : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "RelativePtEC2_uncorrelated_2017"  : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "TimePtEta_uncorrelated_2017"      : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "RelativeSample_uncorrelated_2017" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "RelativeStatEC_uncorrelated_2017" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "RelativeStatFSR_uncorrelated_2017": {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "RelativeStatHF_uncorrelated_2017" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "JER_2017"                         : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "AbsoluteStat_uncorrelated_2018"   : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "RelativeJEREC1_uncorrelated_2018" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "RelativeJEREC2_uncorrelated_2018" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "RelativePtEC1_uncorrelated_2018"  : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "RelativePtEC2_uncorrelated_2018"  : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "TimePtEta_uncorrelated_2018"      : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "RelativeSample_uncorrelated_2018" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "RelativeStatEC_uncorrelated_2018" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "RelativeStatFSR_uncorrelated_2018": {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "RelativeStatHF_uncorrelated_2018" : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "JER_2018"                         : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
    },
    "total_run2" : {
        "btagSFlight_uncorrelated_2016APV" : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "btagSFbc_uncorrelated_2016APV"    : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "btagSFlight_uncorrelated_2016"    : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "btagSFbc_uncorrelated_2016"       : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "btagSFlight_uncorrelated_2017"    : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "btagSFbc_uncorrelated_2017"       : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "btagSFlight_uncorrelated_2018"    : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "btagSFbc_uncorrelated_2018"       : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "JER_2016APV"                      : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "JER_2016"                         : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "JER_2017"                         : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "JER_2018"                         : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},
        "JEC_2016APV"                      : {"yr_rel":"UL16APV", "yr_notrel": ["UL16", "UL17", "UL18"]},
        "JEC_2016"                         : {"yr_rel":"UL16", "yr_notrel": ["UL16APV", "UL17", "UL18"]},
        "JEC_2017"                         : {"yr_rel":"UL17", "yr_notrel": ["UL16APV", "UL16", "UL18"]},
        "JEC_2018"                         : {"yr_rel":"UL18", "yr_notrel": ["UL16APV", "UL16", "UL17"]},

    },

    "run3" : {
        "btagSFbc_uncorrelated_2022"       : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "btagSFbc_uncorrelated_2022EE"     : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "btagSFbc_uncorrelated_2023"       : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "btagSFbc_uncorrelated_2023BPix"   : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "AbsoluteStat_uncorrelated_2022"   : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "RelativeJEREC1_uncorrelated_2022" : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "RelativeJEREC2_uncorrelated_2022" : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "RelativePtEC1_uncorrelated_2022"  : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "RelativePtEC2_uncorrelated_2022"  : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "TimePtEta_uncorrelated_2022"      : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "RelativeSample_uncorrelated_2022" : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "RelativeStatEC_uncorrelated_2022" : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "RelativeStatFSR_uncorrelated_2022": {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "RelativeStatHF_uncorrelated_2022" : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "JER_2022"                         : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "AbsoluteStat_uncorrelated_2022EE"   : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "RelativeJEREC1_uncorrelated_2022EE" : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "RelativeJEREC2_uncorrelated_2022EE" : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "RelativePtEC1_uncorrelated_2022EE"  : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "RelativePtEC2_uncorrelated_2022EE"  : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "TimePtEta_uncorrelated_2022EE"      : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "RelativeSample_uncorrelated_2022EE" : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "RelativeStatEC_uncorrelated_2022EE" : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "RelativeStatFSR_uncorrelated_2022EE": {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "RelativeStatHF_uncorrelated_2022EE" : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "JER_2022EE"                         : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "AbsoluteStat_uncorrelated_2023"   : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "RelativeJEREC1_uncorrelated_2023" : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "RelativeJEREC2_uncorrelated_2023" : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "RelativePtEC1_uncorrelated_2023"  : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "RelativePtEC2_uncorrelated_2023"  : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "TimePtEta_uncorrelated_2023"      : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "RelativeSample_uncorrelated_2023" : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "RelativeStatEC_uncorrelated_2023" : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "RelativeStatFSR_uncorrelated_2023": {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "RelativeStatHF_uncorrelated_2023" : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "JER_2023"                         : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "AbsoluteStat_uncorrelated_2023BPix"   : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "RelativeJEREC1_uncorrelated_2023BPix" : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "RelativeJEREC2_uncorrelated_2023BPix" : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "RelativePtEC1_uncorrelated_2023BPix"  : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "RelativePtEC2_uncorrelated_2023BPix"  : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "TimePtEta_uncorrelated_2023BPix"      : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "RelativeSample_uncorrelated_2023BPix" : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "RelativeStatEC_uncorrelated_2023BPix" : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "RelativeStatFSR_uncorrelated_2023BPix": {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "RelativeStatHF_uncorrelated_2023BPix" : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "JER_2023BPix"                         : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
    },
    "total_run3" : {
        "btagSFbc_uncorrelated_2022"       : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "btagSFbc_uncorrelated_2022EE"     : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "btagSFbc_uncorrelated_2023"       : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "btagSFbc_uncorrelated_2023BPix"   : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "JER_2022"                         : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "JER_2022EE"                         : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "JER_2023"                         : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "JER_2023BPix"                         : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
        "JEC_2022"                         : {"yr_rel":"2022", "yr_notrel": ["2022EE","2023","2023BPix"]},
        "JEC_2022EE"                         : {"yr_rel":"2022EE", "yr_notrel": ["2022","2023","2023BPix"]},
        "JEC_2023"                         : {"yr_rel":"2023", "yr_notrel": ["2022","2022EE","2023BPix"]},
        "JEC_2023BPix"                         : {"yr_rel":"2023BPix", "yr_notrel": ["2022","2022EE","2023"]},
    },

    "y22" : {
        "btagSFbc_uncorrelated_2022"       : {"yr_rel":"2022", "yr_notrel": ["2022EE"]},
        "btagSFbc_uncorrelated_2022EE"     : {"yr_rel":"2022EE", "yr_notrel": ["2022"]},
    },

    "y23" : {
        "btagSFbc_uncorrelated_2023"       : {"yr_rel":"2023", "yr_notrel": ["2023BPix"]},
        "btagSFbc_uncorrelated_2023BPix"   : {"yr_rel":"2023BPix", "yr_notrel": ["2023"]},
    },

}


# Hard code the rateParam lines to put at the end of the card (for background normalization)
RATE_PARAM_LINES = [
    "ZZ_norm rateParam * ZZ 1 [0,5]",
    "txZ_norm rateParam * ttZ 1 [0,5]",
    "txZ_norm rateParam * tWZ 1 [0,5]",
]


########### Writing the datacard ###########

# Make the datacard for a given channel
def make_ch_card(ch,proc_order,year_name,ch_ylds,ch_kappas=None,ch_gmn=None,extra_lines=None,out_dir="."):

    # Building blocks we'll need to build the card formatting
    bin_str = f"bin_{ch}"
    syst_width = 0
    col_width = max(PRECISION*2+5,len(bin_str))
    line_break = "##----------------------------------\n"
    left_width = len(line_break) + 2
    left_width = max(syst_width+len("shape")+1,left_width)

    # The output name, location
    outf_card_name = f"wwz4l_card_{ch}_{year_name}.txt"
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
            if any([x in p for x in sg.SIG_LST]): # Check for if the process is signal or not
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

        # Write any extra lines
        if extra_lines is not None:
            for extra_line in extra_lines:
                f.write(f"{extra_line}\n")



########### Getting and manipulating yields ###########


# Modify the yields dict to properly calculate the per-year systs
#   - Needs all years in input dict, expects input dict with keys: year,cat,syst,proc
#   - Modifies the "FR2" subdict of the yield dict to have proper up and down variations for the per year systematics
#   - Because of how we fill in the processor, the yields for per year systs come _only_ from that year
#   - So this function adds the nominal yields from the other three years to the up/down variation for the relevant year
#   - Note the in_dict is modifed in place (we do not return a copy of the dict)
def handle_per_year_systs_for_fr(in_dict,year_name,do_tot):
    if do_tot:
        systs_special = SYSTS_SPECIAL[f"total_{year_name}"]
    else:
        systs_special = SYSTS_SPECIAL[year_name]
    for cat in in_dict["FR"].keys():
        for sys in systs_special:
            # Find up/down variation for the year relevant to that syst
            yrrel = systs_special[sys]["yr_rel"] # The relevant year for this special syst
            yld_yrrel_up = in_dict[yrrel][cat][f"{sys}Up"]
            yld_yrrel_do = in_dict[yrrel][cat][f"{sys}Down"]
            # Get nominal yld for all years other than the relevant one
            yld_yrrel_nom = in_dict[yrrel][cat]["nominal"]
            yld_yrall_nom = in_dict["FR"][cat]["nominal"]
            yld_allbutyrrel_nom = utils.get_diff_between_dicts(yld_yrall_nom,yld_yrrel_nom,difftype="absolute_diff") # This is: x = sum(nom yld for all years except relevant year)
            # Get the yield with just the relevant year's up/down variation varied (with nominal yld for all other years)
            yld_up = utils.get_diff_between_dicts(yld_allbutyrrel_nom, yld_yrrel_up, difftype="sum") # This is: x + (up   yld for relevant year)
            yld_do = utils.get_diff_between_dicts(yld_allbutyrrel_nom, yld_yrrel_do, difftype="sum") # This is: x + (down yld for relevant year)
            in_dict["FR"][cat][f"{sys}Up"] = yld_up
            in_dict["FR"][cat][f"{sys}Down"] = yld_do


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
            if val <= 0:
                print(f"WARNING: Process \"{proc}\" in cat \"{cat}\" is negative ({val}), replacing with {SMALL} and shifting up/down systematic variations accordingly.")
                out_dict[cat]["nominal"][proc][0] = SMALL
                out_dict[cat]["nominal"][proc][1] = (abs(val) + np.sqrt(var))**2
                for syst in out_dict[cat]:
                    if syst == "nominal": continue # Already handled this one
                    syst_var_orig = out_dict[cat][syst][proc][0] # Dont bother messsing with mc stat error on the syst variation
                    out_dict[cat][syst][proc][0] = (syst_var_orig - val) + SMALL # Center around SMALL

    return out_dict



########### Regarding uncertainties ###########

# Get the rate systs from the intput json, dump into dict (with nested keys: syst, proc)
# Outputs strings, ready to be dumped into datacard
def get_rate_systs(proc_lst,year_tag):
    syst_json = ewkcoffea_path("params/rate_systs.json")
    with open(syst_json) as f_systs: rate_systs_dict = json.load(f_systs)

    # Build up the dictionary
    out_dict = {}

    # Make row for rate uncertainties that impact all processes
    for uncty_name in rate_systs_dict["rate_uncertainties_all_proc"]:
        out_dict[uncty_name] = {}
        for proc in proc_lst:
            if uncty_name == "lumi":
                # The lumi uncty is different depending on year
                out_dict[uncty_name][proc] = str(rate_systs_dict["rate_uncertainties_all_proc"][uncty_name][year_tag])
            else:
                out_dict[uncty_name][proc] = str(rate_systs_dict["rate_uncertainties_all_proc"][uncty_name])

    # Make rows for rate uncertainties that impact a subset of the processes
    for uncty_name in rate_systs_dict["rate_uncertainties_some_proc"]:

        # If the uncty has a decorrelation_dict defined, use that to append a tag to the name of the uncty
        # E.g., for the fake uncertianty, we append the run so that it is fake_run2 or fake_run3
        tag = ""
        if "decorrelation_dict" in rate_systs_dict["rate_uncertainties_some_proc"][uncty_name]:
            tag = "_" + rate_systs_dict["rate_uncertainties_some_proc"][uncty_name]["decorrelation_dict"][year_tag]

        # Loop over procs and put in the uncty or "-" if NA
        for proc_of_interest in rate_systs_dict["rate_uncertainties_some_proc"][uncty_name]["procs"]:
            out_dict[f"{uncty_name}_{proc_of_interest}{tag}"] = {}
            for proc_itr in proc_lst:
                if proc_itr == proc_of_interest:
                    uncty = str(rate_systs_dict["rate_uncertainties_some_proc"][uncty_name]["val"])
                else:
                    uncty = "-"
                out_dict[f"{uncty_name}_{proc_of_interest}{tag}"][proc_itr] = uncty

    return out_dict


#Determines if the up and down variations of a systematic are in the same direction
def determine_updo_same(nom,up,down):
    if nom < 0:
        raise Exception("Negative values should have been fixed by this point!")
    elif ((up > nom) and (down > nom)):
        return True
    elif ((up < nom) and (down < nom)):
        return True
    else: return False

#Fixes the situation when the up and down variation are in the same direction by taking the larger variation and symmetrize
def fix_updown_same(nom,up,down):
    diff_1 = abs(nom - up)
    diff_2 = abs(nom - down)
    diff = max(diff_1,diff_2)
    kappa_up = (nom + diff)/nom
    kappa_down = (nom - diff)/nom
    return kappa_up, kappa_down

# Get kappa dict (e.g. up/nom ratios) from the dict of all histograms
def get_kappa_dict(in_dict_mc,in_dict_data,yrs_lst):

    # Get the list of systematic base names (i.e. without the up and down tags)
    #     - Assumes each syst has a "systnameUp" and a "systnameDown"
    #     - This will drop nominal (since there is no "nominalUp" to tag on)
    def get_syst_base_name_lst(in_lst):
        out_lst = []
        for syst in in_lst:
            if syst.endswith("Up"):
                syst_name_base = syst[:-2]
                if syst_name_base not in out_lst:
                    out_lst.append(syst_name_base)
        return out_lst

    # Get the kappa dict
    kappa_dict = {}
    for cat in in_dict_mc.keys():
        kappa_dict[cat] = {}
        for sys in get_syst_base_name_lst(list(in_dict_mc[cat].keys())):

            # Skip year specific systs that are not relevant
            # E.g. if the pkl file is all r3, but we're only making a 22 datacard, skip btagSFbc_uncorrelated_2023BPix variation
            # Do this by checking if the syst name ends with a year, and then if that year is in our list for this card
            if (sys.split("_")[-1] in ALL_YEARS_LST) and (sys.split("_")[-1] not in yrs_lst): continue

            kappa_dict[cat][sys] = {}
            for proc in in_dict_mc[cat]["nominal"]:

                kappa_dict[cat][sys][proc] = {}
                valvar_up = in_dict_mc[cat][f"{sys}Up"][proc]
                valvar_do = in_dict_mc[cat][f"{sys}Down"][proc]
                valvar_nom = in_dict_mc[cat]["nominal"][proc]
                valvar_kappa_up = yt.valvar_op(valvar_up,valvar_nom,"div")
                valvar_kappa_do = yt.valvar_op(valvar_do,valvar_nom,"div")

                # Handle negative cases
                if (valvar_kappa_up[0]<=0) and (valvar_kappa_do[0]<=0):

#                    if (set(yrs_lst) == set(["2022","2022EE","2023","2023BPix"])) and (cat=="sr_4l_bdt_of_1") and (sys=="FSR") and (proc=="ttZ"):
#                        # One known case of this: FSR for ttZ for sr_4l_bdt_of_1 for r3
#                        #   - Does not matter since we don't use this bin, so skip it
#                        #   - Super super hard coded skip for this case:
#                    else:
#                        # Otherwise raise an error, we should stop and take a look at what's happening
#                        raise Exception(f"Both Kappas Neagtive for process: {proc}, category: {cat}, systematic: {sys}")

#                    raise Exception("Kappas in same direction, but we don't care right now.")

                    valvar_kappa_up[0], valvar_kappa_do[0] = fix_updown_same(valvar_nom[0],valvar_up[0],valvar_do[0])


                if valvar_kappa_up[0] <= 0:
                    print(f"WARNING: Up var for {sys} for {proc} for {cat} is negative, setting to {SMALL}.")
                    valvar_kappa_up[0] = SMALL
                if valvar_kappa_do[0] <= 0:
                    valvar_kappa_do[0] = SMALL
                    print(f"WARNING: Down var for {sys} for {proc} for {cat} is negative, setting to {SMALL}.")

                kappa_dict[cat][sys][proc]["Up"] = valvar_kappa_up
                kappa_dict[cat][sys][proc]["Down"] = valvar_kappa_do


    return kappa_dict


# Get the MC stats from the yld dict and put it into kappa dict
# Can optionally skip processes (e.g. skip data driven ZZ and ttZ since do these different via gmN)
def add_stats_kappas(yld_mc, kappas, skip_procs=[]):

    kappas_out = copy.deepcopy(kappas)

    # Loop over cat and get mc stats
    for cat in yld_mc:

        # Looping over each proc (these will be rows, i.e. one row for each proc's mc stats)
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
                    if do <= 0:
                        print(f"WARNING: For cat \"{cat}\" and proc \"{proc_of_interest}\", the uncertainty {np.sqrt(valvar[1])} is larger than the value {valvar[0]}. Clipping down variation to {SMALL}.")
                        do = SMALL
                    kappas_out[cat][f"stats_{cat}_{proc_of_interest}"][proc_itr]["Up"]   = [up, None] # Rel err up, do not include error on the error (just leave as None)
                    kappas_out[cat][f"stats_{cat}_{proc_of_interest}"][proc_itr]["Down"] = [do, None] # Rel err down, do not include error on the error (just leave as None)

    return kappas_out



########### Put stuff into form to pass to the function to write out cards ###########

# Get just the numbers we want for rate row for datacard
# Also sum all MC rates together into asimov number if not unblind
# Assumes in_dict has nested keys: cat,syst,proc
def get_rate_for_dc(in_dict_mc,in_dict_data,cat,unblind):
    out_dict = {}
    asimov_data = 0
    for proc in in_dict_mc[cat]["nominal"]:
        rate = in_dict_mc[cat]["nominal"][proc][0]
        if rate < 0:
            print(f"\nWarning: Process \"{proc}\" in \"{cat}\" has negative total rate: {rate}.\n")
            raise Exception("This should not be happening. Exiting.")
        out_dict[proc] = str(rate)
        asimov_data += rate

    if not unblind:
        out_dict["data_obs"] = str(asimov_data)
    else:
        out_dict["data_obs"] = str(in_dict_data[cat]["nominal"]["data"][0])

    return out_dict


# Assumes in_dict has nested keys: cat,syst_base,proc,"Up"(or "Down)
# Takes just the val (dropps uncty)
def get_kappa_for_dc(in_dict,cat):
    out_dict = {}
    for systname_base in in_dict[cat]:
        out_dict[systname_base] = {}
        for proc in in_dict[cat][systname_base]:
            u = in_dict[cat][systname_base][proc]['Up'][0]
            d = in_dict[cat][systname_base][proc]['Down'][0]
            if (u is None) and (d is None):
                out_dict[systname_base][proc] = "-"
            else:
                out_dict[systname_base][proc] = f"{d}/{u}"
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
    parser.add_argument("--out-dir","-d",default="./cards_wwz4l",help="Output directory to write root and text datacard files to")
    parser.add_argument("-s","--do-nuisance",action="store_true",help="Include nuisance parameters")
    parser.add_argument("--do-tf",action="store_true",help="Do the TF data-driven background estimation")
    parser.add_argument("--bdt",action="store_true",help="Use BDT SR bins")
    parser.add_argument("--jec-tot",action="store_true",help="Use Total JEC Uncertainty")
    parser.add_argument("--unblind",action="store_true",help="If set, use real data, otherwise use asimov data")
    parser.add_argument('-u', "--run", default='run2', help = "Which years to process", choices=["run2","run3","y22","y23"])

    args = parser.parse_args()
    in_file = args.in_file_name
    out_dir = args.out_dir
    do_nuis = args.do_nuisance
    do_tf   = args.do_tf
    jec_tot = args.jec_tot
    use_bdt_sr = args.bdt
    unblind = args.unblind
    run = args.run

    # Check args
    if out_dir != "." and not os.path.exists(out_dir):
        print(f"Making dir \"{out_dir}\"")
        os.makedirs(out_dir)

    # Set list of years from the run name
    if run == "run2": yrs_lst = ["UL16APV","UL16","UL17","UL18"]
    elif run == "run3": yrs_lst = ["2022","2022EE","2023","2023BPix"]
    elif run == "y22" : yrs_lst = ["2022","2022EE"]
    elif run == "y23" : yrs_lst = ["2023","2023BPix"]
    else: raise Exception("Unknown year")

    # Get the histo
    f = pickle.load(gzip.open(in_file))
    histo = f["njets"] # Let's use njets

    # Get the dictionary defining the mc sample grouping
    sample_names_dict_data = {"FR" : sg.create_data_sample_dict(run)}
    sample_names_dict_mc   = {"FR" : sg.create_mc_sample_dict(run)}
    for year_indiv_name in yrs_lst:
        sample_names_dict_mc[year_indiv_name] = sg.create_mc_sample_dict(year_indiv_name)

    # Get yield dictionary (nested in the order: year,cat,syst,proc)
    yld_dict_mc_allyears = {}
    for year in sample_names_dict_mc:
        yld_dict_mc_allyears[year] = yt.get_yields(histo,sample_names_dict_mc[year])
    if do_nuis:
        handle_per_year_systs_for_fr(yld_dict_mc_allyears,run,jec_tot)

    yld_dict_mc = yld_dict_mc_allyears["FR"]
    yld_dict_data = yt.get_yields(histo,sample_names_dict_data["FR"])

    # Scale yield for any processes (e.g. for testing impacts of small backgrounds)
    scale_dict = {"WZ":1.0}
    yld_dict_mc = yt.scale_yld_dict(yld_dict_mc,scale_dict)

    ####################################################################################
    # Dump some info about a bin (just raw numbers, more or less)
    # This print is before we start messing with the yields (eg to get rid of negatives)
    printinfo = 0
    if printinfo:
        s = "FSR"
        p = "ttZ"
        c = "sr_4l_bdt_of_1"
        cr = "cr_4l_btag_of"
        print(p,c,s,cr)
        print("\nPrinting info:")
        print("mc sr n",yld_dict_mc[c]["nominal"][p])
        print("mc sr u",yld_dict_mc[c][f"{s}Up"][p])
        print("mc sr d",yld_dict_mc[c][f"{s}Down"][p])
        print("mc cr n",yld_dict_mc[cr]["nominal"][p])
        print("mc cr u",yld_dict_mc[cr][f"{s}Up"][p])
        print("da cr n",yld_dict_data[cr]["nominal"]["data"])
        print("da cr u",yld_dict_data[cr][f"{s}Up"]["data"])
        print("\n")
        exit()
    ####################################################################################


    # Get rid of negative yields (and recenter syst variations around SMALL), should happen before computing kappas
    yld_dict_mc = handle_negatives(yld_dict_mc)

    # Get the syst ratios to nominal (i.e. kappas)
    kappa_dict = None
    if do_nuis:
        kappa_dict = get_kappa_dict(yld_dict_mc,yld_dict_data,yrs_lst)
        # Don't do mc stats kappas for data-driven bkg if doing TFs
        if do_tf: skip_stats_kappas_lst = ["ZZ","ttZ"]
        else: skip_stats_kappas_lst = []
        kappa_dict = add_stats_kappas(yld_dict_mc,kappa_dict,skip_procs=skip_stats_kappas_lst)

    # Do the TF calculation
    if do_tf:
        yld_dict_mc, kappa_dict, gmn_dict = yt.do_tf(yld_dict_mc,yld_dict_data,kappa_dict,sg.BKG_TF_MAP)


    #### Make the cards for each channel ####

    # Get list of channels
    cat_lst_cr = ["cr_4l_btag_of_1b", "cr_4l_btag_of_2b", "cr_4l_btag_sf_offZ_met80_1b", "cr_4l_btag_sf_offZ_met80_2b","cr_4l_sf"]
    cat_lst_sr = sg.CAT_LST_CB
    if use_bdt_sr:
        if run in ["run2"]:
            cat_lst_sr = sg.CAT_LST_BDT
        elif run in ["run3", "y22", "y23"]:
            cat_lst_sr = sg.CAT_LST_BDT_COARSE
        else:
            print(run)
            raise Exception("Unknown year")
    cat_lst = cat_lst_sr + cat_lst_cr
    print(f"\nMaking cards for {cat_lst}. \nPutting in {out_dir}.")

    # Loop over channels and make cards
    for ch in cat_lst:

        # Use real data in CRs
        if ch in cat_lst_cr: unblind = True

        # Get just the info we want to put in the card in str form
        rate_for_dc_ch = get_rate_for_dc(yld_dict_mc,yld_dict_data,ch,unblind)

        # Get the kappa and gamma dict for this channel if we are doing systs
        kappa_for_dc_ch = None
        gmn_for_dc_ch   = None
        if do_nuis:
            kappa_for_dc_ch = get_kappa_for_dc(kappa_dict,ch)
            kappa_for_dc_ch.update(get_rate_systs(sg.PROC_LST,run)) # Append in the ones from rate json
        if do_nuis and do_tf and (ch not in cat_lst_cr):
            # TF calculation not meaningful for CRs
            gmn_for_dc_ch = get_gmn_for_dc(gmn_dict[ch],proc_lst=sg.PROC_LST)


        # Make the card for this chan
        make_ch_card(
            ch,
            sg.PROC_LST,
            run,
            rate_for_dc_ch,
            kappa_for_dc_ch,
            gmn_for_dc_ch,
            extra_lines=RATE_PARAM_LINES,
            out_dir=out_dir,
        )


    print("Finished!")


if __name__ == "__main__":
    main()
