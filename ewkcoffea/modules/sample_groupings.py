# Dict to keep track of sample groupings, and associated tools


######################## Lists ########################

PROC_LST = ["WWZ","ZH","ZZ","ttZ","tWZ","WZ","other"]
SIG_LST = ["WWZ","ZH"]
BKG_LST = ["ZZ","ttZ","tWZ","WZ","other"]

SR_SF_CB = ["sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C"]
SR_OF_CB = ["sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]
CAT_LST_CB = SR_SF_CB + SR_OF_CB

SR_SF_BDT = ["sr_4l_bdt_sf_1", "sr_4l_bdt_sf_2", "sr_4l_bdt_sf_3", "sr_4l_bdt_sf_4", "sr_4l_bdt_sf_5", "sr_4l_bdt_sf_6", "sr_4l_bdt_sf_7", "sr_4l_bdt_sf_8"]
SR_OF_BDT = ["sr_4l_bdt_of_1", "sr_4l_bdt_of_2", "sr_4l_bdt_of_3", "sr_4l_bdt_of_4", "sr_4l_bdt_of_5", "sr_4l_bdt_of_6", "sr_4l_bdt_of_7", "sr_4l_bdt_of_8"]
SR_SF_BDT_COARSE = ["sr_4l_bdt_sf_coarse_1", "sr_4l_bdt_sf_coarse_2", "sr_4l_bdt_sf_coarse_3", "sr_4l_bdt_sf_coarse_4"]
SR_OF_BDT_COARSE = ["sr_4l_bdt_of_coarse_1", "sr_4l_bdt_of_coarse_2", "sr_4l_bdt_of_coarse_3", "sr_4l_bdt_of_coarse_4"]

CAT_LST_BDT = SR_SF_BDT + SR_OF_BDT
CAT_LST_BDT_COARSE = SR_SF_BDT_COARSE + SR_OF_BDT_COARSE

CAT_LST_CR = ["cr_4l_btag_of", "cr_4l_btag_sf_offZ_met80", "cr_4l_sf"]

######################## Dictionaries ########################

# Map showing which CR categories are used to calculate NSFs for which SRs for which processes
BKG_TF_MAP = {

    "ZZ" : {

        # Cut based
        "sr_4l_sf_A" : "cr_4l_sf",
        "sr_4l_sf_B" : "cr_4l_sf",
        "sr_4l_sf_C" : "cr_4l_sf",
        "sr_4l_of_1" : "cr_4l_sf",
        "sr_4l_of_2" : "cr_4l_sf",
        "sr_4l_of_3" : "cr_4l_sf",
        "sr_4l_of_4" : "cr_4l_sf",


        # BDT based
        "sr_4l_bdt_sf_1" : "cr_4l_sf",
        "sr_4l_bdt_sf_2" : "cr_4l_sf",
        "sr_4l_bdt_sf_3" : "cr_4l_sf",
        "sr_4l_bdt_sf_4" : "cr_4l_sf",
        "sr_4l_bdt_sf_5" : "cr_4l_sf",
        "sr_4l_bdt_sf_6" : "cr_4l_sf",
        "sr_4l_bdt_sf_7" : "cr_4l_sf",
        "sr_4l_bdt_sf_8" : "cr_4l_sf",

        "sr_4l_bdt_sf_coarse_1" : "cr_4l_sf",
        "sr_4l_bdt_sf_coarse_2" : "cr_4l_sf",
        "sr_4l_bdt_sf_coarse_3" : "cr_4l_sf",
        "sr_4l_bdt_sf_coarse_4" : "cr_4l_sf",

        "sr_4l_bdt_of_1" : "cr_4l_sf",
        "sr_4l_bdt_of_2" : "cr_4l_sf",
        "sr_4l_bdt_of_3" : "cr_4l_sf",
        "sr_4l_bdt_of_4" : "cr_4l_sf",
        "sr_4l_bdt_of_5" : "cr_4l_sf",
        "sr_4l_bdt_of_6" : "cr_4l_sf",
        "sr_4l_bdt_of_7" : "cr_4l_sf",
        "sr_4l_bdt_of_8" : "cr_4l_sf",

        "sr_4l_bdt_of_coarse_1" : "cr_4l_sf",
        "sr_4l_bdt_of_coarse_2" : "cr_4l_sf",
        "sr_4l_bdt_of_coarse_3" : "cr_4l_sf",
        "sr_4l_bdt_of_coarse_4" : "cr_4l_sf",

    },
    "ttZ" : {

        # Cut based
        "sr_4l_sf_A" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_sf_B" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_sf_C" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_of_1" : "cr_4l_btag_of",
        "sr_4l_of_2" : "cr_4l_btag_of",
        "sr_4l_of_3" : "cr_4l_btag_of",
        "sr_4l_of_4" : "cr_4l_btag_of",

        # BDT based
        "sr_4l_bdt_sf_1" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_2" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_3" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_4" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_5" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_6" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_7" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_8" : "cr_4l_btag_sf_offZ_met80",

        "sr_4l_bdt_sf_coarse_1" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_coarse_2" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_coarse_3" : "cr_4l_btag_sf_offZ_met80",
        "sr_4l_bdt_sf_coarse_4" : "cr_4l_btag_sf_offZ_met80",

        "sr_4l_bdt_of_1" : "cr_4l_btag_of",
        "sr_4l_bdt_of_2" : "cr_4l_btag_of",
        "sr_4l_bdt_of_3" : "cr_4l_btag_of",
        "sr_4l_bdt_of_4" : "cr_4l_btag_of",
        "sr_4l_bdt_of_5" : "cr_4l_btag_of",
        "sr_4l_bdt_of_6" : "cr_4l_btag_of",
        "sr_4l_bdt_of_7" : "cr_4l_btag_of",
        "sr_4l_bdt_of_8" : "cr_4l_btag_of",

        "sr_4l_bdt_of_coarse_1" : "cr_4l_btag_of",
        "sr_4l_bdt_of_coarse_2" : "cr_4l_btag_of",
        "sr_4l_bdt_of_coarse_3" : "cr_4l_btag_of",
        "sr_4l_bdt_of_coarse_4" : "cr_4l_btag_of",
    }
}


# The "official" Run 2 groupings
SAMPLE_DICT_BASE_RUN2 = {
    "WWZ" : ["WWZJetsTo4L2Nu"],
    "ZH"  : ["GluGluZH","qqToZHToZTo2L"],

    #"qqZZ": ["ZZTo4l"],
    #"ggZZ": ["ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau"],
    "ZZ"  : ["ZZTo4l", "ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau"],

    "ttZ" : [
        "TTZToLL_M_1to10",
        "TTZToLLNuNu_M_10",
        "TTZToQQ",
    ],

    "tWZ" : ["tWZ4l"], # Old sample tWll

    #"WZ" : ["WZTo3LNu"], # Added to 'other' for now

    "other" : [

        ##"WWZJetsTo4L2Nu",
        ##"GluGluZH","qqToZHToZTo2L",
        ##"ZZTo4l", "ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau",
        ##"TTZToLL_M_1to10","TTZToLLNuNu_M_10","TTZToQQ",
        ##"tWZ4l",

        #"DYJetsToLL_M_10to50_MLM", # No events in SR
        #"DYJetsToLL_M_50_MLM", # No events in SR
        "SSWW",
        "ST_antitop_t-channel",
        "ST_top_s-channel",
        "ST_top_t-channel",
        "tbarW_noFullHad",
        "ttHnobb",
        "TTTo2L2Nu",
        "TTWJetsToLNu",
        "TTWJetsToQQ",
        "tW_noFullHad",
        "tZq",
        "VHnobb",
        ##"WJetsToLNu",
        "WWTo2L2Nu",
        #"WZTo3LNu", # Now by itself

        "WWW",
        "WZZ",
        "ZZZ",

        "ggHToZZ4L",
        "WZTo3LNu",
    ],
}

# The Run3 groupings
SAMPLE_DICT_BASE_RUN3 = {
    "WWZ" : ["WWZJetsTo4L2Nu"],
    "ZH"  : ["GluGluZHTo2WTo2L2Nu","qqToZHTo2WTo2L2Nu"],
    "ZZ"  : ["ZZTo4l", "ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau"],
    "ttZ" : [
        "TTZToLL_M_4to50",
        "TTZToLL_M_50",
        "TTZToQQ",
    ],
    "tWZ" : ["tWZ4l"],
    #"WZ" : ["WZTo3LNu"], # Added to 'other'
    "other" : [
        #"DYJetsToLL_M_10to50_MLM",
        #"DYJetsToLL_M_50_MLM",
        "SSWW_TT",
        "SSWW_TL",
        "SSWW_LL",
        "ST_antitop_t-channel",
        "ST_top_s-channel",
        "ST_top_t-channel",
        "tbarW_leptonic",
        "tbarW_semileptonic",
        "ttHnobb",
        "TTTo2L2Nu",
        "tW_leptonic",
        "tW_semileptonic",
        "tZq",
        "VHnobb",
        "WJetsToLNu",
        "WWTo2L2Nu",
        "WWW",
        "WZZ",
        "ZZZ",
        "ggHToZZ4L",
        "WZTo3LNu", # added to other for now
    ],
}

# Processes indiviudally
SAMPLE_DICT_BASE_INDIV_RUN2 = {
    "WWZJetsTo4L2Nu":            ["WWZJetsTo4L2Nu"],
    "GluGluZH":                  ["GluGluZH"],
    "qqToZHToZTo2L":             ["qqToZHToZTo2L"],
    "ZZTo4l":                    ["ZZTo4l"],
    "ggToZZTo2e2mu":             ["ggToZZTo2e2mu"],
    "ggToZZTo2e2tau":            ["ggToZZTo2e2tau"],
    "ggToZZTo2mu2tau":           ["ggToZZTo2mu2tau"],
    "ggToZZTo4e":                ["ggToZZTo4e"],
    "ggToZZTo4mu":               ["ggToZZTo4mu"],
    "ggToZZTo4tau":              ["ggToZZTo4tau"],
    "TTZToLL_M_1to10":           ["TTZToLL_M_1to10"],
    "TTZToLLNuNu_M_10":          ["TTZToLLNuNu_M_10"],
    "TTZToQQ":                   ["TTZToQQ"],
    "tWZ4l" :                    ["tWZ4l"], # Old sample tWll

    #"DYJetsToLL_M_10to50_MLM": ["DYJetsToLL_M_10to50_MLM"],
    #"DYJetsToLL_M_50_MLM":     ["DYJetsToLL_M_50_MLM"],
    "SSWW":                      ["SSWW"],
    "ST_antitop_t-channel":      ["ST_antitop_t-channel"],
    "ST_top_s-channel":          ["ST_top_s-channel"],
    "ST_top_t-channel":          ["ST_top_t-channel"],
    "tbarW_noFullHad":           ["tbarW_noFullHad"],
    "ttHnobb":                   ["ttHnobb"],
    "TTTo2L2Nu":                 ["TTTo2L2Nu"],
    "TTWJetsToLNu":              ["TTWJetsToLNu"],
    "TTWJetsToQQ":               ["TTWJetsToQQ"],
    "tW_noFullHad":              ["tW_noFullHad"],
    "tZq":                       ["tZq"],
    "VHnobb":                    ["VHnobb"],
    ##"WJetsToLNu":              ["WJetsToLNu"],
    "WWTo2L2Nu":                 ["WWTo2L2Nu"],
    "WZTo3LNu":                  ["WZTo3LNu"],

    "WWW" : ["WWW"],
    "WZZ" : ["WZZ"],
    "ZZZ" : ["ZZZ"],

    "ggHToZZ4L" : ["ggHToZZ4L"],
}

SAMPLE_DICT_BASE_INDIV_RUN3 = {
    "WWZJetsTo4L2Nu"          : ["WWZJetsTo4L2Nu"],
    "GluGluZHTo2WTo2L2Nu"     : ["GluGluZHTo2WTo2L2Nu"],
    "qqToZHTo2WTo2L2Nu"       : ["qqToZHTo2WTo2L2Nu"],
    "ZZTo4l"                  : ["ZZTo4l"],
    "ggToZZTo2e2mu"           : ["ggToZZTo2e2mu"],
    "ggToZZTo2e2tau"          : ["ggToZZTo2e2tau"],
    "ggToZZTo2mu2tau"         : ["ggToZZTo2mu2tau"],
    "ggToZZTo4e"              : ["ggToZZTo4e"],
    "ggToZZTo4mu"             : ["ggToZZTo4mu"],
    "ggToZZTo4tau"            : ["ggToZZTo4tau"],
    "TTZToLL_M_4to50"         : ["TTZToLL_M_4to50"],
    "TTZToLL_M_50"            : ["TTZToLL_M_50"],
    "TTZToQQ"                 : ["TTZToQQ"],
    "tWZ4l"                   : ["tWZ4l"],
    "WZTo3LNu"                : ["WZTo3LNu"],
    #"DYJetsToLL_M_10to50_MLM" : ["DYJetsToLL_M_10to50_MLM"],
    #"DYJetsToLL_M_50_MLM"     : ["DYJetsToLL_M_50_MLM"],
    "SSWW_TT"                 : ["SSWW_TT"],
    "SSWW_TL"                 : ["SSWW_TL"],
    "SSWW_LL"                 : ["SSWW_LL"],
    "ST_antitop_t-channel"    : ["ST_antitop_t-channel"],
    "ST_top_s-channel"        : ["ST_top_s-channel"],
    "ST_top_t-channel"        : ["ST_top_t-channel"],
    "tbarW_leptonic"          : ["tbarW_leptonic"],
    "tbarW_semileptonic"      : ["tbarW_semileptonic"],
    "ttHnobb"                 : ["ttHnobb"],
    "TTTo2L2Nu"               : ["TTTo2L2Nu"],
    "tW_leptonic"             : ["tW_leptonic"],
    "tW_semileptonic"         : ["tW_semileptonic"],
    "tZq"                     : ["tZq"],
    "VHnobb"                  : ["VHnobb"],
    "WJetsToLNu"              : ["WJetsToLNu"],
    "WWTo2L2Nu"               : ["WWTo2L2Nu"],
    "WWW"                     : ["WWW"],
    "WZZ"                     : ["WZZ"],
    "ZZZ"                     : ["ZZZ"],
    "ggHToZZ4L"               : ["ggHToZZ4L"],
}

######################## Tools ########################

# Pass dictionary with the base names for the samples, and return with full list for the years
#     - This is a wrapper around create_mc_sample_dict_single_run
#     - This wrapper handles the stapling together of R2 and R3 dictionaries if necessary)
def create_mc_sample_dict(year,yld_individual=False):
    # If the year is "all" (i.e. run2+run3) get the dict for each, and combine
    if year == "all":
        out_dict = {}
        out_dict_r2 = create_mc_sample_dict_single_run("run2",yld_individual)
        out_dict_r3 = create_mc_sample_dict_single_run("run3",yld_individual)
        if set(out_dict_r2.keys()) != set(out_dict_r3.keys()):
            raise Exception("Run 2 and Run 3 are not combinable, they do not have the same keys.")
        for proc_group_key in out_dict_r2.keys():
            out_dict[proc_group_key] = out_dict_r2[proc_group_key] + out_dict_r3[proc_group_key]
    # Otherwise just call create_mc_sample_dict_single_run directly
    else:
        out_dict = create_mc_sample_dict_single_run(year,yld_individual)
    return out_dict

# Pass dictionary with the base names for the samples, and return with full list for the years
# Does not create combined R2 and R3 dict since the sample names are different
# To make combined, call this once for R2 and once for R2 and glue together the resuls
def create_mc_sample_dict_single_run(year,yld_individual=False):
    out_dict = {}
    all_years = ["UL16APV","UL16","UL17","UL18","2022","2022EE"]
    r2_years = ["UL16APV","UL16","UL17","UL18"]
    #r3_years = ["2022","2022EE","2023","2023BPix"] TODO: USE this when we switch to 2023 MC Samples
    r3_years = ["2022","2022EE", "2023","2023BPix"]
    y22_years = ["2022","2022EE"]
    y23_years = ["2023","2023BPix"]
    if year == "run2":
        years = r2_years
        sample_dict_base = SAMPLE_DICT_BASE_RUN2
        if yld_individual: sample_dict_base = SAMPLE_DICT_BASE_INDIV_RUN2 # If we want individual not grouped yields
    elif year == "run3":
        years = r3_years
        sample_dict_base = SAMPLE_DICT_BASE_RUN3
        if yld_individual: sample_dict_base = SAMPLE_DICT_BASE_INDIV_RUN3 # If we want individual not grouped yields
    elif year == "y22":
        years = y22_years
        sample_dict_base = SAMPLE_DICT_BASE_RUN3
        if yld_individual: sample_dict_base = SAMPLE_DICT_BASE_INDIV_RUN3 # If we want individual not grouped yields
    elif year == "y23":
        years = y23_years
        sample_dict_base = SAMPLE_DICT_BASE_RUN3
        if yld_individual: sample_dict_base = SAMPLE_DICT_BASE_INDIV_RUN3 # If we want individual not grouped yields
    else:
        years = [year]
        if year in r2_years:
            sample_dict_base = SAMPLE_DICT_BASE_RUN2
            if yld_individual: sample_dict_base = SAMPLE_DICT_BASE_INDIV_RUN2 # If we want individual not grouped yields
        elif year in r3_years:
            sample_dict_base = SAMPLE_DICT_BASE_RUN3
            if yld_individual: sample_dict_base = SAMPLE_DICT_BASE_INDIV_RUN3 # If we want individual not grouped yields
        else:
            raise Exception(f"ERROR: Unrecognized year \"{year}\". Exiting.")

    for proc_group in sample_dict_base.keys():
        out_dict[proc_group] = []
        for proc_base_name in sample_dict_base[proc_group]:
            for year_str in years:
                #ST_antitop_t-channel has no yields for 2023BPix. Simple method to skip this sample for particular year
                if (year_str == "2023BPix") and (proc_base_name == "ST_antitop_t-channel"):
                    continue
                else:
                    out_dict[proc_group].append(f"{year_str}_{proc_base_name}")

    return out_dict

# Get data sampel dict
def create_data_sample_dict(year):
    if year == "all":
        grouping_data = {'data': ["UL16APV_data","UL16_data","UL17_data","UL18_data","2022_data","2022EE_data"]}
    elif year == "run2":
        grouping_data = {'data': ["UL16APV_data","UL16_data","UL17_data","UL18_data"]}
    elif year == "run3":
        grouping_data = {'data': ["2022_data","2022EE_data","2023_data","2023BPix_data"]}
    elif year == "y22":
        grouping_data = {'data': ["2022_data","2022EE_data"]}
    elif year == "y23":
        grouping_data = {'data': ["2023_data","2023BPix_data"]}
    else:
        grouping_data = {'data': [f"{year}_data"]}
    return grouping_data

