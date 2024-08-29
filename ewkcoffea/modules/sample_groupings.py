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

all_exclude_dict = {
    "TTZToLL_M_1to10": ["2022","2022EE"],
    "TTZToLLNuNu_M_10": ["2022","2022EE"],
    "SSWW": ["2022","2022EE"],
    "tbarW_noFullHad": ["2022","2022EE"],
    "TTWJetsToLNu": ["2022","2022EE"],
    "TTWJetsToQQ": ["2022","2022EE"],
    "tW_noFullHad": ["2022","2022EE"],
    "GluGluZH": ["2022","2022EE"],
    "qqToZHToZTo2L": ["2022","2022EE"],

    "TTZToLL_M_4to50": ["UL16APV","UL16","UL17","UL18"],
    "TTZToLL_M_50": ["UL16APV","UL16","UL17","UL18"],
    "GluGluZHTo2WTo2L2Nu": ["UL16APV","UL16","UL17","UL18"],
    "qqToZHTo2WTo2L2Nu": ["UL16APV","UL16","UL17","UL18"],
    "DYJetsToLL_M_10to50_MLM": ["UL16APV","UL16","UL17","UL18"],
    "SSWW_TT": ["UL16APV","UL16","UL17","UL18"],
    "SSWW_TL": ["UL16APV","UL16","UL17","UL18"],
    "SSWW_LL": ["UL16APV","UL16","UL17","UL18"],
    "tbarW_leptonic": ["UL16APV","UL16","UL17","UL18"],
    "tbarW_semileptonic": ["UL16APV","UL16","UL17","UL18"],
    "tW_leptonic": ["UL16APV","UL16","UL17","UL18"],
    "tW_semileptonic": ["UL16APV","UL16","UL17","UL18"],
    "WJetsToLNu": ["UL16APV","UL16","UL17","UL18"],
}


# The "official" groupings
SAMPLE_DICT_BASE_ALL = {
    "WWZ" : ["WWZJetsTo4L2Nu"],
    "ZH"  : ["GluGluZH","qqToZHToZTo2L","GluGluZHTo2WTo2L2Nu","qqToZHTo2WTo2L2Nu"],
    "ZZ"  : ["ZZTo4l", "ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau"],

    "ttZ" : [
        "TTZToLL_M_1to10",
        "TTZToLLNuNu_M_10",
        "TTZToQQ",
        "TTZToLL_M_4to50",
        "TTZToLL_M_50",
    ],
    "tWZ" : ["tWZ4l"], 
    "WZ" : ["WZTo3LNu"],

    "other" : [
        "DYJetsToLL_M_50_MLM",
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
        "WWTo2L2Nu",
        "WWW",
        "WZZ",
        "ZZZ",
        "ggHToZZ4L",
        "DYJetsToLL_M_10to50_MLM",
        "SSWW_TT",
        "SSWW_TL",
        "SSWW_LL",
        "tbarW_leptonic",
        "tbarW_semileptonic",
        "tW_leptonic",
        "tW_semileptonic",
        "WJetsToLNu",
    ],
}


# The "official" groupings
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

    "WZ" : ["WZTo3LNu"],

    "other" : [

        ##"WWZJetsTo4L2Nu",
        ##"GluGluZH","qqToZHToZTo2L",
        ##"ZZTo4l", "ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau",
        ##"TTZToLL_M_1to10","TTZToLLNuNu_M_10","TTZToQQ",
        ##"tWZ4l",

        ##"DYJetsToLL_M_10to50_MLM",
        "DYJetsToLL_M_50_MLM",
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
    "WZ" : ["WZTo3LNu"],
    "other" : [
        "DYJetsToLL_M_10to50_MLM",
        "DYJetsToLL_M_50_MLM",
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
        #"TTWJetsToLNu",
        #"TTWJetsToQQ",
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
    "tWZ4l" :                     ["tWZ4l"], # Old sample tWll

    ##"DYJetsToLL_M_10to50_MLM": ["DYJetsToLL_M_10to50_MLM"],
    "DYJetsToLL_M_50_MLM":       ["DYJetsToLL_M_50_MLM"],
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
    "DYJetsToLL_M_10to50_MLM" : ["DYJetsToLL_M_10to50_MLM"],
    "DYJetsToLL_M_50_MLM"     : ["DYJetsToLL_M_50_MLM"],
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
    #"TTWJetsToLNu"            : ["TTWJetsToLNu"],
    #"TTWJetsToQQ"             : ["TTWJetsToQQ"],
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

# Pass dictionary with the base names for the samples, and return with full list for 4 years
def create_mc_sample_dict(year,yld_individual=False):
    out_dict = {}
    all_years = ["UL16APV","UL16","UL17","UL18","2022","2022EE"]
    r2_years = ["UL16APV","UL16","UL17","UL18"]
#    r3_years = ["2022","2022EE","2023","2023BPix"]
    r3_years = ["2022","2022EE"]
    y22_years = ["2022","2022EE"]
    y23_years = ["2023","2023BPix"]
    if year == "all":
        years = all_years
        sample_dict_base = SAMPLE_DICT_BASE_ALL
    elif year == "run2":
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
                if proc_base_name in all_exclude_dict.keys():
                    if year_str in all_exclude_dict[proc_base_name]:
                        continue
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

