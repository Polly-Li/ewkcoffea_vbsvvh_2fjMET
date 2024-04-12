# Dict to keep track of sample groupings, and associated tools


######################## Lists ########################

PROC_LST = ["WWZ","ZH","ZZ","ttZ","tWZ","WZ","other"]
SIG_LST = ["WWZ","ZH"]
BKG_LST = ["ZZ","ttZ","tWZ","WZ","other"]

CAT_LST_CB = ["sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]

SR_SF_BDT = ["sr_4l_bdt_sf_1", "sr_4l_bdt_sf_2", "sr_4l_bdt_sf_3", "sr_4l_bdt_sf_4", "sr_4l_bdt_sf_5", "sr_4l_bdt_sf_6", "sr_4l_bdt_sf_7"]
SR_OF_BDT = ["sr_4l_bdt_of_1", "sr_4l_bdt_of_2", "sr_4l_bdt_of_3", "sr_4l_bdt_of_4", "sr_4l_bdt_of_5", "sr_4l_bdt_of_6", "sr_4l_bdt_of_7", "sr_4l_bdt_of_8"]
CAT_LST_BDT = SR_SF_BDT + SR_OF_BDT


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

        "sr_4l_bdt_of_1" : "cr_4l_sf",
        "sr_4l_bdt_of_2" : "cr_4l_sf",
        "sr_4l_bdt_of_3" : "cr_4l_sf",
        "sr_4l_bdt_of_4" : "cr_4l_sf",
        "sr_4l_bdt_of_5" : "cr_4l_sf",
        "sr_4l_bdt_of_6" : "cr_4l_sf",
        "sr_4l_bdt_of_7" : "cr_4l_sf",
        "sr_4l_bdt_of_8" : "cr_4l_sf",

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

        "sr_4l_bdt_of_1" : "cr_4l_btag_of",
        "sr_4l_bdt_of_2" : "cr_4l_btag_of",
        "sr_4l_bdt_of_3" : "cr_4l_btag_of",
        "sr_4l_bdt_of_4" : "cr_4l_btag_of",
        "sr_4l_bdt_of_5" : "cr_4l_btag_of",
        "sr_4l_bdt_of_6" : "cr_4l_btag_of",
        "sr_4l_bdt_of_7" : "cr_4l_btag_of",
        "sr_4l_bdt_of_8" : "cr_4l_btag_of",
    }
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
        #"TTZToQQ",
    ],

    "tWZ" : ["tWZ4l"],

    "WZ" : ["WZTo3LNu"],

    "other" : [

        "DYJetsToLL_M_10to50_MLM",
        "DYJetsToLL_M_50_MLM",
        #"SSWW",
        #"ST_antitop_t-channel",
        #"ST_top_s-channel",
        #"ST_top_t-channel",
        "tbarW_noFullHad",
        "ttHnobb",
        #"TTTo2L2Nu",
        #"TTWJetsToLNu",
        #"TTWJetsToQQ",
        "tW_noFullHad",
        #"tZq",
        "VHnobb",
        "WJetsToLNu",
        "WWTo2L2Nu",
        "WZTo3LNu",

        "WWW",
        "WZZ",
        #"ZZZ",
    ],
}

# Processes indiviudally
SAMPLE_DICT_BASE_INDIV = {
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
    #"TTZToQQ"                 : ["TTZToQQ"],
    "tWZ4l"                   : ["tWZ4l"],
    "WZTo3LNu"                : ["WZTo3LNu"],

    "DYJetsToLL_M_10to50_MLM" : ["DYJetsToLL_M_10to50_MLM"],
    "DYJetsToLL_M_50_MLM"     : ["DYJetsToLL_M_50_MLM"],
    #"SSWW"                    : ["SSWW"],
    #"ST_antitop_t-channel"    : ["ST_antitop_t-channel"],
    #"ST_top_s-channel"        : ["ST_top_s-channel"],
    #"ST_top_t-channel"        : ["ST_top_t-channel"],
    "tbarW_noFullHad"         : ["tbarW_noFullHad"],
    "ttHnobb"                 : ["ttHnobb"],
    #"TTTo2L2Nu"               : ["TTTo2L2Nu"],
    #"TTWJetsToLNu"            : ["TTWJetsToLNu"],
    #"TTWJetsToQQ"             : ["TTWJetsToQQ"],
    "tW_noFullHad"            : ["tW_noFullHad"],
    #"tZq"                     : ["tZq"],
    "VHnobb"                  : ["VHnobb"],
    "WJetsToLNu"              : ["WJetsToLNu"],
    "WWTo2L2Nu"               : ["WWTo2L2Nu"],
    "WWW"                     : ["WWW"],
    "WZZ"                     : ["WZZ"],
    #"ZZZ"                     : ["ZZZ"],
}


######################## Tools ########################

# Pass dictionary with the base names for the samples, and return with full list for 4 years
def create_mc_sample_dict(year):
    out_dict = {}
    r2_years = ["UL16APV","UL16","UL17","UL18"]
    r3_years = ["2022","2022EE"]
    if year == "all":
        raise Exception("ERROR: We are not ready to sum Run2 and Run3.")
    elif year == "run2":
        years = r2_years
        sample_dict_base = SAMPLE_DICT_BASE_RUN2
    elif year == "run3":
        years = r3_years
        sample_dict_base = SAMPLE_DICT_BASE_RUN3
    else:
        years = [year]
        if year in r2_years:
            sample_dict_base = SAMPLE_DICT_BASE_RUN2
        elif yeaer in r3_years:
            sample_dict_base = SAMPLE_DICT_BASE_RUN3
        else:
            raise Exception(f"ERROR: Unrecognized year \"{year}\". Exiting.")

    for proc_group in sample_dict_base.keys():
        out_dict[proc_group] = []
        for proc_base_name in sample_dict_base[proc_group]:
            for year_str in years:
                out_dict[proc_group].append(f"{year_str}_{proc_base_name}")

    return out_dict

# Get data sampel dict
def create_data_sample_dict(year):
    if year == "all":
        raise Exception("ERROR: We are not ready to run over Run2 and Run3.")
        #grouping_data = {'data': ["UL16APV_data","UL16_data","UL17_data","UL18_data","2022_data","2022EE_data"]}
    elif year == "run2":
        grouping_data = {'data': ["UL16APV_data","UL16_data","UL17_data","UL18_data"]}
    elif year == "run3":
        grouping_data = {'data': ["2022_data","2022EE_data"]}
    else:
        grouping_data = {'data': [f"{year}_data"]}
    return grouping_data

