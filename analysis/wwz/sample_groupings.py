# Dict to keep track of sample groupings, and associated tools

######################## Dictionaries ########################

# The "official" groupings
SAMPLE_DICT_BASE = {
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

    "tWZ" : ["tWll"],
    #"tWZ" : ["tWZ4l"],

    "other" : [

        ##"WWZJetsTo4L2Nu",
        ##"GluGluZH","qqToZHToZTo2L",
        ##"ZZTo4l", "ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau",
        ##"TTZToLL_M_1to10","TTZToLLNuNu_M_10","TTZToQQ",
        ##"tWll",

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
        "WZTo3LNu",

        "WWW",
        "WZZ",
        "ZZZ",
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
    "tWll" :                     ["tWll"],
    #"tWZ4l" :                     ["tWZ4l"],

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


######################## Tools ########################

# Pass dictionary with the base names for the samples, and return with full list for 4 years
def create_mc_sample_dict(in_dict,year):
    out_dict = {}
    if year == "all":
        years = ["UL16APV","UL16","UL17","UL18"]
    else:
        years = [year]
    for proc_group in in_dict.keys():
        out_dict[proc_group] = []
        for proc_base_name in in_dict[proc_group]:
            for year_str in years:
                out_dict[proc_group].append(f"{year_str}_{proc_base_name}")
                #out_dict[proc_group].append(f"{proc_base_name}{year_str}") # TOP22006 format
    return out_dict

# Get data sampel dict
def create_data_sample_dict(year):
    if year == "all":
        grouping_data = {'data': ["UL16APV_data","UL16_data","UL17_data","UL18_data"]}
        #grouping_data = {'data': ["dataUL16APV","dataUL16","dataUL17","dataUL18"]} # TOP22006 format
    else:
        grouping_data = {'data': [f"{year}_data"]}
        #grouping_data = {'data': [f"data{year}"]} # TOP22006 format
    return grouping_data

