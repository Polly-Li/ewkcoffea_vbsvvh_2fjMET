# Dict to keep track of sample groupings, and associated tools


######################## Lists ########################

PROC_LST = ["DY"]
BKG_LST = ["DY"]
SIG_LST = []

CAT_LST_CB = ["mumu", "ee"]

######################## Dictionaries ########################

# The "official" groupings
SAMPLE_DICT_BASE = {
    "DY" : ["DYJetsToLL_M_50_MLM"],
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
        sample_dict_base = SAMPLE_DICT_BASE
    elif year == "run3":
        years = r3_years
        sample_dict_base = SAMPLE_DICT_BASE
    else:
        years = [year]
        if year in r2_years:
            sample_dict_base = SAMPLE_DICT_BASE
        elif year in r3_years:
            sample_dict_base = SAMPLE_DICT_BASE
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

