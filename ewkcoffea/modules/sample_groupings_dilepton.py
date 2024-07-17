# Dict to keep track of sample groupings, and associated tools


######################## Lists ########################

PROC_LST = ["DY","W","WW","TT"]
SIG_LST = ["DY"]
BKG_LST = ["W","WW","TT"]

#CAT_LST_CB = ["sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]

#SR_SF_BDT = ["sr_4l_bdt_sf_1", "sr_4l_bdt_sf_2", "sr_4l_bdt_sf_3", "sr_4l_bdt_sf_4", "sr_4l_bdt_sf_5", "sr_4l_bdt_sf_6", "sr_4l_bdt_sf_7"]
#SR_OF_BDT = ["sr_4l_bdt_of_1", "sr_4l_bdt_of_2", "sr_4l_bdt_of_3", "sr_4l_bdt_of_4", "sr_4l_bdt_of_5", "sr_4l_bdt_of_6", "sr_4l_bdt_of_7", "sr_4l_bdt_of_8"]
#CAT_LST_BDT = SR_SF_BDT + SR_OF_BDT


######################## Dictionaries ########################



# The "official" groupings
SAMPLE_DICT_BASE = {
    "DY" : ["DYJetsToLL_M_50_MLM","DYJetsToLL_M_10to50_MLM"],
    "W" : ["WJetsToLNu"],
    "WW" : ["WWTo2L2Nu"],
    "TT" : ["TTTo2L2Nu"],
}


# Processes indiviudally
SAMPLE_DICT_BASE_INDIV = {
    "DYJetsToLL_M_50_MLM"      : ["DYJetsToLL_M_50_MLM"],
    "DYJetsToLL_M_10to50_MLM"  : ["DYJetsToLL_M_10to50_MLM"],
    "WJetsToLNu"               : ["WJetsToLNu"],
    "WWTo2L2Nu"                : ["WWTo2L2Nu"],
    "TTTo2L2Nu"                : ["TTTo2L2Nu"],
}

######################## Tools ########################

# Pass dictionary with the base names for the samples, and return with full list for 4 years
def create_mc_sample_dict(year):
    out_dict = {}
    r2_years = ["UL16APV","UL16","UL17","UL18"]
    r3_years = ["2022","2022EE","2023","2023BPix"]
    if year == "all":
        raise Exception("ERROR: We are not ready to sum Run2 and Run3.")
    elif year == "run2":
        years = r2_years
        sample_dict_base = SAMPLE_DICT_BASE
    elif year == "run3":
        years = r3_years
        sample_dict_base = SAMPLE_DICT_BASE
    elif year == "yr23":
        years = ["2023","2023BPix"]
        sample_dict_base = SAMPLE_DICT_BASE
    elif year == "yr22":
        years = ["2022","2022EE"]
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
    elif year == "run2":
        grouping_data = {'data': ["UL16APV_data","UL16_data","UL17_data","UL18_data"]}
    elif year == "run3":
        grouping_data = {'data': ["2022_data","2022EE_data","2023_data","2023BPix_data"]}
    elif year == "yr22":
        grouping_data = {'data': ["2022_data","2022EE_data"]}
    elif year == "yr23":
        grouping_data = {'data': ["2023_data","2023BPix_data"]}
    else:
        grouping_data = {'data': [f"{year}_data"]}
    return grouping_data

