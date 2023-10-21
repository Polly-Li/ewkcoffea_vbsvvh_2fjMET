import argparse
import pickle
import json
import gzip
import os
import numpy as np
import matplotlib.pyplot as plt

import hist

from topcoffea.scripts.make_html import make_html

from topcoffea.modules import utils
import topcoffea.modules.MakeLatexTable as mlt

# This script opens a pkl file of histograms produced by wwz processor
# Reads the histograms and dumps out the yields for each group of processes
# Example usage: python get_yld_check.py -f histos/tmp_histo.pkl.gz

# Colors in VVV observation
#ZZ    = (240, 155, 205)  #F09B9B
#ttZ   = (0, 208, 145) #00D091
#WZ    = (163, 155, 47) #A39B2F
#tWZ   = (205, 240, 155) #CDF09B
#Other = (205, 205, 205) #CDCDCD
CLR_LST = ["red","blue","#F09B9B","#00D091","#CDF09B","#CDCDCD"]

#7ad98d5
EWK_REF = {'WWZ': {'sr_4l_sf_A': (2.4242849899619614, 0.008680020463597452), 'sr_4l_sf_B': (2.050828817213187, 0.007941286365426964), 'sr_4l_sf_C': (0.594701078172875, 0.004324531372334359), 'all_events': (143.03255106410325, 0.06621972573333947), '4l_presel': (25.390131740285142, 0.02793500696082887), 'cr_4l_sf': (1.099154420891864, 0.005806987065182123), 'sr_4l_of_1': (0.6709899594807212, 0.004531366935275469), 'sr_4l_of_2': (0.7674646184241283, 0.004852568555205251), 'sr_4l_of_3': (1.538278290983726, 0.006872754942126762), 'sr_4l_of_4': (5.351211327137207, 0.01289551609520589), 'cr_4l_of': (0.31412254933638906, 0.003051399358542798)}, 'ZH': {'sr_4l_sf_A': (0.9588826269718993, 0.008562801130365305), 'sr_4l_sf_B': (1.5290591662496809, 0.009915184511108817), 'sr_4l_sf_C': (0.6829818391088338, 0.005461472864686166), 'all_events': (137.65612493509434, 0.09089195409138587), '4l_presel': (19.875361077707566, 0.03255127455912648), 'cr_4l_sf': (0.07069379744552862, 0.0017785467537088193), 'sr_4l_of_1': (3.0776649291565263, 0.012188953704734596), 'sr_4l_of_2': (1.3634838648740697, 0.008187925676202436), 'sr_4l_of_3': (0.35065870046946657, 0.004282999746471847), 'sr_4l_of_4': (0.15334445705866528, 0.003442573562330441), 'cr_4l_of': (0.26697188716389064, 0.004216783461442276)}, 'ZZ': {'sr_4l_sf_A': (1.3188915546124917, 0.027709709326434746), 'sr_4l_sf_B': (4.589564790327131, 0.051334985548410574), 'sr_4l_sf_C': (2.788085121097538, 0.03945147448378394), 'all_events': (20839.935963596385, 3.5284508076238685), '4l_presel': (2996.3448773944438, 1.3026525735774501), 'cr_4l_sf': (1809.4037161641972, 1.010609187977771), 'sr_4l_of_1': (0.6230402135952318, 0.01913251268840916), 'sr_4l_of_2': (0.6195523298629269, 0.019149313484256587), 'sr_4l_of_3': (0.389804485719651, 0.015214380361572485), 'sr_4l_of_4': (0.503428151627304, 0.017427271906628516), 'cr_4l_of': (0.5690950408679782, 0.018532804854137372)}, 'ttZ': {'sr_4l_sf_A': (1.4632380469702184, 0.07073097766465616), 'sr_4l_sf_B': (1.2456058174138889, 0.06511930763490584), 'sr_4l_sf_C': (0.2973263565218076, 0.03235651540657733), 'all_events': (6400.391729545197, 5.524013760526022), '4l_presel': (155.00430985842831, 0.7354360233527983), 'cr_4l_sf': (0.6504467655904591, 0.047309618956213474), 'sr_4l_of_1': (0.35825034615118057, 0.03331835991644278), 'sr_4l_of_2': (0.4247714438242838, 0.03477282822737832), 'sr_4l_of_3': (0.9259218217339367, 0.05493082335382369), 'sr_4l_of_4': (2.785680097178556, 0.09803585582198722), 'cr_4l_of': (55.61262463382445, 0.4411700479312359)}, 'tWZ': {'sr_4l_sf_A': (0.4517017470789142, 0.01738180013665498), 'sr_4l_sf_B': (0.42845417943317443, 0.016933980145656828), 'sr_4l_sf_C': (0.1038940001744777, 0.00839185864022349), 'all_events': (457.53748542949324, 0.5535594488742449), '4l_presel': (21.40333431190811, 0.11971397281075295), 'cr_4l_sf': (0.2044043393107131, 0.011675895183228412), 'sr_4l_of_1': (0.14237167895771563, 0.009755355520882436), 'sr_4l_of_2': (0.16259158292086795, 0.010390841060134557), 'sr_4l_of_3': (0.2754841863643378, 0.013584184758869055), 'sr_4l_of_4': (0.9119344529462978, 0.024725680873858506), 'cr_4l_of': (6.6433630282990634, 0.06668917741944116)}, 'other': {'sr_4l_sf_A': (0.7716581597924232, 0.2312639955738936), 'sr_4l_sf_B': (1.1314077151473612, 0.38018104654891954), 'sr_4l_sf_C': (0.36015862389467657, 0.18350490893232696), 'all_events': (2676923.753788651, 2177.829454889236), '4l_presel': (52.9465736518614, 2.5140479882052373), 'cr_4l_sf': (12.130903454846703, 0.7388418925835935), 'sr_4l_of_1': (0.7967714633559808, 0.26061527782483923), 'sr_4l_of_2': (0.1370740740094334, 0.12366394284843292), 'sr_4l_of_3': (0.17122043680865318, 0.225799520923806), 'sr_4l_of_4': (1.396175786969252, 0.2706284118286445), 'cr_4l_of': (2.322214526706375, 0.23042204650019732)}, '$S/\\sqrt{B}$': {'sr_4l_sf_A': [1.6904242563854188, None], 'sr_4l_sf_B': [1.3164349201271868, None], 'sr_4l_sf_C': [0.678174872143654, None], 'sr_4l_of_1': [2.7050531360896706, None], 'sr_4l_of_2': [1.8381249239564001, None], 'sr_4l_of_3': [1.4228575364170486, None], 'sr_4l_of_4': [2.326677270959722, None]}, '$S/\\sqrt{S+B}$': {'sr_4l_sf_A': [1.2446314148200222, None], 'sr_4l_sf_B': [1.080609413155987, None], 'sr_4l_sf_C': [0.5815376339963901, None], 'sr_4l_of_1': [1.5744136197796181, None], 'sr_4l_of_2': [1.1431400111439176, None], 'sr_4l_of_3': [0.9885295897591657, None], 'sr_4l_of_4': [1.6520610073253268, None]}, 'Sig': {'sr_4l_sf_A': [3.3831676169338607, None], 'sr_4l_sf_B': [3.579887983462868, None], 'sr_4l_sf_C': [1.2776829172817088, None], 'sr_4l_of_1': [3.7486548886372475, None], 'sr_4l_of_2': [2.130948483298198, None], 'sr_4l_of_3': [1.8889369914531926, None], 'sr_4l_of_4': [5.5045557841958725, None]}, 'Bkg': {'sr_4l_sf_A': [4.005489508454048, None], 'sr_4l_sf_B': [7.3950325023215555, None], 'sr_4l_sf_C': [3.5494641016885, None], 'sr_4l_of_1': [1.9204337020601088, None], 'sr_4l_of_2': [1.343989430617512, None], 'sr_4l_of_3': [1.7624309306265786, None], 'sr_4l_of_4': [5.59721848872141, None]}}

SOVERROOTB = "$S/\sqrt{B}$"
SOVERROOTSPLUSB = "$S/\sqrt{S+B}$"

# Keegan's yields as of Aug 18, 2023
KEEGAN_YIELDS = {
    "WWZ" : {
        "sr_4l_sf_A": 2.33,
        "sr_4l_sf_B": 1.97,
        "sr_4l_sf_C": 0.572,
        "sr_4l_of_1": 0.645,
        "sr_4l_of_2": 0.738,
        "sr_4l_of_3": 1.48,
        "sr_4l_of_4": 5.14,
    },
    "ZH" : {
        "sr_4l_sf_A": 0.952,
        "sr_4l_sf_B": 1.52,
        "sr_4l_sf_C": 0.677,
        "sr_4l_of_1": 3.05,
        "sr_4l_of_2": 1.35,
        "sr_4l_of_3": 0.348,
        "sr_4l_of_4": 0.152,
    },
    "ttZ" : {
        "sr_4l_sf_A": 1.32*(0.281/0.2529),
        "sr_4l_sf_B": 1.12*(0.281/0.2529),
        "sr_4l_sf_C": 0.269*(0.281/0.2529),
        "sr_4l_of_1": 0.322*(0.281/0.2529),
        "sr_4l_of_2": 0.382*(0.281/0.2529),
        "sr_4l_of_3": 0.833*(0.281/0.2529),
        "sr_4l_of_4": 2.51*(0.281/0.2529),
    },
    "ZZ": {
        "sr_4l_sf_A": 1.32,
        "sr_4l_sf_B": 4.58,
        "sr_4l_sf_C": 2.78,
        "sr_4l_of_1": 0.623,
        "sr_4l_of_2": 0.619,
        "sr_4l_of_3": 0.39,
        "sr_4l_of_4": 0.503,
    },
    "other": {
        "sr_4l_sf_A": 0.607 + 0.00855,
        "sr_4l_sf_B": 0.417 + 0.0571,
        "sr_4l_sf_C": 0.0967 + 0.0106,
        "sr_4l_of_1": 0.27 + 0.107 + 0.0069,
        "sr_4l_of_2": -0.055 + 0.0 + 0.0106,
        "sr_4l_of_3": -0.0555 + 0.176 + 0.0446,
        "sr_4l_of_4": 0.452 + 0.12 + 0.0233,
    },
}


sample_dict_base = {
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

    ## Indiviudally
    #"WWZJetsTo4L2Nu":            ["WWZJetsTo4L2Nu"],
    #"GluGluZH":                  ["GluGluZH"],
    #"qqToZHToZTo2L":             ["qqToZHToZTo2L"],
    #"ZZTo4l":                    ["ZZTo4l"],
    #"ggToZZTo2e2mu":             ["ggToZZTo2e2mu"],
    #"ggToZZTo2e2tau":            ["ggToZZTo2e2tau"],
    #"ggToZZTo2mu2tau":           ["ggToZZTo2mu2tau"],
    #"ggToZZTo4e":                ["ggToZZTo4e"],
    #"ggToZZTo4mu":               ["ggToZZTo4mu"],
    #"ggToZZTo4tau":              ["ggToZZTo4tau"],
    #"TTZToLL_M_1to10":           ["TTZToLL_M_1to10"],
    #"TTZToLLNuNu_M_10":          ["TTZToLLNuNu_M_10"],
    #"TTZToQQ":                   ["TTZToQQ"],
    #"tWll" :                     ["tWll"],

    ###"DYJetsToLL_M_10to50_MLM": ["DYJetsToLL_M_10to50_MLM"],
    #"DYJetsToLL_M_50_MLM":       ["DYJetsToLL_M_50_MLM"],
    #"SSWW":                      ["SSWW"],
    #"ST_antitop_t-channel":      ["ST_antitop_t-channel"],
    #"ST_top_s-channel":          ["ST_top_s-channel"],
    #"ST_top_t-channel":          ["ST_top_t-channel"],
    #"tbarW_noFullHad":           ["tbarW_noFullHad"],
    #"ttHnobb":                   ["ttHnobb"],
    #"TTTo2L2Nu":                 ["TTTo2L2Nu"],
    #"TTWJetsToLNu":              ["TTWJetsToLNu"],
    #"TTWJetsToQQ":               ["TTWJetsToQQ"],
    #"tW_noFullHad":              ["tW_noFullHad"],
    #"tZq":                       ["tZq"],
    #"VHnobb":                    ["VHnobb"],
    ###"WJetsToLNu":              ["WJetsToLNu"],
    #"WWTo2L2Nu":                 ["WWTo2L2Nu"],
    #"WZTo3LNu":                  ["WZTo3LNu"],

    #"WWW" : ["WWW"],
    #"WZZ" : ["WZZ"],
    #"ZZZ" : ["ZZZ"],


}



# Pass dictionary with the base names for the samples, and return with full list for 4 years
def create_full_sample_dict(in_dict):
    out_dict = {}
    years = ["UL16APV","UL16","UL17","UL18"]
    for proc_group in in_dict.keys():
        out_dict[proc_group] = []
        for proc_base_name in in_dict[proc_group]:
            for year_str in years:
                out_dict[proc_group].append(f"{year_str}_{proc_base_name}")
    return out_dict


sample_dict = create_full_sample_dict(sample_dict_base)




################### Getting and printing yields ###################

# Get the yields in the SR
def get_yields(histos_dict,raw_counts=False,quiet=False):

    yld_dict = {}

    # Look at the yields in one histo (e.g. njets)
    if raw_counts: dense_axis = "njets_counts"
    else: dense_axis = "njets"
    for proc_name in sample_dict.keys():
        yld_dict[proc_name] = {}
        for cat_name in histos_dict[dense_axis].axes["category"]:
            val = sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name]}].values(flow=True)))
            var = np.sqrt(sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name]}].variances(flow=True))))
            yld_dict[proc_name][cat_name] = (val,var)

    # Print to screen
    if not quiet:
        for proc in yld_dict.keys():
            print(f"\n{proc}:")
            for cat in yld_dict[proc].keys():
                val = yld_dict[proc][cat]
                print(f"\t{cat}: {val}")

    return yld_dict


# Gets the S/sqrt(B) and puts it into the dict
# Hard coded for the summed values (e.g. looking for "ZH" not "GluGluZH","qqToZHToZTo2L")
def put_s_over_root_b(yld_dict):
    sig_lst = ["WWZ","ZH"]
    bkg_lst = ["ZZ","ttZ","tWZ","other"]
    sig_sum = {"sr_4l_sf_A":0, "sr_4l_sf_B":0, "sr_4l_sf_C":0, "sr_4l_of_1":0, "sr_4l_of_2":0, "sr_4l_of_3":0, "sr_4l_of_4":0}
    bkg_sum = {"sr_4l_sf_A":0, "sr_4l_sf_B":0, "sr_4l_sf_C":0, "sr_4l_of_1":0, "sr_4l_of_2":0, "sr_4l_of_3":0, "sr_4l_of_4":0}
    for proc in yld_dict.keys():
        print(proc)
        for cat in yld_dict[proc].keys():
            if cat not in sig_sum: continue
            val,err = yld_dict[proc][cat]
            print("   ",cat,val)
            if proc in sig_lst:
                sig_sum[cat] += val
            if proc in bkg_lst:
                bkg_sum[cat] += val

    yld_dict[SOVERROOTB] = {}
    yld_dict[SOVERROOTSPLUSB] = {}
    yld_dict["Sig"] = {}
    yld_dict["Bkg"] = {}
    for cat in sig_sum.keys():
        yld_dict[SOVERROOTB][cat]      = [sig_sum[cat]/np.sqrt(bkg_sum[cat]) , None]
        yld_dict[SOVERROOTSPLUSB][cat] = [sig_sum[cat]/np.sqrt(sig_sum[cat] + bkg_sum[cat]) , None]
        yld_dict["Sig"][cat] = [sig_sum[cat], None]
        yld_dict["Bkg"][cat] = [bkg_sum[cat], None]


# Print yields
def print_yields(yld_dict,print_fom=True):

    # Dump the yields to dict for latex table
    yld_dict_for_printing = {}
    cats_to_print = [ "sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]
    for proc in yld_dict.keys():
        yld_dict_for_printing[proc] = {}
        for cat in yld_dict[proc].keys():
            if cat not in cats_to_print: continue
            yld_dict_for_printing[proc][cat] = yld_dict[proc][cat]


    # Print the yields directly
    mlt.print_latex_yield_table(
        yld_dict_for_printing,
        tag="All yields",
        key_order=sample_dict_base.keys(),
        subkey_order=cats_to_print,
        print_begin_info=True,
        print_end_info=True,
        print_errs=True,
        column_variable="subkeys",
    )
    #exit()


    # Compare with other yields, print comparison

    #tag1 = "ewkcoffea"
    #tag2 = "VVVNanoLooper"
    tag1 = "btagDeepFlav 0m"
    tag2 = "btagCSV 0l"

    #yld_dict_comp = utils.put_none_errs(KEEGAN_YIELDS)
    yld_dict_comp = EWK_REF

    yld_dict_1 = yld_dict_for_printing
    yld_dict_2 = yld_dict_comp

    pdiff_dict = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="percent_diff",inpercent=True)
    diff_dict  = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="absolute_diff")

    #utils.print_yld_dicts(yld_dict_1,tag1)
    #utils.print_yld_dicts(yld_dict_2,tag2)
    #utils.print_yld_dicts(pdiff_dict,f"Percent diff between {tag1} and {tag2}")
    #utils.print_yld_dicts(diff_dict,f"Diff between {tag1} and {tag2}")

    procs_to_print = list(sample_dict_base.keys())
    if print_fom:
        procs_to_print.append("Sig")
        procs_to_print.append("Bkg")
        procs_to_print.append(SOVERROOTB)
        procs_to_print.append(SOVERROOTSPLUSB)

    mlt.print_begin()
    mlt.print_latex_yield_table(yld_dict_1,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag1)
    mlt.print_latex_yield_table(yld_dict_2,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag2)
    mlt.print_latex_yield_table(pdiff_dict,key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Percent diff between {tag1} and {tag2}")
    mlt.print_latex_yield_table(diff_dict, key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Diff between {tag1} and {tag2}")
    mlt.print_end()


# Dump the counts dict to a latex table
def print_counts(counts_dict):

    cats_to_print = ["all_events", "4l_presel", "sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]

    # Print the yields directly
    mlt.print_latex_yield_table(
        counts_dict,
        tag="Raw MC counts (ewkcoffea)",
        key_order=counts_dict.keys(),
        subkey_order=cats_to_print,
        print_begin_info=True,
        print_end_info=True,
        column_variable="subkeys",
    )


# This should maybe be in a different script
################### Hist manipulation and plotting ###################


# Get the list of categories on the sparese axis
def get_axis_cats(histo,axis_name):
    process_list = [x for x in histo.axes[axis_name]]
    return process_list


# Regroup categories (e.g. processes)
def group(h, oldname, newname, grouping):

    # Build up a grouping dict that drops any proc that is not in our h
    grouping_slim = {}
    proc_lst = get_axis_cats(h,"process")
    for grouping_name in grouping.keys():
        for proc in grouping[grouping_name]:
            if proc in proc_lst:
                if grouping_name not in grouping_slim:
                    grouping_slim[grouping_name] = []
                grouping_slim[grouping_name].append(proc)
            #else:
            #    print(f"WARNING: process {proc} not in this hist")

    # From Nick: https://github.com/CoffeaTeam/coffea/discussions/705#discussioncomment-4604211
    hnew = hist.Hist(
        hist.axis.StrCategory(grouping_slim, name=newname),
        *(ax for ax in h.axes if ax.name != oldname),
        storage=h._storage_type,
    )
    for i, indices in enumerate(grouping_slim.values()):
        hnew.view(flow=True)[i] = h[{oldname: indices}][{oldname: sum}].view(flow=True)

    return hnew


# Takes a mc hist and data hist and plots both
def make_cr_fig(histo_mc,histo_data,title,unit_norm_bool=False):

    # Create the figure
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(7,7),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.07)

    # Plot the mc
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        color=CLR_LST,
        ax=ax,
        flow=None,
    )
    # Plot the data
    histo_data.plot1d(
        stack=False,
        histtype="errorbar",
        color="k",
        ax=ax,
        w2=histo_data.variances(),
        w2method="sqrt",
        flow=None,
    )

    ### Get the err and ratios and plot them by hand ###
    histo_mc_sum = histo_mc[{"process_grp":sum}]
    histo_data_sum = histo_data[{"process_grp":sum}]

    mc_arr = histo_mc_sum.values()
    mc_err_arr = np.sqrt(histo_mc_sum.variances())
    data_arr = histo_data_sum.values()
    data_err_arr = np.sqrt(histo_data_sum.variances())

    err_p = np.append(mc_arr + mc_err_arr, 0)
    err_m = np.append(mc_arr - mc_err_arr, 0)
    err_ratio_p = np.append(1+mc_err_arr/mc_arr,1)
    err_ratio_m = np.append(1-mc_err_arr/mc_arr,1)

    data_ratio_err_p = (data_arr + data_err_arr)/mc_arr
    data_ratio_err_m = (data_arr - data_err_arr)/mc_arr

    bin_edges_arr = histo_mc_sum.axes[0].edges
    bin_centers_arr = histo_mc_sum.axes[0].centers

    ax.fill_between(bin_edges_arr,err_m,err_p, step='post', facecolor='none', edgecolor='gray', alpha=0.5, linewidth=0.0, label='MC stat', hatch='/////')
    rax.fill_between(bin_edges_arr,err_ratio_m,err_ratio_p,step='post', facecolor='none',edgecolor='gray', label='MC stat', linewidth=0.0, hatch='/////',alpha=0.5)
    rax.scatter(bin_centers_arr,data_arr/mc_arr,facecolor='black',edgecolor='black',marker="o")
    rax.vlines(bin_centers_arr,data_ratio_err_p,data_ratio_err_m,color='k')

    # Scale the y axis and labels
    ax.legend()
    ax.set_title(title)
    ax.autoscale(axis='y')
    ax.set_xlabel(None)
    rax.set_ylabel('Ratio')
    rax.set_ylim(0.0,2.0)
    rax.axhline(1.0,linestyle="-",color="k",linewidth=1)

    return fig


# Plots a hist
def make_single_fig(histo_mc,title,unit_norm_bool=False):
    #print("\nPlotting values:",histo.values())
    fig, ax = plt.subplots(1, 1, figsize=(7,7))

    # Guard against accidently looking at data before unblind
    proc_lst = get_axis_cats(histo_mc,"process_grp")
    for proc_name in proc_lst:
        if "data" in proc_name:
            raise Exception(f"CAUTION: Are you sure you want to look at this data in {proc_name}?")

    # Plot the mc
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        #color=["red","yellow","blue","green","orange","grey"],
        color=CLR_LST,
        yerr=True,
    )

    plt.title(title)
    ax.autoscale(axis='y')
    plt.legend()
    return fig



# Main function for making CR plots
def make_plots(histo_dict):

    save_dir_path = "plots"

    for var_name in histo_dict.keys():
        print(f"\n{var_name}")
        if "counts" in var_name: continue
        #if var_name != "njets": continue # TMP
        histo = histo_dict[var_name]

        # Rebin if continous variable
        if var_name not in ["njets","nbtagsl","nleps"]:
            histo = histo[..., ::hist.rebin(4)] # Rebin according to https://github.com/CoffeaTeam/coffea/discussions/705

        # Loop over categories and make plots for each
        for cat_name in histo.axes["category"]:
            if cat_name not in ["cr_4l_sf","cr_4l_of"]: continue # TMP
            print(cat_name)

            histo_cat = histo[{"category":cat_name}]
            #histo_cat.plot1d(overlay="process_grp")

            # Group the mc samples
            grouping_mc = sample_dict
            histo_grouped_mc = group(histo_cat,"process","process_grp",grouping_mc)

            # Group the data samples
            grouping_data = {'data': ["UL16APV_data","UL16_data","UL17_data","UL18_data"]}
            histo_grouped_data = group(histo_cat,"process","process_grp",grouping_data)

            #####
            ##if cat_name == "cr_4l_sf" and var_name == "nleps":
            #if cat_name == "cr_4l_of" and var_name == "nleps":
            #    print("mc\n",histo_grouped_mc)
            #    print("data\n",histo_grouped_data)
            #    print("val mc\n",histo_grouped_mc.values(flow=True))
            #    print("val data\n",histo_grouped_data.values(flow=True))
            #    print("var mc\n",np.sqrt(histo_grouped_mc.variances()))
            #    print("var data\n",np.sqrt(histo_grouped_data.variances()))
            #continue
            #print("\Yields")
            #print(type(histo_cat.values(flow=True)))
            #print(sum(histo_cat.values(flow=True)))
            #print(sum(sum(histo_cat.values(flow=True))))
            #exit()
            ####

            # Make figure
            title = f"{cat_name}_{var_name}"
            if "cr" in title:
                fig = make_cr_fig(histo_grouped_mc,histo_grouped_data,title=title)
            else:
                fig = make_single_fig(histo_grouped_mc,title=title)

            # Save
            save_dir_path_cat = os.path.join(save_dir_path,cat_name)
            if not os.path.exists(save_dir_path_cat): os.mkdir(save_dir_path_cat)
            fig.savefig(os.path.join(save_dir_path_cat,title+".pdf"))
            fig.savefig(os.path.join(save_dir_path_cat,title+".png"))

            make_html(os.path.join(os.getcwd(),save_dir_path_cat))



################### Main ###################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default=".", help = "The path the output files should be saved to")
    parser.add_argument("-n", "--output-name", default="wwz_yields", help = "A name for the output directory")
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))

    # Wrapper around the code for getting the raw counts and dump to latex table
    #counts_dict = get_yields(histo_dict,raw_counts=True)
    #print_counts(counts_dict)
    #exit()

    # Wrapper around the code for getting the yields for sr and bkg samples
    yld_dict = get_yields(histo_dict)
    put_s_over_root_b(yld_dict)
    #print(yld_dict)
    #exit()
    print_yields(yld_dict)

    # Test plotting
    #make_plots(histo_dict)
    #exit()

    # Dump yield dict to json
    if "json" not in args.output_name: output_name = args.output_name + ".json"
    else: output_name = args.output_name
    with open(output_name,"w") as out_file: json.dump(yld_dict, out_file, indent=4)
    print(f"\nSaved json file: {output_name}\n")



if __name__ == "__main__":
    main()

