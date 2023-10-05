import argparse
import pickle
import json
import gzip
import os
import numpy as np
import matplotlib.pyplot as plt

from topcoffea.modules import utils
import topcoffea.modules.MakeLatexTable as mlt

# This script opens a pkl file of histograms produced by wwz processor
# Reads the histograms and dumps out the yields for each group of processes
# Example usage: python get_yld_check.py -f histos/tmp_histo.pkl.gz

EWK_PFMET = {'WWZ': {'4l_wwz_sf_A': (2.4805652549148363, 0.008789550103292692), '4l_wwz_sf_B': (1.8165616188380227, 0.00747052524914959), '4l_wwz_sf_C': (0.5534816082181351, 0.004188046949095953), '4l_wwz_of_1': (0.6708585904361826, 0.004534885362190564), '4l_wwz_of_2': (0.7723650003790681, 0.004866945190715124), '4l_wwz_of_3': (1.5525118519453827, 0.006911280413392028), '4l_wwz_of_4': (5.351180291743731, 0.01289547612370954)}, 'ZH': {'4l_wwz_sf_A': (1.0238004361544881, 0.008788955035482712), '4l_wwz_sf_B': (1.3359813175648014, 0.00934077047862732), '4l_wwz_sf_C': (0.6501900931434648, 0.005294535535514352), '4l_wwz_of_1': (3.083909928617686, 0.012228834559895143), '4l_wwz_of_2': (1.3766249652217084, 0.008169022817831762), '4l_wwz_of_3': (0.3649182918288716, 0.004309977848419291), '4l_wwz_of_4': (0.15344176035887358, 0.0034426940282646663)}, 'ZZ': {'4l_wwz_sf_A': (1.4521438367992232, 0.029205042479336828), '4l_wwz_sf_B': (4.677904277159541, 0.05234694277251101), '4l_wwz_sf_C': (3.389974884881667, 0.04426395340496686), '4l_wwz_of_1': (0.7480083064656355, 0.021154910976938047), '4l_wwz_of_2': (0.9564218484774756, 0.023975422645129762), '4l_wwz_of_3': (0.623115328970016, 0.0195838118311669), '4l_wwz_of_4': (0.503428151627304, 0.017427271906628516)}, 'ttZ': {'4l_wwz_sf_A': (1.5473141528200358, 0.0715998885477053), '4l_wwz_sf_B': (1.0910299081588164, 0.06124966375859094), '4l_wwz_sf_C': (0.24564968887716532, 0.030378417285157187), '4l_wwz_of_1': (0.3568343404913321, 0.0331584923690703), '4l_wwz_of_2': (0.4411137462593615, 0.03545671638186383), '4l_wwz_of_3': (0.9048907430842519, 0.05501831019250101), '4l_wwz_of_4': (2.784085821127519, 0.09804881816154837)}, 'other': {'4l_wwz_sf_A': (1.1430392182664946, 0.4580699456518573), '4l_wwz_sf_B': (2.022533335839398, 0.5401745037312242), '4l_wwz_sf_C': (0.585361891426146, 0.2581708633603538), '4l_wwz_of_1': (3.3825075863860548, 0.7125001691912998), '4l_wwz_of_2': (1.7277471584966406, 0.47367528844190904), '4l_wwz_of_3': (0.12254117743577808, 0.26617433325376494), '4l_wwz_of_4': (1.343750176136382, 0.29745430223262964)}}
EWK_PUPPIMET = {'WWZ': {'4l_wwz_sf_A': (2.424348248217939, 0.008680002565746428), '4l_wwz_sf_C': (0.4677769291865843, 0.003835822331958304), '4l_wwz_of_1': (0.6708247330916493, 0.004530849496833588), '4l_wwz_of_2': (0.7673509698361158, 0.004852220739912817), '4l_wwz_sf_B': (1.8689647813753254, 0.007584448169001832), '4l_wwz_of_3': (1.5385571659608104, 0.006873341626007332), '4l_wwz_of_4': (5.351180291743731, 0.01289547612370954)}, 'ZH': {'4l_wwz_sf_A': (0.9588826269718993, 0.008562801130365305), '4l_wwz_sf_C': (0.5165593981819256, 0.0048819432342064), '4l_wwz_of_1': (3.0776154328386838, 0.012188943059530692), '4l_wwz_of_2': (1.3634339364716652, 0.008187905808049156), '4l_wwz_sf_B': (1.3810158226897329, 0.00942694896955524), '4l_wwz_of_3': (0.35068535435402737, 0.004283041214261076), '4l_wwz_of_4': (0.15344176035887358, 0.0034426940282646663)}, 'ZZ': {'4l_wwz_sf_A': (1.3188915546124917, 0.027709709326434746), '4l_wwz_sf_C': (1.8970280267531052, 0.03264529465008086), '4l_wwz_of_1': (0.6225106355814205, 0.019125182061863494), '4l_wwz_of_2': (0.6200819078767381, 0.01915663487648714), '4l_wwz_sf_B': (3.5339360351426876, 0.04509486030995772), '4l_wwz_of_3': (0.389804485719651, 0.015214380361572485), '4l_wwz_of_4': (0.503428151627304, 0.017427271906628516)}, 'ttZ': {'4l_wwz_sf_A': (1.46149617806077, 0.07070952619060875), '4l_wwz_sf_C': (0.24784770980477333, 0.030016927771481437), '4l_wwz_of_1': (0.35825034615118057, 0.03331835991644278), '4l_wwz_of_2': (0.4247714438242838, 0.03477282822737832), '4l_wwz_sf_B': (1.142400507000275, 0.06225512607803385), '4l_wwz_of_3': (0.9275160977849737, 0.05490768287045149), '4l_wwz_of_4': (2.784085821127519, 0.09804881816154837)}, 'other': {'4l_wwz_sf_A': (0.6555535732768476, 0.4433544142855943), '4l_wwz_sf_C': (0.4607994193211198, 0.23877620130436167), '4l_wwz_of_1': (3.6120765920495614, 0.7101441628058969), '4l_wwz_of_2': (1.3026276694145054, 0.4530368022583915), '4l_wwz_sf_B': (2.212222555070184, 0.5652786715257617), '4l_wwz_of_3': (0.21378649363759905, 0.2694498342344747), '4l_wwz_of_4': (1.343750176136382, 0.29745430223262964)}}

# Keegan's yields as of Aug 18, 2023
KEEGAN_YIELDS = {
    "WWZ" : {
        "4l_wwz_sf_A": 2.33,
        "4l_wwz_sf_B": 1.97,
        "4l_wwz_sf_C": 0.572,
        "4l_wwz_of_1": 0.645,
        "4l_wwz_of_2": 0.738,
        "4l_wwz_of_3": 1.48,
        "4l_wwz_of_4": 5.14,
    },
    "ZH" : {
        "4l_wwz_sf_A": 0.952,
        "4l_wwz_sf_B": 1.52,
        "4l_wwz_sf_C": 0.677,
        "4l_wwz_of_1": 3.05,
        "4l_wwz_of_2": 1.35,
        "4l_wwz_of_3": 0.348,
        "4l_wwz_of_4": 0.152,
    },
    "ttZ" : {
        "4l_wwz_sf_A": 1.32*(0.281/0.2529),
        "4l_wwz_sf_B": 1.12*(0.281/0.2529),
        "4l_wwz_sf_C": 0.269*(0.281/0.2529),
        "4l_wwz_of_1": 0.322*(0.281/0.2529),
        "4l_wwz_of_2": 0.382*(0.281/0.2529),
        "4l_wwz_of_3": 0.833*(0.281/0.2529),
        "4l_wwz_of_4": 2.51*(0.281/0.2529),
    },
    "ZZ": {
        "4l_wwz_sf_A": 1.32,
        "4l_wwz_sf_B": 4.58,
        "4l_wwz_sf_C": 2.78,
        "4l_wwz_of_1": 0.623,
        "4l_wwz_of_2": 0.619,
        "4l_wwz_of_3": 0.39,
        "4l_wwz_of_4": 0.503,
    },
    "other": {
        "4l_wwz_sf_A": 0.607 + 0.00855, 
        "4l_wwz_sf_B": 0.417 + 0.0571,
        "4l_wwz_sf_C": 0.0967 + 0.0106,
        "4l_wwz_of_1": 0.27 + 0.107 + 0.0069,
        "4l_wwz_of_2": -0.055 + 0.0 + 0.0106,
        "4l_wwz_of_3": -0.0555 + 0.176 + 0.0446,
        "4l_wwz_of_4": 0.452 + 0.12 + 0.0233,
    },
}


sample_dict_base = {
    "WWZ" : ["WWZJetsTo4L2Nu"],
    "ZH"  : ["GluGluZH","ggToZHToZTo2L"],

    #"qqZZ": ["ZZTo4l"],
    #"ggZZ": ["ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau"],
    "ZZ"  : ["ZZTo4l", "ggToZZTo2e2mu", "ggToZZTo2e2tau", "ggToZZTo2mu2tau", "ggToZZTo4e", "ggToZZTo4mu", "ggToZZTo4tau"],
    #"ZZ"  : ["ZZTo4l"],

    "ttZ" : [
        "TTZToLL_M_1to10",
        "TTZToLLNuNu_M_10",
        "TTZToQQ",
    ],

    "other" : [
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
    ],

    #"WWZJetsTo4L2Nu":            ["WWZJetsTo4L2Nu"],
    #"GluGluZH":                  ["GluGluZH"],
    #"ggToZHToZTo2L":             ["ggToZHToZTo2L"],
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

# Print yields
def print_yields(yld_dict):

    # Dump the yields to dict for latex table
    yld_dict_for_printing = {}
    cats_to_print = [ "4l_wwz_sf_A", "4l_wwz_sf_B", "4l_wwz_sf_C", "4l_wwz_of_1", "4l_wwz_of_2", "4l_wwz_of_3", "4l_wwz_of_4"]
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
        #print_errs=True,
        column_variable="subkeys",
    )
    #exit()


    # Compare with other yields, print comparison

    #tag1 = "ewkcoffea (pfmet)"
    #tag2 = "ewkcoffea (puppi without 70-\>65)"
    tag1 = "ewkcoffea (puppi)"
    tag2 = "VVVNanoLooper (puppi)"

    yld_dict_comp = utils.put_none_errs(KEEGAN_YIELDS)
    #yld_dict_comp = EWK_PUPPIMET

    yld_dict_1 = yld_dict_for_printing
    yld_dict_2 = yld_dict_comp

    pdiff_dict = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="percent_diff",inpercent=True)
    diff_dict  = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="absolute_diff")

    utils.print_yld_dicts(yld_dict_1,tag1)
    utils.print_yld_dicts(yld_dict_2,tag2)
    utils.print_yld_dicts(pdiff_dict,f"Percent diff between {tag1} and {tag2}")
    utils.print_yld_dicts(diff_dict,f"Diff between {tag1} and {tag2}")

    mlt.print_begin()
    mlt.print_latex_yield_table(yld_dict_1,key_order=sample_dict_base.keys(),subkey_order=cats_to_print,tag=tag1)
    mlt.print_latex_yield_table(yld_dict_2,key_order=sample_dict_base.keys(),subkey_order=cats_to_print,tag=tag2)
    mlt.print_latex_yield_table(pdiff_dict,key_order=sample_dict_base.keys(),subkey_order=cats_to_print,tag=f"Percent diff between {tag1} and {tag2}")
    mlt.print_latex_yield_table(diff_dict, key_order=sample_dict_base.keys(),subkey_order=cats_to_print,tag=f"Diff between {tag1} and {tag2}")
    mlt.print_end()


# Dump the counts dict to a latex table
def print_counts(counts_dict):

    cats_to_print = ["all_events", "4l_presel", "4l_wwz_sf_A", "4l_wwz_sf_B", "4l_wwz_sf_C", "4l_wwz_of_1", "4l_wwz_of_2", "4l_wwz_of_3", "4l_wwz_of_4"]

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


# DRAFTING
################### Plotting ###################

import hist

def group(h, oldname, newname, grouping):
    hnew = hist.Hist(
        hist.axis.StrCategory(grouping, name=newname),
        *(ax for ax in h.axes if ax.name != oldname),
        storage=h._storage_type,
    )
    for i, indices in enumerate(grouping.values()):
        hnew.view(flow=True)[i] = h[{oldname: indices}][{oldname: sum}].view(flow=True)

    return hnew

# Takes a hist with one sparse axis and one dense axis, overlays everything on the sparse axis
def make_single_fig(histo,overlay_name,unit_norm_bool=False):
    #print("\nPlotting values:",histo.values())
    fig, ax = plt.subplots(1, 1, figsize=(7,7))
    histo.plot1d(
        overlay=overlay_name,
        stack=True,
        histtype="fill",
    )
    ax.autoscale(axis='y')
    plt.legend()
    return fig

def make_plots(histo_dict):

    save_dir_path = "plots"

    for var_name in histo_dict.keys():
        print("\n",var_name)
        if var_name == "ptl4": continue

        histo = histo_dict[var_name]

        grouping = sample_dict
        print("grouping",grouping)
        histo_grouped = group(histo,"process","process_grp",grouping)

        for cat_name in histo_grouped.axes["category"]:
            print(cat_name)
            histo_cat = histo_grouped[{"category":cat_name}]
            histo_cat.plot1d(overlay="process_grp")

            fig = make_single_fig(histo_cat,overlay_name="process_grp")
            title = f"{cat_name}_{var_name}"
            fig.savefig(os.path.join(save_dir_path,title))




################### Main ###################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default=".", help = "The path the output files should be saved to")
    parser.add_argument("-n", "--output-name", default="counts_wwz_sync", help = "A name for the output directory")
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))

    # Wrapper around the code for getting the raw counts and dump to latex table
    #counts_dict = get_yields(histo_dict,raw_counts=True)
    #print_counts(counts_dict)
    #exit()

    # Wrapper around the code for getting the yields for sr and bkg samples
    yld_dict = get_yields(histo_dict)
    print_yields(yld_dict)

    # Test plotting
    #make_plots(histo_dict)

    # Dump yield dict to json
    if "json" not in args.output_name: output_name = args.output_name + ".json"
    else: output_name = args.output_name
    with open(output_name,"w") as out_file: json.dump(yld_dict, out_file, indent=4)
    print(f"\nSaved json file: {output_name}\n")



if __name__ == "__main__":
    main()

