#!/bin/env python

import ewkcoffea.modules.sample_groupings as sg
import ewkcoffea.modules.keynote_tools as kt
import ewkcoffea.modules.yield_tools as yt

import argparse
import pickle
import gzip

# Get the yields in the SR
def get_yields(histos_dict,sample_dict,raw_counts=False,quiet=True,blind=True,systematic_name="nominal"):

    yld_dict = {}

    # Look at the yields in one histo (e.g. njets)
    if raw_counts: dense_axis = "njets_counts"
    else: dense_axis = "njets"
    for proc_name in sample_dict.keys():
        yld_dict[proc_name] = {}
        for cat_name in histos_dict[dense_axis].axes["category"]:
            #if "bdt" in cat_name: continue # TMP!!!
            if blind and (("data" in proc_name) and (not cat_name.startswith("cr_"))):
                # If this is data and we're not in a CR category, put placeholder numbers for now
                yld_dict[proc_name][cat_name] = [-999,-999]
            else:
                val = sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name],"systematic":systematic_name}].values(flow=True)))
                var = sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name],"systematic":systematic_name}].variances(flow=True)))
                yld_dict[proc_name][cat_name] = [val,var]

    # Print to screen
    if not quiet:
        for proc in yld_dict.keys():
            print(f"\n{proc}:")
            for cat in yld_dict[proc].keys():
                val = yld_dict[proc][cat]
                print(f"\t{cat}: {val}")

    return yld_dict

# Regroup categories (e.g. processes)
def group(h, oldname, newname, grouping):

    # Build up a grouping dict that drops any proc that is not in our h
    grouping_slim = {}
    proc_lst = get_axis_cats(h,oldname)
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
        storage=h.storage_type(),
    )
    for i, indices in enumerate(grouping_slim.values()):
        hnew.view(flow=True)[i] = h[{oldname: indices}][{oldname: sum}].view(flow=True)

    return hnew


def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default="plots", help = "The path the output files should be saved to")
    parser.add_argument('-y', "--get-yields", action='store_true', help = "Get yields from the pkl file")
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    parser.add_argument('-b', "--get-backgrounds", action='store_true', help = "Get background estimations")
    parser.add_argument('-u', "--ul-year", default='run2', help = "Which year to process", choices=["all","run2","run3","UL16APV","UL16","UL17","UL18","2022","2022EE"])
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))

    sample_dict_mc = sg.create_mc_sample_dict(args.ul_year)
    sample_dict_data = sg.create_data_sample_dict(args.ul_year)
    out_path = args.output_path

    yields_dict = yt.get_yields(histo_dict["njets_counts"],sample_dict_mc)
    kt.print_table(
        yields_dict,
        options={
            "output_name": "raw",
            "yield_prec": 0, # to print raw counts
        },
        sigs=["WWZ", "ZH"],
        bkgs=['ZZ', 'ttZ', 'tWZ', 'WZ', 'other'],
    )

    yields_dict = yt.get_yields(histo_dict["njets"],sample_dict_mc)
    yields_data_dict = yt.get_yields(histo_dict["njets"],sample_dict_data)


    kt.print_table(
        yields_dict,
        options={
            "output_name": "all",
        },
        sigs=["WWZ", "ZH"],
        bkgs=['ZZ', 'ttZ', 'tWZ', 'WZ', 'other'],
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "cutflow_cutbased",
            "regions": ["all_events", "4l_presel", "sr_4l_of_incl", "sr_4l_sf_incl"],
        },
        sigs=["WWZ", "ZH"],
        bkgs=['ZZ', 'ttZ', 'tWZ', 'WZ', 'other'],
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "sr_4l_sf",
            "regions": ["sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C"],
        },
        sigs=["WWZ", "ZH"],
        bkgs=['ZZ', 'ttZ', 'tWZ', 'WZ', 'other'],
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "sr_4l_of",
            "regions": ["sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"],
        },
        sigs=["WWZ", "ZH"],
        bkgs=['ZZ', 'ttZ', 'tWZ', 'WZ', 'other'],
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "cr",
            "regions": ["cr_4l_btag_of", "cr_4l_btag_sf_offZ_met80", "cr_4l_sf"],
        },
        sigs=["WWZ", "ZH"],
        bkgs=['ZZ', 'ttZ', 'tWZ', 'WZ', 'other'],
        yields_dict_data=yields_data_dict,
        data="data",
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "cutflow_bdt",
            "regions": ["all_events", "4l_presel", "sr_4l_bdt_of_presel", "sr_4l_bdt_sf_presel", "sr_4l_bdt_sf_trn"],
        },
        sigs=["WWZ", "ZH"],
        bkgs=['ZZ', 'ttZ', 'tWZ', 'WZ', 'other'],
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "sr_4l_bdt_sf",
            "regions": [f"sr_4l_bdt_sf_{i}" for i in range(1, 8)],
        },
        sigs=["WWZ", "ZH"],
        bkgs=['ZZ', 'ttZ', 'tWZ', 'WZ', 'other'],
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "sr_4l_bdt_of",
            "regions": [f"sr_4l_bdt_of_{i}" for i in range(1, 9)],
        },
        sigs=["WWZ", "ZH"],
        bkgs=['ZZ', 'ttZ', 'tWZ', 'WZ', 'other'],
    )

if __name__ == "__main__":
    main()
