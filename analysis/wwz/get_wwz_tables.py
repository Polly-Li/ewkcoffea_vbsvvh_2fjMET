#!/bin/env python

import ewkcoffea.modules.sample_groupings as sg
import ewkcoffea.modules.keynote_tools as kt
import ewkcoffea.modules.yield_tools as yt

import argparse
import pickle
import gzip

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default="plots", help = "The path the output files should be saved to")
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
        sigs=sg.SIG_LST,
        bkgs=sg.BKG_LST,
    )

    yields_dict = yt.get_yields(histo_dict["njets"],sample_dict_mc)
    yields_data_dict = yt.get_yields(histo_dict["njets"],sample_dict_data)


    kt.print_table(
        yields_dict,
        options={
            "output_name": "all",
        },
        sigs=sg.SIG_LST,
        bkgs=sg.BKG_LST,
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "cutflow_cutbased",
            "regions": ["all_events", "4l_presel", "sr_4l_of_incl", "sr_4l_sf_incl"],
        },
        sigs=sg.SIG_LST,
        bkgs=sg.BKG_LST,
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "sr_4l_sf",
            "regions": ["sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C"],
        },
        sigs=sg.SIG_LST,
        bkgs=sg.BKG_LST,
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "sr_4l_of",
            "regions": ["sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"],
        },
        sigs=sg.SIG_LST,
        bkgs=sg.BKG_LST,
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "cr",
            "regions": ["cr_4l_btag_of", "cr_4l_btag_sf_offZ_met80", "cr_4l_sf"],
        },
        sigs=sg.SIG_LST,
        bkgs=sg.BKG_LST,
        yields_dict_data=yields_data_dict,
        data="data",
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "cutflow_bdt",
            "regions": ["all_events", "4l_presel", "sr_4l_bdt_of_presel", "sr_4l_bdt_sf_presel", "sr_4l_bdt_sf_trn"],
        },
        sigs=sg.SIG_LST,
        bkgs=sg.BKG_LST,
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "sr_4l_bdt_sf",
            "regions": [f"sr_4l_bdt_sf_{i}" for i in range(1, 8)],
        },
        sigs=sg.SIG_LST,
        bkgs=sg.BKG_LST,
    )

    kt.print_table(
        yields_dict,
        options={
            "output_name": "sr_4l_bdt_of",
            "regions": [f"sr_4l_bdt_of_{i}" for i in range(1, 9)],
        },
        sigs=sg.SIG_LST,
        bkgs=sg.BKG_LST,
    )

if __name__ == "__main__":
    main()
