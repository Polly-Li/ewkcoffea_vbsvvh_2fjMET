#!/bin/env python

import argparse
import plottery_wrapper as p

parser = argparse.ArgumentParser(description='Input Root Files and Directory Name')
parser.add_argument('-dir', '--plot_dir', type=str, help='Plot Directory Name')
parser.add_argument('-dy', '--DY_file', type=str, help='DY File Name')
parser.add_argument('-w', '--W_file', type=str, help='W File Name')
parser.add_argument('-ww', '--WW_file', type=str, help='WW File Name')
parser.add_argument('-tt', '--TT_file', type=str, help='TT File Name')
parser.add_argument('-dt', '--Data_file', type=str, help='Data File Name')

args = parser.parse_args()

dir_name = args.plot_dir
dy_file = args.DY_file
w_file = args.W_file
ww_file = args.WW_file
tt_file = args.TT_file
data_file = args.Data_file

p.dump_plot(
    fnames=[
        dy_file,
        tt_file,
        w_file,
        ww_file,
    ],
    data_fname = data_file,

    legend_labels=[
        "DY",
        "TT",
        "W",
        "WW",
    ],

    extraoptions={
        "no_overflow": True,
        "print_yield": True,
        "yield_prec": 4,
        "lumi_value": 35.1,
        #"yaxis_log": True,
        "ratio_range": [0.5, 1.5],
    },

    dirname=dir_name,
)
