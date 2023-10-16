#!/bin/env python
import sys
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
#	sig_fnames=[
#		"MC.root",
#		],
	legend_labels=[
		"DY",
		"TT",
		"W",
                "WW",
		],
#	signal_labels=[
#                "WWZ_SM",
#                "WWZ_NRW",
#		"WZZ_SM",
#                "WZZ_NRW",
#                ],


	extraoptions={
                "no_overflow": True, 
		"print_yield": True, 
                "yield_prec": 4,
		"lumi_value": 3.055, 
#                "signal_scale": 65.22,
#                "fit_bkg": True,
		"nbins": 45,
		"yaxis_log": True,
#                "xaxis_label":"mll",
#                "yaxis_label":"p_{T}  ",
#                "z_axis_range":[0.,1.], 
                "ratio_range": [0., 2.],
#                "draw_option_2d": "colztext",
#                "bin_text_format" : ".3f",
#                "bin_text_size": 1.0,  
#                "yaxis_log": True,
		}, 

	dirname=dir_name,
	)







# #______________________________________________________________________________________________________________________
# def dump_plot(
#         fnames=[],
#         sig_fnames=[],
#         data_fname=None,
#         dirname="plots",
#         legend_labels=[],
#         legend_labels_tex=[],
#         signal_labels=None,
#         signal_labels_tex=None,
#         donorm=False,
#         filter_pattern="",
#         signal_scale=1,
#         extraoptions={},
#         usercolors=None,
#         do_sum=False,
#         output_name=None,
#         dogrep=False,
#         _plotter=plot_hist,
#         doKStest=False,
#         histmodfunc=None,
#         histxaxislabeloptions={},
#         skip2d=False,
#         data_syst=None,
#         bkg_syst=None,
#         sig_syst=None,
#         ):

#        self.recognized_options = {
#
#            # Canvas
#            "canvas_width": {"type": "Int", "desc": "width of TCanvas in pixel", "default": None, "kinds": ["1dratio","graph","2d"], },
#            "canvas_height": {"type": "Int", "desc": "height of TCanvas in pixel", "default": None, "kinds": ["1dratio","graph","2d"], },
#            "canvas_main_y1": {"type": "Float", "desc": "main plot tpad y1", "default": 0.18, "kinds": ["1dratio","graph","2d"], },
#            "canvas_main_topmargin": {"type": "Float", "desc": "ratio plot top margin", "default": None, "kinds": ["1dratio"], },
#            "canvas_main_bottommargin": {"type": "Float", "desc": "ratio plot bottom margin", "default": None, "kinds": ["1dratio"], },
#            "canvas_main_rightmargin": {"type": "Float", "desc": "ratio plot right margin", "default": None, "kinds": ["1dratio"], },
#            "canvas_main_leftmargin": {"type": "Float", "desc": "ratio plot left margin", "default": None, "kinds": ["1dratio"], },
#            "canvas_ratio_y2": {"type": "Float", "desc": "ratio tpad y2", "default": 0.19, "kinds": ["1dratio","graph","2d"], },
#            "canvas_ratio_topmargin": {"type": "Float", "desc": "ratio plot top margin", "default": None, "kinds": ["1dratio"], },
#            "canvas_ratio_bottommargin": {"type": "Float", "desc": "ratio plot bottom margin", "default": None, "kinds": ["1dratio"], },
#            "canvas_ratio_rightmargin": {"type": "Float", "desc": "ratio plot right margin", "default": None, "kinds": ["1dratio"], },
#            "canvas_ratio_leftmargin": {"type": "Float", "desc": "ratio plot left margin", "default": None, "kinds": ["1dratio"], },
#            "canvas_tick_one_side": {"type": "Boolean", "desc": "ratio plot left margin", "default": False, "kinds": ["1dratio"], },
#
#            # Legend
#            "legend_coordinates": { "type": "List", "desc": "4 elements specifying TLegend constructor coordinates", "default": [0.63,0.67,0.93,0.87], "kinds": ["1dratio","graph"], },
#            "legend_alignment": { "type": "String", "desc": "easy alignment of TLegend. String containing two words from: bottom, top, left, right", "default": "", "kinds": ["1dratio","graph"], },
#            "legend_smart": { "type": "Boolean", "desc": "Smart alignment of legend to prevent overlaps", "default": True, "kinds": ["1dratio"], },
#            "legend_border": { "type": "Boolean", "desc": "show legend border?", "default": True, "kinds": ["1dratio","graph"], },
#            "legend_rounded": { "type": "Boolean", "desc": "rounded legend border", "default": True, "kinds": ["1dratio"], },
#            "legend_scalex": { "type": "Float", "desc": "scale width of legend by this factor", "default": 1, "kinds": ["1dratio","graph"], },
#            "legend_scaley": { "type": "Float", "desc": "scale height of legend by this factor", "default": 1, "kinds": ["1dratio","graph"], },
#            "legend_opacity": { "type": "Float", "desc": "from 0 to 1 representing the opacity of the TLegend white background", "default": 0.5, "kinds": ["1dratio","graph"], },
#            "legend_ncolumns": { "type": "Int", "desc": "number of columns in the legend", "default": 1, "kinds": ["1dratio","graph"], },
#            "legend_column_separation": { "type": "Float", "desc": "column separation size", "default": None, "kinds": ["1dratio","graph"], },
#            "legend_percentageinbox": { "type": "Boolean", "desc": "show relative process contributions as %age in the legend thumbnails", "default": True, "kinds": ["1dratio"], },
#            "legend_datalabel": { "type": "String", "desc": "label for the data histogram in the legend", "default": "Data", "kinds": ["1dratio"], },
#
#            # Axes
#            "xaxis_log": { "type": "Boolean", "desc": "log scale x-axis", "default": False, "kinds": ["1dratio","graph","2d"], },
#            "yaxis_log": { "type": "Boolean", "desc": "log scale y-axis", "default": False, "kinds": ["1dratio","graph","2d"], },
#            "zaxis_log": { "type": "Boolean", "desc": "log scale z-axis", "default": False, "kinds": ["2d"], },
#
#            "xaxis_label": { "type": "String", "desc": "label for x axis", "default": "", "kinds": ["1dratio","graph","2d"], },
#            "yaxis_label": { "type": "String", "desc": "label for y axis", "default": "Events", "kinds": ["1dratio","graph","2d"], },
#            "zaxis_label": { "type": "String", "desc": "label for z axis", "default": "", "kinds": ["2d"], },
#
#            "xaxis_label_size_scale": { "type": "Float", "desc": "size of fonts for x axis", "default": 1.0, "kinds": ["1dratio","graph","2d"], },
#            "yaxis_label_size_scale": { "type": "Float", "desc": "size of fonts for y axis", "default": 1.0, "kinds": ["1dratio","graph","2d"], },
#            "zaxis_label_size_scale": { "type": "Float", "desc": "size of fonts for z axis", "default": 1.0, "kinds": ["2d"], },
#
#            "xaxis_title_size": { "type": "Float", "desc": "size of fonts for x axis title", "default": None, "kinds": ["1dratio","graph","2d"], },
#            "yaxis_title_size": { "type": "Float", "desc": "size of fonts for y axis title", "default": None, "kinds": ["1dratio","graph","2d"], },
#
#            "xaxis_title_offset": { "type": "Float", "desc": "offset of x axis title", "default": None, "kinds": ["1dratio","graph","2d"], },
#            "yaxis_title_offset": { "type": "Float", "desc": "offset of y axis title", "default": None, "kinds": ["1dratio","graph","2d"], },
#
#            "xaxis_label_offset_scale": { "type": "Float", "desc": "x axis tickmark labels offset", "default": 1.0, "kinds": ["1dratio","graph","2d"], },
#            "yaxis_label_offset_scale": { "type": "Float", "desc": "y axis tickmark labels offset", "default": 1.0, "kinds": ["1dratio","graph","2d"], },
#
#            "xaxis_tick_length_scale": { "type": "Float", "desc": "x axis tickmark length scale", "default": 1.0, "kinds": ["1dratio","graph","2d"], },
#            "yaxis_tick_length_scale": { "type": "Float", "desc": "y axis tickmark length scale", "default": 1.0, "kinds": ["1dratio","graph","2d"], },
#
#            "xaxis_moreloglabels": { "type": "Boolean", "desc": "show denser labels with logscale for x axis", "default": True, "kinds": ["1dratio","graph","2d"], },
#            "yaxis_moreloglabels": { "type": "Boolean", "desc": "show denser labels with logscale for y axis", "default": True, "kinds": ["1dratio","graph","2d"], },
#            "zaxis_moreloglabels": { "type": "Boolean", "desc": "show denser labels with logscale for z axis", "default": True, "kinds": ["1dratio","graph","2d"], },
#            "xaxis_noexponents": { "type": "Boolean", "desc": "don't show exponents in logscale labels for x axis", "default": False, "kinds": ["1dratio","graph","2d"], },
#            "yaxis_noexponents": { "type": "Boolean", "desc": "don't show exponents in logscale labels for y axis", "default": False, "kinds": ["1dratio","graph","2d"], },
#            "zaxis_noexponents": { "type": "Boolean", "desc": "don't show exponents in logscale labels for z axis", "default": False, "kinds": ["1dratio","graph","2d"], },
#
#            "yaxis_exponent_offset": { "type": "Float", "desc": "offset x10^n left or right", "default": 0.0, "kinds": ["1dratio"], },
#            "yaxis_exponent_vertical_offset": { "type": "Float", "desc": "offset x10^n up or down", "default": 0.0, "kinds": ["1dratio"], },
#
#            "yaxis_ndivisions": { "type": "Int", "desc": "SetNdivisions integer for y-axis", "default": 510, "kinds": ["1dratio", "graph", "2d"], },
#            "xaxis_ndivisions": { "type": "Int", "desc": "SetNdivisions integer for x-axis", "default": 510, "kinds": ["1dratio", "graph", "2d"], },
#
#            "xaxis_range": { "type": "List", "desc": "2 elements to specify x axis range", "default": [], "kinds": ["1dratio","graph","2d"], },
#            "yaxis_range": { "type": "List", "desc": "2 elements to specify y axis range", "default": [], "kinds": ["1dratio","graph","2d"], },
#            "zaxis_range": { "type": "List", "desc": "2 elements to specify z axis range", "default": [], "kinds": ["2d"], },
#
#            # Ratio
#            "ratio_name": { "type": "String", "desc": "name of ratio pad", "default": "Data/MC", "kinds": ["1dratio"], },
#            "ratio_name_size": { "type": "Float", "desc": "size of the name on the ratio pad (e.g. data/MC)", "default": 0.2, "kinds": ["1dratio"], },
#            "ratio_name_offset": { "type": "Float", "desc": "offset to the name of ratio pad", "default": 0.25, "kinds": ["1dratio"], },
#            "ratio_range": { "type": "List", "desc": "pair for min and max y-value for ratio; default auto re-sizes to 3 sigma range", "default": [-1,-1], "kinds": ["1dratio"], },
#            "ratio_horizontal_lines": { "type": "List", "desc": "list of y-values to draw horizontal line", "default": [1.], "kinds": ["1dratio"], },
#            "ratio_chi2prob": { "type": "Boolean", "desc": "show chi2 probability for ratio", "default": False, "kinds": ["1dratio"], },
#            "ratio_pull": { "type": "Boolean", "desc": "show pulls instead of ratios in ratio pad", "default": False, "kinds": ["1dratio"], },
#            "ratio_pull_numbers": { "type": "Boolean", "desc": "show numbers for pulls, and mean/sigma", "default": True, "kinds": ["1dratio"], },
#            "ratio_ndivisions": { "type": "Int", "desc": "SetNdivisions integer for ratio", "default": 505, "kinds": ["1dratio"], },
#            "ratio_numden_indices": { "type": "List", "desc": "Pair of numerator and denominator histogram indices (from `bgs`) for ratio", "default": None, "kinds": ["1dratio"], },
#            "ratio_binomial_errors": { "type": "Boolean", "desc": "Use binomial error propagation when computing ratio eror bars", "default": False, "kinds": ["1dratio"], },
#            "ratio_xaxis_title": { "type": "String", "desc": "X-axis label", "default": "", "kinds": ["1dratio"], },
#            "ratio_xaxis_title_size": { "type": "Float", "desc": "X-axis label size", "default": None, "kinds": ["1dratio"], },
#            "ratio_xaxis_title_offset": { "type": "FLoat", "desc": "X-axis label offset", "default": None, "kinds": ["1dratio"], },
#            "ratio_label_size": { "type": "Float", "desc": "X-axis label size", "default": 0., "kinds": ["1dratio"], },
#            "ratio_xaxis_label_offset": { "type": "Float", "desc": "offset to the x-axis labels (numbers)", "default": None, "kinds": ["1dratio"], },
#            "ratio_yaxis_label_offset": { "type": "Float", "desc": "offset to the y-axis labels (numbers)", "default": None, "kinds": ["1dratio"], },
#            "ratio_tick_length_scale": { "type": "Float", "desc": "Tick length scale of ratio pads", "default": 1.0, "kinds": ["1dratio"], },
#
#            # Overall
#            "title": { "type": "String", "desc": "plot title", "default": "", "kinds": ["1dratio","graph","2d"], },
#            "draw_points": { "type": "Boolean", "desc": "draw points instead of fill", "default": False, "kinds": ["1d","1dratio"], },
#            "draw_option_2d": { "type": "String", "desc": "hist draw option", "default": "colz", "kinds": ["2d"], },
#            "bkg_err_fill_style": { "type": "Int", "desc": "Error shade draw style", "default": 1001, "kinds": ["1d", "1dratio"], },
#            "bkg_err_fill_color": { "type": "Int", "desc": "Error shade color", "default": None, "kinds": ["1d", "1dratio"], },
#
#            # CMS things
#            "cms_label": {"type": "String", "desc": "E.g., 'Preliminary'; default hides label", "default": None, "kinds": ["1dratio","graph","2d"]},
#            "lumi_value": {"type": "String", "desc": "E.g., 35.9; default hides lumi label", "default": "", "kinds": ["1dratio","graph","2d"]},
#            "lumi_unit": {"type": "String", "desc": "Unit for lumi label", "default": "fb", "kinds": ["1dratio","graph","2d"]},
#
#            # Misc
#            "do_stack": { "type": "Boolean", "desc": "stack histograms", "default": True, "kinds": ["1dratio"], },
#            "palette_name": { "type": "String", "desc": "color palette: 'default', 'rainbow', 'susy', etc.", "default": "default", "kinds": ["2d"], },
#            "show_bkg_errors": { "type": "Boolean", "desc": "show error bar for background stack", "default": False, "kinds": ["1dratio"], },
#            "show_bkg_smooth": { "type": "Boolean", "desc": "show smoothed background stack", "default": False, "kinds": ["1dratio"], },
#            "bkg_sort_method": { "type": "Boolean", "desc": "how to sort background stack using integrals: 'unsorted', 'ascending', or 'descending'", "default": 'ascending', "kinds": ["1dratio"], },
#            "no_ratio": { "type": "Boolean", "desc": "do not draw ratio plot", "default": False, "kinds": ["1dratio"], },
#            "no_overflow": { "type": "Boolean", "desc": "do not draw overflow bins", "default": False, "kinds": ["1dratio"], },
#            "stack_signal": { "type": "Boolean", "desc": "stack signal histograms", "default": False, "kinds": ["1dratio"], },
#
#            "max_digits": { "type": "Int", "desc": "integer for max digits", "default": 5, "kinds" : ["1dratio", "graph", "2d"], },
#
#
#            "bin_text_size": { "type": "Float", "desc": "size of text in bins (TH2::SetMarkerSize)", "default": 1.7, "kinds": ["2d"], },
#            "bin_text_format": { "type": "String", "desc": "format string for text in TH2 bins", "default": ".1f", "kinds": ["2d"], },
#            "bin_text_smart": { "type": "Boolean", "desc": "change bin text color for aesthetics", "default": False, "kinds": ["2d"], },
#            "bin_text_format_smart": { "type": "String", "desc": "python-syntax format string for smart text in TH2 bins taking value and bin error", "default": "{0:.0f}#pm{1:.0f}", "kinds": ["2d"], },
#
#            "hist_line_none": { "type": "Boolean", "desc": "No lines for histograms, only fill", "default": False, "kinds": ["1dratio"], },
#            "hist_line_black": { "type": "Boolean", "desc": "Black lines for histograms", "default": False, "kinds": ["1dratio"], },
#            "hist_disable_xerrors": { "type": "Boolean", "desc": "Disable the x-error bars on data for 1D hists", "default": True, "kinds": ["1dratio"], },
#
#            "extra_text": { "type": "List", "desc": "list of strings for textboxes", "default": [], "kinds": [ "1dratio","graph"], },
#            "extra_text_size": { "type": "Float", "desc": "size for extra text", "default": 0.04, "kinds": [ "1dratio","graph"], },
#            "extra_text_xpos": { "type": "Float", "desc": "NDC x position (0 to 1) for extra text", "default": 0.3, "kinds": [ "1dratio","graph"], },
#            "extra_text_ypos": { "type": "Float", "desc": "NDC y position (0 to 1) for extra text", "default": 0.87, "kinds": [ "1dratio","graph"], },
#
#            "extra_lines": { "type": "List", "desc": "list of 4-tuples (x1,y1,x2,y2) for lines", "default": [], "kinds": [ "1dratio","graph"], },
#            "no_overflow": {"type":"Boolean","desc":"Do not plot overflow bins","default": False, "kinds" : ["1dratio"],},
#
#            # Fun
#            "us_flag": { "type": "Boolean", "desc": "show the US flag in the corner", "default": False, "kinds": ["1dratio","graph","2d"], },
#            "us_flag_coordinates": { "type": "List", "desc": "Specify flag location with (x pos, y pos, size)", "default": [0.68,0.96,0.06], "kinds": ["1dratio","graph","2d"], },
#
#            # Output
#            "output_name": { "type": "String", "desc": "output file name/path", "default": "plot.pdf", "kinds": ["1dratio","graph","2d"], },
#            "output_ic": { "type": "Boolean", "desc": "run `ic` (imgcat) on output", "default": False, "kinds": ["1dratio","graph","2d"], },
#            "output_jsroot": { "type": "Boolean", "desc": "output .json for jsroot", "default": False, "kinds": ["1dratio","graph","2d"], },
#            "output_diff_previous": { "type": "Boolean", "desc": "diff the new output file with the previous", "default": False, "kinds": ["1dratio","graph","2d"], },
#
#        }

