import argparse
import pickle
import gzip
import json
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import copy
import warnings
warnings.filterwarnings("ignore", message="List indexing selection is experimental.*")
warnings.filterwarnings("ignore", message="All sumw are zero!*")
from itertools import combinations

#for sample name csv
sample_name_csv = "/data/userdata/pyli/projects/VVHjj/coffea/vvhjj_coffea/analysis/vbsvvh/variables/sample_names.csv"
import csv
from collections import defaultdict

import utils.plotting_tools as plt_tools
from variables.cutflow_config import cutflow_dict

default_cutflow = "MET_objsel"
CAT_LST = cutflow_dict[default_cutflow].keys()

clr_map = {
    'EWK': '#d55e00',
    'Other': '#e377c2',
    'QCD': '#8c564b',
    'ST': '#9467bd',
    'WJets': '#d62728',
    'ZJets': '#2ca02c',
    'ttbar': '#ff7f0e',
    'ttx': '#1f77b4',
}


def parse_sample_csv(csv_path):
    GRP_DICT_FULL = defaultdict(list)
    bkg_color_map = {}

    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["sample_year"] != "2016preVFP": #temp fix: change to get year automatically in the future
                continue
            cat = row["sample_category"].strip().lower()
            sample_type = row["sample_type"].strip()
            sample_name = row["sample_name"].strip()
            color = row["plotting_colour"].strip()

            if cat == "sig":
                GRP_DICT_FULL["Signal"].append(sample_name)
            elif cat == "bkg":
                GRP_DICT_FULL.setdefault(sample_type, []).append(sample_name)
                if sample_type not in bkg_color_map:
                    bkg_color_map[sample_type] = color

    # Get background color list in the same order as background keys
    bkg_CLR_LST = [bkg_color_map[grp] for grp in GRP_DICT_FULL if grp != "Signal"]

    return dict(GRP_DICT_FULL), bkg_CLR_LST

GRP_DICT_FULL, CLR_LST = parse_sample_csv(sample_name_csv)

# Get sig and bkg yield in all categories
def get_yields_per_cat(histo_dict,var_name):
    out_dict = {}

    # Get the initial grouping dict
    grouping_dict = GRP_DICT_FULL

    # Get list of all of the backgrounds together
    bkg_lst = []
    for grp in grouping_dict:
        if grp != "Signal":
            bkg_lst = bkg_lst + grouping_dict[grp]

    # Make the dictionary to get yields for, it includes what's in grouping_dict, plus the backgrounds grouped as one
    groups_to_get_yields_for_dict = copy.deepcopy(grouping_dict)
    groups_to_get_yields_for_dict["Background"] = bkg_lst

    # Loop over cats and fill dict of sig and bkg
    for cat in CAT_LST:
        out_dict[cat] = {}
        
        
        histo_base = histo_dict[var_name][{"systematic":"nominal", "category":cat}]
        histo_base = plt_tools.group(histo_base,"year", "temp_year", {"all_year":['2016preVFP','2016postVFP','2017','2018']})
        histo_base = histo_base[{"temp_year":"all_year"}]

        # Get values per proc
        for group_name,group_lst in groups_to_get_yields_for_dict.items():
            histo = plt_tools.group(histo_base,"process","process",{group_name:group_lst})
            yld = sum(sum(histo.values(flow=True)))
            var = sum(sum(histo.variances(flow=True)))
            out_dict[cat][group_name] = [yld,(var)**0.5]

        # Get the metric
        sig = out_dict[cat]["Signal"][0]
        bkg = out_dict[cat]["Background"][0]
        metric = sig/(bkg)**0.5
        punzi = sig/(1+(bkg)**0.5)
        punzi_p1 = sig/(0.11+(bkg)**0.5)
        out_dict[cat]["metric"] = [metric,None] # Don't bother propagating error
        out_dict[cat]["punzi"] = [punzi,None] # Don't bother propagating error
        out_dict[cat]["punzi_p1"] = [punzi_p1,None] # Don't bother propagating error

    return out_dict

def get_histogram_ABCD(original_hist, cat_arr, axis_name="category"):
    results = {}

    for cat1, cat2 in combinations(cat_arr, 2):
        label = (cat1,cat2)

        # Categories excluding both cat1 and cat2
        cats_except_both = [c for c in cat_arr if c not in [cat1, cat2]]
        # Categories excluding only cat1
        cats_except_cat1 = [c for c in cat_arr if c != cat1]
        # Categories excluding only cat2
        cats_except_cat2 = [c for c in cat_arr if c != cat2]
        # All categories
        cats_all = cat_arr

        results[label] = {
            "all_except_both": original_hist[{axis_name: cats_except_both}],
            "all_except_cat1": original_hist[{axis_name: cats_except_cat1}],
            "all_except_cat2": original_hist[{axis_name: cats_except_cat2}],
            "all": original_hist[{axis_name: cats_all}],
        }

    return results

def plot_hist2d_and_correlation(h, var1, var2, filename=None, logz=False):
    """
    h      : coffea.hist.Hist object with var1 and var2 as dense axes
    var1   : name of x-axis variable
    var2   : name of y-axis variable
    filename : optional output filename, defaults to '{var1}_{var2}_hist2d.png'
    logz   : use log scale on Z axis (colorbar)
    """

    # Project to 2D array
    h2d = h.project(var1, var2)
    values = h2d.values()

    # Bin centers
    x_axis = h2d.axes[var1]
    y_axis = h2d.axes[var2]
    x_centers = 0.5 * (x_axis.edges[:-1] + x_axis.edges[1:])
    y_centers = 0.5 * (y_axis.edges[:-1] + y_axis.edges[1:])
    X, Y = np.meshgrid(x_centers, y_centers, indexing='ij')

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    cmap = "viridis"
    c = ax.pcolormesh(X, Y, values.T, cmap=cmap, shading='auto',
                      norm=(plt.LogNorm() if logz else None))
    fig.colorbar(c, ax=ax, label='Entries')

    ax.set_xlabel(var1)
    ax.set_ylabel(var2)
    ax.set_title(f"2D Histogram: {var1} vs {var2}")

    # Save
    if filename is None:
        filename = f"{var1}_{var2}_hist2d.png"
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

    # Flattened data for correlation (weights)
    xvals = np.repeat(x_centers, len(y_centers))
    yvals = np.tile(y_centers, len(x_centers))
    weights = values.flatten()

    # Mask zero weights
    mask = weights > 0
    xvals, yvals, weights = xvals[mask], yvals[mask], weights[mask]

    # Weighted correlation
    cov = np.cov(xvals, yvals, aweights=weights)
    corr = cov[0, 1] / (np.sqrt(cov[0, 0]) * np.sqrt(cov[1, 1]))

    print(f"Pearson correlation ({var1}, {var2}): {corr:.4f}")
    return corr

# Make the figures for the vvh sudy


def group_hist(histo_dict,var,syst='nominal',group_dict=GRP_DICT_FULL):
    histo_dict = histo_dict[var][{"systematic":syst}]
    histo_dict = plt_tools.group(histo_dict,"year", "temp_year", {"all_year":['2016preVFP','2016postVFP','2017','2018']})
    histo_dict = histo_dict[{"temp_year":"all_year"}]
    histo_dict = plt_tools.group(histo_dict,"process","process_grp",group_dict)
    return histo_dict

### Get the sig and bkg yields and print or dump to json ###
def print_yields(histo_dict,roundat=None,print_counts=False,dump_to_json=True,quiet=False,out_name="yields"):

    # Get ahold of the yields
    yld_dict    = get_yields_per_cat(histo_dict,"nGoodAK4")
    #counts_dict = get_yields_per_cat(histo_dict,"nAK4_counts")

    #group_lst_order = ['Signal', 'Background', 'ttbar', 'VV', 'Vjets', 'QCD', 'single-t', 'ttX', 'VH', 'VVV']
    group_lst_order = list(GRP_DICT_FULL.keys()) + ['Background']
    print("debug",group_lst_order)

    # Print to screen
    if not quiet:

        ### Print readably ###
        print("\n--- Yields ---")
        for cat in yld_dict:
            print(f"\n{cat}")
            for group_name in group_lst_order:
                if group_name == "metric": continue
                yld, err = yld_dict[cat][group_name]
                perr = 100*(err/yld)
                print(f"    {group_name}:  {np.round(yld,roundat)} +- {np.round(perr,2)}%")
            print(f"    -> Metric: {np.round(yld_dict[cat]['metric'][0],3)}")


        ### Print csv, build op as an out string ###

        # Append the header
        out_str = ""
        header = "cat name"
        for proc_name in group_lst_order:
            header = header + f", {proc_name}_val , pm , {proc_name}_pm"
        header = header + ", metric"
        out_str = out_str + header

        # Appead a line for each category, with yields and metric
        for cat in yld_dict:
            line_str = cat
            for group_name in yld_dict[cat]:
                if group_name == "metric": continue
                if group_name == "punzi": continue
                if group_name == "punzi_p1": continue
                yld, err = yld_dict[cat][group_name]
                perr = 100*(err/yld)
                line_str = line_str + f" , {np.round(yld,roundat)} , ± , {np.round(perr,2)}%"
            # And also append the metric
            metric = yld_dict[cat]["metric"][0]
            line_str = line_str + f" , {np.round(metric,3)}"
            # Append the string for this line to the out string
            out_str = out_str + f"\n{line_str}"

        # Print the out string to the screen
        print("\n\n--- Yields CSV formatted ---\n")
        print(out_str)


    # Dump directly to json
    if dump_to_json:
        #out_dict = {"yields":yld_dict, "counts":counts_dict}
        out_dict = {"yields": yld_dict}
        output_name = f"{out_name}.json"  # e.g., 'results/my_output.json'

        # Ensure the parent directory exists
        output_dir = os.path.dirname(output_name)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Now write the file
        with open(output_name, "w") as out_file:
            json.dump(out_dict, out_file, indent=4)
        if not quiet:
            print("\n\n--- Yields json formatted ---")
            print(f"\nSaved json file: {output_name}\n")


def loop_ABCD(histo_dict,output_name):
    grouping_dict = GRP_DICT_FULL

    cat_lst = CAT_LST
    var_lst = histo_dict.keys()

    # Save dir
    save_dir_path = f"ABCD_plots/{output_name}"
    if not os.path.exists(f"./ABCD_plots/{output_name}"): os.makedirs(f"./ABCD_plots/{output_name}", exist_ok=True)
    
    csv_file = f'{save_dir_path}/ABCD_summary.csv'
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['var1','var2','bkg_correlation','pred','err','actual','sig_contam','total_bkg','total_sig','s/sqrt(b)']) #header

    for cat1,cat2 in combinations():
        hist = histo_dict[var]
            #if var not in ["njets","njets_counts","scalarptsum_lepmet"]: continue # TMP

            #histo = copy.deepcopy(histo_dict[var][{"systematic":"nominal", "category":cat}]) 
        histo = copy.deepcopy(histo_dict[var][{"category":cat}]) #keep systematic before group as bkg and sig may use different ones

        #can't sum directly over year
        #histo = histo.sum("year")

        # Clean up a bit (rebin, regroup, and handle overflow)

        # if var not in ["njets","nleps","nbtagsl","nbtagsm","njets_counts","nleps_counts","nfatjets","njets_forward","njets_tot"]:
        #    histo = plt_tools.rebin(histo,6)
        
        histo = plt_tools.group(histo,"process","process_grp",grouping_dict)


        # Get one hist of just sig and one of just bkg
        grp_names_bkg_lst = list(grouping_dict.keys()) # All names, still need to drop signal
        
        #print(grp_names_bkg_lst)

        grp_names_bkg_lst.remove("Signal")

        histo_sig = histo[{"process_grp":["Signal"], "systematic":sig_coupling}] #specify coupling for signal
        #histo = histo[{"systematic":"nominal"}] #poor naming, nominal for bkg is really just weight
        histo_bkg = plt_tools.group(histo[{"systematic":"nominal"}],"process_grp","process_grp",{"Background":grp_names_bkg_lst})
        #print("grp_names_bkg_lst",grp_names_bkg_lst)
        #print("histo key", histo)
        
        available_groups = list(histo.axes["process_grp"])


        # Filter only those group names that exist in both the list and the histo
        valid_bkg_groups = [grp for grp in grp_names_bkg_lst if grp in available_groups]
        # Now safely subset


        histo = plt_tools.merge_overflow(histo)
        hist_bkg_test = histo[{"systematic":"nominal","process_grp": valid_bkg_groups}]
        #same binning for plotting
        print("debug",var,cat)
        histo_bkg,rebin_factor = plt_tools.rebin_to_nonnegative(histo_bkg,"process_grp")
        if rebin_factor >1:
            hist_bkg_test = plt_tools.rebin(hist_bkg_test,rebin_factor)
            histo_sig = plt_tools.rebin(histo_sig,rebin_factor)

##################################### Main #####################################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file. must be a n-1 pkl")
    parser.add_argument('-y', "--get-yields", action='store_true', help = "Get yields from the pkl file")
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    parser.add_argument('-j', "--dump-json", action='store_true', help = "Dump some yield numbers into a json file")
    parser.add_argument('-o', "--output-name", default='vvh', help = "What to name the outputs")
    parser.add_argument('-c', "--cutflow", default=None, help = "cutflow to use")
    parser.add_argument('-s', "--sig_coupling", default=None, help = "sig_coupling (systematic axis)")

    args = parser.parse_args()
    
    if args.cutflow is not None:
        global CAT_LST
        CAT_LST = cutflow_dict[args.cutflow].keys()

    sig_coupling = "nominal" #default setting
    if args.sig_coupling is not None:
        sig_coupling = args.sig_coupling

    # Get the dictionary of histograms from the input pkl file
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))
    histo_dict = copy.deepcopy(histo_dict)
    var_lst = histo_dict.keys()
    for var in var_lst:
        histo_dict[var] = group_hist(histo_dict,var)
    plot_hist2d_and_correlation(histo_dict, 'nAK4', 'nAK8', filename='test')
    
    
    #print(histo_dict)
    #return coffea.processor.accumulator.dict_accumulator (accumulated hists from processor))

    # Print total raw events
    #tot_raw = sum(sum(histo_dict["njets_counts"][{"systematic":"nominal", "category":"all_events"}].values(flow=True)))
    #print("Tot raw events:",tot_raw)

    # Which main functionalities to run
    #print(f"Setting: cutflow-{list(CAT_LST)} weight-{sig_coupling}")
    #loop_ABCD(histo_dict,args.output_name)

main()