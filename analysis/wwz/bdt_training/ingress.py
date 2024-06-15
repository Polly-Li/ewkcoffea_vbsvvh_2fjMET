import argparse
import numpy as np
import pickle
import gzip
import cloudpickle
import ewkcoffea.modules.sample_groupings as sg
from ewkcoffea.modules.paths import ewkcoffea_path as ewkcoffea_path
from topcoffea.modules.get_param_from_jsons import GetParam

def prepare_bdt_training_data(pklfilepath):

    # Open up the pickle file
    d = pickle.load(gzip.open(pklfilepath))

    # Dictionary where we will store our master data
    myd = {}

    # At this point I do wonder whether the following should be in the params.json
    # (N.B. there is also a copy of this in the processor)
    list_output_names = [
        "bdt_of_wwz_list",
        "bdt_of_zh_list",
        "bdt_of_bkg_list",
        "bdt_of_proc_list",
        "bdt_of_wgt_list",
        "bdt_of_evt_list",
        "bdt_sf_wwz_list",
        "bdt_sf_zh_list",
        "bdt_sf_bkg_list",
        "bdt_sf_proc_list",
        "bdt_sf_wgt_list",
        "bdt_sf_evt_list",
    ]

    get_ec_param = GetParam(ewkcoffea_path("params/params.json"))

    # Parse the pickle file and store to myd as numpy arrays
    for list_output_name in list_output_names: myd[list_output_name] = np.array(d[list_output_name])
    for var in get_ec_param("of_bdt_var_lst"): myd["of_bdt_" + var] = np.array(d["of_bdt_" + var])
    for var in get_ec_param("sf_bdt_var_lst"): myd["sf_bdt_" + var] = np.array(d["sf_bdt_" + var])

    # Sample groupings to split the data into signals, and backgrounds
    sample_dict_mc = sg.create_mc_sample_dict("run2")

    # List of categories that will use BDT to label
    categories = ["WWZ", "ZH", "Bkg"]

    # Sample splittings
    sample_splittings = ["train", "test", "all"] # later we can add validate if needed

    # Dictionary where we will store our training and testing data
    bdt_data = {}

    # Creating data structure in bdt_data
    # axis 0: categories
    # axis 1: split
    # axis 2: bdt variables and weight
    for cat in categories:
        bdt_data[cat] = {}
        for split in sample_splittings:
            bdt_data[cat][split] = {}
            # BDT inputs
            for var in get_ec_param("of_bdt_var_lst"): bdt_data[cat][split]["of_bdt_" + var] = np.array([])
            for var in get_ec_param("sf_bdt_var_lst"): bdt_data[cat][split]["sf_bdt_" + var] = np.array([])
            # also weight
            bdt_data[cat][split]["of_bdt_weight"] = np.array([])
            bdt_data[cat][split]["sf_bdt_weight"] = np.array([])
            # also event number
            bdt_data[cat][split]["of_bdt_event"] = np.array([])
            bdt_data[cat][split]["sf_bdt_event"] = np.array([])
            # also TMVA results to compare
            bdt_data[cat][split]["tmva_bdt_of_wwz"] = np.array([])
            bdt_data[cat][split]["tmva_bdt_sf_wwz"] = np.array([])
            bdt_data[cat][split]["tmva_bdt_of_zh"] = np.array([])
            bdt_data[cat][split]["tmva_bdt_sf_zh"] = np.array([])
            bdt_data[cat][split]["tmva_bdt_of_bkg"] = np.array([])
            bdt_data[cat][split]["tmva_bdt_sf_bkg"] = np.array([])

    # Loop over sample groupings
    for proc in sample_dict_mc:

        # Grouping handling for background (all backgrounds are grouped as "Bkg")
        # proc_cat will be either "WWZ", "ZH", or "Bkg"
        proc_cat = proc
        if proc in ["ZZ", "ttZ", "tWZ", "WZ", "other"]:
            proc_cat = "Bkg"

        # Loop over all the groupings and fill the empty numpy array with the actual inputs
        for sample_name in sample_dict_mc[proc]:

            # Processing OF BDT inputs
            for var in get_ec_param("of_bdt_var_lst"):

                bdt_vname = "of_bdt_" + var

                # Get the original numpy array
                a = bdt_data[proc_cat]["all"][bdt_vname]

                # Get more variables for this sample_name
                # masking is done via "bdt_of/sf_proc_list" which holds sample_name (that is how it was saved in wwz4l.py)
                b = myd[bdt_vname][myd["bdt_of_proc_list"] == sample_name]

                # Concatenate
                c = np.concatenate((a, b))

                # Set it back to bdt_data where we will store
                bdt_data[proc_cat]["all"][bdt_vname] = c

            # We need to save weights as well
            bdt_data[proc_cat]["all"]["of_bdt_weight"] = np.concatenate((bdt_data[proc_cat]["all"]["of_bdt_weight"], myd["bdt_of_wgt_list"][myd["bdt_of_proc_list"] == sample_name]))

            # We need to save event numbers as well
            bdt_data[proc_cat]["all"]["of_bdt_event"] = np.concatenate((bdt_data[proc_cat]["all"]["of_bdt_event"], myd["bdt_of_evt_list"][myd["bdt_of_proc_list"] == sample_name]))

            # We need to save tmva variables as well
            bdt_data[proc_cat]["all"]["tmva_bdt_of_wwz"] = np.concatenate((bdt_data[proc_cat]["all"]["tmva_bdt_of_wwz"], myd["bdt_of_wwz_list"][myd["bdt_of_proc_list"] == sample_name]))
            bdt_data[proc_cat]["all"]["tmva_bdt_of_zh"] = np.concatenate((bdt_data[proc_cat]["all"]["tmva_bdt_of_zh"], myd["bdt_of_zh_list"][myd["bdt_of_proc_list"] == sample_name]))
            bdt_data[proc_cat]["all"]["tmva_bdt_of_bkg"] = np.concatenate((bdt_data[proc_cat]["all"]["tmva_bdt_of_bkg"], myd["bdt_of_bkg_list"][myd["bdt_of_proc_list"] == sample_name]))

            # Processing OF BDT inputs
            for var in get_ec_param("sf_bdt_var_lst"):

                bdt_vname = "sf_bdt_" + var

                # Get the original numpy array
                a = bdt_data[proc_cat]["all"][bdt_vname]

                # Get more variables for this sample_name
                # masking is done via "bdt_sf/sf_proc_list" which holds sample_name (that is how it was saved in wwz4l.py)
                b = myd[bdt_vname][myd["bdt_sf_proc_list"] == sample_name]

                # Concatenate
                c = np.concatenate((a, b))

                # Set it back to bdt_data where we will store
                bdt_data[proc_cat]["all"][bdt_vname] = c

            # We need to save weights as well
            bdt_data[proc_cat]["all"]["sf_bdt_weight"] = np.concatenate((bdt_data[proc_cat]["all"]["sf_bdt_weight"], myd["bdt_sf_wgt_list"][myd["bdt_sf_proc_list"] == sample_name]))

            # We need to save event numbers as well
            bdt_data[proc_cat]["all"]["sf_bdt_event"] = np.concatenate((bdt_data[proc_cat]["all"]["sf_bdt_event"], myd["bdt_sf_evt_list"][myd["bdt_sf_proc_list"] == sample_name]))

            # We need to save tmva variables as well
            bdt_data[proc_cat]["all"]["tmva_bdt_sf_wwz"] = np.concatenate((bdt_data[proc_cat]["all"]["tmva_bdt_sf_wwz"], myd["bdt_sf_wwz_list"][myd["bdt_sf_proc_list"] == sample_name]))
            bdt_data[proc_cat]["all"]["tmva_bdt_sf_zh"] = np.concatenate((bdt_data[proc_cat]["all"]["tmva_bdt_sf_zh"], myd["bdt_sf_zh_list"][myd["bdt_sf_proc_list"] == sample_name]))
            bdt_data[proc_cat]["all"]["tmva_bdt_sf_bkg"] = np.concatenate((bdt_data[proc_cat]["all"]["tmva_bdt_sf_bkg"], myd["bdt_sf_bkg_list"][myd["bdt_sf_proc_list"] == sample_name]))

    # Split the testing and training by event number (even event number vs. odd event number)
    for cat in categories:
        for key in bdt_data[cat]["all"].keys():
            if "of" in key:
                test = bdt_data[cat]["all"][key][bdt_data[cat]["all"]["of_bdt_event"] % 2 == 1] # Getting events with event number that is even
                train = bdt_data[cat]["all"][key][bdt_data[cat]["all"]["of_bdt_event"] % 2 == 0] # Getting events with event number that is odd
                bdt_data[cat]["test"][key] = test
                bdt_data[cat]["train"][key] = train
            if "sf" in key:
                test = bdt_data[cat]["all"][key][bdt_data[cat]["all"]["sf_bdt_event"] % 2 == 1] # Getting events with event number that is even
                train = bdt_data[cat]["all"][key][bdt_data[cat]["all"]["sf_bdt_event"] % 2 == 0] # Getting events with event number that is odd
                bdt_data[cat]["test"][key] = test
                bdt_data[cat]["train"][key] = train


    return bdt_data

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    args = parser.parse_args()

    bdt_data = prepare_bdt_training_data(args.pkl_file_path)

    # Gathering categories and variable names
    categories = bdt_data.keys()
    of_variables = []
    sf_variables = []
    for cat in bdt_data.keys():
        for key in bdt_data[cat]["train"].keys():
            if "weight" in key: # skip weight
                continue
            if "event" in key: # skip weight
                continue
            if "of_bdt_" in key:
                of_variables.append(key)
            if "sf_bdt_" in key:
                sf_variables.append(key)

    # Integer labeling (WWZ = 0, ZH = 1, Bkg = 2)
    label_d = {}
    label_d["WWZ"] = 0
    label_d["ZH"] = 1
    label_d["Bkg"] = 2

    # Remove duplicates and get only unique variable names
    # Since we looped over multiple categories same variable names were put in
    of_variables = list(set(of_variables))
    sf_variables = list(set(sf_variables))

    # Building XGBoost input (which needs to be in a flat matrix)
    X_train_of = [] # Input variables
    y_train_of = np.array([]) # Output labels
    w_train_of = np.array([]) # Output labels
    X_train_sf = [] # Input variables
    y_train_sf = np.array([]) # Output labels
    w_train_sf = np.array([]) # Output labels
    X_test_of = [] # Input variables
    y_test_of = np.array([]) # Output labels
    w_test_of = np.array([]) # Output labels
    X_test_sf = [] # Input variables
    y_test_sf = np.array([]) # Output labels
    w_test_sf = np.array([]) # Output labels

    # Allocate some empty arrays where we will store things
    for var in of_variables: X_train_of.append(np.array([])) # Allocate N variable worth of lists
    for var in sf_variables: X_train_sf.append(np.array([])) # Allocate N variable worth of lists
    for var in of_variables: X_test_of.append(np.array([])) # Allocate N variable worth of lists
    for var in sf_variables: X_test_sf.append(np.array([])) # Allocate N variable worth of lists

    # Loop over the variables and store
    for cat in bdt_data.keys():
        # Store the variables
        for ivar, var in enumerate(of_variables):
            X_train_of[ivar] = np.concatenate((X_train_of[ivar], bdt_data[cat]["train"][var]))
            X_test_of[ivar] = np.concatenate((X_test_of[ivar], bdt_data[cat]["test"][var]))
        # Take the last variable and get length and create labels
        y_train_of = np.concatenate((y_train_of, np.full(len(bdt_data[cat]["train"][var]), label_d[cat])))
        y_test_of = np.concatenate((y_test_of, np.full(len(bdt_data[cat]["test"][var]), label_d[cat])))
        w_train_of = np.concatenate((w_train_of, bdt_data[cat]["train"]["of_bdt_weight"]))
        w_test_of = np.concatenate((w_test_of, bdt_data[cat]["test"]["of_bdt_weight"]))
        for ivar, var in enumerate(sf_variables):
            X_train_sf[ivar] = np.concatenate((X_train_sf[ivar], bdt_data[cat]["train"][var]))
            X_test_sf[ivar] = np.concatenate((X_test_sf[ivar], bdt_data[cat]["test"][var]))
        # Take the last variable and get length and create labels
        y_train_sf = np.concatenate((y_train_sf, np.full(len(bdt_data[cat]["train"][var]), label_d[cat])))
        y_test_sf = np.concatenate((y_test_sf, np.full(len(bdt_data[cat]["test"][var]), label_d[cat])))
        w_train_sf = np.concatenate((w_train_sf, bdt_data[cat]["train"]["sf_bdt_weight"]))
        w_test_sf = np.concatenate((w_test_sf, bdt_data[cat]["test"]["sf_bdt_weight"]))

    # Turn them into numpy array
    X_train_of = np.array(X_train_of)
    X_train_sf = np.array(X_train_sf)
    X_test_of = np.array(X_test_of)
    X_test_sf = np.array(X_test_sf)

    # Transpose to have the events in rows and variables in columns
    X_train_of = X_train_of.T
    X_train_sf = X_train_sf.T
    X_test_of = X_test_of.T
    X_test_sf = X_test_sf.T

    print("Printing OF training sample numpy array shapes")
    print(X_train_of.shape)
    print(y_train_of.shape)
    print(w_train_of.shape)
    print("Printing SF training sample numpy array shapes")
    print(X_train_sf.shape)
    print(y_train_sf.shape)
    print(w_train_sf.shape)
    print("Printing OF testing sample numpy array shapes")
    print(X_test_of.shape)
    print(y_test_of.shape)
    print(w_test_of.shape)
    print("Printing SF testing sample numpy array shapes")
    print(X_test_sf.shape)
    print(y_test_sf.shape)
    print(w_test_sf.shape)

    #### TMVA

    # Building TMVA output numpy array (which needs to be in a flat matrix)
    y_train_tmva_of = []
    y_train_tmva_sf = []
    y_test_tmva_of = []
    y_test_tmva_sf = []

    # Allocate some empty arrays where we will store things (there are 3 scores in tmva)
    tmva_scores_categories = ["wwz", "zh", "bkg"]
    for var in range(len(tmva_scores_categories)): y_train_tmva_of.append(np.array([])) # Allocate N variable worth of lists
    for var in range(len(tmva_scores_categories)): y_train_tmva_sf.append(np.array([])) # Allocate N variable worth of lists
    for var in range(len(tmva_scores_categories)): y_test_tmva_of.append(np.array([])) # Allocate N variable worth of lists
    for var in range(len(tmva_scores_categories)): y_test_tmva_sf.append(np.array([])) # Allocate N variable worth of lists

    # Loop over the variables and store
    for cat in bdt_data.keys():
        # Store the variables
        for ivar, var in enumerate(tmva_scores_categories):
            y_train_tmva_of[ivar] = np.concatenate((y_train_tmva_of[ivar], bdt_data[cat]["train"][f"tmva_bdt_of_{var}"]))
            y_test_tmva_of[ivar]  = np.concatenate((y_test_tmva_of[ivar], bdt_data[cat]["test"][f"tmva_bdt_of_{var}"]))
            y_train_tmva_sf[ivar] = np.concatenate((y_train_tmva_sf[ivar], bdt_data[cat]["train"][f"tmva_bdt_sf_{var}"]))
            y_test_tmva_sf[ivar]  = np.concatenate((y_test_tmva_sf[ivar], bdt_data[cat]["test"][f"tmva_bdt_sf_{var}"]))

    y_train_tmva_of = np.array(y_train_tmva_of)
    y_test_tmva_of = np.array(y_test_tmva_of)
    y_train_tmva_sf = np.array(y_train_tmva_sf)
    y_test_tmva_sf = np.array(y_test_tmva_sf)

    y_train_tmva_of = y_train_tmva_of.T
    y_test_tmva_of = y_test_tmva_of.T
    y_train_tmva_sf = y_train_tmva_sf.T
    y_test_tmva_sf = y_test_tmva_sf.T

    dd = {}
    dd["X_train_of"] = X_train_of
    dd["y_train_of"] = y_train_of
    dd["w_train_of"] = w_train_of
    dd["X_train_sf"] = X_train_sf
    dd["y_train_sf"] = y_train_sf
    dd["w_train_sf"] = w_train_sf
    dd["X_test_of"] = X_test_of
    dd["y_test_of"] = y_test_of
    dd["w_test_of"] = w_test_of
    dd["X_test_sf"] = X_test_sf
    dd["y_test_sf"] = y_test_sf
    dd["w_test_sf"] = w_test_sf
    dd["y_train_tmva_of"] = y_train_tmva_of
    dd["y_test_tmva_of"] = y_test_tmva_of
    dd["y_train_tmva_sf"] = y_train_tmva_sf
    dd["y_test_tmva_sf"] = y_test_tmva_sf

    with gzip.open("bdt.pkl.gz", "wb") as fout:
        cloudpickle.dump(dd, fout)

if __name__ == "__main__":

    main()
