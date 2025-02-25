import pickle
import argparse
import gzip
import json

samples_dict = {
    "WWZ": {
        "samples": [
            "WWZJetsTo4L2Nu"
        ],
        "value": 0
    },
    "ZH": {
        "samples": [
            "qqToZHTo2WTo2L2Nu",
            "GluGluZHTo2WTo2L2Nu",
            "qqToZHToZTo2L",
            "GluGluZH"
        ],
        "value": 1
    },
    "ZZ": {
        "samples": [
            "ggToZZTo2e2mu",
            "ggToZZTo4mu",
            "ggToZZTo4e",
            "ggToZZTo2mu2tau",
            "ggToZZTo2e2tau",
            "ZZTo4l",
            "ggToZZTo4tau"
        ],
        "value": 2
    },
    "ttZ": {
        "samples": [
            "TTZToLL_M_50",
            "TTZToLL_M_4to50",
            "TTZToLLNuNu_M_10"
        ],
        "value": 3
    },
    "tWZ": {
        "samples": [
            "tWZ4l"
        ],
        "value": 4
    },
    "WZ": {
        "samples": [
            "WZTo3LNu"
        ],
        "value": 5
    },
    "Other": {
        "samples": [
            "TTTo2L2Nu",
            "ggHToZZ4L",
            "ttHnobb",
            "tZq",
            "ZZZ",
            "WZZ",
            "VHnobb",
            "tW_leptonic",
            "WWW",
            "TTWJetsToLNu",
        ],
        "value": 6
    },
    "Data": {
        "samples": [
            "data"
        ],
        "value": -1
    }
}


def update_dict():
    with open(__file__, "r+") as file:
        lines = file.readlines()

        # Locate the samples_dict definition
        start_idx = None
        end_idx = None
        for i, line in enumerate(lines):
            if line.strip() == "samples_dict = {":
                start_idx = i
            if start_idx is not None and line.strip() == "}":
                end_idx = i+1
                break

        if start_idx is not None and end_idx is not None:
            # Replace the dictionary definition
            new_dict = json.dumps(samples_dict, indent=4)
            lines[start_idx:end_idx + 1] = [f"samples_dict = {new_dict}\n"]

            # Write the updated file
            file.seek(0)
            file.writelines(lines)
            file.truncate()

def search_and_update(histo_dict):

    print("Searching for unrecognized samples")

    big_list = histo_dict["list_bdt_sf_proc"] + histo_dict["list_bdt_of_proc"]

    for i in range(len(big_list)):
        sample_name = big_list[i]
        sample_split_lst = sample_name.split("_",1)
        if len(sample_split_lst) != 2:
            raise Exception(f"Samples should have one underscore! {sample_name} does not follow the convention.")
        sample = sample_split_lst[1]
        # Check to see if the sample is already in the dict and append if not
        if "WWZJetsTo4L2Nu" in sample:
            if not (sample in samples_dict["WWZ"]["samples"]):
                print(f"{sample} is being added to WWZ category")
                samples_dict["WWZ"]["samples"].append(sample)
        elif "ZH" in sample:
            if not (sample in samples_dict["ZH"]["samples"]):
                print(f"{sample} is being added to ZH category")
                samples_dict["ZH"]["samples"].append(sample)
        elif "ZZTo" in sample:
            if not (sample in samples_dict["ZZ"]["samples"]):
                print(f"{sample} is being added to ZZ category")
                samples_dict["ZZ"]["samples"].append(sample)
        elif "TTZ" in sample:
            if not (sample in samples_dict["ttZ"]["samples"]):
                print(f"{sample} is being added to ttZ category")
                samples_dict["ttZ"]["samples"].append(sample)
        elif "tWZ" in sample:
            if not (sample in samples_dict["tWZ"]["samples"]):
                print(f"{sample} is being added to tWZ category")
                samples_dict["tWZ"]["samples"].append(sample)
        elif "WZTo" in sample:
            if not (sample in samples_dict["WZ"]["samples"]):
                print(f"{sample} is being added to WZ category")
                samples_dict["WZ"]["samples"].append(sample)
        elif "data" in sample:
            if not (sample in samples_dict["Data"]["samples"]):
                print(f"{sample} is being added to Data category")
                samples_dict["Data"]["samples"].append(sample)
        else:
            if not (sample in samples_dict["Other"]["samples"]):
                print(f"{sample} is being added to Other category")
                samples_dict["Other"]["samples"].append(sample)

    update_dict()


def get_proc_val(proc):
    # Check the proc is expected format
    sample_split_lst = proc.split("_",1)
    if len(sample_split_lst) != 2:
        raise Exception(f"Samples should have one underscore! {proc} does not follow the convention.")
    sample = sample_split_lst[1]
    for category, category_info in samples_dict.items():
        if sample in category_info["samples"]:
            return category_info["value"]
    raise Exception(f"{proc} not found! Try updating the dictionary!")


def create_txt_file(histo_dict, flavor, out_path):

    # Checking the output name
    if not (".txt" in out_path):
        out_path = out_path + ".txt"

    # Checking the flavor
    if flavor not in ["of","sf"]:
        raise Exception("Flavor must be sf or of!")

    # Getting the appropriate lists
    elif flavor == "sf":
        proc = histo_dict["list_bdt_sf_proc"]
        wwz_score = histo_dict["list_bdt_sf_wwz"]
        zh_score = histo_dict["list_bdt_sf_zh"]
        bkg_score = histo_dict["list_bdt_sf_bkg"]
        wgt = histo_dict["list_bdt_sf_wgt"]
    else:
        proc = histo_dict["list_bdt_of_proc"]
        wwz_score = histo_dict["list_bdt_of_wwz"]
        zh_score = histo_dict["list_bdt_of_zh"]
        bkg_score = histo_dict["list_bdt_of_bkg"]
        wgt = histo_dict["list_bdt_of_wgt"]

    # Checking the list lengths
    if not (len(proc) == len(wwz_score) == len(zh_score) == len(bkg_score) == len(wgt)):
        raise Exception(f"Values for {flavor} are not the same length! This is unexpected.")

    with open(f"{out_path}", "w") as file:
        for i in range(len(proc)):
            weight_val = wgt[i]
            wwz_score_val = wwz_score[i]
            zh_score_val = zh_score[i]
            bkg_score_val = bkg_score[i]
            proc_val = get_proc_val(proc[i])

            line = f"{proc_val} {wwz_score_val} {zh_score_val} {bkg_score_val} {weight_val}\n"
            file.write(line)

    print(f"Output has been written to {out_path}")


def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default = None, help = "The path the output files should be saved to")
    parser.add_argument('-u', "--update-samples", action='store_true', help = "Update Sample dictionary")
    parser.add_argument('-f', "--flavor", default='of', help = "Which flavor to grab", choices=["of","sf"])
    parser.add_argument('-r', "--run", default='run2', help = "Which Run is this", choices=["run2","run3"])
    args = parser.parse_args()

    # Define the arguments
    out_path = args.output_path
    pkl_file = args.pkl_file_path
    up_dict = args.update_samples
    flavor = args.flavor
    run = args.run

    # Make an appropriate output name based on inputs if outname is empty
    if (out_path == None):
        out_path = run + "_" + flavor
    print(f"Output will be saved as {out_path}")

    # Grab the hist
    histo_dict = pickle.load(gzip.open(pkl_file))

    # Check this hist was made with the --siphon option
    if not (("list_bdt_sf_proc" in histo_dict.keys()) and ("list_bdt_of_proc" in histo_dict.keys())):
        raise Exception("PKL file must have been made with --siphon option in processing step!")

    if up_dict:
        search_and_update(histo_dict)

    create_txt_file(histo_dict, flavor, out_path)

if __name__ == "__main__":
    main()
