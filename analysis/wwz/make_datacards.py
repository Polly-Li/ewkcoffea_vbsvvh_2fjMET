import os
import json
import pickle
import gzip
import argparse

from get_wwz_yields import get_yields # Note the fact that we're using this function means it probably belongs in ewkcoffea/ewkcoffea/modules


# Global variables
PRECISION = 6   # Decimal point precision in the text datacard output
PROC_LST = ["WWZ","ZH","ZZ","ttZ","other"]
SIG_LST = ["WWZ","ZH"]
CAT_LST = ["4l_wwz_sf_A", "4l_wwz_sf_B", "4l_wwz_sf_C", "4l_wwz_of_1", "4l_wwz_of_2", "4l_wwz_of_3", "4l_wwz_of_4"]


# Sum the predicted yields over categorires to get asimov data number
def get_fake_data_for_ch(yld_dict,ch):
    data_obs = 0
    for proc in yld_dict.keys():
        data_obs += yld_dict[proc][ch][0]
    return data_obs


# Make the datacard for a given channel
def make_ch_card(ch_ylds,ch,out_dir):

    # Building blocks we'll need to build the card formatting
    bin_str = f"bin_{ch}"
    syst_width = 0
    col_width = max(PRECISION*2+5,len(bin_str))
    line_break = "##----------------------------------\n"
    left_width = len(line_break) + 2
    left_width = max(syst_width+len("shape")+1,left_width)

    # The output name, location
    outf_card_name = f"test_card_{ch}.txt"
    print(f"Generating text file: {outf_card_name}")
    outf_card_name = os.path.join(out_dir,outf_card_name)

    # Create the card for this channel
    with open(outf_card_name,"w") as f:

        # Shapes rows, not sure of the purpose of the shapes lines when we have no shape templates
        #f.write("shapes *        * {ch} FAKE\n")
        f.write("shapes *        * FAKE\n")
        f.write(line_break)
        f.write(f"bin         {bin_str}\n")
        f.write(f"observation {ch_ylds['data_obs']}\n")
        f.write(line_break)
        f.write(line_break)

        # Note: This list is what controls the columns in the text datacard, if a process appears
        proc_order = PROC_LST

        # Bin row
        row = [f"{'bin':<{left_width}}"] # Python string formatting is pretty great!
        for p in proc_order:
            row.append(f"{bin_str:>{col_width}}")
        row = " ".join(row) + "\n"
        f.write(row)

        # 1st process row
        row = [f"{'process':<{left_width}}"]
        for p in proc_order:
            row.append(f"{p:>{col_width}}")
        row = " ".join(row) + "\n"
        f.write(row)

        # 2nd process row
        row = [f"{'process':<{left_width}}"]
        bkgd_count =  1
        sgnl_count = -1
        for p in proc_order:
            if any([x in p for x in SIG_LST]): # Check for if the process is signal or not
                row.append(f"{sgnl_count:>{col_width}}")
                sgnl_count += -1
            else:
                row.append(f"{bkgd_count:>{col_width}}")
                bkgd_count += 1
        row = " ".join(row) + "\n"
        f.write(row)

        # Rate row
        row = [f"{'rate':<{left_width}}"]
        for p in proc_order:
            r = ch_ylds[p][0]
            if r < 0:
                raise Exception(f"\nERROR: Process {p} has negative total rate: {r:.3f}.\n")
            row.append(f"{r:>{col_width}.{PRECISION}f}")
        row = " ".join(row) + "\n"

        # Write to the file
        f.write(row)
        f.write(line_break)


def main():

    # Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file_name",help="Either json file of yields or pickle file with scikit hists")
    parser.add_argument("--out-dir","-d",default="./testcards",help="Output directory to write root and text datacard files to")
    parser.add_argument("--do-nuisance",action="store_true",help="Include nuisance parameters")
    parser.add_argument("--unblind",action="store_true",help="If set, use real data, otherwise use asimov data")

    args = parser.parse_args()
    in_file = args.in_file_name
    out_dir = args.out_dir
    do_nuis = args.do_nuisance
    unblind = args.unblind

    # Check args
    if do_nuis:
        raise Exception("Nuisance params not implimented yet.")
    if out_dir != "." and not os.path.exists(out_dir):
        print(f"Making dir \"{out_dir}\"")
        os.makedirs(out_dir)

    # Get the yields dict from the input file
    #     - We can load a scikit hist (produced by wwz4l.py processor) and get the yields from that, dumpt to a dict
    #     - We can also load a json that contains the yields directly
    if in_file.endswith(".pkl.gz"):
        f = pickle.load(gzip.open(in_file))
        yld_dict = get_yields(f)
    elif in_file.endswith(".json"):
        with open(in_file) as f:
            yld_dict = json.load(f)
    else:
        raise Exception(f"ERROR: This script can only take hists or jsons, not files of type \"{in_file.split('.')[-1]}\".")


    # Make the cards for each channel
    print(f"Making cards for {CAT_LST}. \nPutting in {out_dir}.")
    for ch in CAT_LST:
        # Get yields for just this chan (also make fake asimov data if necessary)
        ch_ylds = {}
        for proc_name in yld_dict.keys():
            ch_ylds[proc_name] = yld_dict[proc_name][ch]
        if not unblind:
            ch_ylds["data_obs"] = get_fake_data_for_ch(yld_dict,ch) # Make asimov data for now
        else:
            raise Exception("We are not unblinded yet.")

        # Make the card for this chan
        make_ch_card(ch_ylds,ch,out_dir)


    print("Finished!")


if __name__ == "__main__":
    main()
