import pickle
import gzip
import argparse
import matplotlib.pyplot as plt

from coffea import lookup_tools

# This script is an example of how to dump some info from the btag eff histos (does not actually currently make a plot)
# Example usage:
#   python btageff_plotter.py ../../ewkcoffea/data/btag_eff/btag_eff_ttZ_srpresel.pkl.gz

def make_2d_fig(histo,title):
    fig, ax = plt.subplots(1, 1, figsize=(7,7))

    histo.plot2d(labels=True)

    plt.title(title)
    #ax.autoscale(axis='y')
    #plt.legend()
    return fig


def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo = pickle.load(gzip.open(args.pkl_file_path))["ptabseta"]

    #pname = "UL17_TTZToLLNuNu_M_10"
    #pname = "UL17_WWZJetsTo4L2Nu"
    #print(histo)
    #print(histo.values())

    for pname in ["UL16APV_TTZToLLNuNu_M_10","UL16_TTZToLLNuNu_M_10","UL17_TTZToLLNuNu_M_10","UL18_TTZToLLNuNu_M_10"]:
        print(f"\n{pname}")

        # Create lookup object and evaluate eff
        # Copy pasted from ewkcoffea corrections.py
        wp = "L"
        histo_proc = histo[{"process":pname}]
        h_eff = histo_proc[{"tag":wp}] / histo_proc[{"tag":"all"}]
        vals = h_eff.values(flow=True)[1:,1:-1,:-1] # Pt (drop underflow), eta (drop under and over flow), flav (drop overflow, there is not underflow)
        h_eff_lookup = lookup_tools.dense_lookup.dense_lookup(vals, [ax.edges for ax in h_eff.axes])

        # Print the histo
        print(h_eff)
        print(vals)

        # Example evaluation
        pt = 50
        abseta = 1.5
        hf = 0
        eff = h_eff_lookup(pt,abseta,hf)
        print(f"eff for pt={pt}, eta={abseta}, flav={hf}: {eff}")


if __name__ == "__main__":
    main()
