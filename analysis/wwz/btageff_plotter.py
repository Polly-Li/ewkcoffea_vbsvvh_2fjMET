import os
import pickle
import gzip
import argparse
import matplotlib.pyplot as plt
import numpy as np
import awkward as ak

from coffea import lookup_tools

# This script is an example of how to dump some info from the btag eff histos (does not actually currently make a plot)
# Example usage:
#   python btageff_plotter.py ../../ewkcoffea/data/btag_eff/btag_eff_ttZ_srpresel.pkl.gz -u run2

def make_2d_plot(h_eff_flav, title, flav):
    xedges = ak.flatten(h_eff_flav.axes.edges[0])
    yedges = ak.flatten(h_eff_flav.axes.edges[1])
    bin_values = h_eff_flav.values()

    # Print bin values for debugging
    print("Bin values:\n", bin_values)

    fig, ax = plt.subplots()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = ax.pcolormesh(X, Y, bin_values.T, cmap='viridis', shading='auto')

    # Determine color bar limits based on flav value
    if flav == 0:
        vmin, vmax = 0, 0.5
    elif flav == 4:
        vmin, vmax = 0.2, 0.9
    elif flav == 5:
        vmin, vmax = 0.7, 1
    else:
        raise ValueError(f"Unsupported flav value: {flav}. Supported values are 0, 4, or 5.")

    # Add color bar with customized range
    cbar = plt.colorbar(mesh, ax=ax, extend='both')  # extend='both' for color bar extensions
    cbar.set_label('Efficiency')
    cbar.set_ticks(np.linspace(vmin, vmax, num=5))  # Example of setting ticks

    # Set color bar limits directly when creating it
    mesh.set_clim(vmin=vmin, vmax=vmax)

    # Set axis labels
    ax.set_xlabel('PT')
    ax.set_ylabel('Eta')
    ax.set_title(title)

    # Annotate each bin with its value
    for i in range(len(yedges)-1):
        for j in range(len(xedges)-1):
            bin_value = f'{bin_values.T[i, j]:.2f}'
            x = (xedges[j] + xedges[j+1]) / 2
            y = (yedges[i] + yedges[i+1]) / 2
            print(f"Annotating bin at ({x}, {y}) with value {bin_value}")
            ax.text(x, y, bin_value, ha='center', va='center', color='white', fontweight='bold', fontsize=8,
                    bbox=dict(facecolor='black', alpha=0.5, edgecolor='none'))

    return fig

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument('-u', "--ul-year", default='run2', help = "Which year to process", choices=["run2","run3"])
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    parser.add_argument('-y', "--get-yields", action='store_true', help = "Make yields from the pkl file")
    args = parser.parse_args()

    year = args.ul_year
    

    # Get the counts from the input hiso
    histo = pickle.load(gzip.open(args.pkl_file_path))["ptabseta"]

    if year == "run2":
        pname_list = ["UL16APV_TTZToLLNuNu_M_10","UL16_TTZToLLNuNu_M_10","UL17_TTZToLLNuNu_M_10","UL18_TTZToLLNuNu_M_10"]
    elif year == "run3":
        pname_list = ["2022EE","2022"]
    else:
        raise Exception("Unkown Input Year!")

    for pname in pname_list:
        print(f"\n{pname}")

        # Create lookup object and evaluate eff
        # Copy pasted from ewkcoffea corrections.py
        wp = "L"
        if year == "run3":
            histo_proc = histo[{"process":f"{pname}_TTZToLL_M_4to50"}] + histo[{"process":f"{pname}_TTZToLL_M_50"}]
        elif year == "run2":
            histo_proc = histo[{"process":pname}]
        else:
            raise Exception("Unknown Input Year!")

        h_eff = histo_proc[{"tag":wp}] / histo_proc[{"tag":"all"}]
        vals = h_eff.values(flow=True)[1:,1:-1,:-1] # Pt (drop underflow), eta (drop under and over flow), flav (drop overflow, there is not underflow)
        h_eff_lookup = lookup_tools.dense_lookup.dense_lookup(vals, [ax.edges for ax in h_eff.axes])

        if args.get_yields:
            # Print the histo
            print(h_eff)
            print(vals)

            # Example evaluation
            pt = 50
            abseta = 1.5
            hf = 0
            eff = h_eff_lookup(pt,abseta,hf)
            print(f"eff for pt={pt}, eta={abseta}, flav={hf}: {eff}")

        if args.make_plots:
            save_dir = "btag_plots"
            os.makedirs(save_dir, exist_ok=True)
            for flav in [0,4,5]:
                h_eff_flav = h_eff[{'flavor':flav}]
                fig = make_2d_plot(h_eff_flav,f"{pname} Flavor:{flav}",flav)
                save_path = os.path.join(save_dir, f"{pname}_Flavor_{flav}.png")
                fig.savefig(save_path)
                plt.close(fig)
                #plt.show()



if __name__ == "__main__":
    main()
