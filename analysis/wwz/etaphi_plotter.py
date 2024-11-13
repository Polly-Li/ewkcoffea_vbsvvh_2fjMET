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

import numpy as np
import matplotlib.pyplot as plt

def make_2d_plot(hist, title, nbins=10):
    # Get the original x and y range
    x_min, x_max = ak.flatten(hist.axes.edges[0])[0], ak.flatten(hist.axes.edges[0])[-1]
    y_min, y_max = ak.flatten(hist.axes.edges[1])[0], ak.flatten(hist.axes.edges[1])[-1]
    
    # Define new bin edges using the same range but with fewer bins
    xedges = np.linspace(x_min, x_max, nbins + 1)
    yedges = np.linspace(y_min, y_max, nbins + 1)
    
    # Aggregate the original bin values to match the new binning
    bin_values = hist.values()
    new_bin_values = np.zeros((nbins, nbins))
    x_bin_width = bin_values.shape[0] // nbins
    y_bin_width = bin_values.shape[1] // nbins
    
    for i in range(nbins):
        for j in range(nbins):
            new_bin_values[i, j] = bin_values[i*x_bin_width:(i+1)*x_bin_width, j*y_bin_width:(j+1)*y_bin_width].sum()
    
    # Plot with the new binning
    fig, ax = plt.subplots()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = ax.pcolormesh(X, Y, new_bin_values.T, cmap='viridis', shading='auto')

    # Add color bar with customized range
    cbar = plt.colorbar(mesh, ax=ax, extend='both')
    cbar.set_label('Jet Counts')

    # Set axis labels and title
    ax.set_xlabel('ETA')
    ax.set_ylabel('PHI')
    ax.set_title(title)

    # Annotate each bin with its value
    #for i in range(nbins):
    #    for j in range(nbins):
    #        bin_value = f'{new_bin_values[i, j]:.2f}'
    #        x = (xedges[j] + xedges[j+1]) / 2
    #        y = (yedges[i] + yedges[i+1]) / 2
    #        ax.text(x, y, bin_value, ha='center', va='center', color='white', fontweight='bold', fontsize=8,
    #                bbox=dict(facecolor='black', alpha=0.5, edgecolor='none'))

    return fig



#def make_2d_plot(hist, title):
#    xedges = ak.flatten(hist.axes.edges[0])
#    yedges = ak.flatten(hist.axes.edges[1])
#    bin_values = hist.values()
#
#    # Print bin values for debugging
#    #print("Bin values:\n", bin_values)
#
#    fig, ax = plt.subplots()
#    X, Y = np.meshgrid(xedges, yedges)
#    mesh = ax.pcolormesh(X, Y, bin_values.T, cmap='viridis', shading='auto')
#
#    # Add color bar with customized range
#    cbar = plt.colorbar(mesh, ax=ax, extend='both')  # extend='both' for color bar extensions
#    cbar.set_label('Jet Counts')
#    #cbar.set_ticks(np.linspace(vmin, vmax, num=5))  # Example of setting ticks
#
#    # Set color bar limits directly when creating it
#    mesh.set_clim(vmin=0, vmax=3)
#
#    # Set axis labels
#    ax.set_xlabel('ETA')
#    ax.set_ylabel('PHI')
#    ax.set_title(title)
#
#    # Annotate each bin with its value
#    #for i in range(len(yedges)-1):
#    #    for j in range(len(xedges)-1):
#    #        bin_value = f'{bin_values.T[i, j]:.2f}'
#    #        x = (xedges[j] + xedges[j+1]) / 2
#    #        y = (yedges[i] + yedges[i+1]) / 2
#    #        print(f"Annotating bin at ({x}, {y}) with value {bin_value}")
#    #        ax.text(x, y, bin_value, ha='center', va='center', color='white', fontweight='bold', fontsize=8,
#    #                bbox=dict(facecolor='black', alpha=0.5, edgecolor='none'))
#
#    return fig


def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument('-u', "--ul-year", default='run2', help = "Which year to process", choices=["UL18","UL17","UL16","UL16APV"])
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    args = parser.parse_args()

    year = args.ul_year


    # Get the counts from the input hiso
    histo = pickle.load(gzip.open(args.pkl_file_path))["etaphi_all"]

    cat_list = [
#        "all_events",
        "cr_4l_sf",
#        "cr_4l_btag_of",
#        "cr_4l_btag_sf_offZ_met80"
    ]

    for cat in cat_list:
        print(f"\n{cat}")

        histo_proc = histo[{"category":cat}]

        if args.make_plots:
            save_dir = "etaphi_plots"
            os.makedirs(save_dir, exist_ok=True)

            fig = make_2d_plot(histo_proc,f"Year:{year} Category:{cat}")
            save_path = os.path.join(save_dir, f"{year}_{cat}_etaphi.png")
            fig.savefig(save_path)
            plt.close(fig)
            plt.show()


if __name__ == "__main__":
    main()
