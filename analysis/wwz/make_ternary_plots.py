import mpltern
import numpy as np
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.lines as mlines

# Avoid import but not used error for mpltern:
#     - We need to import it in order to register the "ternary" projection
#     - But mpltern is not explicitly used anywhere, making flake8 upset
#     - Thus just add an assert line, so that flake8 will not catch it
assert mpltern

###############################################################################
### Plotting Values and Variables. Look here if you want to make some changes!!

# Number of WWZ Events
n_event = 30000

#Size of the figure
figsize_val = 12.0

#Size of the dots
s_val_mc = 0.6
s_val_data = 100

#Rotation of the ternary plot
rotation_val = -60

#Ticksize on the ternary plot
tick_size_val = 18

#Opacity for each process
alpha_wwz = 1.0
alpha_zh = 1.0
alpha_zz = 0.3
alpha_ttz = 0.3
alpha_twz = 0.3
alpha_wz = 0.3
alpha_bkg = 0.3
alpha_data = 1.0

#Colors for each process
wwz_color = 'red'
zh_color = 'blue'
zz_color = 'green'
ttz_color = 'green'
twz_color = 'green'
wz_color = 'green'
bkg_color = 'green'
data_color = 'black'

# Include Data on the ternary plots (data must be present in the txt file)
include_data = True

# Signal Region Bin Lines (False is no lines)
draw_lines = True

#Signal Region Line Width
sr_lw = 3

#CMS Label
CMS_x = -0.07
CMS_y = 0.99
CMS_size = 40

#Supplementary Label
on_off = True # True is on
Supp_x = 0.07
Supp_y = 0.99
Supp_size = 36

#Lumi Stuff
lumi_x = 0.81
lumi_y = 0.99
lumi_size = 35

#Legend Stuff
legend_marker_size = 30
label_font_size = 30
legend_x = 0.90
legend_y = 0.35

#Region and Run Label Stuff
label_size = 36
label_x = 0.15
label_y = 0.31

# Axis Labels
zh_label_x = 0.06
zh_label_y = 0.87
wwz_label_x = 0.99
wwz_label_y = 0.87
bkg_label_x = 0.515
bkg_label_y = 0.03
axis_label_size = 28

###############################################################################

# dictionary for run and flavor dependent values
run_flavor_dict = {
    "run2_of":{
        "zh_fact": 0.898,
        "zz_fact": 2.257,
        "ttz_fact": 0.45,
        "twz_fact": 0.204,
        "wz_fact": 0.221,
        "bkg_fact": 0.181,
        "sf_fact": 1.0,
        "line_1": [[0.92, 0.46], [0.00, 0.46], [0.08, 0.08]],
        "line_2": [[0.79, 0.395], [0.00, 0.395], [0.21, 0.21]],
        "line_3": [[0.59, 0.295], [0.00, 0.295], [0.41, 0.41]],
        "line_4": [[0.00, 0.495], [0.99, 0.495], [0.01, 0.01]],
        "line_5": [[0.00, 0.445], [0.89, 0.445], [0.11, 0.11]],
        "line_6": [[0.00, 0.31], [0.62, 0.31], [0.38, 0.38]],
    },
    "run2_sf":{
        "zh_fact": 0.906,
        "zz_fact": 18.96,
        "ttz_fact": 0.475,
        "twz_fact": 0.221,
        "wz_fact": 0.174,
        "bkg_fact": 0.934,
        "sf_fact": 0.15,
        "line_1": [[0.9, 0.45], [0.00, 0.45], [0.1, 0.1]],
        "line_2": [[0.85, 0.425], [0.00, 0.425], [0.15, 0.15]],
        "line_3": [[0.75, 0.375], [0.00, 0.375], [0.25, 0.25]],
        "line_4": [[0.00, 0.49], [0.98, 0.49], [0.02, 0.02]],
        "line_5": [[0.00, 0.45], [0.9, 0.45], [0.1, 0.1]],
        "line_6": [[0.00, 0.4], [0.8, 0.4], [0.2, 0.2]],
    },
    "run3_of":{
        "zh_fact": 0.942,
        "zz_fact": 2.658,
        "ttz_fact": 0.501,
        "twz_fact": 0.273,
        "wz_fact": 0.154,
        "bkg_fact": 0.165,
        "sf_fact": 1.0,
        "line_1": [[0.84, 0.42], [0.00, 0.42], [0.16, 0.16]],
        "line_2": [[0.00, 0.415], [0.83, 0.415], [0.17, 0.17]],
    },
    "run3_sf":{
        "zh_fact": 0.973,
        "zz_fact": 23.014,
        "ttz_fact": 0.482,
        "twz_fact": 0.277,
        "wz_fact": 0.086,
        "bkg_fact": 1.777,
        "sf_fact": 0.15,
        "line_1": [[0.96, 0.48], [0.00, 0.48], [0.04, 0.04]],
        "line_2": [[0.00, 0.48], [0.96, 0.48], [0.04, 0.04]],
    }
}

###############################################################################

# Scores needed when looping through the txt file

# Scores for the WWZ process
p_wwz_s_wwz = []
p_wwz_s_zh = []
p_wwz_s_bkg = []
# Scores for the ZH process
p_zh_s_wwz = []
p_zh_s_zh  = []
p_zh_s_bkg = []
# Scores for the ZZ process
p_zz_s_wwz = []
p_zz_s_zh  = []
p_zz_s_bkg = []
# Scores for the ttZ process
p_ttz_s_wwz = []
p_ttz_s_zh  = []
p_ttz_s_bkg = []
# Scores for the tWZ process
p_twz_s_wwz = []
p_twz_s_zh  = []
p_twz_s_bkg = []
# Scores for the WZ process
p_wz_s_wwz = []
p_wz_s_zh  = []
p_wz_s_bkg = []
# Scores for the Other process
p_bkg_s_wwz = []
p_bkg_s_zh = []
p_bkg_s_bkg = []
# Scores for the Other process
p_data_s_wwz = []
p_data_s_zh = []
p_data_s_bkg = []
# Event weights
p_wwz_weights = []
p_zh_weights = []
p_zz_weights = []
p_ttz_weights = []
p_twz_weights = []
p_wz_weights = []
p_bkg_weights = []
p_data_weights = []

###############################################################################

def scan_file(txt_file,key):

    n_wwz_events = n_event
    n_zh_events = n_event*run_flavor_dict[key]["zh_fact"]
    n_zz_events = n_event*run_flavor_dict[key]["zz_fact"]*run_flavor_dict[key]["sf_fact"]
    n_ttz_events = n_event*run_flavor_dict[key]["ttz_fact"]
    n_twz_events = n_event*run_flavor_dict[key]["twz_fact"]
    n_wz_events = n_event*run_flavor_dict[key]["wz_fact"]
    n_bkg_events = n_event*run_flavor_dict[key]["bkg_fact"]

    wwz_counter = 0
    zh_counter = 0
    zz_counter = 0
    ttz_counter = 0
    twz_counter = 0
    wz_counter = 0
    bkg_counter = 0

    for x in range(len(txt_file)):
        wwz_score = txt_file[x][1]
        zh_score = txt_file[x][2]
        bkg_score = txt_file[x][3]

        if txt_file[x][0] == 0:  # Process = WWZ
            wwz_counter += 1
            if wwz_counter > n_wwz_events:
                continue
            p_wwz_s_wwz.append(wwz_score)
            p_wwz_s_zh.append(zh_score)
            p_wwz_s_bkg.append(bkg_score)
            p_wwz_weights.append(txt_file[x][4])
        if txt_file[x][0] == -1:  # Process = data
            p_data_s_wwz.append(wwz_score)
            p_data_s_zh.append(zh_score)
            p_data_s_bkg.append(bkg_score)
            p_data_weights.append(txt_file[x][4])
        if txt_file[x][0] == 1:  # Process = ZH
            zh_counter += 1
            if zh_counter > n_zh_events:
                continue
            p_zh_s_wwz.append(wwz_score)
            p_zh_s_zh.append(zh_score)
            p_zh_s_bkg.append(bkg_score)
            p_zh_weights.append(txt_file[x][4])
        if txt_file[x][0] == 2:  # Process = ZZ
            zz_counter += 1
            if zz_counter > n_zz_events:
                continue
            p_zz_s_wwz.append(wwz_score)
            p_zz_s_zh.append(zh_score)
            p_zz_s_bkg.append(bkg_score)
            p_zz_weights.append(txt_file[x][4])
        if txt_file[x][0] == 3:  # Process = ttZ
            ttz_counter += 1
            if ttz_counter > n_ttz_events:
                continue
            p_ttz_s_wwz.append(wwz_score)
            p_ttz_s_zh.append(zh_score)
            p_ttz_s_bkg.append(bkg_score)
            p_ttz_weights.append(txt_file[x][4])
        if txt_file[x][0] == 4:  # Process = twZ
            twz_counter += 1
            if twz_counter > n_twz_events:
                continue
            p_twz_s_wwz.append(wwz_score)
            p_twz_s_zh.append(zh_score)
            p_twz_s_bkg.append(bkg_score)
            p_twz_weights.append(txt_file[x][4])
        if txt_file[x][0] == 5:  # Process = WZ
            wz_counter += 1
            if wz_counter > n_wz_events:
                continue
            p_wz_s_wwz.append(wwz_score)
            p_wz_s_zh.append(zh_score)
            p_wz_s_bkg.append(bkg_score)
            p_wz_weights.append(txt_file[x][4])
        if txt_file[x][0] == 6:  # Process = BKG
            bkg_counter += 1
            if bkg_counter > n_bkg_events:
                continue
            p_bkg_s_wwz.append(wwz_score)
            p_bkg_s_zh.append(zh_score)
            p_bkg_s_bkg.append(bkg_score)
            p_bkg_weights.append(txt_file[x][4])

###############################################################################

def make_plot(run, flavor, key):
    # Lists need to be arrays for mpltern
    arr_p_wwz_s_wwz = np.array(p_wwz_s_wwz)
    arr_p_wwz_s_zh = np.array(p_wwz_s_zh)
    arr_p_wwz_s_bkg = np.array(p_wwz_s_bkg)
    arr_p_zh_s_wwz = np.array(p_zh_s_wwz)
    arr_p_zh_s_zh = np.array(p_zh_s_zh)
    arr_p_zh_s_bkg = np.array(p_zh_s_bkg)
    arr_p_zz_s_wwz = np.array(p_zz_s_wwz)
    arr_p_zz_s_zh = np.array(p_zz_s_zh)
    arr_p_zz_s_bkg = np.array(p_zz_s_bkg)
    arr_p_ttz_s_wwz = np.array(p_ttz_s_wwz)
    arr_p_ttz_s_zh = np.array(p_ttz_s_zh)
    arr_p_ttz_s_bkg = np.array(p_ttz_s_bkg)
    arr_p_twz_s_wwz = np.array(p_twz_s_wwz)
    arr_p_twz_s_zh = np.array(p_twz_s_zh)
    arr_p_twz_s_bkg = np.array(p_twz_s_bkg)
    arr_p_wz_s_wwz = np.array(p_wz_s_wwz)
    arr_p_wz_s_zh = np.array(p_wz_s_zh)
    arr_p_wz_s_bkg = np.array(p_wz_s_bkg)
    arr_p_bkg_s_wwz = np.array(p_bkg_s_wwz)
    arr_p_bkg_s_zh = np.array(p_bkg_s_zh)
    arr_p_bkg_s_bkg = np.array(p_bkg_s_bkg)
    arr_p_data_s_wwz = np.array(p_data_s_wwz)
    arr_p_data_s_zh = np.array(p_data_s_zh)
    arr_p_data_s_bkg = np.array(p_data_s_bkg)
    arr_p_wwz_weights = np.array(p_wwz_weights)
    arr_p_zh_weights = np.array(p_zh_weights)
    arr_p_zz_weights = np.array(p_zz_weights)
    arr_p_ttz_weights = np.array(p_ttz_weights)
    arr_p_twz_weights = np.array(p_twz_weights)
    arr_p_wz_weights = np.array(p_wz_weights)
    arr_p_bkg_weights = np.array(p_bkg_weights)
    arr_p_data_weights = np.array(p_data_weights)

    # Grab the dict
    rel_dict = run_flavor_dict[key]

    fig = plt.figure(figsize=(figsize_val,figsize_val))

    #Create the ternary plot
    ax = fig.add_subplot(111, projection="ternary", rotation=rotation_val)
    ax.tick_params(labelsize=tick_size_val)
    ax.scatter(arr_p_wwz_s_wwz,arr_p_wwz_s_zh,arr_p_wwz_s_bkg, s=s_val_mc, color = wwz_color, alpha = alpha_wwz)
    ax.scatter(arr_p_zh_s_wwz,arr_p_zh_s_zh,arr_p_zh_s_bkg, s=s_val_mc, color = zh_color, alpha = alpha_zh)
    ax.scatter(arr_p_zz_s_wwz,arr_p_zz_s_zh,arr_p_zz_s_bkg, s=s_val_mc, color = zz_color, alpha = alpha_zz)
    ax.scatter(arr_p_ttz_s_wwz,arr_p_ttz_s_zh,arr_p_ttz_s_bkg, s=s_val_mc, color = ttz_color, alpha = alpha_ttz)
    ax.scatter(arr_p_twz_s_wwz,arr_p_twz_s_zh,arr_p_twz_s_bkg, s=s_val_mc, color = twz_color, alpha = alpha_twz)
    ax.scatter(arr_p_wz_s_wwz,arr_p_wz_s_zh,arr_p_wz_s_bkg, s=s_val_mc, color = wz_color, alpha = alpha_wz)
    ax.scatter(arr_p_bkg_s_wwz,arr_p_bkg_s_zh,arr_p_bkg_s_bkg, s=s_val_mc, color = bkg_color, alpha = alpha_bkg)
    if include_data:
        ax.scatter(arr_p_data_s_wwz,arr_p_data_s_zh,arr_p_data_s_bkg, s=s_val_data, color = data_color, alpha = alpha_data)

    if draw_lines:
        # Draw the SR lines
        ax.axline([0.5, 0.5, 0.0], [0.0, 0.0, 1.0], color='black', linewidth = sr_lw)
        if (run == "run2"):
            ax.plot(rel_dict['line_1'][0],rel_dict['line_1'][1],rel_dict['line_1'][2],color='black',linewidth=sr_lw)
            ax.plot(rel_dict['line_2'][0],rel_dict['line_2'][1],rel_dict['line_2'][2],color='black',linewidth=sr_lw)
            ax.plot(rel_dict['line_3'][0],rel_dict['line_3'][1],rel_dict['line_3'][2],color='black',linewidth=sr_lw)
            ax.plot(rel_dict['line_4'][0],rel_dict['line_4'][1],rel_dict['line_4'][2],color='black',linewidth=sr_lw)
            ax.plot(rel_dict['line_5'][0],rel_dict['line_5'][1],rel_dict['line_5'][2],color='black',linewidth=sr_lw)
            ax.plot(rel_dict['line_6'][0],rel_dict['line_6'][1],rel_dict['line_6'][2],color='black',linewidth=sr_lw)
        elif (run == "run3"):
            ax.plot(rel_dict['line_1'][0],rel_dict['line_1'][1],rel_dict['line_1'][2],color='black',linewidth=sr_lw)
            ax.plot(rel_dict['line_2'][0],rel_dict['line_2'][1],rel_dict['line_2'][2],color='black',linewidth=sr_lw)
        else:
            raise Exception("Unrecognized run!")


    # CMS Stuff
    fig.text(CMS_x, CMS_y, "CMS", fontsize=CMS_size, weight="bold", transform=fig.transFigure)
    if on_off:
        fig.text(Supp_x, Supp_y, "$\\it{Supplementary}$", fontsize=Supp_size, transform=fig.transFigure)
    if (run == "run2"):
        fig.text(lumi_x, lumi_y, "138 $\\mathrm{fb^{{-}1}}$ (13 TeV)", fontsize=lumi_size, transform=fig.transFigure)
    elif (run == "run3"):
        fig.text(lumi_x, lumi_y, "62 $\\mathrm{fb^{{-}1}}$ (13.6 TeV)", fontsize=lumi_size, transform=fig.transFigure)

    # Add custom "Opposite/Same Flavor" text
    if (flavor == "of"):
        fig.text(label_x, label_y, "Opposite Flavor \nSignal Region", fontsize=label_size, ha='center', transform=fig.transFigure, weight="bold")
    elif (flavor == "sf"):
        fig.text(label_x, label_y, "Same Flavor \nSignal Region", fontsize=label_size, ha='center', transform=fig.transFigure, weight="bold")

    # Create legend
    legend_dot_data = mlines.Line2D([], [], marker='o', color=data_color, markersize=legend_marker_size, linestyle='None', label='data')
    legend_dot_wwz = mlines.Line2D([], [], marker='o', color=wwz_color, markersize=legend_marker_size, linestyle='None', label='WWZ')
    legend_dot_zh = mlines.Line2D([], [], marker='o', color=zh_color, markersize=legend_marker_size, linestyle='None', label='ZH')
    legend_dot_bkg = mlines.Line2D([], [], marker='o', color=bkg_color, markersize=legend_marker_size, linestyle='None', label='Backgrounds')

    if include_data:
        fig.legend(handles=[legend_dot_data, legend_dot_wwz, legend_dot_zh, legend_dot_bkg], fontsize=label_font_size, bbox_to_anchor=(legend_x, legend_y), loc='center', frameon=False)
    else:
        fig.legend(handles=[legend_dot_wwz, legend_dot_zh, legend_dot_bkg], fontsize=label_font_size, bbox_to_anchor=(legend_x, legend_y), loc='center', frameon=False)

    # Manually put axis labels
    fig.text(zh_label_x, zh_label_y, "ZH BDT \nScore", fontsize=axis_label_size, ha='center', transform=fig.transFigure, weight="normal")
    fig.text(wwz_label_x, wwz_label_y, "WWZ BDT \nScore", fontsize=axis_label_size, ha='center', transform=fig.transFigure, weight="normal")
    fig.text(bkg_label_x, bkg_label_y, "BKG BDT \nScore", fontsize=axis_label_size, ha='center', transform=fig.transFigure, weight="normal")

    # Save the figure as a PNG
    if draw_lines:
        fig.savefig(f'{key}_bins.png', dpi=150, bbox_inches='tight')
    else:
        fig.savefig(f'{key}_nobins.png', dpi=150, bbox_inches='tight')



###############################################################################

def make_key(run, flavor):
    if run not in ["run2","run3"]:
        raise Exception("Run is not recognized")
    if flavor not in ["sf","of"]:
        raise Exception("Flavor not recognized")

    key = (run + "_" + flavor)
    return key

###############################################################################


def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("txt_file_path", help = "The path to the txt file")
    parser.add_argument('-f', "--flavor", default='of', help = "Which flavor to grab", choices=["of","sf"])
    parser.add_argument('-r', "--run", default='run2', help = "Which Run is this", choices=["run2","run3"])
    args = parser.parse_args()

    # Define the arguments
    input_file = args.txt_file_path
    flavor = args.flavor
    run = args.run

    key = make_key(run, flavor)

    # Read the txt file containing the eventlist and scores
    txt_file = pd.read_csv(input_file, sep=' ').values

    # Scan over the txt file
    scan_file(txt_file, key)

    # Make the ternary plot
    make_plot(run, flavor, key)

if __name__ == "__main__":
    main()
