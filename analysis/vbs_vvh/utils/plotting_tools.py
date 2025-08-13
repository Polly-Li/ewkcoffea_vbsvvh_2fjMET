import copy
import hist
import numpy as np

# Get the list of categories on the sparese axis
def get_axis_cats(histo,axis_name):
    process_list = [x for x in histo.axes[axis_name]]
    return process_list


# Rebin according to https://github.com/CoffeaTeam/coffea/discussions/705
def rebin(histo,factor):
    return histo[..., ::hist.rebin(factor)]

def rebin_to_nonnegative(histo, category_axis_to_sum):
    """
    Rebins a 1D hist.Hist such that all bins are non-negative in the summed total histogram.
    - Assumes histo is already sliced to a single dense axis (process: bkg, or process:sig)
    - Sums over the provided categorical axis (e.g., "process") before determining the rebin factor.
    
    Returns: rebinned hist. (or original histo if total sum <0, which is to be avoided during cutflow optimization)
    """
    # Sum over category axis
    total = histo[{category_axis_to_sum: sum}]

    values, edges = total.to_numpy()

    if np.sum(values) < 0:
        print("[rebin_to_nonnegative] Total histogram sum < 0, cannot fix negatives by rebinning.")
        return histo,1

    n_bins = len(values)
    factor = 1

    for fac in range(1, n_bins + 1): #leave it untouch if no negative bins from the start
        if n_bins % fac != 0:
            continue
        rebinned_vals = values.reshape(-1, fac).sum(axis=1)
        if np.all(rebinned_vals >= 0):
            factor = fac
            break

    if factor >1:
        print(f"[rebin_to_nonnegative] Rebinning axis with factor {factor} to eliminate negatives.")
    
    return histo[..., ::hist.rebin(factor)],factor

# Scale the hist, see https://github.com/CoffeaTeam/coffea/discussions/705
# Seems to modify in place
def scale(histo,axis_name,scale_dict):
    for i, name in enumerate(histo.axes[axis_name]):
        histo.view(flow=True)[i] *= scale_dict.get(name,1)
    return histo


# Regroup categories (e.g. processes)
def group(h, oldname, newname, grouping):

    # Build up a grouping dict that drops any proc that is not in our h
    grouping_slim = {}
    proc_lst = get_axis_cats(h,oldname)
    for grouping_name in grouping.keys():
        for proc in grouping[grouping_name]:
            if proc in proc_lst:
                if grouping_name not in grouping_slim:
                    grouping_slim[grouping_name] = []
                grouping_slim[grouping_name].append(proc)
            #else:
            #    print(f"WARNING: process {proc} not in this hist")

    # From Nick: https://github.com/CoffeaTeam/coffea/discussions/705#discussioncomment-4604211
    hnew = hist.Hist(
        hist.axis.StrCategory(grouping_slim, name=newname),
        *(ax for ax in h.axes if ax.name != oldname),
        storage=h.storage_type(),
    )
    for i, indices in enumerate(grouping_slim.values()):
        hnew.view(flow=True)[i] = h[{oldname: indices}][{oldname: sum}].view(flow=True)

    return hnew


# Merges the last bin (overflow) into the second to last bin, zeros the content of the last bin, returns a new hist
# Note assumes just one axis!
def merge_overflow(hin):
    hout = copy.deepcopy(hin)
    for cat_idx,arr in enumerate(hout.values(flow=True)):
        hout.values(flow=True)[cat_idx][-2] += hout.values(flow=True)[cat_idx][-1]
        hout.values(flow=True)[cat_idx][-1] = 0
        hout.variances(flow=True)[cat_idx][-2] += hout.variances(flow=True)[cat_idx][-1]
        hout.variances(flow=True)[cat_idx][-1] = 0
    return hout
