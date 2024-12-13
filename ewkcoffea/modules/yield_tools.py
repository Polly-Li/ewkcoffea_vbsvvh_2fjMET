import numpy as np
import copy

# Yields for getting yields from histograms, and manipulating those yields

########### General ###########

# Takes two pairs [val1,var1] and [val2,var2], returns the product or sum with error propagated
# Note the var is like sumw2 i.e. alreayd squared (in both input and output)
def valvar_op(valvar_1, valvar_2, op):
    val1 = valvar_1[0]
    var1 = valvar_1[1]
    val2 = valvar_2[0]
    var2 = valvar_2[1]

    if op == "prod":
        val = val1*val2
        var = (val**2) * ( (np.sqrt(var1)/val1)**2 + (np.sqrt(var2)/val2)**2 )
    elif op == "div":
        val = val1/val2
        var = (val**2) * ( (np.sqrt(var1)/val1)**2 + (np.sqrt(var2)/val2)**2 )
    elif op == "sum":
        val = val1 + val2
        var = var1 + var2
    elif op == "diff":
        val = val1 - val2
        var = var1 + var2
    else:
        raise Exception("Unknown operatation")

    return [val,var]


############## Getting yields from a histo ##############

# Get the yields (nested in the order: year,cat,syst,proc)
def get_yields(histo,sample_dict,blind=True,systematic_name=None):

    yld_dict = {}

    if systematic_name is None: syst_lst = histo.axes["systematic"]
    else: syst_lst = [systematic_name]

    # Look at the yields in the histo
    for cat_name in histo.axes["category"]:
        yld_dict[cat_name] = {}
        for syst_name in syst_lst:
            yld_dict[cat_name][syst_name ] = {}
            for proc_name in sample_dict.keys():
                if blind and (("data" in proc_name) and (not cat_name.startswith("cr_"))):
                    # If this is data and we're not in a CR category, put placeholder numbers for now
                    yld_dict[cat_name][syst_name][proc_name] = [-999,-999]
                else:
                    val_n = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":"nominal" }].values(flow=True)))
                    var_n = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":"nominal" }].variances(flow=True)))

                    val = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":syst_name }].values(flow=True)))
                    var = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":syst_name }].variances(flow=True)))

                    # Zero out bins that are consistent with zero in the nominal
                    if (np.sqrt(var_n) >= val_n):
                        val = 0
                        var = 0

                    yld_dict[cat_name][syst_name][proc_name] = [val,var]

    return yld_dict


############## Manipulating yields in a yield dict ##############

# Scale processes in yield dictionary by a given factor
# Will scale all systematic variations (not just nominal) by the factor
# Returns a new dictionary, does not modify original
def scale_yld_dict(in_yld_dict,scaling_dict):
    out_dict = copy.deepcopy(in_yld_dict)
    for cat in in_yld_dict:
        for syst in in_yld_dict[cat]:
            for proc in in_yld_dict[cat][syst]:
                if proc in scaling_dict:
                    out_dict[cat][syst][proc][0] = in_yld_dict[cat][syst][proc][0]*scaling_dict[proc]
                    out_dict[cat][syst][proc][1] = in_yld_dict[cat][syst][proc][1]*scaling_dict[proc]
    return out_dict


# Calculate data-driven background estimation from relevant CRs
def do_tf(yld_mc,yld_data,kappas,tf_map,quiet=True):

    yld_mc_out = copy.deepcopy(yld_mc)
    kappas_out = copy.deepcopy(kappas)
    gmn_alpha_out = {}

    nsf_dict = get_nsf_dict(yld_mc,yld_data,tf_map)

    # Loop over cat and do NSF calculation for each relevant proc in each cat
    for cat in nsf_dict:
        for proc in nsf_dict[cat]:

            # Get the NSF = N_CR_with_other_bkg_subtracted / MC_CR, the scaled yld, and the alpha
            valvar_nsf    = nsf_dict[cat][proc]["nsf"]
            valvar_bkg    = nsf_dict[cat][proc]["bkg_scaled"]
            cr_name       = nsf_dict[cat][proc]["cr_name"]
            valvar_alpha  = nsf_dict[cat][proc]["alpha"]
            valvar_cr_obs = nsf_dict[cat][proc]["cr_obs"]

            # Dump NSF value
            if not quiet:
                relerr = 100*(np.sqrt(valvar_nsf[1])/valvar_nsf[0])
                print(f"NSF for {cat} {proc}: {np.round(valvar_nsf[0],2)} +- {np.round(relerr,2)}%")
                #print(f"NSF for {cat} {proc}: {np.round(valvar_nsf[0],2)} +- {np.round(np.sqrt(valvar_nsf[1]),2)}")

            # Put the old yield times the NSF into the out dict
            yld_mc_out[cat]["nominal"][proc] = valvar_bkg

            # Get the gmN numbers
            if cat not in gmn_alpha_out: gmn_alpha_out[cat] = {}
            if cr_name not in gmn_alpha_out[cat]: gmn_alpha_out[cat][cr_name] = {"N": valvar_cr_obs[0], "proc_alpha": {}}
            gmn_alpha_out[cat][cr_name]["proc_alpha"][proc] = valvar_alpha[0]

            ### Handle the kappas ###

            if kappas is not None:

                # Loop over syst and replace the kappas with the e.g. MC_SR_up/MC_CR_up
                for syst_base_name in kappas[cat]:

                    # Skip the stats uncertainties (they are not correlated between CR and SR, so just leave alone)
                    if syst_base_name.startswith("stats_"): continue

                    sr_up = kappas[cat][syst_base_name][proc]["Up"]
                    sr_do = kappas[cat][syst_base_name][proc]["Down"]
                    cr_up = kappas[cr_name][syst_base_name][proc]["Up"]
                    cr_do = kappas[cr_name][syst_base_name][proc]["Down"]
                    new_kappa_up = valvar_op(sr_up,cr_up,"div")
                    new_kappa_do = valvar_op(sr_do,cr_do,"div")

                    kappas_out[cat][syst_base_name][proc]["Up"] = new_kappa_up
                    kappas_out[cat][syst_base_name][proc]["Down"] = new_kappa_do

    return [yld_mc_out, kappas_out, gmn_alpha_out]


# Calculate the NSF for the cats and procs in the given background tf map, return as a dict
def get_nsf_dict(yld_mc, yld_data, tf_map):

    out_dict = {}

    for cat in yld_mc:

        # Just nonimal for now
        for proc_of_interest in yld_mc[cat]["nominal"]:

            # Skip procs that do not get TF calculations
            if proc_of_interest not in tf_map: continue
            elif cat not in tf_map[proc_of_interest]:
                #print(f"Warning, cat \"{cat}\" not defined for this proc.")
                continue

            # Otherwise we go ahead and do the background estimation stuff
            else:

                if cat not in out_dict: out_dict[cat] = {}

                # Get the nominal mc and data yields in the CR
                cr_name = tf_map[proc_of_interest][cat]
                valvar_cr_mc   = yld_mc[cr_name]["nominal"][proc_of_interest]
                valvar_cr_data = yld_data[cr_name]["nominal"]["data"]

                # Sum up all contributions in the CR besides bkg of interest
                valvar_bkg_all_but_bkg_of_interest = [0,0]
                for p in yld_mc[cr_name]["nominal"]:
                    if p != proc_of_interest:
                        valvar_bkg_all_but_bkg_of_interest[0] += yld_mc[cr_name]["nominal"][p][0]
                        valvar_bkg_all_but_bkg_of_interest[1] += yld_mc[cr_name]["nominal"][p][1]

                # Subtract those extra background contributions from the data
                valvar_cr_data_corrected = valvar_op(valvar_cr_data, valvar_bkg_all_but_bkg_of_interest, "diff")

                # Calculate the NSF = N_CR_with_other_bkg_subtracted / MC_CR, use this to scale the yld
                valvar_nsf = valvar_op(valvar_cr_data_corrected, valvar_cr_mc, "div")

                # Scale the yield by this NSF
                valvar_bkg = valvar_op(valvar_nsf, yld_mc[cat]["nominal"][proc_of_interest], "prod")

                # Get the alpha value (needed for datacard)
                valvar_alpha = valvar_op(valvar_bkg,valvar_cr_data, "div") ###

                # Put these outputs (nsf, the scaled bkg, and the alpha) into the out dict
                out_dict[cat][proc_of_interest] = {
                    "nsf"        : valvar_nsf,
                    "bkg_scaled" : valvar_bkg,
                    "cr_name"    : cr_name,
                    "alpha"      : valvar_alpha,
                    "cr_obs"     : valvar_cr_data,
                }

    return out_dict
