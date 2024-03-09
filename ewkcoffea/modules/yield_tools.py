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
                    val = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":syst_name }].values(flow=True)))
                    var = sum(sum(histo[{"category":cat_name,"process":sample_dict[proc_name],"systematic":syst_name }].variances(flow=True)))
                    yld_dict[cat_name][syst_name][proc_name] = [val,var]

    return yld_dict


############## Manipulating yields in a yield dict ##############

# Calculate data-driven background estimation from relevant CRs
def do_tf(yld_mc,yld_data,kappas,tf_map,quiet=True):

    yld_mc_out = copy.deepcopy(yld_mc)
    kappas_out = copy.deepcopy(kappas)
    gmn_alpha_out = {}

    # Loop over cat and do NSF calculation for each relevant proc in each cat
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
                valvar_bkg = valvar_op(valvar_nsf, yld_mc[cat]["nominal"][proc_of_interest], "prod")

                # Dump NSF value
                if not quiet:
                    relerr = 100*(np.sqrt(valvar_nsf[1])/valvar_nsf[0])
                    print(f"NSF for {cat} {proc_of_interest}: {np.round(valvar_nsf[0],2)} +- {np.round(relerr,2)}%")
                    #print(f"NSF for {cat} {proc_of_interest}: {np.round(valvar_nsf[0],2)} +- {np.round(np.sqrt(valvar_nsf[1]),2)}")

                # Put the old yield times the NSF into the out dict
                yld_mc_out[cat]["nominal"][proc_of_interest] = valvar_bkg

                # Get the gmN numbers
                valvar_alpha = valvar_op(valvar_bkg,valvar_cr_data, "div")
                if cat not in gmn_alpha_out: gmn_alpha_out[cat] = {}
                if cr_name not in gmn_alpha_out[cat]: gmn_alpha_out[cat][cr_name] = {"N": valvar_cr_data[0], "proc_alpha": {}}
                gmn_alpha_out[cat][cr_name]["proc_alpha"][proc_of_interest] = valvar_alpha[0]

                ### Handle the kappas ###

                if kappas is not None:

                    # Loop over syst and replace the kappas with the e.g. MC_SR_up/MC_CR_up
                    for syst_base_name in kappas[cat]:

                        # Skip the stats uncertainties (they are not correlated between CR and SR, so just leave alone)
                        if syst_base_name.startswith("stats_"): continue

                        sr_up = kappas[cat][syst_base_name][proc_of_interest]["Up"]
                        sr_do = kappas[cat][syst_base_name][proc_of_interest]["Down"]
                        cr_up = kappas[cr_name][syst_base_name][proc_of_interest]["Up"]
                        cr_do = kappas[cr_name][syst_base_name][proc_of_interest]["Down"]
                        new_kappa_up = valvar_op(sr_up,cr_up,"div")
                        new_kappa_do = valvar_op(sr_do,cr_do,"div")

                        kappas_out[cat][syst_base_name][proc_of_interest]["Up"] = new_kappa_up
                        kappas_out[cat][syst_base_name][proc_of_interest]["Down"] = new_kappa_do

    return [yld_mc_out, kappas_out, gmn_alpha_out]
