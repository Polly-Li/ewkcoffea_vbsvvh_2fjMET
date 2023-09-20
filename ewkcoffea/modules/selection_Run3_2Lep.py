import numpy as np
import awkward as ak
from mt2 import mt2
from coffea.nanoevents.methods import vector

import topcoffea.modules.event_selection as tc_es

from ewkcoffea.modules.paths import ewkcoffea_path
from topcoffea.modules.get_param_from_jsons import GetParam
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))

# The datasets we are using, and the triggers in them
dataset_dict = {


    "2022" : {
        "EGamma" : [
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        ],
        "Muon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        ],
        "MuonEG" : [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
    }

}

trgs_for_matching = {

    "2022" : {
        "m_m" : {
            "trg_lst" : ["Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"],
            "offline_thresholds" : [10.0,10.0],
        },
        "e_e" : {
            "trg_lst" : ["Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"],
            "offline_thresholds" : [10.0,10.0],
        },
        "m_e" : {
            "trg_lst" : ["Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"],
            "offline_thresholds" : [10.0,10.0],
        },
        "e_m" : {
            "trg_lst" : ["Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
            "offline_thresholds" : [10.0,10.0],
        },
    }
}


# Hard coded dictionary for figuring out overlap...
#   - No unique way to do this
#   - Note: In order for this to work properly, you should be processing all of the datastes to be used in the analysis
#   - Otherwise, you may be removing events that show up in other datasets you're not using
exclude_dict = {
    "2022": {
        "Muon"     : [],
        "EGamma"         : dataset_dict["2022"]["Muon"],
        "MuonEG"         : dataset_dict["2022"]["Muon"] + dataset_dict["2022"]["EGamma"],
    },
}


#------------------------------------------------------------------------------------------------------------------------------
# Apply trigger matching requirements to make sure pt is above online thresholds
def trg_matching(events,year):

    # Initialize return array to be True array with same shape as events
    ret_arr = ak.zeros_like(np.array(events.event), dtype=bool)

    # Get the leptons, sort and pad
    el = events.l_Run3_2Lep_t[abs(events.l_Run3_2Lep_t.pdgId)==11]
    el = ak.pad_none(el[ak.argsort(el.pt,axis=-1,ascending=False)],2)
    mu = events.l_Run3_2Lep_t[abs(events.l_Run3_2Lep_t.pdgId)==13]
    mu = ak.pad_none(mu[ak.argsort(mu.pt,axis=-1,ascending=False)],2)

    # Loop over offline cuts, make sure triggers pass the offline cuts for the associated triggers
    for l_l in trgs_for_matching[year]:


        # Check if lep pt passes the offline cuts
        offline_thresholds = trgs_for_matching[year][l_l]["offline_thresholds"]
        if   l_l == "m_m": offline_cut = ak.fill_none(((mu[:,0].pt > offline_thresholds[0]) & (mu[:,1].pt > offline_thresholds[1])),False)
        elif l_l == "e_e": offline_cut = ak.fill_none(((el[:,0].pt > offline_thresholds[0]) & (el[:,1].pt > offline_thresholds[1])),False)
        elif l_l == "m_e": offline_cut = ak.fill_none(((mu[:,0].pt > offline_thresholds[0]) & (el[:,0].pt > offline_thresholds[1]) & (mu[:,0].pt > el[:,0].pt )),False)
        elif l_l == "e_m": offline_cut = ak.fill_none(((el[:,0].pt > offline_thresholds[0]) & (mu[:,0].pt > offline_thresholds[1]) & (mu[:,0].pt < el[:,0].pt )),False)
        else: raise Exception("Unknown offline cut.")

        # Check if trigger passes the associated triggers
        trg_lst = trgs_for_matching[year][l_l]["trg_lst"]
        trg_passes = tc_es.passes_trg_inlst(events,trg_lst)

        # Build the return mask
        # The return mask started from an array of False
        # The way an event becomes True is if it passes a trigger AND passes the offline pt cuts associated with that trg
        false_arr = ak.zeros_like(np.array(events.event), dtype=bool) # False array with same shape as events
        ret_arr = ret_arr | ak.where(trg_passes,offline_cut,false_arr)

    return ret_arr

#------------------------------------------------------------------------------------------------------------------------------
# 2Lep Selection
def add2lmask_Run3_2Lep(events, year, isData):

    # Leptons and padded leptons
    leps = events.l_Run3_2Lep_t

    # Filters
    filter_flags = events.Flag
    filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & filter_flags.ecalBadCalibFilter 

    # Lep multiplicity
    nlep_2 = (ak.num(leps) == 2)

    mask = filters & nlep_2
    events['is2l'] = ak.fill_none(mask,False)


#------------------------------------------------------------------------------------------------------------------------------
# Jet Selection
def addjetmask_Run3_2Lep(events, year, isData):

    jets = events.jets_Run3_2Lep

    njet_2 = (ak.num(jets) >= 2)

    mask = njet_2
    events['has2jets'] = ak.fill_none(mask,False)
#------------------------------------------------------------------------------------------------------------------------------
# Met Mask
def addmetmask_Run3_2Lep(events, year, isData):

    met = events.MET

    met_30 = (met.pt <= 30)

    mask = met_30
    events['metmask'] = ak.fill_none(mask,False)
#------------------------------------------------------------------------------------------------------------------------------
# Do Run3 2Lep pre selection, construct event level mask
# Convenience function around get_Run3_2Lep_candidates() and get_z_candidate_mask()
def attach_Run3_2Lep_preselection_mask(events,lep_collection):

    # Pt requirements (assumes lep_collection is pt sorted and padded)
    pt_mask = ak.fill_none((lep_collection[:,0].pt > 25.0) & (lep_collection[:,1].pt > 15.0),False)

    # SFOS and OPOS masks
    os_mask = ak.any((((lep_collection[:,0:1].pdgId)*(lep_collection[:,1:2].pdgId)) < 0),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    sf_mask = ak.any((((lep_collection[:,0:1].pdgId) + (lep_collection[:,1:2].pdgId)) == 0),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    of_mask = ak.any((abs(lep_collection[:,0:1].pdgId) != abs(lep_collection[:,1:2].pdgId)),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    ee_mask = ak.any((abs(lep_collection[:,0:1].pdgId) == 11),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    mumu_mask = ak.any((abs(lep_collection[:,0:1].pdgId) == 13),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]

    os_mask = ak.fill_none(os_mask,False) # Replace the None with False in the mask just to make it easier to think about
    sf_mask = ak.fill_none(sf_mask,False) # Replace the None with False in the mask just to make it easier to think about
    of_mask = ak.fill_none(of_mask,False) # Replace the None with False in the mask just to make it easier to think about
    ee_mask = ak.fill_none(ee_mask,False) # Replace the None with False in the mask just to make it easier to think about
    mumu_mask = ak.fill_none(mumu_mask,False) # Replace the None with False in the mask just to make it easier to think about

    # mLL masks
    z_mass_mask = ak.any((abs(((lep_collection[:,0:1] + lep_collection[:,1:2]).mass) - 91.2) < 10.0),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    z_mass_mask = ak.fill_none(z_mass_mask,False) # Replace the None with False in the mask just to make it easier to think about

    of_mass_mask = ak.any((((lep_collection[:,0:1] + lep_collection[:,1:2]).mass) > 20.0),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    of_mass_mask = ak.fill_none(of_mass_mask,False) # Replace the None with False in the mask just to make it easier to think about

    # The final preselection mask
    Run3_2Lep_presel_mask = (os_mask & pt_mask)

    # Attach to the lepton objects
    events["Run3_2Lep_presel_sf_ee"] = (Run3_2Lep_presel_mask & sf_mask & ee_mask & z_mass_mask)
    events["Run3_2Lep_presel_sf_mumu"] = (Run3_2Lep_presel_mask & sf_mask & mumu_mask & z_mass_mask)
    events["Run3_2Lep_presel_of"] = (Run3_2Lep_presel_mask & of_mask & of_mass_mask)


