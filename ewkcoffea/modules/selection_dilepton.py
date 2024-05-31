import numpy as np
import awkward as ak
import xgboost as xgb
from mt2 import mt2

from coffea.nanoevents.methods import vector

import topcoffea.modules.event_selection as tc_es

from ewkcoffea.modules.paths import ewkcoffea_path
from topcoffea.modules.get_param_from_jsons import GetParam
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))


# The datasets we are using, and the triggers in them
dataset_dict = {

    "2016" : {
        "DoubleMuon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
            "Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
            "Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
        ],
        "DoubleEG" : [
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        ],
        "MuonEG" : [
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
            "Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
    },

    "2017" : {
        "DoubleMuon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        ],
        "DoubleEG" : [
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        ],
        "MuonEG" : [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
    },

    "2018" : {
        "EGamma" : [
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        ],
        "DoubleMuon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        ],
        "MuonEG" : [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
    },

    "2022" : {
        "EGamma" : [
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Ele30_WPTight_Gsf",
            "Ele32_WPTight_Gsf",
            "Ele32_WPTight_Gsf_L1DoubleEG",
            "Ele35_WPTight_Gsf",
            "Ele115_CaloIdVT_GsfTrkIdT",
            "DoubleEle25_CaloIdL_MW",
        ],
        "Muon" : [
            "IsoMu24",
            "IsoMu27",
            "Mu50",
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        ],
        "SingleMuon" : [
            "IsoMu24",
            "IsoMu27",
            "Mu50",
        ],
        "DoubleMuon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        ],
        "MuonEG" : [
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
    }
}

trgs_for_matching = {

    "2016" : {
        "m_m" : {
            "trg_lst" : dataset_dict["2016"]["DoubleMuon"],
            "offline_thresholds" : [20.0,10.0],
        },
        "e_e" : {
            "trg_lst" : dataset_dict["2016"]["DoubleEG"],
            "offline_thresholds" : [25.0,15.0],
        },
        "m_e" : {
            "trg_lst" : ["Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ"],
            "offline_thresholds" : [25,10],
        },
        "e_m" : {
            "trg_lst" : ["Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
            "offline_thresholds" : [25.0,10.0],
        },
    },
    "2017" : {
        "m_m" : {
            "trg_lst" : dataset_dict["2017"]["DoubleMuon"],
            "offline_thresholds" : [20.0,10.0],
        },
        "e_e" : {
            "trg_lst" : dataset_dict["2017"]["DoubleEG"],
            "offline_thresholds" : [25.0,15.0],
        },
        "m_e" : {
            "trg_lst" : ["Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"],
            "offline_thresholds" : [25.0,15.0],
        },
        "e_m" : {
            "trg_lst" : ["Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
            "offline_thresholds" : [25.0,10.0],
        },
    },
    "2018" : {
        "m_m" : {
            "trg_lst" : dataset_dict["2018"]["DoubleMuon"],
            "offline_thresholds" : [20.0,10.0],
        },
        "e_e" : {
            "trg_lst" : dataset_dict["2018"]["EGamma"],
            "offline_thresholds" : [25.0,15.0],
        },
        "m_e" : {
            "trg_lst" : ["Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"],
            "offline_thresholds" : [25.0,15.0],
        },
        "e_m" : {
            "trg_lst" : ["Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"],
            "offline_thresholds" : [25.0,10.0],
        },
    }
}


# Hard coded dictionary for figuring out overlap...
#   - No unique way to do this
#   - Note: In order for this to work properly, you should be processing all of the datastes to be used in the analysis
#   - Otherwise, you may be removing events that show up in other datasets you're not using
# For Era C which has both, the events in (SingleMuon, DoubleMuon) and (Muon) are exclusive so we do not perform duplicate removal bet
# For Era C, SingleMuon and DoubleMuon fall in the run ranges of [355800,357399] while Muon falls in [356400,357400]
exclude_dict = {
    "2016": {
        "DoubleMuon"     : [],
        "DoubleEG"       : dataset_dict["2016"]["DoubleMuon"],
        "MuonEG"         : dataset_dict["2016"]["DoubleMuon"] + dataset_dict["2016"]["DoubleEG"],
    },
    "2017": {
        "DoubleMuon"     : [],
        "DoubleEG"       : dataset_dict["2017"]["DoubleMuon"],
        "MuonEG"         : dataset_dict["2017"]["DoubleMuon"] + dataset_dict["2017"]["DoubleEG"],
    },
    "2018": {
        "DoubleMuon"     : [],
        "EGamma"         : dataset_dict["2018"]["DoubleMuon"],
        "MuonEG"         : dataset_dict["2018"]["DoubleMuon"] + dataset_dict["2018"]["EGamma"],
    },
    "C": {
        "Muon"           : [],
        "SingleMuon"     : [],
        "DoubleMuon"     : dataset_dict["2022"]["SingleMuon"],
        "EGamma"         : dataset_dict["2022"]["Muon"] + dataset_dict["2022"]["DoubleMuon"] + dataset_dict["2022"]["SingleMuon"],
        "MuonEG"         : dataset_dict["2022"]["Muon"] + dataset_dict["2022"]["DoubleMuon"] + dataset_dict["2022"]["SingleMuon"] + dataset_dict["2022"]["EGamma"],
    },
    "D": {
        "Muon"     : [],
        "EGamma"         : dataset_dict["2022"]["Muon"],
        "MuonEG"         : dataset_dict["2022"]["Muon"] + dataset_dict["2022"]["EGamma"],
    },
    "E": {
        "Muon"     : [],
        "EGamma"         : dataset_dict["2022"]["Muon"],
        "MuonEG"         : dataset_dict["2022"]["Muon"] + dataset_dict["2022"]["EGamma"],
    },
    "F": {
        "Muon"     : [],
        "EGamma"         : dataset_dict["2022"]["Muon"],
        "MuonEG"         : dataset_dict["2022"]["Muon"] + dataset_dict["2022"]["EGamma"],
    },
    "G": {
        "Muon"     : [],
        "EGamma"         : dataset_dict["2022"]["Muon"],
        "MuonEG"         : dataset_dict["2022"]["Muon"] + dataset_dict["2022"]["EGamma"],
    },
}


# Apply trigger matching requirements to make sure pt is above online thresholds
def trg_matching(events,year):

    # The trigger for 2016 and 2016APV are the same
    if year == "2016APV": year = "2016"

    # Initialize return array to be True array with same shape as events
    ret_arr = ak.zeros_like(np.array(events.event), dtype=bool)

    # Get the leptons, sort and pad
    el = events.lep_loose[abs(events.lep_loose.pdgId)==11]
    el = ak.pad_none(el[ak.argsort(el.pt,axis=-1,ascending=False)],2)
    mu = events.lep_loose[abs(events.lep_loose.pdgId)==13]
    mu = ak.pad_none(mu[ak.argsort(mu.pt,axis=-1,ascending=False)],2)

    # Loop over offline cuts, make sure triggers pass the offline cuts for the associated triggers
    for l_l in trgs_for_matching[year]:

        # Check if lep pt passes the offline cuts
        offline_thresholds = trgs_for_matching[year][l_l]["offline_thresholds"]
        if   l_l == "m_m": offline_cut = ak.fill_none(((mu[:,0].pt > offline_thresholds[0]) & (mu[:,1].pt > offline_thresholds[1])),False)
        elif l_l == "e_e": offline_cut = ak.fill_none(((el[:,0].pt > offline_thresholds[0]) & (el[:,1].pt > offline_thresholds[1])),False)
        elif l_l == "m_e": offline_cut = ak.fill_none(((mu[:,0].pt > offline_thresholds[0]) & (el[:,0].pt > offline_thresholds[1])),False)
        elif l_l == "e_m": offline_cut = ak.fill_none(((el[:,0].pt > offline_thresholds[0]) & (mu[:,0].pt > offline_thresholds[1])),False)
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


# 2l selection
def add2lmask(events, year, isData, sample_name,is2022):

    # Leptons and padded leptons
    leps_tight = events.lep_tight
    leps_tight_padded = ak.pad_none(leps_tight,2)
    leps_loose = events.lep_loose
    leps_loose_padded = ak.pad_none(leps_loose,2)

    # Filters
    filter_flags = events.Flag
    if is2022:
        filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & filter_flags.ecalBadCalibFilter & filter_flags.BadPFMuonDzFilter & filter_flags.hfNoisyHitsFilter & filter_flags.eeBadScFilter
    else:
        filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & (((year == "2016")|(year == "2016APV")) | filter_flags.ecalBadCalibFilter) & (isData | filter_flags.eeBadScFilter)

    # Conditions for masks
    n_loose_lep_2 = (ak.num(leps_loose) == 2)
    n_tight_lep_2 = (ak.num(leps_tight) == 2)
    n_tight_lep_al1 = (ak.num(leps_tight) >= 1)
    n_tight_lep_1 = (ak.num(leps_tight) == 1)
    tight_leading = ak.fill_none(((leps_tight_padded[:,0].pt == leps_loose_padded[:,0].pt) & (leps_tight_padded[:,0].eta == leps_loose_padded[:,0].eta)),False)
    tight_subleading = ak.fill_none(((leps_tight_padded[:,0].pt == leps_loose_padded[:,1].pt) & (leps_tight_padded[:,0].eta == leps_loose_padded[:,1].eta)),False)
    tight_leading = ak.fill_none(tight_leading,False)
    tight_subleading = ak.fill_none(tight_subleading,False)

    # Define the masks
    mask_loose      = filters & n_loose_lep_2
    mask_2_tight    = filters & n_loose_lep_2 & n_tight_lep_2
    mask_TnP        = filters & n_loose_lep_2 & n_tight_lep_al1
    mask_LooseTight = filters & n_loose_lep_2 & n_tight_lep_1 & tight_subleading
    mask_TightLoose = filters & n_loose_lep_2 & n_tight_lep_1 & tight_leading

    # Apply the masks
    events['is_2_loose'] = ak.fill_none(mask_loose,False)
    events['is_loosetight'] = ak.fill_none(mask_LooseTight,False)
    events['is_tightloose'] = ak.fill_none(mask_TightLoose,False)
    events['is_2_tight'] = ak.fill_none(mask_2_tight,False)
    events['is_TnP'] = ak.fill_none(mask_TnP,False)

def addlepSF(events, year, isData, sample_name,is2022):

    # Leptons and padded leptons
    leps_tight = events.lep_tight
    leps_tight_padded = ak.pad_none(leps_tight,2)

    # SFs:
    events['sf_2l_muon'] = leps_tight_padded[:,0].sf_nom_muon*leps_tight_padded[:,1].sf_nom_muon
    events['sf_2l_elec'] = leps_tight_padded[:,0].sf_nom_elec*leps_tight_padded[:,1].sf_nom_elec

def attach_dilepton_preselection_mask(events,lep_collection):
    # Pt requirements (assumes lep_collection is pt sorted and padded)
    pt_mask = ak.fill_none((lep_collection[:,0].pt > 25.0) & (lep_collection[:,1].pt > 20.0),False)
    pt_mask = ak.fill_none(pt_mask,False) # Replace the None with False in the mask just to make it easier to think about

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
    z_mass_mask = ak.any((abs(((lep_collection[:,0:1] + lep_collection[:,1:2]).mass) - get_ec_param("zmass")) < 15.0),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    z_mass_mask = ak.fill_none(z_mass_mask,False) # Replace the None with False in the mask just to make it easier to think about

    of_mass_mask = ak.any((((lep_collection[:,0:1] + lep_collection[:,1:2]).mass) > 20.0),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    of_mass_mask = ak.fill_none(of_mass_mask,False) # Replace the None with False in the mask just to make it easier to think about

    # The final preselection mask
    run3_2lep_presel_mask = (os_mask & pt_mask)

    # Attach to the lepton objects
    events["dilep_presel_sf_ee"] = (run3_2lep_presel_mask & sf_mask & ee_mask & z_mass_mask)
    events["dilep_presel_sf_mumu"] = (run3_2lep_presel_mask & sf_mask & mumu_mask & z_mass_mask)
    events["dilep_presel_of"] = (run3_2lep_presel_mask & of_mask & of_mass_mask)
