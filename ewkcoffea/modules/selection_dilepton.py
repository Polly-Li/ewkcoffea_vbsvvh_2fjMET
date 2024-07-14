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

# 2l selection
def add4lmask_wwz(events,year,isData,sample_name,is2022,is2023):

    # Leptons and padded leptons
    leps = events.l_wwz_t
    leps_padded = ak.pad_none(leps,2)

    # Filters
    filter_flags = events.Flag
    if (is2022 or is2023):
        filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & filter_flags.ecalBadCalibFilter & filter_flags.BadPFMuonDzFilter & filter_flags.hfNoisyHitsFilter & filter_flags.eeBadScFilter
    else:
        filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & (((year == "2016")|(year == "2016APV")) | filter_flags.ecalBadCalibFilter) & (isData | filter_flags.eeBadScFilter)

    # Lep multiplicity
    nlep_2 = (ak.num(leps) == 2)

    mask = filters & nlep_2

    # SFs:
    events['sf_4l_muon'] = leps_padded[:,0].sf_nom_muon*leps_padded[:,1].sf_nom_muon
    events['sf_4l_elec'] = leps_padded[:,0].sf_nom_elec*leps_padded[:,1].sf_nom_elec
    events['sf_4l_hi_muon'] = leps_padded[:,0].sf_hi_muon*leps_padded[:,1].sf_hi_muon
    events['sf_4l_hi_elec'] = leps_padded[:,0].sf_hi_elec*leps_padded[:,1].sf_hi_elec
    events['sf_4l_lo_muon'] = leps_padded[:,0].sf_lo_muon*leps_padded[:,1].sf_lo_muon
    events['sf_4l_lo_elec'] = leps_padded[:,0].sf_lo_elec*leps_padded[:,1].sf_lo_elec

    events['is2lWWZ'] = ak.fill_none(mask,False)

# Do Run3 2Lep pre selection, construct event level mask
def attach_dilepton_preselection_mask(events,lep_collection):

    # Pt requirements (assumes lep_collection is pt sorted and padded)
    pt_mask = ak.fill_none(((lep_collection[:,0].pt > 25.0) & (lep_collection[:,1].pt > 20.0)),False)
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
    events["run3_2lep_presel_sf_ee"] = (run3_2lep_presel_mask & sf_mask & ee_mask & z_mass_mask)
    events["run3_2lep_presel_sf_mumu"] = (run3_2lep_presel_mask & sf_mask & mumu_mask & z_mass_mask)
    events["run3_2lep_presel_of"] = (run3_2lep_presel_mask & of_mask & of_mass_mask)
