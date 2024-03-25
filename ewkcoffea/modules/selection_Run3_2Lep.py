import awkward as ak

from ewkcoffea.modules.paths import ewkcoffea_path
from topcoffea.modules.get_param_from_jsons import GetParam
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))

# The datasets we are using, and the triggers in them
dataset_dict = {


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
            "DoublePhoton70",
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
        "Muon" : [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            "IsoMu24",
            "IsoMu27",
            "Mu50",
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
exclude_dict = {
    "B": {
        "SingleMuon"     : [],
        "DoubleMuon"     : dataset_dict["2022"]["SingleMuon"],
        "EGamma"         : dataset_dict["2022"]["SingleMuon"] + dataset_dict["2022"]["DoubleMuon"],
        "MuonEG"         : dataset_dict["2022"]["SingleMuon"] + dataset_dict["2022"]["DoubleMuon"] + dataset_dict["2022"]["EGamma"],
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

#-----------------------------------------------------------------------------------------------------------------------------
# 2Lep Selection
def add2lmask_Run3_2Lep(events, year, isData):

    # Leptons and padded leptons
    leps = events.l_Run3_2Lep_veto
    leps_padded = ak.pad_none(leps,2)

    # Filters
    filter_flags = events.Flag
    filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & filter_flags.ecalBadCalibFilter & filter_flags.BadPFMuonDzFilter & filter_flags.hfNoisyHitsFilter & filter_flags.eeBadScFilter

    # Lep multiplicity
    nlep_2 = (ak.num(leps) == 2)

    mask = filters & nlep_2

    events['muon_sf'] = leps_padded[:,0].sf_nom_muon*leps_padded[:,1].sf_nom_muon
    events['ele_sf'] = leps_padded[:,0].sf_nom_elec*leps_padded[:,1].sf_nom_elec
    events['is2l'] = ak.fill_none(mask,False)

#------------------------------------------------------------------------------------------------------------------------------
# Do Run3 2Lep pre selection, construct event level mask
# Convenience function around get_Run3_2Lep_candidates() and get_z_candidate_mask()
def attach_Run3_2Lep_preselection_mask(events,lep_collection):

    # Pt requirements (assumes lep_collection is pt sorted and padded)
    pt_mask = ak.fill_none((lep_collection[:,0].pt > 25.0),False)
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
    z_mass_mask = ak.any((abs(((lep_collection[:,0:1] + lep_collection[:,1:2]).mass) - 91.2) < 15.0),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    z_mass_mask = ak.fill_none(z_mass_mask,False) # Replace the None with False in the mask just to make it easier to think about

    of_mass_mask = ak.any((((lep_collection[:,0:1] + lep_collection[:,1:2]).mass) > 20.0),axis=1) # Use ak.any() here so that instead of e.g [[None],None,...] we have [False,None,...]
    of_mass_mask = ak.fill_none(of_mass_mask,False) # Replace the None with False in the mask just to make it easier to think about

    # The final preselection mask
    Run3_2Lep_presel_mask = (os_mask & pt_mask)

    # Attach to the lepton objects
    events["Run3_2Lep_presel_sf_ee"] = (Run3_2Lep_presel_mask & sf_mask & ee_mask & z_mass_mask)
    events["Run3_2Lep_presel_sf_mumu"] = (Run3_2Lep_presel_mask & sf_mask & mumu_mask & z_mass_mask)
    events["Run3_2Lep_presel_of"] = (Run3_2Lep_presel_mask & of_mask & of_mass_mask)


