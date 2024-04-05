import awkward as ak

from ewkcoffea.modules.paths import ewkcoffea_path
from topcoffea.modules.get_param_from_jsons import GetParam
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))

#-----------------------------------------------------------------------------------------------------------------------------
# 2Lep Selection
def add2lmask_run3_2lep(events, year, isData):

    # Leptons and padded leptons
    leps = events.l_run3_2lep_veto
    leps_padded = ak.pad_none(leps,2)

    # Filters
    filter_flags = events.Flag
    # Filters are from AN SMP-24-001
    filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & filter_flags.ecalBadCalibFilter & filter_flags.BadPFMuonDzFilter & filter_flags.hfNoisyHitsFilter & filter_flags.eeBadScFilter

    # Lep multiplicity
    nlep_2 = (ak.num(leps) == 2)

    mask = filters & nlep_2

    events['muon_sf'] = leps_padded[:,0].sf_nom_muon*leps_padded[:,1].sf_nom_muon
    events['ele_sf'] = leps_padded[:,0].sf_nom_elec*leps_padded[:,1].sf_nom_elec
    events['is2l'] = ak.fill_none(mask,False)

#------------------------------------------------------------------------------------------------------------------------------
# Do Run3 2Lep pre selection, construct event level mask
def attach_run3_2lep_preselection_mask(events,lep_collection):

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


