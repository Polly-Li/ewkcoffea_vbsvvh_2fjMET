import numpy as np
import awkward as ak
import xgboost as xgb

from topcoffea.modules.paths import topcoffea_path

from ewkcoffea.modules.paths import ewkcoffea_path
from topcoffea.modules.get_param_from_jsons import GetParam
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))

# Clean collection b (e.g. jets) with collection a (e.g. leps)
def get_cleaned_collection(obj_collection_a,obj_collection_b,drcut=0.4):
    obj_a_nearest_to_any_in_b , dr = obj_collection_b.nearest(obj_collection_a,return_metric=True)
    mask = ak.fill_none(dr>drcut,True)
    return obj_collection_b[mask]

######### WWZ 4l analysis object selection #########

# Correctable Jets
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections#Implementation
def get_correctable_jets(jets):
    mask = (
        (jets.pt                       > 15.0) &
        (abs(jets.eta)                 <  5.1)  &
        ((jets.chEmEF + jets.neEmEF)   <  0.9)
    )
    return jets[mask]

# WWZ preselection for electrons
def is_presel_wwz_ele(ele,is2022,is2023):
    mask = (
        (ele.pt               >  get_ec_param("wwz_pres_e_pt"))          &
        (abs(ele.eta)         <  get_ec_param("wwz_pres_e_eta"))         &
        (abs(ele.dxy)         <  get_ec_param("wwz_pres_e_dxy"))         &
        (abs(ele.dz)          <  get_ec_param("wwz_pres_e_dz"))          &
        (abs(ele.sip3d)       <  get_ec_param("wwz_pres_e_sip3d"))       &
        (ele.lostHits         <= get_ec_param("wwz_pres_e_lostHits"))    &
        (ele.tightCharge      == get_ec_param("wwz_pres_e_tightCharge")) &
        (ele.convVeto)
    )

    if (is2022 or is2023):
        mask_year = (
            ele.mvaIso_WP90
        )
    else:
        mask_year = (
            (ele.miniPFRelIso_all < get_ec_param("wwz_pres_e_miniPFRelIso_all"))
        )

    mask_return = mask & mask_year
    return mask_return


# WWZ preselection for muons
def is_presel_wwz_mu(mu,is2022,is2023):
    mask = (
        (mu.pt               >  get_ec_param("wwz_pres_m_pt"))   &
        (abs(mu.eta)         <  get_ec_param("wwz_pres_m_eta"))  &
        (abs(mu.dxy)         <  get_ec_param("wwz_pres_m_dxy"))  &
        (abs(mu.dz)          <  get_ec_param("wwz_pres_m_dz"))   &
        (abs(mu.sip3d)       <  get_ec_param("wwz_pres_m_sip3d"))&
        (mu.mediumId)
    )

    if (is2022 or is2023):
        mask_year = (
            (mu.pfIsoId >= get_ec_param("run3_wwz_pres_m_pfIsoId_Loose"))
        )
    else:
        mask_year = (
            (mu.miniPFRelIso_all < get_ec_param("wwz_pres_m_miniPFRelIso_all"))
        )

    mask_return = mask & mask_year
    return mask_return


# Get MVA score from TOP MVA for electrons
def get_topmva_score_ele(events, year):

    ele = events.Electron

    # Get the model path
    if (year == "2016"):      ulbase = "UL16"
    elif (year == "2016APV"): ulbase = "UL16APV"
    elif (year == "2017"):    ulbase = "UL17"
    elif (year == "2018"):    ulbase = "UL18"
    else: raise Exception(f"Error: Unknown year \"{year}\". Exiting...")
    model_fpath = topcoffea_path(f"data/topmva/lepid_weights/el_TOP{ulbase}_XGB.weights.bin")

    # Get the input data
    ele["btagDeepFlavB"] = ak.fill_none(ele.matched_jet.btagDeepFlavB, 0)
    ele["jetPtRatio"] = 1./(ele.jetRelIso+1.)
    ele["miniPFRelIso_diff_all_chg"] = ele.miniPFRelIso_all - ele.miniPFRelIso_chg
    # The order here comes from https://github.com/cmstas/VVVNanoLooper/blob/8a194165cdbbbee3bcf69f932d837e95a0a265e6/src/ElectronIDHelper.cc#L110-L122
    in_vals = np.array([
        ak.flatten(ele.pt),
        ak.flatten(ele.eta), # Kirill confirms that signed eta was used in the training
        ak.flatten(ele.jetNDauCharged),
        ak.flatten(ele.miniPFRelIso_chg),
        ak.flatten(ele.miniPFRelIso_diff_all_chg),
        ak.flatten(ele.jetPtRelv2),
        ak.flatten(ele.jetPtRatio),
        ak.flatten(ele.pfRelIso03_all),
        ak.flatten(ele.btagDeepFlavB),
        ak.flatten(ele.sip3d),
        ak.flatten(np.log(abs(ele.dxy))),
        ak.flatten(np.log(abs(ele.dz))),
        ak.flatten(ele.mvaFall17V2noIso),
    ])
    in_vals = np.transpose(in_vals) # To go from e.g. [ [pt1,pt1] , [eta1,eta2] ] -> [ [pt1,eta1] , [pt2,eta2] ]
    in_vals = xgb.DMatrix(in_vals) # The format xgb expects

    # Load model and evaluate
    xgb.set_config(verbosity = 0)
    bst = xgb.Booster()
    bst.load_model(model_fpath)
    score = bst.predict(in_vals).reshape(-1)

    # Restore the shape (i.e. unflatten)
    counts = ak.num(ele.pt)
    score = ak.unflatten(score,counts)
    return score


# Get MVA score from TOP MVA for muons
def get_topmva_score_mu(events, year):

    mu = events.Muon

    # Get the model path
    if (year == "2016"):      ulbase = "UL16"
    elif (year == "2016APV"): ulbase = "UL16APV"
    elif (year == "2017"):    ulbase = "UL17"
    elif (year == "2018"):    ulbase = "UL18"
    else: raise Exception(f"Error: Unknown year \"{year}\". Exiting...")
    model_fpath = topcoffea_path(f"data/topmva/lepid_weights/mu_TOP{ulbase}_XGB.weights.bin")

    # Get the input data
    mu["btagDeepFlavB"] = ak.fill_none(mu.matched_jet.btagDeepFlavB, 0)

    mu["jetPtRatio"] = 1./(mu.jetRelIso+1.)
    mu["miniPFRelIso_diff_all_chg"] = mu.miniPFRelIso_all - mu.miniPFRelIso_chg
    in_vals = np.array([
        ak.flatten(mu.pt),
        ak.flatten(mu.eta), # Kirill confirms that signed eta was used in the training
        ak.flatten(mu.jetNDauCharged),
        ak.flatten(mu.miniPFRelIso_chg),
        ak.flatten(mu.miniPFRelIso_diff_all_chg),
        ak.flatten(mu.jetPtRelv2),
        ak.flatten(mu.jetPtRatio),
        ak.flatten(mu.pfRelIso03_all),
        ak.flatten(mu.btagDeepFlavB),
        ak.flatten(mu.sip3d),
        ak.flatten(np.log(abs(mu.dxy))),
        ak.flatten(np.log(abs(mu.dz))),
        ak.flatten(mu.segmentComp),
    ])
    in_vals = np.transpose(in_vals)
    in_vals = xgb.DMatrix(in_vals)

    # Load model and evaluate
    xgb.set_config(verbosity = 0)
    bst = xgb.Booster()
    bst.load_model(model_fpath)
    score = bst.predict(in_vals).reshape(-1)

    # Restore the shape (i.e. unflatten)
    counts = ak.num(mu.pt)
    score = ak.unflatten(score,counts)
    return score

#################### Run 2 and Run 3 TnP IDs #######################

def is_loose_run2_ele_tnp(ele):
    mask = (
        (ele.pt               >  get_ec_param("wwz_pres_e_pt"))  &
        (abs(ele.eta)         <  get_ec_param("wwz_pres_e_eta")) &
        (abs(ele.dxy)         <  get_ec_param("wwz_pres_e_dxy")) &
        (abs(ele.dz)          <  get_ec_param("wwz_pres_e_dz"))  &
        (abs(ele.sip3d)       <  get_ec_param("wwz_pres_e_sip3d")) &
        (ele.lostHits         <= get_ec_param("wwz_pres_e_lostHits"))&
        (ele.tightCharge      == get_ec_param("wwz_pres_e_tightCharge")) &
        (ele.miniPFRelIso_all <  get_ec_param("wwz_pres_e_miniPFRelIso_all")) &
        (ele.convVeto)
    )
    return mask

def is_loose_run2_mu_tnp(mu):
    mask = (
        (mu.pt               >  get_ec_param("wwz_pres_m_pt")) &
        (abs(mu.eta)         <  get_ec_param("wwz_pres_m_eta")) &
        (abs(mu.dxy)         <  get_ec_param("wwz_pres_m_dxy")) &
        (abs(mu.dz)          <  get_ec_param("wwz_pres_m_dz")) &
        (abs(mu.sip3d)       <  get_ec_param("wwz_pres_m_sip3d")) &
        (mu.miniPFRelIso_all <  get_ec_param("wwz_pres_m_miniPFRelIso_all")) &
        (mu.mediumId)
    )
    return mask

def is_loose_run3_ele_tnp(ele):
    mask = (
        (ele.pt               >  get_ec_param("wwz_pres_e_pt"))  &
        (abs(ele.eta)         <  get_ec_param("wwz_pres_e_eta")) &
        (abs(ele.dxy)         <  get_ec_param("wwz_pres_e_dxy")) &
        (abs(ele.dz)          <  get_ec_param("wwz_pres_e_dz"))  &
        (abs(ele.sip3d)       <  get_ec_param("wwz_pres_e_sip3d")) &
        (ele.lostHits         <= get_ec_param("wwz_pres_e_lostHits"))&
        (ele.tightCharge      == get_ec_param("wwz_pres_e_tightCharge")) &
        (ele.cutBased         >= 1) &
        (ele.convVeto)
    )
    return mask

def is_loose_run3_mu_tnp(mu):
    mask = (
        (mu.pt               >  get_ec_param("wwz_pres_m_pt")) &
        (abs(mu.eta)         <  get_ec_param("wwz_pres_m_eta")) &
        (abs(mu.dxy)         <  get_ec_param("wwz_pres_m_dxy")) &
        (abs(mu.dz)          <  get_ec_param("wwz_pres_m_dz")) &
        (abs(mu.sip3d)       <  get_ec_param("wwz_pres_m_sip3d")) &
        #(mu.pfIsoId          >= 1) &
        (mu.mediumId)
    )
    return mask

def is_tight_run3_ele_tnp(ele):
    mask = (
        (ele.pt               >  get_ec_param("wwz_pres_e_pt"))  &
        (abs(ele.eta)         <  get_ec_param("wwz_pres_e_eta")) &
        (abs(ele.dxy)         <  get_ec_param("wwz_pres_e_dxy")) &
        (abs(ele.dz)          <  get_ec_param("wwz_pres_e_dz"))  &
        (abs(ele.sip3d)       <  get_ec_param("wwz_pres_e_sip3d")) &
        (ele.lostHits         <= get_ec_param("wwz_pres_e_lostHits"))&
        (ele.tightCharge      == get_ec_param("wwz_pres_e_tightCharge")) &
        #(ele.cutBased         >= 1) &
        (ele.mvaIso_WP90) &
        (ele.convVeto)
    )
    return mask

def is_tight_run3_mu_tnp(mu):
    mask = (
        (mu.pt               >  get_ec_param("wwz_pres_m_pt")) &
        (abs(mu.eta)         <  get_ec_param("wwz_pres_m_eta")) &
        (abs(mu.dxy)         <  get_ec_param("wwz_pres_m_dxy")) &
        (abs(mu.dz)          <  get_ec_param("wwz_pres_m_dz")) &
        (abs(mu.sip3d)       <  get_ec_param("wwz_pres_m_sip3d")) &
        (mu.pfIsoId          >= 2) &
        (mu.mediumId)
    )
    return mask



