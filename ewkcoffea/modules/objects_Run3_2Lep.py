import numpy as np
import awkward as ak
import xgboost as xgb
from topcoffea.modules.paths import topcoffea_path
from ewkcoffea.modules.paths import ewkcoffea_path
from topcoffea.modules.get_param_from_jsons import GetParam as get_param
get_param = get_param(ewkcoffea_path("params/params.json"))

##############Commented out for now. Will implement later.#########################################
# Clean collection b (e.g. jets) with collection a (e.g. leps)
#def get_cleaned_collection(obj_collection_a,obj_collection_b,drcut=0.4):
#    obj_b_nearest_to_any_in_a , dr = obj_collection_b.nearest(obj_collection_a,return_metric=True)
#    mask = ak.fill_none(dr>drcut,True)
#    return obj_collection_b[mask]
###################################################################################################




######### Run3 2Lep analysis object selection #########

# Run3 2Lep preselection for electrons
def is_presel_Run3_2Lep_ele(ele):
    mask = (
        (ele.pt               >  get_param("Run3_2Lep_pres_e_pt")) &
        (abs(ele.eta)         <  get_param("Run3_2Lep_pres_e_eta")) &
        (abs(ele.dxy)         <  get_param("Run3_2Lep_pres_e_dxy")) &
        (abs(ele.dz)          <  get_param("Run3_2Lep_pres_e_dz")) &
        (ele.pfRelIso03_all <  get_param("Run3_2Lep_pres_e_pfRelIso_all")) &
        (ele.mvaNoIso_WP90) & 
        ((abs(ele.eta) < 1.444) | (abs(ele.eta) > 1.566))
    )
    return mask


# Run3 2Lep preselection for muons
def is_presel_Run3_2Lep_mu(mu):
    mask = (
        (mu.pt               >  get_param("Run3_2Lep_pres_m_pt")) &
        (abs(mu.eta)         <  get_param("Run3_2Lep_pres_m_eta")) &
        (abs(mu.dxy)         <  get_param("Run3_2Lep_pres_m_dxy")) &
        (abs(mu.dz)          <  get_param("Run3_2Lep_pres_m_dz")) &
        (mu.pfRelIso03_all <  get_param("Run3_2Lep_pres_m_pfRelIso_all")) &
        (mu.mediumId)
    )
    return mask

# Run3 2Lep preselection for jets
def is_presel_Run3_2Lep_jets(jets):
    mask = (
        (jets.pt               >  get_param("Run3_2Lep_pres_jets_pt")) 
    )
    return mask



