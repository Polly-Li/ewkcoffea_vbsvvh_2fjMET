import awkward as ak
from ewkcoffea.modules.paths import ewkcoffea_path
from topcoffea.modules.get_param_from_jsons import GetParam as get_param
get_param = get_param(ewkcoffea_path("params/params.json"))

def get_cleaned_collection(obj_collection_a,obj_collection_b,drcut=0.4):
    obj_b_nearest_to_any_in_a , dr = obj_collection_b.nearest(obj_collection_a,return_metric=True)
    mask = ak.fill_none(dr>drcut,True)
    return mask

def is_veto_run3_2lep_ele(ele):
    mask = (
        (ele.pt                 >  get_param("run3_2lep_pres_e_pt")) &
        (abs(ele.eta)           <  get_param("run3_2lep_pres_e_eta")) &
        (ele.cutBased           >= get_param("run3_2lep_pres_e_cutBasedID_veto")) &
        (ele.convVeto)
    )
    return mask

def is_tight_run3_2lep_ele(ele):
    mask = (
        (ele.pt                                >  get_param("run3_2lep_pres_e_pt")) &
        (abs(ele.eta)                          <  get_param("run3_2lep_pres_e_eta")) &
        (abs(ele.dxy)                          <  get_param("run3_2lep_pres_e_dxy")) &
        (abs(ele.dz)                           <  get_param("run3_2lep_pres_e_dz")) &
        (ele.cutBased                          >= get_param("run3_2lep_pres_e_cutBasedID_med"))
    )
    return mask

def is_veto_run3_2lep_mu(mu):
    mask = (
        (mu.pt               >  get_param("run3_2lep_pres_m_pt")) &
        (abs(mu.eta)         <  get_param("run3_2lep_pres_m_eta")) &
        (mu.looseId)
    )
    return mask

def is_tight_run3_2lep_mu(mu):
    mask = (
        (mu.pt               >  get_param("run3_2lep_pres_m_pt")) &
        (abs(mu.eta)         <  get_param("run3_2lep_pres_m_eta")) &
        (abs(mu.dxy)         <  get_param("run3_2lep_pres_m_dxy")) &
        (abs(mu.dz)          <  get_param("run3_2lep_pres_m_dz")) &
        (mu.pfIsoId          >= get_param("run3_2lep_pres_m_pfIsoId_Tight")) &
        (mu.mediumId)
    )
    return mask

def is_presel_run3_2lep_jets(jets):
    mask = (
        (jets.pt               >  get_param("run3_2lep_pres_jets_pt")) &
        (abs(jets.eta)         <  get_param("run3_2lep_pres_jets_eta"))
    )
    return mask
