from hist import axis
from coffea.nanoevents.methods import vector
import awkward as ak
ak.behavior.update(vector.behavior)
import numpy as np

obj = { #reconstruction of H, V1, V2, vbsjs and MET
    "Higgs": lambda events: ak.zip(
        {
            "pt": events.Higgs_pt,
            "mass": events.Higgs_msoftdrop,
            "eta": events.Higgs_eta,
            "phi": events.Higgs_phi,
            "score": events.HiggsScore,
            "second_score": events.V1Score,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "V1": lambda events: ak.zip(
        {
            "pt": events.V1_pt,
            "mass": events.V1_msoftdrop,
            "eta": events.V1_eta,
            "phi": events.V1_phi,
            "score": events.V1Score,
            "second_score": events.HiggsScore,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "vbsj1": lambda events: ak.zip(
        {
            "pt": events.vbsj1_pt,
            "mass": events.vbsj1_m,
            "eta": events.vbsj1_eta,
            "phi": events.vbsj1_phi,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "vbsj2": lambda events: ak.zip(
        {
            "pt": events.vbsj2_pt,
            "mass": events.vbsj2_m,
            "eta": events.vbsj2_eta,
            "phi": events.vbsj2_phi,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "MET": lambda events: ak.zip(
        {
            "pt": events.Met_pt,
            "phi": events.Met_phi,
            "eta": ak.zeros_like(events.Met_pt),  # Set to 0 to ensure valid LorentzVector
            "mass": ak.zeros_like(events.Met_pt), # Same for mass
            "significance": events.Met_significance,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "leadAK8": lambda events: ak.zip(
        {
            "pt": events.leadAK8_pt,
            "mass": events.leadAK8_msoftdrop,
            "eta": events.leadAK8_eta,
            "phi": events.leadAK8_phi,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "goodAK4Jets": lambda events: ak.zip(
        {
            "pt": events.goodAK4Jets_pt,
            "mass": events.goodAK4Jets_mass,
            "eta": events.goodAK4Jets_eta,
            "phi": events.goodAK4Jets_phi,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "centralAK4Jets": lambda events: ak.zip(
        {
            "pt": events.goodAK4Jets_pt[abs(events.goodAK4Jets_eta) < 2.4],
            "mass": events.goodAK4Jets_mass[abs(events.goodAK4Jets_eta) < 2.4],
            "eta": events.goodAK4Jets_eta[abs(events.goodAK4Jets_eta) < 2.4],
            "phi": events.goodAK4Jets_phi[abs(events.goodAK4Jets_eta) < 2.4],
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "centralVBSJets": lambda events: ak.zip(
        {
            "pt": events.goodVBSJets_pt[abs(events.goodVBSJets_eta) < 2.4],
            "mass": events.goodVBSJets_mass[abs(events.goodVBSJets_eta) < 2.4],
            "eta": events.goodVBSJets_eta[abs(events.goodVBSJets_eta) < 2.4],
            "phi": events.goodVBSJets_phi[abs(events.goodVBSJets_eta) < 2.4],
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "leadAK4": lambda events: ak.firsts(
        ak.zip(
            {
                "pt": events.goodAK4Jets_pt,
                "mass": events.goodAK4Jets_mass,
                "eta": events.goodAK4Jets_eta,
                "phi": events.goodAK4Jets_phi,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior
        )[ak.argsort(events.goodAK4Jets_pt, axis=1, ascending=False)]
    ),
}

other_objs = { #reconstruction of all (good) jets and stuffs if needed
    "goodAK8Jets": lambda events: ak.zip(
        {
            "pt": events.goodAK8Jets_pt,
            "mass": events.goodAK8Jets_msoftdrop,
            "eta": events.goodAK8Jets_eta,
            "phi": events.goodAK8Jets_phi,
            "HbbScore":events.goodAK8Jets_HbbScore,
            "WqqScore":events.goodAK8Jets_WqqScore,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "goodVBSJets": lambda events: ak.zip(
        {
            "pt": events.goodAK4Jets_pt,
            "mass": events.goodAK4Jets_mass,
            "eta": events.goodAK4Jets_eta,
            "phi": events.goodAK4Jets_phi,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    "V2": lambda events: ak.zip( #for met2FJ channel, V2 should not exist
        {
            "pt": events.V2_pt,
            "mass": events.V2_msoftdrop,
            "eta": events.V2_eta,
            "phi": events.V2_phi,
            "score": events.V2Score,
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior
    ),
    
}

dense_variables_config = { #name of axis must be same as key    
    "count_plot": {
        "axis": axis.Regular(50, -1, 1, name="count_plot", label="count_plot"),
        "expr": lambda events, objects:ak.ones_like(events.weight), 
    },
    "nGoodAK4": {
        "axis": axis.Regular(25, 0, 25, name="nGoodAK4", label="nGoodAK4"),
        "expr": lambda events, objects: events.nGoodAK4,
    },
    "nGoodAK8": {
        "axis": axis.Regular(6, 0, 6, name="nGoodAK8", label="nGoodAK8"),
        "expr": lambda events, objects: events.nGoodAK8,
    },
    "nAK4": {
        "axis": axis.Regular(25, 0, 25, name="nAK4", label="nAK4"),
        "expr": lambda events, objects: events.nAK4
    },
    "nAK8": {
        "axis": axis.Regular(6, 0, 6, name="nAK8", label="nAK8"),
        "expr": lambda events, objects: events.nAK8,
    },
    "Higgs_pt": {
        "axis": axis.Regular(50, 0, 2000, name="Higgs_pt", label="Higgs pt (GeV)"),
        "expr":  lambda events, objects: objects["Higgs"].pt,
    },
    "Higgs_phi": {
        "axis": axis.Regular(50, -3.5, 3.5, name="Higgs_phi", label="Higgs phi"),
        "expr":  lambda events, objects: objects["Higgs"].phi,
    },
    "Higgs_eta": {
        "axis": axis.Regular(50, 0,6, name="Higgs_eta", label="Higgs eta"),
        "expr":  lambda events, objects: objects["Higgs"].eta,
    },
    "Higgs_mass": {
        "axis": axis.Regular(50, 0, 500, name="Higgs_mass", label="Higgs mass (GeV)"),
        "expr":  lambda events, objects: objects["Higgs"].mass,
    },
    "Higgs_score": {
        "axis": axis.Regular(50, 0, 1, name="Higgs_score", label="Higgs score"),
        "expr":  lambda events, objects: objects["Higgs"].score,
    },
    "Higgs_second_score": {
        "axis": axis.Regular(50, 0, 1, name="Higgs_second_score", label="Higgs' V1 score"),
        "expr":  lambda events, objects: objects["Higgs"].second_score,
    },
    "V1_pt": {
        "axis": axis.Regular(50, 0, 2000, name="V1_pt", label="V1 pt (GeV)"),
        "expr":  lambda events, objects: objects["V1"].pt,
    },
    "V1_phi": {
        "axis": axis.Regular(50, -3.5, 3.5, name="V1_phi", label="V1 phi"),
        "expr":  lambda events, objects: objects["V1"].phi,
    },
    "V1_eta": {
        "axis": axis.Regular(50, 0,6, name="V1_eta", label="V1 eta"),
        "expr":  lambda events, objects: objects["V1"].eta,
    },
    "V1_mass": {
        "axis": axis.Regular(50, 0, 500, name="V1_mass", label="V1 mass (GeV)"),
        "expr":  lambda events, objects: objects["V1"].mass,
    },
    "V1_score": {
        "axis": axis.Regular(50, 0, 1, name="V1_score", label="V1 score"),
        "expr":  lambda events, objects: objects["V1"].score,
    },
    "V1_second_score": {
        "axis": axis.Regular(50, 0, 1, name="V1_second_score", label="V1's Higgs score"),
        "expr":  lambda events, objects: objects["V1"].second_score,
    },
    "Met_pt":{
        "axis": axis.Regular(50, 0, 2000, name="Met_pt", label="MET pt (GeV)"),
        "expr":  lambda events, objects: objects["MET"].pt,
    },
    "Met_significance":{
        "axis": axis.Regular(50, 0, 2000, name="Met_significance", label="MET significance"),
        "expr":  lambda events, objects: objects["MET"].significance,
    },
    "Met_significance_fine":{
        "axis": axis.Regular(50, 0, 200, name="Met_significance_fine", label="MET significance zoomed in"),
        "expr":  lambda events, objects: objects["MET"].significance,
    },
    "Met_phi":{
        "axis": axis.Regular(50, -3.5, 3.5, name="Met_phi", label="MET phi"),
        "expr":  lambda events, objects: objects["MET"].phi,
    },
    "HV1_dR": {
        "axis": axis.Regular(50, 0, 6, name="HV1_dR", label="dR(Higgs, V1)"),
        "expr": lambda events, objects: objects["Higgs"].delta_r(objects["V1"]),
    },
    "HMET_dphi": {
        "axis": axis.Regular(50, 0, 3.5, name="HMET_dphi", label="dphi(Higgs, MET)"),
        "expr": lambda events, objects: deltaPhi(objects["Higgs"],objects["MET"]),
    },
    "V1MET_dphi": {
        "axis": axis.Regular(50, 0, 3.5, name="V1MET_dphi", label="dphi(V1, MET)"),
        "expr": lambda events, objects: deltaPhi(objects["V1"],objects["MET"]),
    },
    "HVMET_sum_pt": {
        "axis": axis.Regular(50, 0, 2000, name="HVMET_sum_pt", label="pt(H + V1 + MET) (GeV)"),
        "expr": lambda events, objects: (objects["Higgs"] + objects["V1"] + objects["MET"]).pt,
    },
    "sum_bosonHT":{
        "axis": axis.Regular(50, 0, 5000, name="sum_bosonHT", label="(Higgs_pt+V1_pt+MET_pt)"),
        "expr": lambda events, objects: objects["Higgs"].pt + objects["V1"].pt +objects["MET"].pt,
    },
    "vbsj_deta":{
        "axis": axis.Regular(50, 0, 10, name="vbsj_deta", label="vbsj_deta"),
        "expr": lambda events, objects: np.abs(objects["vbsj1"].eta - objects["vbsj2"].eta),
    },
    "vbsj_Mjj":{
        "axis": axis.Regular(50, 0, 5000, name="vbsj_Mjj", label="vbsj_Mjj (GeV)"),
        "expr": lambda events, objects: (objects["vbsj1"] + objects["vbsj2"]).mass,
    },
    
    "max_deta":{
        "axis": axis.Regular(50, 0, 10, name="max_deta", label="max_deta"),
        "expr": lambda events, objects: events.goodAK4Jets_maxAbsDeltaEta,
    },
    "leadAK8_MET_dphi":{
        "axis": axis.Regular(50, 0, 3.5, name="leadAK8_MET_dphi", label="dphi(leadAK8, MET)"),
        "expr": lambda events, objects: deltaPhi(objects["leadAK8"],objects["MET"]),
    },
    "sum_dphi":{
        "axis": axis.Regular(50, 0, 6.5, name="sum_dphi", label="sum_dphi_of_HV1MET"),
        "expr": lambda events, objects: deltaPhi(objects["Higgs"],objects["MET"]) + deltaPhi(objects["V1"],objects["MET"]) + deltaPhi(objects["V1"],objects["Higgs"]),
    },
    
    "HV_sum_pt": {
        "axis": axis.Regular(50, 0, 2000, name="HV_sum_pt", label="pt(H + V1) (GeV)"),
        "expr": lambda events, objects: (objects["Higgs"] + objects["V1"]).pt,
    },
    "sum_HV1_HT":{
        "axis": axis.Regular(50, 0, 5000, name="sum_HV1_HT", label="(Higgs_pt+V1_pt)"),
        "expr": lambda events, objects: objects["Higgs"].pt + objects["V1"].pt,
    },
    "nCentralAK4": {
        "axis": axis.Regular(10, 0, 10, name="nCentralAK4", label="nCentralAK4(not excluding by dR(H or V1))"),
        "expr": lambda events, objects: ak.num(objects["centralAK4Jets"]),
    },
    "nCentralvbs": {
        "axis": axis.Regular(10, 0, 10, name="nCentralvbs", label="nCentralvbs(excluding dR(H or V1)<0.8)"),
        "expr": lambda events, objects: ak.num(objects["centralVBSJets"]),
    },
    "leadAK4_MET_dphi":{
        "axis": axis.Regular(50, 0, 3.5, name="leadAK4_MET_dphi", label="dphi(leadAK4, MET)"),
        "expr": lambda events, objects: deltaPhi(objects["leadAK4"],objects["MET"]),
    },
    "vbsj1MET_dphi": {
        "axis": axis.Regular(50, 0, 3.5, name="vbsj1MET_dphi", label="dphi(vbsj1, MET)"),
        "expr": lambda events, objects: deltaPhi(objects["vbsj1"],objects["MET"]),
    },
    "vbsj2MET_dphi": {
        "axis": axis.Regular(50, 0, 3.5, name="vbsj2MET_dphi", label="dphi(vbsj2, MET)"),
        "expr": lambda events, objects: deltaPhi(objects["vbsj2"],objects["MET"]),
    },
    
    "Hvbsj1_dR": {
        "axis": axis.Regular(50, 0, 6, name="Hvbsj1_dR", label="dR(Higgs, vbsj1)"),
        "expr": lambda events, objects: objects["Higgs"].delta_r(objects["vbsj1"]),
    },
    "Hvbsj2_dR": {
        "axis": axis.Regular(50, 0, 6, name="Hvbsj2_dR", label="dR(Higgs, vbsj2)"),
        "expr": lambda events, objects: objects["Higgs"].delta_r(objects["vbsj2"]),
    },
    "V1vbsj1_dR": {
        "axis": axis.Regular(50, 0, 6, name="V1vbsj1_dR", label="dR(V1, vbsj1)"),
        "expr": lambda events, objects: objects["V1"].delta_r(objects["vbsj1"]),
    },
    "V1vbsj2_dR": {
        "axis": axis.Regular(50, 0, 6, name="V1vbsj2_dR", label="dR(V1, vbsj2)"),
        "expr": lambda events, objects: objects["V1"].delta_r(objects["vbsj2"]),
    },
    "vbsj_dR": {
        "axis": axis.Regular(50, 0, 6, name="vbsj_dR", label="dR(vbsj1, vbsj2)"),
        "expr": lambda events, objects: objects["vbsj2"].delta_r(objects["vbsj1"]),
    },
}

# seems delta_r is the only available function for coffea?  
def deltaR(v1, v2): 
    return v1.delta_r(v2)

def deltaPhi(v1, v2):
    phi1 = v1.phi
    phi2 = v2.phi
    abs_diff = np.abs(phi1 - phi2)
    dphi = ak.where(
        abs_diff < np.pi,
        abs_diff,
        2 * np.pi - abs_diff
    ) #compare element-wise
    return dphi

def deltaEta(v1,v2):
    eta1 = v1.eta
    eta2 = v2.eta
    return np.abs(eta1-eta2)