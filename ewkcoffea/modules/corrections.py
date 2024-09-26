import numpy as np
import pickle
import gzip
import awkward as ak

import correctionlib
from coffea import lookup_tools

from topcoffea.modules.paths import topcoffea_path
from ewkcoffea.modules.paths import ewkcoffea_path
from topcoffea.modules.CorrectedJetsFactory import CorrectedJetsFactory

extLepSF = lookup_tools.extractor()

#Map for accessing topcoffea clib files TODO Update methods to use these files and remove from ewkcoffea
clib_year_map = {
    "2016APV": "2016preVFP_UL",
    "2016preVFP": "2016preVFP_UL",
    "2016": "2016postVFP_UL",
    "2017": "2017_UL",
    "2018": "2018_UL",
    "2022": "2022_Summer22",
    "2022EE": "2022_Summer22EE",
    "2023": "2022_Summer23",
    "2023BPix": "2022_Summer23BPix",
}

#JERC dictionary for various keys
jerc_dict = {
    "2016": {
        "jec_mc"  : "Summer19UL16_V7_MC",
        "jec_data": "Summer19UL16_RunFGH_V7_DATA",
        "jer"     : "Summer20UL16_JRV3"
    },
    "2016APV": {
        "jec_mc": "Summer19UL16APV_V7_MC",
        "jec_data": {
            "B": "Summer19UL16APV_RunBCD_V7_DATA",
            "C": "Summer19UL16APV_RunBCD_V7_DATA",
            "D": "Summer19UL16APV_RunBCD_V7_DATA",
            "E": "Summer19UL16APV_RunEF_V7_DATA",
            "F": "Summer19UL16APV_RunEF_V7_DATA",
        },
        "jer": "Summer20UL16APV_JRV3"             
    },
    "2017": {
        "jec_mc": "Summer19UL17_V5_MC",
        "jec_data": {
            "B": "Summer19UL17_RunB_V5_DATA",
            "C": "Summer19UL17_RunC_V5_DATA",
            "D": "Summer19UL17_RunD_V5_DATA",
            "E": "Summer19UL17_RunE_V5_DATA",
            "F": "Summer19UL17_RunF_V5_DATA",
        },
        "jer": "Summer19UL17_JRV2"
    },
    "2018": {
        "jec_mc": "Summer19UL18_V5_MC",
        "jec_data": {
            "A": "Summer19UL18_RunA_V5_DATA",
            "B": "Summer19UL18_RunB_V5_DATA",
            "C": "Summer19UL18_RunC_V5_DATA",
            "D": "Summer19UL18_RunD_V5_DATA",
        },
        "jer": "Summer19UL18_JRV2"
    },
    "2022": {
        "jec_mc"  : "Summer22_22Sep2023_V2_MC",
        "jec_data": "Summer22_22Sep2023_RunCD_V2_DATA",
        "jer"     : "Summer22_22Sep2023_JRV1"
    },
    "2022EE": {
        "jec_mc": "Summer22EE_22Sep2023_V2_MC",
        "jec_data": {
            "E": "Summer22EE_22Sep2023_RunE_V2_DATA",
            "F": "Summer22EE_22Sep2023_RunF_V2_DATA",
            "G": "Summer22EE_22Sep2023_RunG_V2_DATA",
        },
        "jer": "Summer22EE_22Sep2023_JRV1"             
    },
    "2023": {
        "jec_mc": "Summer23Prompt23_V1_MC",
        "jec_data": {
            "C1": "Summer23Prompt23_RunCv123_V1_DATA",
            "C2": "Summer23Prompt23_RunCv123_V1_DATA",
            "C3": "Summer23Prompt23_RunCv123_V1_DATA",
            "C4": "Summer23Prompt23_RunCv4_V1_DATA",
        },
        "jer": "Summer23Prompt23_RunCv1234_JRV1"             
    },
    "2023BPix": {
        "jec_mc"  : "Summer23BPixPrompt23_V1_MC",
        "jec_data": "Summer23BPixPrompt23_RunD_V1_DATA",
        "jer"     : "Summer23BPixPrompt23_RunD_JRV1"
    }
}

# Method for extracting the jerc keys and jet algorithm from the dictionary above
def get_jerc_keys(year,era,isdata):

    # Jet Algorithm
    if year in ['2022','2022EE','2023','2023BPix']:
        jet_algo = 'AK4PFPuppi'
    elif year in ['2016','2016APV','2017','2018']:
        jet_algo = 'AK4PFchs'
    # jerc keys
    if not isdata:
        jec_key   = jerc_dict[year]['jec_mc']
        jer_key      = jerc_dict[year]['jer']
    else:
        if year in ['2016','2022','2023BPix']:
            jec_key = jerc_dict[year]['jec_data']
            jer_key      = None
        else:
            jec_key = jerc_dict[year]['jec_data'][era]
            jer_key      = None
    return jet_algo,jec_key,jer_key

###### Muon: tight (topmva) ######
extLepSF.add_weight_sets(["MuonTightSF_2016 NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_value %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL16.json')])
extLepSF.add_weight_sets(["MuonTightSF_2016APV NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_value %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL16APV.json')])
extLepSF.add_weight_sets(["MuonTightSF_2017 NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_value %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL17.json')])
extLepSF.add_weight_sets(["MuonTightSF_2018 NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_value %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL18.json')])
# Syst uncertainty
extLepSF.add_weight_sets(["MuonTightSF_2016_syst NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_syst %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL16.json')])
extLepSF.add_weight_sets(["MuonTightSF_2016APV_syst NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_syst %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL16APV.json')])
extLepSF.add_weight_sets(["MuonTightSF_2017_syst NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_syst %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL17.json')])
extLepSF.add_weight_sets(["MuonTightSF_2018_syst NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_syst %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL18.json')])
# Stat uncertainty
extLepSF.add_weight_sets(["MuonTightSF_2016_stat NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_stat %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL16.json')])
extLepSF.add_weight_sets(["MuonTightSF_2016APV_stat NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_stat %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL16APV.json')])
extLepSF.add_weight_sets(["MuonTightSF_2017_stat NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_stat %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL17.json')])
extLepSF.add_weight_sets(["MuonTightSF_2018_stat NUM_LeptonMvaTight_DEN_TrackerMuons/abseta_pt_stat %s" % ewkcoffea_path('data/topmva_lep_sf/NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt_UL18.json')])


###### Electron: tight (topmva) ######
extLepSF.add_weight_sets(["EleTightSF_2016 EGamma_SF2D %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL16.root')])
extLepSF.add_weight_sets(["EleTightSF_2016APV EGamma_SF2D %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL16APV.root')])
extLepSF.add_weight_sets(["EleTightSF_2017 EGamma_SF2D %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL17.root')])
extLepSF.add_weight_sets(["EleTightSF_2018 EGamma_SF2D %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL18.root')])
# Syst uncertainty
extLepSF.add_weight_sets(["EleTightSF_2016_syst sys %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL16.root')])
extLepSF.add_weight_sets(["EleTightSF_2016APV_syst sys %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL16APV.root')])
extLepSF.add_weight_sets(["EleTightSF_2017_syst sys %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL17.root')])
extLepSF.add_weight_sets(["EleTightSF_2018_syst sys %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL18.root')])
# Stat uncertainty
extLepSF.add_weight_sets(["EleTightSF_2016_stat stat %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL16.root')])
extLepSF.add_weight_sets(["EleTightSF_2016APV_stat stat %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL16APV.root')])
extLepSF.add_weight_sets(["EleTightSF_2017_stat stat %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL17.root')])
extLepSF.add_weight_sets(["EleTightSF_2018_stat stat %s" % ewkcoffea_path('data/topmva_lep_sf/egammaEffi_txt_EGM2D_UL18.root')])


extLepSF.finalize()
SFevaluator = extLepSF.make_evaluator()


def AttachMuonSF(muons, year):
    '''
      Description:
          Inserts 'sf_nom', 'sf_hi', and 'sf_lo' into the muons array passed to this function. These
          values correspond to the nominal, up, and down muon scalefactor values respectively.
    '''

    if year not in ['2016','2016APV','2017','2018']: raise Exception(f"Error: Unknown year \"{year}\".")
    eta = np.abs(muons.eta)
    pt = muons.pt

    tight_sf   = SFevaluator[f'MuonTightSF_{year}'](eta,pt)
    tight_syst = SFevaluator[f'MuonTightSF_{year}_syst'](eta,pt)
    tight_stat = SFevaluator[f'MuonTightSF_{year}_stat'](eta,pt)
    tight_err  = np.sqrt(tight_syst*tight_syst + tight_stat*tight_stat)

    muons['sf_nom_muon'] = tight_sf
    muons['sf_hi_muon']  = (tight_sf + tight_err)
    muons['sf_lo_muon']  = (tight_sf - tight_err)
    muons['sf_nom_elec'] = ak.ones_like(tight_sf)
    muons['sf_hi_elec']  = ak.ones_like(tight_sf)
    muons['sf_lo_elec']  = ak.ones_like(tight_sf)


def AttachElectronSF(electrons, year):
    '''
      Description:
          Inserts 'sf_nom', 'sf_hi', and 'sf_lo' into the electrons array passed to this function. These
          values correspond to the nominal, up, and down electron scalefactor values respectively.
    '''

    if year not in ['2016','2016APV','2017','2018']:
        raise Exception(f"Error: Unknown year \"{year}\".")
    eta = electrons.eta
    pt = electrons.pt

    tight_sf   = SFevaluator[f'EleTightSF_{year}'](eta,pt)
    tight_syst = SFevaluator[f'EleTightSF_{year}_syst'](eta,pt)
    tight_stat = SFevaluator[f'EleTightSF_{year}_stat'](eta,pt)
    tight_err  = np.sqrt(tight_syst*tight_syst + tight_stat*tight_stat)

    electrons['sf_nom_elec'] = tight_sf
    electrons['sf_hi_elec']  = (tight_sf + tight_err)
    electrons['sf_lo_elec']  = (tight_sf - tight_err)
    electrons['sf_nom_muon'] = ak.ones_like(tight_sf)
    electrons['sf_hi_muon']  = ak.ones_like(tight_sf)
    electrons['sf_lo_muon']  = ak.ones_like(tight_sf)


# Evaluate the btag eff
def btag_eff_eval(jets,wp,year):

    # Get the right process name for the given year and read in the histo
    pname_base = "TTZToLLNuNu_M_10"
    if year == "2016APV":
        pname = f"UL16APV_{pname_base}"
    elif year == "2016":
        pname = f"UL16_{pname_base}"
    elif year == "2017":
        pname = f"UL17_{pname_base}"
    elif year == "2018":
        pname = f"UL18_{pname_base}"
    elif year in ["2022","2022EE"]:
        pname_22_1   = "2022_TTZToLL_M_4to50"
        pname_22_2   = "2022_TTZToLL_M_50"
        pname_22EE_1 = "2022EE_TTZToLL_M_4to50"
        pname_22EE_2 = "2022EE_TTZToLL_M_50"
    elif year in ["2023","2023BPix"]:
        pname_23_1   = "2023_TTZToLL_M_4to50"
        pname_23_2   = "2023_TTZToLL_M_50"
        pname_23BPix_1 = "2023BPix_TTZToLL_M_4to50"
        pname_23BPix_2 = "2023BPix_TTZToLL_M_50"
    else:
        raise Exception(f"Not a known year: {year}")

    if year in ["2022","2022EE"]:
        pkl_file_path = ewkcoffea_path("data/btag_eff/btag_eff_2022_ttZ_srpresel.pkl.gz")
        histo = pickle.load(gzip.open(pkl_file_path))["ptabseta"]
        histo_proc = histo[{"process":pname_22_1}] + histo[{"process":pname_22_2}] + histo[{"process":pname_22EE_1}] + histo[{"process":pname_22EE_2}] #Adding 2022EE and 2022 together to fix low-stat bug
    elif year in ["2023","2023BPix"]:
        pkl_file_path = ewkcoffea_path("data/btag_eff/btag_eff_2023_ttZ_srpresel.pkl.gz")
        histo = pickle.load(gzip.open(pkl_file_path))["ptabseta"]
        histo_proc = histo[{"process":pname_23_1}] + histo[{"process":pname_23_2}] + histo[{"process":pname_23BPix_1}] + histo[{"process":pname_23BPix_2}] #Adding 2023 and 2023BPix together to be parallel with 2022
    else:
        pkl_file_path = ewkcoffea_path("data/btag_eff/btag_eff_ttZ_srpresel.pkl.gz")
        histo = pickle.load(gzip.open(pkl_file_path))["ptabseta"]
        histo_proc = histo[{"process":pname}]

    # Make sure wp is known
    if (wp != "L") and (wp != "M"):
        raise Exception(f"Not a known WP: {wp}")

    # Create lookup object and evaluate eff
    h_eff = histo_proc[{"tag":wp}] / histo_proc[{"tag":"all"}]
    vals = h_eff.values(flow=True)[1:,1:-1,:-1] # Pt (drop underflow), eta (drop under and over flow), flav (drop overflow, there is not underflow)
    h_eff_lookup = lookup_tools.dense_lookup.dense_lookup(vals, [ax.edges for ax in h_eff.axes])
    eff = h_eff_lookup(jets.pt,abs(jets.eta),jets.hadronFlavour)

    return eff

def run3_muons_sf_attach(muons,year,id_method,iso_method):

    # Get the right sf json for the given campaign
    if year == "2022EE":
        fname = ewkcoffea_path("data/run3_lep_sf/muon_sf/ScaleFactors_Muon_Z_ID_ISO_2022_EE_schemaV2.json")
    elif year == "2022":
        fname = ewkcoffea_path("data/run3_lep_sf/muon_sf/ScaleFactors_Muon_Z_ID_ISO_2022_schemaV2.json")
    elif year == "2023":
        fname = ewkcoffea_path("data/run3_lep_sf/muon_sf/ScaleFactors_Muon_Z_ID_ISO_2023_schemaV2.json")
    elif year == "2023BPix":
        fname = ewkcoffea_path("data/run3_lep_sf/muon_sf/ScaleFactors_Muon_Z_ID_ISO_2023_BPix_schemaV2.json")
    else:
        raise Exception("Trying to apply Run3 Muon SF where they shouldn't be!")

    # Flatten the input
    abseta_flat = ak.flatten(abs(muons.eta))
    pt_flat = ak.flatten(muons.pt)

    # For now, cap all pt at 199.9 and min of 15.0 TODO: Update so that we have all pt ranges
    pt_flat = ak.where(pt_flat>199.9,199.9,pt_flat)
    pt_flat = ak.where(pt_flat<15.0,15.0,pt_flat)

    # Evaluate the ID SF
    ceval = correctionlib.CorrectionSet.from_file(fname)
    sf_id_flat_nom = ceval[id_method].evaluate(abseta_flat,pt_flat,"nominal")
    sf_id_flat_hi  = ceval[id_method].evaluate(abseta_flat,pt_flat,"systup")
    sf_id_flat_lo  = ceval[id_method].evaluate(abseta_flat,pt_flat,"systdown")

    # Evaluate the Iso SF
    sf_iso_flat_nom = ceval[iso_method].evaluate(abseta_flat,pt_flat,"nominal")
    sf_iso_flat_hi  = ceval[iso_method].evaluate(abseta_flat,pt_flat,"systup")
    sf_iso_flat_lo = ceval[iso_method].evaluate(abseta_flat,pt_flat,"systdown")

    # Getting the overall SF (ID * Iso)
    sf_flat_nom = sf_id_flat_nom * sf_iso_flat_nom
    sf_flat_hi = sf_id_flat_hi * sf_iso_flat_hi
    sf_flat_lo = sf_id_flat_lo * sf_iso_flat_lo

    # Unflatten to match the original structure
    sf_nom = ak.unflatten(sf_flat_nom,ak.num(muons.pt))
    sf_hi = ak.unflatten(sf_flat_hi,ak.num(muons.pt))
    sf_lo = ak.unflatten(sf_flat_lo,ak.num(muons.pt))

    muons['sf_nom_muon'] = sf_nom
    muons['sf_hi_muon']  = sf_hi
    muons['sf_lo_muon']  = sf_lo
    muons['sf_nom_elec'] = ak.ones_like(sf_nom)
    muons['sf_hi_elec']  = ak.ones_like(sf_nom)
    muons['sf_lo_elec']  = ak.ones_like(sf_nom)

def run3_electrons_sf_attach(electrons,year,wp):

    # Get the right sf json for the given campaign
    if year == "2022EE":
        n_year = "2022Re-recoE+PromptFG" # key for accessing the 2022EE SFs
        fname = ewkcoffea_path("data/run3_lep_sf/electron_sf/2022EE_ele/electron.json")
    elif year == "2022":
        n_year = "2022Re-recoBCD" # key for accessing the 2022 SFs
        fname = ewkcoffea_path("data/run3_lep_sf/electron_sf/2022_ele/electron.json")
    elif year == "2023":
        n_year = "2023PromptC" # key for accessing the 2023 SFs
        fname = ewkcoffea_path("data/run3_lep_sf/electron_sf/2023_ele/electron.json")
    elif year == "2023BPix":
        n_year = "2023PromptD" # key for accessing the 2023BPix SFs
        fname = ewkcoffea_path("data/run3_lep_sf/electron_sf/2023BPix_ele/electron.json")
    else:
        raise Exception("Trying to apply run3 SF where they shouldn't be!")

    # Flatten the input
    phi_flat = ak.flatten(electrons.phi)
    eta_flat = ak.flatten(electrons.eta)
    pt_flat = ak.flatten(electrons.pt)
    ceval = correctionlib.CorrectionSet.from_file(fname)

    # Create three pT regions (based on the json pT regions for the working points)
    pt_flat_20 = ak.where(pt_flat >= 20.0,19.9,pt_flat)
    pt_flat_2075 = ak.where(pt_flat < 20.0, 20.0, ak.where(pt_flat >= 75.0, 74.9, pt_flat))
    pt_flat_75 = ak.where(pt_flat < 75.0, 75.0,pt_flat)

    #Get the Reco SF for all three region lists
    if year in ["2023","2023BPix"]:
        sf_flat_20 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","RecoBelow20",eta_flat,pt_flat_20,phi_flat)
        sf_flat_20_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","RecoBelow20",eta_flat,pt_flat_20,phi_flat)
        sf_flat_20_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","RecoBelow20",eta_flat,pt_flat_20,phi_flat)
        sf_flat_2075 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","Reco20to75",eta_flat,pt_flat_2075,phi_flat)
        sf_flat_2075_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","Reco20to75",eta_flat,pt_flat_2075,phi_flat)
        sf_flat_2075_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","Reco20to75",eta_flat,pt_flat_2075,phi_flat)
        sf_flat_75 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","RecoAbove75",eta_flat,pt_flat_75,phi_flat)
        sf_flat_75_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","RecoAbove75",eta_flat,pt_flat_75,phi_flat)
        sf_flat_75_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","RecoAbove75",eta_flat,pt_flat_75,phi_flat)
    elif year in ["2022","2022EE"]:
        sf_flat_20 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","RecoBelow20",eta_flat,pt_flat_20)
        sf_flat_20_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","RecoBelow20",eta_flat,pt_flat_20)
        sf_flat_20_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","RecoBelow20",eta_flat,pt_flat_20)
        sf_flat_2075 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","Reco20to75",eta_flat,pt_flat_2075)
        sf_flat_2075_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","Reco20to75",eta_flat,pt_flat_2075)
        sf_flat_2075_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","Reco20to75",eta_flat,pt_flat_2075)
        sf_flat_75 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","RecoAbove75",eta_flat,pt_flat_75)
        sf_flat_75_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","RecoAbove75",eta_flat,pt_flat_75)
        sf_flat_75_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","RecoAbove75",eta_flat,pt_flat_75)
    else:
        raise Exception("Unidentified Run 3 year. Exciting!")

    # Remove the unwanted Reco SF
    # We assigned values in the correct pT range in order to obtain the SF. We now need to remove the unwanted SF based on the original pt_flat
    reco_sf_flat_20 = ak.where(pt_flat >= 20.0, 0, sf_flat_20)
    reco_sf_flat_20_hi = ak.where(pt_flat >= 20.0, 0, sf_flat_20_hi)
    reco_sf_flat_20_lo = ak.where(pt_flat >= 20.0, 0, sf_flat_20_lo)

    reco_sf_flat_2075 = ak.where(pt_flat < 20.0, 0, ak.where(pt_flat >= 75.0, 0, sf_flat_2075))
    reco_sf_flat_2075_hi = ak.where(pt_flat < 20.0, 0, ak.where(pt_flat >= 75.0, 0, sf_flat_2075_hi))
    reco_sf_flat_2075_lo = ak.where(pt_flat < 20.0, 0, ak.where(pt_flat >= 75.0, 0, sf_flat_2075_lo))

    reco_sf_flat_75 = ak.where(pt_flat < 75.0, 0, sf_flat_75)
    reco_sf_flat_75_hi = ak.where(pt_flat < 75.0, 0, sf_flat_75_hi)
    reco_sf_flat_75_lo = ak.where(pt_flat < 75.0, 0, sf_flat_75_lo)

    #Add up the sf lists
    sf_reco = reco_sf_flat_20 + reco_sf_flat_2075 + reco_sf_flat_75
    sf_reco_hi = reco_sf_flat_20_hi + reco_sf_flat_2075_hi + reco_sf_flat_75_hi
    sf_reco_lo = reco_sf_flat_20_lo + reco_sf_flat_2075_lo + reco_sf_flat_75_lo

    # Evaluate the ID SF
    if "2023" in year:
        sf_id_flat = ceval["Electron-ID-SF"].evaluate(n_year,"sf",wp,eta_flat,pt_flat,phi_flat)
        hi_id_flat = ceval["Electron-ID-SF"].evaluate(n_year,"sfup",wp,eta_flat,pt_flat,phi_flat)
        lo_id_flat = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown",wp,eta_flat,pt_flat,phi_flat)
    else:
        sf_id_flat = ceval["Electron-ID-SF"].evaluate(n_year,"sf",wp,eta_flat,pt_flat)
        hi_id_flat = ceval["Electron-ID-SF"].evaluate(n_year,"sfup",wp,eta_flat,pt_flat)
        lo_id_flat = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown",wp,eta_flat,pt_flat)

    sf_return = sf_id_flat * sf_reco
    sf_return_hi = hi_id_flat * sf_reco_hi
    sf_return_lo = lo_id_flat * sf_reco_lo

    sf = ak.unflatten(sf_return,ak.num(electrons.pt))
    hi = ak.unflatten(sf_return_hi,ak.num(electrons.pt))
    lo = ak.unflatten(sf_return_lo,ak.num(electrons.pt))

    electrons['sf_nom_muon'] = ak.ones_like(sf)
    electrons['sf_hi_muon']  = ak.ones_like(sf)
    electrons['sf_lo_muon']  = ak.ones_like(sf)
    electrons['sf_nom_elec'] = sf
    electrons['sf_hi_elec']  = hi
    electrons['sf_lo_elec']  = lo

def run3_pu_attach(pileup,year,sys):

    # Get the right sf json for the given campaign
    if year == "2022EE":
        fname = ewkcoffea_path("data/run3_pu/pu_2022EE/puWeights.json")
    elif year == "2022":
        fname = ewkcoffea_path("data/run3_pu/pu_2022/puWeights.json")
    elif year == "2023":
        fname = ewkcoffea_path("data/run3_pu/pu_2023/puWeights.json")
    elif year == "2023BPix":
        fname = ewkcoffea_path("data/run3_pu/pu_2023BPix/puWeights.json")
    else:
        raise Exception("Trying to apply run3 SF where they shouldn't be!")

    # Evaluate the SF
    ceval = correctionlib.CorrectionSet.from_file(fname)
    if year == "2022EE":
        pu_corr = ceval["Collisions2022_359022_362760_eraEFG_GoldenJson"].evaluate(pileup.nTrueInt,"nominal")
        pu_corr_hi = ceval["Collisions2022_359022_362760_eraEFG_GoldenJson"].evaluate(pileup.nTrueInt,"up")
        pu_corr_lo = ceval["Collisions2022_359022_362760_eraEFG_GoldenJson"].evaluate(pileup.nTrueInt,"down")
    if year == "2022":
        pu_corr = ceval["Collisions2022_355100_357900_eraBCD_GoldenJson"].evaluate(pileup.nTrueInt,"nominal")
        pu_corr_hi = ceval["Collisions2022_355100_357900_eraBCD_GoldenJson"].evaluate(pileup.nTrueInt,"up")
        pu_corr_lo = ceval["Collisions2022_355100_357900_eraBCD_GoldenJson"].evaluate(pileup.nTrueInt,"down")
    if year == "2023":
        pu_corr = ceval["Collisions2023_366403_369802_eraBC_GoldenJson"].evaluate(pileup.nTrueInt,"nominal")
        pu_corr_hi = ceval["Collisions2023_366403_369802_eraBC_GoldenJson"].evaluate(pileup.nTrueInt,"up")
        pu_corr_lo = ceval["Collisions2023_366403_369802_eraBC_GoldenJson"].evaluate(pileup.nTrueInt,"down")
    if year == "2023BPix":
        pu_corr = ceval["Collisions2023_369803_370790_eraD_GoldenJson"].evaluate(pileup.nTrueInt,"nominal")
        pu_corr_hi = ceval["Collisions2023_369803_370790_eraD_GoldenJson"].evaluate(pileup.nTrueInt,"up")
        pu_corr_lo = ceval["Collisions2023_369803_370790_eraD_GoldenJson"].evaluate(pileup.nTrueInt,"down")
    if sys == "nominal":
        return pu_corr
    if sys == "hi":
        return pu_corr_hi
    if sys == "lo":
        return pu_corr_lo
    if sys not in ["nominal","hi","lo"]:
        raise Exception("ERROR: Not a recognized parameter.")

#def ApplyJetCorrections(year, isdata, era, corr_type):
#
#    if year not in clib_year_map.keys():
#        raise Exception(f"Error: Unknown year \"{year}\".")
#
#    jet_algo,jec_data_key,jec_mc_key,jer_key = get_jerc_keys(year,era,isdata)
#
#    jec_year = year
#    if year.startswith("2016"):
#        jec_year = "2016preVFP" if year == "2016APV" else "2016postVFP"
#    if int(year.replace("APV", "")) < 2020:
#        jec_year += "_UL"
#    json_path = topcoffea_path(f"data/POG/JME/{jec_year}/jet_jerc.json.gz")
#
#    jer_tag = jer_key
#    if isdata:
#        jec_tag = jec_data_key
#    else:
#        jec_tag = jec_mc_key
#
#    jec_names_clib = [
#        f"{jec_tag}_L1FastJet_{jet_algo}",
#        f"{jec_tag}_L2Relative_{jet_algo}",
#        f"{jec_tag}_L3Absolute_{jet_algo}",
#        f"{jec_tag}_L2L3Residual_{jet_algo}",
#    ]
#    jer_names_clib = [
#        f"{jer_tag}_MC_SF_{jet_algo}",
#        f"{jer_tag}_MC_PtResolution_{jet_algo}",
#    ]
#
#    if not isdata:
#        jec_names_clib += jer_names_clib
#
#    jec_types_clib = [
#        'FlavorQCD', 'FlavorPureBottom', 'FlavorPureQuark', 'FlavorPureGluon', 'FlavorPureCharm',
#        'Regrouped_BBEC1', 'Regrouped_Absolute', 'Regrouped_RelativeBal', 'RelativeSample'
#    ]
#    jec_regroup_clib = [f"Quad_Summer19UL{jec_tag}_MC_UncertaintySources_{jet_algo}_{jec_type}" for jec_type in jec_types_clib]
#
#    jec_names_clib.append(json_path)
#    jec_names_clib.append(False) # This boolean will be used to realize if the user wants to save the different level corrections or not
#    jec_stack = jec_names_clib
#
#    name_map = {}
#    name_map['JetPt'] = 'pt'
#    name_map['JetMass'] = 'mass'
#    name_map['JetEta'] = 'eta'
#    name_map['JetPhi'] = 'phi'
#    name_map['JetA'] = 'area'
#    name_map['ptGenJet'] = 'pt_gen'
#    name_map['ptRaw'] = 'pt_raw'
#    name_map['massRaw'] = 'mass_raw'
#    name_map['Rho'] = 'rho'
#    name_map['METpt'] = 'pt'
#    name_map['METphi'] = 'phi'
#    name_map['UnClusteredEnergyDeltaX'] = 'MetUnclustEnUpDeltaX'
#    name_map['UnClusteredEnergyDeltaY'] = 'MetUnclustEnUpDeltaY'
#    return CorrectedJetsFactory(name_map, jec_stack)
#
#
#def jerc_systematics(year,era,syst,jet_collection,isdata):
#
#    if syst == 'nominal':
#        jet_collection["pt"] = jet_collection.pt
#    elif syst.startswith('JER') and syst.endswith('Up'):
#        jet_collection["pt"] = jet_collection.JER.up.pt
#    elif syst.startswith('JER') and syst.endswith('Down'):
#        jet_collection["pt"] = jet_collection.JER.down.pt
#    else:
#        #Get the correct files
#        if year == "2016APV":
#            fname = ewkcoffea_path("data/wwz_jerc/2016APV_jerc/jet_jerc.json")
#        elif year == "2016":
#            fname = ewkcoffea_path("data/wwz_jerc/2016_jerc/jet_jerc.json")
#        elif year == "2017":
#            fname = ewkcoffea_path("data/wwz_jerc/2017_jerc/jet_jerc.json")
#        elif year == "2018":
#            fname = ewkcoffea_path("data/wwz_jerc/2018_jerc/jet_jerc.json")
#        elif year == "2022":
#            fname = ewkcoffea_path("data/wwz_jerc/2022_jerc/jet_jerc.json")
#        elif year == "2022EE":
#            fname = ewkcoffea_path("data/wwz_jerc/2022EE_jerc/jet_jerc.json")
#        elif year == "2023":
#            fname = ewkcoffea_path("data/wwz_jerc/2023_jerc/jet_jerc.json")
#        elif year == "2023BPix":
#            fname = ewkcoffea_path("data/wwz_jerc/2023BPix_jerc/jet_jerc.json")
#        else:
#            raise Exception("Unrecognized year for jerc systematics. Exciting!")
#
#        # Initiliaze files we need
#        jerc_ceval = correctionlib.CorrectionSet.from_file(fname)
#
#        #Inputs we will need
#        jet_pt = jet_collection.pt
#        jet_eta = jet_collection.eta
#
#        #We need to flatten the arrays
#        jet_pt_flat = ak.flatten(jet_pt)
#        jet_eta_flat = ak.flatten(jet_eta)
#        
#        #Get the appropriate key and jet type
#        jet_algo,jec_data_key,jec_mc_key,jer_key = get_jerc_keys(year,era,isdata)
#        #Figure out the syst and up/down
#        syst_type = syst.split("_")[0]
#        updown = syst.split("_")[-1]
#        if updown.endswith("Up"):
#            updown = "Up"
#        elif updown.endswith("Down"):
#            updown = "Down"
#        else:
#            raise Exception("JEC uncertainty should be up or down!")
#
#        #Get the factor for the up/down variation for JEC systematics
#        if syst.startswith("JEC"):
#            factor = jet_pt_flat * jerc_ceval[f"{jec_mc_key}_Total_{jet_algo}"].evaluate(jet_eta_flat,jet_pt_flat)
#        else:
#            factor = jet_pt_flat * jerc_ceval[f"{jec_mc_key}_{syst_type}_{jet_algo}"].evaluate(jet_eta_flat,jet_pt_flat)
#
#        if updown == "Up":
#            jet_pt_final = jet_pt_flat + factor
#        elif updown == "Down":
#            jet_pt_final = jet_pt_flat - factor
#
#        jet_collection["pt"] = ak.unflatten(jet_pt_final,ak.num(jet_collection.pt))

def jerc_corrections(year,era,isdata,correction,jet_collection,events):

    #Get the correct files
    smear_fname = ewkcoffea_path("data/wwz_jerc/jer_smear.json")
    if year == "2016APV":
        fname = ewkcoffea_path("data/wwz_jerc/2016APV_jerc/jet_jerc.json")
    elif year == "2016":
        fname = ewkcoffea_path("data/wwz_jerc/2016_jerc/jet_jerc.json")
    elif year == "2017":
        fname = ewkcoffea_path("data/wwz_jerc/2017_jerc/jet_jerc.json")
    elif year == "2018":
        fname = ewkcoffea_path("data/wwz_jerc/2018_jerc/jet_jerc.json")
    elif year == "2022":
        fname = ewkcoffea_path("data/wwz_jerc/2022_jerc/jet_jerc.json")
    elif year == "2022EE":
        fname = ewkcoffea_path("data/wwz_jerc/2022EE_jerc/jet_jerc.json")
    elif year == "2023":
        fname = ewkcoffea_path("data/wwz_jerc/2023_jerc/jet_jerc.json")
    elif year == "2023BPix":
        fname = ewkcoffea_path("data/wwz_jerc/2023BPix_jerc/jet_jerc.json")
    else:
        raise Exception("Unrecognized year for jerc_corrections. Exciting!")

    # Initiliaze files we need
    jerc_ceval = correctionlib.CorrectionSet.from_file(fname)
    smear_ceval = correctionlib.CorrectionSet.from_file(smear_fname)


    jet_raw_pt = (1 - jet_collection.rawFactor)*jet_collection.pt
    jet_eta = jet_collection.eta
    jet_area = jet_collection.area
    if year in ["2022","2022EE","2023","2023BPix"]:
        rho = events.Rho.fixedGridRhoFastjetAll
    elif year in ["2016APV","2016","2017","2018"]:
        rho = events.fixedGridRhoFastjetAll
    event = events.event
    jet_genpt = jet_collection.pt_gen
    jet_phi = jet_collection.phi

    #We need to broadcast Rho  and Event Number to have the same shape as jet pt
    fixed_rho, _ = ak.broadcast_arrays(rho, jet_raw_pt)
    fixed_event, _ = ak.broadcast_arrays(event, jet_raw_pt)

    #We need to flatten the arrays
    jet_rho_flat = ak.flatten(fixed_rho)
    jet_raw_pt_flat = ak.flatten(jet_raw_pt)
    jet_eta_flat = ak.flatten(jet_eta)
    jet_area_flat = ak.flatten(jet_area)
    jet_event_flat = ak.flatten(fixed_event)
    jet_genpt_flat = ak.flatten(jet_genpt)
    jet_phi_flat = ak.flatten(jet_phi)
    
    #Get the appropriate key and jet type
    jet_algo,jec_key,jer_key = get_jerc_keys(year,era,isdata)

    #JEC
    #L1FastJet
    l1_fj = jerc_ceval[f"{jec_key}_L1FastJet_{jet_algo}"].evaluate(jet_area_flat,jet_eta_flat,jet_raw_pt_flat,jet_rho_flat)
    pt_l1 = l1_fj * jet_raw_pt_flat
    #L2Relative
    if year == "2023BPix":
        l2_rel = jerc_ceval[f"{jec_key}_L2Relative_{jet_algo}"].evaluate(jet_eta_flat,jet_phi_flat,pt_l1)
        pt_l2 = l2_rel * pt_l1
    else:
        l2_rel = jerc_ceval[f"{jec_key}_L2Relative_{jet_algo}"].evaluate(jet_eta_flat,pt_l1)
        pt_l2 = l2_rel * pt_l1
    #L3Absolute 
    l3_abs = jerc_ceval[f"{jec_key}_L3Absolute_{jet_algo}"].evaluate(jet_eta_flat,pt_l2)
    pt_l3 = l3_abs * pt_l2
    #L2L3Residual (Note that this is not needed for MC but will always return 1)
    l2l3_res = jerc_ceval[f"{jec_key}_L2L3Residual_{jet_algo}"].evaluate(jet_eta_flat,pt_l3)
    pt_l4 = l2l3_res * pt_l3

    jec_pt_final = pt_l4

    #Data does not have uncertaintity or JER
    if isdata:
        jet_collection["pt"] = ak.unflatten(jec_pt_final,ak.num(jet_collection.pt))
    else:
        #JER smearing time
        if (correction.startswith("JER") and correction.endswith("Up")):
            up_down_nom = "up"
        elif (correction.startswith("JER") and correction.endswith("Down")):
            up_down_nom = "down"
        else:
            up_down_nom = "nom"

        #Get the JER and JER SF
        jer = jerc_ceval[f"{jer_key}_MC_PtResolution_{jet_algo}"].evaluate(jet_eta_flat,jec_pt_final,jet_rho_flat)

        if year in ["2016","2016APV","2017","2018"]:
            jer_sf = jerc_ceval[f"{jer_key}_MC_ScaleFactor_{jet_algo}"].evaluate(jet_eta_flat,up_down_nom)

        elif year in ["2022","2022EE","2023","2023BPix"]:
            jer_sf = jerc_ceval[f"{jer_key}_MC_ScaleFactor_{jet_algo}"].evaluate(jet_eta_flat,jec_pt_final,up_down_nom)

        #JER Smearing
        smear_factor = smear_ceval["JERSmear"].evaluate(jec_pt_final,jet_eta_flat,jet_genpt_flat,jet_rho_flat,jet_event_flat,jer,jer_sf)

        #Get new nominal pt
        pt_jerc_final = smear_factor*jec_pt_final

        #If nominal or JER correction, we stop here
        if ((correction == 'nominal') or correction.startswith('JER')):
            jet_collection["pt"] = ak.unflatten(pt_jerc_final,ak.num(jet_collection.pt))

        else:

            #Figure out the syst and up/down
            syst = correction.split("_")[0]
            updown = correction.split("_")[-1]
            if updown.endswith("Up"):
                updown = "Up"
            elif updown.endswith("Down"):
                updown = "Down"
            else:
                raise Exception("JEC uncertainty should be up or down!")

            if syst == "JEC":
                factor = pt_jerc_final * jerc_ceval[f"{jec_key}_Total_{jet_algo}"].evaluate(jet_eta_flat,pt_jerc_final)
            else:
                factor = pt_jerc_final * jerc_ceval[f"{jec_key}_{syst}_{jet_algo}"].evaluate(jet_eta_flat,pt_jerc_final)

            if updown == "Up":
                pt_jerc_final = pt_jerc_final + factor
            elif updown == "Down":
                pt_jerc_final = pt_jerc_final - factor

            jet_collection["pt"] = ak.unflatten(pt_jerc_final,ak.num(jet_collection.pt))
