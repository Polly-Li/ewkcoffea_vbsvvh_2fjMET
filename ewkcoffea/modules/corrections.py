import numpy as np
import pickle
import gzip
import awkward as ak

import correctionlib
from coffea import lookup_tools

from ewkcoffea.modules.paths import ewkcoffea_path


extLepSF = lookup_tools.extractor()

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
    elif year == "2022":
        pname = f"UL18_{pname_base}" #TODO Update with 2022 efficiency when available
    elif year == "2022EE":
        pname = f"UL18_{pname_base}" #TODO Update with 2022EE efficiency when available
    elif year == "2023":
        pname = f"UL18_{pname_base}" #TODO Update with 2023 efficiency when available
    elif year == "2023BPix":
        pname = f"UL18_{pname_base}" #TODO Update with 2023BPix efficiency when available
    else:
        raise Exception(f"Not a known year: {year}")

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
    if "2023" in year:
        sf_flat_20 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","RecoBelow20",eta_flat,pt_flat_20,phi_flat)
        sf_flat_20_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","RecoBelow20",eta_flat,pt_flat_20,phi_flat)
        sf_flat_20_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","RecoBelow20",eta_flat,pt_flat_20,phi_flat)
        sf_flat_2075 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","Reco20to75",eta_flat,pt_flat_2075,phi_flat)
        sf_flat_2075_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","Reco20to75",eta_flat,pt_flat_2075,phi_flat)
        sf_flat_2075_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","Reco20to75",eta_flat,pt_flat_2075,phi_flat)
        sf_flat_75 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","RecoAbove75",eta_flat,pt_flat_75,phi_flat)
        sf_flat_75_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","RecoAbove75",eta_flat,pt_flat_75,phi_flat)
        sf_flat_75_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","RecoAbove75",eta_flat,pt_flat_75,phi_flat)
    else:
        sf_flat_20 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","RecoBelow20",eta_flat,pt_flat_20)
        sf_flat_20_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","RecoBelow20",eta_flat,pt_flat_20)
        sf_flat_20_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","RecoBelow20",eta_flat,pt_flat_20)
        sf_flat_2075 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","Reco20to75",eta_flat,pt_flat_2075)
        sf_flat_2075_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","Reco20to75",eta_flat,pt_flat_2075)
        sf_flat_2075_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","Reco20to75",eta_flat,pt_flat_2075)
        sf_flat_75 = ceval["Electron-ID-SF"].evaluate(n_year,"sf","RecoAbove75",eta_flat,pt_flat_75)
        sf_flat_75_hi = ceval["Electron-ID-SF"].evaluate(n_year,"sfup","RecoAbove75",eta_flat,pt_flat_75)
        sf_flat_75_lo = ceval["Electron-ID-SF"].evaluate(n_year,"sfdown","RecoAbove75",eta_flat,pt_flat_75)

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
