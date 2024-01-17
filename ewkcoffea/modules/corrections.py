import numpy as np
import pickle
import gzip
import awkward as ak

from coffea import lookup_tools

from ewkcoffea.modules.paths import ewkcoffea_path
from topcoffea.modules.paths import topcoffea_path


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


