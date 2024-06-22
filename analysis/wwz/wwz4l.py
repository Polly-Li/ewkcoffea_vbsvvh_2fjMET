#!/usr/bin/env python
#import sys
import coffea
import numpy as np
import awkward as ak
import copy
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import processor
import hist
from hist import axis
from coffea.analysis_tools import PackedSelection
from coffea.lumi_tools import LumiMask

from topcoffea.modules.paths import topcoffea_path
import topcoffea.modules.event_selection as es_tc
import topcoffea.modules.object_selection as os_tc
import topcoffea.modules.corrections as cor_tc

from ewkcoffea.modules.paths import ewkcoffea_path as ewkcoffea_path
import ewkcoffea.modules.selection_wwz as es_ec
import ewkcoffea.modules.objects_wwz as os_ec
import ewkcoffea.modules.corrections as cor_ec

from topcoffea.modules.get_param_from_jsons import GetParam
get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))


# Small helper function for creating the list of systematics
# Append "Up" and "Down" to all base strings in a given syst list
def append_up_down_to_sys_base(sys_lst_in):
    sys_lst_out = []
    for s in sys_lst_in:
        sys_lst_out.append(f"{s}Up")
        sys_lst_out.append(f"{s}Down")
    return sys_lst_out

# Takes a list of variable names, and a dictionary mapping those to the ak arrays
# Loops through each variable, and replaces None values with whatever you specify with none_val
# I.e. this is really just a wrapper around ak.fill_none for a list of variables
# Currently just used to get rid of Nones in the variable list we pass to sr bdt evaluation
def fill_none_in_list(var_names,var_names_vals_map,none_val):
    out_lst = []
    for var_name in var_names:
        vals = var_names_vals_map[var_name]
        out_lst.append(ak.fill_none(vals,none_val))
    return out_lst

class AnalysisProcessor(processor.ProcessorABC):

    def __init__(self, samples, wc_names_lst=[], hist_lst=None, ecut_threshold=None, do_errors=False, do_systematics=False, split_by_lepton_flavor=False, skip_signal_regions=False, skip_control_regions=False, muonSyst='nominal', dtype=np.float32, siphon_bdt_data=False):

        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype

        # Create the dense axes for the histograms
        self._dense_axes_dict = {
            "mt2"   : axis.Regular(180, 0, 250, name="mt2",  label="mt2"),
            "met"   : axis.Regular(180, 0, 300, name="met",  label="met"),
            "metphi": axis.Regular(180, -3.1416, 3.1416, name="metphi", label="met phi"),
            "ptl4"  : axis.Regular(180, 0, 500, name="ptl4", label="ptl4"),
            "scalarptsum_lep" : axis.Regular(180, 0, 600, name="scalarptsum_lep", label="S_T"),
            "scalarptsum_lepmet" : axis.Regular(180, 0, 600, name="scalarptsum_lepmet", label="S_T + metpt"),
            "scalarptsum_lepmetjet" : axis.Regular(180, 0, 1100, name="scalarptsum_lepmetjet", label="S_T + metpt + H_T"),
            "scalarptsum_jet" : axis.Regular(180, 0, 500, name="scalarptsum_jet", label="H_T"),
            "mll_01": axis.Regular(180, 0, 200, name="mll_01",  label="mll_l0_l1"),
            "mllll": axis.Regular(180, 0, 600, name="mllll",  label="mllll"),
            "l0pt"  : axis.Regular(180, 0, 500, name="l0pt", label="l0pt"),
            "j0pt"  : axis.Regular(180, 0, 500, name="j0pt", label="j0pt"),

            "w_lep0_pt"  : axis.Regular(180, 0, 300, name="w_lep0_pt", label="Leading W lep pt"),
            "w_lep1_pt"  : axis.Regular(180, 0, 300, name="w_lep1_pt", label="Subleading W lep pt"),
            "z_lep0_pt"  : axis.Regular(180, 0, 300, name="z_lep0_pt", label="Leading Z lep pt"),
            "z_lep1_pt"  : axis.Regular(180, 0, 300, name="z_lep1_pt", label="Subleading Z lep pt"),
            "w_lep0_eta" : axis.Regular(180, -3, 3, name="w_lep0_eta", label="Leading W lep eta"),
            "w_lep1_eta" : axis.Regular(180, -3, 3, name="w_lep1_eta", label="Subleading W lep eta"),
            "z_lep0_eta" : axis.Regular(180, -3, 3, name="z_lep0_eta", label="Leading Z lep eta"),
            "z_lep1_eta" : axis.Regular(180, -3, 3, name="z_lep1_eta", label="Subleading Z lep eta"),
            "w_lep0_phi" : axis.Regular(180, -3.1416, 3.1416, name="w_lep0_phi", label="Leading W lep phi"),
            "w_lep1_phi" : axis.Regular(180, -3.1416, 3.1416, name="w_lep1_phi", label="Subleading W lep phi"),
            "z_lep0_phi" : axis.Regular(180, -3.1416, 3.1416, name="z_lep0_phi", label="Leading Z lep phi"),
            "z_lep1_phi" : axis.Regular(180, -3.1416, 3.1416, name="z_lep1_phi", label="Subleading Z lep phi"),
            "mll_wl0_wl1" : axis.Regular(180, 0, 350, name="mll_wl0_wl1", label="mll(W lep0, W lep1)"),
            "mll_zl0_zl1" : axis.Regular(180, 0, 200, name="mll_zl0_zl1", label="mll(Z lep0, Z lep1)"),

            "pt_zl0_zl1" : axis.Regular(180, 0, 300, name="pt_zl0_zl1", label="pt(Zl0 + Zl1)"),
            "pt_wl0_wl1" : axis.Regular(180, 0, 300, name="pt_wl0_wl1", label="pt(Wl0 + Wl1)"),
            "dr_zl0_zl1" : axis.Regular(180, 0, 5, name="dr_zl0_zl1", label="dr(Zl0,Zl1)"),
            "dr_wl0_wl1" : axis.Regular(180, 0, 5, name="dr_wl0_wl1", label="dr(Wl0,Wl1)"),
            "dr_wleps_zleps" : axis.Regular(180, 0, 5, name="dr_wleps_zleps", label="dr((Wl0+Wl1),(Zl0,Zl1))"),
            "dr_wl0_j_min" : axis.Regular(180, 0, 5, name="dr_wl0_j_min",  label="min dr(Wl0,j)"),
            "dr_wl1_j_min" : axis.Regular(180, 0, 5, name="dr_wl1_j_min",  label="min dr(Wl1,j)"),

            "mt_4l_met"   : axis.Regular(180, 0, 500, name="mt_4l_met", label="mT of 4l system and met"),
            "mt_wleps_met": axis.Regular(180, 0, 300, name="mt_wleps_met", label="mT of W leptons system and met"),
            "mt_wl0_met"  : axis.Regular(180, 0, 300, name="mt_wl0_met", label="mT of W lep0 and met"),
            "mt_wl1_met"  : axis.Regular(180, 0, 300, name="mt_wl1_met", label="mT of W lep1 and met"),

            "absdphi_zl0_zl1": axis.Regular(180, 0, 3.1416, name="absdphi_zl0_zl1", label="abs dphi(Zl0,Zl1)"),
            "absdphi_wl0_wl1": axis.Regular(180, 0, 3.1416, name="absdphi_wl0_wl1", label="abs dphi(Wl0,Wl1)"),
            "absdphi_z_ww"   : axis.Regular(180, 0, 3.1416, name="absdphi_z_ww", label="abs dphi((Zl0+Zl1),(Wl0+Wl1+met))"),
            "absdphi_4l_met" : axis.Regular(180, 0, 3.1416, name="absdphi_4l_met", label="abs dphi((Zl0+Zl1+Wl0+Wl1),met)"),
            "absdphi_zleps_met" : axis.Regular(180, 0, 3.1416, name="absdphi_zleps_met", label="absdphi((Zl0+Zl1),met)"),
            "absdphi_wleps_met" : axis.Regular(180, 0, 3.1416, name="absdphi_wleps_met", label="abs dphi((Wl0+Wl1),met)"),
            "absdphi_wl0_met": axis.Regular(180, 0, 3.1416, name="absdphi_wl0_met", label="abs dphi(Wl0,met)"),
            "absdphi_wl1_met": axis.Regular(180, 0, 3.1416, name="absdphi_wl1_met", label="abs dphi(Wl1,met)"),

            "absdphi_min_afas" : axis.Regular(180, 0, 3.1416, name="absdphi_min_afas",  label="min(abs(delta phi of all pairs))"),
            "absdphi_min_afos" : axis.Regular(180, 0, 3.1416, name="absdphi_min_afos",  label="min(abs(delta phi of OS pairs))"),
            "absdphi_min_sfos" : axis.Regular(180, 0, 3.1416, name="absdphi_min_sfos",  label="min(abs(delta phi of SFOS pairs))"),
            "mll_min_afas" : axis.Regular(180, 0, 150, name="mll_min_afas",  label="min mll of all pairs"),
            "mll_min_afos" : axis.Regular(180, 0, 150, name="mll_min_afos",  label="min mll of OF pairs"),
            "mll_min_sfos" : axis.Regular(180, 0, 150, name="mll_min_sfos",  label="min mll of SFOF pairs"),

            "cos_helicity_x" : axis.Regular(180, 0, 1, name="cos_helicity_x",  label="cos_helicity_x"),

            "mlb_min" : axis.Regular(180, 0, 300, name="mlb_min",  label="min mass(b+l)"),
            "mlb_max" : axis.Regular(180, 0, 1000, name="mlb_max",  label="max mass(b+l)"),

            "njets"   : axis.Regular(8, 0, 8, name="njets",   label="Jet multiplicity"),
            "nleps"   : axis.Regular(5, 0, 5, name="nleps",   label="Lep multiplicity"),
            "nbtagsl" : axis.Regular(6, 0, 6, name="nbtagsl", label="Loose btag multiplicity"),
            "nbtagsm" : axis.Regular(4, 0, 4, name="nbtagsm", label="Medium btag multiplicity"),

            "njets_counts"   : axis.Regular(30, 0, 30, name="njets_counts",   label="Jet multiplicity counts"),
            "nleps_counts"   : axis.Regular(30, 0, 30, name="nleps_counts",   label="Lep multiplicity counts"),
            "nbtagsl_counts" : axis.Regular(30, 0, 30, name="nbtagsl_counts", label="Loose btag multiplicity counts"),

            "bdt_of_wwz": axis.Regular(180, 0, 1, name="bdt_of_wwz", label="Score bdt_of_wwz"),
            "bdt_sf_wwz": axis.Regular(180, 0, 1, name="bdt_sf_wwz", label="Score bdt_sf_wwz"),
            "bdt_of_zh" : axis.Regular(180, 0, 1, name="bdt_of_zh", label="Score bdt_of_zh"),
            "bdt_sf_zh" : axis.Regular(180, 0, 1, name="bdt_sf_zh", label="Score bdt_sf_zh"),
            "bdt_of_bkg" : axis.Regular(180, 0, 1, name="bdt_of_bkg", label="Score bdt_of_bkg"),
            "bdt_sf_bkg" : axis.Regular(180, 0, 1, name="bdt_sf_bkg", label="Score bdt_sf_bkg"),
            "bdt_of_wwz_m_zh" : axis.Regular(180, -1, 1, name="bdt_of_wwz_m_zh", label="Score bdt_of_wwz - bdt_of_zh"),
            "bdt_sf_wwz_m_zh" : axis.Regular(180, -1, 1, name="bdt_sf_wwz_m_zh", label="Score bdt_sf_wwz - bdt_sf_zh"),
            "bdt_of_bin" : axis.Regular(8, 0, 8, name="bdt_of_bin", label="Binned bdt_of"),
            "bdt_sf_bin" : axis.Regular(8, 0, 8, name="bdt_sf_bin", label="Binned bdt_sf"),

        }

        # Add histograms to dictionary that will be passed on to dict_accumulator
        dout = {}
        for dense_axis_name in self._dense_axes_dict.keys():
            dout[dense_axis_name] = hist.Hist(
                hist.axis.StrCategory([], growth=True, name="process", label="process"),
                hist.axis.StrCategory([], growth=True, name="category", label="category"),
                hist.axis.StrCategory([], growth=True, name="systematic", label="systematic"),
                self._dense_axes_dict[dense_axis_name],
                storage="weight", # Keeps track of sumw2
                name="Counts",
            )

        if siphon_bdt_data:
            # Adding list accumulators for BDT output variables and weights
            list_output_names = [
                "bdt_of_wwz_list",
                "bdt_of_zh_list",
                "bdt_of_bkg_list",
                "bdt_of_proc_list",
                "bdt_of_wgt_list",
                "bdt_of_evt_list",
                "bdt_sf_wwz_list",
                "bdt_sf_zh_list",
                "bdt_sf_bkg_list",
                "bdt_sf_proc_list",
                "bdt_sf_wgt_list",
                "bdt_sf_evt_list",
            ]
            for list_output_name in list_output_names:
                dout[list_output_name] = processor.list_accumulator([])

            # Also adding BDT input variable list accumulators
            for of_bdt_var_name in get_ec_param("of_bdt_var_lst"):
                dout["of_bdt_" + of_bdt_var_name] = processor.list_accumulator([])
            for sf_bdt_var_name in get_ec_param("sf_bdt_var_lst"):
                dout["sf_bdt_" + sf_bdt_var_name] = processor.list_accumulator([])

        # Set the accumulator
        self._accumulator = processor.dict_accumulator(dout)

        # Set the list of hists to fill
        if hist_lst is None:
            # If the hist list is none, assume we want to fill all hists
            self._hist_lst = list(self._accumulator.keys())
        else:
            # Otherwise, just fill the specified subset of hists
            for hist_to_include in hist_lst:
                if hist_to_include not in self._accumulator.keys():
                    raise Exception(f"Error: Cannot specify hist \"{hist_to_include}\", it is not defined in the processor.")
            self._hist_lst = hist_lst # Which hists to fill

        # Set the energy threshold to cut on
        self._ecut_threshold = ecut_threshold

        # Set the booleans
        self._do_errors = do_errors # Whether to calculate and store the w**2 coefficients
        self._do_systematics = do_systematics # Whether to process systematic samples
        self._split_by_lepton_flavor = split_by_lepton_flavor # Whether to keep track of lepton flavors individually
        self._skip_signal_regions = skip_signal_regions # Whether to skip the SR categories
        self._skip_control_regions = skip_control_regions # Whether to skip the CR categories
        self._siphon_bdt_data = siphon_bdt_data # Whether to write out bdt data or not

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    # Main function: run on a given dataset
    def process(self, events):

        # Dataset parameters
        json_name = events.metadata["dataset"]

        isData             = self._samples[json_name]["isData"]
        histAxisName       = self._samples[json_name]["histAxisName"]
        year               = self._samples[json_name]["year"]
        xsec               = self._samples[json_name]["xsec"]
        sow                = self._samples[json_name]["nSumOfWeights"]

        # Set a flag if this is a 2022 year
        is2022 = year in ["2022","2022EE"]

        # If this is a 2022 sample, get the era info
        if isData and is2022:
            era = self._samples[json_name]["era"]
        else:
            era = None


        # Get up down weights from input dict
        if (self._do_systematics and not isData):
            lhe_sow = self._samples[json_name]["nSumOfLheWeights"]
            # This assumes we have an NLO xsec, so for these systs we will have e.g. xsec_NLO*(N_pass_up/N_gen_up)
            # Thus these systs should only affect acceptance and effeciency and shape
            # The uncty on xsec comes from NLO and is applied as a rate uncty in the text datacard
            if lhe_sow == []:
                sow_renormDown     = sow
                sow_factDown       = sow
                sow_factUp         = sow
                sow_renormUp       = sow
            elif len(lhe_sow) == 9:
                sow_renormDown     = lhe_sow[1]
                sow_factDown       = lhe_sow[3]
                sow_factUp         = lhe_sow[5]
                sow_renormUp       = lhe_sow[7]
            elif len(lhe_sow) == 8:
                sow_renormDown     = lhe_sow[1]
                sow_factDown       = lhe_sow[3]
                sow_factUp         = lhe_sow[4]
                sow_renormUp       = lhe_sow[6]
            else: raise Exception("ERROR: Unknown LHE weights length {len(lhe_sow)}")
        else:
            sow_renormUp       = -1
            sow_renormDown     = -1
            sow_factUp         = -1
            sow_factDown       = -1


        # Get the dataset name (used for duplicate removal) and check to make sure it is an expected name
        # Get name for MC cases too, since "dataset" is passed to overlap removal function in all cases (though it's not actually used in the MC case)
        dataset = json_name.split('_')[0]
        if isData:
            datasets = ["SingleMuon", "SingleElectron", "EGamma", "MuonEG", "DoubleMuon", "DoubleElectron", "DoubleEG","Muon"]
            if dataset not in datasets:
                raise Exception("ERROR: Unexpected dataset name for data file.")

        # Initialize objects
        #met  = events.MET
        met  = events.PuppiMET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet

        # An array of lenght events that is just 1 for each event
        # Probably there's a better way to do this, but we use this method elsewhere so I guess why not..
        events.nom = ak.ones_like(met.pt)

        # Get the lumi mask for data
        if year == "2016" or year == "2016APV":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")
        elif year == "2017":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt")
        elif year == "2018":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt")
        elif year == "2022" or year == "2022EE":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_Collisions2022_355100_362760_Golden.txt")
        else:
            raise ValueError(f"Error: Unknown year \"{year}\".")
        lumi_mask = LumiMask(golden_json_path)(events.run,events.luminosityBlock)

        ################### Lepton selection ####################

        # Do the object selection for the WWZ eleectrons
        ele_presl_mask = os_ec.is_presel_wwz_ele(ele,is2022)
        if not is2022:
            ele["topmva"] = os_ec.get_topmva_score_ele(events, year)
            ele["is_tight_lep_for_wwz"] = ((ele.topmva > get_tc_param("topmva_wp_t_e")) & ele_presl_mask)
        else:
            ele["is_tight_lep_for_wwz"] = (ele_presl_mask)

        # Do the object selection for the WWZ muons
        mu_presl_mask = os_ec.is_presel_wwz_mu(mu, is2022)
        if not is2022:
            mu["topmva"] = os_ec.get_topmva_score_mu(events, year)
            mu["is_tight_lep_for_wwz"] = ((mu.topmva > get_tc_param("topmva_wp_t_m")) & mu_presl_mask)
        else:
            mu["is_tight_lep_for_wwz"] = (mu_presl_mask)

        # Get tight leptons for WWZ selection
        ele_wwz_t = ele[ele.is_tight_lep_for_wwz]
        mu_wwz_t = mu[mu.is_tight_lep_for_wwz]

        # Attach the lepton SFs to the electron and muons collections
        if is2022:
            cor_ec.run3_muons_sf_attach(mu_wwz_t,year,"NUM_MediumID_DEN_TrackerMuons","NUM_LoosePFIso_DEN_MediumID")
            cor_ec.run3_electrons_sf_attach(ele_wwz_t,year,"wp90iso")
        else:
            cor_ec.AttachElectronSF(ele_wwz_t,year=year)
            cor_ec.AttachMuonSF(mu_wwz_t,year=year)

        l_wwz_t = ak.with_name(ak.concatenate([ele_wwz_t,mu_wwz_t],axis=1),'PtEtaPhiMCandidate')
        l_wwz_t = l_wwz_t[ak.argsort(l_wwz_t.pt, axis=-1,ascending=False)] # Sort by pt


        # For WWZ: Compute pair invariant masses
        llpairs_wwz = ak.combinations(l_wwz_t, 2, fields=["l0","l1"])
        os_pairs_mask = (llpairs_wwz.l0.pdgId*llpairs_wwz.l1.pdgId < 0)   # Maks for opposite-sign pairs
        sfos_pairs_mask = (llpairs_wwz.l0.pdgId == -llpairs_wwz.l1.pdgId) # Mask for same-flavor-opposite-sign pairs
        ll_absdphi_pairs = abs(llpairs_wwz.l0.delta_phi(llpairs_wwz.l1))
        ll_mass_pairs = (llpairs_wwz.l0+llpairs_wwz.l1).mass            # The mll for each ll pair
        absdphi_min_afas = ak.min(ll_absdphi_pairs,axis=-1)
        absdphi_min_afos = ak.min(ll_absdphi_pairs[os_pairs_mask],axis=-1)
        absdphi_min_sfos = ak.min(ll_absdphi_pairs[sfos_pairs_mask],axis=-1)
        mll_min_afas = ak.min(ll_mass_pairs,axis=-1)
        mll_min_afos = ak.min(ll_mass_pairs[os_pairs_mask],axis=-1)
        mll_min_sfos = ak.min(ll_mass_pairs[sfos_pairs_mask],axis=-1)
        events["min_mll_afos"] = mll_min_afos # Attach this one to the event info since we need it for selection

        # For WWZ
        l_wwz_t_padded = ak.pad_none(l_wwz_t, 4)
        l0 = l_wwz_t_padded[:,0]
        l1 = l_wwz_t_padded[:,1]
        l2 = l_wwz_t_padded[:,2]
        l3 = l_wwz_t_padded[:,3]

        nleps = ak.num(l_wwz_t)

        # Put njets and l_fo_conept_sorted into events and get 4l event selection mask
        events["l_wwz_t"] = l_wwz_t
        es_ec.add4lmask_wwz(events, year, isData, histAxisName, is2022)


        ######### Normalization and weights ###########


        # These weights can go outside of the outside sys loop since they do not depend on pt of mu or jets
        # We only calculate these values if not isData
        # Note: add() will generally modify up/down weights, so if these are needed for any reason after this point, we should instead pass copies to add()
        # Note: Here we will to the weights object the SFs that do not depend on any of the forthcoming loops
        weights_obj_base = coffea.analysis_tools.Weights(len(events),storeIndividual=True)
        if not isData:
            genw = events["genWeight"]

            # If it's an EFT sample, just take SM piece
            sm_wgt = 1.0
            eft_coeffs = ak.to_numpy(events["EFTfitCoefficients"]) if hasattr(events, "EFTfitCoefficients") else None
            if eft_coeffs is not None:
                sm_wgt = eft_coeffs[:,0]

            # Normalize by (xsec/sow)*genw where genw is 1 for EFT samples
            # Note that for theory systs, will need to multiply by sow/sow_wgtUP to get (xsec/sow_wgtUp)*genw and same for Down
            lumi = 1000.0*get_tc_param(f"lumi_{year}")
            weights_obj_base.add("norm",(xsec/sow)*genw*lumi*sm_wgt)


            # Scale weights
            cor_tc.AttachPSWeights(events)
            cor_tc.AttachScaleWeights(events)
            # FSR/ISR weights
            # For now only consider variations in the numerator
            weights_obj_base.add('ISR', events.nom, events.ISRUp, events.ISRDown)
            weights_obj_base.add('FSR', events.nom, events.FSRUp, events.FSRDown)
            # Renorm/fact scale
            weights_obj_base.add('renorm', events.nom, events.renormUp*(sow/sow_renormUp), events.renormDown*(sow/sow_renormDown))
            weights_obj_base.add('fact', events.nom, events.factUp*(sow/sow_factUp), events.factDown*(sow/sow_factDown))
            if not is2022:
                # Misc other experimental SFs and systs
                weights_obj_base.add('PreFiring', events.L1PreFiringWeight.Nom,  events.L1PreFiringWeight.Up,  events.L1PreFiringWeight.Dn)
                weights_obj_base.add('PU', cor_tc.GetPUSF((events.Pileup.nTrueInt), year), cor_tc.GetPUSF(events.Pileup.nTrueInt, year, 'up'), cor_tc.GetPUSF(events.Pileup.nTrueInt, year, 'down'))
            else:
                weights_obj_base.add("PU", cor_ec.run3_pu_attach(events.Pileup,year,"nominal"), cor_ec.run3_pu_attach(events.Pileup,year,"hi"), cor_ec.run3_pu_attach(events.Pileup,year,"lo"))

            # Lepton SFs and systs
            weights_obj_base.add("lepSF_muon", events.sf_4l_muon, copy.deepcopy(events.sf_4l_hi_muon), copy.deepcopy(events.sf_4l_lo_muon))
            weights_obj_base.add("lepSF_elec", events.sf_4l_elec, copy.deepcopy(events.sf_4l_hi_elec), copy.deepcopy(events.sf_4l_lo_elec))


        # Set up the list of systematics that are handled via event weight variations
        wgt_correction_syst_lst_common = [
            "btagSFbc_correlated", f"btagSFbc_uncorrelated_{year}",
            "lepSF_elec", "lepSF_muon", "PU",
            "renorm", "fact", "ISR", "FSR",
        ]
        if not is2022:
            wgt_correction_syst_lst = wgt_correction_syst_lst_common + ["PreFiring","btagSFlight_correlated",f"btagSFlight_uncorrelated_{year}"]
        else:
            wgt_correction_syst_lst = wgt_correction_syst_lst_common + ["btagSFlight"]
        wgt_correction_syst_lst = append_up_down_to_sys_base(wgt_correction_syst_lst)


        ######### The rest of the processor is inside this loop over systs that affect object kinematics  ###########

        obj_correction_systs = [] # Will have e.g. jes etc

        # If we're doing systematics and this isn't data, we will loop over the obj correction syst lst list
        if self._do_systematics and not isData: obj_corr_syst_var_list = ["nominal"] + obj_correction_systs
        # Otherwise loop juse once, for nominal
        else: obj_corr_syst_var_list = ['nominal']

        # Loop over the list of systematic variations (that impact object kinematics) that we've constructed
        for obj_corr_syst_var in obj_corr_syst_var_list:
            # Make a copy of the base weights object, so that each time through the loop we do not double count systs
            # In this loop over systs that impact kinematics, we will add to the weights objects the SFs that depend on the object kinematics
            weights_obj_base_for_kinematic_syst = copy.deepcopy(weights_obj_base)


            #################### Jets ####################

            # Clean with dr (though another option is to use jetIdx)
            cleanedJets = os_ec.get_cleaned_collection(l_wwz_t,jets)

            # Selecting jets and cleaning them
            jetptname = "pt_nom" if hasattr(cleanedJets, "pt_nom") else "pt"
            cleanedJets["is_good"] = os_tc.is_tight_jet(getattr(cleanedJets, jetptname), cleanedJets.eta, cleanedJets.jetId, pt_cut=20., eta_cut=get_ec_param("wwz_eta_j_cut"), id_cut=get_ec_param("wwz_jet_id_cut"))
            goodJets = cleanedJets[cleanedJets.is_good]

            # Count jets
            njets = ak.num(goodJets)
            ht = ak.sum(goodJets.pt,axis=-1)
            j0 = goodJets[ak.argmax(goodJets.pt,axis=-1,keepdims=True)]


            # Loose DeepJet WP
            btagger = "btag" # For deep flavor WPs
            #btagger = "btagcsv" # For deep CSV WPs
            if year == "2017":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_UL17")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_UL17")
            elif year == "2018":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_UL18")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_UL18")
            elif year=="2016":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_UL16")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_UL16")
            elif year=="2016APV":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_UL16APV")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_UL16APV")
            elif year=="2022":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_2022")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_2022")
            elif year=="2022EE":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_2022EE")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_2022EE")
            else:
                raise ValueError(f"Error: Unknown year \"{year}\".")

            if btagger == "btag":
                isBtagJetsLoose = (goodJets.btagDeepFlavB > btagwpl)
                isBtagJetsMedium = (goodJets.btagDeepFlavB > btagwpm)
            if btagger == "btagcsv":
                isBtagJetsLoose = (goodJets.btagDeepB > btagwpl)
                isBtagJetsMedium = (goodJets.btagDeepB > btagwpm)

            isNotBtagJetsLoose = np.invert(isBtagJetsLoose)
            nbtagsl = ak.num(goodJets[isBtagJetsLoose])

            isNotBtagJetsMedium = np.invert(isBtagJetsMedium)
            nbtagsm = ak.num(goodJets[isBtagJetsMedium])


            ######### Apply SFs #########

            if not isData:

                ### Evaluate btag weights ###
                jets_light = goodJets[goodJets.hadronFlavour==0]
                jets_bc    = goodJets[goodJets.hadronFlavour>0]

                # Workaround to use UL16APV SFs for UL16 for light jets
                year_light = year
                if year == "2016": year_light = "2016APV"

                if not is2022:
                    btag_sf_light = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_incl","central")
                else:
                    btag_sf_light = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_light","central")
                btag_sf_bc = cor_tc.btag_sf_eval(jets_bc,"L",year,"deepJet_comb","central")

                btag_eff_light = cor_ec.btag_eff_eval(jets_light,"L",year)
                btag_eff_bc = cor_ec.btag_eff_eval(jets_bc,"L",year)

                wgt_light = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light, jets_light.btagDeepFlavB>btagwpl)
                wgt_bc    = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc,    jets_bc.btagDeepFlavB>btagwpl)

                wgt_btag_nom = wgt_light*wgt_bc
                weights_obj_base_for_kinematic_syst.add("btagSF", wgt_btag_nom)

                # Put the btagging up and down weight variations into the weights object
                if self._do_systematics:

                    # Run3 2022 btagging systematics stuff
                    # Note light correlated and uncorrelated are missing, so just using total, as suggested by the pog
                    # See this for more info: https://cms-talk.web.cern.ch/t/2022-btag-sf-recommendations/42262
                    if is2022:
                        for btag_sys in ["correlated", "uncorrelated"]:
                            year_tag = f"_{year}"
                            if btag_sys == "correlated": year_tag = ""
                            btag_sf_bc_up      = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"up_{btag_sys}")
                            btag_sf_bc_down    = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"down_{btag_sys}")
                            wgt_bc_up      = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_up,    jets_bc.btagDeepFlavB>btagwpl)
                            wgt_bc_down    = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_down,    jets_bc.btagDeepFlavB>btagwpl)
                            # Note, up and down weights scaled by 1/wgt_btag_nom so that don't double count the central btag correction (i.e. don't apply it also in the case of up and down variations)
                            weights_obj_base_for_kinematic_syst.add(f"btagSFbc_{btag_sys}{year_tag}",    events.nom, wgt_light*wgt_bc_up/wgt_btag_nom, wgt_light*wgt_bc_down/wgt_btag_nom)

                        # Light have no correlated/uncorrelated so just use total:
                        btag_sf_light_up   = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_light","up")
                        btag_sf_light_down = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_light","down")
                        wgt_light_up   = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_up, jets_light.btagDeepFlavB>btagwpl)
                        wgt_light_down = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_down, jets_light.btagDeepFlavB>btagwpl)
                        # Note, up and down weights scaled by 1/wgt_btag_nom so that don't double count the central btag correction (i.e. don't apply it also in the case of up and down variations)
                        weights_obj_base_for_kinematic_syst.add("btagSFlight", events.nom, wgt_light_up*wgt_bc/wgt_btag_nom, wgt_light_down*wgt_bc/wgt_btag_nom)

                    # Run2 btagging systematics stuff
                    else:
                        for btag_sys in ["correlated", "uncorrelated"]:
                            year_tag = f"_{year}"
                            if btag_sys == "correlated": year_tag = ""

                            btag_sf_light_up   = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_incl",f"up_{btag_sys}")
                            btag_sf_light_down = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_incl",f"down_{btag_sys}")
                            btag_sf_bc_up      = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"up_{btag_sys}")
                            btag_sf_bc_down    = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"down_{btag_sys}")

                            wgt_light_up   = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_up, jets_light.btagDeepFlavB>btagwpl)
                            wgt_bc_up      = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_up,    jets_bc.btagDeepFlavB>btagwpl)
                            wgt_light_down = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_down, jets_light.btagDeepFlavB>btagwpl)
                            wgt_bc_down    = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_down,    jets_bc.btagDeepFlavB>btagwpl)

                            # Note, up and down weights scaled by 1/wgt_btag_nom so that don't double count the central btag correction (i.e. don't apply it also in the case of up and down variations)
                            weights_obj_base_for_kinematic_syst.add(f"btagSFlight_{btag_sys}{year_tag}", events.nom, wgt_light_up*wgt_bc/wgt_btag_nom, wgt_light_down*wgt_bc/wgt_btag_nom)
                            weights_obj_base_for_kinematic_syst.add(f"btagSFbc_{btag_sys}{year_tag}",    events.nom, wgt_light*wgt_bc_up/wgt_btag_nom, wgt_light*wgt_bc_down/wgt_btag_nom)


            ######### Masks we need for the selection ##########

            # Pass trigger mask
            pass_trg = es_tc.trg_pass_no_overlap(events,isData,dataset,str(year),dataset_dict=es_ec.dataset_dict,exclude_dict=es_ec.exclude_dict,era=era)
            # Skip trigger matching requirementmes for R3 for now
            if not is2022:
                pass_trg = (pass_trg & es_ec.trg_matching(events,year))

            # b jet masks
            bmask_atleast1med_atleast2loose = ((nbtagsm>=1)&(nbtagsl>=2)) # Used for 2lss and 4l
            bmask_exactly0loose = (nbtagsl==0) # Used for 4l WWZ SR
            bmask_exactly0med = (nbtagsm==0) # Used for 3l CR and 2los Z CR
            bmask_exactly1med = (nbtagsm==1) # Used for 3l SR and 2lss CR
            bmask_exactly2med = (nbtagsm==2) # Used for CRtt
            bmask_atleast2med = (nbtagsm>=2) # Used for 3l SR
            bmask_atmost2med  = (nbtagsm< 3) # Used to make 2lss mutually exclusive from tttt enriched
            bmask_atleast3med = (nbtagsm>=3) # Used for tttt enriched
            bmask_atleast1med = (nbtagsm>=1)
            bmask_atleast1loose = (nbtagsl>=1)
            bmask_atleast2loose = (nbtagsl>=2)


            ######### WWZ event selection stuff #########

            # Get some preliminary things we'll need
            es_ec.attach_wwz_preselection_mask(events,l_wwz_t_padded[:,0:4]) # Attach preselection sf and of flags to the events
            leps_from_z_candidate_ptordered, leps_not_z_candidate_ptordered = es_ec.get_wwz_candidates(l_wwz_t_padded[:,0:4]) # Get ahold of the leptons from the Z and from the W

            w_lep0 = leps_not_z_candidate_ptordered[:,0]
            w_lep1 = leps_not_z_candidate_ptordered[:,1]
            mll_wl0_wl1 = (w_lep0 + w_lep1).mass

            # Make masks for the SF regions
            w_candidates_mll_far_from_z = ak.fill_none(abs(mll_wl0_wl1 - get_ec_param("zmass")) > 10.0,False) # Will enforce this for SF in the PackedSelection
            ptl4 = (l0+l1+l2+l3).pt
            sf_A = ak.fill_none(met.pt >= 120.0,False) # This should never be None, but just keep syntax same as other categories
            sf_B = ak.fill_none((met.pt >= 65.0) & (met.pt < 120.0) & (ptl4 >= 70.0),False)
            sf_C = ak.fill_none((met.pt >= 65.0) & (met.pt < 120.0) & (ptl4 >= 40.0) & (ptl4 < 70.0),False)

            # Make masks for the OF regions
            of_1 = ak.fill_none((mll_wl0_wl1 >= 0.0)  & (mll_wl0_wl1 < 40.0),False)
            of_2 = ak.fill_none((mll_wl0_wl1 >= 40.0) & (mll_wl0_wl1 < 60.0),False)
            of_3 = ak.fill_none((mll_wl0_wl1 >= 60.0) & (mll_wl0_wl1 < 100.0),False)
            of_4 = ak.fill_none((mll_wl0_wl1 >= 100.0),False)

            # Mask for mt2 cut
            mt2_val = es_ec.get_mt2(w_lep0,w_lep1,met)
            mt2_mask = ak.fill_none(mt2_val>25.0,False)




            ######### Get variables we haven't already calculated #########

            l0pt = l0.pt
            j0pt = ak.flatten(j0.pt) # Flatten to go from [[j0pt],[j0pt],...] -> [j0pt,j0pt,...]
            mll_01 = (l0+l1).mass
            mllll = (l0+l1+l2+l3).mass
            scalarptsum_lep = l0.pt + l1.pt + l2.pt + l3.pt
            scalarptsum_lepmet = l0.pt + l1.pt + l2.pt + l3.pt + met.pt
            scalarptsum_lepmetjet = l0.pt + l1.pt + l2.pt + l3.pt + met.pt + ak.sum(goodJets.pt,axis=-1)
            scalarptsum_jet = ak.sum(goodJets.pt,axis=-1)

            # Get lep from Z
            z_lep0 = leps_from_z_candidate_ptordered[:,0]
            z_lep1 = leps_from_z_candidate_ptordered[:,1]

            mll_zl0_zl1 = (z_lep0 + z_lep1).mass

            pt_zl0_zl1 = (z_lep0 + z_lep1).pt
            pt_wl0_wl1 = (w_lep0 + w_lep1).pt

            dr_zl0_zl1 = z_lep0.delta_r(z_lep1)
            dr_wl0_wl1 = w_lep0.delta_r(w_lep1)
            dr_wleps_zleps = (w_lep0 + w_lep1).delta_r(z_lep0 + z_lep1)

            dr_wl0_j_min = ak.min(w_lep0.delta_r(goodJets),axis=-1)
            dr_wl1_j_min = ak.min(w_lep1.delta_r(goodJets),axis=-1)
            dr_wl0_j_min = ak.where(njets>0,dr_wl0_j_min,0)
            dr_wl1_j_min = ak.where(njets>0,dr_wl1_j_min,0)

            absdphi_zl0_zl1 = abs(z_lep0.delta_phi(z_lep1))
            absdphi_wl0_wl1 = abs(w_lep0.delta_phi(w_lep1))
            absdphi_z_ww = abs((z_lep0 + z_lep1).delta_phi(w_lep0 + w_lep1 + met))
            absdphi_4l_met = abs((z_lep0 + z_lep1 + w_lep0 + w_lep1).delta_phi(met))
            absdphi_wleps_met = abs((w_lep0 + w_lep1).delta_phi(met))
            absdphi_wl0_met = abs(w_lep0.delta_phi(met))
            absdphi_wl1_met = abs(w_lep1.delta_phi(met))
            absdphi_zleps_met = abs((z_lep0 + z_lep1).delta_phi(met))

            # Transverse mass
            mt_4l_met    = es_ec.get_mt((l0+l1+l2+l3),met)
            mt_wleps_met = es_ec.get_mt((w_lep0+w_lep1),met)
            mt_wl0_met   = es_ec.get_mt((w_lep0),met)
            mt_wl1_met   = es_ec.get_mt((w_lep1),met)

            # The helicity w0 w1 variable
            cos_helicity_x = es_ec.helicity(w_lep0,w_lep1)

            # lb pairs (i.e. always one lep, one bjet)
            bjets = goodJets[isBtagJetsLoose]
            lb_pairs = ak.cartesian({"l":l_wwz_t,"j":bjets})
            mlb_min = ak.min((lb_pairs["l"] + lb_pairs["j"]).mass,axis=-1)
            mlb_max = ak.max((lb_pairs["l"] + lb_pairs["j"]).mass,axis=-1)

            # Put the variables into a dictionary for easy access later
            dense_variables_dict = {
                "mt2" : mt2_val,
                "met" : met.pt,
                "metphi" : met.phi,
                "ptl4" : ptl4,
                "scalarptsum_lep" : scalarptsum_lep,
                "scalarptsum_lepmet" : scalarptsum_lepmet,
                "scalarptsum_lepmetjet" : scalarptsum_lepmetjet,
                "scalarptsum_jet" : scalarptsum_jet,
                "mll_01" : mll_01,
                "mllll" : mllll,
                "l0pt" : l0pt,
                "j0pt" : j0pt,

                "z_lep0_pt" : z_lep0.pt,
                "z_lep1_pt" : z_lep1.pt,
                "w_lep0_pt" : w_lep0.pt,
                "w_lep1_pt" : w_lep1.pt,
                "z_lep0_eta" : z_lep0.eta,
                "z_lep1_eta" : z_lep1.eta,
                "w_lep0_eta" : w_lep0.eta,
                "w_lep1_eta" : w_lep1.eta,
                "z_lep0_phi" : z_lep0.phi,
                "z_lep1_phi" : z_lep1.phi,
                "w_lep0_phi" : w_lep0.phi,
                "w_lep1_phi" : w_lep1.phi,

                "mll_wl0_wl1" : mll_wl0_wl1,
                "mll_zl0_zl1" : mll_zl0_zl1,

                "pt_zl0_zl1" : pt_zl0_zl1,
                "pt_wl0_wl1" : pt_wl0_wl1,
                "dr_zl0_zl1" : dr_zl0_zl1,
                "dr_wl0_wl1" : dr_wl0_wl1,
                "dr_wleps_zleps" : dr_wleps_zleps,
                "dr_wl0_j_min" : dr_wl0_j_min,
                "dr_wl1_j_min" : dr_wl1_j_min,
                "mt_4l_met" : mt_4l_met,
                "mt_wleps_met" : mt_wleps_met,
                "mt_wl0_met" : mt_wl0_met,
                "mt_wl1_met" : mt_wl1_met,
                "absdphi_zl0_zl1" : absdphi_zl0_zl1,
                "absdphi_wl0_wl1" : absdphi_wl0_wl1,
                "absdphi_z_ww" : absdphi_z_ww,
                "absdphi_4l_met" : absdphi_4l_met,
                "absdphi_zleps_met" : absdphi_zleps_met,
                "absdphi_wleps_met" : absdphi_wleps_met,
                "absdphi_wl0_met" : absdphi_wl0_met,
                "absdphi_wl1_met" : absdphi_wl1_met,

                "cos_helicity_x" : cos_helicity_x,

                "nleps" : nleps,
                "njets" : njets,
                "nbtagsl" : nbtagsl,

                "nleps_counts" : nleps,
                "njets_counts" : njets,
                "nbtagsl_counts" : nbtagsl,

                "absdphi_min_afas" : absdphi_min_afas,
                "absdphi_min_afos" : absdphi_min_afos,
                "absdphi_min_sfos" : absdphi_min_sfos,
                "mll_min_afas" : mll_min_afas,
                "mll_min_afos" : mll_min_afos,
                "mll_min_sfos" : mll_min_sfos,

                "mlb_min" : mlb_min,
                "mlb_max" : mlb_max,

            }


            ######### Evaluate the BDTs (get WWZ, ZH, and WZ scores for SF and OF) #########

            if not is2022:

                # Get the list of variables for the BDTs (and fill None with -9999 to not cause problems), and eval
                bdt_vars_sf_wwz = fill_none_in_list(get_ec_param("sf_bdt_var_lst"),dense_variables_dict,-9999)
                bdt_vars_of_wwz = fill_none_in_list(get_ec_param("of_bdt_var_lst"),dense_variables_dict,-9999)

                # Evaluate BDT v7
                bdt_of_tern = ak.Array(es_ec.eval_of_tern_bdt(bdt_vars_of_wwz))
                bdt_of_wwz = bdt_of_tern[:, 0]
                bdt_of_zh = bdt_of_tern[:, 1]
                bdt_of_bkg = bdt_of_tern[:, 2]
                bdt_sf_tern = ak.Array(es_ec.eval_sf_tern_bdt(bdt_vars_sf_wwz))
                bdt_sf_wwz = bdt_sf_tern[:, 0]
                bdt_sf_zh = bdt_sf_tern[:, 1]
                bdt_sf_bkg = bdt_sf_tern[:, 2]
                bdt_of_wwz_m_zh = bdt_of_wwz - bdt_of_zh
                bdt_sf_wwz_m_zh = bdt_sf_wwz - bdt_sf_zh

                # Philip's version of the v7 binning

                of_thr_zh_1 = 0.08
                of_thr_zh_2 = 0.14
                of_thr_zh_3 = 0.39

                of_thr_wwz_1 = 0.04
                of_thr_wwz_2 = 0.13
                of_thr_wwz_3 = 0.29

                sf_thr_zh_1 = 0.04
                sf_thr_zh_2 = 0.11
                sf_thr_zh_3 = 0.30

                sf_thr_wwz_1 = 0.03
                sf_thr_wwz_2 = 0.05
                sf_thr_wwz_3 = 0.18

                bdt_of_wwz_vs_zh_divider = bdt_of_wwz_m_zh
                bdt_sf_wwz_vs_zh_divider = bdt_sf_wwz_m_zh

                bdt_of_wwz_vs_zh_divider_threshold = 0
                bdt_sf_wwz_vs_zh_divider_threshold = 0

                # Calculating the masks for OF bins
                bdt_of_1 =                                (bdt_of_bkg < of_thr_wwz_1) & (bdt_of_wwz_vs_zh_divider > bdt_of_wwz_vs_zh_divider_threshold)
                bdt_of_2 = (bdt_of_bkg >= of_thr_wwz_1) & (bdt_of_bkg < of_thr_wwz_2) & (bdt_of_wwz_vs_zh_divider > bdt_of_wwz_vs_zh_divider_threshold)
                bdt_of_3 = (bdt_of_bkg >= of_thr_wwz_2) & (bdt_of_bkg < of_thr_wwz_3) & (bdt_of_wwz_vs_zh_divider > bdt_of_wwz_vs_zh_divider_threshold)
                bdt_of_4 = (bdt_of_bkg >= of_thr_wwz_3)                               & (bdt_of_wwz_vs_zh_divider > bdt_of_wwz_vs_zh_divider_threshold)
                bdt_of_5 =                                (bdt_of_bkg < of_thr_zh_1)  & (bdt_of_wwz_vs_zh_divider <= bdt_of_wwz_vs_zh_divider_threshold)
                bdt_of_6 = (bdt_of_bkg >= of_thr_zh_1)  & (bdt_of_bkg < of_thr_zh_2)  & (bdt_of_wwz_vs_zh_divider <= bdt_of_wwz_vs_zh_divider_threshold)
                bdt_of_7 = (bdt_of_bkg >= of_thr_zh_2)  & (bdt_of_bkg < of_thr_zh_3)  & (bdt_of_wwz_vs_zh_divider <= bdt_of_wwz_vs_zh_divider_threshold)
                bdt_of_8 = (bdt_of_bkg >= of_thr_zh_3)                                & (bdt_of_wwz_vs_zh_divider <= bdt_of_wwz_vs_zh_divider_threshold)
                bdt_of_bin = ak.full_like(events.nom,-999)
                bdt_of_bin = ak.where(bdt_of_1, 0, bdt_of_bin)
                bdt_of_bin = ak.where(bdt_of_2, 1, bdt_of_bin)
                bdt_of_bin = ak.where(bdt_of_3, 2, bdt_of_bin)
                bdt_of_bin = ak.where(bdt_of_4, 3, bdt_of_bin)
                bdt_of_bin = ak.where(bdt_of_5, 4, bdt_of_bin)
                bdt_of_bin = ak.where(bdt_of_6, 5, bdt_of_bin)
                bdt_of_bin = ak.where(bdt_of_7, 6, bdt_of_bin)
                bdt_of_bin = ak.where(bdt_of_8, 7, bdt_of_bin)

                # Calculating the bins and computing bin index for each event
                bdt_sf_1 =                                (bdt_sf_bkg < sf_thr_wwz_1) & (bdt_sf_wwz_vs_zh_divider > bdt_sf_wwz_vs_zh_divider_threshold)
                bdt_sf_2 = (bdt_sf_bkg >= sf_thr_wwz_1) & (bdt_sf_bkg < sf_thr_wwz_2) & (bdt_sf_wwz_vs_zh_divider > bdt_sf_wwz_vs_zh_divider_threshold)
                bdt_sf_3 = (bdt_sf_bkg >= sf_thr_wwz_2) & (bdt_sf_bkg < sf_thr_wwz_3) & (bdt_sf_wwz_vs_zh_divider > bdt_sf_wwz_vs_zh_divider_threshold)
                bdt_sf_4 = (bdt_sf_bkg >= sf_thr_wwz_3)                               & (bdt_sf_wwz_vs_zh_divider > bdt_sf_wwz_vs_zh_divider_threshold)
                bdt_sf_5 =                                (bdt_sf_bkg < sf_thr_zh_1)  & (bdt_sf_wwz_vs_zh_divider <= bdt_sf_wwz_vs_zh_divider_threshold)
                bdt_sf_6 = (bdt_sf_bkg >= sf_thr_zh_1)  & (bdt_sf_bkg < sf_thr_zh_2)  & (bdt_sf_wwz_vs_zh_divider <= bdt_sf_wwz_vs_zh_divider_threshold)
                bdt_sf_7 = (bdt_sf_bkg >= sf_thr_zh_2)  & (bdt_sf_bkg < sf_thr_zh_3)  & (bdt_sf_wwz_vs_zh_divider <= bdt_sf_wwz_vs_zh_divider_threshold)
                bdt_sf_8 = (bdt_sf_bkg >= sf_thr_zh_3)                                & (bdt_sf_wwz_vs_zh_divider <= bdt_sf_wwz_vs_zh_divider_threshold)
                bdt_sf_bin = ak.full_like(events.nom,-999)
                bdt_sf_bin = ak.where(bdt_sf_1, 0, bdt_sf_bin)
                bdt_sf_bin = ak.where(bdt_sf_2, 1, bdt_sf_bin)
                bdt_sf_bin = ak.where(bdt_sf_3, 2, bdt_sf_bin)
                bdt_sf_bin = ak.where(bdt_sf_4, 3, bdt_sf_bin)
                bdt_sf_bin = ak.where(bdt_sf_5, 4, bdt_sf_bin)
                bdt_sf_bin = ak.where(bdt_sf_6, 5, bdt_sf_bin)
                bdt_sf_bin = ak.where(bdt_sf_7, 6, bdt_sf_bin)
                bdt_sf_bin = ak.where(bdt_sf_8, 7, bdt_sf_bin)

                # Creating the event mask for BDT regions when split between WWZ vs. ZH
                bdt_of_bin_wwz = (bdt_of_wwz_m_zh > 0)
                bdt_of_bin_zh  = (bdt_of_wwz_m_zh <= 0)
                bdt_sf_bin_wwz = (bdt_sf_wwz_m_zh > 0)
                bdt_sf_bin_zh  = (bdt_sf_wwz_m_zh <= 0)

                # Put the bdt variables into the dict of variables too
                dense_variables_dict["bdt_of_wwz"]      = bdt_of_wwz
                dense_variables_dict["bdt_sf_wwz"]      = bdt_sf_wwz
                dense_variables_dict["bdt_of_zh"]       = bdt_of_zh
                dense_variables_dict["bdt_sf_zh"]       = bdt_sf_zh
                dense_variables_dict["bdt_of_bkg"]      = bdt_of_bkg
                dense_variables_dict["bdt_sf_bkg"]      = bdt_sf_bkg
                dense_variables_dict["bdt_of_wwz_m_zh"] = bdt_of_wwz_m_zh
                dense_variables_dict["bdt_sf_wwz_m_zh"] = bdt_sf_wwz_m_zh
                dense_variables_dict["bdt_of_bin"]      = bdt_of_bin
                dense_variables_dict["bdt_sf_bin"]      = bdt_sf_bin




            ######### Store boolean masks with PackedSelection ##########

            selections = PackedSelection(dtype='uint64')

            # Lumi mask (for data)
            selections.add("is_good_lumi",lumi_mask)

            zeroj = (njets==0)

            # Fill the packed slection object

            # CRs
            ww_ee = ((abs(w_lep0.pdgId) == 11) & (abs(w_lep1.pdgId) == 11))
            ww_mm = ((abs(w_lep0.pdgId) == 13) & (abs(w_lep1.pdgId) == 13))
            ww_em = ((abs(w_lep0.pdgId) == 11) & (abs(w_lep1.pdgId) == 13))
            ww_me = ((abs(w_lep0.pdgId) == 13) & (abs(w_lep1.pdgId) == 11))
            selections.add("cr_4l_btag_of",            (pass_trg & events.is4lWWZ & bmask_atleast1loose & events.wwz_presel_of))
            selections.add("cr_4l_btag_sf_offZ_met80", (pass_trg & events.is4lWWZ & bmask_atleast1loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & (met.pt > 80.0)))
            selections.add("cr_4l_sf", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & (~w_candidates_mll_far_from_z)))

            selections.add("sr_4l_sf_A", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & sf_A))
            selections.add("sr_4l_sf_B", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & sf_B))
            selections.add("sr_4l_sf_C", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & sf_C))
            selections.add("sr_4l_of_1", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of & of_1 & mt2_mask))
            selections.add("sr_4l_of_2", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of & of_2 & mt2_mask))
            selections.add("sr_4l_of_3", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of & of_3 & mt2_mask))
            selections.add("sr_4l_of_4", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of & of_4))

            selections.add("all_events", (events.is4lWWZ | (~events.is4lWWZ))) # All events.. this logic is a bit roundabout to just get an array of True
            selections.add("4l_presel", (events.is4lWWZ)) # This matches the VVV looper selection (object selection and event selection)

            selections.add("sr_4l_sf_incl", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & (met.pt >= 65.0))) # Inclusive over SF sr (only applying cuts that are applied to all SF SRs), just use for visualization
            selections.add("sr_4l_of_incl", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of)) # Inclusive over OF sr (only applying cuts that are applied to all OF SRs), just use for visualization

            # For BDT SRs

            if not is2022:
                sr_4l_bdt_sf_presel = (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z)
                sr_4l_bdt_sf_trn    = (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & mt2_mask)
                sr_4l_bdt_of_presel = (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of)
                sr_4l_bdt_of_trn    = sr_4l_bdt_of_presel # For OF, presel and trn regions are the same
                selections.add("sr_4l_bdt_sf_presel", sr_4l_bdt_sf_presel)
                selections.add("sr_4l_bdt_sf_trn"   , sr_4l_bdt_sf_trn)
                selections.add("sr_4l_bdt_of_presel", sr_4l_bdt_of_presel)
                selections.add("sr_4l_bdt_of_trn"   , sr_4l_bdt_of_trn)

                selections.add("sr_4l_bdt_sf_1", (sr_4l_bdt_sf_trn & bdt_sf_1))
                selections.add("sr_4l_bdt_sf_2", (sr_4l_bdt_sf_trn & bdt_sf_2))
                selections.add("sr_4l_bdt_sf_3", (sr_4l_bdt_sf_trn & bdt_sf_3))
                selections.add("sr_4l_bdt_sf_4", (sr_4l_bdt_sf_trn & bdt_sf_4))
                selections.add("sr_4l_bdt_sf_5", (sr_4l_bdt_sf_trn & bdt_sf_5))
                selections.add("sr_4l_bdt_sf_6", (sr_4l_bdt_sf_trn & bdt_sf_6))
                selections.add("sr_4l_bdt_sf_7", (sr_4l_bdt_sf_trn & bdt_sf_7))
                selections.add("sr_4l_bdt_sf_8", (sr_4l_bdt_sf_trn & bdt_sf_8))

                selections.add("sr_4l_bdt_of_1", (sr_4l_bdt_of_trn & bdt_of_1))
                selections.add("sr_4l_bdt_of_2", (sr_4l_bdt_of_trn & bdt_of_2))
                selections.add("sr_4l_bdt_of_3", (sr_4l_bdt_of_trn & bdt_of_3))
                selections.add("sr_4l_bdt_of_4", (sr_4l_bdt_of_trn & bdt_of_4))
                selections.add("sr_4l_bdt_of_5", (sr_4l_bdt_of_trn & bdt_of_5))
                selections.add("sr_4l_bdt_of_6", (sr_4l_bdt_of_trn & bdt_of_6))
                selections.add("sr_4l_bdt_of_7", (sr_4l_bdt_of_trn & bdt_of_7))
                selections.add("sr_4l_bdt_of_8", (sr_4l_bdt_of_trn & bdt_of_8))

                selections.add("sr_4l_bdt_of_wwz", (sr_4l_bdt_of_trn & bdt_of_bin_wwz))
                selections.add("sr_4l_bdt_of_zh" , (sr_4l_bdt_of_trn & bdt_of_bin_zh))
                selections.add("sr_4l_bdt_sf_wwz", (sr_4l_bdt_sf_trn & bdt_sf_bin_wwz))
                selections.add("sr_4l_bdt_sf_zh" , (sr_4l_bdt_sf_trn & bdt_sf_bin_zh))

                bdt_sr_names = [
                    "sr_4l_bdt_sf_1",
                    "sr_4l_bdt_sf_2",
                    "sr_4l_bdt_sf_3",
                    "sr_4l_bdt_sf_4",
                    "sr_4l_bdt_sf_5",
                    "sr_4l_bdt_sf_6",
                    "sr_4l_bdt_sf_7",
                    "sr_4l_bdt_sf_8",

                    "sr_4l_bdt_of_1",
                    "sr_4l_bdt_of_2",
                    "sr_4l_bdt_of_3",
                    "sr_4l_bdt_of_4",
                    "sr_4l_bdt_of_5",
                    "sr_4l_bdt_of_6",
                    "sr_4l_bdt_of_7",
                    "sr_4l_bdt_of_8",

                    "sr_4l_bdt_sf_wwz",
                    "sr_4l_bdt_sf_zh",

                    "sr_4l_bdt_of_wwz",
                    "sr_4l_bdt_of_zh",
                ]
                bdt_misc_names = [
                    "sr_4l_bdt_sf_presel",
                    "sr_4l_bdt_sf_trn",
                    "sr_4l_bdt_of_presel",
                    "sr_4l_bdt_of_trn",
                ]

            cat_dict = {
                "lep_chan_lst" : [
                    "sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C","sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4",
                    "all_events","4l_presel", "sr_4l_sf_incl", "sr_4l_of_incl",
                    "cr_4l_btag_of", "cr_4l_btag_sf_offZ_met80", "cr_4l_sf",
                ]
            }

            if not is2022:
                cat_dict["lep_chan_lst"] = cat_dict["lep_chan_lst"] + bdt_sr_names + bdt_misc_names

            ######### Fill histos #########

            # List the hists that are only defined for some categories
            analysis_cats = ["sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C","sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4"]
            if not is2022:
                analysis_cats = analysis_cats + bdt_sr_names
            exclude_var_dict = {
                "mt2" : ["all_events"],
                "ptl4" : ["all_events"],
                "j0pt" : ["all_events", "4l_presel", "sr_4l_sf_incl", "sr_4l_of_incl", "cr_4l_sf"] + analysis_cats,
                "l0pt" : ["all_events"],
                "mll_01" : ["all_events"],
                "mllll" : ["all_events"],
                "scalarptsum_lep" : ["all_events"],
                "scalarptsum_lepmet" : ["all_events"],
                "scalarptsum_lepmetjet" : ["all_events"],
                "scalarptsum_jet" : ["all_events"],
                "w_lep0_pt"  : ["all_events"],
                "w_lep1_pt"  : ["all_events"],
                "z_lep0_pt"  : ["all_events"],
                "z_lep1_pt"  : ["all_events"],
                "w_lep0_eta" : ["all_events"],
                "w_lep1_eta" : ["all_events"],
                "z_lep0_eta" : ["all_events"],
                "z_lep1_eta" : ["all_events"],
                "w_lep0_phi" : ["all_events"],
                "w_lep1_phi" : ["all_events"],
                "z_lep0_phi" : ["all_events"],
                "z_lep1_phi" : ["all_events"],
                "mll_wl0_wl1" : ["all_events"],
                "mll_zl0_zl1" : ["all_events"],

                "pt_zl0_zl1" : ["all_events"],
                "pt_wl0_wl1" : ["all_events"],
                "dr_zl0_zl1" : ["all_events"],
                "dr_wl0_wl1" : ["all_events"],
                "dr_wleps_zleps" : ["all_events"],
                "dr_wl0_j_min" : ["all_events"],
                "dr_wl1_j_min" : ["all_events"],
                "mt_4l_met" : ["all_events"],
                "mt_wleps_met" : ["all_events"],
                "mt_wl0_met" : ["all_events"],
                "mt_wl1_met" : ["all_events"],
                "absdphi_zl0_zl1" : ["all_events"],
                "absdphi_wl0_wl1" : ["all_events"],
                "absdphi_z_ww" : ["all_events"],
                "absdphi_4l_met" : ["all_events"],
                "absdphi_zleps_met" : ["all_events"],
                "absdphi_wleps_met" : ["all_events"],
                "absdphi_wl0_met" : ["all_events"],
                "absdphi_wl1_met" : ["all_events"],

                "cos_helicity_x" : ["all_events"],

                "absdphi_min_afas" : ["all_events"],
                "absdphi_min_afos" : ["all_events"],
                "absdphi_min_sfos" : ["all_events"],
                "mll_min_afas" : ["all_events"],
                "mll_min_afos" : ["all_events"],
                "mll_min_sfos" : ["all_events"],

                "mlb_min" : ["all_events","4l_presel", "sr_4l_sf_incl", "sr_4l_of_incl", "cr_4l_sf"] + analysis_cats,
                "mlb_max" : ["all_events","4l_presel", "sr_4l_sf_incl", "sr_4l_of_incl", "cr_4l_sf"] + analysis_cats,
            }

            if not is2022:
                exclude_var_dict["bdt_of_wwz"] = ["all_events"]
                exclude_var_dict["bdt_sf_wwz"] = ["all_events"]
                exclude_var_dict["bdt_of_zh"] = ["all_events"]
                exclude_var_dict["bdt_sf_zh"] = ["all_events"]
                exclude_var_dict["bdt_of_bkg"] = ["all_events"]
                exclude_var_dict["bdt_sf_bkg"] = ["all_events"]
                exclude_var_dict["bdt_of_wwz_m_zh"] = ["all_events"]
                exclude_var_dict["bdt_sf_wwz_m_zh"] = ["all_events"]
                exclude_var_dict["j0pt"] = exclude_var_dict["j0pt"] + bdt_misc_names
                exclude_var_dict["mlb_max"] = exclude_var_dict["mlb_max"] + bdt_misc_names
                exclude_var_dict["mlb_min"] = exclude_var_dict["mlb_min"] + bdt_misc_names

            # Set up the list of weight fluctuations to loop over
            # For now the syst do not depend on the category, so we can figure this out outside of the filling loop
            wgt_var_lst = ["nominal"]
            if self._do_systematics:
                if not isData:
                    if (obj_corr_syst_var != "nominal"):
                        # In this case, we are dealing with systs that change the kinematics of the objs (e.g. JES)
                        # So we don't want to loop over up/down weight variations here
                        wgt_var_lst = [obj_corr_syst_var]
                    else:
                        # Otherwise we want to loop over the up/down weight variations
                        wgt_var_lst = wgt_var_lst + wgt_correction_syst_lst



            # Loop over the hists we want to fill
            for dense_axis_name, dense_axis_vals in dense_variables_dict.items():
                if dense_axis_name not in self._hist_lst:
                    #print(f"Skipping \"{dense_axis_name}\", it is not in the list of hists to include.")
                    continue
                #print("\ndense_axis_name,vals",dense_axis_name)
                #print("\ndense_axis_name,vals",vals)

                # Loop over weight fluctuations
                for wgt_fluct in wgt_var_lst:

                    # Get the appropriate weight fluctuation
                    if (wgt_fluct == "nominal") or (wgt_fluct in obj_corr_syst_var_list):
                        # In the case of "nominal", no weight systematic variation is used
                        weight = weights_obj_base_for_kinematic_syst.weight(None)
                    else:
                        # Otherwise get the weight from the Weights object
                        weight = weights_obj_base_for_kinematic_syst.weight(wgt_fluct)


                    # Loop over categories
                    for sr_cat in cat_dict["lep_chan_lst"]:

                        # Skip filling if this variable is not relevant for this selection
                        if (dense_axis_name in exclude_var_dict) and (sr_cat in exclude_var_dict[dense_axis_name]): continue

                        # If this is a counts hist, forget the weights and just fill with unit weights
                        if dense_axis_name.endswith("_counts"): weight = events.nom
                        #else: weights = weights_obj_base_for_kinematic_syst.partial_weight(include=["norm"]) # For testing
                        #else: weights = weights_obj_base_for_kinematic_syst.weight(None) # For testing

                        # Make the cuts mask
                        cuts_lst = [sr_cat]
                        if isData: cuts_lst.append("is_good_lumi") # Apply golden json requirements if this is data
                        all_cuts_mask = selections.all(*cuts_lst)

                        #run = events.run[all_cuts_mask]
                        #luminosityBlock = events.luminosityBlock[all_cuts_mask]
                        #event = events.event[all_cuts_mask]
                        #w = weights[all_cuts_mask]
                        #if dense_axis_name == "njets":
                        #    print("\nSTARTPRINT")
                        #    for i,j in enumerate(w):
                        #        out_str = f"PRINTTAG {i} {dense_axis_name} {year} {sr_cat} {event[i]} {run[i]} {luminosityBlock[i]} {w[i]}"
                        #        print(out_str,file=sys.stderr,flush=True)
                        #    print("ENDPRINT\n")
                        #print("\ndense_axis_name",dense_axis_name)
                        #print("sr_cat",sr_cat)
                        #print("dense_axis_vals[all_cuts_mask]",dense_axis_vals[all_cuts_mask])
                        #print("end")

                        # Fill the histos
                        axes_fill_info_dict = {
                            dense_axis_name : dense_axis_vals[all_cuts_mask],
                            "weight"        : weight[all_cuts_mask],
                            "process"       : histAxisName,
                            "category"      : sr_cat,
                            "systematic"    : wgt_fluct,
                        }
                        self.accumulator[dense_axis_name].fill(**axes_fill_info_dict)

            # Fill the list accumulator
            if self._siphon_bdt_data:
                for chan,mask in {"of": sr_4l_bdt_of_trn, "sf": sr_4l_bdt_sf_trn}.items():
                    self.accumulator[f"bdt_{chan}_wwz_list"]  += bdt_of_wwz[mask].to_list()
                    self.accumulator[f"bdt_{chan}_zh_list"]   += bdt_of_zh[mask].to_list()
                    self.accumulator[f"bdt_{chan}_bkg_list"]  += bdt_of_bkg[mask].to_list()
                    self.accumulator[f"bdt_{chan}_evt_list"]  += events.event[mask].to_list()
                    self.accumulator[f"bdt_{chan}_proc_list"] += [histAxisName] * len(bdt_of_bkg[mask])
                    self.accumulator[f"bdt_{chan}_wgt_list"]  += weights_obj_base_for_kinematic_syst.weight(None)[mask]
                    for ivar, var in enumerate(get_ec_param(f"{chan}_bdt_var_lst")):
                        if chan == "of": self.accumulator[f"{chan}_bdt_{var}"] += bdt_vars_of_wwz[ivar][mask]
                        if chan == "sf": self.accumulator[f"{chan}_bdt_{var}"] += bdt_vars_sf_wwz[ivar][mask]

        return self.accumulator

    def postprocess(self, accumulator):
        return accumulator
