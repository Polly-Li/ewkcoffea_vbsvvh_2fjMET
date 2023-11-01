#!/usr/bin/env python
#import sys
import math
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

from ewkcoffea.modules.paths import ewkcoffea_path as ewkcoffea_path
import ewkcoffea.modules.selection_wwz as es_ec
import ewkcoffea.modules.objects_wwz as os_ec
import ewkcoffea.modules.corrections as cor_ec

from topcoffea.modules.get_param_from_jsons import GetParam
get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))


class AnalysisProcessor(processor.ProcessorABC):

    def __init__(self, samples, wc_names_lst=[], hist_lst=None, ecut_threshold=None, do_errors=False, do_systematics=False, split_by_lepton_flavor=False, skip_signal_regions=False, skip_control_regions=False, muonSyst='nominal', dtype=np.float32):

        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype

        # Create the dense axes for the histograms
        self._dense_axes_dict = {
            "mt2"   : axis.Regular(180, 0, 100, name="mt2",  label="mt2"),
            "met"   : axis.Regular(180, 0, 500, name="met",  label="met"),
            "metphi": axis.Regular(180, -4, 4, name="metphi", label="met phi"),
            "ptl4"  : axis.Regular(180, 0, 500, name="ptl4", label="ptl4"),
            "scalarptsum_lep" : axis.Regular(180, 0, 500, name="scalarptsum_lep", label="S_T"),
            "scalarptsum_lepmet" : axis.Regular(180, 0, 600, name="scalarptsum_lepmet", label="S_T + metpt"),
            "scalarptsum_lepmetjet" : axis.Regular(180, 0, 600, name="scalarptsum_lepmetjet", label="S_T + metpt + H_T"),
            "mll_01": axis.Regular(180, 0, 200, name="mll_01",  label="mll_l0_l1"),
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
            "w_lep0_phi" : axis.Regular(180, -4, 4, name="w_lep0_phi", label="Leading W lep phi"),
            "w_lep1_phi" : axis.Regular(180, -4, 4, name="w_lep1_phi", label="Subleading W lep phi"),
            "z_lep0_phi" : axis.Regular(180, -4, 4, name="z_lep0_phi", label="Leading Z lep phi"),
            "z_lep1_phi" : axis.Regular(180, -4, 4, name="z_lep1_phi", label="Subleading Z lep phi"),
            "mll_wl0_wl1" : axis.Regular(180, 0, 200, name="mll_wl0_wl1", label="mll(W lep0, W lep1)"),
            "mll_zl0_zl1" : axis.Regular(180, 0, 200, name="mll_zl0_zl1", label="mll(Z lep0, Z lep1)"),

            "pt_zl0_zl1" : axis.Regular(180, 0, 300, name="pt_zl0_zl1", label="pt(Zl0 + Zl1)"),
            "pt_wl0_wl1" : axis.Regular(180, 0, 300, name="pt_wl0_wl1", label="pt(Wl0 + Wl1)"),
            "dr_zl0_zl1" : axis.Regular(180, 0, 5, name="dr_zl0_zl1", label="dr(Zl0,Zl1)"),
            "dr_wl0_wl1" : axis.Regular(180, 0, 5, name="dr_wl0_wl1", label="dr(Wl0,Wl1)"),
            "dr_wleps_zleps" : axis.Regular(180, 0, 5, name="dr_wleps_zleps", label="dr((Wl0+Wl1),(Zl0,Zl1))"),

            "absdphi_zl0_zl1" : axis.Regular(180, 0, 4, name="absdphi_zl0_zl1", label="abs dphi(Zl0,Zl1)"),
            "absdphi_wl0_wl1" : axis.Regular(180, 0, 4, name="absdphi_wl0_wl1", label="abs dphi(Wl0,Wl1)"),
            "absdphi_z_ww"    : axis.Regular(180, 0, 4, name="absdphi_z_ww", label="abs dphi((Zl0+Zl1),(Wl0+Wl1+met))"),
            "dphi_4l_met"    : axis.Regular(180, -4, 4, name="dphi_4l_met", label="dphi((Zl0+Zl1+Wl0+Wl1),met)"),
            "dphi_zleps_met" : axis.Regular(180, -4, 4, name="dphi_zleps_met", label="dphi((Zl0+Zl1),met)"),
            "dphi_wleps_met" : axis.Regular(180, -4, 4, name="dphi_wleps_met", label="dphi((Wl0+Wl1),met)"),

            "absdphi_min_afas" : axis.Regular(180, 0, 4, name="absdphi_min_afas",  label="min(abs(delta phi of all pairs))"),
            "absdphi_min_afos" : axis.Regular(180, 0, 4, name="absdphi_min_afos",  label="min(abs(delta phi of OS pairs))"),
            "absdphi_min_sfos" : axis.Regular(180, 0, 4, name="absdphi_min_sfos",  label="min(abs(delta phi of SFOS pairs))"),
            "mll_min_afas" : axis.Regular(180, 0, 150, name="mll_min_afas",  label="min mll of all pairs"),
            "mll_min_afos" : axis.Regular(180, 0, 150, name="mll_min_afos",  label="min mll of OF pairs"),
            "mll_min_sfos" : axis.Regular(180, 0, 150, name="mll_min_sfos",  label="min mll of SFOF pairs"),

            "mlb_min" : axis.Regular(180, 0, 300, name="mlb_min",  label="min mass(b+l)"),
            "mlb_max" : axis.Regular(180, 0, 500, name="mlb_max",  label="max mass(b+l)"),

            "njets"   : axis.Regular(8, 0, 8, name="njets",   label="Jet multiplicity"),
            "nleps"   : axis.Regular(5, 0, 5, name="nleps",   label="Lep multiplicity"),
            "nbtagsl" : axis.Regular(4, 0, 4, name="nbtagsl", label="Loose btag multiplicity"),

            "njets_counts"   : axis.Regular(30, 0, 30, name="njets_counts",   label="Jet multiplicity counts"),
            "nleps_counts"   : axis.Regular(30, 0, 30, name="nleps_counts",   label="Lep multiplicity counts"),
            "nbtagsl_counts" : axis.Regular(30, 0, 30, name="nbtagsl_counts", label="Loose btag multiplicity counts"),

            "bdt_of_wwz_raw": axis.Regular(180, -3.5, 3.5, name="bdt_of_wwz_raw", label="Raw score bdt_of_wwz"),
            "bdt_sf_wwz_raw": axis.Regular(180, -3.5, 3.5, name="bdt_sf_wwz_raw", label="Raw score bdt_sf_wwz"),
            "bdt_of_zh_raw" : axis.Regular(180, -3.5, 3.5, name="bdt_of_zh_raw", label="Raw score bdt_of_zh"),
            "bdt_sf_zh_raw" : axis.Regular(180, -3.5, 3.5, name="bdt_sf_zh_raw", label="Raw score bdt_sf_zh"),

            "bdt_of_wwz": axis.Regular(180, 0, 1, name="bdt_of_wwz", label="Score bdt_of_wwz"),
            "bdt_sf_wwz": axis.Regular(180, 0, 1, name="bdt_sf_wwz", label="Score bdt_sf_wwz"),
            "bdt_of_zh" : axis.Regular(180, 0, 1, name="bdt_of_zh", label="Score bdt_of_zh"),
            "bdt_sf_zh" : axis.Regular(180, 0, 1, name="bdt_sf_zh", label="Score bdt_sf_zh"),

        }

        # Set the list of hists to fill
        if hist_lst is None:
            # If the hist list is none, assume we want to fill all hists
            self._hist_lst = list(self._dense_axes_dict.keys())
        else:
            # Otherwise, just fill the specified subset of hists
            for hist_to_include in hist_lst:
                if hist_to_include not in self._dense_axes_dict.keys():
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


    @property
    def columns(self):
        return self._columns

    # Main function: run on a given dataset
    def process(self, events):

        # Dataset parameters
        dataset = events.metadata["dataset"]

        isData             = self._samples[dataset]["isData"]
        histAxisName       = self._samples[dataset]["histAxisName"]
        year               = self._samples[dataset]["year"]
        xsec               = self._samples[dataset]["xsec"]
        sow                = self._samples[dataset]["nSumOfWeights"]

        # Get up down weights from input dict
        if (self._do_systematics and not isData):
            # Otherwise we have an NLO xsec, so for these systs we will have e.g. xsec_NLO*(N_pass_up/N_gen_up)
            # Thus these systs should only affect acceptance and effeciency and shape
            # The uncty on xsec comes from NLO and is applied as a rate uncty in the text datacard
            sow_ISRUp          = self._samples[dataset]["nSumOfWeights_ISRUp"          ]
            sow_ISRDown        = self._samples[dataset]["nSumOfWeights_ISRDown"        ]
            sow_FSRUp          = self._samples[dataset]["nSumOfWeights_FSRUp"          ]
            sow_FSRDown        = self._samples[dataset]["nSumOfWeights_FSRDown"        ]
            sow_renormUp       = self._samples[dataset]["nSumOfWeights_renormUp"       ]
            sow_renormDown     = self._samples[dataset]["nSumOfWeights_renormDown"     ]
            sow_factUp         = self._samples[dataset]["nSumOfWeights_factUp"         ]
            sow_factDown       = self._samples[dataset]["nSumOfWeights_factDown"       ]
            sow_renormfactUp   = self._samples[dataset]["nSumOfWeights_renormfactUp"   ]
            sow_renormfactDown = self._samples[dataset]["nSumOfWeights_renormfactDown" ]
        else:
            sow_ISRUp          = -1
            sow_ISRDown        = -1
            sow_FSRUp          = -1
            sow_FSRDown        = -1
            sow_renormUp       = -1
            sow_renormDown     = -1
            sow_factUp         = -1
            sow_factDown       = -1
            sow_renormfactUp   = -1
            sow_renormfactDown = -1

        datasets = ["SingleMuon", "SingleElectron", "EGamma", "MuonEG", "DoubleMuon", "DoubleElectron", "DoubleEG"]
        for d in datasets:
            if d in dataset: dataset = dataset.split('_')[0]

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
        else:
            raise ValueError(f"Error: Unknown year \"{year}\".")
        lumi_mask = LumiMask(golden_json_path)(events.run,events.luminosityBlock)


        ################### Lepton selection ####################

        # Do the object selection for the WWZ eleectrons
        ele_presl_mask = os_ec.is_presel_wwz_ele(ele,tight=True)
        ele["topmva"] = os_ec.get_topmva_score_ele(events, year)
        ele["is_tight_lep_for_wwz"] = ((ele.topmva > get_tc_param("topmva_wp_t_e")) & ele_presl_mask)

        # Do the object selection for the WWZ muons
        mu_presl_mask = os_ec.is_presel_wwz_mu(mu)
        mu["topmva"] = os_ec.get_topmva_score_mu(events, year)
        mu["is_tight_lep_for_wwz"] = ((mu.topmva > get_tc_param("topmva_wp_t_m")) & mu_presl_mask)

        # Get tight leptons for WWZ selection
        ele_wwz_t = ele[ele.is_tight_lep_for_wwz]
        mu_wwz_t = mu[mu.is_tight_lep_for_wwz]

        # Attach the lepton SFs to the electron and muons collections
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


        ######### Systematics ###########


        # These weights can go outside of the outside sys loop since they do not depend on pt of mu or jets
        # We only calculate these values if not isData
        # Note: add() will generally modify up/down weights, so if these are needed for any reason after this point, we should instead pass copies to add()
        # Note: Here we will to the weights object the SFs that do not depend on any of the forthcoming loops
        weights_obj_base = coffea.analysis_tools.Weights(len(events),storeIndividual=True)
        if not isData:
            genw = events["genWeight"]

            # Normalize by (xsec/sow)*genw where genw is 1 for EFT samples
            # Note that for theory systs, will need to multiply by sow/sow_wgtUP to get (xsec/sow_wgtUp)*genw and same for Down
            lumi = 1000.0*get_tc_param(f"lumi_{year}")
            weights_obj_base.add("norm",(xsec/sow)*genw*lumi)

        # We do not have systematics yet
        syst_var_list = ['nominal']

        # Loop over the list of systematic variations we've constructed
        for syst_var in syst_var_list:
            # Make a copy of the base weights object, so that each time through the loop we do not double count systs
            # In this loop over systs that impact kinematics, we will add to the weights objects the SFs that depend on the object kinematics
            weights_obj_base_for_kinematic_syst = copy.deepcopy(weights_obj_base)


            #################### Jets ####################

            # Jet cleaning, before any jet selection
            ##vetos_tocleanjets = ak.with_name( ak.concatenate([tau, l_fo], axis=1), "PtEtaPhiMCandidate")
            #vetos_tocleanjets = ak.with_name( l_wwz_t, "PtEtaPhiMCandidate")
            #tmp = ak.cartesian([ak.local_index(jets.pt), vetos_tocleanjets.jetIdx], nested=True)
            #cleanedJets = jets[~ak.any(tmp.slot0 == tmp.slot1, axis=-1)] # this line should go before *any selection*, otherwise lep.jetIdx is not aligned with the jet index

            # Clean with dr for now
            cleanedJets = os_ec.get_cleaned_collection(l_wwz_t,jets)

            # Selecting jets and cleaning them
            # NOTE: The jet id cut is commented for now in objects.py for the sync
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

            #################### Add variables into event object so that they persist ####################

            # Put njets and l_fo_conept_sorted into events
            events["njets"] = njets
            events["l_wwz_t"] = l_wwz_t

            es_ec.add4lmask_wwz(events, year, isData, histAxisName)

            ######### Apply SFs #########

            if not isData:
                weights_obj_base_for_kinematic_syst.add("lepSF_muon", events.sf_4l_muon, copy.deepcopy(events.sf_4l_hi_muon), copy.deepcopy(events.sf_4l_lo_muon))
                weights_obj_base_for_kinematic_syst.add("lepSF_elec", events.sf_4l_elec, copy.deepcopy(events.sf_4l_hi_elec), copy.deepcopy(events.sf_4l_lo_elec))


            ######### Masks we need for the selection ##########

            # Pass trigger mask
            pass_trg = es_tc.trg_pass_no_overlap(events,isData,dataset,str(year),dataset_dict=es_ec.dataset_dict,exclude_dict=es_ec.exclude_dict)
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
            sf_A = (met.pt > 120.0)
            sf_B = ((met.pt > 65.0) & (met.pt < 120.0) & (ptl4 > 70.0))
            sf_C = ((met.pt > 65.0) & (met.pt < 120.0) & (ptl4 > 40.0) & (ptl4 < 70.0))

            # Make masks for the OF regions
            of_1 = ak.fill_none((mll_wl0_wl1 > 0.0) & (mll_wl0_wl1 < 40.0),False)
            of_2 = ak.fill_none((mll_wl0_wl1 > 40.0) & (mll_wl0_wl1 < 60.0),False)
            of_3 = ak.fill_none((mll_wl0_wl1 > 60.0) & (mll_wl0_wl1 < 100.0),False)
            of_4 = ak.fill_none((mll_wl0_wl1 > 100.0),False)

            # Mask for mt2 cut
            mt2_val = es_ec.get_mt2(w_lep0,w_lep1,met)
            mt2_mask = ak.fill_none(mt2_val>25.0,False)



            ######### Store boolean masks with PackedSelection ##########


            selections = PackedSelection(dtype='uint64')

            # Lumi mask (for data)
            selections.add("is_good_lumi",lumi_mask)

            zeroj = (njets==0)

            # For WWZ selection
            selections.add("sr_4l_sf_A", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & sf_A))
            selections.add("sr_4l_sf_B", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & sf_B))
            selections.add("sr_4l_sf_C", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & sf_C))
            selections.add("sr_4l_of_1", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of & of_1 & mt2_mask))
            selections.add("sr_4l_of_2", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of & of_2 & mt2_mask))
            selections.add("sr_4l_of_3", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of & of_3 & mt2_mask))
            selections.add("sr_4l_of_4", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of & of_4))

            selections.add("all_events", (events.is4lWWZ | (~events.is4lWWZ))) # All events.. this logic is a bit roundabout to just get an array of True
            selections.add("4l_presel", (events.is4lWWZ)) # This matches the VVV looper selection (object selection and event selection)

            selections.add("sr_4l_sf_presel", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & (met.pt > 65.0)))
            selections.add("sr_4l_of_presel", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of))

            # CRs
            selections.add("cr_4l_of", (pass_trg & events.is4lWWZ & bmask_atleast1med & events.wwz_presel_of))
            selections.add("cr_4l_sf", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & (~w_candidates_mll_far_from_z)))

            cat_dict = {
                "lep_chan_lst" : [
                    "sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C","sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4",
                    "sr_4l_sf_presel", "sr_4l_of_presel",
                    "all_events","4l_presel",
                    "cr_4l_of","cr_4l_sf",
                ],
            }

            ######### Get variables #########

            l0pt = l0.pt
            j0pt = ak.flatten(j0.pt) # Flatten to go from [[j0pt],[j0pt],...] -> [j0pt,j0pt,...]
            mll_01 = (l0+l1).mass
            scalarptsum_lep = l0.pt + l1.pt + l2.pt + l3.pt
            scalarptsum_lepmet = l0.pt + l1.pt + l2.pt + l3.pt + met.pt
            scalarptsum_lepmetjet = l0.pt + l1.pt + l2.pt + l3.pt + met.pt + ak.sum(goodJets.pt,axis=-1)

            # Get lep from Z
            z_lep0 = leps_from_z_candidate_ptordered[:,0]
            z_lep1 = leps_from_z_candidate_ptordered[:,1]

            mll_zl0_zl1 = (z_lep0 + z_lep1).mass

            pt_zl0_zl1 = (z_lep0 + z_lep1).pt
            pt_wl0_wl1 = (w_lep0 + w_lep1).pt

            dr_zl0_zl1 = z_lep0.delta_r(z_lep1)
            dr_wl0_wl1 = w_lep0.delta_r(w_lep1)
            dr_wleps_zleps = (w_lep0 + w_lep1).delta_r(z_lep0 + z_lep1)

            absdphi_zl0_zl1 = abs(z_lep0.delta_phi(z_lep1))
            absdphi_wl0_wl1 = abs(w_lep0.delta_phi(w_lep1))
            absdphi_z_ww = abs((z_lep0 + z_lep1).delta_phi(w_lep0 + w_lep1 + met))
            dphi_4l_met = (z_lep0 + z_lep1 + w_lep0 + w_lep1).delta_phi(met)
            dphi_wleps_met = (w_lep0 + w_lep1).delta_phi(met)
            dphi_zleps_met = (z_lep0 + z_lep1).delta_phi(met)

            # lb pairs (i.e. always one lep, one bjet)
            bjets = goodJets[isBtagJetsLoose]
            lb_pairs = ak.cartesian({"l":l_wwz_t,"j":bjets})
            mlb_min = ak.min((lb_pairs["l"] + lb_pairs["j"]).mass,axis=-1)
            mlb_max = ak.max((lb_pairs["l"] + lb_pairs["j"]).mass,axis=-1)

            # BDT values
            #print("mll_wl0_wl1",mll_wl0_wl1 )
            #print("dphi_4l_met", dphi_4l_met)
            #print("dphi_zleps_met", dphi_zleps_met)
            #print("dphi_wleps_met", dphi_wleps_met)
            #print("dr_wl0_wl1", dr_wl0_wl1)
            #print("dr_zl0_zl1", dr_zl0_zl1)
            #print("dr_wleps_zleps", dr_wleps_zleps)
            #print("met.pt", met.pt)
            #print("mt2_val", mt2_val)
            #print("ptl4", ptl4)
            #print("scalarptsum_lepmet", scalarptsum_lepmet)
            #print("scalarptsum_lepmetjet", scalarptsum_lepmetjet)
            #print("z_lep0_pt", z_lep0.pt)
            #print("z_lep1_pt", z_lep1.pt)
            #print("w_lep0_pt", w_lep0.pt)
            #print("w_lep1_pt", w_lep1.pt)
            bdt_vars = [
                ak.fill_none(mll_wl0_wl1,-9999),
                ak.fill_none(dphi_4l_met,-9999),
                ak.fill_none(dphi_zleps_met,-9999),
                ak.fill_none(dphi_wleps_met,-9999),
                ak.fill_none(dr_wl0_wl1,-9999),
                ak.fill_none(dr_zl0_zl1,-9999),
                ak.fill_none(dr_wleps_zleps,-9999),
                ak.fill_none(met.pt,-9999),
                ak.fill_none(mt2_val,-9999),
                ak.fill_none(ptl4,-9999),
                ak.fill_none(scalarptsum_lepmet,-9999),
                ak.fill_none(scalarptsum_lepmetjet,-9999),
                ak.fill_none(z_lep0.pt,-9999),
                ak.fill_none(z_lep1.pt,-9999),
                ak.fill_none(w_lep0.pt,-9999),
                ak.fill_none(w_lep1.pt,-9999),
                #mll_wl0_wl1,
                #dphi_4l_met,
                #dphi_zleps_met,
                #dphi_wleps_met,
                #dr_wl0_wl1,
                #dr_zl0_zl1,
                #dr_wleps_zleps,
                #met.pt,
                #mt2_val,
                #ptl4,
                #scalarptsum_lepmet,
                #scalarptsum_lepmetjet,
                #z_lep0.pt,
                #z_lep1.pt,
                #w_lep0.pt,
                #w_lep1.pt,
            ]

            bdt_of_wwz_raw = es_ec.eval_sig_bdt(events,bdt_vars,ewkcoffea_path(f"data/wwz_zh_bdt/of_WWZ.json"))
            bdt_sf_wwz_raw = es_ec.eval_sig_bdt(events,bdt_vars,ewkcoffea_path(f"data/wwz_zh_bdt/sf_WWZ.json"))
            bdt_of_zh_raw  = es_ec.eval_sig_bdt(events,bdt_vars,ewkcoffea_path(f"data/wwz_zh_bdt/of_ZH.json"))
            bdt_sf_zh_raw  = es_ec.eval_sig_bdt(events,bdt_vars,ewkcoffea_path(f"data/wwz_zh_bdt/sf_ZH.json"))
            bdt_of_wwz = (1.0+math.e**(-bdt_of_wwz_raw))**(-1)
            bdt_sf_wwz = (1.0+math.e**(-bdt_sf_wwz_raw))**(-1)
            bdt_of_zh  = (1.0+math.e**(-bdt_of_zh_raw))**(-1)
            bdt_sf_zh  = (1.0+math.e**(-bdt_sf_zh_raw))**(-1)


            ######### Fill histos #########

            hout = {}

            dense_variables_dict = {
                "mt2" : mt2_val,
                "met" : met.pt,
                "metphi" : met.phi,
                "ptl4" : ptl4,
                "scalarptsum_lep" : scalarptsum_lep,
                "scalarptsum_lepmet" : scalarptsum_lepmet,
                "scalarptsum_lepmetjet" : scalarptsum_lepmetjet,
                "mll_01" : mll_01,
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
                "absdphi_zl0_zl1" : absdphi_zl0_zl1,
                "absdphi_wl0_wl1" : absdphi_wl0_wl1,
                "absdphi_z_ww" : absdphi_z_ww,
                "dphi_4l_met" : dphi_4l_met,
                "dphi_zleps_met" : dphi_zleps_met,
                "dphi_wleps_met" : dphi_wleps_met,

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

                "bdt_of_wwz_raw": bdt_of_wwz_raw,
                "bdt_sf_wwz_raw": bdt_sf_wwz_raw,
                "bdt_of_zh_raw" : bdt_of_zh_raw,
                "bdt_sf_zh_raw" : bdt_sf_zh_raw,
                "bdt_of_wwz": bdt_of_wwz,
                "bdt_sf_wwz": bdt_sf_wwz,
                "bdt_of_zh" : bdt_of_zh,
                "bdt_sf_zh" : bdt_sf_zh,
            }

            # List the hists that are only defined for some categories
            analysis_cats = ["sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C","sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4"]
            exclude_var_dict = {
                "mt2" : ["all_events"],
                "ptl4" : ["all_events"],
                "j0pt" : ["all_events", "4l_presel", "sr_4l_sf_presel", "sr_4l_of_presel", "cr_4l_sf"] + analysis_cats,
                "l0pt" : ["all_events"],
                "mll_01" : ["all_events"],
                "scalarptsum_lep" : ["all_events"],
                "scalarptsum_lepmet" : ["all_events"],
                "scalarptsum_lepmetjet" : ["all_events"],
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
                "absdphi_zl0_zl1" : ["all_events"],
                "absdphi_wl0_wl1" : ["all_events"],
                "absdphi_z_ww" : ["all_events"],
                "dphi_4l_met" : ["all_events"],
                "dphi_zleps_met" : ["all_events"],
                "dphi_wleps_met" : ["all_events"],

                "absdphi_min_afas" : ["all_events"],
                "absdphi_min_afos" : ["all_events"],
                "absdphi_min_sfos" : ["all_events"],
                "mll_min_afas" : ["all_events"],
                "mll_min_afos" : ["all_events"],
                "mll_min_sfos" : ["all_events"],

                "mlb_min" : ["all_events","4l_presel", "sr_4l_sf_presel", "sr_4l_of_presel", "cr_4l_sf"] + analysis_cats,
                "mlb_max" : ["all_events","4l_presel", "sr_4l_sf_presel", "sr_4l_of_presel", "cr_4l_sf"] + analysis_cats,

                "bdt_of_wwz_raw": ["all_events"],
                "bdt_sf_wwz_raw": ["all_events"],
                "bdt_of_zh_raw" : ["all_events"],
                "bdt_sf_zh_raw" : ["all_events"],
                "bdt_of_wwz": ["all_events"],
                "bdt_sf_wwz": ["all_events"],
                "bdt_of_zh" : ["all_events"],
                "bdt_sf_zh" : ["all_events"],
            }


            # Loop over the hists we want to fill
            for dense_axis_name, dense_axis_vals in dense_variables_dict.items():
                #print("\ndense_axis_name,vals",dense_axis_name)
                #print("dense_axis_name,vals",dense_axis_vals)

                # Create the hist for this dense axis variable
                hout[dense_axis_name] = hist.Hist(
                    hist.axis.StrCategory([], growth=True, name="process", label="process"),
                    hist.axis.StrCategory([], growth=True, name="category", label="category"),
                    self._dense_axes_dict[dense_axis_name],
                    storage="weight", # Keeps track of sumw2
                    name="Counts",
                )

                # Decide if we are filling this hist with weight or raw event counts
                if dense_axis_name.endswith("_counts"): weights = events.nom
                #else: weights = weights_obj_base_for_kinematic_syst.partial_weight(include=["norm"])
                else: weights = weights_obj_base_for_kinematic_syst.weight(None)

                # Loop over categories
                for sr_cat in cat_dict["lep_chan_lst"]:

                    # Skip filling if this variable is not relevant for this selection
                    if (dense_axis_name in exclude_var_dict) and (sr_cat in exclude_var_dict[dense_axis_name]): continue

                    # Make the cuts mask
                    cuts_lst = [sr_cat]
                    if isData: cuts_lst.append("is_good_lumi") # Apply golden json requirements if this is data
                    all_cuts_mask = selections.all(*cuts_lst)

                    #run = events.run[all_cuts_mask]
                    #luminosityBlock = events.luminosityBlock[all_cuts_mask]
                    #event = events.event[all_cuts_mask]
                    #w = weights[all_cuts_mask]
                    #if dense_axis_name == "njets":
                    #    print("STARTPRINT")
                    #    for i,j in enumerate(w):
                    #        out_str = f"PRINTTAG {i} {dense_axis_name} {year} {sr_cat} {event[i]} {run[i]} {luminosityBlock[i]} {w[i]}"
                    #        print(out_str,file=sys.stderr,flush=True)
                    #    print("ENDPRINT")
                    #    print("ENDPRINT")
                    #print("\ndense_axis_name",dense_axis_name)
                    #print("sr_cat",sr_cat)
                    #print("dense_axis_vals[all_cuts_mask]",dense_axis_vals[all_cuts_mask])
                    #print("this")

                    # Fill the histos
                    axes_fill_info_dict = {
                        dense_axis_name : dense_axis_vals[all_cuts_mask],
                        "weight"        : weights[all_cuts_mask],
                        "process"       : histAxisName,
                        "category"      : sr_cat,
                        #"systematic"    : "nominal",
                    }
                    hout[dense_axis_name].fill(**axes_fill_info_dict)

        return hout

    def postprocess(self, accumulator):
        return accumulator

