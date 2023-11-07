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

        # Create the hist for this dense axis variable
        self._histo_dict = {
            "ptabseta": hist.Hist(
                hist.axis.StrCategory([], growth=True, name="process", label="process"),
                hist.axis.StrCategory([], growth=True, name="flavor", label="flavor"),
                hist.axis.StrCategory([], growth=True, name="tag", label="tag"),
                axis.Variable([0, 20, 30, 60, 120], name="pt",  label="pt"),
                axis.Variable([0, 1, 1.8, 2.4], name="abseta",  label="abseta"),
                storage="weight", # Keeps track of sumw2
                name="Counts",
            )
        }


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

        datasets = ["SingleMuon", "SingleElectron", "EGamma", "MuonEG", "DoubleMuon", "DoubleElectron", "DoubleEG"]
        for d in datasets:
            if d in dataset: dataset = dataset.split('_')[0]

        # Initialize objects
        met  = events.PuppiMET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet


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
        ll_mass_pairs = (llpairs_wwz.l0+llpairs_wwz.l1).mass            # The mll for each ll pair
        mll_min_afos = ak.min(ll_mass_pairs[os_pairs_mask],axis=-1)
        events["min_mll_afos"] = mll_min_afos # Attach this one to the event info since we need it for selection

        # For WWZ
        l_wwz_t_padded = ak.pad_none(l_wwz_t, 4)
        l0 = l_wwz_t_padded[:,0]
        l1 = l_wwz_t_padded[:,1]
        l2 = l_wwz_t_padded[:,2]
        l3 = l_wwz_t_padded[:,3]

        events["l_wwz_t"] = l_wwz_t
        es_ec.add4lmask_wwz(events, year, isData, histAxisName)

        #################### Jets ####################

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

        nbtagsl = ak.num(goodJets[isBtagJetsLoose])


        ######### Masks we need for the event selection ##########

        # Pass trigger mask
        pass_trg = es_tc.trg_pass_no_overlap(events,isData,dataset,str(year),dataset_dict=es_ec.dataset_dict,exclude_dict=es_ec.exclude_dict)
        pass_trg = (pass_trg & es_ec.trg_matching(events,year))

        # Get some preliminary things we'll need
        es_ec.attach_wwz_preselection_mask(events,l_wwz_t_padded[:,0:4]) # Attach preselection sf and of flags to the events
        leps_from_z_candidate_ptordered, leps_not_z_candidate_ptordered = es_ec.get_wwz_candidates(l_wwz_t_padded[:,0:4]) # Get ahold of the leptons from the Z and from the W

        w_lep0 = leps_not_z_candidate_ptordered[:,0]
        w_lep1 = leps_not_z_candidate_ptordered[:,1]
        mll_wl0_wl1 = (w_lep0 + w_lep1).mass

        w_candidates_mll_far_from_z = ak.fill_none(abs(mll_wl0_wl1 - get_ec_param("zmass")) > 10.0,False) # Will enforce this for SF in the PackedSelection

        bmask_exactly0loose = (nbtagsl==0)

        selections = PackedSelection(dtype='uint64')
        selections.add("all_events", (events.is4lWWZ | (~events.is4lWWZ))) # All events.. this logic is a bit roundabout to just get an array of True
        selections.add("4l_presel", (events.is4lWWZ)) # This matches the VVV looper selection (object selection and event selection)
        selections.add("sr_4l_sf_presel", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & (met.pt > 65.0)))
        selections.add("sr_4l_of_presel", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_of))

        selection_mask = selections.all("sr_4l_sf_presel") | selections.all("sr_4l_of_presel")


        ######### Fill histos #########

        hout = self._histo_dict

        dense_axis_name = "ptabseta"

        # Get flat jets for all selected events (lose event level info, we no longer care about event level info)
        #all_cuts_mask = selections.all("sr_4l_sf_presel") | selections.all("sr_4l_of_presel")
        all_cuts_mask = selections.all("all_events")
        jets_sel = ak.flatten(goodJets[all_cuts_mask])
        weights = ak.ones_like(jets_sel.pt)

        pt = jets_sel.pt
        abseta = abs(jets_sel.eta)

        flav_l_mask = np.abs(jets_sel.hadronFlavour) <= 3
        flav_c_mask = np.abs(jets_sel.hadronFlavour) == 4
        flav_b_mask = np.abs(jets_sel.hadronFlavour) == 5
        tag_a_mask = jets_sel.btagDeepFlavB > -999
        tag_l_mask = jets_sel.btagDeepFlavB > btagwpl
        tag_m_mask = jets_sel.btagDeepFlavB > btagwpm

        # Fill the histos
        axes_fill_info_dict = {
            "process"  : histAxisName,
            "flavor"   : "light",
            "tag"      : "medium",
            "weight"   : weights[tag_l_mask],
            "pt"       : pt[tag_l_mask],
            "abseta"   : abseta[tag_l_mask],
        }
        hout[dense_axis_name].fill(**axes_fill_info_dict)

        return hout

    def postprocess(self, accumulator):
        return accumulator

