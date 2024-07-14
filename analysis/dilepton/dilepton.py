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
import topcoffea.modules.corrections as cor_tc

from ewkcoffea.modules.paths import ewkcoffea_path as ewkcoffea_path
import ewkcoffea.modules.selection_dilepton as es_ec
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
            "met"   : axis.Regular(180, 0, 300, name="met",  label="met"),
            "metphi": axis.Regular(180, -3.1416, 3.1416, name="metphi", label="met phi"),
            "mll"   : axis.Regular(180, 0, 200, name="mll",  label="mll"),

            "l0pt"  : axis.Regular(180, 0, 500, name="l0pt", label="l0pt"),
            "l1pt"  : axis.Regular(180, 0, 500, name="l1pt", label="l1pt"),

            "l0_eta" : axis.Regular(180, -3, 3, name="lep0_eta", label="Leading Z lep eta"),
            "l1_eta" : axis.Regular(180, -3, 3, name="lep1_eta", label="Subleading Z lep eta"),

            "l0_phi" : axis.Regular(180, -3.1416, 3.1416, name="lep0_phi", label="Leading Z lep phi"),
            "l1_phi" : axis.Regular(180, -3.1416, 3.1416, name="lep1_phi", label="Subleading Z lep phi"),

            "njets"   : axis.Regular(8, 0, 8, name="njets",   label="Jet multiplicity"),
            "nleps"   : axis.Regular(5, 0, 5, name="nleps",   label="Lep multiplicity"),
            "nbtagsl" : axis.Regular(6, 0, 6, name="nbtagsl", label="Loose btag multiplicity"),
            "nbtagsm" : axis.Regular(4, 0, 4, name="nbtagsm", label="Medium btag multiplicity"),

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
        json_name = events.metadata["dataset"]

        isData             = self._samples[json_name]["isData"]
        histAxisName       = self._samples[json_name]["histAxisName"]
        year               = self._samples[json_name]["year"]
        xsec               = self._samples[json_name]["xsec"]
        sow                = self._samples[json_name]["nSumOfWeights"]

        # Set a flag if this is a 2022 year
        is2022 = year in ["2022","2022EE"]
        is2023 = year in ["2023","2023BPix"]

        # If this is a 2022 sample, get the era info
        if isData and is2022:
            era = self._samples[json_name]["era"]
        else:
            era = None

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
        elif year == "2023" or year == "2023BPix":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_Collisions2023_366442_370790_Golden.txt")
        else:
            raise ValueError(f"Error: Unknown year \"{year}\".")
        lumi_mask = LumiMask(golden_json_path)(events.run,events.luminosityBlock)

        ################### Lepton selection ####################

        # Do the object selection for the dilepton electrons (WWZ ID)
        ele_presl_mask = os_ec.is_presel_wwz_ele(ele,is2022,is2023)
        if not (is2022 or is2023):
            ele["topmva"] = os_ec.get_topmva_score_ele(events, year)
            ele["is_tight_lep_for_dilepton"] = ((ele.topmva > get_tc_param("topmva_wp_t_e")) & ele_presl_mask)
        else:
            ele["is_tight_lep_for_dilepton"] = (ele_presl_mask)

        # Do the object selection for the dilepton muons (WWZ ID)
        mu_presl_mask = os_ec.is_presel_wwz_mu(mu,is2022,is2023)
        if not (is2022 or is2023):
            mu["topmva"] = os_ec.get_topmva_score_mu(events, year)
            mu["is_tight_lep_for_wwz"] = ((mu.topmva > get_tc_param("topmva_wp_t_m")) & mu_presl_mask)
        else:
            mu["is_tight_lep_for_wwz"] = (mu_presl_mask)

        # Get tight leptons for WWZ selection
        ele_wwz_t = ele[ele.is_tight_lep_for_wwz]
        mu_wwz_t = mu[mu.is_tight_lep_for_wwz]

        # Attach the lepton SFs to the electron and muons collections
        if (is2022 or is2023):
            cor_ec.run3_muons_sf_attach(mu_wwz_t,year,"NUM_MediumID_DEN_TrackerMuons","NUM_LoosePFIso_DEN_MediumID")
            cor_ec.run3_electrons_sf_attach(ele_wwz_t,year,"wp90iso")
        else:
            cor_ec.AttachElectronSF(ele_wwz_t,year=year)
            cor_ec.AttachMuonSF(mu_wwz_t,year=year)

        l_wwz_t = ak.with_name(ak.concatenate([ele_wwz_t,mu_wwz_t],axis=1),'PtEtaPhiMCandidate')
        l_wwz_t = l_wwz_t[ak.argsort(l_wwz_t.pt, axis=-1,ascending=False)] # Sort by pt

        # For WWZ
        l_wwz_t_padded = ak.pad_none(l_wwz_t, 2)
        l0 = l_wwz_t_padded[:,0]
        l1 = l_wwz_t_padded[:,1]

        nleps = ak.num(l_wwz_t)

        # Put njets and l_fo_conept_sorted into events and get 4l event selection mask
        events["l_wwz_t"] = l_wwz_t
        es_ec.add2lmask_wwz(events,year,isData,histAxisName,is2022,is2023)


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
        wgt_correction_syst_lst = [
            "btagSFlight_correlated", "btagSFbc_correlated", f"btagSFlight_uncorrelated_{year}", f"btagSFbc_uncorrelated_{year}",
            "lepSF_elec", "lepSF_muon", "PU",
            "renorm", "fact", "ISR", "FSR",
        ]
        if not (is2022 or is2023):
            wgt_correction_syst_lst = wgt_correction_syst_lst + ["PreFiring"]
        wgt_correction_syst_lst = append_up_down_to_sys_base(wgt_correction_syst_lst)


        ######### The rest of the processor is inside this loop over systs that affect object kinematics  ###########

        obj_correction_systs = [] # Will have e.g. jes etc

        # Otherwise loop juse once, for nominal
        obj_corr_syst_var_list = ['nominal']

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
            elif year=="2023":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_2023")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_2023")
            elif year=="2023BPix":
                btagwpl = get_tc_param(f"{btagger}_wp_loose_2023BPix")
                btagwpm = get_tc_param(f"{btagger}_wp_medium_2023BPix")
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

            ######### Masks we need for the selection ##########

            # Pass trigger mask
            pass_trg = es_tc.trg_pass_no_overlap(events,isData,dataset,str(year),dataset_dict=es_ec.dataset_dict,exclude_dict=es_ec.exclude_dict,era=era)


            ######### WWZ event selection stuff #########

            # Get some preliminary things we'll need
            es_ec.attach_wwz_preselection_mask(events,l_wwz_t_padded[:,0:4]) # Attach preselection sf and of flags to the events

            # Put the variables into a dictionary for easy access later
            dense_variables_dict = {
                "met" : met.pt,
                "metphi" : met.phi,
                "nleps" : nleps,
                "njets" : njets,
                "nbtagsl" : nbtagsl,
            }


            ######### Store boolean masks with PackedSelection ##########

            selections = PackedSelection(dtype='uint64')

            # Lumi mask (for data)
            selections.add("is_good_lumi",lumi_mask)

            # Fill the packed slection object

            # SRs
            selections.add("mumu_2l_sf",            (pass_trg & events.is4lWWZ & bmask_atleast1loose & events.wwz_presel_of))
            selections.add("ee_2l_sf", (pass_trg & events.is4lWWZ & bmask_atleast1loose & events.wwz_presel_sf & w_candidates_mll_far_from_z & (met.pt > 80.0)))
            selections.add("emu_2l_of", (pass_trg & events.is4lWWZ & bmask_exactly0loose & events.wwz_presel_sf & (~w_candidates_mll_far_from_z)))

            cat_dict = {
                "lep_chan_lst" : [
                    "mumu_2l_sf","ee_2l_sf","emu_2l_of",
                ]
            }


            ######### Fill histos #########

            hout = {}


            # List the hists that are only defined for some categories
            analysis_cats = ["mumu_2l_sf","ee_2l_sf","emu_2l_of"]
            exclude_var_dict = {}

            # Set up the list of weight fluctuations to loop over
            # For now the syst do not depend on the category, so we can figure this out outside of the filling loop
            wgt_var_lst = ["nominal"]

            # Loop over the hists we want to fill
            for dense_axis_name, dense_axis_vals in dense_variables_dict.items():
                if dense_axis_name not in self._hist_lst:
                    continue
                # Create the hist for this dense axis variable
                hout[dense_axis_name] = hist.Hist(
                    hist.axis.StrCategory([], growth=True, name="process", label="process"),
                    hist.axis.StrCategory([], growth=True, name="category", label="category"),
                    hist.axis.StrCategory([], growth=True, name="systematic", label="systematic"),
                    self._dense_axes_dict[dense_axis_name],
                    storage="weight", # Keeps track of sumw2
                    name="Counts",
                )

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

                        # Make the cuts mask
                        cuts_lst = [sr_cat]
                        if isData: cuts_lst.append("is_good_lumi") # Apply golden json requirements if this is data
                        all_cuts_mask = selections.all(*cuts_lst)

                        # Fill the histos
                        axes_fill_info_dict = {
                            dense_axis_name : dense_axis_vals[all_cuts_mask],
                            "weight"        : weight[all_cuts_mask],
                            "process"       : histAxisName,
                            "category"      : sr_cat,
                            "systematic"    : wgt_fluct,
                        }
                        hout[dense_axis_name].fill(**axes_fill_info_dict)

        return hout

    def postprocess(self, accumulator):
        return accumulator
