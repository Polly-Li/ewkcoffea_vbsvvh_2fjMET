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
#import ewkcoffea.modules.selection_wwz as es_ec
import ewkcoffea.modules.selection_vvh as es_ec
import ewkcoffea.modules.objects_wwz as os_ec
import ewkcoffea.modules.corrections as cor_ec

from topcoffea.modules.get_param_from_jsons import GetParam
get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))


class AnalysisProcessor(processor.ProcessorABC):

    def __init__(self, samples, wc_names_lst=[], hist_lst=None, ecut_threshold=None, do_errors=False, do_systematics=False, skip_obj_systematics=False, split_by_lepton_flavor=False, skip_signal_regions=False, skip_control_regions=False, muonSyst='nominal', dtype=np.float32, siphon_bdt_data=False):

        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype

        # Create the dense axes for the histograms
        self._dense_axes_dict = {
            "met"   : axis.Regular(180, 0, 1000, name="met",  label="met"),
            "metphi": axis.Regular(180, -3.1416, 3.1416, name="metphi", label="met phi"),
            "scalarptsum_jetCentFwd" : axis.Regular(180, 0, 1500, name="scalarptsum_jetCentFwd", label="H_T small radius"),
            "scalarptsum_lep" : axis.Regular(180, 0, 800, name="scalarptsum_lep", label="S_T"),
            "scalarptsum_lepmet" : axis.Regular(180, 0, 1200, name="scalarptsum_lepmet", label="S_T + metpt"),
            "scalarptsum_lepmetFJ" : axis.Regular(180, 0, 1900, name="scalarptsum_lepmetFJ", label="S_T + metpt + FJ pt"),
            "l0_pt"  : axis.Regular(180, 0, 1000, name="l0_pt", label="l0_pt"),
            "l0_eta"  : axis.Regular(180, -3,3, name="l0_eta", label="l0 eta"),

            #"mlb_min" : axis.Regular(180, 0, 300, name="mlb_min",  label="min mass(b+l)"),
            #"mlb_max" : axis.Regular(180, 0, 1000, name="mlb_max",  label="max mass(b+l)"),

            "njets"   : axis.Regular(8, 0, 8, name="njets",   label="Jet multiplicity"),
            "nleps"   : axis.Regular(5, 0, 5, name="nleps",   label="Lep multiplicity"),
            "nbtagsl" : axis.Regular(4, 0, 4, name="nbtagsl", label="Loose btag multiplicity"),
            "nbtagsm" : axis.Regular(4, 0, 4, name="nbtagsm", label="Medium btag multiplicity"),

            "njets_counts"   : axis.Regular(30, 0, 30, name="njets_counts",   label="Jet multiplicity counts (central)"),
            "nleps_counts"   : axis.Regular(30, 0, 30, name="nleps_counts",   label="Lep multiplicity counts (central)"),

            #
            "nfatjets"   : axis.Regular(8, 0, 8, name="nfatjets",   label="Fat jet multiplicity"),
            "njets_forward"   : axis.Regular(8, 0, 8, name="njets_forward",   label="Jet multiplicity (forward)"),
            "njets_tot"   : axis.Regular(8, 0, 8, name="njets_tot",   label="Jet multiplicity (central and forward)"),

            "fj0_pt"  : axis.Regular(180, 0, 2000, name="fj0_pt", label="fj0 pt"),
            "fj0_mass"  : axis.Regular(180, 0, 400, name="fj0_mass", label="fj0 mass"),
            "fj0_eta" : axis.Regular(180, -5, 5, name="fj0_eta", label="fj0 eta"),
            "fj0_phi" : axis.Regular(180, -3.1416, 3.1416, name="fj0_phi", label="j0 phi"),

            "j0central_pt"  : axis.Regular(180, 0, 800, name="j0central_pt", label="j0 pt (central jets)"), # Naming
            "j0central_eta" : axis.Regular(180, -5, 5, name="j0central_eta", label="j0 eta (central jets)"), # Naming
            "j0central_phi" : axis.Regular(180, -3.1416, 3.1416, name="j0central_phi", label="j0 phi (central jets)"), # Naming


            "j0forward_pt"  : axis.Regular(180, 0, 800, name="j0forward_pt", label="j0 pt (forward jets)"),
            "j0forward_eta" : axis.Regular(180, -5, 5, name="j0forward_eta", label="j0 eta (forward jets)"),
            "j0forward_phi" : axis.Regular(180, -3.1416, 3.1416, name="j0forward_phi", label="j0 phi (forward jets)"),

            "j0any_pt"  : axis.Regular(180, 0, 800, name="j0any_pt", label="j0 pt (all regular jets)"),
            "j0any_eta" : axis.Regular(180, -5, 5, name="j0any_eta", label="j0 eta (all regular jets)"),
            "j0any_phi" : axis.Regular(180, -3.1416, 3.1416, name="j0any_phi", label="j0 phi (all regular jets)"),

            "dr_fj0l0" : axis.Regular(180, 0, 6, name="dr_fj0l0", label="dr between FJ and lepton"),
            "dr_j0fwdj1fwd" : axis.Regular(180, 0, 6, name="dr_j0fwdj1fwd", label="dr between leading two forward jets"),
            "dr_j0centj1cent" : axis.Regular(180, 0, 6, name="dr_j0centj1cent", label="dr between leading two central jets"),
            "dr_j0anyj1any" : axis.Regular(180, 0, 6, name="dr_j0anyj1any", label="dr between leading two jets"),

            "mass_j0centj1cent" : axis.Regular(180, 0, 1500, name="mass_j0centj1cent", label="mjj of two leading non-forward jets"),
            "mass_j0fwdj1fwd" : axis.Regular(180, 0, 1500, name="mass_j0fwdj1fwd", label="mjj of two leading forward jets"),
            "mass_j0anyj1any" : axis.Regular(180, 0, 1500, name="mass_j0anyj1any", label="mjj of two leading jets"),

            "mass_b0b1" : axis.Regular(180, 0, 500, name="mass_b0b1", label="mjj of two leading b jets"),

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

        # Adding list accumulators for BDT output variables and weights
        if siphon_bdt_data:
            list_output_names = [
                "list_bdt_of_wwz",
                "list_bdt_of_zh",
                "list_bdt_of_bkg",
                "list_bdt_of_proc",
                "list_bdt_of_wgt",
                "list_bdt_of_evt",
                "list_bdt_sf_wwz",
                "list_bdt_sf_zh",
                "list_bdt_sf_bkg",
                "list_bdt_sf_proc",
                "list_bdt_sf_wgt",
                "list_bdt_sf_evt",
            ]
            for list_output_name in list_output_names:
                dout[list_output_name] = processor.list_accumulator([])
            for of_bdt_var_name in get_ec_param("of_bdt_var_lst"):
                dout[f"list_of_bdt_{of_bdt_var_name}"] = processor.list_accumulator([])
            for sf_bdt_var_name in get_ec_param("sf_bdt_var_lst"):
                dout[f"list_sf_bdt_{sf_bdt_var_name}"] = processor.list_accumulator([])

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

        # Set the booleans
        self._do_errors = do_errors # Whether to calculate and store the w**2 coefficients
        self._do_systematics = do_systematics # Whether to process systematic samples
        self._skip_obj_systematics = skip_obj_systematics # Skip the JEC/JER/MET systematics (even if running with do_systematics on)
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

        # Set a flag for Run3 years
        is2022 = year in ["2022","2022EE"]
        is2023 = year in ["2023","2023BPix"]

        if is2022 or is2023:
            run_tag = "run3"
            com_tag = "13p6TeV"
        elif year in ["2016","2016APV","2017","2018"]:
            run_tag = "run2"
            com_tag = "13TeV"
        else:
            raise Exception(f"ERROR: Unknown year {year}.")

        # Era Needed for all samples
        if isData:
            era = self._samples[json_name]["era"]
        else:
            era = None


        # Get the dataset name (used for duplicate removal) and check to make sure it is an expected name
        # Get name for MC cases too, since "dataset" is passed to overlap removal function in all cases (though it's not actually used in the MC case)
        dataset = json_name.split('_')[0]
        if isData:
            datasets = ["SingleElectron", "EGamma", "MuonEG", "DoubleMuon", "DoubleElectron", "DoubleEG","Muon"]
            if dataset not in datasets:
                raise Exception(f"ERROR: Unexpected dataset name for data file: {dataset}")

        # Initialize objects
        #met  = events.MET
        met  = events.PuppiMET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet
        fatjets = events.FatJet
        if (is2022 or is2023):
            rho = events.Rho.fixedGridRhoFastjetAll
        else:
            rho = events.fixedGridRhoFastjetAll

        # Assigns some original values that will be changed via kinematic corrections
        met["pt_original"] = met.pt
        met["phi_original"] = met.phi
        jets["pt_original"] = jets.pt
        jets["mass_original"] = jets.mass


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

        # Do the object selection for the VVH eleectrons
        ele["is_tight_lep_for_vvh"] = os_ec.is_loose_vvh_ele(ele) & os_ec.is_tight_vvh_ele(ele)

        # Do the object selection for the WWZ muons
        mu["is_tight_lep_for_vvh"] = os_ec.is_loose_vvh_muo(mu) & os_ec.is_tight_vvh_muo(mu)

        # Get tight leptons for WWZ selection
        ele_vvh_t = ele[ele.is_tight_lep_for_vvh]
        mu_vvh_t = mu[mu.is_tight_lep_for_vvh]

        # Attach the lepton SFs to the electron and muons collections
        #cor_ec.AttachElectronSF(ele_wwz_t,year=year)
        #cor_ec.AttachMuonSF(mu_wwz_t,year=year)

        l_vvh_t = ak.with_name(ak.concatenate([ele_vvh_t,mu_vvh_t],axis=1),'PtEtaPhiMCandidate')
        l_vvh_t = l_vvh_t[ak.argsort(l_vvh_t.pt, axis=-1,ascending=False)] # Sort by pt


        # For WWZ
        l_vvh_t_padded = ak.pad_none(l_vvh_t, 4)
        l0 = l_vvh_t_padded[:,0]
        l1 = l_vvh_t_padded[:,1]

        nleps = ak.num(l_vvh_t)

        # Put njets and l_fo_conept_sorted into events and get 4l event selection mask
        events["l_vvh_t"] = l_vvh_t
        #es_ec.add4lmask_vvh(events, year, isData, histAxisName, is2022,is2023)


        ######### Normalization and weights ###########


        # These weights can go outside of the outside sys loop since they do not depend on pt of mu or jets
        # We only calculate these values if not isData
        # Note: add() will generally modify up/down weights, so if these are needed for any reason after this point, we should instead pass copies to add()
        # Note: Here we will to the weights object the SFs that do not depend on any of the forthcoming loops
        weights_obj_base = coffea.analysis_tools.Weights(len(events),storeIndividual=True)
        if not isData:
            if ak.any(events["LHEReweightingWeight"]):
                genw = events["LHEReweightingWeight"][:,60]
            else:
                genw = events["genWeight"]

            # If it's an EFT sample, just take SM piece
            sm_wgt = 1.0
            eft_coeffs = ak.to_numpy(events["EFTfitCoefficients"]) if hasattr(events, "EFTfitCoefficients") else None
            if eft_coeffs is not None:
                sm_wgt = eft_coeffs[:,0]

            # Normalize by (xsec/sow)*genw where genw is 1 for EFT samples
            # Note that for theory systs, will need to multiply by sow/sow_wgtUP to get (xsec/sow_wgtUp)*genw and same for Down
            #lumi = 1000.0*get_tc_param(f"lumi_{year}")
            lumi = 1000.0*59.7 # Just 18
            weights_obj_base.add("norm",(xsec/sow)*genw*lumi*sm_wgt)


            # Scale weights
            cor_tc.AttachPSWeights(events)
            cor_tc.AttachScaleWeights(events)
            # FSR/ISR weights
            # For now only consider variations in the numerator
            #weights_obj_base.add('ps_isr', events.nom, events.ISRUp, events.ISRDown)
            #weights_obj_base.add('ps_fsr', events.nom, events.FSRUp, events.FSRDown)
            # Renorm/fact scale
            #weights_obj_base.add('QCDscale_ren', events.nom, events.renormUp*(sow/sow_renormUp), events.renormDown*(sow/sow_renormDown))
            #weights_obj_base.add('QCDscale_fac', events.nom, events.factUp*(sow/sow_factUp), events.factDown*(sow/sow_factDown))
            if not (is2022 or is2023):
                # Misc other experimental SFs and systs
                weights_obj_base.add('CMS_l1_ecal_prefiring', events.L1PreFiringWeight.Nom,  events.L1PreFiringWeight.Up,  events.L1PreFiringWeight.Dn)
                weights_obj_base.add('CMS_pileup', cor_tc.GetPUSF((events.Pileup.nTrueInt), year), cor_tc.GetPUSF(events.Pileup.nTrueInt, year, 'up'), cor_tc.GetPUSF(events.Pileup.nTrueInt, year, 'down'))
            else:
                weights_obj_base.add("CMS_pileup", cor_ec.run3_pu_attach(events.Pileup,year,"nominal"), cor_ec.run3_pu_attach(events.Pileup,year,"hi"), cor_ec.run3_pu_attach(events.Pileup,year,"lo"))

            # Lepton SFs and systs
            #weights_obj_base.add(f"CMS_eff_m_{com_tag}", events.sf_4l_muon, copy.deepcopy(events.sf_4l_hi_muon), copy.deepcopy(events.sf_4l_lo_muon))
            #weights_obj_base.add(f"CMS_eff_e_{com_tag}", events.sf_4l_elec, copy.deepcopy(events.sf_4l_hi_elec), copy.deepcopy(events.sf_4l_lo_elec))


        # Set up the list of systematics that are handled via event weight variations
        wgt_correction_syst_lst_common = []
        wgt_correction_syst_lst = wgt_correction_syst_lst_common
        #wgt_correction_syst_lst = append_up_down_to_sys_base(wgt_correction_syst_lst)


        ######### The rest of the processor is inside this loop over systs that affect object kinematics  ###########

        obj_correction_systs = [
            #f"CMS_scale_j_{year}",
            #f"CMS_res_j_{year}",
            #f"CMS_scale_met_unclustered_energy_{year}",
        ]
        #obj_correction_systs = append_up_down_to_sys_base(obj_correction_systs)

        # If we're doing systematics and this isn't data, we will loop over the obj correction syst lst list
        if self._do_systematics and not isData and not self._skip_obj_systematics: obj_corr_syst_var_list = ["nominal"] + obj_correction_systs
        # Otherwise loop juse once, for nominal
        else: obj_corr_syst_var_list = ['nominal']

        # Loop over the list of systematic variations (that impact object kinematics) that we've constructed
        for obj_corr_syst_var in obj_corr_syst_var_list:
            # Make a copy of the base weights object, so that each time through the loop we do not double count systs
            # In this loop over systs that impact kinematics, we will add to the weights objects the SFs that depend on the object kinematics
            weights_obj_base_for_kinematic_syst = copy.deepcopy(weights_obj_base)


            #################### Jets ####################

            # Fat jets
            goodfatjets = fatjets[os_ec.is_good_fatjet(fatjets)]
            goodfatjets = os_ec.get_cleaned_collection(l_vvh_t,goodfatjets) # Clean against leps # TODO Should clean lets against FJs?

            # Clean with dr (though another option is to use jetIdx)
            cleanedJets = os_ec.get_cleaned_collection(l_vvh_t,jets) # Clean against leps
            cleanedJets = os_ec.get_cleaned_collection(goodfatjets,cleanedJets) # Clean against fat jets
            jetptname = "pt_nom" if hasattr(cleanedJets, "pt_nom") else "pt"

            # Jet Veto Maps
            # Removes events that have ANY jet in a specific eta-phi space (not required for Run 2)
            # Zero is passing the veto map, so Run 2 will be assigned an array of length events with all zeros
            veto_map_array = cor_ec.ApplyJetVetoMaps(cleanedJets, year) if (is2022 or is2023) else ak.zeros_like(met.pt)
            veto_map_mask = (veto_map_array == 0)

            ##### JME Stuff #####

            # Getting the raw pT and raw mass for jets
            cleanedJets["pt_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.pt_original
            cleanedJets["mass_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.mass_original

            # Getting the generated pT (zeros for unmatched jets)
            # Note this is not used for data, so we use ak.ones_like to create a dummy object
            if not isData:
                cleanedJets["pt_gen"] =ak.values_astype(ak.fill_none(cleanedJets.matched_gen.pt, 0), np.float32)
            else:
                cleanedJets["pt_gen"] =ak.ones_like(cleanedJets.pt)

            # Need to broadcast Rho to have same structure as cleanedJets
            cleanedJets["rho"] = ak.broadcast_arrays(rho, cleanedJets.pt)[0]

            events_cache = events.caches[0] # used for storing intermediary values for corrections
            cleanedJets = cor_ec.ApplyJetCorrections(year,isData, era).build(cleanedJets,lazy_cache=events_cache)
            cleanedJets = cor_ec.ApplyJetSystematics(year,cleanedJets,obj_corr_syst_var)

            # Grab the correctable jets
            correctionJets = os_ec.get_correctable_jets(cleanedJets)

            met = cor_ec.CorrectedMETFactory(correctionJets,year,met,obj_corr_syst_var,isData)
            ##### End of JME #####

            # Selecting jets and cleaning them
            ##cleanedJets["is_good"] = os_tc.is_tight_jet(getattr(cleanedJets, jetptname), cleanedJets.eta, cleanedJets.jetId, pt_cut=20., eta_cut=get_ec_param("wwz_eta_j_cut"), id_cut=get_ec_param("wwz_jet_id_cut"))
            cleanedJets["is_good"] = os_ec.is_good_vbs_jet(cleanedJets,year)
            goodJets = cleanedJets[cleanedJets.is_good & (abs(cleanedJets.eta) <= 2.4)]
            goodJets_forward = cleanedJets[cleanedJets.is_good & (abs(cleanedJets.eta) > 2.4)] # TODO probably not corrected properly

            # Count jets
            njets = ak.num(goodJets)
            njets_forward = ak.num(goodJets_forward)
            njets_tot = njets + njets_forward
            nfatjets = ak.num(goodfatjets)
            ht = ak.sum(goodJets.pt,axis=-1)
            #j0 = goodJets[ak.argmax(goodJets.pt,axis=-1,keepdims=True)]
            #fj0 = goodfatjets[ak.argmax(goodfatjets.pt,axis=-1,keepdims=True)]

            goodJets_ptordered = goodJets[ak.argsort(goodJets.pt,axis=-1,ascending=False)]
            goodJets_ptordered_padded = ak.pad_none(goodJets_ptordered, 2)
            j0 = goodJets_ptordered_padded[:,0]
            j1 = goodJets_ptordered_padded[:,1]

            goodJets_forward_ptordered = goodJets_forward[ak.argsort(goodJets_forward.pt,axis=-1,ascending=False)]
            goodJets_forward_ptordered_padded = ak.pad_none(goodJets_forward_ptordered, 2)
            j0forward = goodJets_forward_ptordered_padded[:,0]
            j1forward = goodJets_forward_ptordered_padded[:,1]

            goodJetsCentFwd = ak.with_name(ak.concatenate([goodJets,goodJets_forward],axis=1),'PtEtaPhiMLorentzVector')
            goodJetsCentFwd_ptordered = goodJetsCentFwd[ak.argsort(goodJetsCentFwd.pt,axis=-1,ascending=False)]
            goodJetsCentFwd_ptordered_padded = ak.pad_none(goodJetsCentFwd_ptordered, 2)
            j0any = goodJetsCentFwd_ptordered_padded[:,0]
            j1any = goodJetsCentFwd_ptordered_padded[:,1]

            goodfatjets_ptordered = goodfatjets[ak.argsort(goodfatjets.pt,axis=-1,ascending=False)]
            goodfatjets_ptordered_padded = ak.pad_none(goodfatjets_ptordered, 2)
            fj0 = goodfatjets_ptordered_padded[:,0]
            fj1 = goodfatjets_ptordered_padded[:,1]

            scalarptsum_jetCentFwd = ak.sum(goodJetsCentFwd.pt,axis=-1)

            # Loose DeepJet WP
            btagger = "btag" # For deep flavor WPs
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

                if not (is2022 or is2023):
                    btag_sf_light = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_incl","central")
                else:
                    btag_sf_light = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_light","central")
                btag_sf_bc = cor_tc.btag_sf_eval(jets_bc,"L",year,"deepJet_comb","central")

                btag_eff_light = cor_ec.btag_eff_eval(jets_light,"L",year)
                btag_eff_bc = cor_ec.btag_eff_eval(jets_bc,"L",year)

                wgt_light = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light, jets_light.btagDeepFlavB>btagwpl)
                wgt_bc    = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc,    jets_bc.btagDeepFlavB>btagwpl)

                wgt_btag_nom = wgt_light*wgt_bc
                ##weights_obj_base_for_kinematic_syst.add("btagSF", wgt_btag_nom)

                # Put the btagging up and down weight variations into the weights object
                if self._do_systematics:

                    # Run3 2022 btagging systematics stuff
                    # Note light correlated and uncorrelated are missing, so just using total, as suggested by the pog
                    # See this for more info: https://cms-talk.web.cern.ch/t/2022-btag-sf-recommendations/42262
                    if (is2022 or is2023):
                        for corr_str in ["correlated", "uncorrelated"]:
                            year_tag = f"_{year}"
                            if corr_str == "correlated": year_tag = ""
                            btag_sf_bc_up      = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"up_{corr_str}")
                            btag_sf_bc_down    = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"down_{corr_str}")
                            wgt_bc_up      = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_up,    jets_bc.btagDeepFlavB>btagwpl)
                            wgt_bc_down    = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_down,    jets_bc.btagDeepFlavB>btagwpl)
                            # Note, up and down weights scaled by 1/wgt_btag_nom so that don't double count the central btag correction (i.e. don't apply it also in the case of up and down variations)
                            ##weights_obj_base_for_kinematic_syst.add(f"CMS_btag_fixedWP_comb_bc_{corr_str}{year_tag}",    events.nom, wgt_light*wgt_bc_up/wgt_btag_nom, wgt_light*wgt_bc_down/wgt_btag_nom)

                        # Light have no correlated/uncorrelated so just use total:
                        btag_sf_light_up   = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_light","up")
                        btag_sf_light_down = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_light","down")
                        wgt_light_up   = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_up, jets_light.btagDeepFlavB>btagwpl)
                        wgt_light_down = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_down, jets_light.btagDeepFlavB>btagwpl)
                        # Note, up and down weights scaled by 1/wgt_btag_nom so that don't double count the central btag correction (i.e. don't apply it also in the case of up and down variations)
                        ##weights_obj_base_for_kinematic_syst.add("CMS_btag_fixedWP_incl_light_correlated", events.nom, wgt_light_up*wgt_bc/wgt_btag_nom, wgt_light_down*wgt_bc/wgt_btag_nom)

                    # Run2 btagging systematics stuff
                    else:
                        for corr_str in ["correlated", "uncorrelated"]:
                            year_tag = f"_{year}"
                            if corr_str == "correlated": year_tag = ""

                            btag_sf_light_up   = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_incl",f"up_{corr_str}")
                            btag_sf_light_down = cor_tc.btag_sf_eval(jets_light, "L",year_light,"deepJet_incl",f"down_{corr_str}")
                            btag_sf_bc_up      = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"up_{corr_str}")
                            btag_sf_bc_down    = cor_tc.btag_sf_eval(jets_bc,    "L",year,      "deepJet_comb",f"down_{corr_str}")

                            wgt_light_up   = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_up, jets_light.btagDeepFlavB>btagwpl)
                            wgt_bc_up      = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_up,    jets_bc.btagDeepFlavB>btagwpl)
                            wgt_light_down = cor_tc.get_method1a_wgt_singlewp(btag_eff_light,btag_sf_light_down, jets_light.btagDeepFlavB>btagwpl)
                            wgt_bc_down    = cor_tc.get_method1a_wgt_singlewp(btag_eff_bc,   btag_sf_bc_down,    jets_bc.btagDeepFlavB>btagwpl)

                            # Note, up and down weights scaled by 1/wgt_btag_nom so that don't double count the central btag correction (i.e. don't apply it also in the case of up and down variations)
                            ##weights_obj_base_for_kinematic_syst.add(f"CMS_btag_fixedWP_incl_light_{corr_str}{year_tag}", events.nom, wgt_light_up*wgt_bc/wgt_btag_nom, wgt_light_down*wgt_bc/wgt_btag_nom)
                            ##weights_obj_base_for_kinematic_syst.add(f"CMS_btag_fixedWP_comb_bc_{corr_str}{year_tag}",    events.nom, wgt_light*wgt_bc_up/wgt_btag_nom, wgt_light*wgt_bc_down/wgt_btag_nom)


            ######### Masks we need for the selection ##########

            # Pass trigger mask
            era_for_trg_check = era
            if not (is2022 or is2023):
                # Era not used for R2
                era_for_trg_check = None
            #pass_trg = es_tc.trg_pass_no_overlap(events,isData,dataset,str(year),dataset_dict=es_ec.dataset_dict,exclude_dict=es_ec.exclude_dict,era=era_for_trg_check)
            #pass_trg = (pass_trg & es_ec.trg_matching(events,year))

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


            ######### VVH event selection stuff #########

            filter_mask = es_ec.get_filter_flag_mask(events,year,is2022,is2023)

            ######### Get variables we haven't already calculated #########

            l0_pt = l0.pt
            l0_eta = l0.eta
            j0central_pt = j0.pt
            j0central_eta = j0.eta
            j0central_phi = j0.phi
            mll_01 = (l0+l1).mass
            scalarptsum_lep = ak.sum(l_vvh_t.pt,axis=-1)
            ##scalarptsum_lep = l0.pt + l1.pt 
            scalarptsum_lepmet = scalarptsum_lep + met.pt
            scalarptsum_lepmetFJ = scalarptsum_lep + met.pt + fj0.pt
            #scalarptsum_lepmetjet = l0.pt + l1.pt + met.pt + ak.sum(goodJets.pt,axis=-1)
            #scalarptsum_jet = ak.sum(goodJets.pt,axis=-1)

            # lb pairs (i.e. always one lep, one bjet)
            bjets = goodJets[isBtagJetsLoose]
            lb_pairs = ak.cartesian({"l":l_vvh_t,"j":bjets})
            mlb_min = ak.min((lb_pairs["l"] + lb_pairs["j"]).mass,axis=-1)
            mlb_max = ak.max((lb_pairs["l"] + lb_pairs["j"]).mass,axis=-1)

            bjets_ptordered = bjets[ak.argsort(bjets.pt,axis=-1,ascending=False)]
            bjets_ptordered_padded = ak.pad_none(bjets_ptordered, 2)
            b0 = bjets_ptordered_padded[:,0]
            b1 = bjets_ptordered_padded[:,1]
            mass_b0b1 = (b0+b1).mass

            # Put the variables into a dictionary for easy access later
            dense_variables_dict = {
                "met" : met.pt,
                "metphi" : met.phi,
                "scalarptsum_lep" : scalarptsum_lep,
                "scalarptsum_jetCentFwd" : scalarptsum_jetCentFwd,
                "scalarptsum_lepmet" : scalarptsum_lepmet,
                "scalarptsum_lepmetFJ" : scalarptsum_lepmetFJ,
                #"scalarptsum_lepmetjet" : scalarptsum_lepmetjet,
                #"scalarptsum_jet" : scalarptsum_jet,
                "l0_pt" : l0_pt,
                "l0_eta" : l0_eta,

                "j0central_pt" : j0central_pt,
                "j0central_eta" : j0central_eta,
                "j0central_phi" : j0central_phi,

                "j0forward_pt" : j0forward.pt,
                "j0forward_eta" : j0forward.eta,
                "j0forward_phi" : j0forward.phi,

                "j0any_pt" : j0any.pt,
                "j0any_eta" : j0any.eta,
                "j0any_phi" : j0any.phi,

                "nleps" : nleps,
                "njets" : njets,
                "nbtagsl" : nbtagsl,

                "nleps_counts" : nleps,
                "njets_counts" : njets,
                "nbtagsl_counts" : nbtagsl,

                ###

                "nbtagsm" : nbtagsm,
                "nbtagsl" : nbtagsl,

                "nfatjets" : nfatjets,
                "njets_forward" : njets_forward,
                "njets_tot" : njets_tot,
                "fj0_pt" : fj0.pt,
                "fj0_mass" : fj0.mass,
                "fj0_eta" : fj0.eta,
                "fj0_phi" : fj0.phi,

                "j0_pt" : j0.pt,
                "j0_eta" : j0.eta,
                "j0_phi" : j0.phi,

                "dr_fj0l0" : fj0.delta_r(l0),
                "dr_j0fwdj1fwd" : j0forward.delta_r(j1forward),
                "dr_j0centj1cent" : j0.delta_r(j1),
                "dr_j0anyj1any" : j0any.delta_r(j1any),

                "mass_j0centj1cent" : (j0+j1).mass,
                "mass_j0fwdj1fwd" : (j0forward+j1forward).mass,
                "mass_j0anyj1any" : (j0any+j1any).mass,

                "mass_b0b1" : mass_b0b1,

            }


            ######### Store boolean masks with PackedSelection ##########

            selections = PackedSelection(dtype='uint64')

            # Lumi mask (for data)
            selections.add("is_good_lumi",lumi_mask)

            zeroj = (njets==0)

            # For Cut Based SRs

            selections.add("all_events", (veto_map_mask | (~veto_map_mask))) # All events.. this logic is a bit roundabout to just get an array of True
            selections.add("filters"                      , veto_map_mask & filter_mask)
            selections.add("exactly1lep"                  , veto_map_mask & filter_mask & (nleps==1))
            selections.add("exactly1lep_exactly1fj"       , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1))
            selections.add("exactly1lep_exactly1fj550"    , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (fj0.pt>550))
            selections.add("exactly1lep_exactly1fj550_2j" , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (fj0.pt>550) & (njets_tot>=2))

            selections.add("exactly1lep_exactly1fj_2j"       , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (njets_tot>=2))

            selections.add("exactly1lep_exactly1fj1100"             , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (fj0.pt>1100))
            selections.add("exactly1lep_exactly1fj800_0j1j"         , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (fj0.pt>800) & (njets_tot<=1))
            selections.add("exactly1lep_exactly1fj700_0jcent1jcent" , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (fj0.pt>700) & (njets<=1))
            selections.add("exactly1lep_exactly1fj700_0j"           , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (fj0.pt>700) & (njets_tot==0))

            selections.add("exactly1lep_exactly1fj_STmet900"    , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (scalarptsum_lepmet>900))
            selections.add("exactly1lep_exactly1fj_STmet1100"    , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (scalarptsum_lepmet>1100))
            selections.add("exactly1lep_exactly1fj_ST600"    , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (scalarptsum_lep>600))

            selections.add("exactly1lep_exactly1fj_STmetFjpt1000"    , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (scalarptsum_lepmetFJ>1000))
            selections.add("exactly1lep_exactly1fj_STmetFjpt1500"    , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (scalarptsum_lepmetFJ>1500))
            selections.add("exactly1lep_exactly1fj_STmetFjpt1700"    , veto_map_mask & filter_mask & (nleps==1) & (nfatjets==1) & (scalarptsum_lepmetFJ>1700))

            cat_dict = {
                "lep_chan_lst" : [
                    "all_events",
                    "filters",
                    "exactly1lep",
                    "exactly1lep_exactly1fj",
                    "exactly1lep_exactly1fj550",
                    "exactly1lep_exactly1fj550_2j",

                    "exactly1lep_exactly1fj_2j",

                    "exactly1lep_exactly1fj1100",
                    "exactly1lep_exactly1fj800_0j1j",
                    "exactly1lep_exactly1fj700_0jcent1jcent",
                    "exactly1lep_exactly1fj700_0j",

                    "exactly1lep_exactly1fj_STmet900",
                    "exactly1lep_exactly1fj_STmet1100",
                    "exactly1lep_exactly1fj_ST600",

                    "exactly1lep_exactly1fj_STmetFjpt1000",
                    "exactly1lep_exactly1fj_STmetFjpt1500",
                    "exactly1lep_exactly1fj_STmetFjpt1700",
                ]
            }


            ######### Fill histos #########

            exclude_var_dict = {
                "mt2" : ["all_events"],
            }

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
                        if (dense_axis_name not in ["njets","njets_counts"]) and (sr_cat in ["all_events","filters","exactly1lep"]): continue # TMP TODO, only fill njets for the cats that dont have all objects

                        # If this is a counts hist, forget the weights and just fill with unit weights
                        if dense_axis_name.endswith("_counts"): weight = events.nom
                        #else: weights = weights_obj_base_for_kinematic_syst.partial_weight(include=["norm"]) # For testing
                        #else: weights = weights_obj_base_for_kinematic_syst.weight(None) # For testing

                        # Make the cuts mask
                        cuts_lst = [sr_cat]
                        if isData: cuts_lst.append("is_good_lumi") # Apply golden json requirements if this is data
                        all_cuts_mask = selections.all(*cuts_lst)

                        # Print info about the events
                        #import sys
                        #run = events.run[all_cuts_mask]
                        #luminosityBlock = events.luminosityBlock[all_cuts_mask]
                        #event = events.event[all_cuts_mask]
                        #w = weight[all_cuts_mask]
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

                        #print("this")
                        #print("\ndense_axis_name",dense_axis_name)
                        #print("\nvals",dense_axis_vals[all_cuts_mask])
                        #print("\ncat",sr_cat)
                        #print("this")

                        # Fill the histos
                        axes_fill_info_dict = {
                            dense_axis_name : ak.fill_none(dense_axis_vals[all_cuts_mask],0),
                            "weight"        : ak.fill_none(weight[all_cuts_mask],0),
                            "process"       : histAxisName,
                            "category"      : sr_cat,
                            "systematic"    : wgt_fluct,
                        }
                        self.accumulator[dense_axis_name].fill(**axes_fill_info_dict)

            # Fill the list accumulator
            if self._siphon_bdt_data:
                for chan,mask in {"of": sr_4l_bdt_of_trn, "sf": sr_4l_bdt_sf_trn}.items():
                    if isData: mask = mask & lumi_mask # Apply golden json to data
                    self.accumulator[f"list_bdt_{chan}_wwz"]  += dense_variables_dict[f"bdt_{chan}_wwz"][mask].to_list()
                    self.accumulator[f"list_bdt_{chan}_zh"]   += dense_variables_dict[f"bdt_{chan}_zh"][mask].to_list()
                    self.accumulator[f"list_bdt_{chan}_bkg"]  += dense_variables_dict[f"bdt_{chan}_bkg"][mask].to_list()
                    self.accumulator[f"list_bdt_{chan}_evt"]  += events.event[mask].to_list()
                    self.accumulator[f"list_bdt_{chan}_wgt"]  += weights_obj_base_for_kinematic_syst.weight(None)[mask]
                    self.accumulator[f"list_bdt_{chan}_proc"] += [histAxisName] * len(dense_variables_dict[f"bdt_{chan}_bkg"][mask])
                    for ivar, var in enumerate(get_ec_param(f"{chan}_bdt_var_lst")):
                        if chan == "of": self.accumulator[f"list_{chan}_bdt_{var}"] += bdt_vars_of_wwz[ivar][mask]
                        if chan == "sf": self.accumulator[f"list_{chan}_bdt_{var}"] += bdt_vars_sf_wwz[ivar][mask]

        return self.accumulator

    def postprocess(self, accumulator):
        return accumulator
