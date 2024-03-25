import coffea
import numpy as np
import awkward as ak
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import processor
import hist
import copy
from hist import axis
from coffea.analysis_tools import PackedSelection
from coffea.lumi_tools import LumiMask
from topcoffea.modules.paths import topcoffea_path
import topcoffea.modules.event_selection as es_tc
from ewkcoffea.modules.paths import ewkcoffea_path as ewkcoffea_path
import ewkcoffea.modules.selection_Run3_2Lep as selRun3_2Lep
import ewkcoffea.modules.corrections as ewk_corrections
import ewkcoffea.modules.objects_Run3_2Lep as objRun3_2Lep
from topcoffea.modules.get_param_from_jsons import GetParam
get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_ec_param = GetParam(ewkcoffea_path("params/params.json"))

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], hist_lst=None, ecut_threshold=None, do_errors=False, do_systematics=False, split_by_lepton_flavor=False, skip_signal_regions=False, skip_control_regions=False, muonSyst='nominal', dtype=np.float32):
        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype
        # Create the histograms
        self._dense_axes_dict = {
            "run": axis.Regular(40, 355800, 356399, name="run", label="run number"),
            "mLL": axis.Regular(180, 0, 160, name="mLL", label="mLL for sf Channel"),
            "pt0": axis.Regular(39, 25, 220, name="pt0", label="Leading Lep pt for sf Channel"),
            "pt1": axis.Regular(180, 0, 160, name="pt1", label="SubLeading Lep pt for sf Channel"),
            "eta0": axis.Regular(180, -2.4, 2.4, name="eta0", label="Leading Lep eta for sf Channel"),
            "eta1": axis.Regular(180, -2.4, 2.4, name="eta1", label="SubLeading Lep eta for sf Channel"),
            "nleps": axis.Regular(20, 0, 20, name="nleps",   label="Lep multiplicity"),
            "njets": axis.Regular(15, 0, 15, name="njets",   label="Jet multiplicity"),
            "nBjets_loose": axis.Regular(15, 0, 15, name="nBjets_loose",   label="Loose B Jet multiplicity"),
            "nBjets_medium": axis.Regular(15, 0, 15, name="nBjets_medium",   label="Medium B Jet multiplicity"),
            "pt_mu": axis.Regular(180,0, 160, name="pt_mu", label="Muon pt in OF Channel"),
            "pt_e": axis.Regular(180,0, 160, name="pt_e", label="Electron pt in OF Channel"),
            "pt_jet0": axis.Regular(180,0, 160, name="pt_jet0", label="Leading Jet pt in OF Channel"),
            "pt_jet1": axis.Regular(180,0, 160, name="pt_jet1", label="Subleading Jet pt in OF Channel"),
            "eta_mu": axis.Regular(180,-2.4, 2.4, name="eta_mu", label="Muon eta in OF Channel"),
            "eta_e": axis.Regular(180,-2.4, 2.4, name="eta_e", label="Electron eta in OF Channel"),
            "eta_jet0": axis.Regular(180,-5, 5, name="eta_jet0", label="Leading Jet eta in OF Channel"),
            "eta_jet1": axis.Regular(180,-5, 5, name="eta_jet1", label="Subleading Jet eta in OF Channel"),
            "reliso0": axis.Regular(180, 0, 0.4, name="reliso0", label="Leading Lep reliso for sf Channel"),
            "reliso1": axis.Regular(180, 0, 0.4, name="reliso1", label="SubLeading Lep reliso for sf Channel"),
            "reliso_mu": axis.Regular(180, 0, 0.4, name="reliso_mu", label="Muon reliso for of Channel"),
            "reliso_e": axis.Regular(180, 0, 0.4, name="reliso_e", label="Electron reliso for of Channel"),
            "dz0": axis.Regular(180, 0, 0.1, name="dz0", label="Leading Lep dz for sf Channel"),
            "dz1": axis.Regular(180, 0, 0.1, name="dz1", label="SubLeading Lep dz for sf Channel"),
            "dz_mu": axis.Regular(180, 0, 0.1, name="dz_mu", label="Muon dz for of Channel"),
            "dz_e": axis.Regular(180, 0, 0.1, name="dz_e", label="Electron dz for of Channel"),
            "dxy0": axis.Regular(180, 0, 0.05, name="dxy0", label="Leading Lep dxy for sf Channel"),
            "dxy1": axis.Regular(180, 0, 0.05, name="dxy1", label="SubLeading Lep dxy for sf Channel"),
            "dxy_mu": axis.Regular(180, 0, 0.05, name="dxy_mu", label="Muon dxy for of Channel"),
            "dxy_e": axis.Regular(180, 0, 0.05, name="dxy_e", label="Electron dxy for of Channel"),
            "met": axis.Regular(180, 0, 160, name="met", label="Missing transverse energy"),
            "phi_met": axis.Regular(180, -3.1416, 3.1416, name="phi_met", label="phi for met"),
            "phi_jet0": axis.Regular(180, -3.1416, 3.1416, name="phi_jet0", label="phi for leading jet"),
            "phi_jet1": axis.Regular(180, -3.1416, 3.1416, name="phi_jet1", label="phi for subleading jet"),
            "phi_mu": axis.Regular(180, -3.1416, 3.1416, name="phi_mu", label="phi for muon"),
            "phi_e": axis.Regular(180, -3.1416, 3.1416, name="phi_e", label="phi for electron"),
            "phi_0": axis.Regular(180, -3.1416, 3.1416, name="phi_0", label="phi for leading sf lep"),
            "phi_1": axis.Regular(180, -3.1416, 3.1416, name="phi_1", label="phi for subleading sf lep"),
            "npvs": axis.Regular(70, 0.5, 70.5, name="npvs",   label="NPVS"),
            "npvsGood": axis.Regular(70, 0.5, 70.5, name="npvsGood",   label="NPVS Good"),
        }

        # Set the list of hists to fill
        if hist_lst is None:
            self._hist_lst = list(self._dense_axes_dict.keys())
        else:
            # Otherwise, just fill the specified subset of hists
            for hist_to_include in hist_lst:
                if hist_to_include not in self._dense_hists_dict.keys():
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
        if isData and ("2022" in year):
            era            = self._samples[dataset]["era"]
        # Get up down weights from input dict
        if (self._do_systematics and not isData):
            if histAxisName in get_ec_param("lo_xsec_samples"):
                # We have a LO xsec for these samples, so for these systs we will have e.g. xsec_LO*(N_pass_up/N_gen_nom)
                # Thus these systs will cover the cross section uncty and the acceptance and effeciency and shape
                # So no NLO rate uncty for xsec should be applied in the text data card
                sow_ISRUp          = self._samples[dataset]["nSumOfWeights"]
                sow_ISRDown        = self._samples[dataset]["nSumOfWeights"]
                sow_FSRUp          = self._samples[dataset]["nSumOfWeights"]
                sow_FSRDown        = self._samples[dataset]["nSumOfWeights"]
                sow_renormUp       = self._samples[dataset]["nSumOfWeights"]
                sow_renormDown     = self._samples[dataset]["nSumOfWeights"]
                sow_factUp         = self._samples[dataset]["nSumOfWeights"]
                sow_factDown       = self._samples[dataset]["nSumOfWeights"]
                sow_renormfactUp   = self._samples[dataset]["nSumOfWeights"]
                sow_renormfactDown = self._samples[dataset]["nSumOfWeights"]
            else:
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
        datasets = ["EGamma", "MuonEG", "DoubleMuon", "SingleMuon","Muon"]
        for d in datasets:
            if d in dataset:
                if d == "Muon" and (("MuonEG" in dataset) or ("SingleMuon" in dataset) or ("DoubleMuon" in dataset)):
                    continue  #
                else:
                    dataset = d

        # Initialize objects
        run  = events.run
        met  = events.MET
        pv  = events.PV
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet
        if not isData:
            pileup = events.Pileup

        # An array of lenght events that is just 1 for each event
        # Probably there's a better way to do this, but we use this method elsewhere so I guess why not..
        events.nom = ak.ones_like(events.MET.pt)

        # Get the lumi mask for data
        if year == "2016" or year == "2016APV":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")
        elif year == "2017":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt")
        elif year == "2018":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt")
        elif year == "2022" or year == "2022EE":
            golden_json_path = topcoffea_path("data/goldenJsons/Cert_Collisions2022_355100_362760_Golden.json")
        else:
            raise ValueError(f"Error: Unknown year \"{year}\".")
        lumi_mask = LumiMask(golden_json_path)(events.run,events.luminosityBlock)

        ################### Object Selections ####################

        # Get the pre-selected electrons and sort by PT
        ele_veto_mask = objRun3_2Lep.is_veto_Run3_2Lep_ele(ele)
        ele["is_veto_lep_for_Run3_2Lep"] = (ele_veto_mask)
        ele_Run3_2Lep_veto = ele[ele.is_veto_lep_for_Run3_2Lep]
        ele_Run3_2Lep_veto = ele_Run3_2Lep_veto[ak.argsort(ele_Run3_2Lep_veto.pt, axis=-1,ascending=False)] # Sort by pt

        # Grab the pre-selected muons and sort by PT
        mu_veto_mask = objRun3_2Lep.is_veto_Run3_2Lep_mu(mu)
        mu["is_veto_lep_for_Run3_2Lep"] = (mu_veto_mask)
        mu_Run3_2Lep_veto = mu[mu.is_veto_lep_for_Run3_2Lep]
        mu_Run3_2Lep_veto = mu_Run3_2Lep_veto[ak.argsort(mu_Run3_2Lep_veto.pt, axis=-1,ascending=False)] # Sort by pt

        #Attach the SF
        ewk_corrections.run3_muons_sf_Attach(mu_Run3_2Lep_veto,year,"NUM_MediumID_DEN_TrackerMuons","NUM_TightPFIso_DEN_MediumID")
        ewk_corrections.run3_electrons_sf_Attach(ele_Run3_2Lep_veto,year,"Medium")

        # Create a List of Leptons from the Muons and Electrons
        l_Run3_2Lep_veto = ak.with_name(ak.concatenate([ele_Run3_2Lep_veto,mu_Run3_2Lep_veto],axis=1),'PtEtaPhiMCandidate')
        l_Run3_2Lep_veto = l_Run3_2Lep_veto[ak.argsort(l_Run3_2Lep_veto.pt, axis=-1,ascending=False)] # Sort by pt

        events["l_Run3_2Lep_veto"] = l_Run3_2Lep_veto
        selRun3_2Lep.add2lmask_Run3_2Lep(events, year, isData)

        #Tighter Ele Selection
        ele_tight_mask = objRun3_2Lep.is_tight_Run3_2Lep_ele(ele_Run3_2Lep_veto)
        ele_Run3_2Lep_veto["is_tight_lep_for_Run3_2Lep"] = (ele_tight_mask)
        ele_Run3_2Lep_tight = ele_Run3_2Lep_veto[ele_Run3_2Lep_veto.is_tight_lep_for_Run3_2Lep]
        ele_Run3_2Lep_tight = ele_Run3_2Lep_tight[ak.argsort(ele_Run3_2Lep_tight.pt, axis=-1,ascending=False)] # Sort by pt

        #Tighter Muon Selection
        mu_tight_mask = objRun3_2Lep.is_tight_Run3_2Lep_mu(mu_Run3_2Lep_veto)
        mu_Run3_2Lep_veto["is_tight_lep_for_Run3_2Lep"] = (mu_tight_mask)
        mu_Run3_2Lep_tight = mu_Run3_2Lep_veto[mu_Run3_2Lep_veto.is_tight_lep_for_Run3_2Lep]
        mu_Run3_2Lep_tight = mu_Run3_2Lep_tight[ak.argsort(mu_Run3_2Lep_tight.pt, axis=-1,ascending=False)] # Sort by pt

        # Create a List of Leptons from the Tight Muons and Electrons
        l_Run3_2Lep_tight = ak.with_name(ak.concatenate([ele_Run3_2Lep_tight,mu_Run3_2Lep_tight],axis=1),'PtEtaPhiMCandidate')
        l_Run3_2Lep_tight = l_Run3_2Lep_tight[ak.argsort(l_Run3_2Lep_tight.pt, axis=-1,ascending=False)] # Sort by pt
        events["l_Run3_2Lep_tight"] = l_Run3_2Lep_tight

        # Get Leps and important values
        l_Run3_2Lep_tight_padded = ak.pad_none(l_Run3_2Lep_tight, 2)
        l0 = l_Run3_2Lep_tight_padded[:,0]
        l1 = l_Run3_2Lep_tight_padded[:,1]
        mll = (l0+l1).mass
        nleps = ak.num(l_Run3_2Lep_tight)



        #################### Jet selection ######################

        # Do the object selection for the Run3 jets
        jets_cleaned_mask = objRun3_2Lep.get_cleaned_collection(l_Run3_2Lep_veto,jets)
        jets["cleaned_jets"] = (jets_cleaned_mask)
        jets_cleaned = jets[jets.cleaned_jets]

        jets_presl_mask = objRun3_2Lep.is_presel_Run3_2Lep_jets(jets_cleaned)
        jets_cleaned["is_jets_for_Run3_2Lep"] = (jets_presl_mask)
        jets_Run3_2Lep = jets_cleaned[jets_cleaned.is_jets_for_Run3_2Lep]
        jets_Run3_2Lep = jets_Run3_2Lep[ak.argsort(jets_Run3_2Lep.pt, axis=-1,ascending=False)] # Sort by pt
        njets = ak.num(jets_Run3_2Lep)

        #Do Some B-Tagging

        #btagwpl = get_tc_param("btag_wp_loose_22EE")
        #btagwpm = get_tc_param("btag_wp_medium_22EE")

        #isBtagJetsLoose = (jets_Run3_2Lep.btagDeepFlavB > btagwpl)
        #isBtagJetsMedium = (jets_Run3_2Lep.btagDeepFlavB > btagwpm)

        #nbtagsl = ak.num(jets_Run3_2Lep[isBtagJetsLoose])
        #nbtagsm = ak.num(jets_Run3_2Lep[isBtagJetsMedium])

        #jets_Run3_2Lep_padded = ak.pad_none(jets_Run3_2Lep, 2)
        #jet0 = jets_Run3_2Lep_padded[:,0]
        #jet1 = jets_Run3_2Lep_padded[:,1]
        #selRun3_2Lep.addjetispresent_Run3_2Lep(jet0, jet1, jets_Run3_2Lep)



        ######### Systematics ###########

        # These weights can go outside of the outside sys loop since they do not depend on pt of mu or jets
        # We only calculate these values if not isData

        weights_obj_base = coffea.analysis_tools.Weights(len(events),storeIndividual=True)
        if not isData:
            sm_wgt = 1.0   #Keep as 1.0 for now
            genw = events["genWeight"]
            lumi = 1000.0*get_tc_param(f"lumi_{year}")
            weights_obj_base.add("norm",((xsec/sow)*genw*lumi*sm_wgt))

        # We do not have systematics yet
        syst_var_list = ['nominal']

        # Loop over the list of systematic variations we've constructed
        for syst_var in syst_var_list:

            events["jets_Run3_2Lep"] = jets_Run3_2Lep
            weights_obj_base_for_kinematics_syst = copy.deepcopy(weights_obj_base)
            if not isData:
                ewk_corrections.run3_pu_Attach(pileup,year)
                weights_obj_base_for_kinematics_syst.add("pu_corr", pileup.pileup_corr)
                weights_obj_base_for_kinematics_syst.add("lepSF_muon", events.muon_sf)
                weights_obj_base_for_kinematics_syst.add("lepSF_ele", events.ele_sf)

            #################### Add variables into event object so that they persist ####################
            #selRun3_2Lep.addjetmask_Run3_2Lep(events, year, isData)
            #selRun3_2Lep.addmetmask_Run3_2Lep(events, year, isData)

            ######### Masks we need for the selection ##########
            # Pass trigger mask
            if isData:
                pass_trg = es_tc.trg_pass_no_overlap(events,isData,dataset,str(year),dataset_dict=selRun3_2Lep.dataset_dict,exclude_dict=selRun3_2Lep.exclude_dict,era=str(era))
            else:
                pass_trg = es_tc.trg_pass_no_overlap(events,isData,dataset,str(year),dataset_dict=selRun3_2Lep.dataset_dict,exclude_dict=selRun3_2Lep.exclude_dict)
            #pass_trg = (pass_trg & selRun3_2Lep.trg_matching(events,year))

            #BTag Mask
            #bmask_atleast1med = (nbtagsm>=1)

            ######### Run3 2Lep event selection stuff #########

            selRun3_2Lep.attach_Run3_2Lep_preselection_mask(events,l_Run3_2Lep_tight_padded[:,0:2])                                              # Attach preselection sf and of flags to the events
            selections = PackedSelection(dtype='uint64')

            # Lumi mask (for data)
            selections.add("is_good_lumi",lumi_mask)

            # For Run3 2Lep selection
            selections.add("2l_sf_ee", (pass_trg & events.is2l & events.Run3_2Lep_presel_sf_ee))
            selections.add("2l_sf_mumu", (pass_trg & events.is2l & events.Run3_2Lep_presel_sf_mumu))
            #selections.add("2l_of", (pass_trg & events.is2l & events.has2jets & events.Run3_2Lep_presel_of))
            #selections.add("2l_of_btag", (bmask_atleast1med & pass_trg & events.is2l & events.Run3_2Lep_presel_of))

            sr_cat_dict = {
                "lep_chan_lst" : ["2l_sf_mumu", "2l_sf_ee"]#, "2l_of", "2l_of_btag"],
            }

            ######### Fill histos #########
            hout = {}

            dense_variables_dict = {
                "2l_sf_mumu" : {
                    "run" : run,
                    "nleps" : nleps,
                    "njets" : njets,
                    #"nBjets_loose" : nbtagsl,
                    #"nBjets_medium" : nbtagsm,
                    "mLL" : mll,
                    "pt0" : l0.pt,
                    "pt1" : l1.pt,
                    "eta0": l0.eta,
                    "eta1": l1.eta,
                    "reliso0": l0.pfRelIso03_all,
                    "reliso1": l1.pfRelIso03_all,
                    "dxy0": l0.dxy,
                    "dxy1": l1.dxy,
                    "dz0": l0.dz,
                    "dz1": l1.dz,
                    "met": met.pt,
                    "phi_met": met.phi,
                    "phi_0": l0.phi,
                    "phi_1": l1.phi,
                    "npvs": pv.npvs,
                    "npvsGood": pv.npvsGood,
                },
                "2l_sf_ee" : {
                    "run" : run,
                    "nleps" : nleps,
                    "njets" : njets,
                    #"nBjets_loose" : nbtagsl,
                    #"nBjets_medium" : nbtagsm,
                    "mLL" : mll,
                    "pt0" : l0.pt,
                    "pt1" : l1.pt,
                    "eta0": l0.eta,
                    "eta1": l1.eta,
                    "reliso0": l0.pfRelIso03_all,
                    "reliso1": l1.pfRelIso03_all,
                    "dxy0": l0.dxy,
                    "dxy1": l1.dxy,
                    "dz0": l0.dz,
                    "dz1": l1.dz,
                    "met": met.pt,
                    "phi_met": met.phi,
                    "phi_0": l0.phi,
                    "phi_1": l1.phi,
                    "npvs": pv.npvs,
                    "npvsGood": pv.npvsGood,
                },
            }

            #First Loop over the SR_CAT_DICT and loop through the relevant SRs
            #Then loop through dense_variables_dict to find correct histos
            for sr_cat_name, sr_cat_lst in sr_cat_dict.items():
                for sr_name in sr_cat_lst:
                    hist_dict = dense_variables_dict[sr_name]
                    for dense_axis_name, dense_axis_vals in hist_dict.items():
                        hist_name = sr_name + "_"+ dense_axis_name
                        hout[hist_name] = hist.Hist(
                            self._dense_axes_dict[dense_axis_name],
                            storage="weight", # Keeps track of sumw2
                        )

                        # Decide if we are filling this hist with weight or raw event counts
                        #if dense_axis_name.endswith("_counts"): weights = events.nom
                        weights = weights_obj_base_for_kinematics_syst.weight(None)

                        # Make the cuts mask
                        cuts_lst = [sr_name]
                        if isData: cuts_lst.append("is_good_lumi") # Apply golden json requirements if this is data
                        all_cuts_mask = selections.all(*cuts_lst)
                        # Fill the histos
                        axes_fill_info_dict = {
                            dense_axis_name : dense_axis_vals[all_cuts_mask],
                            "weight"        : weights[all_cuts_mask],
                        }

                        hout[hist_name].fill(**axes_fill_info_dict)

        return hout

    def postprocess(self, accumulator):
        return accumulator
