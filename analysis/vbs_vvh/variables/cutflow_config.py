import awkward as ak

cutflow_dict = {
    "MET_objsel": {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "Met_trigger": lambda events, var, obj: events.Pass_MetTriggers == 1,
        "MET_pt": lambda events, var, obj: events.Met_pt > 100,
        "MET_significance": lambda events, var, obj: events.Met_significance > 20,
        "ak4_bveto": lambda events, var, obj: events.Pass_NoAK4BJet == 1,
        "loose_H": lambda events, var, obj: events.HiggsScore > 0.1,
        "loose_V1": lambda events, var, obj: events.V1Score > 0.1,
        "has_vbs": lambda events, var, obj: ak.all(events.vbs_idx_max_Mjj > -1, axis=1),
    },
    "MET_loose": {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "Met_trigger": lambda events, var, obj: events.Pass_MetTriggers == 1,
        "has2GoodFJ": lambda events, var, obj: events.nGoodAK8 >= 2,
        "has_vbs": lambda events, var, obj: ak.all(events.vbs_idx_max_Mjj > -1, axis=1),
        "MET_pt": lambda events, var, obj: events.Met_pt > 100,
        "MET_significance": lambda events, var, obj: events.Met_significance > 100,
        "ak4_bveto": lambda events, var, obj: events.Pass_NoAK4BJet == 1,
        "loose_H": lambda events, var, obj: events.HiggsScore > 0.1,
        "loose_V1": lambda events, var, obj: events.V1Score > 0.1,
    },
    "MET_chan1" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "Met_trigger": lambda events, var, obj: events.Pass_MetTriggers == 1,
        "has2GoodFJ": lambda events, var, obj: events.nGoodAK8 >= 2,
        "has_vbs": lambda events, var, obj: ak.all(events.vbs_idx_max_Mjj > -1, axis=1),
        "MET_pt": lambda events, var, obj: events.Met_pt > 100,
        "MET_significance": lambda events, var, obj: events.Met_significance > 100,
        "ak4_bveto": lambda events, var, obj: events.Pass_NoAK4BJet == 1,
        "loose_H": lambda events, var, obj: events.HiggsScore > 0.1,
        "loose_V1": lambda events, var, obj: events.V1Score > 0.1,
        "medium_H": lambda events, var, obj: events.HiggsScore > 0.5,
        "medium_V1": lambda events, var, obj: events.V1Score > 0.5,
    },
    "MET_chan2" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "Met_trigger": lambda events, var, obj: events.Pass_MetTriggers == 1,
        "has2GoodFJ": lambda events, var, obj: events.nGoodAK8 >= 2,
        "has_vbs": lambda events, var, obj: ak.all(events.vbs_idx_max_Mjj > -1, axis=1),
        "MET_pt": lambda events, var, obj: events.Met_pt > 500,
        "MET_significance": lambda events, var, obj: events.Met_significance > 100,
        "ak4_bveto": lambda events, var, obj: events.Pass_NoAK4BJet == 1,
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        'vbsj_deta_3':lambda events, var, obj: (var["vbsj_deta"]>3),
        "vbsj_Mjj_1000":lambda events, var, obj: (var["vbsj_Mjj"]>1000),
        "V1_m":lambda events, var, obj: ((events.V1_msoftdrop > 50) & (events.V1_msoftdrop<200)),
        "H_m":lambda events, var, obj: ((events.Higgs_msoftdrop > 50) & (events.Higgs_msoftdrop<200)),
    },
    
    "negative_check" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "Met_trigger": lambda events, var, obj: events.Pass_MetTriggers == 1,
        "has2GoodFJ": lambda events, var, obj: events.nGoodAK8 >= 2,
        "has_vbs": lambda events, var, obj: ak.all(events.vbs_idx_max_Mjj > -1, axis=1),
        "MET_pt": lambda events, var, obj: events.Met_pt > 500,
        "MET_significance": lambda events, var, obj: events.Met_significance > 100,
        "ak4_bveto": lambda events, var, obj: events.Pass_NoAK4BJet == 1,
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        'vbsj_deta_3':lambda events, var, obj: (var["vbsj_deta"]>3),
        "vbsj_Mjj_1000":lambda events, var, obj: (var["vbsj_Mjj"]>1000),
        "V1_m":lambda events, var, obj: ((events.V1_msoftdrop > 50) & (events.V1_msoftdrop<200)),
        "H_m":lambda events, var, obj: ((events.Higgs_msoftdrop > 50) & (events.Higgs_msoftdrop<200)),
        "negative_weight":lambda events, var, obj: (events.weight<0),
    },

    "MET_med" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "Met_trigger": lambda events, var, obj: events.Pass_MetTriggers == 1,
        "MET_pt": lambda events, var, obj: events.Met_pt > 100,
        "MET_significance": lambda events, var, obj: events.Met_significance > 20,
        "ak4_bveto": lambda events, var, obj: events.Pass_NoAK4BJet == 1,
        "loose_H": lambda events, var, obj: events.HiggsScore > 0.1,
        "loose_V1": lambda events, var, obj: events.V1Score > 0.1,
        
        "MET_pt_2": lambda events, var, obj: events.Met_pt > 600,
        "MET_significance_2": lambda events, var, obj: events.Met_significance > 60,
        
        'vbsj_deta':lambda events, var, obj: (var["vbsj_deta"]>6),
        "vbsj_Mjj":lambda events, var, obj: (var["vbsj_Mjj"]>1000),

        "H_pt": lambda events, var, obj: events.Higgs_pt > 200,
        "V1_pt": lambda events, var, obj: events.V1_pt > 200,
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        "V1_m":lambda events, var, obj: ((events.V1_msoftdrop > 75)),
        "H_m":lambda events, var, obj: ((events.Higgs_msoftdrop > 85)),
    },
    
    "MET_max" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "Met_trigger": lambda events, var, obj: events.Pass_MetTriggers == 1,
        "MET_pt": lambda events, var, obj: events.Met_pt > 100,
        "MET_significance": lambda events, var, obj: events.Met_significance > 20,
        "ak4_bveto": lambda events, var, obj: events.Pass_NoAK4BJet == 1,
        "loose_H": lambda events, var, obj: events.HiggsScore > 0.1,
        "loose_V1": lambda events, var, obj: events.V1Score > 0.1,
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 1000,
        "MET_significance_tight": lambda events, var, obj: events.Met_significance > 60,
        "H_pt": lambda events, var, obj: events.Higgs_pt > 600,
        "V1_pt": lambda events, var, obj: events.V1_pt > 600,
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,

        'vbsj_deta':lambda events, var, obj: (var["vbsj_deta"]>6),
        "vbsj_Mjj":lambda events, var, obj: (var["vbsj_Mjj"]>1000),

        "V1_m":lambda events, var, obj: ((events.V1_msoftdrop > 75)),
        "H_m":lambda events, var, obj: ((events.Higgs_msoftdrop > 85)),
    },
    "MET_temp" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "Met_trigger": lambda events, var, obj: events.Pass_MetTriggers == 1,
        "MET_pt": lambda events, var, obj: events.Met_pt > 100,
        "MET_significance": lambda events, var, obj: events.Met_significance > 20,
        "ak4_bveto": lambda events, var, obj: events.Pass_NoAK4BJet == 1,
        "loose_H": lambda events, var, obj: events.HiggsScore > 0.1,
        "loose_V1": lambda events, var, obj: events.V1Score > 0.1,
        
        "MET_pt_2": lambda events, var, obj: events.Met_pt > 600,
        "MET_significance_2": lambda events, var, obj: events.Met_significance > 60,

        'vbsj_deta':lambda events, var, obj: (var["vbsj_deta"]>6),
        "vbsj_Mjj":lambda events, var, obj: (var["vbsj_Mjj"]>1000),

        "H_pt": lambda events, var, obj: events.Higgs_pt > 200,
        "V1_pt": lambda events, var, obj: events.V1_pt > 200,
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        "V1_m":lambda events, var, obj: ((events.V1_msoftdrop > 75)),
        "H_m":lambda events, var, obj: ((events.Higgs_msoftdrop > 85)),
        
        #template
        "max_deta_cut":lambda events, var, obj: var["vbsj_deta"] > 8,
        "HV_dR": lambda events, var, obj: var["HV1_dR"] > 8,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 8,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 8,

        "MET_leadAK8_dphi":lambda events, var, obj: var["vbsj_deta"] > 8,
    },

    
    "MET_opti" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
    },
    
    "MET_opti2" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
    },

    "MET_opti3" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
    },

    "MET_opti4" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
        "tight_V1": lambda events, var, obj: events.HiggsScore > 0.9,
    },
    "MET_opti5" : {
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
    },
    "MET_opti5a" : { #because removing tight H increases metric
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,

        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.35,
    },
    "MET_opti6" : { #make more sense to use H but Mjj is better
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.35,

        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
    },
    "MET_opti6a" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.35,

        "vbsj_Mjj":lambda events, var, obj: (var["vbsj_Mjj"]>1000),
    },
    "MET_opti7" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.35,

        "vbsj_Mjj":lambda events, var, obj: (var["vbsj_Mjj"]>1000),
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
    },
    "MET_opti_f" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 700,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.35,

        "vbsj_Mjj":lambda events, var, obj: (var["vbsj_Mjj"]>1000),
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
    },
    "MET_opti_met300" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.35,

        "vbsj_Mjj":lambda events, var, obj: (var["vbsj_Mjj"]>1000),
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
    },
    "MET_opti_met300_2" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),
        
        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.35,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.35,

        "vbsj_Mjj":lambda events, var, obj: (var["vbsj_Mjj"]>1000),
        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
        "sumHT": lambda events, var, obj: var["sum_bosonHT"] > 1350,
    },
    "MET_met300_1" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
    },
    "MET_met300_2" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
    },
    "MET_met300_3" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.6,
    },
    "MET_met300_4" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.6,

        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
    },
    "MET_met300_5a" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.6,

        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
    },
    "MET_met300_5b" : { 
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.6,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.6,

        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
    },
    "MET_met300_6" : {  #from 5b
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        "max_deta_cut":lambda events, var, obj: var["max_deta"] > 6,
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.6,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.6,

        "tight_H": lambda events, var, obj: events.HiggsScore > 0.9,
        "tight_V1": lambda events, var, obj: events.V1Score > 0.9,
    },
    "MET_met300_6a" : {  #adjustment on 6
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        'vbsj_deta':lambda events, var, obj: (var["vbsj_deta"]>6),
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.6,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.6,
    },
    "MET_met300_6b" : {  #adjustment on 6
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        'vbsj_deta':lambda events, var, obj: (var["vbsj_deta"]>6),
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.6,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.6,
        "cut_H": lambda events, var, obj: events.HiggsScore > 0.6,
        "cut_V1": lambda events, var, obj: events.V1Score > 0.6,

    },
    "MET_met300_7" : {  #from 6b
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1) & ak.all(events.vbs_idx_max_Mjj > -1, axis=1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        'vbsj_deta':lambda events, var, obj: (var["vbsj_deta"]>6),
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.6,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.6,
        "cut_H": lambda events, var, obj: events.HiggsScore > 0.6,
        "cut_V1": lambda events, var, obj: events.V1Score > 0.6,

        "vbsj_Mjj_1250":lambda events, var, obj: (var["vbsj_Mjj"]>1250),
    },
    
    "MET_met300_7a" : {  #from 6b
        "all_events": lambda events, var, obj: (events.Pass_MetTriggers == 1) | ~(events.Pass_MetTriggers == 1),
        "objsel": lambda events, var, obj: (events.Pass_MetTriggers == 1) & (events.Met_pt > 100) & (events.Met_significance > 20) & (events.Pass_NoAK4BJet == 1) & (events.HiggsScore > 0.1) & (events.V1Score > 0.1) & ak.all(events.vbs_idx_max_Mjj > -1, axis=1),

        "MET_pt_tight": lambda events, var, obj: events.Met_pt > 300,
        'vbsj_deta':lambda events, var, obj: (var["vbsj_deta"]>6),
        "H_MET_dphi":lambda events, var, obj: var["HMET_dphi"] > 0.6,
        "V1_MET_dphi":lambda events, var, obj: var["V1MET_dphi"] > 0.6,

        "vbsj_Mjj_1250":lambda events, var, obj: (var["vbsj_Mjj"]>1250),
        "cut_V1": lambda events, var, obj: events.V1Score > 0.6,
    },
}