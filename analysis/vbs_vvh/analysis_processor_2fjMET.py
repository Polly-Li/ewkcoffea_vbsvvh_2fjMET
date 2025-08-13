#!/usr/bin/env python
#import sys
import coffea
import numpy as np
import awkward as ak
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import processor
import hist
from hist import axis
from coffea.analysis_tools import PackedSelection

from coffea.nanoevents.methods import vector
from variables.var_config import obj, dense_variables_config

from utils.config_handling import get_cutflow
from variables.paths import default_output_dir,cutflow_yamls_dir,default_cutflow_yaml 

from paths import objsel_cf as default_objsel_cf
objsel_cf = default_objsel_cf

class AnalysisProcessor(processor.ProcessorABC):

    def __init__(self, samples, wc_names_lst=[], hist_lst=None, ecut_threshold=None, do_errors=False, do_systematics=False, skip_obj_systematics=False, split_by_lepton_flavor=False, skip_signal_regions=False, skip_control_regions=False, siphon_bdt_data=False, n_minus_1=False, one_cut = False, last_cut = False,project=None, cutflow_name=None):
        
        if cutflow_name is None:
            print(f"default cutflow is used.")
            self.cutflow = objsel_cf
        else: 
            print(f'{cutflow_yamls_dir}/{project}.yaml',cutflow_name)
            try: #try project, then try all.yaml, finally just use basic
                self.cutflow = get_cutflow(f'{cutflow_yamls_dir}/{project}.yaml',cutflow_name)
            except:
                try:
                    self.cutflow = get_cutflow(default_cutflow_yaml,cutflow_name)
                    print(f"Cannot get {cutflow_name} in {str(project)}.yaml. Found cutflow in defualt yaml and use it instead")
                except:
                    try:
                        print(get_cutflow(f'{cutflow_yamls_dir}/{project}.yaml'))
                        print("can get project yaml but not the cutflow. Use default cutflow instead.")
                        self.cutflow = objsel_cf
                    except:
                        self.cutflow = objsel_cf
                        print(f"Warning: Cutflow {cutflow_name} not found. Use default cutflow instead.")
            else:
                print(f"Cutflow {cutflow_name} is used.")
        self._n_minus_1 = n_minus_1 #sel mask will be n-1 cut
        self._one_cut = one_cut #sel mask will be one cut only
        self._last_cut = last_cut #only 1 sel mask: all  (but the name will still be the name of last cut only)

        self._samples = samples
        self._wc_names_lst = wc_names_lst #wilson coefficients; not used
        self._dtype = np.float32 #note that this has been changed from function input to this fixed defintion
        

        # Create the dense axes for the histograms
        self._dense_axes_dict = {
            var_name: cfg["axis"]
            for var_name, cfg in dense_variables_config.items()
        }

        # Add histograms to dictionary that will be passed on to dict_accumulator
        dout = {}
        for dense_axis_name in self._dense_axes_dict.keys():
            dout[dense_axis_name] = hist.Hist(
                hist.axis.StrCategory([], growth=True, name="process", label="process"),
                hist.axis.StrCategory([], growth=True, name="category", label="category"),
                hist.axis.StrCategory([], growth=True, name="systematic", label="systematic"),
                hist.axis.StrCategory([], growth=True, name="year", label="year"),
                self._dense_axes_dict[dense_axis_name],
                storage="weight", # Keeps track of sumw2
                name="Counts",
            )

        # Adding list accumulators for BDT output variables and weights
        if siphon_bdt_data:
            list_output_names = []
            for list_output_name in list_output_names:
                dout[list_output_name] = processor.list_accumulator([])

        # Set the accumulator
        self._accumulator = processor.dict_accumulator(dout)
        '''
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
        '''

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

        ak.behavior.update(vector.behavior)

        #move triggers to sparse axes? later?
        Flag_NoAK4BJet = events.Pass_NoAK4BJet

        objects = {
            name: builder(events)
            for name, builder in obj.items()
        }

        var = {
            var_name: cfg["expr"](events, objects)
            for var_name, cfg in dense_variables_config.items()
        }
        ''' #print the firs 5 events' var
        for name, values in var.items():
            print(f"\n{name}:")
            print(ak.to_list(values[:5]))
        '''

        # Dataset parameters
        json_name = events.metadata["dataset"]
        isSig             = self._samples[json_name]["isSig"]
        '''
        isData             = self._samples[json_name]["isData"]
        histAxisName       = self._samples[json_name]["histAxisName"]
        year               = self._samples[json_name]["year"]
        xsec               = self._samples[json_name]["xsec"]
        sow                = self._samples[json_name]["nSumOfWeights"]
        '''

        # An array of lenght events that is just 1 for each event
        # Probably there's a better way to do this, but we use this method elsewhere so I guess why not..
        events.nom = ak.ones_like(var['nGoodAK4'])

        ################### Lepton selection ####################

        ######### Normalization and weights ###########


        # These weights can go outside of the outside sys loop since they do not depend on pt of mu or jets
        # We only calculate these values if not isData
        # Note: add() will generally modify up/down weights, so if these are needed for any reason after this point, we should instead pass copies to add()
        # Note: Here we will to the weights object the SFs that do not depend on any of the forthcoming loops
        # If we're doing systematics and this isn't data, we will loop over the obj correction syst lst list
        obj_corr_syst_var_list = ['nominal']

        # Loop over the list of systematic variations (that impact object kinematics) that we've constructed
        # Make a copy of the base weights object, so that each time through the loop we do not double count systs
        # In this loop over systs that impact kinematics, we will add to the weights objects the SFs that depend on the object kinematics

        #################### Jets ####################
        
        # Put the variables we'll plot into a dictionary for easy access later
        dense_variables_dict  = var

        ######### Store boolean masks with PackedSelection ##########

        selections = PackedSelection(dtype='uint64')
        '''
        for sel, crit_fn_tuple in self.cutflow.items():
            crit_fn = crit_fn_tuple[0]   # get the function
            selections.add(sel, crit_fn(events, var, obj))

        for sel, crit_fn in objsel_cf.items():
            crit_fn = crit_fn_tuple[0]   # get the function
            selections.add(sel, crit_fn(events, var, obj))'''

        for sel, crit_fn in self.cutflow.items():
            mask = crit_fn(events, var, obj)
            selections.add(sel, mask)
            
        for sel, crit_fn in objsel_cf.items():
            mask = crit_fn(events, var, obj)
            selections.add(sel, mask)
        '''
        for sel, crit_fn_tuple in objsel_cf.items():
            crit_fn = crit_fn_tuple[0]                     # extract the function
            mask = crit_fn(events, var, obj)               # evaluate → bool array
            selections.add(sel, mask)
            '''
        ######### Fill histos #########

        exclude_var_dict = {} # Any particular ones to skip

        # Set up the list of weight fluctuations to loop over
        # For now the syst do not depend on the category, so we can figure this out outside of the filling loop
        # Define all weight variations
        wgt_var_dict = {
            "nominal": events.weight,
            #"count"  : events.nom,
        }

        if isSig:                 # or use the ak.firsts trick
            print("Add other coupling points")
            wgt_var_dict["SM"] = events.weight * ak.fill_none(events.LHEReweightingWeight[:,60],0)
            wgt_var_dict["c2v_1p5"] = events.weight * ak.fill_none(events.LHEReweightingWeight[:,72],0)
            # add other cp
        # -----------------------------------------------------------


        # Build cut mask dictionary for each selection
        cut_mask_dict = {}
        if self._n_minus_1:
            # n-1 cutflow: exclude the current cut from the mask
            for sel in self.cutflow:
                exclusive_sels = [k for k in self.cutflow if k != sel] + list(objsel_cf)
                cut_mask_dict[sel] = selections.all(*exclusive_sels)
        elif self._one_cut:
            # one cut only (+ basic cuts)
            for sel in self.cutflow:
                exclusive_sels = [sel] + list(objsel_cf)
                cut_mask_dict[sel] = selections.all(*exclusive_sels)
        elif self._last_cut:
            sels = list(self.cutflow)
            cut_mask_dict[sels[-1]] = selections.all(*sels)
        else:
            # normal cutflow: cumulative cuts
            cumulative_cuts = []
            for sel in self.cutflow:
                cumulative_cuts.append(sel)
                cut_mask_dict[sel] = selections.all(*cumulative_cuts)

        # Produce histograms for each weight variation and each selection
        for wgt_key, wgt in wgt_var_dict.items():
            for sel, sel_mask in cut_mask_dict.items():
                for dense_axis_name, dense_axis_vals in dense_variables_dict.items():
                    axes_fill_info_dict = {
                        dense_axis_name : ak.fill_none(dense_axis_vals[sel_mask], 0),
                        "weight"        : ak.fill_none(wgt[sel_mask], 0),
                        "process"       : ak.fill_none(events.sample_name[sel_mask], "unknown"),
                        "category"      : sel,
                        "systematic"    : wgt_key,
                    }

                    # Include year if needed (can be safely added here)
                    if "year" in self.accumulator[dense_axis_name].axes.name:
                        axes_fill_info_dict["year"] = ak.fill_none(events.sample_year[sel_mask], "unknown")

                    self.accumulator[dense_axis_name].fill(**axes_fill_info_dict)
        return self.accumulator

    def postprocess(self, accumulator):
        return accumulator
