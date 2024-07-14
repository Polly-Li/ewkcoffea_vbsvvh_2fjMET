import os
import subprocess
import topcoffea.modules.sample_lst_jsons_tools as sjt

# Example
#   "UL16APV_"  : { "histAxisName" : "UL16APV_"  , "path" : "" , "xsecName" : ""  } ,
#   "UL16_"    : { "histAxisName" : "UL16_"  , "path" : "" , "xsecName" : ""  } ,
#   "UL17_"    : { "histAxisName" : "UL17_"  , "path" : "" , "xsecName" : ""  } ,
#   "UL18_"    : { "histAxisName" : "UL18_"  , "path" : "" , "xsecName" : ""  } ,


############################ Data samples ############################





data_2022 = {
    "SingleMuon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/SingleMuon_Run2022C-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/DoubleMuon_Run2022C-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,

    "Muon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/Muon_Run2022C-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "Muon_Run2022D-22Sep2023" : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/Muon_Run2022D-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,

    "EGamma_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/EGamma_Run2022C-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "EGamma_Run2022D-22Sep2023" : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/EGamma_Run2022D-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,

    "MuonEG_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/MuonEG_Run2022C-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "MuonEG_Run2022D-22Sep2023" : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/MuonEG_Run2022D-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
}

data_2022EE = {
    "Muon_Run2022E-22Sep2023" : { "Era" : "E" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/Muon_Run2022E-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "Muon_Run2022F-22Sep2023" : { "Era" : "F" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/Muon_Run2022F-22Sep2023-v2_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "Muon_Run2022G-22Sep2023" : { "Era" : "G" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/Muon_Run2022G-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,

    "EGamma_Run2022E-22Sep2023" : { "Era" : "E" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/EGamma_Run2022E-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "EGamma_Run2022F-22Sep2023" : { "Era" : "F" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/EGamma_Run2022F-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "EGamma_Run2022G-22Sep2023" : { "Era" : "G" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/EGamma_Run2022G-22Sep2023-v2_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,

    "MuonEG_Run2022E-22Sep2023" : { "Era" : "E" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/MuonEG_Run2022E-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "MuonEG_Run2022F-22Sep2023" : { "Era" : "F" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/MuonEG_Run2022F-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "MuonEG_Run2022G-22Sep2023" : { "Era" : "G" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/MuonEG_Run2022G-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
}
############################ Bkg samples ############################
central_2022_bkg_dict = {
    "2022_TTTo2L2Nu"               : { "histAxisName" : "2022_TTTo2L2Nu"               , "path" : "", "xsecName" : "TTTo2L2Nu" , } ,
    "2022_WJetsToLNu"              : { "histAxisName" : "2022_WJetsToLNu"              , "path" : "", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2022_WWTo2L2Nu"               : { "histAxisName" : "2022_WWTo2L2Nu"               , "path" : "", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}
central_2022EE_bkg_dict = {
    "2022EE_TTTo2L2Nu"               : { "histAxisName" : "2022EE_TTTo2L2Nu"               , "path" : "", "xsecName" : "TTTo2L2Nu" , } ,
    "2022EE_WJetsToLNu"              : { "histAxisName" : "2022EE_WJetsToLNu"              , "path" : "", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2022EE_WWTo2L2Nu"               : { "histAxisName" : "2022EE_WWTo2L2Nu"               , "path" : "", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}
central_2023_bkg_dict = {
    "2023_TTTo2L2Nu"               : { "histAxisName" : "2023_TTTo2L2Nu"               , "path" : "", "xsecName" : "TTTo2L2Nu" , } ,
    "2023_WJetsToLNu"              : { "histAxisName" : "2023_WJetsToLNu"              , "path" : "", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2023_WWTo2L2Nu"               : { "histAxisName" : "2023_WWTo2L2Nu"               , "path" : "", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}
central_2022EE_bkg_dict = {
    "2023BPix_TTTo2L2Nu"               : { "histAxisName" : "2023BPix_TTTo2L2Nu"               , "path" : "", "xsecName" : "TTTo2L2Nu" , } ,
    "2023BPix_WJetsToLNu"              : { "histAxisName" : "2023BPix_WJetsToLNu"              , "path" : "", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2023BPix_WWTo2L2Nu"               : { "histAxisName" : "2023BPix_WWTo2L2Nu"               , "path" : "", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}

############################ Signal samples ############################

central_2022_sig_dict = {
    "2022_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2022_DYJetsToLL_M_10to50_MLM" , "path" : "", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2022_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2022_DYJetsToLL_M_50_MLM"     , "path" : "", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}
central_2022EE_sig_dict = {
    "2022EE_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2022EE_DYJetsToLL_M_10to50_MLM" , "path" : "", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2022EE_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2022EE_DYJetsToLL_M_50_MLM"     , "path" : "", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}
central_2023_sig_dict = {
    "2023_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2023_DYJetsToLL_M_10to50_MLM" , "path" : "", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2023_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2023_DYJetsToLL_M_50_MLM"     , "path" : "", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}
central_2022EE_sig_dict = {
    "2023BPix_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2023BPix_DYJetsToLL_M_10to50_MLM" , "path" : "", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2023BPix_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2023BPix_DYJetsToLL_M_50_MLM"     , "path" : "", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}


############################ Convenience function ############################

# Convenience function for running sample_lst_jsons_tools make_json() on all entries in a dictionary of samples, and moving the results to out_dir
def make_jsons_for_dict_of_samples(samples_dict,prefix,year,out_dir,era_op = None,on_das=False):
    failed = []
    for sample_name,sample_info in sorted(samples_dict.items()):
        print(f"\n\nMaking JSON for {sample_name}...")
        path = sample_info["path"]
        hist_axis_name = sample_info["histAxisName"]
        xsec_name = sample_info["xsecName"]
        if era_op is not None:
            era = sample_info["Era"]
        else:
            era = None
        print("we are here")
        sjt.make_json(
            sample_dir = path,
            sample_name = sample_name,
            prefix = prefix,
            sample_yr = year,
            xsec_name = xsec_name,
            hist_axis_name = hist_axis_name,
            era = era,
            on_das = on_das,
            include_lhe_wgts_arr = True,
        )
        out_name = sample_name+".json"
        if not os.path.exists(out_name):
            failed.append(sample_name)

        subprocess.run(["mv",out_name,out_dir])

    if len(failed):
        print("Failed:")
        for l in failed:
            print(f"\t{l}")
    else:
        print("Failed: None")


############################ Main ############################

# Uncomment the make_jsons_for_dict_of_samples() lines for the jsons you want to make/remake
def main():

    # Specify output paths
    jsons_path = "../../input_samples/sample_jsons/"
    out_dir_data_22 = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2022")
    out_dir_data_22EE = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2022EE")
    out_dir_data_23 = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2023")
    out_dir_data_23BPix = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2023BPix")
    out_dir_bkg = os.path.join(jsons_path,"dilepton_analysis_skims/bkg_samples/")
    out_dir_sig = os.path.join(jsons_path,"dilepton_analysis_skims/sig_samples/")

    # Make configs for data samples
    make_jsons_for_dict_of_samples(data_2022, "/cmsuf/data/","2022", out_dir_data_22,era_op=1,on_das=False)
    make_jsons_for_dict_of_samples(data_2022EE, "/cmsuf/data/","2022EE", out_dir_data_22EE,era_op=1,on_das=False)
    make_jsons_for_dict_of_samples(data_2023, "/cmsuf/data/","2023", out_dir_data_23,on_das=False)
    make_jsons_for_dict_of_samples(data_2023BPix, "/cmsuf/data/","2023BPix", out_dir_data_23BPix,on_das=False)

    # Make configs for bkg samples
    make_jsons_for_dict_of_samples(central_2022_bkg_dict, "/cmsuf/data/","2022", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2022EE_bkg_dict, "/cmsuf/data/","2022EE", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2023_bkg_dict, "/cmsuf/data/","2023", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2023BPix_bkg_dict, "/cmsuf/data/","2023BPix", out_dir_bkg,on_das=False)

    # Make configs for sig samples
    make_jsons_for_dict_of_samples(central_2022_sig_dict, "/cmsuf/data/","2022", out_dir_sig,on_das=False)
    make_jsons_for_dict_of_samples(central_2022EE_sig_dict, "/cmsuf/data/","2022EE", out_dir_sig,on_das=False)
    make_jsons_for_dict_of_samples(central_2023_sig_dict, "/cmsuf/data/","2023", out_dir_sig,on_das=False)
    make_jsons_for_dict_of_samples(central_2023BPix_sig_dict, "/cmsuf/data/","2023BPix", out_dir_sig,on_das=False)

if __name__ == "__main__":
    main()
