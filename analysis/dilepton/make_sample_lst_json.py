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
    "SingleMuon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022/SingleMuon_Run2022C-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022/DoubleMuon_Run2022C-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,

    "Muon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022/Muon_Run2022C-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "Muon_Run2022D-22Sep2023" : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022/Muon_Run2022D-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,

    "EGamma_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022/EGamma_Run2022C-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "EGamma_Run2022D-22Sep2023" : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022/EGamma_Run2022D-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,

    "MuonEG_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022/MuonEG_Run2022C-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "MuonEG_Run2022D-22Sep2023" : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022/MuonEG_Run2022D-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
}
data_2022EE = {
    "Muon_Run2022E-22Sep2023" : { "Era" : "E" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022EE/Muon_Run2022E-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "Muon_Run2022F-22Sep2023" : { "Era" : "F" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022EE/Muon_Run2022F-22Sep2023-v2_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "Muon_Run2022G-22Sep2023" : { "Era" : "G" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022EE/Muon_Run2022G-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,

    "EGamma_Run2022E-22Sep2023" : { "Era" : "E" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022EE/EGamma_Run2022E-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "EGamma_Run2022F-22Sep2023" : { "Era" : "F" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022EE/EGamma_Run2022F-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "EGamma_Run2022G-22Sep2023" : { "Era" : "G" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022EE/EGamma_Run2022G-22Sep2023-v2_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,

    "MuonEG_Run2022E-22Sep2023" : { "Era" : "E" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022EE/MuonEG_Run2022E-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "MuonEG_Run2022F-22Sep2023" : { "Era" : "F" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022EE/MuonEG_Run2022F-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
    "MuonEG_Run2022G-22Sep2023" : { "Era" : "G" , "histAxisName" : "2022EE_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2022EE/MuonEG_Run2022G-22Sep2023-v1_NANOAOD_Data_Run3_DiLep_20240712" , "xsecName" : "data" , } ,
}
data_2023 = {
    "Muon_Run2023C-22Sep2023"   : {"histAxisName" : "2023_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2023/Muon" , "xsecName" : "data" , } ,
    "EGamma_Run2023C-22Sep2023" : {"histAxisName" : "2023_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2023/EGamma" , "xsecName" : "data" , } ,
    "MuonEG_Run2023C-22Sep2023" : {"histAxisName" : "2023_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2023/MuonEG" , "xsecName" : "data" , } ,
}
data_2023BPix = {
    "Muon_Run2023D-22Sep2023"   : {"histAxisName" : "2023BPix_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2023BPix/Muon" , "xsecName" : "data" , } ,
    "EGamma_Run2023D-22Sep2023" : {"histAxisName" : "2023BPix_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2023BPix/EGamma" , "xsecName" : "data" , } ,
    "MuonEG_Run2023D-22Sep2023" : {"histAxisName" : "2023BPix_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/Data_Run3_DiLep_20240712/2023BPix/MuonEG" , "xsecName" : "data" , } ,
}
############################ Bkg samples ############################
central_2022_bkg_dict = {
    "2022_TTTo2L2Nu"               : { "histAxisName" : "2022_TTTo2L2Nu"               , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "TTTo2L2Nu" , } ,
    "2022_WJetsToLNu"              : { "histAxisName" : "2022_WJetsToLNu"              , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2022_WWTo2L2Nu"               : { "histAxisName" : "2022_WWTo2L2Nu"               , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}
central_2022EE_bkg_dict = {
    "2022EE_TTTo2L2Nu"               : { "histAxisName" : "2022EE_TTTo2L2Nu"           , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022EE/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "TTTo2L2Nu" , } ,
    "2022EE_WJetsToLNu"              : { "histAxisName" : "2022EE_WJetsToLNu"          , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022EE/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2022EE_WWTo2L2Nu"               : { "histAxisName" : "2022EE_WWTo2L2Nu"           , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022EE/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}
central_2023_bkg_dict = {
    "2023_TTTo2L2Nu"               : { "histAxisName" : "2023_TTTo2L2Nu"               , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "TTTo2L2Nu" , } ,
    "2023_WJetsToLNu"              : { "histAxisName" : "2023_WJetsToLNu"              , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2023_WWTo2L2Nu"               : { "histAxisName" : "2023_WWTo2L2Nu"               , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v4_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}
central_2023BPix_bkg_dict = {
    "2023BPix_TTTo2L2Nu"               : { "histAxisName" : "2023BPix_TTTo2L2Nu"       , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023BPix/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "TTTo2L2Nu" , } ,
    "2023BPix_WJetsToLNu"              : { "histAxisName" : "2023BPix_WJetsToLNu"      , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023BPix/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2023BPix_WWTo2L2Nu"               : { "histAxisName" : "2023BPix_WWTo2L2Nu"       , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023BPix/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}

############################ Signal samples ############################

central_2022_sig_dict = {
    "2022_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2022_DYJetsToLL_M_10to50_MLM"         , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2022_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2022_DYJetsToLL_M_50_MLM"             , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}
central_2022EE_sig_dict = {
    "2022EE_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2022EE_DYJetsToLL_M_10to50_MLM"     , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022EE/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2022EE_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2022EE_DYJetsToLL_M_50_MLM"         , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2022EE/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}
central_2023_sig_dict = {
    "2023_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2023_DYJetsToLL_M_10to50_MLM"         , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14_ext1-v3_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2023_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2023_DYJetsToLL_M_50_MLM"             , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v1_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}
central_2023BPix_sig_dict = {
    "2023BPix_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2023BPix_DYJetsToLL_M_10to50_MLM" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023BPix/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2_ext1-v3_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2023BPix_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2023BPix_DYJetsToLL_M_50_MLM"     , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_2L_Skims/Run3_Dilepton_2L_Skims_20240714/MC_Run3_DiLep_20240712/2023BPix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3_NANOAODSIM_MC_Run3_DiLep_20240712", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
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
    out_dir_data_22 = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2022/")
    out_dir_data_22EE = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2022EE/")
    out_dir_data_23 = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2023/")
    out_dir_data_23BPix = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2023BPix/")
    out_dir_bkg = os.path.join(jsons_path,"dilepton_analysis_skims/bkg_samples/")
    out_dir_sig = os.path.join(jsons_path,"dilepton_analysis_skims/sig_samples/")

    # Make configs for data samples
    #make_jsons_for_dict_of_samples(data_2022, "/cmsuf/data/","2022", out_dir_data_22,era_op=1,on_das=False)
    #make_jsons_for_dict_of_samples(data_2022EE, "/cmsuf/data/","2022EE", out_dir_data_22EE,era_op=1,on_das=False)
    #make_jsons_for_dict_of_samples(data_2023, "/cmsuf/data/","2023", out_dir_data_23,on_das=False)
    #make_jsons_for_dict_of_samples(data_2023BPix, "/cmsuf/data/","2023BPix", out_dir_data_23BPix,on_das=False)

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
