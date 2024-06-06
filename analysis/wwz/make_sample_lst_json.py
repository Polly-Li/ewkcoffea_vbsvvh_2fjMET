import os
import subprocess
import topcoffea.modules.sample_lst_jsons_tools as sjt

# Example
#   "UL16APV_"  : { "histAxisName" : "UL16APV_"  , "path" : "" , "xsecName" : ""  } ,
#   "UL16_"    : { "histAxisName" : "UL16_"  , "path" : "" , "xsecName" : ""  } ,
#   "UL17_"    : { "histAxisName" : "UL17_"  , "path" : "" , "xsecName" : ""  } ,
#   "UL18_"    : { "histAxisName" : "UL18_"  , "path" : "" , "xsecName" : ""  } ,


############################ Data samples ############################

data_UL16APV = {
    "DoubleEG_Run2016B-ver1_HIPM" : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016B-ver2_HIPM" : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v3_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016C-HIPM"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
    "DoubleEG_Run2016D-HIPM"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
    "DoubleEG_Run2016E-HIPM"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
    "DoubleEG_Run2016F-HIPM"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,

    "MuonEG_Run2016B-ver1_HIPM_UL2016" : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2016B-ver2_HIPM_UL2016" : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2016C-HIPM_UL2016"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
    "MuonEG_Run2016D-HIPM_UL2016"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
    "MuonEG_Run2016E-HIPM_UL2016"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
    "MuonEG_Run2016F-HIPM_UL2016"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,

    "DoubleMuon_Run2016B-ver1_HIPM_UL2016" : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016B-ver2_HIPM_UL2016" : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016C-HIPM_UL2016"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016D-HIPM_UL2016"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016E-HIPM_UL2016"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016F-HIPM_UL2016"      : { "histAxisName" : "UL16APV_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep"      , "xsecName" : "data" , } ,
}

data_UL16 = {
    "DoubleMuon_Run2016F-UL2016" : { "histAxisName" : "UL16_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016G-UL2016" : { "histAxisName" : "UL16_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2016G-UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016H-UL2016" : { "histAxisName" : "UL16_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,

    "DoubleEG_Run2016F-UL2016" : { "histAxisName" : "UL16_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016G-UL2016" : { "histAxisName" : "UL16_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2016G-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016H-UL2016" : { "histAxisName" : "UL16_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,

    "MuonEG_Run2016F-UL2016" : { "histAxisName" : "UL16_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2016G-UL2016" : { "histAxisName" : "UL16_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2016G-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2016H-UL2016" : { "histAxisName" : "UL16_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
}

data_UL17 = {
    "DoubleMuon_Run2017B-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017C-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017D-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017E-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017F-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,

    "DoubleEG_Run2017B-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017C-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017D-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017E-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017F-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleEG_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,

    "MuonEG_Run2017B-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2017C-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2017D-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2017E-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2017F-UL2017" : { "histAxisName" : "UL17_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
}

data_UL18 = {
    "MuonEG_Run2018A-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2018B-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2018C-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "MuonEG_Run2018D-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/MuonEG_Run2018D-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,

    "EGamma_Run2018A-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/EGamma_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "EGamma_Run2018B-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/EGamma_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "EGamma_Run2018C-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/EGamma_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "EGamma_Run2018D-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/EGamma_Run2018D-UL2018_MiniAODv2_NanoAODv9-v3_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,

    "DoubleMuon_Run2018A-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018B-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018C-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018D-UL2018" : { "histAxisName" : "UL18_data" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DoubleMuon_Run2018D-UL2018_MiniAODv2_NanoAODv9-v2_NANOAOD_3LepTau_4Lep" , "xsecName" : "data" , } ,
}

data_2022 = {
    "SingleMuon_Run2022B-22Sep2023" : { "Era" : "B" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/SingleMuon_Run2022B-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "SingleMuon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/SingleMuon_Run2022C-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,

    "DoubleMuon_Run2022B-22Sep2023" : { "Era" : "B" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/DoubleMuon_Run2022B-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/DoubleMuon_Run2022C-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,

    "Muon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/Muon_Run2022C-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "Muon_Run2022D-22Sep2023" : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/Muon_Run2022D-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,

    "EGamma_Run2022B-22Sep2023" : { "Era" : "B" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/EGamma_Run2022B-22Sep2023-v2_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "EGamma_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/EGamma_Run2022C-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
    "EGamma_Run2022D-22Sep2023" : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/EGamma_Run2022D-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,

    "MuonEG_Run2022B-22Sep2023" : { "Era" : "B" , "histAxisName" : "2022_data" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_Data_4L_02162024_Loose/MuonEG_Run2022B-22Sep2023-v1_NANOAOD_Run3_Data_4L_02162024_Loose" , "xsecName" : "data" , } ,
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


# Get some from TOP-22-006 skims #
central_UL16APV_22006_dict = {
    "UL16APV_tWZ4l" : {
        "path" : "/store/user/rucio/kmohrman/skims/mc/new-lepMVA-v2/central_bkgd_p7/TWZToLL/v1/UL16APV_TWZToLL_tlept_Wlept/",
        "histAxisName": "UL16APV_tWZ4l",
        "xsecName": "TWZToLL_tlept_Wlept",
    }
}
central_UL16_22006_dict = {
    "UL16_tWZ4l" : {
        "path" : "/store/user/rucio/kmohrman/skims/mc/new-lepMVA-v2/central_bkgd_p7/TWZToLL/v1/UL16_TWZToLL_tlept_Wlept/",
        "histAxisName": "UL16_tWZ4l",
        "xsecName": "TWZToLL_tlept_Wlept",
    }
}
central_UL17_22006_dict = {
    "UL17_tWZ4l" : {
        "path" : "/store/user/rucio/kmohrman/skims/mc/new-lepMVA-v2/central_bkgd_p7/TWZToLL/v1/UL17_TWZToLL_tlept_Wlept/",
        "histAxisName": "UL17_tWZ4l",
        "xsecName": "TWZToLL_tlept_Wlept",
    }
}
central_UL18_22006_dict = {
    "UL18_tWZ4l" : {
        "path" : "/store/user/rucio/kmohrman/skims/mc/new-lepMVA-v2/central_bkgd_p7/TWZToLL/v1/UL18_TWZToLL_tlept_Wlept/",
        "histAxisName": "UL18_tWZ4l",
        "xsecName": "TWZToLL_tlept_Wlept",
    }
}
central_2022_dict = {
    "2022_tWZ4l" : {
        #"path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TWZ_TtoLNu_WtoLNu_Zto2L_DR2_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3_NANOAODSIM_Run3_MC_4L_04022024_Loose/",
        "path" : "/store/user/t2/users/matthew.dittrich/skims/2022_TWZ_DR1_Loose_v2/TWZ_TtoLNu_WtoLNu_Zto2L_DR1_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v4_NANOAODSIM_2022_TWZ_DR1_Loose_v2/",
        "histAxisName": "2022_tWZ4l",
        "xsecName": "TWZToLL_tlept_Wlept_13p6TeV",
    }
}
central_2022EE_dict = {
    "2022EE_tWZ4l" : {
        #"path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TWZ_TtoLNu_WtoLNu_Zto2L_DR2_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3_NANOAODSIM_Run3_MC_4L_04022024_Loose/",
        "path" : "/store/user/t2/users/matthew.dittrich/skims/2022_TWZ_DR1_Loose_v2/TWZ_TtoLNu_WtoLNu_Zto2L_DR1_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3_NANOAODSIM_2022_TWZ_DR1_Loose_v2/",
        "histAxisName": "2022EE_tWZ4l",
        "xsecName": "TWZToLL_tlept_Wlept_13p6TeV",
    }
}


# All the rest of the backgrounds (located at ucsd)
central_UL16APV_bkg_dict = {

    "UL16APV_ZZTo4L" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ZZTo4L_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL16APV_ZZTo4l",
        "xsecName": "ZZTo4LK",
    },

    "UL16APV_ggToZZTo2e2mu"   : { "histAxisName" : "UL16APV_ggToZZTo2e2mu"    , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"    , "xsecName" : "ggToZZTo2e2muK" , } ,
    "UL16APV_ggToZZTo2e2tau"  : { "histAxisName" : "UL16APV_ggToZZTo2e2tau"   , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"   , "xsecName" : "ggToZZTo2e2tauK" , } ,
    "UL16APV_ggToZZTo2mu2tau" : { "histAxisName" : "UL16APV_ggToZZTo2mu2tau"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"  , "xsecName" : "ggToZZTo2mu2tauK" , } ,
    "UL16APV_ggToZZTo4e"      : { "histAxisName" : "UL16APV_ggToZZTo4e"       , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"       , "xsecName" : "ggToZZTo4eK" , } ,
    "UL16APV_ggToZZTo4mu"     : { "histAxisName" : "UL16APV_ggToZZTo4mu"      , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"      , "xsecName" : "ggToZZTo4muK" , } ,
    "UL16APV_ggToZZTo4tau"    : { "histAxisName" : "UL16APV_ggToZZTo4tau"     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"     , "xsecName" : "ggToZZTo4tauK" , } ,

    "UL16APV_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "UL16APV_DYJetsToLL_M_10to50_MLM" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                        , "xsecName" : "DYJetsToLL_M_10to50_MLM" , } ,
    "UL16APV_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "UL16APV_DYJetsToLL_M_50_MLM"     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                            , "xsecName" : "DYJetsToLL_M_50_MLM" , } ,
    "UL16APV_SSWW"                    : { "histAxisName" : "UL16APV_SSWW"                    , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                                          , "xsecName" : "SSWW" , } ,
    "UL16APV_ST_antitop_t-channel"    : { "histAxisName" : "UL16APV_ST_antitop_t-channel"    , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "ST_antitop_t-channel" , } ,
    "UL16APV_ST_top_s-channel"        : { "histAxisName" : "UL16APV_ST_top_s-channel"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                  , "xsecName" : "ST_top_s-channel" , } ,
    "UL16APV_ST_top_t-channel"        : { "histAxisName" : "UL16APV_ST_top_t-channel"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"     , "xsecName" : "ST_top_t-channel" , } ,
    "UL16APV_TTTo2L2Nu"               : { "histAxisName" : "UL16APV_TTTo2L2Nu"               , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "TTTo2L2Nu" , } ,
    "UL16APV_TTWJetsToLNu"            : { "histAxisName" : "UL16APV_TTWJetsToLNu"            , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"                      , "xsecName" : "TTWJetsToLNu" , } ,
    "UL16APV_TTWJetsToQQ"             : { "histAxisName" : "UL16APV_TTWJetsToQQ"             , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"                       , "xsecName" : "TTWJetsToQQ" , } ,
    "UL16APV_TTZToLLNuNu_M_10"        : { "histAxisName" : "UL16APV_TTZToLLNuNu_M_10"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                              , "xsecName" : "TTZToLLNuNu_M_10" , } ,
    "UL16APV_TTZToLL_M_1to10"         : { "histAxisName" : "UL16APV_TTZToLL_M_1to10"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                               , "xsecName" : "TTZToLL_M_1to10" , } ,
    "UL16APV_TTZToQQ"                 : { "histAxisName" : "UL16APV_TTZToQQ"                 , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "TTZToQQ" , } ,
    "UL16APV_VHnobb"                  : { "histAxisName" : "UL16APV_VHnobb"                  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"                    , "xsecName" : "VHnobb" , } ,
    "UL16APV_WJetsToLNu"              : { "histAxisName" : "UL16APV_WJetsToLNu"              , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep"                                , "xsecName" : "WJetsToLNu" , } ,
    "UL16APV_WWTo2L2Nu"               : { "histAxisName" : "UL16APV_WWTo2L2Nu"               , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "WWTo2L2Nu" , } ,
    "UL16APV_WZTo3LNu"                : { "histAxisName" : "UL16APV_WZTo3LNu"                , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                                  , "xsecName" : "WZTo3LNu" , } ,
    "UL16APV_tW_noFullHad"            : { "histAxisName" : "UL16APV_tW_noFullHad"            , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"              , "xsecName" : "tW_noFullHad" , } ,
    "UL16APV_tZq"                     : { "histAxisName" : "UL16APV_tZq"                     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                             , "xsecName" : "tZq" , } ,
    "UL16APV_tbarW_noFullHad"         : { "histAxisName" : "UL16APV_tbarW_noFullHad"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"          , "xsecName" : "tbarW_noFullHad" , } ,
    "UL16APV_ttHnobb"                 : { "histAxisName" : "UL16APV_ttHnobb"                 , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep"                , "xsecName" : "ttHnobb" , } ,

    "UL16APV_tWll" : { "histAxisName" : "UL16APV_tWll", "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tWll_5f_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "tWll" , } ,
    "UL16APV_WWW"  : { "histAxisName" : "UL16APV_WWW"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "WWW"  } ,
    "UL16APV_WZZ"  : { "histAxisName" : "UL16APV_WZZ"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "WZZ"  } ,
    "UL16APV_ZZZ"  : { "histAxisName" : "UL16APV_ZZZ"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "ZZZ"  } ,
}

central_UL16_bkg_dict = {

    "UL16_ZZTo4L" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ZZTo4L_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL16_ZZTo4l",
        "xsecName": "ZZTo4LK",
    },

    "UL16_ggToZZTo2e2mu"      : { "histAxisName" : "UL16_ggToZZTo2e2mu"       , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"              , "xsecName" : "ggToZZTo2e2muK" , } ,
    "UL16_ggToZZTo2e2tau"     : { "histAxisName" : "UL16_ggToZZTo2e2tau"      , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"             , "xsecName" : "ggToZZTo2e2tauK" , } ,
    "UL16_ggToZZTo2mu2tau"    : { "histAxisName" : "UL16_ggToZZTo2mu2tau"     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"            , "xsecName" : "ggToZZTo2mu2tauK" , } ,
    "UL16_ggToZZTo4e"         : { "histAxisName" : "UL16_ggToZZTo4e"          , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_3LepTau_4Lep"                 , "xsecName" : "ggToZZTo4eK" , } ,
    "UL16_ggToZZTo4mu"        : { "histAxisName" : "UL16_ggToZZTo4mu"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_3LepTau_4Lep"                , "xsecName" : "ggToZZTo4muK" , } ,
    "UL16_ggToZZTo4tau"       : { "histAxisName" : "UL16_ggToZZTo4tau"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"               , "xsecName" : "ggToZZTo4tauK" , } ,

    "UL16_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "UL16_DYJetsToLL_M_10to50_MLM" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                        , "xsecName" : "DYJetsToLL_M_10to50_MLM" , } ,
    "UL16_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "UL16_DYJetsToLL_M_50_MLM"     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                            , "xsecName" : "DYJetsToLL_M_50_MLM" , } ,
    "UL16_SSWW"                    : { "histAxisName" : "UL16_SSWW"                    , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                                          , "xsecName" : "SSWW" , } ,
    "UL16_ST_antitop_t-channel"    : { "histAxisName" : "UL16_ST_antitop_t-channel"    , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "ST_antitop_t-channel" , } ,
    "UL16_ST_top_s-channel"        : { "histAxisName" : "UL16_ST_top_s-channel"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                  , "xsecName" : "ST_top_s-channel" , } ,
    "UL16_ST_top_t-channel"        : { "histAxisName" : "UL16_ST_top_t-channel"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"     , "xsecName" : "ST_top_t-channel" , } ,
    "UL16_TTTo2L2Nu"               : { "histAxisName" : "UL16_TTTo2L2Nu"               , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "TTTo2L2Nu" , } ,
    "UL16_TTWJetsToLNu"            : { "histAxisName" : "UL16_TTWJetsToLNu"            , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                      , "xsecName" : "TTWJetsToLNu" , } ,
    "UL16_TTWJetsToQQ"             : { "histAxisName" : "UL16_TTWJetsToQQ"             , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                       , "xsecName" : "TTWJetsToQQ" , } ,
    "UL16_TTZToLLNuNu_M_10"        : { "histAxisName" : "UL16_TTZToLLNuNu_M_10"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                              , "xsecName" : "TTZToLLNuNu_M_10" , } ,
    "UL16_TTZToLL_M_1to10"         : { "histAxisName" : "UL16_TTZToLL_M_1to10"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                               , "xsecName" : "TTZToLL_M_1to10" , } ,
    "UL16_TTZToQQ"                 : { "histAxisName" : "UL16_TTZToQQ"                 , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "TTZToQQ" , } ,
    "UL16_VHnobb"                  : { "histAxisName" : "UL16_VHnobb"                  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_3LepTau_4Lep"                    , "xsecName" : "VHnobb" , } ,
    "UL16_WJetsToLNu"              : { "histAxisName" : "UL16_WJetsToLNu"              , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_3LepTau_4Lep"                                , "xsecName" : "WJetsToLNu" , } ,
    "UL16_WWTo2L2Nu"               : { "histAxisName" : "UL16_WWTo2L2Nu"               , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "WWTo2L2Nu" , } ,
    "UL16_WZTo3LNu"                : { "histAxisName" : "UL16_WZTo3LNu"                , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                                  , "xsecName" : "WZTo3LNu" , } ,
    "UL16_tW_noFullHad"            : { "histAxisName" : "UL16_tW_noFullHad"            , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"              , "xsecName" : "tW_noFullHad" , } ,
    "UL16_tZq"                     : { "histAxisName" : "UL16_tZq"                     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                             , "xsecName" : "tZq" , } ,
    "UL16_tbarW_noFullHad"         : { "histAxisName" : "UL16_tbarW_noFullHad"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"          , "xsecName" : "tbarW_noFullHad" , } ,
    "UL16_ttHnobb"                 : { "histAxisName" : "UL16_ttHnobb"                 , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep"                , "xsecName" : "ttHnobb" , } ,

    "UL16_tWll"   : { "histAxisName" : "UL16_tWll"   , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tWll_5f_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "tWll" , } ,
    "UL16_WWW"    : { "histAxisName" : "UL16_WWW"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "WWW"  } ,
    "UL16_WZZ"    : { "histAxisName" : "UL16_WZZ"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "WZZ"  } ,
    "UL16_ZZZ"    : { "histAxisName" : "UL16_ZZZ"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "ZZZ"  } ,
}

central_UL17_bkg_dict = {

    "UL17_ZZTo4L" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ZZTo4L_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL17_ZZTo4l",
        "xsecName": "ZZTo4LK",
    },

    "UL17_ggToZZTo2e2mu"      : { "histAxisName" : "UL17_ggToZZTo2e2mu"       , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"                , "xsecName" : "ggToZZTo2e2muK" , } ,
    "UL17_ggToZZTo2e2tau"     : { "histAxisName" : "UL17_ggToZZTo2e2tau"      , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"               , "xsecName" : "ggToZZTo2e2tauK" , } ,
    "UL17_ggToZZTo2mu2tau"    : { "histAxisName" : "UL17_ggToZZTo2mu2tau"     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"              , "xsecName" : "ggToZZTo2mu2tauK" , } ,
    "UL17_ggToZZTo4e"         : { "histAxisName" : "UL17_ggToZZTo4e"          , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"                   , "xsecName" : "ggToZZTo4eK" , } ,
    "UL17_ggToZZTo4mu"        : { "histAxisName" : "UL17_ggToZZTo4mu"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"                  , "xsecName" : "ggToZZTo4muK" , } ,
    "UL17_ggToZZTo4tau"       : { "histAxisName" : "UL17_ggToZZTo4tau"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"                 , "xsecName" : "ggToZZTo4tauK" , } ,

    "UL17_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "UL17_DYJetsToLL_M_10to50_MLM" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                        , "xsecName" : "DYJetsToLL_M_10to50_MLM" , } ,
    "UL17_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "UL17_DYJetsToLL_M_50_MLM"     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                            , "xsecName" : "DYJetsToLL_M_50_MLM" , } ,
    "UL17_SSWW"                    : { "histAxisName" : "UL17_SSWW"                    , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                                          , "xsecName" : "SSWW" , } ,
    "UL17_ST_antitop_t-channel"    : { "histAxisName" : "UL17_ST_antitop_t-channel"    , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "ST_antitop_t-channel" , } ,
    "UL17_ST_top_s-channel"        : { "histAxisName" : "UL17_ST_top_s-channel"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                  , "xsecName" : "ST_top_s-channel" , } ,
    "UL17_ST_top_t-channel"        : { "histAxisName" : "UL17_ST_top_t-channel"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"     , "xsecName" : "ST_top_t-channel" , } ,
    "UL17_TTTo2L2Nu"               : { "histAxisName" : "UL17_TTTo2L2Nu"               , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "TTTo2L2Nu" , } ,
    "UL17_TTWJetsToLNu"            : { "histAxisName" : "UL17_TTWJetsToLNu"            , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                      , "xsecName" : "TTWJetsToLNu" , } ,
    "UL17_TTWJetsToQQ"             : { "histAxisName" : "UL17_TTWJetsToQQ"             , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                       , "xsecName" : "TTWJetsToQQ" , } ,
    "UL17_TTZToLLNuNu_M_10"        : { "histAxisName" : "UL17_TTZToLLNuNu_M_10"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                              , "xsecName" : "TTZToLLNuNu_M_10" , } ,
    "UL17_TTZToLL_M_1to10"         : { "histAxisName" : "UL17_TTZToLL_M_1to10"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                               , "xsecName" : "TTZToLL_M_1to10" , } ,
    "UL17_TTZToQQ"                 : { "histAxisName" : "UL17_TTZToQQ"                 , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "TTZToQQ" , } ,
    "UL17_VHnobb"                  : { "histAxisName" : "UL17_VHnobb"                  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"                    , "xsecName" : "VHnobb" , } ,
    "UL17_WJetsToLNu"              : { "histAxisName" : "UL17_WJetsToLNu"              , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"                                , "xsecName" : "WJetsToLNu" , } ,
    "UL17_WWTo2L2Nu"               : { "histAxisName" : "UL17_WWTo2L2Nu"               , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "WWTo2L2Nu" , } ,
    "UL17_WZTo3LNu"                : { "histAxisName" : "UL17_WZTo3LNu"                , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep"                                  , "xsecName" : "WZTo3LNu" , } ,
    "UL17_tW_noFullHad"            : { "histAxisName" : "UL17_tW_noFullHad"            , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"              , "xsecName" : "tW_noFullHad" , } ,
    "UL17_tZq"                     : { "histAxisName" : "UL17_tZq"                     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                             , "xsecName" : "tZq" , } ,
    "UL17_tbarW_noFullHad"         : { "histAxisName" : "UL17_tbarW_noFullHad"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"          , "xsecName" : "tbarW_noFullHad" , } ,
    "UL17_ttHnobb"                 : { "histAxisName" : "UL17_ttHnobb"                 , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep"                , "xsecName" : "ttHnobb" , } ,

    "UL17_tWll"   : { "histAxisName" : "UL17_tWll"   , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tWll_5f_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "tWll" , } ,
    "UL17_WWW"    : { "histAxisName" : "UL17_WWW"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "WWW"  } ,
    "UL17_WZZ"    : { "histAxisName" : "UL17_WZZ"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "WZZ"  } ,
    "UL17_ZZZ"    : { "histAxisName" : "UL17_ZZZ"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "ZZZ"  } ,
}

central_UL18_bkg_dict = {

    "UL18_ZZTo4L" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ZZTo4L_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL18_ZZTo4l",
        "xsecName": "ZZTo4LK",
    },

    "UL18_ggToZZTo2e2mu"      : { "histAxisName" : "UL18_ggToZZTo2e2mu"       , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"     , "xsecName" : "ggToZZTo2e2muK" , } ,
    "UL18_ggToZZTo2e2tau"     : { "histAxisName" : "UL18_ggToZZTo2e2tau"      , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"    , "xsecName" : "ggToZZTo2e2tauK" , } ,
    "UL18_ggToZZTo2mu2tau"    : { "histAxisName" : "UL18_ggToZZTo2mu2tau"     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"   , "xsecName" : "ggToZZTo2mu2tauK" , } ,
    "UL18_ggToZZTo4e"         : { "histAxisName" : "UL18_ggToZZTo4e"          , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"        , "xsecName" : "ggToZZTo4eK" , } ,
    "UL18_ggToZZTo4mu"        : { "histAxisName" : "UL18_ggToZZTo4mu"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"       , "xsecName" : "ggToZZTo4muK" , } ,
    "UL18_ggToZZTo4tau"       : { "histAxisName" : "UL18_ggToZZTo4tau"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"      , "xsecName" : "ggToZZTo4tauK" , } ,

    "UL18_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "UL18_DYJetsToLL_M_10to50_MLM" , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                        , "xsecName" : "DYJetsToLL_M_10to50_MLM" , } ,
    "UL18_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "UL18_DYJetsToLL_M_50_MLM"     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                            , "xsecName" : "DYJetsToLL_M_50_MLM" , } ,
    "UL18_SSWW"                    : { "histAxisName" : "UL18_SSWW"                    , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                                          , "xsecName" : "SSWW" , } ,
    "UL18_ST_antitop_t-channel"    : { "histAxisName" : "UL18_ST_antitop_t-channel"    , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "ST_antitop_t-channel" , } ,
    "UL18_ST_top_s-channel"        : { "histAxisName" : "UL18_ST_top_s-channel"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                  , "xsecName" : "ST_top_s-channel" , } ,
    "UL18_ST_top_t-channel"        : { "histAxisName" : "UL18_ST_top_t-channel"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"     , "xsecName" : "ST_top_t-channel" , } ,
    "UL18_TTTo2L2Nu"               : { "histAxisName" : "UL18_TTTo2L2Nu"               , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "TTTo2L2Nu" , } ,
    "UL18_TTWJetsToLNu"            : { "histAxisName" : "UL18_TTWJetsToLNu"            , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                      , "xsecName" : "TTWJetsToLNu" , } ,
    "UL18_TTWJetsToQQ"             : { "histAxisName" : "UL18_TTWJetsToQQ"             , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                       , "xsecName" : "TTWJetsToQQ" , } ,
    "UL18_TTZToLLNuNu_M_10"        : { "histAxisName" : "UL18_TTZToLLNuNu_M_10"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                              , "xsecName" : "TTZToLLNuNu_M_10" , } ,
    "UL18_TTZToLL_M_1to10"         : { "histAxisName" : "UL18_TTZToLL_M_1to10"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                               , "xsecName" : "TTZToLL_M_1to10" , } ,
    "UL18_TTZToQQ"                 : { "histAxisName" : "UL18_TTZToQQ"                 , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "TTZToQQ" , } ,
    "UL18_VHnobb"                  : { "histAxisName" : "UL18_VHnobb"                  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"                    , "xsecName" : "VHnobb" , } ,
    "UL18_WJetsToLNu"              : { "histAxisName" : "UL18_WJetsToLNu"              , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"                                , "xsecName" : "WJetsToLNu" , } ,
    "UL18_WWTo2L2Nu"               : { "histAxisName" : "UL18_WWTo2L2Nu"               , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"                                       , "xsecName" : "WWTo2L2Nu" , } ,
    "UL18_WZTo3LNu"                : { "histAxisName" : "UL18_WZTo3LNu"                , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep"                                  , "xsecName" : "WZTo3LNu" , } ,
    "UL18_tW_noFullHad"            : { "histAxisName" : "UL18_tW_noFullHad"            , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"              , "xsecName" : "tW_noFullHad" , } ,
    "UL18_tZq"                     : { "histAxisName" : "UL18_tZq"                     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                             , "xsecName" : "tZq" , } ,
    "UL18_tbarW_noFullHad"         : { "histAxisName" : "UL18_tbarW_noFullHad"         , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"          , "xsecName" : "tbarW_noFullHad" , } ,
    "UL18_ttHnobb"                 : { "histAxisName" : "UL18_ttHnobb"                 , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep"                , "xsecName" : "ttHnobb" , } ,

    "UL18_tWll"   : { "histAxisName" : "UL18_tWll"   , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ST_tWll_5f_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "tWll" , } ,
    "UL18_WWW"    : { "histAxisName" : "UL18_WWW"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "WWW"  } ,
    "UL18_WZZ"    : { "histAxisName" : "UL18_WZZ"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "WZZ"  } ,
    "UL18_ZZZ"    : { "histAxisName" : "UL18_ZZZ"  , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "ZZZ"  } ,
}

central_2022_bkg_dict = {

    "2022_ZZTo4L" : {
        "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose",
        "histAxisName": "2022_ZZTo4l",
        "xsecName": "ZZto4L_13p6TeV",
    },

    "2022_ggToZZTo2e2mu"      : { "histAxisName" : "2022_ggToZZTo2e2mu"       , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGlutoContinto2Zto2E2Mu_TuneCP5_13p6TeV_mcfm701-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v1_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGluToContinto2Zto2E2Mu_13p6TeV", } ,
    "2022_ggToZZTo2e2tau"     : { "histAxisName" : "2022_ggToZZTo2e2tau"      , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGluToContinto2Zto2E2Tau_TuneCP5_13p6TeV_mcfm701-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGluToContinto2Zto2E2Tau_13p6TeV", } ,
    "2022_ggToZZTo2mu2tau"    : { "histAxisName" : "2022_ggToZZTo2mu2tau"     , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGluToContinto2Zto2Mu2Tau_TuneCP5_13p6TeV_mcfm701-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGluToContinto2Zto2Mu2Tau_13p6TeV", } ,
    "2022_ggToZZTo4e"         : { "histAxisName" : "2022_ggToZZTo4e"          , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGlutoContinto2Zto4E_TuneCP5_13p6TeV_mcfm-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGlutoContinto2Zto4E_13p6TeV" , } ,
    "2022_ggToZZTo4mu"        : { "histAxisName" : "2022_ggToZZTo4mu"         , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGlutoContinto2Zto4Mu_TuneCP5_13p6TeV_mcfm-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGlutoContinto2Zto4Mu_13p6TeV" , } ,
    "2022_ggToZZTo4tau"       : { "histAxisName" : "2022_ggToZZTo4tau"        , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGlutoContinto2Zto4Tau_TuneCP5_13p6TeV_mcfm-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGlutoContinto2Zto4Tau_13p6TeV" , } ,

    "2022_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2022_DYJetsToLL_M_10to50_MLM" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2022_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2022_DYJetsToLL_M_50_MLM"     , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v1_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
    #"2022_SSWW"                    : { "histAxisName" : "2022_SSWW"                    , "path" : "", "xsecName" : "" , } ,
    #"2022_ST_antitop_t-channel"    : { "histAxisName" : "2022_ST_antitop_t-channel"    , "path" : "", "xsecName" : "" , } ,
    #"2022_ST_top_s-channel"        : { "histAxisName" : "2022_ST_top_s-channel"        , "path" : "", "xsecName" : "" , } ,
    #"2022_ST_top_t-channel"        : { "histAxisName" : "2022_ST_top_t-channel"        , "path" : "", "xsecName" : "" , } ,
    #"2022_TTTo2L2Nu"               : { "histAxisName" : "2022_TTTo2L2Nu"               , "path" : "", "xsecName" : "" , } ,
    #"2022_TTWJetsToLNu"            : { "histAxisName" : "2022_TTWJetsToLNu"            , "path" : "", "xsecName" : "" , } ,
    #"2022_TTWJetsToQQ"             : { "histAxisName" : "2022_TTWJetsToQQ"             , "path" : "", "xsecName" : "" , } ,
    "2022_TTZToLL_M_4to50"         : { "histAxisName" : "2022_TTZToLL_M_4to50"         , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TTLL_MLL-4to50_13p6TeV" , } ,
    "2022_TTZToLL_M_50"            : { "histAxisName" : "2022_TTZToLL_M_50"            , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TTLL_MLL-50_13p6TeV" , } ,
    "2022_TTZToQQ"                 : { "histAxisName" : "2022_TTZToQQ"                 , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TTZToQQ_13p6TeV" , } ,
    "2022_VHnobb"                  : { "histAxisName" : "2022_VHnobb"                  , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/VH_HtoNonbb_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v1_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "VH_HtoNonbb_13p6TeV" , } ,
    "2022_WJetsToLNu"              : { "histAxisName" : "2022_WJetsToLNu"              , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2022_WWTo2L2Nu"               : { "histAxisName" : "2022_WWTo2L2Nu"               , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
    "2022_WZTo3LNu"                : { "histAxisName" : "2022_WZTo3LNu"                , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "WZto3LNu_13p6TeV" , } ,
    "2022_tW_noFullHad"            : { "histAxisName" : "2022_tW_noFullHad"            , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TWminus_DR_AtLeastOneLepton_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TWminus_DR_AtLeastOneLepton_13p6TeV" , } ,
    #"2022_tZq"                     : { "histAxisName" : "2022_tZq"                     , "path" : "", "xsecName" : "" , } ,
    "2022_tbarW_noFullHad"         : { "histAxisName" : "2022_tbarW_noFullHad"         , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TbarWplus_DR_AtLeastOneLepton_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TbarWplus_DR_AtLeastOneLepton_13p6TeV" , } ,
    "2022_ttHnobb"                 : { "histAxisName" : "2022_ttHnobb"                 , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TTHtoNon2B_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v4_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TTHtoNon2B_13p6TeV" , } ,

    #"2022_tWll"   : { "histAxisName" : "2022_tWll"   , "path" : "" , "xsecName" : "" , } ,
    "2022_WWW"    : { "histAxisName" : "2022_WWW"  , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose" , "xsecName" : "WWW_13p6TeV"  } ,
    "2022_WZZ"    : { "histAxisName" : "2022_WZZ"  , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose" , "xsecName" : "WZZ_13p6TeV"  } ,
    #"2022_ZZZ"    : { "histAxisName" : "2022_ZZZ"  , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose" , "xsecName" : "ZZZ_13p6TeV"  } ,
}
central_2022EE_bkg_dict = {

    "2022EE_ZZTo4L" : {
        "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose",
        "histAxisName": "2022EE_ZZTo4l",
        "xsecName": "ZZto4L_13p6TeV",
    },

    "2022EE_ggToZZTo2e2mu"      : { "histAxisName" : "2022EE_ggToZZTo2e2mu"       , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGluToContinto2Zto2E2Mu_TuneCP5_13p6TeV_mcfm701-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGluToContinto2Zto2E2Mu_13p6TeV", } ,
    "2022EE_ggToZZTo2e2tau"     : { "histAxisName" : "2022EE_ggToZZTo2e2tau"      , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGluToContinto2Zto2E2Tau_TuneCP5_13p6TeV_mcfm701-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v4_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGluToContinto2Zto2E2Tau_13p6TeV", } ,
    "2022EE_ggToZZTo2mu2tau"    : { "histAxisName" : "2022EE_ggToZZTo2mu2tau"     , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGluToContinto2Zto2Mu2Tau_TuneCP5_13p6TeV_mcfm701-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v4_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGluToContinto2Zto2Mu2Tau_13p6TeV", } ,
    "2022EE_ggToZZTo4e"         : { "histAxisName" : "2022EE_ggToZZTo4e"          , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGlutoContinto2Zto4E_TuneCP5_13p6TeV_mcfm-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGlutoContinto2Zto4E_13p6TeV" , } ,
    "2022EE_ggToZZTo4mu"        : { "histAxisName" : "2022EE_ggToZZTo4mu"         , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGlutoContinto2Zto4Mu_TuneCP5_13p6TeV_mcfm-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGlutoContinto2Zto4Mu_13p6TeV" , } ,
    "2022EE_ggToZZTo4tau"       : { "histAxisName" : "2022EE_ggToZZTo4tau"        , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGlutoContinto2Zto4Tau_TuneCP5_13p6TeV_mcfm-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "GluGlutoContinto2Zto4Tau_13p6TeV" , } ,

    "2022EE_DYJetsToLL_M_10to50_MLM" : { "histAxisName" : "2022EE_DYJetsToLL_M_10to50_MLM" , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "DYto2L-2Jets_MLL-10to50_13p6TeV", } ,
    "2022EE_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2022EE_DYJetsToLL_M_50_MLM"     , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
    #"2022_SSWW"                    : { "histAxisName" : "2022EE_SSWW"                    , "path" : "", "xsecName" : "" , } ,
    #"2022_ST_antitop_t-channel"    : { "histAxisName" : "2022EE_ST_antitop_t-channel"    , "path" : "", "xsecName" : "" , } ,
    #"2022_ST_top_s-channel"        : { "histAxisName" : "2022EE_ST_top_s-channel"        , "path" : "", "xsecName" : "" , } ,
    #"2022_ST_top_t-channel"        : { "histAxisName" : "2022EE_ST_top_t-channel"        , "path" : "", "xsecName" : "" , } ,
    #"2022_TTTo2L2Nu"               : { "histAxisName" : "2022EE_TTTo2L2Nu"               , "path" : "", "xsecName" : "" , } ,
    #"2022_TTWJetsToLNu"            : { "histAxisName" : "2022EE_TTWJetsToLNu"            , "path" : "", "xsecName" : "" , } ,
    #"2022_TTWJetsToQQ"             : { "histAxisName" : "2022EE_TTWJetsToQQ"             , "path" : "", "xsecName" : "" , } ,
    "2022EE_TTZToLL_M_4to50"         : { "histAxisName" : "2022EE_TTZToLL_M_4to50"         , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TTLL_MLL-4to50_13p6TeV" , } ,
    "2022EE_TTZToLL_M_50"            : { "histAxisName" : "2022EE_TTZToLL_M_50"            , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TTLL_MLL-50_13p6TeV" , } ,
    "2022EE_TTZToQQ"                 : { "histAxisName" : "2022EE_TTZToQQ"                 , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TTZToQQ_13p6TeV" , } ,
    "2022EE_VHnobb"                  : { "histAxisName" : "2022EE_VHnobb"                  , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/VH_HtoNonbb_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v1_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "VH_HtoNonbb_13p6TeV" , } ,
    "2022EE_WJetsToLNu"              : { "histAxisName" : "2022EE_WJetsToLNu"              , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2022EE_WWTo2L2Nu"               : { "histAxisName" : "2022EE_WWTo2L2Nu"               , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
    "2022EE_WZTo3LNu"                : { "histAxisName" : "2022EE_WZTo3LNu"                , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "WZto3LNu_13p6TeV" , } ,
    "2022EE_tW_noFullHad"            : { "histAxisName" : "2022EE_tW_noFullHad"            , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TWminus_DR_AtLeastOneLepton_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TWminus_DR_AtLeastOneLepton_13p6TeV" , } ,
    #"2022_tZq"                     : { "histAxisName" : "2022EE_tZq"                     , "path" : "", "xsecName" : "" , } ,
    "2022EE_tbarW_noFullHad"         : { "histAxisName" : "2022EE_tbarW_noFullHad"         , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TbarWplus_DR_AtLeastOneLepton_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TbarWplus_DR_AtLeastOneLepton_13p6TeV" , } ,
    "2022EE_ttHnobb"                 : { "histAxisName" : "2022EE_ttHnobb"                 , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/TTHtoNon2B_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose", "xsecName" : "TTHtoNon2B_13p6TeV" , } ,

    #"2022_tWll"   : { "histAxisName" : "2022_tWll"   , "path" : "" , "xsecName" : "" , } ,
    "2022EE_WWW"    : { "histAxisName" : "2022EE_WWW"  , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose" , "xsecName" : "WWW_13p6TeV"  } ,
    "2022EE_WZZ"    : { "histAxisName" : "2022EE_WZZ"  , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose" , "xsecName" : "WZZ_13p6TeV"  } ,
    #"2022EE_ZZZ"    : { "histAxisName" : "2022EE_ZZZ"  , "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose" , "xsecName" : "ZZZ_13p6TeV"  } ,
}


############################ Signal samples ############################



central_UL16APV_sig_dict = {
    "UL16APV_WWZJetsTo4L2Nu" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL16APV_WWZJetsTo4L2Nu",
        "xsecName": "WWZ4l",
    },
    "UL16APV_GluGluZH" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL16APV_GluGluZH",
        "xsecName": "ggToZHToHToWWTo2L2Nu",
    },
    "UL16APV_qqToZHToZTo2L"     : { "histAxisName" : "UL16APV_qqToZHToZTo2L"     , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "qqToZHToZTo2L" , } ,
}

central_UL16_sig_dict = {
    "UL16_WWZJetsTo4L2Nu" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL16_WWZJetsTo4L2Nu",
        "xsecName": "WWZ4l",
    },
    "UL16_GluGluZH" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL16_GluGluZH",
        "xsecName": "ggToZHToHToWWTo2L2Nu",
    },
    "UL16_qqToZHToZTo2L"        : { "histAxisName" : "UL16_qqToZHToZTo2L"        , "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "qqToZHToZTo2L" , } ,
}

central_UL17_sig_dict = {
    "UL17_WWZJetsTo4L2Nu" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL17_WWZJetsTo4L2Nu",
        "xsecName": "WWZ4l",
    },
    "UL17_GluGluZH" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL17_GluGluZH",
        "xsecName": "ggToZHToHToWWTo2L2Nu",
    },
    "UL17_qqToZHToZTo2L" : { "histAxisName" : "UL17_qqToZHToZTo2L", "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "qqToZHToZTo2L" , } ,
}

central_UL18_sig_dict = {
    "UL18_WWZJetsTo4L2Nu" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL18_WWZJetsTo4L2Nu",
        "xsecName": "WWZ4l",
    },
    "UL18_GluGluZH" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL18_GluGluZH",
        "xsecName": "ggToZHToHToWWTo2L2Nu",
    },
    "UL18_qqToZHToZTo2L" : { "histAxisName" : "UL18_qqToZHToZTo2L", "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_3LepTau_4Lep" , "xsecName" : "qqToZHToZTo2L" , } ,
}

central_2022_sig_dict = {
    "2022_WWZJetsTo4L2Nu" : {
        "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WWZto4L2Nu_4F_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose",
        "histAxisName": "2022_WWZJetsTo4L2Nu",
        "xsecName": "WWZ4l_13p6TeV",
    },
    "2022_qqToZHTo2WTo2L2Nu" : {
        "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/ZH_ZtoAll_Hto2Wto2L2Nu_M-125_TuneCP5_13p6TeV_powheg-minlo-HZJ-jhugenv752-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose",
        "histAxisName": "2022_qqToZHTo2WTo2L2Nu",
        "xsecName": "ZHToHtoWWto2L2Nu_13p6TeV",
    },
    "2022_GluGluZHTo2WTo2L2Nu" : {
        "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGluZH_ZtoAll_Hto2Wto2L2Nu_M-125_TuneCP5_13p6TeV_powheg-jhugenv752-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose/",
        "histAxisName": "2022_GluGluZHTo2WTo2L2Nu",
        "xsecName": "GluGluZHToHtoWWto2L2Nu_13p6TeV",
    },
}

central_2022EE_sig_dict = {
    "2022EE_WWZJetsTo4L2Nu" : {
        "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/WWZto4L2Nu_4F_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose",
        "histAxisName": "2022EE_WWZJetsTo4L2Nu",
        "xsecName": "WWZ4l_13p6TeV",
    },
    "2022EE_qqToZHTo2WTo2L2Nu" : {
        "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/ZH_ZtoAll_Hto2Wto2L2Nu_M-125_TuneCP5_13p6TeV_powheg-minlo-HZJ-jhugenv752-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose",
        "histAxisName": "2022EE_qqToZHTo2WTo2L2Nu",
        "xsecName": "ZHToHtoWWto2L2Nu_13p6TeV",
    },
    "2022EE_GluGluZHTo2WTo2L2Nu" : {
        "path" : "/store/user/t2/users/matthew.dittrich/skims/Run3_MC_4L_04022024_Loose/GluGluZH_ZtoAll_Hto2Wto2L2Nu_M-125_TuneCP5_13p6TeV_powheg-jhugenv752-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_Run3_MC_4L_04022024_Loose/",
        "histAxisName": "2022EE_GluGluZHTo2WTo2L2Nu",
        "xsecName": "GluGluZHToHtoWWto2L2Nu_13p6TeV",
    },
}


########################### 4L Skim #######################################

# Data
data_UL16APV_4l_skim_dict = {
    "DoubleEG_Run2016B-ver1_HIPM"          : {"path" : "/skim4l_v5/DoubleEG_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016B-ver2_HIPM"          : {"path" : "/skim4l_v5/DoubleEG_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v3_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016C-HIPM"               : {"path" : "/skim4l_v5/DoubleEG_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"        , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016D-HIPM"               : {"path" : "/skim4l_v5/DoubleEG_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"        , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016E-HIPM"               : {"path" : "/skim4l_v5/DoubleEG_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"        , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016F-HIPM"               : {"path" : "/skim4l_v5/DoubleEG_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"        , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2016B-ver1_HIPM_UL2016"     : {"path" : "/skim4l_v5/MuonEG_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2016B-ver2_HIPM_UL2016"     : {"path" : "/skim4l_v5/MuonEG_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2016C-HIPM_UL2016"          : {"path" : "/skim4l_v5/MuonEG_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"          , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2016D-HIPM_UL2016"          : {"path" : "/skim4l_v5/MuonEG_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"          , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2016E-HIPM_UL2016"          : {"path" : "/skim4l_v5/MuonEG_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"          , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2016F-HIPM_UL2016"          : {"path" : "/skim4l_v5/MuonEG_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"          , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016B-ver1_HIPM_UL2016" : {"path" : "/skim4l_v5/DoubleMuon_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016B-ver2_HIPM_UL2016" : {"path" : "/skim4l_v5/DoubleMuon_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016C-HIPM_UL2016"      : {"path" : "/skim4l_v5/DoubleMuon_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"      , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016D-HIPM_UL2016"      : {"path" : "/skim4l_v5/DoubleMuon_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"      , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016E-HIPM_UL2016"      : {"path" : "/skim4l_v5/DoubleMuon_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"      , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016F-HIPM_UL2016"      : {"path" : "/skim4l_v5/DoubleMuon_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/"      , "histAxisName" : "UL16APV_data" , "xsecName" : "data" , } ,
}

data_UL16_4l_skim_dict = {
    "DoubleMuon_Run2016F-UL2016" : {"path" : "/skim4l_v5/DoubleMuon_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL16_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016G-UL2016" : {"path" : "/skim4l_v5/DoubleMuon_Run2016G-UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL16_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016H-UL2016" : {"path" : "/skim4l_v5/DoubleMuon_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL16_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016F-UL2016"   : {"path" : "/skim4l_v5/DoubleEG_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL16_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016G-UL2016"   : {"path" : "/skim4l_v5/DoubleEG_Run2016G-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL16_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016H-UL2016"   : {"path" : "/skim4l_v5/DoubleEG_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL16_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2016F-UL2016"     : {"path" : "/skim4l_v5/MuonEG_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL16_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2016G-UL2016"     : {"path" : "/skim4l_v5/MuonEG_Run2016G-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL16_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2016H-UL2016"     : {"path" : "/skim4l_v5/MuonEG_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL16_data" , "xsecName" : "data" , } ,
}

data_UL17_4l_skim_dict = {
    "DoubleMuon_Run2017B-UL2017" : {"path" : "/skim4l_v5/DoubleMuon_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017C-UL2017" : {"path" : "/skim4l_v5/DoubleMuon_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017D-UL2017" : {"path" : "/skim4l_v5/DoubleMuon_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017E-UL2017" : {"path" : "/skim4l_v5/DoubleMuon_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017F-UL2017" : {"path" : "/skim4l_v5/DoubleMuon_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017B-UL2017"   : {"path" : "/skim4l_v5/DoubleEG_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017C-UL2017"   : {"path" : "/skim4l_v5/DoubleEG_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017D-UL2017"   : {"path" : "/skim4l_v5/DoubleEG_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017E-UL2017"   : {"path" : "/skim4l_v5/DoubleEG_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017F-UL2017"   : {"path" : "/skim4l_v5/DoubleEG_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"   , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2017B-UL2017"     : {"path" : "/skim4l_v5/MuonEG_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2017C-UL2017"     : {"path" : "/skim4l_v5/MuonEG_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2017D-UL2017"     : {"path" : "/skim4l_v5/MuonEG_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2017E-UL2017"     : {"path" : "/skim4l_v5/MuonEG_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2017F-UL2017"     : {"path" : "/skim4l_v5/MuonEG_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL17_data" , "xsecName" : "data" , } ,
}

data_UL18_4l_skim_dict = {
    "MuonEG_Run2018A-UL2018"     : {"path" : "/skim4l_v5/MuonEG_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2018B-UL2018"     : {"path" : "/skim4l_v5/MuonEG_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2018C-UL2018"     : {"path" : "/skim4l_v5/MuonEG_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "MuonEG_Run2018D-UL2018"     : {"path" : "/skim4l_v5/MuonEG_Run2018D-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "EGamma_Run2018A-UL2018"     : {"path" : "/skim4l_v5/EGamma_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "EGamma_Run2018B-UL2018"     : {"path" : "/skim4l_v5/EGamma_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "EGamma_Run2018C-UL2018"     : {"path" : "/skim4l_v5/EGamma_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "EGamma_Run2018D-UL2018"     : {"path" : "/skim4l_v5/EGamma_Run2018D-UL2018_MiniAODv2_NanoAODv9-v3_NANOAOD_skim4l_v5/merged/"     , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018A-UL2018" : {"path" : "/skim4l_v5/DoubleMuon_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018B-UL2018" : {"path" : "/skim4l_v5/DoubleMuon_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018C-UL2018" : {"path" : "/skim4l_v5/DoubleMuon_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018D-UL2018" : {"path" : "/skim4l_v5/DoubleMuon_Run2018D-UL2018_MiniAODv2_NanoAODv9-v2_NANOAOD_skim4l_v5/merged/" , "histAxisName" : "UL18_data" , "xsecName" : "data" , } ,
}

## Signal
central_UL16APV_sig_4l_skim_dict = {
    "UL16APV_WWZJetsTo4L2Nu" : {"path" : "/skim4l_v5/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL16APV_WWZJetsTo4L2Nu" , "xsecName" : "WWZ4l"                , } ,
    # PSWeight sample that has more stats but isn't the one we use.
    # "UL16APV_GluGluZH"       : {"path" : "/skim4l_v5/GluGluZH_HToWWTo2L2Nu_M125_13TeV_powheg_pythia8_TuneCP5_PSweights_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"   , "histAxisName" : "UL16APV_GluGluZH"       , "xsecName" : "ggToZHToHToWWTo2L2Nu" , } ,
    "UL16APV_GluGluZH"       : {"path" : "/skim4l_v5/GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"            , "histAxisName" : "UL16APV_GluGluZH"       , "xsecName" : "ggToZHToHToWWTo2L2Nu" , } ,
    "UL16APV_qqToZHToZTo2L"  : {"path" : "/skim4l_v5/HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/" , "histAxisName" : "UL16APV_qqToZHToZTo2L"  , "xsecName" : "qqToZHToZTo2L"        , } ,
}

central_UL16_sig_4l_skim_dict = {
    "UL16_WWZJetsTo4L2Nu" : {"path" : "/skim4l_v5/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL16_WWZJetsTo4L2Nu" , "xsecName" : "WWZ4l"                , } ,
    # PSWeight sample that has more stats but isn't the one we use.
    # "UL16_GluGluZH"       : {"path" : "/skim4l_v5/GluGluZH_HToWWTo2L2Nu_M125_13TeV_powheg_pythia8_TuneCP5_PSweights_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_skim4l_v5/merged/"   , "histAxisName" : "UL16_GluGluZH"       , "xsecName" : "ggToZHToHToWWTo2L2Nu" , } ,
    "UL16_GluGluZH"       : {"path" : "/skim4l_v5/GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_skim4l_v5/merged/"            , "histAxisName" : "UL16_GluGluZH"       , "xsecName" : "ggToZHToHToWWTo2L2Nu" , } ,
    "UL16_qqToZHToZTo2L"  : {"path" : "/skim4l_v5/HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_skim4l_v5/merged/" , "histAxisName" : "UL16_qqToZHToZTo2L"  , "xsecName" : "qqToZHToZTo2L"        , } ,
}

central_UL17_sig_4l_skim_dict = {
    "UL17_WWZJetsTo4L2Nu" : {"path" : "/skim4l_v5/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL17_WWZJetsTo4L2Nu" , "xsecName" : "WWZ4l"                , } ,
    # PSWeight sample that has more stats but isn't the one we use.
    # "UL17_GluGluZH"       : {"path" : "/skim4l_v5/GluGluZH_HToWWTo2L2Nu_M125_13TeV_powheg_pythia8_TuneCP5_PSweights_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"   , "histAxisName" : "UL17_GluGluZH"       , "xsecName" : "ggToZHToHToWWTo2L2Nu" , } ,
    "UL17_GluGluZH"       : {"path" : "/skim4l_v5/GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"            , "histAxisName" : "UL17_GluGluZH"       , "xsecName" : "ggToZHToHToWWTo2L2Nu" , } ,
    "UL17_qqToZHToZTo2L"  : {"path" : "/skim4l_v5/HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/" , "histAxisName" : "UL17_qqToZHToZTo2L"  , "xsecName" : "qqToZHToZTo2L"        , } ,
}

central_UL18_sig_4l_skim_dict = {
    "UL18_WWZJetsTo4L2Nu" : {"path" : "/skim4l_v5/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL18_WWZJetsTo4L2Nu" , "xsecName" : "WWZ4l"                , } ,
    # PSWeight sample that has more stats but isn't the one we use.
    # "UL18_GluGluZH"       : {"path" : "/skim4l_v5/GluGluZH_HToWWTo2L2Nu_M125_13TeV_powheg_pythia8_TuneCP5_PSweights_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"   , "histAxisName" : "UL18_GluGluZH"       , "xsecName" : "ggToZHToHToWWTo2L2Nu" , } ,
    "UL18_GluGluZH"       : {"path" : "/skim4l_v5/GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"            , "histAxisName" : "UL18_GluGluZH"       , "xsecName" : "ggToZHToHToWWTo2L2Nu" , } ,
    "UL18_qqToZHToZTo2L"  : {"path" : "/skim4l_v5/HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/" , "histAxisName" : "UL18_qqToZHToZTo2L"  , "xsecName" : "qqToZHToZTo2L"        , } ,
}

## Background
central_UL16APV_bkg_4l_skim_dict = {
    "UL16APV_ZZTo4L"                  : {"path" : "/skim4l_v5/ZZTo4L_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                          , "histAxisName" : "UL16APV_ZZTo4l"                  , "xsecName" : "ZZTo4LK"                  , } ,
    "UL16APV_ggToZZTo2e2mu"           : {"path" : "/skim4l_v5/GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                      , "histAxisName" : "UL16APV_ggToZZTo2e2mu"           , "xsecName" : "ggToZZTo2e2muK"          , } ,
    "UL16APV_ggToZZTo2e2tau"          : {"path" : "/skim4l_v5/GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                     , "histAxisName" : "UL16APV_ggToZZTo2e2tau"          , "xsecName" : "ggToZZTo2e2tauK"         , } ,
    "UL16APV_ggToZZTo2mu2tau"         : {"path" : "/skim4l_v5/GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL16APV_ggToZZTo2mu2tau"         , "xsecName" : "ggToZZTo2mu2tauK"        , } ,
    "UL16APV_ggToZZTo4e"              : {"path" : "/skim4l_v5/GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                         , "histAxisName" : "UL16APV_ggToZZTo4e"              , "xsecName" : "ggToZZTo4eK"             , } ,
    "UL16APV_ggToZZTo4mu"             : {"path" : "/skim4l_v5/GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                        , "histAxisName" : "UL16APV_ggToZZTo4mu"             , "xsecName" : "ggToZZTo4muK"            , } ,
    "UL16APV_ggToZZTo4tau"            : {"path" : "/skim4l_v5/GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                       , "histAxisName" : "UL16APV_ggToZZTo4tau"            , "xsecName" : "ggToZZTo4tauK"           , } ,
    "UL16APV_DYJetsToLL_M_10to50_MLM" : {"path" : "/skim4l_v5/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                        , "histAxisName" : "UL16APV_DYJetsToLL_M_10to50_MLM" , "xsecName" : "DYJetsToLL_M_10to50_MLM" , } ,
    "UL16APV_DYJetsToLL_M_50_MLM"     : {"path" : "/skim4l_v5/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                            , "histAxisName" : "UL16APV_DYJetsToLL_M_50_MLM"     , "xsecName" : "DYJetsToLL_M_50_MLM"     , } ,
    "UL16APV_SSWW"                    : {"path" : "/skim4l_v5/SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                          , "histAxisName" : "UL16APV_SSWW"                    , "xsecName" : "SSWW"                    , } ,
    "UL16APV_ST_antitop_t-channel"    : {"path" : "/skim4l_v5/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/" , "histAxisName" : "UL16APV_ST_antitop_t-channel"    , "xsecName" : "ST_antitop_t-channel"    , } ,
    "UL16APV_ST_top_s-channel"        : {"path" : "/skim4l_v5/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                  , "histAxisName" : "UL16APV_ST_top_s-channel"        , "xsecName" : "ST_top_s-channel"        , } ,
    "UL16APV_ST_top_t-channel"        : {"path" : "/skim4l_v5/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"     , "histAxisName" : "UL16APV_ST_top_t-channel"        , "xsecName" : "ST_top_t-channel"        , } ,
    "UL16APV_TTTo2L2Nu"               : {"path" : "/skim4l_v5/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL16APV_TTTo2L2Nu"               , "xsecName" : "TTTo2L2Nu"               , } ,
    "UL16APV_TTWJetsToLNu"            : {"path" : "/skim4l_v5/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                      , "histAxisName" : "UL16APV_TTWJetsToLNu"            , "xsecName" : "TTWJetsToLNu"            , } ,
    "UL16APV_TTWJetsToQQ"             : {"path" : "/skim4l_v5/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                       , "histAxisName" : "UL16APV_TTWJetsToQQ"             , "xsecName" : "TTWJetsToQQ"             , } ,
    "UL16APV_TTZToLLNuNu_M_10"        : {"path" : "/skim4l_v5/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                              , "histAxisName" : "UL16APV_TTZToLLNuNu_M_10"        , "xsecName" : "TTZToLLNuNu_M_10"        , } ,
    "UL16APV_TTZToLL_M_1to10"         : {"path" : "/skim4l_v5/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                               , "histAxisName" : "UL16APV_TTZToLL_M_1to10"         , "xsecName" : "TTZToLL_M_1to10"         , } ,
    "UL16APV_TTZToQQ"                 : {"path" : "/skim4l_v5/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL16APV_TTZToQQ"                 , "xsecName" : "TTZToQQ"                 , } ,
    "UL16APV_VHnobb"                  : {"path" : "/skim4l_v5/VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL16APV_VHnobb"                  , "xsecName" : "VHnobb"                  , } ,
    "UL16APV_WJetsToLNu"              : {"path" : "/skim4l_v5/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_skim4l_v5/merged/"                                , "histAxisName" : "UL16APV_WJetsToLNu"              , "xsecName" : "WJetsToLNu"              , } ,
    "UL16APV_WWTo2L2Nu"               : {"path" : "/skim4l_v5/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL16APV_WWTo2L2Nu"               , "xsecName" : "WWTo2L2Nu"               , } ,
    "UL16APV_WZTo3LNu"                : {"path" : "/skim4l_v5/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                  , "histAxisName" : "UL16APV_WZTo3LNu"                , "xsecName" : "WZTo3LNu"                , } ,
    "UL16APV_tW_noFullHad"            : {"path" : "/skim4l_v5/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"              , "histAxisName" : "UL16APV_tW_noFullHad"            , "xsecName" : "tW_noFullHad"            , } ,
    "UL16APV_tZq"                     : {"path" : "/skim4l_v5/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                             , "histAxisName" : "UL16APV_tZq"                     , "xsecName" : "tZq"                     , } ,
    "UL16APV_tbarW_noFullHad"         : {"path" : "/skim4l_v5/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"          , "histAxisName" : "UL16APV_tbarW_noFullHad"         , "xsecName" : "tbarW_noFullHad"         , } ,
    "UL16APV_ttHnobb"                 : {"path" : "/skim4l_v5/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                , "histAxisName" : "UL16APV_ttHnobb"                 , "xsecName" : "ttHnobb"                 , } ,
    "UL16APV_tWll"                    : {"path" : "/skim4l_v5/ST_tWll_5f_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                    , "histAxisName" : "UL16APV_tWll"                    , "xsecName" : "tWll"                    , } ,
    "UL16APV_WWW"                     : {"path" : "/skim4l_v5/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                        , "histAxisName" : "UL16APV_WWW"                     , "xsecName" : "WWW"                     , } ,
    "UL16APV_WZZ"                     : {"path" : "/skim4l_v5/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                           , "histAxisName" : "UL16APV_WZZ"                     , "xsecName" : "WZZ"                     , } ,
    "UL16APV_ZZZ"                     : {"path" : "/skim4l_v5/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                                           , "histAxisName" : "UL16APV_ZZZ"                     , "xsecName" : "ZZZ"                     , } ,
    "UL16APV_tWZ4l"                   : {"path" : "/skim4l_v5/TWZToLL_tlept_Wlept_5f_DR_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_skim4l_v5/merged/"                     , "histAxisName" : "UL16APV_tWZ4l"                   , "xsecName" : "TWZToLL_tlept_Wlept"     , } ,
}

central_UL16_bkg_4l_skim_dict = {
    "UL16_ZZTo4L"                  : {"path" : "/skim4l_v5/ZZTo4L_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                          , "histAxisName" : "UL16_ZZTo4l"                  , "xsecName" : "ZZTo4LK"                  , } ,
    "UL16_ggToZZTo2e2mu"           : {"path" : "/skim4l_v5/GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                      , "histAxisName" : "UL16_ggToZZTo2e2mu"           , "xsecName" : "ggToZZTo2e2muK"          , } ,
    "UL16_ggToZZTo2e2tau"          : {"path" : "/skim4l_v5/GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                     , "histAxisName" : "UL16_ggToZZTo2e2tau"          , "xsecName" : "ggToZZTo2e2tauK"         , } ,
    "UL16_ggToZZTo2mu2tau"         : {"path" : "/skim4l_v5/GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL16_ggToZZTo2mu2tau"         , "xsecName" : "ggToZZTo2mu2tauK"        , } ,
    "UL16_ggToZZTo4e"              : {"path" : "/skim4l_v5/GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_skim4l_v5/merged/"                         , "histAxisName" : "UL16_ggToZZTo4e"              , "xsecName" : "ggToZZTo4eK"             , } ,
    "UL16_ggToZZTo4mu"             : {"path" : "/skim4l_v5/GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_skim4l_v5/merged/"                        , "histAxisName" : "UL16_ggToZZTo4mu"             , "xsecName" : "ggToZZTo4muK"            , } ,
    "UL16_ggToZZTo4tau"            : {"path" : "/skim4l_v5/GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                       , "histAxisName" : "UL16_ggToZZTo4tau"            , "xsecName" : "ggToZZTo4tauK"           , } ,
    "UL16_DYJetsToLL_M_10to50_MLM" : {"path" : "/skim4l_v5/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                        , "histAxisName" : "UL16_DYJetsToLL_M_10to50_MLM" , "xsecName" : "DYJetsToLL_M_10to50_MLM" , } ,
    "UL16_DYJetsToLL_M_50_MLM"     : {"path" : "/skim4l_v5/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                            , "histAxisName" : "UL16_DYJetsToLL_M_50_MLM"     , "xsecName" : "DYJetsToLL_M_50_MLM"     , } ,
    "UL16_SSWW"                    : {"path" : "/skim4l_v5/SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                          , "histAxisName" : "UL16_SSWW"                    , "xsecName" : "SSWW"                    , } ,
    "UL16_ST_antitop_t-channel"    : {"path" : "/skim4l_v5/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/" , "histAxisName" : "UL16_ST_antitop_t-channel"    , "xsecName" : "ST_antitop_t-channel"    , } ,
    "UL16_ST_top_s-channel"        : {"path" : "/skim4l_v5/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                  , "histAxisName" : "UL16_ST_top_s-channel"        , "xsecName" : "ST_top_s-channel"        , } ,
    "UL16_ST_top_t-channel"        : {"path" : "/skim4l_v5/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"     , "histAxisName" : "UL16_ST_top_t-channel"        , "xsecName" : "ST_top_t-channel"        , } ,
    "UL16_TTTo2L2Nu"               : {"path" : "/skim4l_v5/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL16_TTTo2L2Nu"               , "xsecName" : "TTTo2L2Nu"               , } ,
    "UL16_TTWJetsToLNu"            : {"path" : "/skim4l_v5/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                      , "histAxisName" : "UL16_TTWJetsToLNu"            , "xsecName" : "TTWJetsToLNu"            , } ,
    "UL16_TTWJetsToQQ"             : {"path" : "/skim4l_v5/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                       , "histAxisName" : "UL16_TTWJetsToQQ"             , "xsecName" : "TTWJetsToQQ"             , } ,
    "UL16_TTZToLLNuNu_M_10"        : {"path" : "/skim4l_v5/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                              , "histAxisName" : "UL16_TTZToLLNuNu_M_10"        , "xsecName" : "TTZToLLNuNu_M_10"        , } ,
    "UL16_TTZToLL_M_1to10"         : {"path" : "/skim4l_v5/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                               , "histAxisName" : "UL16_TTZToLL_M_1to10"         , "xsecName" : "TTZToLL_M_1to10"         , } ,
    "UL16_TTZToQQ"                 : {"path" : "/skim4l_v5/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL16_TTZToQQ"                 , "xsecName" : "TTZToQQ"                 , } ,
    "UL16_VHnobb"                  : {"path" : "/skim4l_v5/VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL16_VHnobb"                  , "xsecName" : "VHnobb"                  , } ,
    "UL16_WJetsToLNu"              : {"path" : "/skim4l_v5/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_skim4l_v5/merged/"                                , "histAxisName" : "UL16_WJetsToLNu"              , "xsecName" : "WJetsToLNu"              , } ,
    "UL16_WWTo2L2Nu"               : {"path" : "/skim4l_v5/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL16_WWTo2L2Nu"               , "xsecName" : "WWTo2L2Nu"               , } ,
    "UL16_WZTo3LNu"                : {"path" : "/skim4l_v5/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                  , "histAxisName" : "UL16_WZTo3LNu"                , "xsecName" : "WZTo3LNu"                , } ,
    "UL16_tW_noFullHad"            : {"path" : "/skim4l_v5/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"              , "histAxisName" : "UL16_tW_noFullHad"            , "xsecName" : "tW_noFullHad"            , } ,
    "UL16_tZq"                     : {"path" : "/skim4l_v5/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                             , "histAxisName" : "UL16_tZq"                     , "xsecName" : "tZq"                     , } ,
    "UL16_tbarW_noFullHad"         : {"path" : "/skim4l_v5/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"          , "histAxisName" : "UL16_tbarW_noFullHad"         , "xsecName" : "tbarW_noFullHad"         , } ,
    "UL16_ttHnobb"                 : {"path" : "/skim4l_v5/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                , "histAxisName" : "UL16_ttHnobb"                 , "xsecName" : "ttHnobb"                 , } ,
    "UL16_tWll"                    : {"path" : "/skim4l_v5/ST_tWll_5f_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                    , "histAxisName" : "UL16_tWll"                    , "xsecName" : "tWll"                    , } ,
    "UL16_WWW"                     : {"path" : "/skim4l_v5/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                        , "histAxisName" : "UL16_WWW"                     , "xsecName" : "WWW"                     , } ,
    "UL16_WZZ"                     : {"path" : "/skim4l_v5/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                           , "histAxisName" : "UL16_WZZ"                     , "xsecName" : "WZZ"                     , } ,
    "UL16_ZZZ"                     : {"path" : "/skim4l_v5/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                                           , "histAxisName" : "UL16_ZZZ"                     , "xsecName" : "ZZZ"                     , } ,
    "UL16_tWZ4l"                   : {"path" : "/skim4l_v5/TWZToLL_tlept_Wlept_5f_DR_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_skim4l_v5/merged/"                     , "histAxisName" : "UL16_tWZ4l"                   , "xsecName" : "TWZToLL_tlept_Wlept"     , } ,
}

central_UL17_bkg_4l_skim_dict = {
    "UL17_ZZTo4L"                  : {"path" : "/skim4l_v5/ZZTo4L_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                                          , "histAxisName" : "UL17_ZZTo4l"                  , "xsecName" : "ZZTo4LK"                  , } ,
    "UL17_ggToZZTo2e2mu"           : {"path" : "/skim4l_v5/GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                      , "histAxisName" : "UL17_ggToZZTo2e2mu"           , "xsecName" : "ggToZZTo2e2muK"          , } ,
    "UL17_ggToZZTo2e2tau"          : {"path" : "/skim4l_v5/GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                     , "histAxisName" : "UL17_ggToZZTo2e2tau"          , "xsecName" : "ggToZZTo2e2tauK"         , } ,
    "UL17_ggToZZTo2mu2tau"         : {"path" : "/skim4l_v5/GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL17_ggToZZTo2mu2tau"         , "xsecName" : "ggToZZTo2mu2tauK"        , } ,
    "UL17_ggToZZTo4e"              : {"path" : "/skim4l_v5/GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                         , "histAxisName" : "UL17_ggToZZTo4e"              , "xsecName" : "ggToZZTo4eK"             , } ,
    "UL17_ggToZZTo4mu"             : {"path" : "/skim4l_v5/GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                        , "histAxisName" : "UL17_ggToZZTo4mu"             , "xsecName" : "ggToZZTo4muK"            , } ,
    "UL17_ggToZZTo4tau"            : {"path" : "/skim4l_v5/GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                       , "histAxisName" : "UL17_ggToZZTo4tau"            , "xsecName" : "ggToZZTo4tauK"           , } ,
    "UL17_DYJetsToLL_M_10to50_MLM" : {"path" : "/skim4l_v5/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                        , "histAxisName" : "UL17_DYJetsToLL_M_10to50_MLM" , "xsecName" : "DYJetsToLL_M_10to50_MLM" , } ,
    "UL17_DYJetsToLL_M_50_MLM"     : {"path" : "/skim4l_v5/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                            , "histAxisName" : "UL17_DYJetsToLL_M_50_MLM"     , "xsecName" : "DYJetsToLL_M_50_MLM"     , } ,
    "UL17_SSWW"                    : {"path" : "/skim4l_v5/SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                                          , "histAxisName" : "UL17_SSWW"                    , "xsecName" : "SSWW"                    , } ,
    "UL17_ST_antitop_t-channel"    : {"path" : "/skim4l_v5/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/" , "histAxisName" : "UL17_ST_antitop_t-channel"    , "xsecName" : "ST_antitop_t-channel"    , } ,
    "UL17_ST_top_s-channel"        : {"path" : "/skim4l_v5/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                  , "histAxisName" : "UL17_ST_top_s-channel"        , "xsecName" : "ST_top_s-channel"        , } ,
    "UL17_ST_top_t-channel"        : {"path" : "/skim4l_v5/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"     , "histAxisName" : "UL17_ST_top_t-channel"        , "xsecName" : "ST_top_t-channel"        , } ,
    "UL17_TTTo2L2Nu"               : {"path" : "/skim4l_v5/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL17_TTTo2L2Nu"               , "xsecName" : "TTTo2L2Nu"               , } ,
    "UL17_TTWJetsToLNu"            : {"path" : "/skim4l_v5/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                      , "histAxisName" : "UL17_TTWJetsToLNu"            , "xsecName" : "TTWJetsToLNu"            , } ,
    "UL17_TTWJetsToQQ"             : {"path" : "/skim4l_v5/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                       , "histAxisName" : "UL17_TTWJetsToQQ"             , "xsecName" : "TTWJetsToQQ"             , } ,
    "UL17_TTZToLLNuNu_M_10"        : {"path" : "/skim4l_v5/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                              , "histAxisName" : "UL17_TTZToLLNuNu_M_10"        , "xsecName" : "TTZToLLNuNu_M_10"        , } ,
    "UL17_TTZToLL_M_1to10"         : {"path" : "/skim4l_v5/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                               , "histAxisName" : "UL17_TTZToLL_M_1to10"         , "xsecName" : "TTZToLL_M_1to10"         , } ,
    "UL17_TTZToQQ"                 : {"path" : "/skim4l_v5/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL17_TTZToQQ"                 , "xsecName" : "TTZToQQ"                 , } ,
    "UL17_VHnobb"                  : {"path" : "/skim4l_v5/VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL17_VHnobb"                  , "xsecName" : "VHnobb"                  , } ,
    "UL17_WJetsToLNu"              : {"path" : "/skim4l_v5/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                                , "histAxisName" : "UL17_WJetsToLNu"              , "xsecName" : "WJetsToLNu"              , } ,
    "UL17_WWTo2L2Nu"               : {"path" : "/skim4l_v5/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL17_WWTo2L2Nu"               , "xsecName" : "WWTo2L2Nu"               , } ,
    "UL17_WZTo3LNu"                : {"path" : "/skim4l_v5/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_skim4l_v5/merged/"                                  , "histAxisName" : "UL17_WZTo3LNu"                , "xsecName" : "WZTo3LNu"                , } ,
    "UL17_tW_noFullHad"            : {"path" : "/skim4l_v5/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"              , "histAxisName" : "UL17_tW_noFullHad"            , "xsecName" : "tW_noFullHad"            , } ,
    "UL17_tZq"                     : {"path" : "/skim4l_v5/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                             , "histAxisName" : "UL17_tZq"                     , "xsecName" : "tZq"                     , } ,
    "UL17_tbarW_noFullHad"         : {"path" : "/skim4l_v5/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"          , "histAxisName" : "UL17_tbarW_noFullHad"         , "xsecName" : "tbarW_noFullHad"         , } ,
    "UL17_ttHnobb"                 : {"path" : "/skim4l_v5/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                , "histAxisName" : "UL17_ttHnobb"                 , "xsecName" : "ttHnobb"                 , } ,
    "UL17_tWll"                    : {"path" : "/skim4l_v5/ST_tWll_5f_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                                    , "histAxisName" : "UL17_tWll"                    , "xsecName" : "tWll"                    , } ,
    "UL17_WWW"                     : {"path" : "/skim4l_v5/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                                        , "histAxisName" : "UL17_WWW"                     , "xsecName" : "WWW"                     , } ,
    "UL17_WZZ"                     : {"path" : "/skim4l_v5/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                                           , "histAxisName" : "UL17_WZZ"                     , "xsecName" : "WZZ"                     , } ,
    "UL17_ZZZ"                     : {"path" : "/skim4l_v5/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                                           , "histAxisName" : "UL17_ZZZ"                     , "xsecName" : "ZZZ"                     , } ,
    "UL17_tWZ4l"                   : {"path" : "/skim4l_v5/TWZToLL_tlept_Wlept_5f_DR_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_skim4l_v5/merged/"                     , "histAxisName" : "UL17_tWZ4l"                   , "xsecName" : "TWZToLL_tlept_Wlept"     , } ,
}

central_UL18_bkg_4l_skim_dict = {
    "UL18_ZZTo4L"                  : {"path" : "/skim4l_v5/ZZTo4L_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                                          , "histAxisName" : "UL18_ZZTo4l"                  , "xsecName" : "ZZTo4LK"                  , } ,
    "UL18_ggToZZTo2e2mu"           : {"path" : "/skim4l_v5/GluGluToContinToZZTo2e2mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                      , "histAxisName" : "UL18_ggToZZTo2e2mu"           , "xsecName" : "ggToZZTo2e2muK"          , } ,
    "UL18_ggToZZTo2e2tau"          : {"path" : "/skim4l_v5/GluGluToContinToZZTo2e2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                     , "histAxisName" : "UL18_ggToZZTo2e2tau"          , "xsecName" : "ggToZZTo2e2tauK"         , } ,
    "UL18_ggToZZTo2mu2tau"         : {"path" : "/skim4l_v5/GluGluToContinToZZTo2mu2tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL18_ggToZZTo2mu2tau"         , "xsecName" : "ggToZZTo2mu2tauK"        , } ,
    "UL18_ggToZZTo4e"              : {"path" : "/skim4l_v5/GluGluToContinToZZTo4e_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                         , "histAxisName" : "UL18_ggToZZTo4e"              , "xsecName" : "ggToZZTo4eK"             , } ,
    "UL18_ggToZZTo4mu"             : {"path" : "/skim4l_v5/GluGluToContinToZZTo4mu_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                        , "histAxisName" : "UL18_ggToZZTo4mu"             , "xsecName" : "ggToZZTo4muK"            , } ,
    "UL18_ggToZZTo4tau"            : {"path" : "/skim4l_v5/GluGluToContinToZZTo4tau_TuneCP5_13TeV-mcfm701-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                       , "histAxisName" : "UL18_ggToZZTo4tau"            , "xsecName" : "ggToZZTo4tauK"           , } ,
    "UL18_DYJetsToLL_M_10to50_MLM" : {"path" : "/skim4l_v5/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                        , "histAxisName" : "UL18_DYJetsToLL_M_10to50_MLM" , "xsecName" : "DYJetsToLL_M_10to50_MLM" , } ,
    "UL18_DYJetsToLL_M_50_MLM"     : {"path" : "/skim4l_v5/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                            , "histAxisName" : "UL18_DYJetsToLL_M_50_MLM"     , "xsecName" : "DYJetsToLL_M_50_MLM"     , } ,
    "UL18_SSWW"                    : {"path" : "/skim4l_v5/SSWW_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                                          , "histAxisName" : "UL18_SSWW"                    , "xsecName" : "SSWW"                    , } ,
    "UL18_ST_antitop_t-channel"    : {"path" : "/skim4l_v5/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/" , "histAxisName" : "UL18_ST_antitop_t-channel"    , "xsecName" : "ST_antitop_t-channel"    , } ,
    "UL18_ST_top_s-channel"        : {"path" : "/skim4l_v5/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                  , "histAxisName" : "UL18_ST_top_s-channel"        , "xsecName" : "ST_top_s-channel"        , } ,
    "UL18_ST_top_t-channel"        : {"path" : "/skim4l_v5/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"     , "histAxisName" : "UL18_ST_top_t-channel"        , "xsecName" : "ST_top_t-channel"        , } ,
    "UL18_TTTo2L2Nu"               : {"path" : "/skim4l_v5/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL18_TTTo2L2Nu"               , "xsecName" : "TTTo2L2Nu"               , } ,
    "UL18_TTWJetsToLNu"            : {"path" : "/skim4l_v5/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                      , "histAxisName" : "UL18_TTWJetsToLNu"            , "xsecName" : "TTWJetsToLNu"            , } ,
    "UL18_TTWJetsToQQ"             : {"path" : "/skim4l_v5/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                       , "histAxisName" : "UL18_TTWJetsToQQ"             , "xsecName" : "TTWJetsToQQ"             , } ,
    "UL18_TTZToLLNuNu_M_10"        : {"path" : "/skim4l_v5/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                              , "histAxisName" : "UL18_TTZToLLNuNu_M_10"        , "xsecName" : "TTZToLLNuNu_M_10"        , } ,
    "UL18_TTZToLL_M_1to10"         : {"path" : "/skim4l_v5/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                               , "histAxisName" : "UL18_TTZToLL_M_1to10"         , "xsecName" : "TTZToLL_M_1to10"         , } ,
    "UL18_TTZToQQ"                 : {"path" : "/skim4l_v5/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL18_TTZToQQ"                 , "xsecName" : "TTZToQQ"                 , } ,
    "UL18_VHnobb"                  : {"path" : "/skim4l_v5/VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                    , "histAxisName" : "UL18_VHnobb"                  , "xsecName" : "VHnobb"                  , } ,
    "UL18_WJetsToLNu"              : {"path" : "/skim4l_v5/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                                , "histAxisName" : "UL18_WJetsToLNu"              , "xsecName" : "WJetsToLNu"              , } ,
    "UL18_WWTo2L2Nu"               : {"path" : "/skim4l_v5/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                                       , "histAxisName" : "UL18_WWTo2L2Nu"               , "xsecName" : "WWTo2L2Nu"               , } ,
    "UL18_WZTo3LNu"                : {"path" : "/skim4l_v5/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_skim4l_v5/merged/"                                  , "histAxisName" : "UL18_WZTo3LNu"                , "xsecName" : "WZTo3LNu"                , } ,
    "UL18_tW_noFullHad"            : {"path" : "/skim4l_v5/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"              , "histAxisName" : "UL18_tW_noFullHad"            , "xsecName" : "tW_noFullHad"            , } ,
    "UL18_tZq"                     : {"path" : "/skim4l_v5/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                             , "histAxisName" : "UL18_tZq"                     , "xsecName" : "tZq"                     , } ,
    "UL18_tbarW_noFullHad"         : {"path" : "/skim4l_v5/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"          , "histAxisName" : "UL18_tbarW_noFullHad"         , "xsecName" : "tbarW_noFullHad"         , } ,
    "UL18_ttHnobb"                 : {"path" : "/skim4l_v5/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                , "histAxisName" : "UL18_ttHnobb"                 , "xsecName" : "ttHnobb"                 , } ,
    "UL18_tWll"                    : {"path" : "/skim4l_v5/ST_tWll_5f_TuneCP5_13TeV-madgraph-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                                    , "histAxisName" : "UL18_tWll"                    , "xsecName" : "tWll"                    , } ,
    "UL18_WWW"                     : {"path" : "/skim4l_v5/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                                        , "histAxisName" : "UL18_WWW"                     , "xsecName" : "WWW"                     , } ,
    "UL18_WZZ"                     : {"path" : "/skim4l_v5/WZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                                           , "histAxisName" : "UL18_WZZ"                     , "xsecName" : "WZZ"                     , } ,
    "UL18_ZZZ"                     : {"path" : "/skim4l_v5/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                                           , "histAxisName" : "UL18_ZZZ"                     , "xsecName" : "ZZZ"                     , } ,
    "UL18_tWZ4l"                   : {"path" : "/skim4l_v5/TWZToLL_tlept_Wlept_5f_DR_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_skim4l_v5/merged/"                     , "histAxisName" : "UL18_tWZ4l"                   , "xsecName" : "TWZToLL_tlept_Wlept"     , } ,
}



############################ Test example samples ############################

# Test dict
test_wwz_dict = {
    "UL17_WWZJetsTo4L2Nu" : {
        "path" : "/store/user/kmohrman/samples/from_keegan_skims_3LepTau_4Lep/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL17_WWZJetsTo4L2Nu",
        "xsecName": "WWZ4l",
    },
}

# CI example
# Note, if using this to remake the json for the CI, should replace the file in the "files" with just output_1.root (no path) since assumes will be downloaded locally
ci_dict = {
    "UL17_WWZJetsTo4L2Nu_forCI" : {
        "path" : "/cmsuf/data/store/user/kmohrman/test/ci_for_wwz",
        "histAxisName": "UL17_WWZJetsTo4L2Nu",
        "xsecName": "WWZ4l",
    },
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

    # A simple example
    #make_jsons_for_dict_of_samples(test_wwz_dict, "/ceph/cms/","2017",".",on_das=False) # An example
    #make_jsons_for_dict_of_samples(ci_dict, "","2017","../../input_samples/sample_jsons/test_samples/",on_das=False) # For CI json

    # Specify output paths
    jsons_path = "../../input_samples/sample_jsons/"
    #out_dir_data_16APV = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/UL16APV")
    #out_dir_data_16 = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/UL16")
    #out_dir_data_17 = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/UL17")
    #out_dir_data_18 = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/UL18")
    out_dir_data_22 = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/2022")
    out_dir_data_22EE = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/2022EE")
    out_dir_bkg = os.path.join(jsons_path,"wwz_analysis_skims_v0/bkg_samples/")
    out_dir_sig = os.path.join(jsons_path,"wwz_analysis_skims_v0/sig_samples/")

    # Make configs for data samples
    #make_jsons_for_dict_of_samples(data_UL16APV, "/ceph/cms/","2016APV", out_dir_data_16APV,on_das=False)
    #make_jsons_for_dict_of_samples(data_UL16, "/ceph/cms/","2016", out_dir_data_16,on_das=False)
    #make_jsons_for_dict_of_samples(data_UL17, "/ceph/cms/","2017", out_dir_data_17,on_das=False)
    #make_jsons_for_dict_of_samples(data_UL18, "/ceph/cms/","2018", out_dir_data_18,on_das=False)
    make_jsons_for_dict_of_samples(data_2022, "/cmsuf/data/","2022", out_dir_data_22,era_op=1,on_das=False)
    make_jsons_for_dict_of_samples(data_2022EE, "/cmsuf/data/","2022EE", out_dir_data_22EE,era_op=1,on_das=False)

    # Make configs for bkg samples
    #make_jsons_for_dict_of_samples(central_UL16APV_22006_dict, "/cmsuf/data/","2016APV", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL16_22006_dict, "/cmsuf/data/","2016", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL17_22006_dict, "/cmsuf/data/","2017", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL18_22006_dict, "/cmsuf/data/","2018", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2022_dict, "/cmsuf/data/","2022", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2022EE_dict, "/cmsuf/data/","2022EE", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL16APV_bkg_dict, "/ceph/cms/","2016APV", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL16_bkg_dict, "/ceph/cms/","2016", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL17_bkg_dict, "/ceph/cms/","2017", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL18_bkg_dict, "/ceph/cms/","2018", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2022_bkg_dict, "/cmsuf/data/","2022", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2022EE_bkg_dict, "/cmsuf/data/","2022EE", out_dir_bkg,on_das=False)

    # Make configs for sig samples
    #make_jsons_for_dict_of_samples(central_UL16APV_sig_dict, "/ceph/cms/","2016APV", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL16_sig_dict, "/ceph/cms/","2016", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL17_sig_dict, "/ceph/cms/","2017", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL18_sig_dict, "/ceph/cms/","2018", out_dir_sig,on_das=False)
    make_jsons_for_dict_of_samples(central_2022_sig_dict, "/cmsuf/data/","2022", out_dir_sig,on_das=False)
    make_jsons_for_dict_of_samples(central_2022EE_sig_dict, "/cmsuf/data/","2022EE", out_dir_sig,on_das=False)

    # Replace xsec numbers
    #replace_xsec_for_dict_of_samples(central_UL16APV_bkg_dict,out_dir_bkg)
    #replace_xsec_for_dict_of_samples(central_UL16_bkg_dict,out_dir_bkg)
    #replace_xsec_for_dict_of_samples(central_UL17_bkg_dict,out_dir_bkg)
    #replace_xsec_for_dict_of_samples(central_UL18_bkg_dict,out_dir_bkg)

    # 4L Skim
    #jsons_path = "../../input_samples/sample_jsons/"
    #out_dir_data = os.path.join(jsons_path,"wwz_analysis_4l_skims_v0/data_samples")
    #os.makedirs(out_dir_data, exist_ok=True)
    #out_dir_bkg = os.path.join(jsons_path,"wwz_analysis_4l_skims_v0/bkg_samples")
    #os.makedirs(out_dir_bkg, exist_ok=True)
    #out_dir_sig = os.path.join(jsons_path,"wwz_analysis_4l_skims_v0/sig_samples")
    #os.makedirs(out_dir_sig, exist_ok=True)

    #make_jsons_for_dict_of_samples(central_UL16APV_sig_4l_skim_dict, "/data/userdata/phchang/", "2016APV", out_dir_sig, on_das=False)
    #make_jsons_for_dict_of_samples(central_UL16_sig_4l_skim_dict, "/data/userdata/phchang/", "2016", out_dir_sig, on_das=False)
    #make_jsons_for_dict_of_samples(central_UL17_sig_4l_skim_dict, "/data/userdata/phchang/", "2017", out_dir_sig, on_das=False)
    #make_jsons_for_dict_of_samples(central_UL18_sig_4l_skim_dict, "/data/userdata/phchang/", "2018", out_dir_sig, on_das=False)

    #make_jsons_for_dict_of_samples(central_UL16APV_bkg_4l_skim_dict, "/data/userdata/phchang/", "2016APV", out_dir_bkg, on_das=False)
    #make_jsons_for_dict_of_samples(central_UL16_bkg_4l_skim_dict, "/data/userdata/phchang/", "2016", out_dir_bkg, on_das=False)
    #make_jsons_for_dict_of_samples(central_UL17_bkg_4l_skim_dict, "/data/userdata/phchang/", "2017", out_dir_bkg, on_das=False)
    #make_jsons_for_dict_of_samples(central_UL18_bkg_4l_skim_dict, "/data/userdata/phchang/", "2018", out_dir_bkg, on_das=False)

    #make_jsons_for_dict_of_samples(data_UL16APV_4l_skim_dict , "/data/userdata/phchang/", "2016APV", out_dir_data, on_das=False)
    #make_jsons_for_dict_of_samples(data_UL16_4l_skim_dict , "/data/userdata/phchang/", "2016", out_dir_data, on_das=False)
    #make_jsons_for_dict_of_samples(data_UL17_4l_skim_dict , "/data/userdata/phchang/", "2017", out_dir_data, on_das=False)
    #make_jsons_for_dict_of_samples(data_UL18_4l_skim_dict , "/data/userdata/phchang/", "2018", out_dir_data, on_das=False)


if __name__ == "__main__":
    main()
