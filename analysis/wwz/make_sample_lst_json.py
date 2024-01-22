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



############################ Bkg samples ############################


central_UL16APV_bkg_dict = {

    "UL16APV_ZZTo4L" : {
        "path" : "/store/user/kdownham/skimOutput/3LepTau_4Lep/ZZTo4L_TuneCP5_13TeV_powheg_pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_3LepTau_4Lep",
        "histAxisName": "UL16APV_ZZTo4l",
        "xsecName": "ZZTo4L",
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
        "xsecName": "ZZTo4L",
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
        "xsecName": "ZZTo4L",
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
        "xsecName": "ZZTo4L",
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
def make_jsons_for_dict_of_samples(samples_dict,prefix,year,out_dir,on_das=False):
    failed = []
    for sample_name,sample_info in sorted(samples_dict.items()):
        print(f"\n\nMaking JSON for {sample_name}...")
        path = sample_info["path"]
        hist_axis_name = sample_info["histAxisName"]
        xsec_name = sample_info["xsecName"]
        print("we are here")
        sjt.make_json(
            sample_dir = path,
            sample_name = sample_name,
            prefix = prefix,
            sample_yr = year,
            xsec_name = xsec_name,
            hist_axis_name = hist_axis_name,
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
    make_jsons_for_dict_of_samples(ci_dict, "","2017","../../input_samples/sample_jsons/test_samples/",on_das=False) # For CI json

    # Specify output paths
    jsons_path = "../../input_samples/sample_jsons/"
    out_dir_data_16APV = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/UL16APV")
    out_dir_data_16 = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/UL16")
    out_dir_data_17 = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/UL17")
    out_dir_data_18 = os.path.join(jsons_path,"wwz_analysis_skims_v0/data_samples/UL18")
    out_dir_bkg = os.path.join(jsons_path,"wwz_analysis_skims_v0/bkg_samples/")
    out_dir_sig = os.path.join(jsons_path,"wwz_analysis_skims_v0/sig_samples/")

    # Make configs for data samples
    #make_jsons_for_dict_of_samples(data_UL16APV, "/ceph/cms/","2016APV", out_dir_data_16APV,on_das=False)
    #make_jsons_for_dict_of_samples(data_UL16, "/ceph/cms/","2016", out_dir_data_16,on_das=False)
    #make_jsons_for_dict_of_samples(data_UL17, "/ceph/cms/","2017", out_dir_data_17,on_das=False)
    #make_jsons_for_dict_of_samples(data_UL18, "/ceph/cms/","2018", out_dir_data_18,on_das=False)

    # Make configs for bkg samples
    make_jsons_for_dict_of_samples(central_UL16APV_bkg_dict, "/ceph/cms/","2016APV", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_UL16_bkg_dict, "/ceph/cms/","2016", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_UL17_bkg_dict, "/ceph/cms/","2017", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_UL18_bkg_dict, "/ceph/cms/","2018", out_dir_bkg,on_das=False)

    # Make configs for sig samples
    #make_jsons_for_dict_of_samples(central_UL16APV_sig_dict, "/ceph/cms/","2016APV", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL16_sig_dict, "/ceph/cms/","2016", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL17_sig_dict, "/ceph/cms/","2017", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_UL18_sig_dict, "/ceph/cms/","2018", out_dir_sig,on_das=False)

    # Replace xsec numbers
    #replace_xsec_for_dict_of_samples(central_UL16APV_bkg_dict,out_dir_bkg)
    #replace_xsec_for_dict_of_samples(central_UL16_bkg_dict,out_dir_bkg)
    #replace_xsec_for_dict_of_samples(central_UL17_bkg_dict,out_dir_bkg)
    #replace_xsec_for_dict_of_samples(central_UL18_bkg_dict,out_dir_bkg)


if __name__ == "__main__":
    main()
