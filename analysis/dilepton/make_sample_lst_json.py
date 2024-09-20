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
    "DoubleEG_Run2016B-ver1_HIPM"          : { "Era" : "B" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleEG_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016B-ver2_HIPM"          : { "Era" : "B" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleEG_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v3_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016C-HIPM"               : { "Era" : "C" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleEG_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016D-HIPM"               : { "Era" : "D" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleEG_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016E-HIPM"               : { "Era" : "E" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleEG_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016F-HIPM"               : { "Era" : "F" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleEG_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,

    "MuonEG_Run2016B-ver1_HIPM_UL2016"     : { "Era" : "B" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/MuonEG_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2016B-ver2_HIPM_UL2016"     : { "Era" : "B" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/MuonEG_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2016C-HIPM_UL2016"          : { "Era" : "C" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/MuonEG_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2016D-HIPM_UL2016"          : { "Era" : "D" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/MuonEG_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2016E-HIPM_UL2016"          : { "Era" : "E" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/MuonEG_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2016F-HIPM_UL2016"          : { "Era" : "F" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/MuonEG_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,

    "DoubleMuon_Run2016B-ver1_HIPM_UL2016" : { "Era" : "B" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleMuon_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016B-ver2_HIPM_UL2016" : { "Era" : "B" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleMuon_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016C-HIPM_UL2016"      : { "Era" : "C" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleMuon_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016D-HIPM_UL2016"      : { "Era" : "D" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleMuon_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016E-HIPM_UL2016"      : { "Era" : "E" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleMuon_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016F-HIPM_UL2016"      : { "Era" : "F" , "histAxisName" : "UL16APV_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016APV/DoubleMuon_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
}

data_UL16 = {
    "DoubleMuon_Run2016F-UL2016" : { "Era" : "F" , "histAxisName" : "UL16_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016/DoubleMuon_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016G-UL2016" : { "Era" : "G" , "histAxisName" : "UL16_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016/DoubleMuon_Run2016G-UL2016_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2016H-UL2016" : { "Era" : "H" , "histAxisName" : "UL16_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016/DoubleMuon_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,

    "DoubleEG_Run2016F-UL2016"   : { "Era" : "F" , "histAxisName" : "UL16_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016/DoubleEG_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016G-UL2016"   : { "Era" : "G" , "histAxisName" : "UL16_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016/DoubleEG_Run2016G-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2016H-UL2016"   : { "Era" : "H" , "histAxisName" : "UL16_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016/DoubleEG_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,

    "MuonEG_Run2016F-UL2016"     : { "Era" : "F" , "histAxisName" : "UL16_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016/MuonEG_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2016G-UL2016"     : { "Era" : "G" , "histAxisName" : "UL16_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016/MuonEG_Run2016G-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2016H-UL2016"     : { "Era" : "H" , "histAxisName" : "UL16_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2016/MuonEG_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
}

data_UL17 = {
    "DoubleMuon_Run2017B-UL2017" : { "Era" : "B" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleMuon_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017C-UL2017" : { "Era" : "C" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleMuon_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017D-UL2017" : { "Era" : "D" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleMuon_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017E-UL2017" : { "Era" : "E" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleMuon_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2017F-UL2017" : { "Era" : "F" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleMuon_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,

    "DoubleEG_Run2017B-UL2017"   : { "Era" : "B" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleEG_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017C-UL2017"   : { "Era" : "C" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleEG_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017D-UL2017"   : { "Era" : "D" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleEG_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017E-UL2017"   : { "Era" : "E" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleEG_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleEG_Run2017F-UL2017"   : { "Era" : "F" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/DoubleEG_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,

    "MuonEG_Run2017B-UL2017"     : { "Era" : "B" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/MuonEG_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2017C-UL2017"     : { "Era" : "C" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/MuonEG_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2017D-UL2017"     : { "Era" : "D" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/MuonEG_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2017E-UL2017"     : { "Era" : "E" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/MuonEG_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2017F-UL2017"     : { "Era" : "F" , "histAxisName" : "UL17_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2017/MuonEG_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
}

data_UL18 = {
    "MuonEG_Run2018A-UL2018"     : { "Era" : "A" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/MuonEG_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2018B-UL2018"     : { "Era" : "B" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/MuonEG_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2018C-UL2018"     : { "Era" : "C" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/MuonEG_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2018D-UL2018"     : { "Era" : "D" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/MuonEG_Run2018D-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,

    "EGamma_Run2018A-UL2018"     : { "Era" : "A" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/EGamma_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "EGamma_Run2018B-UL2018"     : { "Era" : "B" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/EGamma_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "EGamma_Run2018C-UL2018"     : { "Era" : "C" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/EGamma_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "EGamma_Run2018D-UL2018"     : { "Era" : "D" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/EGamma_Run2018D-UL2018_MiniAODv2_NanoAODv9-v3_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,

    "DoubleMuon_Run2018A-UL2018" : { "Era" : "A" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/DoubleMuon_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018B-UL2018" : { "Era" : "B" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/DoubleMuon_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018C-UL2018" : { "Era" : "C" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/DoubleMuon_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
    "DoubleMuon_Run2018D-UL2018" : { "Era" : "D" , "histAxisName" : "UL18_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/data/2018/DoubleMuon_Run2018D-UL2018_MiniAODv2_NanoAODv9-v2_NANOAOD_WWZ_Dilepton_Run2_09112024/" , "xsecName" : "data" , } ,
}

data_2022 = {
    "DoubleMuon_Run2022C-22Sep2023" : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022/DoubleMuon_Run2022C-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,

    "Muon_Run2022C-22Sep2023"       : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022/Muon_Run2022C-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_Run2022D-22Sep2023"       : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022/Muon_Run2022D-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,

    "EGamma_Run2022C-22Sep2023"     : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022/EGamma_Run2022C-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_Run2022D-22Sep2023"     : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022/EGamma_Run2022D-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,

    "MuonEG_Run2022C-22Sep2023"     : { "Era" : "C" , "histAxisName" : "2022_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022/MuonEG_Run2022C-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2022D-22Sep2023"     : { "Era" : "D" , "histAxisName" : "2022_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022/MuonEG_Run2022D-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
}

data_2022EE = {
    "Muon_Run2022E-22Sep2023"   : { "Era" : "E" , "histAxisName" : "2022EE_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022EE/Muon_Run2022E-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_Run2022F-22Sep2023"   : { "Era" : "F" , "histAxisName" : "2022EE_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022EE/Muon_Run2022F-22Sep2023-v2_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_Run2022G-22Sep2023"   : { "Era" : "G" , "histAxisName" : "2022EE_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022EE/Muon_Run2022G-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,

    "EGamma_Run2022E-22Sep2023" : { "Era" : "E" , "histAxisName" : "2022EE_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022EE/EGamma_Run2022E-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_Run2022F-22Sep2023" : { "Era" : "F" , "histAxisName" : "2022EE_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022EE/EGamma_Run2022F-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_Run2022G-22Sep2023" : { "Era" : "G" , "histAxisName" : "2022EE_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022EE/EGamma_Run2022G-22Sep2023-v2_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,

    "MuonEG_Run2022E-22Sep2023" : { "Era" : "E" , "histAxisName" : "2022EE_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022EE/MuonEG_Run2022E-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2022F-22Sep2023" : { "Era" : "F" , "histAxisName" : "2022EE_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022EE/MuonEG_Run2022F-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "MuonEG_Run2022G-22Sep2023" : { "Era" : "G" , "histAxisName" : "2022EE_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2022EE/MuonEG_Run2022G-22Sep2023-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
}

data_2023 = {
    "Muon_0_V1_Run2023C-22Sep2023"   : {"Era" : "C1" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/Muon0_Run2023C-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_0_V2_Run2023C-22Sep2023"   : {"Era" : "C2" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/Muon0_Run2023C-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_0_V3_Run2023C-22Sep2023"   : {"Era" : "C3" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/Muon0_Run2023C-22Sep2023_v3-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_0_V4_Run2023C-22Sep2023"   : {"Era" : "C4" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/Muon0_Run2023C-22Sep2023_v4-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_1_V1_Run2023C-22Sep2023"   : {"Era" : "C1" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/Muon1_Run2023C-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_1_V2_Run2023C-22Sep2023"   : {"Era" : "C2" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/Muon1_Run2023C-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_1_V3_Run2023C-22Sep2023"   : {"Era" : "C3" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/Muon1_Run2023C-22Sep2023_v3-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_1_V4_Run2023C-22Sep2023"   : {"Era" : "C4" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/Muon1_Run2023C-22Sep2023_v4-v2_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,

    "EGamma_0_V1_Run2023C-22Sep2023" : {"Era" : "C1" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/EGamma0_Run2023C-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_0_V2_Run2023C-22Sep2023" : {"Era" : "C2" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/EGamma0_Run2023C-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_0_V3_Run2023C-22Sep2023" : {"Era" : "C3" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/EGamma0_Run2023C-22Sep2023_v3-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_0_V4_Run2023C-22Sep2023" : {"Era" : "C4" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/EGamma0_Run2023C-22Sep2023_v4-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_1_V1_Run2023C-22Sep2023" : {"Era" : "C1" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/EGamma1_Run2023C-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_1_V2_Run2023C-22Sep2023" : {"Era" : "C2" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/EGamma1_Run2023C-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_1_V3_Run2023C-22Sep2023" : {"Era" : "C3" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/EGamma1_Run2023C-22Sep2023_v3-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_1_V4_Run2023C-22Sep2023" : {"Era" : "C4" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/EGamma1_Run2023C-22Sep2023_v4-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,

    "MuonEG_V1_Run2023C-22Sep2023"   : {"Era" : "C1" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/MuonEG_Run2023C-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "MuonEG_V2_Run2023C-22Sep2023"   : {"Era" : "C2" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/MuonEG_Run2023C-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "MuonEG_V3_Run2023C-22Sep2023"   : {"Era" : "C3" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/MuonEG_Run2023C-22Sep2023_v3-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "MuonEG_V4_Run2023C-22Sep2023"   : {"Era" : "C4" ,"histAxisName" : "2023_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023/MuonEG_Run2023C-22Sep2023_v4-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
}

data_2023BPix = {
    "Muon_0_V1_Run2023D-22Sep2023"   : {"Era" : "D1" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/Muon0_Run2023D-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_0_V2_Run2023D-22Sep2023"   : {"Era" : "D2" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/Muon0_Run2023D-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_1_V1_Run2023D-22Sep2023"   : {"Era" : "D1" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/Muon1_Run2023D-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "Muon_1_V2_Run2023D-22Sep2023"   : {"Era" : "D2" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/Muon1_Run2023D-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,

    "EGamma_0_V1_Run2023D-22Sep2023" : {"Era" : "D1" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/EGamma0_Run2023D-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_0_V2_Run2023D-22Sep2023" : {"Era" : "D2" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/EGamma0_Run2023D-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_1_V1_Run2023D-22Sep2023" : {"Era" : "D1" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/EGamma1_Run2023D-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "EGamma_1_V2_Run2023D-22Sep2023" : {"Era" : "D2" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/EGamma1_Run2023D-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,

    "MuonEG_V1_Run2023D-22Sep2023"   : {"Era" : "D1" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/MuonEG_Run2023D-22Sep2023_v1-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
    "MuonEG_V2_Run2023D-22Sep2023"   : {"Era" : "D2" ,"histAxisName" : "2023BPix_data" , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/data/2023BPix/MuonEG_Run2023D-22Sep2023_v2-v1_NANOAOD_WWZ_Dilepton_Run3_09102024/" , "xsecName" : "data" , } ,
}

############################ Bkg samples ############################
central_2016APV_bkg_dict = {
    "UL16APV_TTTo2L2Nu"            : { "histAxisName" : "UL16APV_TTTo2L2Nu"            , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "TTTo2L2Nu" , } ,
    "UL16APV_WJetsToLNu"           : { "histAxisName" : "UL16APV_WJetsToLNu"           , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v2_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "WJetsToLNu" , } ,
    "UL16APV_WWTo2L2Nu"            : { "histAxisName" : "UL16APV_WWTo2L2Nu"            , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024", "xsecName" : "WWTo2L2Nu" , } ,
}
central_2016_bkg_dict = {
    "UL16_TTTo2L2Nu"               : { "histAxisName" : "UL16_TTTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "TTTo2L2Nu" , } ,
    "UL16_WJetsToLNu"              : { "histAxisName" : "UL16_WJetsToLNu"              , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v2_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "WJetsToLNu" , } ,
    "UL16_WWTo2L2Nu"               : { "histAxisName" : "UL16_WWTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "WWTo2L2Nu" , } ,
}
central_2017_bkg_dict = {
    "UL17_TTTo2L2Nu"               : { "histAxisName" : "UL17_TTTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "TTTo2L2Nu" , } ,
    "UL17_WJetsToLNu"              : { "histAxisName" : "UL17_WJetsToLNu"              , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "WJetsToLNu" , } ,
    "UL17_WWTo2L2Nu"               : { "histAxisName" : "UL17_WWTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "WWTo2L2Nu" , } ,
}
central_2018_bkg_dict = {
    "UL18_TTTo2L2Nu"               : { "histAxisName" : "UL18_TTTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "TTTo2L2Nu" , } ,
    "UL18_WJetsToLNu"              : { "histAxisName" : "UL18_WJetsToLNu"              , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "WJetsToLNu" , } ,
    "UL18_WWTo2L2Nu"               : { "histAxisName" : "UL18_WWTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "WWTo2L2Nu" , } ,
}
central_2022_bkg_dict = {
    "2022_TTTo2L2Nu"               : { "histAxisName" : "2022_TTTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "TTto2L2Nu_13p6TeV" , } ,
    "2022_WJetsToLNu"              : { "histAxisName" : "2022_WJetsToLNu"              , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2022_WWTo2L2Nu"               : { "histAxisName" : "2022_WWTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}
central_2022EE_bkg_dict = {
    "2022EE_TTTo2L2Nu"             : { "histAxisName" : "2022EE_TTTo2L2Nu"             , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "TTto2L2Nu_13p6TeV" , } ,
    "2022EE_WJetsToLNu"            : { "histAxisName" : "2022EE_WJetsToLNu"            , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2022EE_WWTo2L2Nu"             : { "histAxisName" : "2022EE_WWTo2L2Nu"             , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}
central_2023_bkg_dict = {
    "2023_TTTo2L2Nu"               : { "histAxisName" : "2023_TTTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "TTto2L2Nu_13p6TeV" , } ,
    "2023_WJetsToLNu"              : { "histAxisName" : "2023_WJetsToLNu"              , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2023_WWTo2L2Nu"               : { "histAxisName" : "2023_WWTo2L2Nu"               , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v4_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}
central_2023BPix_bkg_dict = {
    "2023BPix_TTTo2L2Nu"           : { "histAxisName" : "2023BPix_TTTo2L2Nu"           , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "TTto2L2Nu_13p6TeV" , } ,
    "2023BPix_WJetsToLNu"          : { "histAxisName" : "2023BPix_WJetsToLNu"          , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "WtoLNu-2Jets_13p6TeV" , } ,
    "2023BPix_WWTo2L2Nu"           : { "histAxisName" : "2023BPix_WWTo2L2Nu"           , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "WWto2L2Nu_13p6TeV" , } ,
}

############################ Signal samples ############################

central_2016APV_sig_dict = {
    "UL16APV_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "UL16APV_DYJetsToLL_M_50_MLM"       , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "DYJetsToLL_M_50_MLM", } ,
}
central_2016_sig_dict = {
    "UL16_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "UL16_DYJetsToLL_M_50_MLM"             , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "DYJetsToLL_M_50_MLM", } ,
}
central_2017_sig_dict = {
    "UL17_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "UL17_DYJetsToLL_M_50_MLM"             , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "DYJetsToLL_M_50_MLM", } ,
}
central_2018_sig_dict = {
    "UL18_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "UL18_DYJetsToLL_M_50_MLM"             , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run2_09112024/mc/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_WWZ_Dilepton_Run2_09112024/", "xsecName" : "DYJetsToLL_M_50_MLM", } ,
}
central_2022_sig_dict = {
    "2022_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2022_DYJetsToLL_M_50_MLM"             , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}
central_2022EE_sig_dict = {
    "2022EE_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2022EE_DYJetsToLL_M_50_MLM"         , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}
central_2023_sig_dict = {
    "2023_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2023_DYJetsToLL_M_50_MLM"             , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v1_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
}
central_2023BPix_sig_dict = {
    "2023BPix_DYJetsToLL_M_50_MLM"     : { "histAxisName" : "2023BPix_DYJetsToLL_M_50_MLM"     , "path" : "/store/user/mdittric/skim/WWZ_Dilepton_Run3_09102024/mc/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3_NANOAODSIM_WWZ_Dilepton_Run3_09102024/", "xsecName" : "DYto2L-2Jets_MLL-50_13p6TeV", } ,
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
            exit()
    else:
        print("Failed: None")


############################ Main ############################

# Uncomment the make_jsons_for_dict_of_samples() lines for the jsons you want to make/remake
def main():

    # Specify output paths
    jsons_path = "../../input_samples/sample_jsons/"
    out_dir_data_16APV = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2016APV/")
    out_dir_data_16 = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2016/")
    out_dir_data_17 = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2017/")
    out_dir_data_18 = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2018/")
    out_dir_data_22 = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2022/")
    out_dir_data_22EE = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2022EE/")
    out_dir_data_23 = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2023/")
    out_dir_data_23BPix = os.path.join(jsons_path,"dilepton_analysis_skims/data_samples/2023BPix/")
    out_dir_bkg = os.path.join(jsons_path,"dilepton_analysis_skims/bkg_samples/")
    out_dir_sig = os.path.join(jsons_path,"dilepton_analysis_skims/sig_samples/")

    # Make configs for data samples
    make_jsons_for_dict_of_samples(data_UL16APV, "/ceph/cms/","2016APV", out_dir_data_16APV,era_op=1,on_das=False)
    make_jsons_for_dict_of_samples(data_UL16, "/ceph/cms/","2016", out_dir_data_16,era_op=1,on_das=False)
    make_jsons_for_dict_of_samples(data_UL17, "/ceph/cms/","2017", out_dir_data_17,era_op=1,on_das=False)
    make_jsons_for_dict_of_samples(data_UL18, "/ceph/cms/","2018", out_dir_data_18,era_op=1,on_das=False)
    #make_jsons_for_dict_of_samples(data_2022, "/ceph/cms/","2022", out_dir_data_22,era_op=1,on_das=False)
    #make_jsons_for_dict_of_samples(data_2022EE, "/ceph/cms/","2022EE", out_dir_data_22EE,era_op=1,on_das=False)
    make_jsons_for_dict_of_samples(data_2023, "/ceph/cms/","2023", out_dir_data_23,era_op=1,on_das=False)
    make_jsons_for_dict_of_samples(data_2023BPix, "/ceph/cms/","2023BPix", out_dir_data_23BPix,era_op=1,on_das=False)

    # Make configs for bkg samples
    make_jsons_for_dict_of_samples(central_2016APV_bkg_dict, "/ceph/cms/","2016APV", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2016_bkg_dict, "/ceph/cms/","2016", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2017_bkg_dict, "/ceph/cms/","2017", out_dir_bkg,on_das=False)
    make_jsons_for_dict_of_samples(central_2018_bkg_dict, "/ceph/cms/","2018", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_2022_bkg_dict, "/ceph/cms/","2022", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_2022EE_bkg_dict, "/ceph/cms/","2022EE", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_2023_bkg_dict, "/ceph/cms/","2023", out_dir_bkg,on_das=False)
    #make_jsons_for_dict_of_samples(central_2023BPix_bkg_dict, "/ceph/cms/","2023BPix", out_dir_bkg,on_das=False)

    # Make configs for sig samples
    make_jsons_for_dict_of_samples(central_2016APV_sig_dict, "/ceph/cms/","2016APV", out_dir_sig,on_das=False)
    make_jsons_for_dict_of_samples(central_2016_sig_dict, "/ceph/cms/","2016", out_dir_sig,on_das=False)
    make_jsons_for_dict_of_samples(central_2017_sig_dict, "/ceph/cms/","2017", out_dir_sig,on_das=False)
    make_jsons_for_dict_of_samples(central_2018_sig_dict, "/ceph/cms/","2018", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_2022_sig_dict, "/ceph/cms/","2022", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_2022EE_sig_dict, "/ceph/cms/","2022EE", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_2023_sig_dict, "/ceph/cms/","2023", out_dir_sig,on_das=False)
    #make_jsons_for_dict_of_samples(central_2023BPix_sig_dict, "/ceph/cms/","2023BPix", out_dir_sig,on_das=False)

if __name__ == "__main__":
    main()
