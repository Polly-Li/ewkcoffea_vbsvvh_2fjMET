#Run all samples

trap "kill 0" SIGINT

#-r root://redirector.t2.ucsd.edu:1095/

#Data
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/DoubleMuon_B_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/DoubleMuon_C_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/EGamma_B_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/EGamma_C_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/EGamma_D_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/EGamma_E_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/EGamma_F_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/EGamma_G_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/Muon_C_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/Muon_D_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/Muon_E_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/MuonEG_B_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/MuonEG_C_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/MuonEG_D_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/MuonEG_E_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/MuonEG_F_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/MuonEG_G_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/Muon_F_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/Muon_G_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/SingleMuon_B_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/data_skims_hpg/SingleMuon_C_2L_Loose_Skim.json  -x futures -o Data_Out -n 125 &
wait

#TT
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_dilepton_skims/mc_skims_hpg/TTto2L2Nu_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/mc_skims_hpg/TTto2L2Nu_EE_2L_Loose_Skim.json -x futures -o TT_Out -n 125 &
wait

#WW
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_dilepton_skims/mc_skims_hpg/WW_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/mc_skims_hpg/WW_EE_2L_Loose_Skim.json -x futures -o WW_Out -n 125 &
wait

#W
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_dilepton_skims/mc_skims_hpg/WtoLNu_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/mc_skims_hpg/WtoLNu_EE_2L_Loose_Skim.json -x futures -o W_Out -n 125 &
wait

#DY
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_dilepton_skims/mc_skims_hpg/DY_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_dilepton_skims/mc_skims_hpg/DY_EE_2L_Loose_Skim.json -x futures -o DY_Out -n 125 &
wait

#Make Root Files
python pickleroot_plotter.py Data_Out.pkl.gz DY_Out.pkl.gz TT_Out.pkl.gz WW_Out.pkl.gz W_Out.pkl.gz -r histos/
wait

#Make the Plots (Change directory name as needed)
python plot.py -dir DiLep_Ele_Tight -tt TT_Out.root -ww WW_Out.root -w W_Out.root -dy DY_Out.root -dt Data_Out.root
wait

#Publish them to UAF
publish DiLep_Ele_Tight/

