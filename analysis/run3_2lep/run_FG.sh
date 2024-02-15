#Run all samples

trap "kill 0" SIGINT

#Data
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/EGamma_EraF_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/EGamma_EraG_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/MuonEG_EraF_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/MuonEG_EraG_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/Muon_EraF_2L_Loose_Skim.json,../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/Muon_EraG_2L_Loose_Skim.json -x futures -o Data_Out -n 125 &
wait

#TT
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/TT_EE_2L_Loose_Skim.json -x futures -o TT_Out -n 125 &
wait

#WW
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/WW_EE_2L_Loose_Skim.json -x futures -o WW_Out -n 125 &
wait

#W
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/WtoLNu_EE_2L_Loose_Skim.json -x futures -o W_Out -n 125 &
wait

#DY
python run_run3_2lep.py ../../input_samples/sample_jsons/run3_2lep_skims/FG_Skims/DY_EE_2L_Loose_Skim.json -x futures -o DY_Out -n 125 &
wait

#Make Root Files
python pickleroot_plotter.py Data_Out.pkl.gz DY_Out.pkl.gz TT_Out.pkl.gz WW_Out.pkl.gz W_Out.pkl.gz -r histos/
wait

#Make the Plots (Change directory name as needed)
python plot.py -dir DiLep_EraFG_2Lep_LIN -tt TT_Out.root -ww WW_Out.root -w W_Out.root -dy DY_Out.root -dt Data_Out.root
wait

#Publish them to UAF
publish DiLep_EraFG_2Lep_LIN/

