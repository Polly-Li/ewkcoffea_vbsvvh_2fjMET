#Run all samples

trap "kill 0" SIGINT

#Data
python run_Run3_2Lep.py Run3_JSONS/2022EE_EGamma_E_2l.json,Run3_JSONS/2022EE_EGamma_F_2l.json,Run3_JSONS/2022EE_EGamma_G_2l.json,Run3_JSONS/2022EE_Muon_E_2l.json,Run3_JSONS/2022EE_MuonEG_E_2l.json,Run3_JSONS/2022EE_MuonEG_F_2l.json,Run3_JSONS/2022EE_MuonEG_G_2l.json,Run3_JSONS/2022EE_Muon_F_2l.json,Run3_JSONS/2022EE_Muon_G_2l.json,Run3_JSONS/2022_EGamma_C_2l.json,Run3_JSONS/2022_EGamma_D_2l.json,Run3_JSONS/2022_Muon_C_2l.json,Run3_JSONS/2022_Muon_D_2l.json,Run3_JSONS/2022_MuonEG_C_2l.json,Run3_JSONS/2022_MuonEG_D_2l.json -x futures -o Data_Out -n 60 &
wait

#TT
python run_Run3_2Lep.py Run3_JSONS/2022EE_TT_2l.json,Run3_JSONS/2022_TT_2l.json -x futures -o TT_Out -n 60 &
wait

#WW
python run_Run3_2Lep.py Run3_JSONS/2022EE_WW_2l.json,Run3_JSONS/2022_WW_2l.json -x futures -o WW_Out -n 60 &
wait

#W
python run_Run3_2Lep.py Run3_JSONS/2022EE_W_2l.json,Run3_JSONS/2022_W_2l.json -x futures -o W_Out -n 60 &
wait

#DY
python run_Run3_2Lep.py Run3_JSONS/2022EE_DY_2l.json,Run3_JSONS/2022_DY_2l.json -x futures -o DY_Out -n 60 &
wait

#Make Root Files
python pickleroot_plotter.py Data_Out.pkl.gz DY_Out.pkl.gz TT_Out.pkl.gz WW_Out.pkl.gz W_Out.pkl.gz -r histos/
wait

#Make the Plots (Change directory name as needed)
python plot.py -dir 2022_ALL_DiLep_03122023_2 -tt TT_Out.root -ww WW_Out.root -w W_Out.root -dy DY_Out.root -dt Data_Out.root
wait

#Publish them to UAF
publish 2022_ALL_DiLep_03122023_2/

