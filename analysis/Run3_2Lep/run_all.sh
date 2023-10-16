#Run all samples

trap "kill 0" SIGINT

#Data
python run_Run3_2Lep.py Run3_JSONS/Muon_2022G.json,Run3_JSONS/MuonEG_2022G.json,Run3_JSONS/EGamma_2022G.json -x futures -o Data_Out -n 127 &
wait

#TT
python run_Run3_2Lep.py Run3_JSONS/Run3_TTto2L2Nu.json -x futures -o TT_Out -n 127 &
wait

#WW
python run_Run3_2Lep.py Run3_JSONS/Run3_WW.json -x futures -o WW_Out -n 127 &
wait

#W
python run_Run3_2Lep.py Run3_JSONS/Run3_WtoLNu_2Jets.json -x futures -o W_Out -n 127 &
wait

#DY
python run_Run3_2Lep.py Run3_JSONS/Run3_DYto2L_2Jets_MLL50.json -x futures -o DY_Out -n 127 &
wait

#Make Root Files
python pickleroot_plotter.py Data_Out.pkl.gz DY_Out.pkl.gz TT_Out.pkl.gz WW_Out.pkl.gz W_Out.pkl.gz -r histos/
wait

#Make the Plots (Change directory name as needed)
python plot.py -dir DiLep_PU_15102023_1 -tt TT_Out.root -ww WW_Out.root -w W_Out.root -dy DY_Out.root -dt Data_Out.root
wait

#Publish them to UAF
publish DiLep_PU_15102023_1/

