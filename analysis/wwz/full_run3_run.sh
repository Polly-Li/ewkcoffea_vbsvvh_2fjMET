# Some example run commands for R3 configs
# ALL RUN SCRIPTS ARE TO BE RUN LOCALLY AT UAF!

#################################################################################################################################
# Run at scale (with futures) 
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run3_v1.cfg -o wwz_run3_histos -x futures -n 200
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run3_v1.cfg -o wwz_run3_histos_withsyst -x futures -n 200 --do-systs

# Run with the siphon on 
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run3_v1.cfg -o wwz_run3_histos_siphon -x futures -n 200 --siphon

# Run at scale (with wq)
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run3_v1.cfg -o wwz_run3_histos --hist-list few
#################################################################################################################################

#Run at scale (with futures) for JUST 2022
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_2022_v1.cfg -o y22_wwz_histos_noSyst -x futures -n 200 --hist-list bdt
time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_2022_v1.cfg -o y22_wwz_histos_withSyst -x futures -n 200 --do-systs --hist-list njets njets_counts

# Run at scale (with futures) for JUST 2023
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_2023_v1.cfg -o y23_wwz_histos_noSyst -x futures -n 200

#Run at scale (with futures) for scaled 2022 (This is 2022 MC + 2022 & 2023 Data)
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_scaled_run3.cfg -o y22scaled_wwz_histos_noSyst -x futures -n 200

#Run at scale (with futures) for Run 2 + Run 3 (This is scaled 2022 as explained above)
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4l_run2run3.cfg -o wwz_run2run3_histos -x futures -n 200
