# Some example run commands for R3 configs
# ALL RUN SCRIPTS ARE TO BE RUN LOCALLY AT UAF!

# Run at scale (with futures) 
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run3_v1.cfg -o wwz_run3_histos -x futures -n 200
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run3_v1.cfg -o wwz_run3_histos_withsyst -x futures -n 200 --do-systs

# Run with the siphon on 
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run3_v1.cfg -o wwz_run3_histos_siphon -x futures -n 200 --siphon

# Run at scale (with wq)
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run3_v1.cfg -o wwz_run3_histos --hist-list few

#Run at scale (with futures) for JUST 2022
time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_2022_v1.cfg -o wwz_2022_histos -x futures -n 200

# Run at scale (with futures) for JUST 2023
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_2023_v1.cfg -o wwz_2023_histos -x futures -n 200
