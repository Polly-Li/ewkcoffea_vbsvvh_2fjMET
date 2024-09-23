# Some example commands for running at scale

### Examples with 4l skims configs ### (ONLY WORKS AT UAF. SAMPLES HAVE NOT BEEN MOVED TO HPG!)

# Run at scale with futures
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run2_v1.cfg -x futures -n 200 -s 100000000 -o wwz_histos --hist-list bdt
time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run2_v1.cfg,../../input_samples/cfgs/wwz_analysis/samples_4lskim_2022_v1.cfg,../../input_samples/cfgs/wwz_analysis/samples_4lskim_2023_v1.cfg -x futures -n 200 -s 100000000 -o run2run3_JERC_Total_syst.pkl.gz --hist-list njets njets_counts --do-systs

