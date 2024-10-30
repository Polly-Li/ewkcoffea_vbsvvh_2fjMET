# Some example commands for running at scale


### Examples with 3l skim configs (OLD) ###

# Run at scale with futures DO NOT RUN THIS ON LOGIN NODE!
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples.cfg,../../input_samples/cfgs/wwz_analysis/data_samples.cfg -o wwz_histos -x futures -n 128

# Run at scale with wq
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples.cfg,../../input_samples/cfgs/wwz_analysis/data_samples.cfg -o wwz_histos --do-systs --hist-list bdt

# Run at scale with wq, with hpg copy
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples_hpg.cfg,../../input_samples/cfgs/wwz_analysis/data_samples_hpg.cfg -o wwz_histos_hpg --do-systs --hist-list njets njets_counts


### Examples with 4l skims configs ### (ONLY WORKS AT UAF. SAMPLES HAVE NOT BEEN MOVED TO HPG!)

# Run at scale with futures
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run2_v1.cfg -x futures -n 200 -s 100000000 -o r2_wwz_histos_noSyst --hist-list bdt
time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run2_v1.cfg -x futures -n 200 -s 100000000 -o r2_wwz_histos_withSyst --hist-list njets njets_counts --do-systs

# Run with the siphon turned on (might want to comment out data in the input cfg)
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run2_v1.cfg -x futures -n 200 -o r2_wwz_histos_noSyst_siphon --siphon

# Run at scale with wq
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run2_v1.cfg -o r2_wwz_histos_noSyst --hist-list bdt
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_run2_v1.cfg -o r2_wwz_histos_withSyst --hist-list njets njets_counts --do-systs
