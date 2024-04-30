# Some example commands for running at scale

# Run at scale (with wq)
time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples.cfg,../../input_samples/cfgs/wwz_analysis/data_samples.cfg -o wwz_histos --do-systs --hist-list few

# Run at scale (with wq) with hpg copy
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples_hpg.cfg,../../input_samples/cfgs/wwz_analysis/data_samples_hpg.cfg -o wwz_histos_hpg --do-systs --hist-list few
