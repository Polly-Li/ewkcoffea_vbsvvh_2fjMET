# Run at scale (with wq)
time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples.cfg,../../input_samples/cfgs/wwz_analysis/data_samples.cfg -o wwz_histos

# Run at scale (with futures) DO NOT RUN THIS ON LOGIN NODE!
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples.cfg,../../input_samples/cfgs/wwz_analysis/data_samples.cfg -o wwz_histos -x futures -n 128
