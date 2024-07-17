# Some example run commands for R3 configs

# Run at scale (with futures) DO NOT RUN THIS ON LOGIN NODE!
time python run_dilepton.py ../../input_samples/cfgs/dilepton_skims/mc_sig_bkg_samples_2023.cfg,../../input_samples/cfgs/dilepton_skims/data_samples_2023.cfg -o dilepton_2023 -x futures -n 125

# Run at scale (with wq)
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples_run3.cfg,../../input_samples/cfgs/wwz_analysis/data_samples_run3.cfg -o wwz_run3_histos --hist-list few
