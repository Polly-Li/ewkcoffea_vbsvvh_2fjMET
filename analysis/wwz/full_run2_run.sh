# Some example commands for running at scale

# Run at scale (with futures) DO NOT RUN THIS ON LOGIN NODE!
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples.cfg,../../input_samples/cfgs/wwz_analysis/data_samples.cfg -o wwz_histos -x futures -n 128

# Run at scale (with wq)
time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples.cfg,../../input_samples/cfgs/wwz_analysis/data_samples.cfg -o wwz_histos --do-systs --hist-list few

# Run at scale (with wq) with hpg copy
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples_hpg.cfg,../../input_samples/cfgs/wwz_analysis/data_samples_hpg.cfg -o wwz_histos_hpg --do-systs --hist-list few

# Run at scale (with 4l skim) with uaf-2 copy
# Run only nominal with many histograms
# time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim.cfg -x futures -n 200 -s 100000000 -o wwz_histos_4l_skim -p /data/userdata/phchang/histos
# Run all systematics with few histograms
# time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim.cfg -x futures -n 200 -s 100000000 --do-systs --hist-list few -o wwz_histos_4l_skim_syst
