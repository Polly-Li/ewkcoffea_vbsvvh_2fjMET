# Some example run commands for R3 configs

# Run at scale (with futures) DO NOT RUN THIS ON LOGIN NODE AT UF!
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_run3.cfg -o wwz_run3_histos_nosyst -x futures -n 200
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_run3.cfg -o wwz_run3_histos_withsyst -x futures -n 200 --hist-list njets njets_counts --do-systs

# Run with the siphon on DO NOT RUN THIS ON LOGIN NODE AT UF!
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_run3.cfg -o wwz_run3_histos_siphon -x futures -n 200 --siphon

# Run at scale (with wq)
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_run3.cfg -o wwz_run3_histos --hist-list few


time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/samples_4lskim_2022_v1.cfg -o wwz_2022_JECTEST -x futures -n 200 --do-systs --hist-list njets
