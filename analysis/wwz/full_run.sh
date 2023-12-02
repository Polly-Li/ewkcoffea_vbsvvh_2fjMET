# An example of running at scale with either WQ or futures

# Run at scale (with wq)
# To submit workers with wq, e.g.:
#   slurm_submit_workers --cores 64 --memory 500000 -M ${USER}-workqueue-coffea -p  "--account avery --qos avery-b --time 1:00:00" 1
time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples.cfg,../../input_samples/cfgs/wwz_analysis/data_samples.cfg -o wwz_histos

# Run at scale (with futures)
# Make you have srun before running at scale with futures (DO NOT RUN THIS ON LOGIN NODE!), e.g. like this:
#    srun -t 600 --qos=avery --account=avery --cpus-per-task=128 --mem=512gb --pty bash -i
#time python run_wwz4l.py ../../input_samples/cfgs/wwz_analysis/mc_sig_bkg_samples.cfg,../../input_samples/cfgs/wwz_analysis/data_samples.cfg -o wwz_histos -x futures -n 128
