# This script reproduces the reference yields file that the CI compares against
# Run this script when you want to update the reference file for 2022

# Get the file the CI uses, and move it to the directory the JSON expects
printf "\nDownloading root file...\n"
wget -nc http://uaf-10.t2.ucsd.edu/~mdittric/for_ci/for_wwz/WWZto4L2Nu_4F_TuneCP5_13p6TeV_amcatnlo-pythia8_Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2_NANOAODSIM_WWZ_MC_2024_0811/output_9.root

# Run the processor
printf "\nRunning the processor...\n"
time python run_wwz4l.py ../../input_samples/sample_jsons/test_samples/2022_WWZJetsTo4L2Nu_forCI.json -x iterative -o new_ref_histos_2022

# Make the JSON file of the yields
printf "\nMaking the yields JSON file...\n"
python get_wwz_counts.py -f histos/new_ref_histos_2022.pkl.gz -n new_ref_yields_2022 -s 2022_WWZJetsTo4L2Nu

# Compare the JSON file of the yields
printf "\nCompare the new yields JSON file to old ref...\n"
python comp_json_yields.py new_ref_yields_2022.json ref_for_ci/counts_wwz_ref_2022.json -t1 "New yields" -t2 "Old ref yields"

# Replace the reference yields with the new reference yields
printf "\nReplacing 2022 ref yields JSON with new file...\n"
mv new_ref_yields_2022.json ref_for_ci/counts_wwz_ref_2022.json
printf "\n\nDone.\n\n"
