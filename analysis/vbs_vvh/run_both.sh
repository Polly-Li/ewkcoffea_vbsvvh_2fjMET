project_name=$1
cutflow_name=$2

python3 "$HOME/analysis/vvh/coffea_scripts/run_analysis.py" "$HOME/analysis/vvh/coffea_scripts/configs/all_year.cfg" --cutflow "${cutflow_name}" --project "${project_name}" -x futures -n 64 
python3 "$HOME/analysis/vvh/coffea_scripts/make_plots.py" --cutflow "${cutflow_name}" --project "${project_name}" -p -y