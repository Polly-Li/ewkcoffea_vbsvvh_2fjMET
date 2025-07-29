# Example run commands

# Example for the from-Nano processor
#time python run_analysis.py input_cfg_r2.cfg -x futures -n 64 -o vvh_hists_from_nano -p 1lep1jf_nano

# Example for the from-RDF processor
time python run_analysis.py ../../input_samples/sample_jsons/vbs_vvh/rdf_input_1lep1fj_mc.json -x futures -n 32 -o vvh_hists --hist-list njets njets_counts fj0_mparticlenet -p 1lep1jf
