# Wrapper around making full analysis test data cards

R2_PKL="r2_wwz_histos_withSyst.pkl.gz"
R3_PKL="r3_wwz_histos_withSyst.pkl.gz"

rm cards_wwz4l/*

# For unblind BDT results
python make_datacards.py histos/$R2_PKL -u run2 -s --bdt --unblind --zero-low-mc
python make_datacards.py histos/$R3_PKL -u run3 -s --bdt --unblind --zero-low-mc

# For blind BDT results
#python make_datacards.py histos/$R2_PKL -u run2 -s --bdt --zero-low-mc
#python make_datacards.py histos/$R3_PKL -u run3 -s --bdt --zero-low-mc

# For Cut Based results
#python make_datacards.py histos/$R2_PKL -u run2 -s
#python make_datacards.py histos/$R3_PKL -u run3 -s
