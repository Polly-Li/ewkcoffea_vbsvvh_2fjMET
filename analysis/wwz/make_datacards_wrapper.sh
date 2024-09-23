# Wrapper around making full analysis test data cards

R2_PKL="r2_wwz_histos_withSyst.pkl.gz"
R3_PKL="r3_wwz_histos_withSyst.pkl.gz"

rm cards_wwz4l/*

python make_datacards.py histos/$R3_PKL -u run3 -s --bdt
python make_datacards.py histos/$R2_PKL -u run2 -s --bdt
