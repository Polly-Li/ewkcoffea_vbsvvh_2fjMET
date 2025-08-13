import csv

sample_categories = ['sig', 'bkg', 'data']

sample_types = {
    'sig': ['WWH_OS', 'WWH_SS', 'WZH', 'ZZH'],
    'bkg': ["ttbar", "ttx", "ST", "WJets", "ZJets", "EWK", "QCD", "Other"],
    'data': ['data']
}

#choose buckets according to sample_types
sample_names = {
    'WWH_OS': ["VBSWWH_OS_VBSCuts"],
    'WWH_SS': ["VBSWWH_SS_VBSCuts"],
    'WZH': ["VBSWZH_VBSCuts"],
    'ZZH': ["VBSZZH_VBSCuts"],
    "DY": ["DY"],
    "Other": ["Other"],
    "ttx": ["ttx"],
    "WJets": ["WJets"],
    "ZJets": ["ZJets"],
    'QCD': ["QCD_HT50to100", "QCD_HT100to200", "QCD_HT200to300", "QCD_HT300to500", "QCD_HT500to700", "QCD_HT700to1000", "QCD_HT1000to1500", "QCD_HT1500to2000", "QCD_HT2000toInf", "QCD"],
    'ttbar_had': ["TTToHadronic"],
    'ttbar_SL': ["TTToSemiLeptonic"],
    'ttbar': ["TTToHadronic", "TTToSemiLeptonic", "ttbar"],
    'ST': ["ST_t-channel_antitop_4f_InclusiveDecays", "ST_t-channel_top_4f_InclusiveDecays", "ST_tW_antitop_5f_inclusiveDecays", "ST_tW_top_5f_inclusiveDecays", "ST"],
    'WJet': ["WJetsToQQ_HT-200to400", "WJetsToQQ_HT-400to600", "WJetsToQQ_HT-600to800", "WJetsToQQ_HT-800toInf", "WJets"],
    'ZJet': ["ZJetsToQQ_HT-200to400", "ZJetsToQQ_HT-400to600", "ZJetsToQQ_HT-600to800", "ZJetsToQQ_HT-800toInf", "ZJets"],
    'EWKV': ["EWKWminus2Jets_WToQQ_dipoleRecoilOn", "EWKWplus2Jets_WToQQ_dipoleRecoilOn", "EWKZ2Jets_ZToLL_M-50", "EWKZ2Jets_ZToNuNu_M-50", "EWKZ2Jets_ZToQQ_dipoleRecoilOn"],
    "EWK": ["EWK"],
    'ttX': ["TTWW", "TTWZ", "TTWJetsToQQ", "ttHToNonbb_M125", "ttHTobb_M125", "TTbb_4f_TTToHadronic"],
    'bosons': ["VHToNonbb_M125", "WWTo1L1Nu2Q_4f", "WWTo4Q_4f", "WWW_4F", "WWZ_4F", "WZJJ_EWK_InclusivePolarization", "WZTo1L1Nu2Q_4f", "WZTo2Q2L_mllmin4p0", "WZZ", "WminusH_HToBB_WToLNu_M-125", "WplusH_HToBB_WToLNu_M-125", "ZH_HToBB_ZToQQ_M-125", "ZZTo2Nu2Q_5f", "ZZTo2Q2L_mllmin4p0", "ZZTo4Q_5f", "ZZZ"],
    "data": ['MET'],
}

sample_colour_mpl = {
    "data": "blue",
    "sig": "red",
    "WWH_OS": "crimson",
    "WWH_SS": "lightcoral",
    "WZH": "orangered",
    "ZZH": "tomato",
    "bkg_total": "dimgray",
    "DY": "palegreen",
    "ttbar": "mediumturquoise",
    "ttbar_had": "powderblue",
    "ttbar_SL": "deepskyblue",
    "ttx": "cadetblue",
    "ttX": "cadetblue",
    "ST": "lightsteelblue",
    "WJets": "darkorange",
    "WJet": "darkorange",
    "ZJets": "lightyellow",
    "ZJet": "lightyellow",
    "EWKV": "skyblue",
    "EWK": "skyblue",
    "QCD": "lightcyan",
    "bosons": "sandybrown",
    "Other": "orchid",
}

# Assign and store rows
rows = []

for cat_index, sample_cat in enumerate(sample_categories):
    type_list = sample_types[sample_cat]
    for type_index, sample_type in enumerate(type_list):
        name_list = sample_names.get(sample_type, [])
        for name_index, sample_name in enumerate(name_list):
            sample_code = f"{cat_index}{type_index}{name_index}"
            plotting_colour = sample_colour_mpl.get(sample_type, "black")  # fallback to black if missing
            rows.append([sample_cat, sample_type, sample_name, sample_code, plotting_colour])

# Write to CSV
with open('sample_names.csv', mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['sample_category', 'sample_type', 'sample_name', 'sample_code', 'plotting_colour'])
    for row in rows:
        writer.writerow(row)