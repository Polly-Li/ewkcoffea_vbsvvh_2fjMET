# Filter flags
def get_filter_flag_mask(events, year, is2022,is2023):

    # Filters
    filter_flags = events.Flag
    if (is2022 or is2023):
        filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & filter_flags.ecalBadCalibFilter & filter_flags.BadPFMuonDzFilter & filter_flags.hfNoisyHitsFilter & filter_flags.eeBadScFilter
    elif year in ["2016","2016APV"]:
        filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & filter_flags.BadPFMuonDzFilter & filter_flags.eeBadScFilter
    else:
        filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & filter_flags.BadPFMuonDzFilter & filter_flags.eeBadScFilter & filter_flags.ecalBadCalibFilter

    return filters
