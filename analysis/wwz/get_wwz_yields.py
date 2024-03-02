import argparse
import pickle
import json
import gzip
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import copy
import hist

from topcoffea.scripts.make_html import make_html
from topcoffea.modules import utils
import topcoffea.modules.MakeLatexTable as mlt

import ewkcoffea.modules.yield_tools as yt
import ewkcoffea.modules.sample_groupings as sg


# This script opens a pkl file of histograms produced by wwz processor
# Reads the histograms and dumps out the yields for each group of processes
# Example usage: python get_yld_check.py histos/tmp_histo.pkl.gz -y

# Colors in VVV observation
#ZZ    = (240, 155, 205)  #F09B9B
#ttZ   = (0, 208, 145) #00D091
#WZ    = (163, 155, 47) #A39B2F
#tWZ   = (205, 240, 155) #CDF09B
#Other = (205, 205, 205) #CDCDCD
CLR_LST = ["red","blue","#F09B9B","#00D091","#CDF09B","#A39B2F","#CDCDCD"]
#CLR_LST = ["#F09B9B","#00D091","#CDF09B"]

SAMPLE_DICT_BASE = sg.SAMPLE_DICT_BASE

# Names of the cut-based and BDT SRs
SR_SF_CB = ["sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C"]
SR_OF_CB = ["sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4"]
SR_SF_BDT = ["sr_4l_bdt_sf_wwz_sr1", "sr_4l_bdt_sf_wwz_sr2", "sr_4l_bdt_sf_wwz_sr3", "sr_4l_bdt_sf_wwz_sr4", "sr_4l_bdt_sf_zh_sr1", "sr_4l_bdt_sf_zh_sr2", "sr_4l_bdt_sf_zh_sr3"]
SR_OF_BDT = ["sr_4l_bdt_of_wwz_sr1", "sr_4l_bdt_of_wwz_sr2", "sr_4l_bdt_of_wwz_sr3", "sr_4l_bdt_of_wwz_sr4", "sr_4l_bdt_of_zh_sr1", "sr_4l_bdt_of_zh_sr2", "sr_4l_bdt_of_zh_sr3", "sr_4l_bdt_of_zh_sr4",]

BDT_INPUT_LST = [
    "mll_wl0_wl1",
    "dphi_4l_met",
    "dphi_zleps_met",
    "dphi_wleps_met",
    "dr_wl0_wl1",
    "dr_zl0_zl1",
    "dr_wleps_zleps",
    "met",
    "mt2",
    "ptl4",
    "scalarptsum_lepmet",
    "scalarptsum_lepmetjet",
    "z_lep0_pt",
    "z_lep1_pt",
    "w_lep0_pt",
    "w_lep1_pt",
]

TMP_VAR_LST = [
    "j0pt",
    "njets",
    "nbtagsl",
    "nleps",
    "met",
    "l0pt",
]

# Ref with no  SFs applied (note naming change to presel categories since the time this was created)
# feb27_baseline_noSFnoSys.pkl.gz
#EWK_REF = {'WWZ': {'sr_4l_sf_B': [1.9092815210624394, 5.868219968367227e-05], 'sr_4l_sf_C': [0.5507712043909123, 1.7292573108649608e-05], 'sr_4l_of_3': [1.433722815978399, 4.4002888617222127e-05], 'sr_4l_of_4': [4.9765614088337315, 0.00015460696855228978], 'sr_4l_sf_presel': [8.218720214481436, 0.0002543393204535014], 'sr_4l_sf_trn': [5.38486802954867, 0.00016650880857958835], 'sr_4l_of_presel': [8.977618631754012, 0.00027742868980258986], 'all_events': [143.03255106410325, 0.004385052076200205], '4l_presel': [25.390131740285142, 0.000780364613901527], 'cr_4l_btag_sf_offZ_met80': [1.1789255509174836, 3.543312373031215e-05], 'cr_4l_sf': [1.0234111444715381, 3.1406010820999584e-05], 'sr_4l_sf_A': [2.248098852023759, 6.983801460398903e-05], 'sr_4l_of_1': [0.6239818002213724, 1.9074980915948606e-05], 'sr_4l_of_2': [0.7141527369931282, 2.195676707741774e-05], 'cr_4l_btag_of': [2.332873423842102, 7.052031927015096e-05], 'sr_sf_all_cutbased': [4.7081515774771105, 0.0001458127873963109], 'sr_of_all_cutbased': [7.748418762026631, 0.00023964160516287824], 'sr_all_cutbased': [12.456570339503742, 0.00038545439255918923]}, 'ZH': {'sr_4l_sf_B': [1.4252621539799293, 9.199921721242454e-05], 'sr_4l_sf_C': [0.6304575557069256, 2.7550956936472543e-05], 'sr_4l_of_3': [0.32845933757653256, 1.773571780234943e-05], 'sr_4l_of_4': [0.14063188295949658, 1.0594258701507093e-05], 'sr_4l_sf_presel': [7.500271647765658, 0.00038186652465116205], 'sr_4l_sf_trn': [4.522856767529447, 0.00022087671808088043], 'sr_4l_of_presel': [7.463838904360273, 0.0003795543171914625], 'all_events': [137.65612493509434, 0.008261347318551384], '4l_presel': [19.875361077707566, 0.0010595854754236367], 'cr_4l_btag_sf_offZ_met80': [0.6045969923652592, 4.986340568511177e-05], 'cr_4l_sf': [0.06583264991149917, 2.963363355469417e-06], 'sr_4l_sf_A': [0.8808398754063091, 6.696664634914041e-05], 'sr_4l_of_1': [2.8618885583600786, 0.00013797226190703833], 'sr_4l_of_2': [1.2643676003663131, 6.131047557490622e-05], 'cr_4l_btag_of': [1.6822630345304788, 0.00010680937268688129], 'sr_sf_all_cutbased': [2.936559585093164, 0.0001865168204980375], 'sr_of_all_cutbased': [4.595347379262421, 0.00022761271398580107], 'sr_all_cutbased': [7.531906964355585, 0.00041412953448383865]}, 'ZZ': {'sr_4l_sf_B': [4.06302963554117, 0.002343257632045054], 'sr_4l_sf_C': [2.4055145616730442, 0.0013489537268409197], 'sr_4l_of_3': [0.3522159432868648, 0.00020917560040396212], 'sr_4l_of_4': [0.4632858518671128, 0.00027953866293970814], 'sr_4l_sf_presel': [625.3544098984094, 0.3562166451167309], 'sr_4l_sf_trn': [90.69184155950279, 0.05143268771989755], 'sr_4l_of_presel': [19.06702512973061, 0.011167948980671838], 'all_events': [20839.935963596385, 12.44996510182031], '4l_presel': [2996.3448773944438, 1.6969037274479712], 'cr_4l_btag_sf_offZ_met80': [2.358815064686496, 0.0013679280585803325], 'cr_4l_sf': [1683.1491276815104, 0.9502221381958484], 'sr_4l_sf_A': [1.19626292139219, 0.0006944572452723625], 'sr_4l_of_1': [0.5727219310065266, 0.0003361000024790477], 'sr_4l_of_2': [0.562517372865841, 0.00033284359074704046], 'cr_4l_btag_of': [3.7302420161504415, 0.0021904190096882813], 'sr_sf_all_cutbased': [7.664807118606404, 0.004386668604158336], 'sr_of_all_cutbased': [1.9507410990263452, 0.0011576578565697584], 'sr_all_cutbased': [9.61554821763275, 0.005544326460728095]}, 'ttZ': {'sr_4l_sf_B': [0.9408668631222099, 0.0030715775891913986], 'sr_4l_sf_C': [0.23643477854784578, 0.0008312534314345273], 'sr_4l_of_3': [0.7172879233257845, 0.0020985737304703837], 'sr_4l_of_4': [2.0460971421562135, 0.006911022257176301], 'sr_4l_sf_presel': [3.5659334138035774, 0.01180990449602128], 'sr_4l_sf_trn': [2.4085794143611565, 0.00810546769571046], 'sr_4l_of_presel': [3.787581167416647, 0.012204073027489382], 'all_events': [6400.391729545197, 30.514728026480668], '4l_presel': [155.00430985842831, 0.5408661444449773], 'cr_4l_btag_sf_offZ_met80': [31.04904586810153, 0.10758473272086123], 'cr_4l_sf': [0.4853159709600732, 0.0016168707289103397], 'sr_4l_sf_A': [1.0880661539267749, 0.003780754297758552], 'sr_4l_of_1': [0.26424446457531303, 0.0007768689629915706], 'sr_4l_of_2': [0.2673212867230177, 0.0008240181824348471], 'cr_4l_btag_of': [63.414810959715396, 0.22234178470495916], 'sr_sf_all_cutbased': [2.2653677955968305, 0.007683585318384478], 'sr_of_all_cutbased': [3.2949508167803288, 0.010610483133073103], 'sr_all_cutbased': [5.560318612377159, 0.018294068451457583]}, 'tWZ': {'sr_4l_sf_B': [0.4198845353530487, 0.00010717563775915021], 'sr_4l_sf_C': [0.10839767166180536, 2.8695803200881038e-05], 'sr_4l_of_3': [0.3047727922530612, 7.527883120245118e-05], 'sr_4l_of_4': [0.9282776923646452, 0.00024224380419460127], 'sr_4l_sf_presel': [1.5975971029547509, 0.00040201487686921143], 'sr_4l_sf_trn': [1.123730185412569, 0.0002804039804139507], 'sr_4l_of_presel': [1.730163810585509, 0.0004416391065870169], 'all_events': [178.15963437006576, 0.04486192040271183], '4l_presel': [26.544821587594924, 0.0066209831302612675], 'cr_4l_btag_sf_offZ_met80': [4.7047833817778155, 0.0011600529197435352], 'cr_4l_sf': [0.20248869292845484, 5.108321546004097e-05], 'sr_4l_sf_A': [0.46284899866441265, 0.00011133764734468921], 'sr_4l_of_1': [0.1342859514261363, 3.176125450297847e-05], 'sr_4l_of_2': [0.1488241527986247, 3.757876565006425e-05], 'cr_4l_btag_of': [10.202555451760418, 0.002522789947753728], 'sr_sf_all_cutbased': [0.9911312056792667, 0.00024720908830472044], 'sr_of_all_cutbased': [1.5161605888424674, 0.0003868626555500952], 'sr_all_cutbased': [2.507291794521734, 0.0006340717438548156]}, 'WZ': {'sr_4l_sf_B': [0.4695464987307787, 0.01739613314260788], 'sr_4l_sf_C': [0.18745648488402367, 0.005562775915829661], 'sr_4l_of_3': [0.01870139129459858, 0.01571616535855198], 'sr_4l_of_4': [0.81333789229393, 0.04587245091006516], 'sr_4l_sf_presel': [1.9591068103909492, 0.10288223760163802], 'sr_4l_sf_trn': [0.9570152945816517, 0.04763624418418125], 'sr_4l_of_presel': [1.9694814532995224, 0.11593636356868098], 'all_events': [64295.440694248304, 3785.8715240238676], '4l_presel': [9.901992252096534, 0.5879718428319408], 'cr_4l_btag_sf_offZ_met80': [0.14544222131371498, 0.01005621476728371], 'cr_4l_sf': [0.19009067304432392, 0.012810969442493516], 'sr_4l_sf_A': [0.05745987221598625, 0.006227639598605526], 'sr_4l_of_1': [0.4658583365380764, 0.027649146792590645], 'sr_4l_of_2': [0.13211736641824245, 0.004049024052467726], 'cr_4l_btag_of': [0.7625117842108011, 0.04328537125233069], 'sr_sf_all_cutbased': [0.7144628558307886, 0.02918654865704307], 'sr_of_all_cutbased': [1.4300149865448475, 0.09328678711367551], 'sr_all_cutbased': [2.144477842375636, 0.12247333577071857]}, 'other': {'sr_4l_sf_B': [0.7585755839245394, 0.12648709826160898], 'sr_4l_sf_C': [0.17917197081260383, 0.029325680594073936], 'sr_4l_of_3': [0.06966101261787117, 0.01827217095315412], 'sr_4l_of_4': [0.5390066548716277, 0.027797331479117317], 'sr_4l_sf_presel': [6.328916352707893, 0.7409668954355337], 'sr_4l_sf_trn': [2.1822817949578166, 0.23177169086467067], 'sr_4l_of_presel': [1.8504821928218007, 0.20673090259240398], 'all_events': [2612628.3130944027, 4739155.263059123], '4l_presel': [43.044581399764866, 5.732465444166863], 'cr_4l_btag_sf_offZ_met80': [1.6483099256874993, 0.051794245569599134], 'cr_4l_sf': [11.429296468617395, 0.5261506111547267], 'sr_4l_sf_A': [0.58405411161948, 0.034532871510194546], 'sr_4l_of_1': [0.2247555028880015, 0.03622101953838245], 'sr_4l_of_2': [0.0004370792303234339, 0.010489449262681468], 'cr_4l_btag_of': [2.6049697652924806, 0.07841745866607606], 'sr_sf_all_cutbased': [1.5218016663566232, 0.19034565036587744], 'sr_of_all_cutbased': [0.8338602496078238, 0.09277997123333537], 'sr_all_cutbased': [2.355661915964447, 0.2831256215992128]}, '$S/\\sqrt{B}$': {'sr_4l_sf_A': [1.6997341741567213, None], 'sr_4l_sf_B': [1.2928955851584796, None], 'sr_4l_sf_C': [0.6690634727354225, None], 'sr_4l_of_1': [2.704040585334127, None], 'sr_4l_of_2': [1.8768995488452491, None], 'sr_4l_of_3': [1.4570760328884798, None], 'sr_4l_of_4': [2.338104016392618, None], 'sr_sf_all_cutbased': [2.2379278780754537, None], 'sr_of_all_cutbased': [4.2923639594842085, None], 'sr_all_cutbased': [4.840734401735572, None]}, '$S/\\sqrt{S+B}$': {'sr_4l_sf_A': [1.22560967692322, None], 'sr_4l_sf_B': [1.0551906003922045, None], 'sr_4l_sf_C': [0.5697580749248254, None], 'sr_4l_of_1': [1.5363956729406407, None], 'sr_4l_of_2': [1.1255886886700301, None], 'sr_4l_of_3': [0.9812913654271258, None], 'sr_4l_of_4': [1.6257598316487298, None], 'sr_sf_all_cutbased': [1.714692551819697, None], 'sr_of_all_cutbased': [2.6895147023128505, None], 'sr_all_cutbased': [3.18961751331142, None]}, 'Sig': {'sr_4l_sf_A': [3.128938727430068, 0.00013680466095312944], 'sr_4l_sf_B': [3.3345436750423687, 0.0001506814168960968], 'sr_4l_sf_C': [1.181228760097838, 4.484353004512215e-05], 'sr_4l_of_1': [3.485870358581451, 0.00015704724282298693], 'sr_4l_of_2': [1.9785203373594413, 8.326724265232396e-05], 'sr_4l_of_3': [1.7621821535549316, 6.173860641957155e-05], 'sr_4l_of_4': [5.117193291793228, 0.00016520122725379688], 'sr_sf_all_cutbased': [7.644711162570275, 0.0003323296078943484], 'sr_of_all_cutbased': [12.343766141289052, 0.00046725431914867934], 'sr_all_cutbased': [19.988477303859327, 0.0007995839270430276]}, 'Bkg': {'sr_4l_sf_A': [3.3886920578188438, 0.04534706029917568], 'sr_4l_sf_B': [6.651903116671747, 0.14940524226321245], 'sr_4l_sf_C': [3.116975467579323, 0.037097359471379925], 'sr_4l_of_1': [1.6618661864340538, 0.0650148965509467], 'sr_4l_of_2': [1.1112172580360493, 0.015732913853981147], 'sr_4l_of_3': [1.4626390627781802, 0.03637136447378289], 'sr_4l_of_4': [4.790005233553529, 0.08110258711349308], 'sr_sf_all_cutbased': [13.157570642069913, 0.23184966203376806], 'sr_of_all_cutbased': [9.025727740801813, 0.19822176199220382], 'sr_all_cutbased': [22.183298382871726, 0.43007142402597176]}, 'Zmetric': {'sr_4l_sf_A': [1.506006574717711, None], 'sr_4l_sf_B': [1.2026702524278647, None], 'sr_4l_sf_C': [0.6323467198919773, None], 'sr_4l_of_1': [2.1606677653842192, None], 'sr_4l_of_2': [1.536970006909602, None], 'sr_4l_of_3': [1.254969951949225, None], 'sr_4l_of_4': [2.04091460317212, None], 'sr_sf_all_cutbased': [2.028382092537306, None], 'sr_of_all_cutbased': [3.5736597476881857, None], 'sr_all_cutbased': [4.109182145582356, None]}}

# Ref with all current SFs applied (btag, lep, PU, prefire) Mar 02, 2024
# feb27_baseline_withSFnoSys_withCBpresel.pkl.gz
EWK_REF = {'WWZ': {'sr_4l_sf_A': [1.8796742427699058, 5.0011966352230034e-05], 'sr_4l_sf_B': [1.596064584710653, 4.1913361034359516e-05], 'sr_4l_sf_C': [0.46329822942371474, 1.2478240476784822e-05], 'sr_4l_of_1': [0.5089704111430158, 1.2928098662250936e-05], 'sr_4l_of_2': [0.5863582014957153, 1.5105571655629712e-05], 'sr_4l_of_3': [1.1904609448245147, 3.094890749661624e-05], 'sr_4l_of_4': [4.221321537855876, 0.00011366124456835628], 'all_events': [150.97044577809945, 0.005508468309953707], '4l_presel': [22.942228933056704, 0.0006758208952004465], 'sr_4l_sf_incl': [4.086939376681597, 0.0001082900240082549], 'sr_4l_of_incl': [7.524496275892502, 0.0001990142032474757], 'cr_4l_btag_of': [2.6844044528430437, 0.00010151390499609684], 'cr_4l_btag_sf_offZ_met80': [1.352840130589533, 5.099584303344667e-05], 'cr_4l_sf': [0.8578601477161325, 2.251552161108323e-05], 'sr_sf_all_cutbased': [3.939037056904273, 0.00010440356786337436], 'sr_of_all_cutbased': [6.507111095319122, 0.00017264382238285316], 'sr_all_cutbased': [10.446148152223394, 0.0002770473902462275]}, 'ZH': {'sr_4l_sf_A': [0.729182685851617, 4.514358815906539e-05], 'sr_4l_sf_B': [1.1802983133547118, 6.439208199854825e-05], 'sr_4l_sf_C': [0.5238591034915036, 1.9353931412792887e-05], 'sr_4l_of_1': [2.3741376736563735, 9.673793250202969e-05], 'sr_4l_of_2': [1.0512964047749698, 4.274738875297518e-05], 'sr_4l_of_3': [0.27487046285644734, 1.2516520867463701e-05], 'sr_4l_of_4': [0.11932959533821874, 7.591287106296967e-06], 'all_events': [145.62539717282453, 0.01004702745166159], '4l_presel': [17.639006952903085, 0.000897724478318867], 'sr_4l_sf_incl': [2.585684086912913, 0.00013520743152816494], 'sr_4l_of_incl': [6.207629455186265, 0.0002664887798313454], 'cr_4l_btag_of': [1.8996142891682295, 0.00014267977972230728], 'cr_4l_btag_sf_offZ_met80': [0.686375081593778, 6.94158621939911e-05], 'cr_4l_sf': [0.055719411432755335, 2.139009834632954e-06], 'sr_sf_all_cutbased': [2.4333401026978323, 0.00012888960157040653], 'sr_of_all_cutbased': [3.8196341366260094, 0.00015959312922876555], 'sr_all_cutbased': [6.2529742393238426, 0.00028848273079917205]}, 'ZZ': {'sr_4l_sf_A': [0.9897744828273218, 0.0004867788096212909], 'sr_4l_sf_B': [3.3126167722389104, 0.0015901957641724686], 'sr_4l_sf_C': [1.971583247515798, 0.0009209113661524344], 'sr_4l_of_1': [0.4613732783213163, 0.00022237833938276952], 'sr_4l_of_2': [0.46108370261056303, 0.00022786268728842662], 'sr_4l_of_3': [0.2935959769568094, 0.00014762270992032994], 'sr_4l_of_4': [0.3891131801194942, 0.0002007563039401573], 'all_events': [22132.43491159535, 15.642873844819166], '4l_presel': [2691.1766010171436, 1.4391414532425766], 'sr_4l_sf_incl': [9.809160716833057, 0.00478333333872781], 'sr_4l_of_incl': [15.519336306919715, 0.007525034060702233], 'cr_4l_btag_of': [4.1592982373747365, 0.002939994893574839], 'cr_4l_btag_sf_offZ_met80': [2.7554293416441475, 0.0020865159266828443], 'cr_4l_sf': [1431.1449772069113, 0.698969642521474], 'sr_sf_all_cutbased': [6.27397450258203, 0.002997885939946194], 'sr_of_all_cutbased': [1.605166138008183, 0.0007986200405316834], 'sr_all_cutbased': [7.879140640590213, 0.003796505980477877]}, 'ttZ': {'sr_4l_sf_A': [0.9744323818401107, 0.002988213693670242], 'sr_4l_sf_B': [0.8260914480482815, 0.0024411443158970365], 'sr_4l_sf_C': [0.2117777372878653, 0.0006793733688719802], 'sr_4l_of_1': [0.2274315031538946, 0.0006141271088714311], 'sr_4l_of_2': [0.24264839413990152, 0.0006803625592641666], 'sr_4l_of_3': [0.6497391620795682, 0.0017203709264879009], 'sr_4l_of_4': [1.8370650087275844, 0.005608680086156667], 'all_events': [6584.347460577734, 36.5765832029854], '4l_presel': [139.51959476092196, 0.46530609600697137], 'sr_4l_sf_incl': [2.1481188650717815, 0.006441869102479813], 'sr_4l_of_incl': [3.3987615017713977, 0.009902575175096328], 'cr_4l_btag_of': [57.02843920704273, 0.19120787997139224], 'cr_4l_btag_sf_offZ_met80': [27.96411533276142, 0.09254164193739435], 'cr_4l_sf': [0.42989697675843075, 0.001306420936671004], 'sr_sf_all_cutbased': [2.012301567176258, 0.006108731378439258], 'sr_of_all_cutbased': [2.9568840681009485, 0.008623540680780166], 'sr_all_cutbased': [4.969185635277206, 0.014732272059219422]}, 'tWZ': {'sr_4l_sf_A': [0.40417646308070293, 8.703750596671525e-05], 'sr_4l_sf_B': [0.37410199380609116, 8.532189290243701e-05], 'sr_4l_sf_C': [0.0957249409079119, 2.258685361984745e-05], 'sr_4l_of_1': [0.11518523443695693, 2.4047369970902342e-05], 'sr_4l_of_2': [0.1296259288057434, 2.883456646287243e-05], 'sr_4l_of_3': [0.26653774477742836, 5.913066647446086e-05], 'sr_4l_of_4': [0.8301279967178373, 0.00019654830703977157], 'all_events': [192.49430732626087, 0.05881629917636074], '4l_presel': [23.97344544819388, 0.005715484631594204], 'sr_4l_sf_incl': [0.9235819068409694, 0.00020599139869891832], 'sr_4l_of_incl': [1.5294180689195565, 0.0003516027067161033], 'cr_4l_btag_of': [9.237700282558725, 0.0022026516549542774], 'cr_4l_btag_sf_offZ_met80': [4.255080503872744, 0.0010152383797496227], 'cr_4l_sf': [0.17782139062034538, 4.003134982851281e-05], 'sr_sf_all_cutbased': [0.874003397794706, 0.0001949462524889997], 'sr_of_all_cutbased': [1.341476904737966, 0.0003085609099480072], 'sr_all_cutbased': [2.215480302532672, 0.0005035071624370069]}, 'WZ': {'sr_4l_sf_A': [0.07096471161892892, 0.0031669096309620373], 'sr_4l_sf_B': [0.38406512417674726, 0.012442842468890386], 'sr_4l_sf_C': [0.16260528877019054, 0.004210805556691592], 'sr_4l_of_1': [0.3724353881761476, 0.0178115542873584], 'sr_4l_of_2': [0.10258944561447007, 0.0023701419327491835], 'sr_4l_of_3': [0.011818492291620748, 0.010013353958452284], 'sr_4l_of_4': [0.6973773005085421, 0.03337936876254605], 'all_events': [67206.13135117064, 4553.624485563815], '4l_presel': [8.859950688690056, 0.4746705553140165], 'sr_4l_sf_incl': [0.7129665173515559, 0.02258580279217667], 'sr_4l_of_incl': [1.6520052483661258, 0.0804926867847766], 'cr_4l_btag_of': [0.8860196415778555, 0.059575910521280286], 'cr_4l_btag_sf_offZ_met80': [0.11606954799352177, 0.010279649434932512], 'cr_4l_sf': [0.14670822278034434, 0.009247351563810902], 'sr_sf_all_cutbased': [0.6176351245658667, 0.019820557656544016], 'sr_of_all_cutbased': [1.1842206265907804, 0.06357441894110591], 'sr_all_cutbased': [1.8018557511566469, 0.08339497659764994]}, 'other': {'sr_4l_sf_A': [0.4993868268217272, 0.027154017210982366], 'sr_4l_sf_B': [0.6033849836294176, 0.07499426824703304], 'sr_4l_sf_C': [0.09591669267138649, 0.02215997619394746], 'sr_4l_of_1': [0.1393089487979516, 0.019281127575627503], 'sr_4l_of_2': [-0.006085849094745107, 0.007863615554409087], 'sr_4l_of_3': [0.05654106816261955, 0.011564137349614385], 'sr_4l_of_4': [0.4417725773375794, 0.019546908278451205], 'all_events': [2747882.281747165, 5947477.843846047], '4l_presel': [37.71068049919267, 3.5641995992582793], 'sr_4l_sf_incl': [1.2401878850076171, 0.12465785077497647], 'sr_4l_of_incl': [1.4477455619548447, 0.12790510023309565], 'cr_4l_btag_of': [2.494962082595951, 0.06885015229405128], 'cr_4l_btag_sf_offZ_met80': [1.5760296544906491, 0.052401971126985615], 'cr_4l_sf': [9.464967629072424, 0.37356337925725397], 'sr_sf_all_cutbased': [1.1986885031225314, 0.12430826165196288], 'sr_of_all_cutbased': [0.6315367452034055, 0.05825578875810218], 'sr_all_cutbased': [1.8302252483259367, 0.18256405041006507]}, '$S/\\sqrt{B}$': {'sr_4l_sf_A': [1.5218437321425393, None], 'sr_4l_sf_B': [1.1838171020186103, None], 'sr_4l_sf_C': [0.6196894699730852, None], 'sr_4l_of_1': [2.5134880048170003, None], 'sr_4l_of_2': [1.6982950716613776, None], 'sr_4l_of_3': [1.2960774065043659, None], 'sr_4l_of_4': [2.119167541718631, None], 'sr_sf_all_cutbased': [2.025202783745074, None], 'sr_of_all_cutbased': [3.9207800007540765, None], 'sr_all_cutbased': [4.412931240071834, None]}, '$S/\\sqrt{S+B}$': {'sr_4l_sf_A': [1.1076384270013098, None], 'sr_4l_sf_B': [0.9650496099402631, None], 'sr_4l_sf_C': [0.5258008422662838, None], 'sr_4l_of_1': [1.4070066010959612, None], 'sr_4l_of_2': [1.022034867751365, None], 'sr_4l_of_3': [0.8846639338623044, None], 'sr_4l_of_4': [1.4856788881128868, None], 'sr_sf_all_cutbased': [1.5603365535498177, None], 'sr_of_all_cutbased': [2.452365160592572, None], 'sr_all_cutbased': [2.906672503264163, None]}, 'Sig': {'sr_4l_sf_A': [2.608856928621523, 9.515555451129542e-05], 'sr_4l_sf_B': [2.776362898065365, 0.00010630544303290777], 'sr_4l_sf_C': [0.9871573329152183, 3.183217188957771e-05], 'sr_4l_of_1': [2.8831080847993893, 0.00010966603116428062], 'sr_4l_of_2': [1.6376546062706852, 5.785296040860489e-05], 'sr_4l_of_3': [1.465331407680962, 4.3465428364079944e-05], 'sr_4l_of_4': [4.340651133194095, 0.00012125253167465325], 'sr_sf_all_cutbased': [6.372377159602106, 0.00023329316943378092], 'sr_of_all_cutbased': [10.326745231945132, 0.0003322369516116187], 'sr_all_cutbased': [16.699122391547235, 0.0005655301210453997]}, 'Bkg': {'sr_4l_sf_A': [2.938734866188791, 0.03388295685120265], 'sr_4l_sf_B': [5.500260321899448, 0.09155377268889536], 'sr_4l_sf_C': [2.5376079071531517, 0.027993653339283314], 'sr_4l_of_1': [1.315734352886267, 0.03795323468121101], 'sr_4l_of_2': [0.9298616220759329, 0.011170817300173736], 'sr_4l_of_3': [1.2782324442680462, 0.02350461561094936], 'sr_4l_of_4': [4.195456063411037, 0.058932261738133845], 'sr_sf_all_cutbased': [10.97660309524139, 0.15343038287938132], 'sr_of_all_cutbased': [7.719284482641283, 0.13156092933046795], 'sr_all_cutbased': [18.695887577882672, 0.2849913122098493]}, 'Zmetric': {'sr_4l_sf_A': [1.3535126937574338, None], 'sr_4l_sf_B': [1.1007210618985692, None], 'sr_4l_sf_C': [0.5848969322063289, None], 'sr_4l_of_1': [1.9946363806225782, None], 'sr_4l_of_2': [1.3928855127238835, None], 'sr_4l_of_3': [1.122627723636731, None], 'sr_4l_of_4': [1.8561096981857161, None], 'sr_sf_all_cutbased': [1.8400238828814135, None], 'sr_of_all_cutbased': [3.2594693678628235, None], 'sr_all_cutbased': [3.7429705381701934, None]}}

SOVERROOTB = "$S/\sqrt{B}$"
SOVERROOTSPLUSB = "$S/\sqrt{S+B}$"



################### Getting and printing yields ###################

# Get the yields in the SR
def get_yields(histos_dict,sample_dict,raw_counts=False,quiet=True,blind=True,systematic_name="nominal"):

    yld_dict = {}

    # Look at the yields in one histo (e.g. njets)
    if raw_counts: dense_axis = "njets_counts"
    else: dense_axis = "njets"
    for proc_name in sample_dict.keys():
        yld_dict[proc_name] = {}
        for cat_name in histos_dict[dense_axis].axes["category"]:
            #if "bdt" in cat_name: continue # TMP!!!
            if blind and (("data" in proc_name) and (not cat_name.startswith("cr_"))):
                # If this is data and we're not in a CR category, put placeholder numbers for now
                yld_dict[proc_name][cat_name] = [-999,-999]
            else:
                val = sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name],"systematic":systematic_name}].values(flow=True)))
                var = sum(sum(histos_dict[dense_axis][{"category":cat_name,"process":sample_dict[proc_name],"systematic":systematic_name}].variances(flow=True)))
                yld_dict[proc_name][cat_name] = [val,var]

    # Print to screen
    if not quiet:
        for proc in yld_dict.keys():
            print(f"\n{proc}:")
            for cat in yld_dict[proc].keys():
                val = yld_dict[proc][cat]
                print(f"\t{cat}: {val}")

    return yld_dict


# Gets the process sums for S and B and gets metrics e.g. S/sqrt(B) and puts it into the dict
# Hard coded for the summed values (e.g. looking for "ZH" not "GluGluZH","qqToZHToZTo2L")
def put_proc_row_sums(yld_dict,sr_cat_lst, sig_lst=sg.SIG_LST,bkg_lst=sg.BKG_LST):
    # Build up the empty dicts for sig and bkg that we can later fill and then put into the yld dict
    # Will look something like this: {"sr_4l_sf_A":[0,0], "sr_4l_sf_B":[0,0], "sr_4l_sf_C":[0,0], "sr_4l_of_1":[0,0], "sr_4l_of_2":[0,0], "sr_4l_of_3":[0,0], "sr_4l_of_4":[0,0]}
    sig_sum = {}
    bkg_sum = {}
    for sr_cat_name in sr_cat_lst:
        sig_sum[sr_cat_name] = [0,0]
        bkg_sum[sr_cat_name] = [0,0]

    # Finding sums
    for proc in yld_dict.keys():
        print(proc)
        for cat in yld_dict[proc].keys():
            if cat not in sig_sum: continue
            val,var = yld_dict[proc][cat]
            print("   ",cat,val)
            if proc in sig_lst:
                sig_sum[cat][0] += val
                sig_sum[cat][1] += var
            if proc in bkg_lst:
                bkg_sum[cat][0] += val
                bkg_sum[cat][1] += var

    # Finding metrics, and putting sums and metrics into the yld dict
    if SOVERROOTB not in yld_dict:      yld_dict[SOVERROOTB] = {}
    if SOVERROOTSPLUSB not in yld_dict: yld_dict[SOVERROOTSPLUSB] = {}
    if "Sig" not in yld_dict: yld_dict["Sig"] = {}
    if "Bkg" not in yld_dict: yld_dict["Bkg"] = {}
    if "Zmetric" not in yld_dict: yld_dict["Zmetric"] = {}
    for cat in sig_sum.keys():
        s = sig_sum[cat][0]
        b = bkg_sum[cat][0]
        s_var = sig_sum[cat][1]
        b_var = bkg_sum[cat][1]
        yld_dict[SOVERROOTB][cat]      = [s/math.sqrt(b) , None]
        yld_dict[SOVERROOTSPLUSB][cat] = [s/math.sqrt(s+b) , None]
        yld_dict["Zmetric"][cat] = [math.sqrt(2 * ((s + b) * math.log(1 + s / b) - s)), None] # Eq 18 https://cds.cern.ch/record/2203244/files/1087459_109-114.pdf
        yld_dict["Sig"][cat] = [s, s_var]
        yld_dict["Bkg"][cat] = [b, b_var]


# Gets the sums of categoreis (assumed to be columns in the input dict) and puts them into the dict
# Special handling for rows that are metrics (e.g. s/sqrt(b)), sums these in quadrature
def put_cat_col_sums(yld_dict,sr_sf_lst,sr_of_lst,metrics_names_lst=["Zmetric",SOVERROOTB,SOVERROOTSPLUSB],tag=""):

    # The full SR list should be the sf and of together
    sr_lst = sr_sf_lst + sr_of_lst

    # Loop over rows (processes) and sum columns together, fill the result into new_dict
    new_dict = {}
    for proc in yld_dict:
        sr_sf_val = 0
        sr_sf_var = 0
        sr_of_val = 0
        sr_of_var = 0
        sr_val = 0
        sr_var = 0
        for cat in yld_dict[proc]:
            val = yld_dict[proc][cat][0]
            var = yld_dict[proc][cat][1]
            if cat in sr_sf_lst:
                if proc in metrics_names_lst:
                    sr_sf_val += val*val
                else:
                    sr_sf_val += val
                    sr_sf_var += var
            if cat in sr_of_lst:
                if proc in metrics_names_lst:
                    sr_of_val += val*val
                else:
                    sr_of_val += val
                    sr_of_var += var
            if cat in sr_lst:
                if proc in metrics_names_lst:
                    sr_val += val*val
                else:
                    sr_val += val
                    sr_var += var

        # Fill our new_dict with what we've computed
        new_dict[proc] = {}
        if proc in metrics_names_lst:
            new_dict[proc]["sr_sf_all"] = [np.sqrt(sr_sf_val),None]
            new_dict[proc]["sr_of_all"] = [np.sqrt(sr_of_val),None]
            new_dict[proc]["sr_all"]    = [np.sqrt(sr_val),None]
        else:
            new_dict[proc]["sr_sf_all"] = [sr_sf_val, sr_sf_var]
            new_dict[proc]["sr_of_all"] = [sr_of_val, sr_of_var]
            new_dict[proc]["sr_all"]    = [sr_val, sr_var]

    # Put the columns into the yld_dict
    for proc in new_dict:
        yld_dict[proc][f"sr_sf_all{tag}"] = new_dict[proc]["sr_sf_all"]
        yld_dict[proc][f"sr_of_all{tag}"] = new_dict[proc]["sr_of_all"]
        yld_dict[proc][f"sr_all{tag}"] = new_dict[proc]["sr_all"]


# Print yields
def print_yields(yld_dict_in,cats_to_print,procs_to_print,print_fom=True,hlines=[]):

    # Get err from var
    def get_err_from_var(in_dict):
        out_dict = {}
        for proc in in_dict:
            out_dict[proc] = {}
            for cat in in_dict[proc]:
                if in_dict[proc][cat][1] is None: var = None
                else: var = np.sqrt(in_dict[proc][cat][1])
                out_dict[proc][cat] = [in_dict[proc][cat][0],var]
        return out_dict

    yld_dict = get_err_from_var(yld_dict_in)

    # Print the yields directly
    mlt.print_latex_yield_table(
        yld_dict,
        tag="All yields",
        key_order=SAMPLE_DICT_BASE.keys(),
        subkey_order=cats_to_print,
        print_begin_info=True,
        print_end_info=True,
        print_errs=True,
        column_variable="subkeys",
        size="tiny",
        hz_line_lst=[6],
    )
    #exit()

    ### Compare with other yields, print comparison ###

    tag1 = "New"
    tag2 = "Ref"

    yld_dict_comp = get_err_from_var(EWK_REF)

    yld_dict_1 = copy.deepcopy(yld_dict)
    yld_dict_2 = copy.deepcopy(yld_dict_comp)

    pdiff_dict = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="percent_diff",inpercent=True)
    diff_dict  = utils.get_diff_between_nested_dicts(yld_dict_1,yld_dict_2,difftype="absolute_diff")

    mlt.print_begin()
    mlt.print_latex_yield_table(yld_dict_1,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag1,hz_line_lst=hlines,print_errs=True,size="tiny",column_variable="keys")
    mlt.print_latex_yield_table(yld_dict_2,key_order=procs_to_print,subkey_order=cats_to_print,tag=tag2,hz_line_lst=hlines,print_errs=True,size="tiny",column_variable="keys")
    mlt.print_latex_yield_table(pdiff_dict,key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Percent diff between {tag1} and {tag2}",hz_line_lst=hlines,size="tiny",column_variable="keys")
    mlt.print_latex_yield_table(diff_dict, key_order=procs_to_print,subkey_order=cats_to_print,tag=f"Diff between {tag1} and {tag2}",hz_line_lst=hlines,size="tiny",column_variable="keys")
    mlt.print_end()


# Dump the counts dict to a latex table
def print_counts(counts_dict):

    cats_to_print = ["all_events", "4l_presel", "sr_4l_sf_A", "sr_4l_sf_B", "sr_4l_sf_C", "sr_4l_of_1", "sr_4l_of_2", "sr_4l_of_3", "sr_4l_of_4"]

    # Print the yields directly
    mlt.print_latex_yield_table(
        counts_dict,
        tag="Raw MC counts (ewkcoffea)",
        key_order=counts_dict.keys(),
        subkey_order=cats_to_print,
        print_begin_info=True,
        print_end_info=True,
        column_variable="subkeys",
    )


# This should maybe be in a different script
################### Hist manipulation and plotting ###################


# Get the list of categories on the sparese axis
def get_axis_cats(histo,axis_name):
    process_list = [x for x in histo.axes[axis_name]]
    return process_list


# Merges the last bin (overflow) into the second to last bin, zeros the content of the last bin, returns a new hist
# Note assumes just one axis!
def merge_overflow(hin):
    hout = copy.deepcopy(hin)
    for cat_idx,arr in enumerate(hout.values(flow=True)):
        hout.values(flow=True)[cat_idx][-2] += hout.values(flow=True)[cat_idx][-1]
        hout.values(flow=True)[cat_idx][-1] = 0
        hout.variances(flow=True)[cat_idx][-2] += hout.variances(flow=True)[cat_idx][-1]
        hout.variances(flow=True)[cat_idx][-1] = 0
    return hout


# Rebin according to https://github.com/CoffeaTeam/coffea/discussions/705
def rebin(histo,factor):
    return histo[..., ::hist.rebin(factor)]


# Regroup categories (e.g. processes)
def group(h, oldname, newname, grouping):

    # Build up a grouping dict that drops any proc that is not in our h
    grouping_slim = {}
    proc_lst = get_axis_cats(h,oldname)
    for grouping_name in grouping.keys():
        for proc in grouping[grouping_name]:
            if proc in proc_lst:
                if grouping_name not in grouping_slim:
                    grouping_slim[grouping_name] = []
                grouping_slim[grouping_name].append(proc)
            #else:
            #    print(f"WARNING: process {proc} not in this hist")

    # From Nick: https://github.com/CoffeaTeam/coffea/discussions/705#discussioncomment-4604211
    hnew = hist.Hist(
        hist.axis.StrCategory(grouping_slim, name=newname),
        *(ax for ax in h.axes if ax.name != oldname),
        storage=h.storage_type(),
    )
    for i, indices in enumerate(grouping_slim.values()):
        hnew.view(flow=True)[i] = h[{oldname: indices}][{oldname: sum}].view(flow=True)

    return hnew


# Takes a mc hist and data hist and plots both
def make_cr_fig(histo_mc,histo_data=None,title="test",unit_norm_bool=False):

    # Create the figure
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(7,7),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.07)

    # Plot the mc
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        color=CLR_LST,
        ax=ax,
    )
    # Plot the data
    if histo_data is not None:
        histo_data.plot1d(
            stack=False,
            histtype="errorbar",
            color="k",
            ax=ax,
            w2=histo_data.variances(),
            w2method="sqrt",
        )
    # Plot a dummy hist on rax to get the label to show up
    histo_mc.plot1d(alpha=0, ax=rax)

    ### Get the errs on MC and plot them by hand ###
    histo_mc_sum = histo_mc[{"process_grp":sum}]
    mc_arr = histo_mc_sum.values()
    mc_err_arr = np.sqrt(histo_mc_sum.variances())
    err_p = np.append(mc_arr + mc_err_arr, 0)
    err_m = np.append(mc_arr - mc_err_arr, 0)
    bin_edges_arr = histo_mc_sum.axes[0].edges
    bin_centers_arr = histo_mc_sum.axes[0].centers
    ax.fill_between(bin_edges_arr,err_m,err_p, step='post', facecolor='none', edgecolor='gray', alpha=0.5, linewidth=0.0, label='MC stat', hatch='/////')

    ### Get the errs on data and ratios and plot them by hand ###
    if histo_data is not None:
        histo_data_sum = histo_data[{"process_grp":sum}]

        data_arr = histo_data_sum.values()
        data_err_arr = np.sqrt(histo_data_sum.variances())

        err_ratio_p = np.append(1+mc_err_arr/mc_arr,1)
        err_ratio_m = np.append(1-mc_err_arr/mc_arr,1)

        data_ratio_err_p = (data_arr + data_err_arr)/mc_arr
        data_ratio_err_m = (data_arr - data_err_arr)/mc_arr

        rax.fill_between(bin_edges_arr,err_ratio_m,err_ratio_p,step='post', facecolor='none',edgecolor='gray', label='MC stat', linewidth=0.0, hatch='/////',alpha=0.5)
        rax.scatter(bin_centers_arr,data_arr/mc_arr,facecolor='black',edgecolor='black',marker="o")
        rax.vlines(bin_centers_arr,data_ratio_err_p,data_ratio_err_m,color='k')

    # Scale the y axis and labels
    ax.legend(fontsize="12")
    ax.set_title(title)
    ax.autoscale(axis='y')
    ax.set_xlabel(None)
    rax.set_ylabel('Ratio')
    rax.set_ylim(0.0,2.0)
    rax.axhline(1.0,linestyle="-",color="k",linewidth=1)
    ax.tick_params(axis='y', labelsize=16)
    rax.tick_params(axis='x', labelsize=16)
    #ax.set_yscale('log')

    return fig


# Plots a hist
def make_single_fig(histo_mc,ax_to_overlay="process",err_p=None,err_m=None,ylims=None,unit_norm_bool=False,title=None):
    #print("\nPlotting values:",histo.values())
    fig, ax = plt.subplots(1, 1, figsize=(12,7))

    # Plot the mc
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        color=CLR_LST,
        yerr=True,
        overlay=ax_to_overlay
    )

    # Set title and y lims
    if title is not None: plt.title(title)
    if ylims is None:
        ax.autoscale(axis='y')
    else:
        ax.set_ylim(ylims)

    # Draw errors
    if (err_p is not None) and (err_m is not None):
        bin_edges_arr = histo_mc.axes[0].edges
        bin_centers_arr = histo_mc.axes[0].centers
        ax.fill_between(bin_edges_arr,err_m,err_p, step='post', facecolor='none', edgecolor='gray', alpha=0.5, linewidth=0.0, label='MC stat', hatch='/////')

    plt.legend()
    return fig


# Takes a mc hist and data hist and plots both
def make_syst_fig(histo_mc,mc_up_arr,mc_do_arr,syst,histo_data=None,title="test",unit_norm_bool=False):

    # Create the figure
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(7,7),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.07)

    # Plot the mc
    histo_mc.plot1d(
        stack=True,
        histtype="fill",
        color=CLR_LST,
        ax=ax,
    )

    # Plot the syst
    histo_mc_sum = histo_mc[{"process_grp":sum}]
    bin_edges_arr = histo_mc_sum.axes[0].edges
    bin_centers_arr = histo_mc_sum.axes[0].centers
    ax.stairs(mc_up_arr, bin_edges_arr, color="cyan", linestyle="--", label=f'{syst} up')
    ax.stairs(mc_do_arr, bin_edges_arr, color="magenta", linestyle="--", label=f'{syst} down')

    # Plot the syst on ratio plots
    mc_arr = histo_mc_sum.values()
    rax.scatter(bin_centers_arr,mc_up_arr/mc_arr,facecolor='cyan',edgecolor='cyan',marker="o")
    rax.scatter(bin_centers_arr,mc_do_arr/mc_arr,facecolor='magenta',edgecolor='magenta',marker="o")

    # Scale the y axis and labels
    ax.legend(fontsize="12")
    ax.set_title(title)
    ax.autoscale(axis='y')
    ax.set_xlabel(None)
    rax.set_ylabel('Ratio')
    rax.set_ylim(0.9,1.1)
    rax.axhline(1.0,linestyle="-",color="k",linewidth=1)
    ax.tick_params(axis='y', labelsize=16)
    rax.tick_params(axis='x', labelsize=16)

    return fig


# IN PROGRESS
# Main function for checking individual systematics
def make_syst_plots(histo_dict,grouping_mc,grouping_data,save_dir_path,year):

    for var_name in histo_dict.keys():
        #print(f"\n{var_name}")
        if var_name not in TMP_VAR_LST: continue
        histo = histo_dict[var_name]

        cat_lst = [
            "sr_4l_sf_A",
            "sr_4l_sf_B",
            "sr_4l_sf_C",
            "sr_4l_of_1",
            "sr_4l_of_2",
            "sr_4l_of_3",
            "sr_4l_of_4",
            #"cr_4l_sf",
            #"cr_4l_btag_sf_offZ_met80",
            #"cr_4l_btag_of",
            #"sr_4l_of_incl",
            #"sr_4l_sf_incl",
        ]

        # Rebin if continous variable
        if var_name not in ["njets","nbtagsl","nleps"]:
            histo = rebin(histo,6)

        # Get the list of systematic base names (i.e. without the up and down tags)
        # Assumes each syst has a "systnameUp" and a "systnameDown" category on the systematic axis
        syst_var_lst = []
        all_syst_var_lst = histo.axes["systematic"]
        for syst_var_name in all_syst_var_lst:
            if syst_var_name.endswith("Up"):
                syst_name_base = syst_var_name.replace("Up","")
                if syst_name_base not in syst_var_lst:
                    syst_var_lst.append(syst_name_base)

        for cat in cat_lst:
            print("\n",cat)
            if "cr_4l_of" not in cat and var_name == "j0pt": continue
            histo_cat = histo[{"category":cat}]
            histo_grouped_mc = group(histo_cat,"process","process_grp",grouping_mc)
            histo_grouped_data = group(histo_cat,"process","process_grp",grouping_data)

            mc_nom   = merge_overflow(histo_grouped_mc[{"systematic":"nominal"}])
            data_nom = merge_overflow(histo_grouped_data[{"systematic":"nominal"}])

            for syst in syst_var_lst:
                #if "btag" not in syst: continue
                #if "uncorrelated" not in syst: continue
                #if "lepSF" not in syst: continue
                #if "PreFiring" not in syst: continue
                #if "PU" not in syst: continue
                #if "ISR" not in syst and "FSR" not in syst: continue
                if "renorm" not in syst and "fact" not in syst: continue

                # Skip the variations that don't apply (TODO: why are these in the hist to begin with??)
                if year == "UL16APV": blacklist_years = ["2016","2017","2018"]
                if year == "UL16": blacklist_years = ["2016APV","2017","2018"]
                if year == "UL17": blacklist_years = ["2016APV","2016","2018"]
                if year == "UL18": blacklist_years = ["2016APV","2016","2017"]
                if year == "all": blacklist_years = []
                skip = False
                for y in blacklist_years:
                    if syst.endswith(y):
                        skip = True
                if skip: continue

                if year == "all": yeartag = "FullR2"
                else: yeartag = year

                mc_up     = merge_overflow(histo_grouped_mc[{"systematic":f"{syst}Up"}])
                mc_down   = merge_overflow(histo_grouped_mc[{"systematic":f"{syst}Down"}])
                data_up   = merge_overflow(histo_grouped_data[{"systematic":f"{syst}Up"}])
                data_down = merge_overflow(histo_grouped_data[{"systematic":f"{syst}Down"}])

                mc_up_arr = mc_up[{"process_grp":sum}].values()
                mc_down_arr = mc_down[{"process_grp":sum}].values()

                # Print individual syst numbers
                #if var_name != "nleps": continue
                #n = sum(sum(mc_nom.values()))
                #u = sum(mc_up_arr)
                #d = sum(mc_down_arr)
                #print("\n",syst)
                #print("nom",n)
                #print("up",u)
                #print("do",d)
                #r_up = abs((n-u)/n)
                #r_do = abs((n-d)/n)
                #r = (r_up+r_do)/2
                #print("err",np.round(100*abs(n-u)/n,1),"%")
                #print("err up",np.round(100*r_up,1),"%")
                #print("err do",np.round(100*r_do,1),"%")
                #print("err do",np.round(100*r,1),"%")
                #continue

                fig = make_syst_fig(mc_nom,mc_up_arr,mc_down_arr,syst,title=f"{var_name}_{yeartag}_{cat}_{syst}")

                out_path_for_this_cat = os.path.join(save_dir_path,os.path.join(yeartag,cat))
                if not os.path.exists(out_path_for_this_cat): os.makedirs(out_path_for_this_cat)
                fig.savefig(f"{out_path_for_this_cat}/{var_name}_{yeartag}_{cat}_{syst}.png")

            make_html(os.path.join(os.getcwd(),out_path_for_this_cat))


# A function for making a summary plot of SR yields
def make_sr_comb_plot(histo_dict,grouping_mc,grouping_data):

    # Declare the hist we'll be filling
    sr_lst  = ["sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C" , "sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4"]
    proc_lst  = ["WWZ", "ZH", "ZZ", "ttZ", "tWZ", "WZ", "other"]
    histo_comb = hist.Hist(
        hist.axis.StrCategory(proc_lst, name="process", label="process"),
        hist.axis.StrCategory(sr_lst,   name="cat",     label="Cut-based SRs"),
    )

    # Get the yield dict
    year = "all"
    histo = histo_dict["njets"]
    sample_names_dict_mc   = sg.create_mc_sample_dict(sg.SAMPLE_DICT_BASE,"all")
    sample_names_dict_data = sg.create_data_sample_dict(year)
    yld_dict_mc   = yt.get_yields(histo,sample_names_dict_mc)
    yld_dict_data = yt.get_yields(histo,sample_names_dict_data)

    # Apply the data-driven normalization for ZZ, ttZ
    yld_dict_mc, _, _ = yt.do_tf(yld_dict_mc,yld_dict_data,None,sg.BKG_TF_MAP,quiet=False)

    # Get the values and fill the combined hist
    histo = histo_dict["nleps"][{"systematic":"nominal"}]
    err_lst_p = []
    err_lst_m = []
    for cat_name in sr_lst:
        val_sum = 0
        var_sum = 0
        # Loop over processes, get yields to fill hist with (along with summing errs)
        for proc_name in proc_lst:
            val = yld_dict_mc[cat_name]["nominal"][proc_name][0]
            histo_comb[{"process": proc_name, "cat": cat_name}] = val
            # Sum the variances so we can have stat uncertainty (also need sum of nom)
            var_sum += yld_dict_mc[cat_name]["nominal"][proc_name][1]
            val_sum += val
        err_lst_p.append(val_sum + np.sqrt(var_sum))
        err_lst_m.append(val_sum - np.sqrt(var_sum))

    # Append a 0 err for overflow bin
    err_lst_p.append(0)
    err_lst_m.append(0)

    # Make plot
    fig = make_single_fig(histo_comb,err_p=err_lst_p,err_m=err_lst_m,ylims=[0,9])
    fig.savefig("sr_comb_plot.png")


# Main function for making CR plots
def make_plots(histo_dict,grouping_mc,grouping_data,save_dir_path):

    for var_name in histo_dict.keys():
        #print(f"\n{var_name}")
        if "counts" in var_name: continue
        #if var_name != "njets": continue # TMP
        if var_name == "nbtagsm": continue # TMP
        #if var_name not in BDT_INPUT_LST and "bdt" not in var_name: continue # TMP
        #if var_name not in TMP_VAR_LST: continue # TMP
        histo = histo_dict[var_name]

        # Just plot nominal syst for now
        histo = histo[{"systematic":"nominal"}]

        # Rebin if continous variable
        if var_name not in ["njets","nbtagsl","nleps"]:
            histo = rebin(histo,6)

        # Group SR procs together
        #grouping_sr_procs = {"sr_4l_sf":["sr_4l_sf_A","sr_4l_sf_B","sr_4l_sf_C"],"sr_4l_of":["sr_4l_of_1","sr_4l_of_2","sr_4l_of_3","sr_4l_of_4"]}
        #histo = group(histo,"category","category",grouping_sr_procs)


        # Loop over categories and make plots for each
        for cat_name in histo.axes["category"]:
            #if cat_name not in ["cr_4l_sf","cr_4l_btag_of"]: continue # TMP
            #if "cr" not in cat_name: continue # TMP
            #if "bdt" in cat_name: continue # TMP
            if cat_name not in ["cr_4l_sf","sr_4l_of_incl","sr_4l_sf_incl","sr_4l_bdt_sf_presel","sr_4l_bdt_sf_trn","sr_4l_bdt_of_presel"]: continue # TMP
            #print(cat_name)

            histo_cat = histo[{"category":cat_name}]

            # Group the mc and data samples
            histo_grouped_mc = group(histo_cat,"process","process_grp",grouping_mc)
            histo_grouped_data = group(histo_cat,"process","process_grp",grouping_data)

            ######
            #if (cat_name == "cr_4l_btag_sf_offZ_met80") and var_name == "nleps":
            #if ("met80" in cat_name and ("ee" in cat_name or "mm" in cat_name)) and var_name == "nleps":
            #if ("ee" in cat_name or "mm" in cat_name) or ("cutflow" in cat_name and var_name == "nleps"):
            #    #print("mc\n",histo_grouped_mc)
            #    #print("data\n",histo_grouped_data)
            #    #print("val mc\n",histo_grouped_mc.values(flow=True))
            #    #print("val data\n",histo_grouped_data.values(flow=True))
            #    ##print("var mc\n",(histo_grouped_mc.variances(flow=True)))
            #    ##print("var data\n",(histo_grouped_data.variances(flow=True)))
            #    print("val mc\n",sum(histo_grouped_mc.values(flow=True)))
            #    print("val data\n",sum(histo_grouped_data.values(flow=True)))
            #    #print("var mc\n",(histo_grouped_mc.variances(flow=True)))
            #    #print("var data\n",(histo_grouped_data.variances(flow=True)))
            #continue
            #if (cat_name == "cr_4l_sf" or cat_name == "cr_4l_btag_of" or cat_name=="cr_4l_btag_sf_offZ_met80") and var_name == "nleps":
            #    print(f"\n{cat_name} {var_name}:")
            #    print("Yields")
            #    print("mc:",sum(sum(histo_grouped_mc.values(flow=True))))
            #    print("data:",sum(sum(histo_grouped_data.values(flow=True))))
            #    print("data/mc:",sum(sum(histo_grouped_data.values(flow=True))) / sum(sum(histo_grouped_mc.values(flow=True))) )
            #continue
            #####

            # Merge overflow into last bin (so it shows up in the plot)
            histo_grouped_data = merge_overflow(histo_grouped_data)
            histo_grouped_mc = merge_overflow(histo_grouped_mc)

            # Make figure
            title = f"{cat_name}_{var_name}"
            if "cr" in title:
                fig = make_cr_fig(histo_grouped_mc,histo_grouped_data,title=title)
            else:
                fig = make_cr_fig(histo_grouped_mc,title=title)

            # Save
            save_dir_path_cat = os.path.join(save_dir_path,cat_name)
            if not os.path.exists(save_dir_path_cat): os.mkdir(save_dir_path_cat)
            fig.savefig(os.path.join(save_dir_path_cat,title+".pdf"))
            fig.savefig(os.path.join(save_dir_path_cat,title+".png"))

            make_html(os.path.join(os.getcwd(),save_dir_path_cat))



###### Transfer factors for background ######

# TODO move to same function as used in datacard maker
# Function for getting a dict with NSF and TF etc
def get_background_dict(yld_dict_mc,yld_dict_data,bkg_proc,cr_name,sr_name):

    # Get the sum of all other contributions
    bkg_all_but_bkg_of_interest = [0,0]
    for proc in yld_dict_mc.keys():
        if proc != bkg_proc:
            bkg_all_but_bkg_of_interest[0] += yld_dict_mc[proc][cr_name][0]
            bkg_all_but_bkg_of_interest[1] += yld_dict_mc[proc][cr_name][1]

    n_cr = yld_dict_data["data"][cr_name][0] - bkg_all_but_bkg_of_interest[0]
    m_cr = yld_dict_mc[bkg_proc][cr_name][0]
    m_sr = yld_dict_mc[bkg_proc][sr_name][0]

    n_cr_err = np.sqrt(yld_dict_data["data"][cr_name][1] + bkg_all_but_bkg_of_interest[1])
    m_cr_err = np.sqrt(yld_dict_mc[bkg_proc][cr_name][1])
    m_sr_err = np.sqrt(yld_dict_mc[bkg_proc][sr_name][1])

    out_dict = {
        "n_sr_est" : [n_cr*(m_sr/m_cr) , (n_cr*(m_sr/m_cr))*np.sqrt((n_cr_err/n_cr)**2 + (m_sr_err/m_sr)**2 + (m_cr_err/m_cr)**2)],
        "m_sr"     : [m_sr , m_sr_err],
        "n_cr"     : [n_cr , n_cr_err],
        "m_cr"     : [m_cr , m_cr_err],
        "tf"       : [m_sr/m_cr , (m_sr/m_cr)*np.sqrt((m_sr_err/m_sr)**2+(m_cr_err/m_cr)**2)],
        "nsf"      : [n_cr/m_cr , (n_cr/m_cr)*np.sqrt((n_cr_err/n_cr)**2+(m_cr_err/m_cr)**2)],
    }
    return out_dict

# TODO move to same function as used in datacard maker
# Wrapper around the background estimation of TFs and yields
def do_background_estimation(yld_dict_mc,yld_dict_data):

    # Map between short name and name to display in table
    kname_dict = {
        "n_sr_est" : "$N_{SR \\rm \\; est} = TF \\cdot N_{CR}$",
        "m_sr"     : "$MC_{SR}$",
        "n_cr"     : "$N_{CR}$",
        "m_cr"     : "$MC_{CR}$",
        "tf"       : "TF",
        "nsf"      : "NSF",
    }


    print_dict = {}

    # Do the ttZ and ZZ estimation for cut-based SRs
    print_dict["ttZ SR_OF"] = get_background_dict(yld_dict_mc,yld_dict_data,"ttZ","cr_4l_btag_of","sr_of_all")
    print_dict["ZZ SR_OF"]  = get_background_dict(yld_dict_mc,yld_dict_data,"ZZ","cr_4l_sf","sr_of_all")
    print_dict["ttZ SR_SF"] = get_background_dict(yld_dict_mc,yld_dict_data,"ttZ","cr_4l_btag_sf_offZ_met80","sr_sf_all")
    print_dict["ZZ SR_SF"]  = get_background_dict(yld_dict_mc,yld_dict_data,"ZZ","cr_4l_sf","sr_sf_all")

    # Do the ttZ and ZZ estimation for BDT SRs
    for bdt_sr in SR_OF_BDT:
        print_dict[f"ttZ {bdt_sr}"] = get_background_dict(yld_dict_mc,yld_dict_data,"ttZ","cr_4l_btag_of",bdt_sr)
        print_dict[f"ZZ {bdt_sr}"] = get_background_dict(yld_dict_mc,yld_dict_data,"ZZ","cr_4l_sf",bdt_sr)
    for bdt_sr in SR_SF_BDT:
        print_dict[f"ttZ {bdt_sr}"] = get_background_dict(yld_dict_mc,yld_dict_data,"ttZ","cr_4l_btag_sf_offZ_met80",bdt_sr)
        print_dict[f"ZZ {bdt_sr}"] = get_background_dict(yld_dict_mc,yld_dict_data,"ZZ", "cr_4l_sf",bdt_sr)

    for k in print_dict:
        print(k,"\n",print_dict[k])

    # Print the dicts so we can look at the values
    # First replace key names with more descriptive names for printing
    for tf_dict_name in print_dict.keys():
        for kname in kname_dict.keys():
            print_dict[tf_dict_name][kname_dict[kname]] = print_dict[tf_dict_name].pop(kname)
    mlt.print_latex_yield_table(
        print_dict,
        tag="NSFs and TFs for ttZ and ZZ SR estimations",
        print_begin_info=True,
        print_end_info=True,
        roundat=3,
        print_errs=True,
        size="footnotesize",
        hz_line_lst=[1,3,3,19]
    )



################### Main ###################

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("pkl_file_path", help = "The path to the pkl file")
    parser.add_argument("-o", "--output-path", default=".", help = "The path the output files should be saved to")
    parser.add_argument('-y', "--get-yields", action='store_true', help = "Get yields from the pkl file")
    parser.add_argument('-p', "--make-plots", action='store_true', help = "Make plots from the pkl file")
    parser.add_argument('-b', "--get-backgrounds", action='store_true', help = "Get background estimations")
    parser.add_argument('-u', "--ul-year", default='all', help = "Which year to process", choices=["all","UL16APV","UL16","UL17","UL18"])
    args = parser.parse_args()

    # Get the counts from the input hiso
    histo_dict = pickle.load(gzip.open(args.pkl_file_path))

    sample_dict_mc = sg.create_mc_sample_dict(SAMPLE_DICT_BASE,args.ul_year)
    sample_dict_data = sg.create_data_sample_dict(args.ul_year)
    out_path = "plots" # Could make this an argument


    # Wrapper around the code for getting the raw counts and dump to latex table
    #counts_dict = get_yields(histo_dict,sample_dict_mc,raw_counts=True)
    #print_counts(counts_dict)
    #exit()


    # Wrapper around the code for getting the TFs and background estimation factors
    # TODO move to same function as used in datacard maker
    if args.get_backgrounds:
        yld_dict_data = get_yields(histo_dict,sample_dict_data,quiet=True,blind=True)
        yld_dict_mc   = get_yields(histo_dict,sample_dict_mc,quiet=True)
        put_cat_col_sums(yld_dict_mc, sr_sf_lst=SR_SF_CB, sr_of_lst=SR_OF_CB)
        do_background_estimation(yld_dict_mc,yld_dict_data)


    # Wrapper around the code for getting the yields for sr and bkg samples
    if args.get_yields:

        # Get the yield dict and put the extra columns and rows into it
        yld_dict = get_yields(histo_dict,sample_dict_mc)
        put_proc_row_sums(yld_dict, SR_SF_CB+SR_OF_CB)
        put_proc_row_sums(yld_dict, SR_SF_BDT+SR_OF_BDT)
        put_cat_col_sums(yld_dict, sr_sf_lst=SR_SF_CB, sr_of_lst=SR_OF_CB, tag="_cutbased")
        put_cat_col_sums(yld_dict, sr_sf_lst=SR_SF_BDT, sr_of_lst=SR_OF_BDT, tag="_bdt")
        #print(yld_dict)
        #exit()

        # Dump latex table for cut based
        hlines = [2,3,7,8]
        sr_cats_to_print = SR_SF_CB + ["sr_sf_all_cutbased"] + SR_OF_CB + ["sr_of_all_cutbased","sr_all_cutbased"]
        #sr_cats_to_print = ["sr_sf_all_cutbased" , "sr_of_all_cutbased" , "sr_all_cutbased" , "sr_4l_sf_presel" , "sr_4l_sf_trn" , "sr_4l_of_presel"] # Preselection SR categories
        procs_to_print = ["WWZ","ZH","Sig","ZZ","ttZ","tWZ","WZ","other","Bkg",SOVERROOTB,SOVERROOTSPLUSB,"Zmetric"]
        print_yields(yld_dict,sr_cats_to_print,procs_to_print,hlines=hlines)

        # Dump latex table for BDT
        #hlines = [6,7,15,16]
        #sr_cats_to_print = SR_SF_BDT + ["sr_sf_all_bdt"] + SR_OF_BDT + ["sr_of_all_bdt","sr_all_bdt"]
        #procs_to_print = ["WWZ","ZH","Sig","ZZ","ttZ","tWZ","WZ","other","Bkg",SOVERROOTB,SOVERROOTSPLUSB,"Zmetric"]
        #print_yields(yld_dict,sr_cats_to_print,procs_to_print,hlines=hlines)

        # Dump yield dict to json
        json_name = "process_yields.json" # Could be an argument
        json_name = os.path.join(out_path,json_name)
        with open(json_name,"w") as out_file: json.dump(yld_dict, out_file, indent=4)
        print(f"\nSaved json file: {json_name}\n")


    # Make plots
    if args.make_plots:
        make_plots(histo_dict,sample_dict_mc,sample_dict_data,save_dir_path=out_path)
        #make_syst_plots(histo_dict,sample_dict_mc,sample_dict_data,out_path,args.ul_year) # Check on individual systematics
        #make_sr_comb_plot(histo_dict,sample_dict_mc,sample_dict_data) # Make plot of all SR yields in one plot




if __name__ == "__main__":
    main()

