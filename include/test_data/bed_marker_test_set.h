#pragma once

#define BMT_NUM_INDIVIDUALS 10
#define BMT_NUM_MARKERS 3
#define BMT_NUM_PHEN 2

// That should be this matrix in padded .bed format, column-major,
// without magic numbers
// [
//     [1, 2, 1],
//     [1, 1, 1],
//     [1, 0, 1],
//     [0, 2, 0],
//     [0, 1, 1],
//     [0, 0, 1],
//     [2, 2, 2],
//     [1, 1, 1],
//     [0, 0, 0],
//     [0, 2, 0],
// ]
// Column-major
const unsigned char bmt_marker_vals[9] = {0x2a, 0xb0, 0x00, 0xcb, 0xb2, 0x0c, 0x2a, 0xba, 0x00};
const float bmt_marker_mean[3] = {1.4, 0.9, 1.2};
const float bmt_marker_std[3] = {0.66332496, 0.83066239, 0.6};

// Column-major
const float bmt_phen_vals[20] = {0.72031609,  0.59822862, 0.47614114,  -0.62264611, -1.110996,
                                 -1.23308348, 2.06327829, 0.59822862,  -0.86682106, -0.62264611,
                                 0.71351355,  0.2937997,  -0.12591416, -0.82543724, 0.15389508,
                                 -0.26581877, 2.25246434, 0.2937997,   -1.66486495, -0.82543724};

#define BMT2_NUM_INDIVIDUALS 100
#define BMT2_NUM_MARKERS 7
#define BMT2_NUM_PHEN 5

// Column-major
const unsigned char bmt2_marker_vals[175] = {
    0xff, 0xff, 0xff, 0xff, 0xef, 0xbf, 0xff, 0xff, 0xfe, 0xef, 0xbb,
    0xff, 0xbe, 0xee, 0xef, 0xff, 0xff, 0xfb, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xfe, 0xff, 0xfe, 0xfb, 0xeb, 0xff, 0xff, 0xff, 0xff, 0xee,
    0xbf, 0xfb, 0xbf, 0xff, 0xee, 0xef, 0xfb, 0xfb, 0xff, 0xee, 0xbf,
    0xff, 0xff, 0xff, 0xee, 0xff, 0xff, 0xff, 0xef, 0xcf, 0xbf, 0xff,
    0xff, 0xbe, 0xef, 0xef, 0xbf, 0xbb, 0xbf, 0xff, 0xfa, 0x3f, 0xbf,
    0xfb, 0xff, 0xff, 0xef, 0xfe, 0xff, 0xee, 0xfe, 0xfa, 0xff, 0xbf,
    0xef, 0xea, 0xef, 0xff, 0xfe, 0xef, 0xff, 0xeb, 0xff, 0xff, 0xfe,
    0xbf, 0xff, 0xef, 0xff, 0xff, 0xf2, 0xff, 0xfc, 0xfb, 0xff, 0xcf,
    0xff, 0xeb, 0xef, 0xbf, 0xfe, 0xce, 0xef, 0xff, 0xee, 0xbe, 0xff,
    0xfe, 0xee, 0xfb, 0x3f, 0xf3, 0xff, 0xbe, 0xfa, 0xff, 0xf2, 0xeb,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xfa, 0xff, 0xeb, 0xff, 0xbb, 0xfb,
    0xfb, 0xef, 0xff, 0xff, 0xbf, 0xfb, 0xff, 0xfb, 0xff, 0xff, 0xff,
    0xff, 0xeb, 0xff, 0xef, 0xef, 0xff, 0xff, 0xef, 0xff, 0xff, 0xfb,
    0xef, 0xff, 0xff, 0xff, 0xff, 0xef, 0xff, 0xbf, 0xfe, 0xff, 0xbf,
    0xbf, 0xff, 0xff, 0xff, 0xfa, 0xfb, 0xff, 0xfc, 0xbb, 0xff};
const float bmt2_marker_mean[7] = {0.13, 0.19, 0.25, 0.21, 0.3 , 0.16, 0.15};
const float bmt2_marker_std[7] = {
    0.33630343, 0.3923009, 0.4769696,  0.47528939, 0.53851648, 0.36660606, 0.38405729};
// Column-major
const float bmt2_phen_vals[500] = {
        0.45658955,  0.3883214 , -0.5023145 ,  0.5896911 , -0.1572791 ,
       -0.20246696, -1.70189   ,  0.46559504, -1.3703188 ,  0.5033582 ,
       -1.8437334 ,  0.4927548 , -0.39656088, -0.6715747 ,  0.7769599 ,
        0.35416725, -1.5182887 ,  0.50020564,  0.41317415,  1.2763122 ,
        1.8268765 , -0.09960622,  0.39863086, -1.2270334 ,  0.52653354,
       -1.6354271 , -1.0266266 ,  1.588925  , -0.11973821,  1.8420416 ,
        0.31255826,  0.7847544 , -0.5858978 ,  1.1462014 ,  0.856273  ,
        0.5267046 ,  0.5942773 ,  0.79937804,  0.46572468, -0.56663716,
       -1.4691442 ,  0.7136877 ,  1.058059  ,  0.01829539,  0.10607949,
        1.6434177 , -0.72756946, -1.9334239 ,  0.72790927, -0.48926145,
       -0.57130367,  0.4632781 ,  0.18706506, -1.1178232 ,  0.50549203,
       -0.73647875,  1.7630601 , -0.8845268 , -0.87803054,  0.7258894 ,
        1.6400113 ,  1.2780308 , -0.14943528,  0.25826108,  1.5058503 ,
        0.6786161 , -0.01728863, -0.55536526,  1.5956221 ,  1.1676278 ,
       -1.5130728 , -0.00383193, -0.42370775, -1.8786397 , -0.4110004 ,
       -1.8136702 ,  0.95197314, -0.7397573 ,  1.0483606 , -0.92710245,
        0.16646247,  0.4920646 , -0.8828411 , -2.0339873 ,  1.3607923 ,
       -1.7498811 , -0.53693557, -0.1133754 ,  0.40260845, -0.60086757,
        0.5618609 ,  0.5142982 , -1.8389893 , -1.1593816 ,  0.43318668,
       -0.37545028,  0.6168837 , -1.6369767 ,  1.1433955 ,  0.18036513,
        2.2580855 , -0.5046296 ,  0.2938357 , -1.1919513 , -1.0142725 ,
        0.11889374, -0.7107633 , -0.28773588,  0.4101044 , -0.78050065,
        0.68900144, -0.04065982, -1.4815598 , -0.01066885, -0.7068034 ,
        0.08698451,  0.16290063, -1.0705953 , -1.2466321 , -0.48474786,
        0.784336  ,  1.6895221 ,  1.0675449 , -0.758196  , -1.0945497 ,
        0.44669908,  0.7096018 , -0.9631225 ,  0.7808227 ,  0.41920725,
       -1.3547165 , -0.5190844 ,  0.03649549, -1.1841527 ,  0.46805525,
        0.89476985,  0.34273314, -0.5053757 , -1.2420778 , -1.8532791 ,
        0.8997586 ,  0.1228796 ,  1.5214083 , -0.59616494,  0.1829966 ,
        1.9599913 ,  1.407658  , -1.2797905 ,  1.8402249 , -1.2024289 ,
       -0.81788915,  0.7523436 ,  0.49915326,  0.8862532 ,  0.31349263,
       -0.78440714,  1.0187364 , -0.85157233,  0.5500754 , -0.18916497,
       -0.4890473 , -0.24577557, -2.121697  , -0.08310501,  0.56547606,
       -1.29487   , -0.787411  ,  0.4977027 , -0.71691173, -0.83039516,
       -1.3482047 , -0.08951432, -1.7232354 ,  1.5401218 ,  0.01104911,
       -1.4662101 ,  0.6197949 , -1.3010374 ,  0.30373433, -0.18916668,
       -1.1788627 , -0.01648137,  0.73190445, -0.46560314,  0.1780113 ,
       -0.06641519,  0.936053  ,  0.36897007,  2.3915753 , -1.2132392 ,
        1.3585634 ,  0.86266226,  2.1107314 , -0.10230985,  1.8798674 ,
        0.7422043 ,  0.479838  , -0.39222977,  1.6442095 ,  0.01217979,
       -0.5298833 ,  1.6245551 , -0.08780789, -0.08638468,  0.54999405,
       -1.0701356 ,  1.5549483 , -1.811252  , -1.2412096 ,  1.2534842 ,
       -1.3033094 , -0.3373329 , -1.2716349 , -0.8907215 ,  0.26870698,
       -0.3321176 ,  2.20515   , -0.43507272,  1.1403048 , -1.9120461 ,
       -0.58175325,  0.84618604, -1.4213018 , -0.6493653 ,  0.7155052 ,
       -0.01343828, -0.6421576 ,  0.73595876, -0.17374666,  2.1337962 ,
        0.7809203 , -0.95852834,  0.03378114,  1.0589321 , -1.3520716 ,
       -0.05985617,  0.83466613, -1.096371  , -2.130412  , -0.72432446,
        0.9253744 ,  1.532807  ,  0.4205967 , -0.20600763, -0.68768734,
       -0.89174885,  2.4784164 , -0.22087133, -0.17965627, -2.0950341 ,
        0.79119086,  0.9070774 , -0.92291063, -0.8214471 , -0.4657141 ,
        0.4310035 ,  0.6301061 ,  1.3129928 ,  0.01265131,  0.99804807,
       -0.18399341, -1.7864938 , -0.35720152,  0.5111877 , -0.86216545,
        0.28013274,  0.550326  ,  1.1044638 ,  0.01592322, -0.7875195 ,
       -1.2198281 ,  1.2198263 , -1.9464952 ,  0.5729823 , -0.48069063,
        0.1910828 ,  1.8581097 ,  1.2720052 , -0.65184855,  0.31179485,
       -0.44200203,  1.0732551 , -0.55425185,  0.5977572 ,  0.95671344,
       -0.1646186 ,  0.22978196, -0.19046202,  1.3025788 , -0.35673317,
       -0.18123756,  1.4133941 , -0.42553818, -0.0054767 , -0.96519387,
       -0.46942663,  0.6691835 ,  0.08573594, -1.119411  ,  0.36050984,
        0.700852  , -0.68545586, -1.613617  ,  1.392482  , -1.1386515 ,
        0.37463915,  1.155209  , -0.00567964, -0.19691576,  0.15177512,
       -0.6743645 , -0.50141215, -0.64028645,  1.621953  ,  0.37040174,
       -0.3921195 ,  0.74400395, -1.480079  ,  0.34752983, -1.257002  ,
        0.2759558 ,  0.8192143 , -0.11538097,  1.0298897 , -1.414216  ,
        1.9180714 , -0.91320264, -0.89417315,  1.5806994 , -1.0828593 ,
        0.40723482,  0.07115416, -0.12848042, -0.20507547,  0.6897468 ,
        1.5994344 ,  0.35512742,  1.690201  , -0.00625354, -0.3474178 ,
        0.412956  , -0.39506474,  0.9805077 , -1.2264913 ,  1.8469428 ,
        1.0736976 ,  1.2120395 , -1.3332183 ,  1.3730906 , -0.14134209,
       -1.5187352 , -1.0168145 , -0.27253482,  1.6613712 ,  0.11610639,
       -0.52950984, -1.1711963 , -0.1189286 ,  1.5767605 , -0.08938554,
        0.69765115, -1.5611292 ,  2.0434506 ,  1.3230833 , -0.02023116,
       -0.88674885, -0.79863006, -0.27632037, -0.47907817, -0.18460535,
        0.3330314 , -0.6959863 , -0.38247654,  1.658424  , -0.9838214 ,
       -0.6602707 , -1.541014  , -0.8410523 ,  0.7292709 ,  0.18367825,
       -2.3203933 , -1.2743138 ,  0.05918577,  0.21470164, -0.05583009,
       -0.27221662, -1.3448744 ,  1.1326247 , -0.5770459 , -0.15310074,
        0.22652227, -2.1116087 ,  1.0921068 ,  0.83355653,  0.7156834 ,
        0.7622796 ,  0.85508376,  0.7545153 , -0.8777326 , -1.3595562 ,
       -0.19399792, -0.4408626 , -1.0062122 , -1.1163563 , -0.00298953,
        0.74644744, -1.2899733 , -0.03070637, -0.24752176,  1.4492426 ,
       -0.24281104, -0.7776911 ,  0.06985169, -1.1746643 ,  0.7199208 ,
        0.5014228 , -1.2628534 , -1.349086  , -1.0496694 , -1.845117  ,
        1.5734488 ,  0.19872044, -1.2104807 , -0.9119749 , -0.09603857,
       -0.7364486 ,  1.7580259 ,  0.6737085 ,  2.4134464 , -0.20675385,
       -0.5094084 ,  1.1291335 ,  1.4282926 ,  1.5967138 , -0.24831814,
        0.67067015, -0.4512341 ,  2.2629101 , -0.23888125,  0.05547969,
       -0.803195  ,  1.1218097 ,  0.6150953 ,  0.06858319,  0.9025978 ,
       -0.0332593 , -0.94318223,  1.4506996 ,  0.55292827,  0.84083205,
        1.0900208 , -1.0873573 ,  0.23655602, -2.12738   ,  0.41418839,
       -0.94166684, -0.4674001 , -0.42933   ,  2.7442896 ,  1.2791594 ,
       -0.8730737 ,  0.28909925,  0.00631015,  0.06376325,  0.08346462,
       -1.0781054 , -1.9223868 , -1.1361057 , -0.43886518,  0.59315455,
       -0.76074505, -0.04301741, -0.6964569 ,  0.64186275, -1.6838744 ,
        0.66421044,  0.4681458 , -0.57594043,  0.11822977, -0.72844297,
        1.9315511 , -1.2177409 , -0.8063162 ,  0.23980655,  0.3991668 ,
       -0.09732963, -0.5790043 ,  0.3560067 , -2.2197993 , -0.38286266,
        0.34048212,  0.5714715 ,  0.8551161 , -1.2240363 ,  1.1241771 ,
        1.0077208 ,  0.64568007,  0.6809951 ,  0.14803135,  0.14428207
};
// Row-major
const float bmt2_marker_corrs[21] = {
        0.06306042,  0.10184015,  0.06821211, -0.02498852, -0.1371644,
        0.15336412,  0.07140777,  0.04720982,  0.01071173, -0.00436877,
       -0.06193996,  0.00567552, -0.33084848,  0.02414851,  0.1076548 ,
       -0.03949853, -0.10308874,  0.15181038, -0.07253344,  0.16289531,
        0.08953973};
// Row-major
const float bmt2_marker_phen_corrs_kendall[35] = {
        0.00863026,  0.04711764, -0.07097485,  0.08156585,  0.16324314,
        0.03925866,  0.00739839, -0.07107905,  0.058585  ,  0.20846228,
       -0.01625403, -0.06131007,  0.03722   , -0.04979232,  0.12395052,
       -0.09545691, -0.26968717, -0.2553251 ,  0.13979574,  0.06350141,
       -0.14876796, -0.09046499,  0.27322491, -0.0119455 , -0.26747445,
       -0.12512528, -0.04748418, -0.02557523,  0.06693996, -0.001218,
       -0.04286607, -0.02623651,  0.18516187, -0.03647124, -0.14603322};
// Row-major
const float bmt2_marker_phen_corrs_pearson[35] = {
        0.00393405,  0.04614673, -0.06364793,  0.06854644,  0.14055866,
        0.01886041,  0.05001539, -0.06688912,  0.0442911 ,  0.1572077 ,
       -0.05746425, -0.02293673,  0.0267954 , -0.05645598,  0.06809136,
       -0.08503284, -0.11812583, -0.19953439,  0.08492274,  0.12475495,
       -0.09109957, -0.13154999,  0.25308736, -0.01611377, -0.19248518,
       -0.09524657, -0.03775284, -0.01394576,  0.0485203 , -0.01227726,
       -0.01410438,  0.0731126 ,  0.1425229 , -0.0414379 , -0.1589431 
};
// Row-major
const float bmt2_phen_corrs[10] = {
        0.0910091794065144, -0.015190814711046488, -0.17576266783016972, 0.10610964506206991,
        0.07019207988092789, 0.24396844719170566, 0.03809660888209857, -0.10084166118126799,
        -0.0884855163147392, 0.10069993744517854};

