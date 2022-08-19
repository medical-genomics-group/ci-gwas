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
        0.45658958,  0.38832142, -0.50231455,  0.58969111, -0.15727911,
       -0.20246698, -1.70189016,  0.46559508, -1.37031878,  0.50335817,
       -1.84373354,  0.49275481, -0.39656091, -0.67157479,  0.77695995,
        0.35416727, -1.51828881,  0.50020567,  0.41317419,  1.27631234,
        1.82687669, -0.09960622,  0.39863089, -1.22703348,  0.52653355,
       -1.63542732, -1.02662668,  1.58892522, -0.11973822,  1.84204181,
        0.31255828,  0.78475446, -0.58589786,  1.14620151,  0.85627307,
        0.52670463,  0.59427733,  0.79937805,  0.46572471, -0.56663725,
       -1.46914434,  0.71368775,  1.05805913,  0.01829539,  0.10607949,
        1.64341788, -0.72756952, -1.9334239 ,  0.7279093 , -0.48926149,
       -0.57130376,  0.46327817,  0.18706507, -1.11782331,  0.50549205,
       -0.73647885,  1.76306029, -0.88452687, -0.87803062,  0.72588945,
        1.64001152,  1.27803087, -0.14943529,  0.25826108,  1.50585054,
        0.67861612, -0.01728863, -0.55536535,  1.59562221,  1.167628  ,
       -1.51307295, -0.00383193, -0.42370779, -1.87863986, -0.41100043,
       -1.81367018,  0.95197318, -0.73975743,  1.04836064, -0.92710255,
        0.16646248,  0.49206463, -0.88284117, -2.03398747,  1.36079246,
       -1.74988122, -0.53693566, -0.1133754 ,  0.4026085 , -0.60086764,
        0.56186095,  0.51429819, -1.83898934, -1.15938166,  0.43318672,
       -0.3754503 ,  0.61688372, -1.63697682,  1.14339566,  0.18036514,
        2.25593458, -0.00428075,  0.2891537 , -1.19854779, -1.0206399 ,
        0.11398634, -0.71673966, -0.29316723,  0.40557223, -0.78656693,
        0.68482867, -0.04577281, -1.48852947, -0.01574319, -0.71277468,
        0.08203599,  0.15804994, -1.07703532, -1.25329912, -0.49043303,
        0.78028608,  1.68663851,  1.06385981, -0.76423353, -1.10102063,
        0.44221407,  0.70545561, -0.9694241 ,  0.77676825,  0.41468681,
       -1.36152278, -0.5248138 ,  0.03148192, -1.19073915,  0.46359776,
        0.89086225,  0.33811416, -0.51108746, -1.24873889, -1.86072756,
        0.89585741,  0.11797734,  1.51830811, -0.60199368,  0.1781718 ,
        1.95745632,  1.4044111 , -1.28650017,  1.8375355 , -1.20903902,
       -0.82400362,  0.74825244,  0.49473588,  0.88233462,  0.30883597,
       -0.79047848,  1.01498844, -0.85773016,  0.54572364, -0.19446931,
       -0.494738  , -0.25115286, -2.12949141, -0.08827268,  0.5611441 ,
       -1.30159917, -0.79348615,  0.49328344, -0.72289605, -0.83652571,
       -1.35500246, -0.09469026, -1.73051649,  1.53704568,  0.00600275,
       -1.47316004,  0.61553299, -1.3077745 ,  0.2990651 , -0.19447101,
       -1.18544226, -0.02156319,  0.72778696, -0.47126363,  0.17318007,
       -0.07156136,  0.93219842,  0.36438488,  2.38959635, -1.21986312,
        1.35525337,  0.85871329,  2.10839056, -0.10750228,  1.87722918,
        0.73810008,  0.47539577, -0.39779577,  1.6412675 ,  0.00713489,
       -0.52988327,  1.62455502, -0.08780786, -0.08638467,  0.54999403,
       -1.0701355 ,  1.55494828, -1.81125184, -1.24120952,  1.25348427,
       -1.30330926, -0.33733287, -1.27163493, -0.89072149,  0.26870699,
       -0.33211756,  2.20514981, -0.4350727 ,  1.1403047 , -1.91204604,
       -0.58175319,  0.84618602, -1.42130181, -0.64936529,  0.71550518,
       -0.01343827, -0.64215755,  0.73595875, -0.17374664,  2.13379609,
        0.78092029, -0.9585283 ,  0.03378115,  1.05893212, -1.35207156,
       -0.05985615,  0.83466611, -1.096371  , -2.13041177, -0.72432445,
        0.92537438,  1.53280691,  0.4205967 , -0.20600761, -0.68768729,
       -0.89174878,  2.47841646, -0.2208713 , -0.17965625, -2.09503381,
        0.79119078,  0.90707734, -0.92291062, -0.82144702, -0.46571406,
        0.4310035 ,  0.63010606,  1.31299282,  0.01265133,  0.99804803,
       -0.1839934 , -1.78649373, -0.35720151,  0.51118765, -0.86216539,
        0.28013275,  0.55032593,  1.10446371,  0.01592324, -0.78751946,
       -1.21982803,  1.21982634, -1.94649502,  0.57298228, -0.4806906 ,
        0.19108281,  1.85810955,  1.27200517, -0.65184852,  0.31179482,
       -0.44200199,  1.07325508, -0.55425182,  0.59775725,  0.95671338,
       -0.16461858,  0.22978197, -0.190462  ,  1.30257874, -0.35673313,
       -0.18123754,  1.41339405, -0.42553817, -0.00547668, -0.96519379,
       -0.46942658,  0.66918347,  0.08573596, -1.11941095,  0.36050985,
        0.70085202, -0.68545591, -1.61361695,  1.39248221, -1.13865147,
        0.37463919,  1.15520908, -0.00567961, -0.19691573,  0.15177515,
       -0.67436446, -0.50141211, -0.64028645,  1.62195319,  0.37040179,
       -0.3921195 ,  0.74400398, -1.48007893,  0.34752985, -1.25700191,
        0.27595582,  0.8192143 , -0.11538094,  1.02988986, -1.41421594,
        1.91807142, -0.91320262, -0.89417313,  1.58069957, -1.0828592 ,
        0.40723484,  0.0711542 , -0.12848039, -0.20507546,  0.68974684,
        1.59943441,  0.35512747,  1.69020121, -0.00625351, -0.34741776,
        0.41295605, -0.39506472,  0.98050783, -1.22649129,  1.84694272,
        1.07369759,  1.2120396 , -1.33321832,  1.37309081, -0.14134206,
       -1.51873521, -1.01681443, -0.27253478,  1.66137131,  0.11610643,
       -0.52950979, -1.17119629, -0.11892857,  1.57676073, -0.08938551,
        0.69765119, -1.56112927,  2.04345058,  1.32308342, -0.02023113,
       -0.88674883, -0.79863008, -0.27632035, -0.47907814, -0.18460532,
        0.33303142, -0.69598629, -0.38247653,  1.6584241 , -0.98382136,
       -0.66027068, -1.54101403, -0.84105229,  0.72927097,  0.18367827,
       -2.32039345, -1.27431378,  0.05918579,  0.21470167, -0.05583006,
       -0.27221659, -1.34487442,  1.13262482, -0.57704591, -0.15310072,
        0.22652229, -2.11160874,  1.09210692,  0.83355654,  0.71568342,
        0.76227965,  0.85508381,  0.75451529, -0.87773258, -1.35955612,
       -0.19399791, -0.44086259, -1.00621221, -1.11635624, -0.00298953,
        0.74644743, -1.28997331, -0.03070637, -0.24752176,  1.44924259,
       -0.24281103, -0.77769108,  0.06985169, -1.17466425,  0.71992085,
        0.50142285, -1.26285335, -1.34908609, -1.04966928, -1.84511699,
        1.57344877,  0.19872045, -1.21048057, -0.91197489, -0.09603856,
       -0.73644853,  1.75802596,  0.67370848,  2.41344638, -0.20675385,
       -0.50940841,  1.1291334 ,  1.42829262,  1.59671384, -0.24831814,
        0.67067019, -0.45123411,  2.26291015, -0.23888125,  0.05547969,
       -0.80319503,  1.12180975,  0.61509528,  0.06858319,  0.9025978 ,
       -0.03325929, -0.94318225,  1.45069951,  0.55292827,  0.84083205,
        1.09002077, -1.08735719,  0.23655603, -2.12737997,  0.41418838,
       -0.94166683, -0.46740007, -0.42932998,  2.7442897 ,  1.27915947,
       -0.87307368,  0.28909928,  0.00631015,  0.06376326,  0.08346463,
       -1.07810543, -1.92238681, -1.13610569, -0.43886517,  0.59315457,
       -0.76074506, -0.04301741, -0.69645688,  0.64186276, -1.68387437,
        0.66421043,  0.4681458 , -0.5759404 ,  0.11822977, -0.72844295,
        1.93155107, -1.2177408 , -0.80631618,  0.23980655,  0.39916679,
       -0.09732962, -0.57900428,  0.3560067 , -2.21979922, -0.38286266,
        0.34048212,  0.57147153,  0.85511611, -1.22403629,  1.12417713,
        1.00772088,  0.64568013,  0.68099514,  0.14803135,  0.14428207
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
        0.09309158657114662, -0.015190826156410004, -0.17576266786050795, 0.10610965160829801,
        0.07850373641853728, 0.24081399072948173, 0.03591467883664174, -0.10084165780456639,
        -0.08848550805624861, 0.10069993242360148};
