import numpy as np
import matplotlib as mpl
import matplotlib
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
from dataclasses import dataclass
import json
import pandas as pd
import seaborn as sns
import queue
from scipy.io import mmread
import scipy
from glob import glob

# import scienceplots
# plt.style.use('nature')


from adjustText import adjust_text

from sklearn.linear_model import LinearRegression
import itertools

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

rng = np.random.default_rng()

BASE_INDEX = 1
XLABEL_ROTATION = 60

chr_lengths = {
    1: 248956422,
    2: 242193529,
    3: 198295559,
    4: 190214555,
    5: 181538259,
    6: 170805979,
    7: 159345973,
    8: 145138636,
    9: 138394717,
    10: 133797422,
    11: 135086622,
    12: 133275309,
    13: 114364328,
    14: 107043718,
    15: 101991189,
    16: 90338345,
    17: 83257441,
    18: 80373285,
    19: 58617616,
    20: 64444167,
    21: 46709983,
    22: 50818468,
}

global_chr_starts = {i: sum(chr_lengths[j] for j in range(1, i)) for i in range(1, 24)}

cause_gamma = pd.DataFrame(
    [
        {
            "y1": "bw".upper(),
            "y2": "t2d".upper(),
            "gamma": -0.274768438468858,
            "p_value": 0.0639138595336843,
        },
        {
            "y1": "bw".upper(),
            "y2": "ST".upper(),
            "gamma": -0.11133390989017,
            "p_value": 0.0877269352676329,
        },
        {
            "y1": "bw".upper(),
            "y2": "AT".upper(),
            "gamma": 0.0569668051462395,
            "p_value": 0.725407139890355,
        },
        {
            "y1": "bw".upper(),
            "y2": "cad".upper(),
            "gamma": -0.142129232069397,
            "p_value": 0.00142540782923592,
        },
        {
            "y1": "bmi".upper(),
            "y2": "t2d".upper(),
            "gamma": 0.753001009486836,
            "p_value": 0.00483103642108524,
        },
        {
            "y1": "bmi".upper(),
            "y2": "ST".upper(),
            "gamma": 0.0724579259674975,
            "p_value": 0.229567550936659,
        },
        {
            "y1": "bmi".upper(),
            "y2": "AT".upper(),
            "gamma": 0.127570936884617,
            "p_value": 0.164044172388618,
        },
        {
            "y1": "bmi".upper(),
            "y2": "cad".upper(),
            "gamma": 0.254304518656026,
            "p_value": 1.24e-04,
        },
        {
            "y1": "HT".upper(),
            "y2": "t2d".upper(),
            "gamma": 0.0133993365533514,
            "p_value": 0.964876784670837,
        },
        {
            "y1": "HT".upper(),
            "y2": "ST".upper(),
            "gamma": -0.0158921479966716,
            "p_value": 0.54766997929668,
        },
        {
            "y1": "HT".upper(),
            "y2": "AT".upper(),
            "gamma": -0.00266884948627959,
            "p_value": 1,
        },
        {
            "y1": "HT".upper(),
            "y2": "cad".upper(),
            "gamma": -0.0646821638302025,
            "p_value": 1.80e-04,
        },
        {
            "y1": "hdl".upper(),
            "y2": "t2d".upper(),
            "gamma": -0.157545489184377,
            "p_value": 0.0563712638516308,
        },
        {
            "y1": "hdl".upper(),
            "y2": "ST".upper(),
            "gamma": -0.0346983958224145,
            "p_value": 0.269303011987186,
        },
        {
            "y1": "hdl".upper(),
            "y2": "AT".upper(),
            "gamma": 0.00954960164699658,
            "p_value": 0.99836256615099,
        },
        {
            "y1": "hdl".upper(),
            "y2": "cad".upper(),
            "gamma": -0.203750511894064,
            "p_value": 4.42e-04,
        },
        {
            "y1": "ldl".upper(),
            "y2": "t2d".upper(),
            "gamma": -0.125971830503692,
            "p_value": 0.3242828409554,
        },
        {
            "y1": "ldl".upper(),
            "y2": "ST".upper(),
            "gamma": 0.0618108230339071,
            "p_value": 0.0457418341958065,
        },
        {
            "y1": "ldl".upper(),
            "y2": "AT".upper(),
            "gamma": -0.0178475072825962,
            "p_value": 0.869422530345408,
        },
        {
            "y1": "ldl".upper(),
            "y2": "cad".upper(),
            "gamma": 0.363195509520469,
            "p_value": 6.29e-12,
        },
        {
            "y1": "TRIG".upper(),
            "y2": "t2d".upper(),
            "gamma": 0.181573783879022,
            "p_value": 0.15119741082514,
        },
        {
            "y1": "TRIG".upper(),
            "y2": "ST".upper(),
            "gamma": 0.0081104586416291,
            "p_value": 0.98962812259515,
        },
        {
            "y1": "TRIG".upper(),
            "y2": "AT".upper(),
            "gamma": -0.0917176821192589,
            "p_value": 0.167567891536276,
        },
        {
            "y1": "TRIG".upper(),
            "y2": "cad".upper(),
            "gamma": 0.281762788578154,
            "p_value": 0.0853378458891273,
        },
        {
            "y1": "ALC".upper(),
            "y2": "t2d".upper(),
            "gamma": 0.0692860734762807,
            "p_value": 0.993922398394544,
        },
        {
            "y1": "ALC".upper(),
            "y2": "ST".upper(),
            "gamma": 0.125010784103481,
            "p_value": 0.253475957376267,
        },
        {
            "y1": "ALC".upper(),
            "y2": "AT".upper(),
            "gamma": -0.070909459789924,
            "p_value": 0.807946004128504,
        },
        {
            "y1": "ALC".upper(),
            "y2": "cad".upper(),
            "gamma": 0.0369981351799616,
            "p_value": 0.815004206171256,
        },
        {
            "y1": "SMK".upper(),
            "y2": "t2d".upper(),
            "gamma": 0.153669125333255,
            "p_value": 0.671514364328843,
        },
        {
            "y1": "SMK".upper(),
            "y2": "ST".upper(),
            "gamma": 0.287713602378084,
            "p_value": 0.0230622332612941,
        },
        {
            "y1": "SMK".upper(),
            "y2": "AT".upper(),
            "gamma": 0.129522157700426,
            "p_value": 0.501793441537144,
        },
        {
            "y1": "SMK".upper(),
            "y2": "cad".upper(),
            "gamma": 0.479817266067469,
            "p_value": 9.21e-08,
        },
        {
            "y1": "bfp".upper(),
            "y2": "t2d".upper(),
            "gamma": 0.0647735280855574,
            "p_value": 0.999860686969604,
        },
        {
            "y1": "bfp".upper(),
            "y2": "ST".upper(),
            "gamma": 0.0044428460841386,
            "p_value": 0.999999416083672,
        },
        {
            "y1": "bfp".upper(),
            "y2": "AT".upper(),
            "gamma": 0.0957086524513934,
            "p_value": 0.569696461982255,
        },
        {
            "y1": "bfp".upper(),
            "y2": "cad".upper(),
            "gamma": 0.133460312029649,
            "p_value": 0.328145432358914,
        },
        {
            "y1": "fg".upper(),
            "y2": "t2d".upper(),
            "gamma": 1.32033326715937,
            "p_value": 0.0133860198514224,
        },
        {
            "y1": "fg".upper(),
            "y2": "ST".upper(),
            "gamma": 0.00957326536858055,
            "p_value": 0.992341254004135,
        },
        {
            "y1": "fg".upper(),
            "y2": "AT".upper(),
            "gamma": -0.200412032138845,
            "p_value": 0.34591203996851,
        },
        {
            "y1": "fg".upper(),
            "y2": "cad".upper(),
            "gamma": 0.111579686635068,
            "p_value": 0.247491416509988,
        },
        {
            "y1": "dbp".upper(),
            "y2": "t2d".upper(),
            "gamma": 0.0204781118822844,
            "p_value": 0.020117240438985,
        },
        {
            "y1": "dbp".upper(),
            "y2": "ST".upper(),
            "gamma": 0.0318219908469929,
            "p_value": 0.00113343221355482,
        },
        {
            "y1": "dbp".upper(),
            "y2": "AT".upper(),
            "gamma": 0.00387409082195085,
            "p_value": 0.593576411806745,
        },
        {
            "y1": "dbp".upper(),
            "y2": "cad".upper(),
            "gamma": 0.0371147380330953,
            "p_value": 1.32e-26,
        },
        {
            "y1": "sbp".upper(),
            "y2": "t2d".upper(),
            "gamma": 0.0142748064016252,
            "p_value": 0.0101128632377118,
        },  # changed from sbp,
        {
            "y1": "sbp".upper(),
            "y2": "ST".upper(),
            "gamma": 0.0216220295177384,
            "p_value": 4.12e-09,
        },  # changed from sbp,
        {
            "y1": "sbp".upper(),
            "y2": "AT".upper(),
            "gamma": 0.00246673566316583,
            "p_value": 0.755589637124656,
        },  # changed from sbp,
        {
            "y1": "sbp".upper(),
            "y2": "cad".upper(),
            "gamma": 0.024807681904989,
            "p_value": 6.59e-31,
        },
    ]
)

diseases = set(["AT", "ST", "T2D", "CAD"])

risk_factors = set(
    [
        "BMI",
        "HT",
        "ALC",
        "SMK",
        "CHOL",
        "DBP",
        "GLU",
        "HDL",
        "HbA1c",
        "LDL",
        "SBP",
        "TRIG",
        "WHR",
    ]
)


def get_pheno_codes(phen_path):
    with open(phen_path, "r") as fin:
        header = fin.readline()
    return header.strip().split("\t")[2:]


def get_block_out_stems(blockpath):
    res = []
    with open(blockpath, "r") as fin:
        for line in fin:
            line = line.strip()
            fields = line.split("\t")
            res.append(f"{fields[0]}_{fields[1]}_{fields[2]}")
    return res


def load_mdim(basepath: str):
    with open(basepath + ".mdim", "r") as fin:
        fields = fin.readline().strip().split("\t")
    return [int(f) for f in fields]


def make_dm_ix_to_sm_ix(num_m: int, num_p: int, marker_offset: int) -> np.array:
    ixs = np.arange(num_m + num_p)
    new_ixs = np.zeros_like(ixs)
    new_ixs[np.where(ixs < num_m)] = (
        ixs[np.where(ixs < num_m)] + marker_offset + num_p + BASE_INDEX
    )
    new_ixs[np.where(ixs >= num_m)] = ixs[np.where(ixs >= num_m)] - num_m + BASE_INDEX
    return new_ixs


def load_mat_sparse(
    basepath: str, num_m: int, num_p: int, marker_offset, dtype, suffix
):
    res = {}
    n = num_m + num_p
    dm = np.fromfile(basepath + suffix, dtype=dtype)
    dm = dm.reshape(n, n)
    dm2sm = make_dm_ix_to_sm_ix(num_m, num_p, marker_offset)
    nz = np.where(dm != 0)
    nz_sm_i = dm2sm[nz[0]]
    nz_sm_j = dm2sm[nz[1]]
    for tix in range(len(nz_sm_i)):
        res[(nz_sm_i[tix], nz_sm_j[tix])] = dm[nz[0][tix], nz[1][tix]]
    return res


def load_sepset_sparse(
    basepath: str, num_m: int, num_p: int, max_level: int, marker_offset
):
    res = {}
    n = num_m + num_p
    dm = np.fromfile(basepath + ".sep", dtype=np.int32)
    dm2sm = make_dm_ix_to_sm_ix(num_m, num_p, marker_offset)
    for i in range(n):
        for j in range(n):
            new_s = []
            for six in range(max_level):
                e = dm[(i * n + j) * max_level + six]
                if e == -1:
                    break
                new_s.append(dm2sm[e])
            if new_s:
                # print(i, j, [dm[(i * n + j) * max_level + six] for six in range(max_level)])
                if dm2sm[i] in new_s or dm2sm[j] in new_s:
                    raise ValueError("SepSet(x, y) contains x or y")
                res[(dm2sm[i], dm2sm[j])] = set(new_s)

    return res


def load_corr_sparse(basepath: str, num_m: int, num_p: int, marker_offset):
    return load_mat_sparse(basepath, num_m, num_p, marker_offset, np.float32, ".corr")


def load_skeleton(basepath: str, num_m: int, num_p: int, marker_offset):
    dtype = np.int32
    suffix = ".adj"
    res = {}
    n = num_m + num_p
    dm = np.fromfile(basepath + suffix, dtype=dtype)
    dm = dm.reshape(n, n)
    dm2sm = make_dm_ix_to_sm_ix(num_m, num_p, marker_offset)
    not_zero = np.where(dm != 0)

    for i, j in zip(dm2sm[not_zero[0]], dm2sm[not_zero[1]]):
        if i not in res:
            res[i] = set()
        if j not in res:
            res[j] = set()
        res[i].add(j)
        res[j].add(j)
    return res


def load_adj_sparse(basepath: str, num_m: int, num_p: int, marker_offset=0):
    return load_mat_sparse(basepath, num_m, num_p, marker_offset, np.int32, ".adj")


def load_global_marker_indices(
    basepath: str,
    num_m: int,
    num_p: int,
    selected_marker_offset=0,
    global_marker_offset=0,
):
    global_marker_indices = {}
    rel_to_block = np.fromfile(basepath + ".ixs", dtype=np.int32)
    dm2sm = make_dm_ix_to_sm_ix(num_m, num_p, selected_marker_offset)
    for dm_ix in range(len(dm2sm)):
        sm_ix = dm2sm[dm_ix]
        if sm_ix >= (num_p + BASE_INDEX):
            global_marker_indices[sm_ix] = rel_to_block[dm_ix] + global_marker_offset
    return global_marker_indices


def add_gmi(a, b):
    a.update(b)


def add_scm(a, b):
    a.update(b)


def add_sam(a, b, num_p: int):
    # remove non-intersection p x p links
    for i in range(num_p):
        for j in range(num_p):
            if (i, j) in a and (i, j) not in b:
                del a[(i, j)]

    for (i, j), k in b.items():
        if i >= num_p or j >= num_p:
            a[(i, j)] = k


def add_ssm(a, b):
    for k, v in b.items():
        if k in a:
            a[k] = a[k] | v
        else:
            a[k] = v


def merge_block_outputs(blockfile: str, outdir: str):
    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    try:
        bo = BlockOutput(basepaths[0])
        marker_offset = bo.num_markers()
        global_marker_offset = bo.block_size()
        sam = bo.sam()
        scm = bo.scm()
        ssm = bo.ssm()
        gmi = bo.gmi()
    except FileNotFoundError:
        path = basepaths[0]
        print(f"Missing: {path}")
        global_marker_offset = block_size(path)
        marker_offset = 0
        sam = {}
        scm = {}
        ssm = {}
        gmi = {}

    for path in basepaths[1:]:
        try:
            bo = BlockOutput(path, marker_offset, global_marker_offset)
        except FileNotFoundError:
            print(f"Missing: {path}")
            global_marker_offset += block_size(path)
            continue
        add_sam(sam, bo.sam(), bo.num_phen())
        add_scm(scm, bo.scm())
        add_ssm(ssm, bo.ssm())
        add_gmi(gmi, bo.gmi())
        marker_offset += bo.num_markers()
        global_marker_offset += bo.block_size()

    return GlobalBdpcResult(
        sam, scm, ssm, gmi, marker_offset + bo.num_phen(), bo.num_phen(), bo.max_level()
    )


def global_upm(blockfile: str, outdir: str, max_depth=np.inf):
    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    bo = BlockOutput(basepaths[0])
    marker_offset = bo.num_markers()

    upm = bo.union_pleiotropy_mat(max_depth=max_depth)

    for path in basepaths[1:]:
        try:
            bo = BlockOutput(path, marker_offset)
        except FileNotFoundError:
            continue
        print("processing block: ", bo)
        marker_offset += bo.num_markers()
        for k, v in bo.union_pleiotropy_mat(max_depth=max_depth).items():
            if k in upm:
                upm[k] += v
            else:
                upm[k] = v
    return upm


def global_epm(blockfile: str, outdir: str, max_depth=np.inf):
    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    bo = BlockOutput(basepaths[0])
    marker_offset = bo.num_markers()

    epm = bo.exclusive_pleiotropy_mat(max_depth=max_depth)

    for path in basepaths[1:]:
        try:
            bo = BlockOutput(path, marker_offset)
        except FileNotFoundError:
            continue
        print("processing block: ", bo)
        marker_offset += bo.num_markers()
        for k, v in bo.exclusive_pleiotropy_mat(max_depth=max_depth).items():
            if k in epm:
                epm[k] += v
            else:
                epm[k] = v
    return epm


def global_eps(blockfile: str, outdir: str, max_depth=np.inf):
    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    bo = BlockOutput(basepaths[0])
    marker_offset = bo.num_markers()

    eps = bo.exclusive_pleiotropy_sets(max_depth=max_depth)

    for path in basepaths[1:]:
        bo = BlockOutput(path, marker_offset)
        marker_offset += bo.num_markers()
        for k, v in bo.exclusive_pleiotropy_sets(max_depth=max_depth).items():
            if k in eps:
                eps[k].update(v)
            else:
                eps[k] = v
    return eps


def global_ancestor_sets(blockfile: str, outdir: str, reduced_indices=False, depth=1):
    if not reduced_indices:
        gr = merge_block_outputs(blockfile, outdir)

    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    bo = BlockOutput(basepaths[0])
    marker_offset = bo.num_markers()

    eps = bo.pheno_ancestor_sets(depth)

    for path in basepaths[1:]:
        try:
            bo = BlockOutput(path, marker_offset)
        except FileNotFoundError as e:
            print(e)
            continue
        marker_offset += bo.num_markers()
        for k, v in bo.pheno_ancestor_sets(depth).items():
            if not reduced_indices:
                v = set(gr.gmi[ix] for ix in v)
            if k in eps:
                eps[k].update(v)
            else:
                eps[k] = v
    return eps


def global_parent_sets(blockfile: str, outdir: str, reduced_indices=False):
    if not reduced_indices:
        gr = merge_block_outputs(blockfile, outdir)

    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    bo = BlockOutput(basepaths[0])
    marker_offset = bo.num_markers()

    eps = bo.pheno_direct_parents()

    for path in basepaths[1:]:
        bo = BlockOutput(path, marker_offset)
        marker_offset += bo.num_markers()
        for k, v in bo.pheno_direct_parents().items():
            if not reduced_indices:
                v = set([int(gr.gmi[ix]) for ix in v])
            if k in eps:
                eps[k].update(v)
            else:
                eps[k] = v
    return eps


def pag_pheno_parent_sets(pag, num_phen, neighbor_fn, depth=1):
    """Compute upper bound of markers that could affect each phenotype or combination of phenotypes"""
    res = {}
    plr = pag.tolil().rows
    phens = set(np.arange(num_phen))
    for pix in range(num_phen):
        visited = set()
        nq = queue.Queue()
        q = queue.Queue()
        q.put(pix)
        for _ in range(depth):
            while not q.empty():
                v1 = q.get()
                for v2 in plr[v1]:
                    if (
                        not v2 in phens
                        and not v2 in visited
                        and neighbor_fn(pag, v1, v2)
                    ):
                        nq.put(v2)
                        visited.add(v2)
            q = nq
            nq = queue.Queue()
        res[pix] = visited
    return res


def is_child(pag, v1, v2):
    return pag[v2, v1] == 2 and pag[v1, v2] == 3


def is_possible_child(pag, v1, v2):
    return is_child(pag, v1, v2) or (pag[v2, v1] == 2 and pag[v1, v2] == 1)


def pag_exclusive_pleiotropy_sets(pag_path: str, pheno_path: str, neighbor_fn, depth=1):
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    pag = mmread(pag_path).tocsr()

    pm = pag_pheno_parent_sets(pag, num_phen, neighbor_fn, depth)
    pleiotropic_markers = set()
    res = {}
    for i in range(num_phen):
        for j in range(i + 1, num_phen):
            s = set.intersection(pm[i], pm[j])
            res[(i, j)] = s
            res[(j, i)] = s
            pleiotropic_markers.update(s)
    for i in range(num_phen):
        res[(i, i)] = pm[i] - pleiotropic_markers
    return res


def block_size(basepath: str) -> int:
    first, last = basepath.split("_")[-2:]
    return int(last) - int(first) + 1


class BlockOutput:
    def __init__(self, basepath: str, marker_offset=0, global_marker_offset=0):
        self.basepath = basepath
        self.mdim = load_mdim(basepath)
        # number of selected markers in all previous blocks
        self.marker_offset = marker_offset
        # .bim row index of first marker in block definition
        self.global_marker_offset = global_marker_offset

    def __str__(self) -> str:
        chr, first, last = self.basepath.split("/")[-1].split("_")
        return f"{chr}:{first}-{last}"

    def block_size(self) -> int:
        first, last = self.basepath.split("_")[-2:]
        return int(last) - int(first) + 1

    def skeleton(self):
        return load_skeleton(
            self.basepath, self.num_markers(), self.num_phen(), self.marker_offset
        )

    def max_level(self) -> int:
        return self.mdim[2]

    def num_markers(self) -> int:
        return self.mdim[0] - self.mdim[1]

    def num_phen(self) -> int:
        return self.mdim[1]

    def has_markers(self) -> bool:
        return self.num_markers() > 0

    def sam(self):
        return load_adj_sparse(
            self.basepath, self.num_markers(), self.num_phen(), self.marker_offset
        )

    def scm(self):
        return load_corr_sparse(
            self.basepath, self.num_markers(), self.num_phen(), self.marker_offset
        )

    def ssm(self):
        return load_sepset_sparse(
            self.basepath,
            self.num_markers(),
            self.num_phen(),
            self.max_level(),
            self.marker_offset,
        )

    def gmi(self):
        return load_global_marker_indices(
            self.basepath,
            self.num_markers(),
            self.num_phen(),
            self.marker_offset,
            self.global_marker_offset,
        )

    def marker_indices(self):
        first = self.num_phen() + self.marker_offset
        last = first + self.num_markers()
        return np.arange(first, last) + BASE_INDEX

    def pheno_indices(self):
        return np.arange(0, self.num_phen()) + BASE_INDEX

    def pheno_parents(self, max_depth=np.inf):
        """Compute upper bound of markers that could affect each phenotype or combination of phenotypes"""
        res = {}
        adj = self.skeleton()
        phens = set(self.pheno_indices())
        for pix in self.pheno_indices():
            visited = set()
            q = queue.Queue()
            next_q = queue.Queue()
            q.put(pix)
            depth = 0
            while depth < max_depth:
                while not q.empty():
                    v1 = q.get()
                    for v2 in adj[v1]:
                        if not v2 in phens and not v2 in visited:
                            next_q.put(v2)
                            visited.add(v2)
                if next_q.empty():
                    break
                q = next_q
                next_q = queue.Queue()
                depth += 1
            res[pix] = visited
        return res

    def pheno_ancestor_sets(self, depth: int):
        res = {}
        adj = self.sam()
        for pix in self.pheno_indices():
            res[pix] = set()
            for parent in self.marker_indices():
                if (pix, parent) in adj:
                    res[pix].add(parent)
                    q = queue.Queue()
                    q.put(parent)
                    next_q = queue.Queue()
                    for _ in range(depth - 1):
                        while not q.empty():
                            v1 = q.get()
                            for v2 in self.marker_indices():
                                if (v1, v2) in adj and not v2 in res[pix]:
                                    next_q.put(v2)
                                    res[pix].add(v2)
        return res

    def pheno_direct_parents(self):
        """For each phenotype, find all markers that are direct parents to it."""
        res = {}
        adj = self.sam()
        for pix in self.pheno_indices():
            res[pix] = set()
            for parent_candidate in self.marker_indices():
                if (pix, parent_candidate) in adj:
                    res[pix].add(parent_candidate)
        return res

    def exclusive_pleiotropy_mat(self, max_depth=np.inf):
        pm = self.pheno_parents(max_depth=max_depth)
        pleiotropic_markers = set()
        res = {}
        for i in self.pheno_indices():
            for j in range(i + 1, self.num_phen() + BASE_INDEX):
                s = set.intersection(pm[i], pm[j])
                res[(i, j)] = len(s)
                res[(j, i)] = len(s)
                pleiotropic_markers.update(s)
        for i in self.pheno_indices():
            res[(i, i)] = len(pm[i] - pleiotropic_markers)
        return res

    def union_pleiotropy_mat(self, max_depth=np.inf):
        pm = self.pheno_parents(max_depth=max_depth)
        pleiotropic_markers = set()
        res = {(i, i): set(pm[i]) for i in self.pheno_indices()}
        for i in self.pheno_indices():
            for j in range(i + 1, self.num_phen() + BASE_INDEX):
                s = set.intersection(pm[i], pm[j])
                res[(i, i)] = s.union(res[(i, i)])
                res[(j, j)] = s.union(res[(j, j)])
                res[(i, j)] = len(s)
                res[(j, i)] = len(s)
                pleiotropic_markers.update(s)
        for i in self.pheno_indices():
            res[(i, i)] = len(res[(i, i)])
        return res

    def exclusive_pleiotropy_sets(self, max_depth=np.inf):
        pm = self.pheno_parents(max_depth=max_depth)
        pleiotropic_markers = set()
        res = {}
        for i in self.pheno_indices():
            for j in range(i + 1, self.num_phen() + BASE_INDEX):
                s = set.intersection(pm[i], pm[j])
                res[(i, j)] = s
                res[(j, i)] = s
                pleiotropic_markers.update(s)
        for i in self.pheno_indices():
            res[(i, i)] = pm[i] - pleiotropic_markers
        return res


@dataclass
class GlobalBdpcResult:
    sam: dict
    scm: dict
    ssm: dict
    gmi: dict
    num_var: int
    num_phen: int
    max_level: int

    def write_mm(self, basepath: str):
        # only sam and scm at the moment
        with open(basepath + "_sam" + ".mtx", "w") as fout:
            L = len(self.sam)
            N = M = max(t[0] for t in self.sam.keys())
            # N = M = len(set([t[0] for t in self.sam.keys()]))
            fout.write("%%MatrixMarket matrix coordinate integer general\n")
            fout.write(f"{N}\t{M}\t{L}\n")
            for (t1, t2), v in self.sam.items():
                fout.write(f"{t1}\t{t2}\t{v}\n")

        with open(basepath + "_scm" + ".mtx", "w") as fout:
            L = len(self.scm)
            N = M = max(t[0] for t in self.sam.keys())
            # N = M = len(set([t[0] for t in self.scm.keys()]))
            fout.write("%%MatrixMarket matrix coordinate real general\n")
            fout.write(f"{N}\t{M}\t{L}\n")
            for (t1, t2), v in self.scm.items():
                fout.write(f"{t1}\t{t2}\t{v}\n")

        with open(basepath + ".ssm", "w") as fout:
            for (t1, t2), v in self.ssm.items():
                outline = " ".join([str(e) for e in [t1, t2] + sorted(list(v))])
                fout.write(outline + "\n")

        with open(basepath + ".mdim", "w") as fout:
            fout.write(f"{self.num_var}\t{self.num_phen}\t{self.max_level}\n")


def heatmap(
    data,
    row_labels,
    col_labels,
    ax=None,
    cbar=True,
    cbar_kw=None,
    cbarlabel="",
    xlabel=None,
    ylabel=None,
    title=None,
    bad_color=None,
    title_kw=dict(),
    cbarlabel_rotation=0,
    label_interspacing=False,
    **kwargs,
):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    if kwargs["cmap"]:
        cm = plt.get_cmap(kwargs["cmap"])
        if bad_color is None:
            cm.set_bad("white")
        else:
            cm.set_bad(bad_color)
        kwargs["cmap"] = cm

    # ax.grid(which="major", color="gray", linestyle=":", linewidth=0.5)

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    if cbar:
        # Create colorbar
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
        if cbarlabel_rotation != 0:
            labels = cbar.ax.set_yticklabels(
                cbar.ax.get_yticklabels(),
                rotation=cbarlabel_rotation,
                rotation_mode="anchor",
                ha="left",
            )
        if label_interspacing:
            for i, label in enumerate(labels):
                label.set_x(label.get_position()[0] + 0.3 + (i % 2) * 1.3)
    else:
        cbar = None

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)
    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    if title:
        ax.set_title(title, **title_kw)

    return im, cbar


def get_skeleton_pleiotropy_mat(
    outdir: str,
    blockfile: str,
    pheno_path: str,
    max_depth=np.inf,
    mat_type="exclusive",
    num_phen=None,
):
    if num_phen is None:
        p_names = get_pheno_codes(pheno_path)
        num_phen = len(p_names)

    if mat_type == "exclusive":
        gepm = global_epm(blockfile, outdir, max_depth=max_depth)
    elif mat_type == "union":
        gepm = global_upm(blockfile, outdir, max_depth=max_depth)
    else:
        raise ValueError(f"Invalid mat_type: {mat_type}")

    z = np.zeros(shape=(num_phen, num_phen))
    for (i, j), c in gepm.items():
        z[i - 1, j - 1] = c

    return z


def plot_skeleton_pleiotropy_mat_z(
    z: np.array,
    pheno_path: str,
    ax=None,
    cbar_kw=None,
    title=None,
    title_kw=dict(),
    cmap="BuPu",
    norm=None,
    cbar=True,
    **kwargs,
):
    p_names = get_pheno_codes(pheno_path)
    mask = ~np.tri(z.shape[0], k=-1, dtype=bool)
    z = np.ma.array(z, mask=mask)  # mask out the lower triangle
    cmap = plt.get_cmap(cmap)
    cmap.set_bad("w")  # default value is 'k'
    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    im, _ = heatmap(
        z,
        p_names,
        p_names,
        cmap=cmap,
        cbar_kw=cbar_kw,
        cbarlabel=r"# shared ancestral markers",
        # vmin=-max_z,
        # vmax=max_z,
        # xlabel=r"$y_2$",
        # ylabel=r"$y_1$",
        title=title,
        title_kw=title_kw,
        ax=ax,
        norm=norm,
        cbar=cbar,
        **kwargs,
    )


def plot_skeleton_pleiotropy_mat(
    outdir: str,
    blockfile: str,
    pheno_path: str,
    max_depth=np.inf,
    ax=None,
    cbar_kw=None,
    title=None,
    title_kw=dict(),
    cmap="BuPu",
):
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)

    gepm = global_epm(blockfile, outdir, max_depth=max_depth)

    z = np.zeros(shape=(num_phen, num_phen))
    for (i, j), c in gepm.items():
        z[i - 1, j - 1] = c

    mask = ~np.tri(z.shape[0], k=-1, dtype=bool)
    z = np.ma.array(z, mask=mask)  # mask out the lower triangle
    cmap = plt.get_cmap(cmap)
    cmap.set_bad("w")  # default value is 'k'
    im, _ = heatmap(
        z,
        p_names,
        p_names,
        cmap=cmap,
        cbar_kw=cbar_kw,
        cbarlabel=r"# shared parent markers",
        # vmin=-max_z,
        # vmax=max_z,
        xlabel=r"$y_2$",
        ylabel=r"$y_1$",
        title=title,
        title_kw=title_kw,
        ax=ax,
    )


def plot_pleiotropy_mat(
    pag_path: str,
    pheno_path: str,
    neighbor_fn=is_possible_child,
    depth=1,
    ax=None,
    cbar_kw=None,
    title=None,
    title_kw=dict(),
    cmap="BuPu",
):
    poss_parents = pag_exclusive_pleiotropy_sets(
        pag_path, pheno_path, neighbor_fn, depth
    )
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    z = [[0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(i + 1):
            if i != j:
                z[i][j] = len(poss_parents[(i, j)])
    z = np.array(z)
    mask = ~np.tri(z.shape[0], k=-1, dtype=bool)
    z = np.ma.array(z, mask=mask)
    # z = np.ma.masked_array(z, z == 0.0)
    cmap = plt.get_cmap(cmap)
    cmap.set_bad("w")  # default value is 'k'
    im, _ = heatmap(
        z,
        p_names,
        p_names,
        cmap=cmap,
        cbar_kw=cbar_kw,
        cbarlabel=r"# shared parent markers",
        # vmin=-max_z,
        # vmax=max_z,
        xlabel=r"$y_2$",
        ylabel=r"$y_1$",
        title=title,
        title_kw=title_kw,
        ax=ax,
    )


def load_ace_directed_only(ace_path: str, pag_path: str, pheno_path: str) -> np.array:
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    ace = mmread(ace_path).tocsr()
    pag = mmread(pag_path).tocsr()
    directed = [[False for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(num_phen):
            v1 = pag[i, j]
            v2 = pag[j, i]
            if (v1, v2) == (2, 3):
                directed[i][j] = True

    z = [[0.0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(num_phen):
            if directed[i][j]:
                z[i][j] = ace[i, j]

    return np.array(z)


def load_ace(ace_path: str, pheno_path: str) -> np.array:
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    ace = mmread(ace_path).tocsr()
    z = [[0.0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(num_phen):
            z[i][j] = ace[i, j]
    return np.array(z)


def plot_ace_directed_only(
    ace_path: str,
    pag_path: str,
    pheno_path: str,
    title=None,
    title_kw=dict(),
    cmap="bwr",
    cbarlabel=r"$ACE \: (y_1 \rightarrow y_2)$",
    cbar_kw=None,
    cbar=True,
    ax=None,
    norm=None,
    xlabel=r"$y_2$",
    ylabel=r"$y_1$",
):
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    ace = mmread(ace_path).tocsr()

    pag = mmread(pag_path).tocsr()

    directed = [[False for _ in range(num_phen)] for _ in range(num_phen)]

    for i in range(num_phen):
        for j in range(num_phen):
            v1 = pag[i, j]
            v2 = pag[j, i]
            if (v1, v2) == (2, 3):
                directed[i][j] = True

    z = [[0.0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(num_phen):
            if directed[i][j]:
                z[i][j] = ace[i, j]

    max_z = np.max(np.abs(z))
    z = np.array(z)
    z = np.ma.masked_array(z, z == 0.0)
    im, _ = heatmap(
        z,
        p_names,
        p_names,
        cmap=cmap,
        cbarlabel=cbarlabel,
        cbar_kw=cbar_kw,
        cbar=cbar,
        # vmin=-max_z,
        # vmax=max_z,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        title_kw=title_kw,
        ax=ax,
        norm=norm,
    )
    return im


def plot_ace(
    ace_path: str,
    pheno_path: str,
    title=None,
    title_kw=dict(),
    cmap="bwr",
    cbarlabel=r"$ACE \: (y_1 \rightarrow y_2)$",
    cbar_kw=None,
    cbar=True,
    ax=None,
    norm=None,
    xlabel=r"$y_2$",
    ylabel=r"$y_1$",
):
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    ace = mmread(ace_path).tocsr()
    z = [[0.0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(num_phen):
            z[i][j] = ace[i, j]
    max_z = np.max(np.abs(z))
    z = np.array(z)
    z = np.ma.masked_array(z, z == 0.0)
    im, _ = heatmap(
        z,
        p_names,
        p_names,
        cmap=cmap,
        cbarlabel=cbarlabel,
        cbar_kw=cbar_kw,
        cbar=cbar,
        # vmin=-max_z,
        # vmax=max_z,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        title_kw=title_kw,
        ax=ax,
        norm=norm,
    )
    return im


@dataclass
class EdgeEncoding:
    str_rep: list[str]
    int_rep: dict[tuple[int, int], int]
    cmap: mpl.colors.ListedColormap


all_edge_types = EdgeEncoding(
    [
        r"$y_1 \; \; \; y_2$",
        r"$y_1$ o-o $y_2$",
        r"$y_1$ <-o $y_2$",
        r"$y_1$ o-> $y_2$",
        r"$y_1$ -o $y_2$",
        r"$y_1$ o- $y_2$",
        r"$y_1$ <-> $y_2$",
        r"$y_1$ -> $y_2$",
        r"$y_1$ <- $y_2$",
        r"$y_1$ - $y_2$",
    ],
    {
        (0, 0): 0,
        (1, 1): 1,
        (1, 2): 2,
        (2, 1): 3,
        (1, 3): 4,
        (3, 1): 5,
        (2, 2): 6,
        (2, 3): 7,
        (3, 2): 8,
        (3, 3): 9,
    },
    mpl.colors.ListedColormap(
        np.array(
            [
                "#ffffff",  # white
                "#003f5c",
                "#2f4b7c",
                "#665191",
                "#a05195",
                "#d45087",
                "#f95d6a",
                "#ff7c43",
                "#ffa600",
                "#ffe300",
            ]
        )
    ),
)

simulation_edge_types = EdgeEncoding(
    [
        r"$y_1 \; \; \; y_2$",
        r"$y_1$ -> $y_2$",
        r"$y_1$ <- $y_2$",
        r"$y_1$ - $y_2$",
    ],
    {
        (0, 0): 0,
        (2, 3): 1,
        (3, 2): 2,
        (3, 3): 3,
    },
    mpl.colors.ListedColormap(
        np.array(
            [
                "#ffffff",  # white
                "#fcc006",  # marigold
                # "#a6cee3", # light blue
                "#1f78b4",  # darker blue
                "#10a674",  # blueish green
            ]
        )
    ),
)


seven_edge_types = EdgeEncoding(
    [
        r"$y_1 \; \; \; y_2$",
        r"$y_1$ <-> $y_2$",
        r"$y_1$ -> $y_2$",
        r"$y_1$ <- $y_2$",
        r"$y_1$ <-o $y_2$",
        r"$y_1$ o-> $y_2$",
        r"$y_1$ - $y_2$",
        r"$y_1$ o-o $y_2$",
    ],
    {
        (0, 0): 0,
        (2, 2): 1,
        (2, 3): 2,
        (3, 2): 3,
        (1, 2): 4,
        (2, 1): 5,
        (3, 3): 6,
        (1, 1): 7,
    },
    mpl.colors.ListedColormap(
        np.array(
            [
                "#ffffff",  # white
                "#b2df8a",  # green
                "#fcc006",  # marigold
                # "#a6cee3", # light blue
                "#1f78b4",  # darker blue
                "#510ac9",  # violet blue
                "#fd411e",  # orange red
                "#41fdfe",  # bright cyan
                "#d8dcd6",  # light grey
            ]
        )
    ),
)

six_edge_types = EdgeEncoding(
    [
        r"$y_1 \; \; \; y_2$",
        r"$y_1$ <-> $y_2$",
        r"$y_1$ -> $y_2$",
        r"$y_1$ <- $y_2$",
        r"$y_1$ <-o $y_2$",
        r"$y_1$ o-> $y_2$",
        r"$y_1$ o-o $y_2$",
    ],
    {
        (0, 0): 0,
        (2, 2): 1,
        (2, 3): 2,
        (3, 2): 3,
        (1, 2): 4,
        (2, 1): 5,
        (1, 1): 6,
    },
    mpl.colors.ListedColormap(
        np.array(
            [
                "#ffffff",  # white
                "#b2df8a",  # green
                "#fcc006",  # marigold
                # "#a6cee3", # light blue
                "#1f78b4",  # darker blue
                "#510ac9",  # violet blue
                "#fd411e",  # orange red
                "#d8dcd6",  # light grey
            ]
        )
    ),
)

five_common_edge_types = EdgeEncoding(
    [
        r"$y_1 \; \; \; y_2$",
        r"$y_1$ <-> $y_2$",
        r"$y_1$ -> $y_2$",
        r"$y_1$ <- $y_2$",
        r"$y_1$ <-o $y_2$",
        r"$y_1$ o-> $y_2$",
    ],
    {
        (0, 0): 0,
        (2, 2): 1,
        (2, 3): 2,
        (3, 2): 3,
        (1, 2): 4,
        (2, 1): 5,
    },
    mpl.colors.ListedColormap(
        np.array(
            [
                "#ffffff",  # white
                "#b2df8a",  # green
                "#fcc006",  # marigold
                # "#a6cee3", # light blue
                "#1f78b4",  # darker blue
                "#510ac9",  # violet blue
                "#fd411e",  # orange red
            ]
        )
    ),
)

two_common_edge_types = EdgeEncoding(
    [
        r"$y_1 \; \; \; y_2$",
        r"$y_1$ <-> $y_2$",
        r"$y_1$ -> $y_2$",
        r"$y_1$ <- $y_2$",
    ],
    {
        (0, 0): 0,
        (2, 2): 1,
        (2, 3): 2,
        (3, 2): 3,
    },
    mpl.colors.ListedColormap(
        np.array(
            [
                "#ffffff",  # white
                "#b2df8a",  # green
                "#fcc006",  # marigold
                # "#a6cee3", # light blue
                "#1f78b4",  # darker blue
            ]
        )
    ),
)


def plot_pag(
    pag_path: str,
    pheno_path: str,
    title=None,
    title_kw=dict(),
    cbar_kw=None,
    edge_encoding=all_edge_types,
    ax=None,
    cbar=True,
    pheno_names=None,
    pheno_subset=None,
    pheno_offset=0,
    pag=None,
):
    if pheno_names is None:
        pheno_names = get_pheno_codes(pheno_path)
    num_phen = len(pheno_names)

    if pheno_subset is None:
        pheno_indices = list(range(num_phen))
    else:
        pheno_indices = [pheno_names.index(e) for e in pheno_subset]
        num_phen = len(pheno_indices)
        pheno_names = pheno_subset

    if pag is None:
        pag = mmread(pag_path).tocsr()

    z = [[0 for _ in range(num_phen)] for _ in range(num_phen)]

    for i in range(num_phen):
        for j in range(i):
            v1 = pag[pheno_offset + pheno_indices[i], pheno_offset + pheno_indices[j]]
            v2 = pag[pheno_offset + pheno_indices[j], pheno_offset + pheno_indices[i]]
            z[i][j] = edge_encoding.int_rep[(v1, v2)]

    ne = len(edge_encoding.int_rep)

    norm = mpl.colors.BoundaryNorm(np.linspace(0, ne, ne + 1), ne)
    fmt = mpl.ticker.FuncFormatter(lambda x, pos: edge_encoding.str_rep[norm(x)])

    if cbar_kw is None:
        cbar_kw = {}
    cbar_kw["ticks"] = np.arange(ne) + 0.5
    cbar_kw["format"] = fmt

    im, _ = heatmap(
        np.array(z),
        pheno_names,
        pheno_names,
        cmap=edge_encoding.cmap,
        norm=norm,
        cbar_kw=cbar_kw,
        # cbarlabel="Edge Type",
        xlabel=r"$y_2$",
        ylabel=r"$y_1$",
        title=title,
        title_kw=title_kw,
        ax=ax,
        cbar=cbar,
        cbarlabel_rotation=-50,
        label_interspacing=False,
    )

    return im


def get_possibly_causal_paths(pag_path: str, pheno_path: str, pheno_names=None):
    if pheno_names is None:
        p_names = get_pheno_codes(pheno_path)
        num_phen = len(p_names)
    else:
        num_phen = len(pheno_names)

    pag = mmread(pag_path).tocsr()

    causal_paths = np.zeros(shape=(num_phen, num_phen))

    for start_node in range(num_phen):
        q = queue.Queue()
        q.put(start_node)
        descendants = set()
        while not q.empty():
            current_node = q.get()
            for neighbor in range(num_phen):
                link = (pag[current_node, neighbor], pag[neighbor, current_node])
                if neighbor not in descendants and (link == (2, 3) or link == (2, 1)):
                    descendants.add(neighbor)
                    q.put(neighbor)
        for j in descendants:
            causal_paths[start_node, j] = 1

    return causal_paths


def get_causal_paths(
    pag_path: str, pheno_path: str, pheno_names=None, max_path_len=np.inf
):
    if pheno_names is None:
        p_names = get_pheno_codes(pheno_path)
        num_phen = len(p_names)
    else:
        num_phen = len(pheno_names)
    pag = mmread(pag_path).tocsr()

    causal_paths = np.zeros(shape=(num_phen, num_phen))

    for start_node in range(num_phen):
        q = queue.Queue()
        q.put(start_node)
        descendants = set()
        path_len = 0
        while path_len < max_path_len:
            next_q = queue.Queue()
            while not q.empty():
                current_node = q.get()
                for neighbor in range(num_phen):
                    if neighbor not in descendants and (
                        pag[current_node, neighbor],
                        pag[neighbor, current_node],
                    ) == (2, 3):
                        descendants.add(neighbor)
                        next_q.put(neighbor)
            path_len += 1
            q = next_q
        for j in descendants:
            causal_paths[start_node, j] = 1

    return causal_paths


def plot_ace_rf_to_d(
    ace_path: str,
    pag_path: str,
    pheno_path: str,
    title=None,
    title_kw=dict(),
    ax=None,
    cbarlabel=r"$ACE \: (y_1 \rightarrow y_2)$",
    cmap="bwr",
    cbar=True,
    xlabel=None,
    ylabel=None,
    cbar_kw=None,
    norm=None,
    **kwargs,
):
    if ax is None:
        plt.figure(figsize=(3, 10))
        ax = plt.gca()

    disease_list = sorted(list(diseases))

    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)

    ace = mmread(ace_path).tocsr()
    z = [[0.0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(num_phen):
            z[i][j] = ace[i, j]

    rf_to_d_ace = {}

    cause_ys = set(cause_gamma["y1"].values)
    cause_ys.update(set(cause_gamma["y2"].values))
    reg_pnames = cause_ys.intersection(set(p_names))

    for i in range(num_phen):
        for j in range(num_phen):
            name_i = p_names[i]
            name_j = p_names[j]
            if name_i in risk_factors and name_j in diseases:
                rf_to_d_ace[(name_i, name_j)] = z[i][j]

    rf_intersection = sorted(list(risk_factors.intersection(reg_pnames)))
    rf_cg_only = sorted(list(risk_factors - set(rf_intersection)))

    risk_factor_list = rf_intersection + rf_cg_only

    rf_to_d_ace_mat = [
        [0 for _ in range(len(diseases))] for _ in range(len(risk_factors))
    ]

    for i, rf in enumerate(risk_factor_list):
        for j, ds in enumerate(disease_list):
            if (rf, ds) in rf_to_d_ace:
                rf_to_d_ace_mat[i][j] = rf_to_d_ace[(rf, ds)]

    col_labels = [r"$\bf{{{0}}}$".format(e) for e in disease_list]
    row_labels = [r"$\bf{{{0}}}$".format(e) for e in rf_intersection] + rf_cg_only

    links = get_rf_to_d_causal_path_mat(pag_path, pheno_path)

    z = np.array(rf_to_d_ace_mat)

    z[(z == 0.0) & (links != 0)] = np.nan

    im, _ = heatmap(
        z,
        risk_factor_list,
        disease_list,
        cmap=cmap,
        cbarlabel=cbarlabel,
        cbar_kw=cbar_kw,
        cbar=cbar,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        title_kw=title_kw,
        ax=ax,
        norm=norm,
        bad_color="#d8dcd6",
        **kwargs,
    )
    return im


def get_rf_to_d_causal_path_mat(
    pag_path: str,
    pheno_path: str,
):
    disease_list = sorted(list(diseases))
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    causal_paths = get_causal_paths(pag_path, pheno_path)
    ci_gwas_links = {}
    cause_ys = set(cause_gamma["y1"].values)
    cause_ys.update(set(cause_gamma["y2"].values))
    reg_pnames = cause_ys.intersection(set(p_names))
    for i in range(num_phen):
        for j in range(num_phen):
            name_i = p_names[i]
            name_j = p_names[j]
            if name_i in risk_factors and name_j in diseases:
                if causal_paths[i, j] == 1:
                    ci_gwas_links[(name_i, name_j)] = 1

    rf_intersection = sorted(list(risk_factors.intersection(reg_pnames)))
    rf_cg_only = sorted(list(risk_factors - set(rf_intersection)))
    risk_factor_list = rf_intersection + rf_cg_only
    ci_gwas_links_mat = [
        [0 for _ in range(len(diseases))] for _ in range(len(risk_factors))
    ]
    for i, rf in enumerate(risk_factor_list):
        for j, ds in enumerate(disease_list):
            if (rf, ds) in ci_gwas_links:
                ci_gwas_links_mat[i][j] = 1
    return np.array(ci_gwas_links_mat)


def plot_causal_paths(
    pag_path: str,
    pheno_path: str,
    title=None,
    title_kw=dict(),
    ax=None,
    pheno_names=None,
):
    if ax is None:
        plt.figure(figsize=(10, 10))
        ax = plt.gca()

    if pheno_names is None:
        p_names = get_pheno_codes(pheno_path)
    else:
        p_names = pheno_names

    ci_gwas_link_mat = get_causal_paths(pag_path, pheno_path, pheno_names)

    return heatmap(
        np.array(ci_gwas_link_mat),
        p_names,
        p_names,
        cbar=False,
        ax=ax,
        title=title,
        title_kw=title_kw,
        cmap="binary",
        vmax=1.5,
        xlabel=r"$y_2$",
        ylabel=r"$y_1$",
    )


def plot_causal_paths_rf_to_d(
    pag_path: str,
    pheno_path: str,
    title=None,
    title_kw=dict(),
    ax=None,
):
    if ax is None:
        plt.figure(figsize=(3, 10))
        ax = plt.gca()

    disease_list = sorted(list(diseases))

    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)

    cause_ys = set(cause_gamma["y1"].values)
    cause_ys.update(set(cause_gamma["y2"].values))
    reg_pnames = cause_ys.intersection(set(p_names))

    rf_intersection = sorted(list(risk_factors.intersection(reg_pnames)))
    rf_cg_only = sorted(list(risk_factors - set(rf_intersection)))

    risk_factor_list = rf_intersection + rf_cg_only

    ci_gwas_link_mat = get_rf_to_d_causal_path_mat(pag_path, pheno_path)

    col_labels = [r"$\bf{{{0}}}$".format(e) for e in disease_list]
    row_labels = [r"$\bf{{{0}}}$".format(e) for e in rf_intersection] + rf_cg_only

    return heatmap(
        np.array(ci_gwas_links_mat),
        risk_factor_list,
        disease_list,
        cbar=False,
        ax=ax,
        title=title,
        title_kw=title_kw,
        cmap="Greys",
        vmax=1.5,
    )


def plot_ci_gwas_cause_ace_comparison_tri(
    pag_path: str,
    ace_path: str,
    pheno_path: str,
    p_thr=0.05,
    title=None,
    title_kw=dict(),
    ax=None,
    cbar_kw=None,
    cbarlabel=r"$ACE \: (y_1 \rightarrow y_2), \ CAUSE \ \gamma$",
    max_path_len=np.inf,
):
    if ax is None:
        plt.figure(figsize=(3, 8))
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    disease_list = sorted(list(diseases))

    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)

    causal_paths = get_causal_paths(pag_path, pheno_path, max_path_len=max_path_len)
    ci_gwas_links = {}

    for i in range(num_phen):
        for j in range(num_phen):
            name_i = p_names[i]
            name_j = p_names[j]
            if name_i in risk_factors and name_j in diseases:
                if causal_paths[i, j] == 1:
                    ci_gwas_links[(name_i, name_j)] = 1

    cause_ys = set(cause_gamma["y1"].values)
    cause_ys.update(set(cause_gamma["y2"].values))
    reg_pnames = cause_ys.intersection(set(p_names))

    rf_intersection = sorted(list(risk_factors.intersection(reg_pnames)))
    rf_cg_only = sorted(list(risk_factors - set(rf_intersection)))

    risk_factor_list = rf_intersection + rf_cg_only

    ci_gwas_links_arr = np.array([False] * (len(diseases) * len(risk_factors)))
    cause_sig_arr = np.array([False] * (len(diseases) * len(risk_factors)))
    ci_gwas_ace_mat = [
        [0 for _ in range(len(diseases))] for _ in range(len(risk_factors))
    ]
    cause_ace_mat = [
        [0 for _ in range(len(diseases))] for _ in range(len(risk_factors))
    ]

    ace = load_ace(ace_path, pheno_path)

    for i, rf in enumerate(risk_factor_list):
        for j, ds in enumerate(disease_list):
            glob_rf = p_names.index(rf)
            glob_ds = p_names.index(ds)
            ci_gwas_ace_mat[i][j] = ace[glob_rf, glob_ds]
            cg = cause_gamma[
                (cause_gamma["y1"] == rf) & (cause_gamma["y2"] == ds)
            ].gamma.values
            cp = cause_gamma[
                (cause_gamma["y1"] == rf) & (cause_gamma["y2"] == ds)
            ].p_value.values
            if len(cg) == 1:
                cause_ace_mat[i][j] = cg[0]
                if cp[0] <= p_thr:
                    cause_sig_arr[i * len(diseases) + j] = True
            if (rf, ds) in ci_gwas_links:
                ci_gwas_links_arr[i * len(diseases) + j] = True

    col_labels = [r"$\bf{{{0}}}$".format(e) for e in disease_list]
    row_labels = [r"$\bf{{{0}}}$".format(e) for e in rf_intersection] + rf_cg_only

    M = len(diseases)
    N = len(risk_factors)
    x = np.arange(M + 1)
    y = np.arange(N + 1)
    xs, ys = np.meshgrid(x, y)

    triangles1 = [
        (i + j * (M + 1), i + 1 + j * (M + 1), i + (j + 1) * (M + 1))
        for j in range(N)
        for i in range(M)
    ]
    triangles2 = [
        (i + 1 + j * (M + 1), i + 1 + (j + 1) * (M + 1), i + (j + 1) * (M + 1))
        for j in range(N)
        for i in range(M)
    ]
    tri_ci_gwas_link_mask = np.array([1] * len(triangles1))
    tri_ci_gwas_link_mask[np.where(ci_gwas_links_arr)] = False
    tri_ci_gwas_link = Triangulation(
        xs.ravel() - 0.5,
        ys.ravel() - 0.5,
        triangles=np.array(triangles1),
        mask=tri_ci_gwas_link_mask,
    )

    tri_cause_sig_mask = np.array([1] * len(triangles1))
    tri_cause_sig_mask[np.where(cause_sig_arr)] = False
    tri_cause_sig = Triangulation(
        xs.ravel() - 0.5,
        ys.ravel() - 0.5,
        triangles=np.array(triangles2),
        mask=tri_cause_sig_mask,
    )

    triang1 = Triangulation(xs.ravel() - 0.5, ys.ravel() - 0.5, triangles1)
    triang2 = Triangulation(xs.ravel() - 0.5, ys.ravel() - 0.5, triangles2)

    tri_color = "#1f0954"

    img1 = ax.tripcolor(
        triang1,
        np.array(ci_gwas_ace_mat).ravel(),
        cmap="RdBu",
        # vmin=-vm,
        # vmax=vm,
        norm=mpl.colors.SymLogNorm(vmin=-1.0, vmax=1.0, linthresh=0.01),
    )
    img2 = ax.tripcolor(
        triang2,
        np.array(cause_ace_mat).ravel(),
        cmap="RdBu",
        # vmin=-vm,
        # vmax=vm,
        norm=mpl.colors.SymLogNorm(vmin=-1.0, vmax=1.0, linthresh=0.01),
    )
    _ = ax.tripcolor(
        tri_cause_sig,
        np.array(cause_ace_mat).ravel(),
        cmap="RdBu",
        # vmin=-vm,
        # vmax=vm,
        edgecolor="#df4ec8",  # purpleish pink
        linewidth=4,
        norm=mpl.colors.SymLogNorm(vmin=-1.0, vmax=1.0, linthresh=0.01),
    )
    _ = ax.tripcolor(
        tri_ci_gwas_link,
        np.array(ci_gwas_ace_mat).ravel(),
        # cmap=mpl.colors.ListedColormap(["w", "#fcb001"]),
        cmap="RdBu",
        # vmin=-vm,
        # vmax=vm,
        edgecolor="#019529",  # irish green
        linewidth=4,
        linestyle=":",
        norm=mpl.colors.SymLogNorm(vmin=-1.0, vmax=1.0, linthresh=0.01),
    )

    # plt.colorbar(img2, ticks=range(10), pad=-0.05)
    # ax.figure.colorbar(img2, cax=ax, **cbar_kw)
    cbar_kw = {"fraction": 0.1, "pad": 0.04}
    cbar = plt.colorbar(img1, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(M))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(N), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)
    ax.set_xticks(np.arange(M + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="gray", linestyle="--", linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_xlabel("disease")
    ax.set_ylabel("risk factor")

    ci_gwas_patch = mpatches.Patch(color="#fcb001", label="CI-GWAS")
    cause_patch = mpatches.Patch(color="#448ee4", label="CAUSE")

    # ax.legend(
    #     handles=[ci_gwas_patch, cause_patch], bbox_to_anchor=(1, 0.5), loc="lower left"
    # )

    if title:
        ax.set_title(title, **title_kw)

    return img1, img2


def plot_direct_link_cause_comparison_wide(
    pag_path: str,
    pheno_path: str,
    p_thr=0.05,
    title=None,
    title_kw=dict(),
    ax=None,
    max_path_len=np.inf,
):
    if ax is None:
        plt.figure(figsize=(3, 10))
        ax = plt.gca()

    disease_list = sorted(list(diseases))

    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)

    causal_paths = get_causal_paths(pag_path, pheno_path, max_path_len=max_path_len)

    ci_gwas_links = {}
    cause_links = {}

    cause_ys = set(cause_gamma["y1"].values)
    cause_ys.update(set(cause_gamma["y2"].values))
    reg_pnames = cause_ys.intersection(set(p_names))

    for i in range(num_phen):
        for j in range(num_phen):
            name_i = p_names[i]
            name_j = p_names[j]
            if name_i in risk_factors and name_j in diseases:
                cp = cause_gamma[
                    (cause_gamma["y1"] == name_i) & (cause_gamma["y2"] == name_j)
                ].p_value.values
                if causal_paths[i, j] == 1:
                    ci_gwas_links[(name_i, name_j)] = 1
                if len(cp) >= 1 and cp[0] <= p_thr:
                    cause_links[(name_i, name_j)] = 1

    rf_intersection = sorted(list(risk_factors.intersection(reg_pnames)))
    rf_cg_only = sorted(list(risk_factors - set(rf_intersection)))

    risk_factor_list = rf_intersection + rf_cg_only

    ci_gwas_links_mat = [
        [0 for _ in range(len(diseases))] for _ in range(len(risk_factors))
    ]
    cause_links_mat = [
        [0 for _ in range(len(diseases))] for _ in range(len(risk_factors))
    ]

    for i, rf in enumerate(risk_factor_list):
        for j, ds in enumerate(disease_list):
            if (rf, ds) in ci_gwas_links:
                ci_gwas_links_mat[i][j] = 1
            if (rf, ds) in cause_links:
                cause_links_mat[i][j] = 1

    row_labels = [r"$\bf{{{0}}}$".format(e) for e in disease_list]
    col_labels = [r"$\bf{{{0}}}$".format(e) for e in rf_intersection] + rf_cg_only

    N = len(diseases)
    M = len(risk_factors)
    x = np.arange(M + 1)
    y = np.arange(N + 1)
    xs, ys = np.meshgrid(x, y)

    triangles1 = [
        (i + j * (M + 1), i + 1 + j * (M + 1), i + (j + 1) * (M + 1))
        for j in range(N)
        for i in range(M)
    ]
    triangles2 = [
        (i + 1 + j * (M + 1), i + 1 + (j + 1) * (M + 1), i + (j + 1) * (M + 1))
        for j in range(N)
        for i in range(M)
    ]
    triang1 = Triangulation(xs.ravel() - 0.5, ys.ravel() - 0.5, triangles1)
    triang2 = Triangulation(xs.ravel() - 0.5, ys.ravel() - 0.5, triangles2)
    img1 = ax.tripcolor(
        triang1,
        np.array(ci_gwas_links_mat).T.ravel(),
        cmap=mpl.colors.ListedColormap(["w", "#fcb001"]),
    )
    img2 = ax.tripcolor(
        triang2,
        np.array(cause_links_mat).T.ravel(),
        cmap=mpl.colors.ListedColormap(["w", "#448ee4"]),
    )

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(M))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(N), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)
    ax.set_xticks(np.arange(M + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="gray", linestyle="--", linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_ylabel("disease")
    ax.set_xlabel("risk factor")

    ci_gwas_patch = mpatches.Patch(color="#fcb001", label="CI-GWAS")
    cause_patch = mpatches.Patch(color="#448ee4", label="CAUSE")

    ax.legend(
        handles=[ci_gwas_patch, cause_patch], bbox_to_anchor=(1, 0.5), loc="lower left"
    )

    if title:
        ax.set_title(title, **title_kw)

    return img1, img2


def plot_direct_link_cause_comparison(
    pag_path: str,
    pheno_path: str,
    p_thr=0.05,
    title=None,
    title_kw=dict(),
    ax=None,
):
    if ax is None:
        plt.figure(figsize=(3, 10))
        ax = plt.gca()

    disease_list = sorted(list(diseases))

    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)

    causal_paths = get_causal_paths(pag_path, pheno_path)

    ci_gwas_links = {}
    cause_links = {}

    cause_ys = set(cause_gamma["y1"].values)
    cause_ys.update(set(cause_gamma["y2"].values))
    reg_pnames = cause_ys.intersection(set(p_names))

    for i in range(num_phen):
        for j in range(num_phen):
            name_i = p_names[i]
            name_j = p_names[j]
            if name_i in risk_factors and name_j in diseases:
                cp = cause_gamma[
                    (cause_gamma["y1"] == name_i) & (cause_gamma["y2"] == name_j)
                ].p_value.values
                if causal_paths[i, j] == 1:
                    ci_gwas_links[(name_i, name_j)] = 1
                if len(cp) >= 1 and cp[0] <= p_thr:
                    cause_links[(name_i, name_j)] = 1

    rf_intersection = sorted(list(risk_factors.intersection(reg_pnames)))
    rf_cg_only = sorted(list(risk_factors - set(rf_intersection)))

    risk_factor_list = rf_intersection + rf_cg_only

    ci_gwas_links_mat = [
        [0 for _ in range(len(diseases))] for _ in range(len(risk_factors))
    ]
    cause_links_mat = [
        [0 for _ in range(len(diseases))] for _ in range(len(risk_factors))
    ]

    for i, rf in enumerate(risk_factor_list):
        for j, ds in enumerate(disease_list):
            if (rf, ds) in ci_gwas_links:
                ci_gwas_links_mat[i][j] = 1
            if (rf, ds) in cause_links:
                cause_links_mat[i][j] = 1

    col_labels = [r"$\bf{{{0}}}$".format(e) for e in disease_list]
    row_labels = [r"$\bf{{{0}}}$".format(e) for e in rf_intersection] + rf_cg_only

    M = len(diseases)
    N = len(risk_factors)
    x = np.arange(M + 1)
    y = np.arange(N + 1)
    xs, ys = np.meshgrid(x, y)

    triangles1 = [
        (i + j * (M + 1), i + 1 + j * (M + 1), i + (j + 1) * (M + 1))
        for j in range(N)
        for i in range(M)
    ]
    triangles2 = [
        (i + 1 + j * (M + 1), i + 1 + (j + 1) * (M + 1), i + (j + 1) * (M + 1))
        for j in range(N)
        for i in range(M)
    ]
    triang1 = Triangulation(xs.ravel() - 0.5, ys.ravel() - 0.5, triangles1)
    triang2 = Triangulation(xs.ravel() - 0.5, ys.ravel() - 0.5, triangles2)
    img1 = ax.tripcolor(
        triang1,
        np.array(ci_gwas_links_mat).ravel(),
        cmap=mpl.colors.ListedColormap(["w", "#fcb001"]),
    )
    img2 = ax.tripcolor(
        triang2,
        np.array(cause_links_mat).ravel(),
        cmap=mpl.colors.ListedColormap(["w", "#448ee4"]),
    )

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(M))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(N), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)
    ax.set_xticks(np.arange(M + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="gray", linestyle="--", linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_ylabel("risk factor")
    ax.set_xlabel("disease")

    ci_gwas_patch = mpatches.Patch(color="#fcb001", label="CI-GWAS")
    cause_patch = mpatches.Patch(color="#448ee4", label="CAUSE")

    ax.legend(
        handles=[ci_gwas_patch, cause_patch], bbox_to_anchor=(1, 0.5), loc="lower left"
    )

    if title:
        ax.set_title(title, **title_kw)

    return img1, img2


def marker_pheno_associations_with_pnames(
    blockfile: str,
    outdir: str,
    p_names: list[str],
    bim_path: str,
    depth=1,
):
    bim_df = pd.read_csv(bim_path, sep="\t", header=None)

    rs_ids = bim_df[1].values
    chrs = bim_df[0].values
    on_chr_positions = bim_df[3].values

    num_phen = len(p_names)
    assoc_markers = []

    anc_sets = global_ancestor_sets(
        blockfile, outdir, reduced_indices=False, depth=depth
    )

    for pix in np.arange(0, num_phen) + BASE_INDEX:
        for bim_line in anc_sets[pix]:
            try:
                assoc_markers.append(
                    {
                        "phenotype": p_names[pix - BASE_INDEX],
                        "rsID": rs_ids[bim_line],
                        "bim_line_ix": bim_line,
                        "chr": chrs[bim_line],
                        "bp": on_chr_positions[bim_line],
                    }
                )
            except IndexError:
                print("pix: ", pix, "bim_line: ", bim_line)

    return pd.DataFrame(assoc_markers)


def marker_pheno_associations(
    blockfile: str,
    outdir: str,
    bim_path: str,
    corr_path: str,
    adj_path: str,
    num_phen=None,
    pheno_path=None,
):
    if num_phen is None and pheno_path is None:
        raise RuntimeError("Either num_phen or pheno_path have to specified")

    if pheno_path is None:
        p_names = list(range(1, num_phen + 1))
    else:
        p_names = get_pheno_codes(pheno_path)
        num_phen = len(p_names)

    bim_df = pd.read_csv(bim_path, sep="\t", header=None)

    rs_ids = bim_df[1].values
    chrs = bim_df[0].values
    on_chr_positions = bim_df[3].values

    assoc_markers = []

    adj = mmread(adj_path).toarray()
    corr = mmread(corr_path).toarray()

    gr = merge_block_outputs(blockfile, outdir)
    glob_ixs = np.array(sorted(list(gr.gmi.values())))

    for pix in np.arange(0, num_phen):
        bim_lines = glob_ixs[np.where(adj[pix, num_phen:])]
        corrs = corr[pix, num_phen:][np.where(adj[pix, num_phen:])]
        for bim_line, c in zip(bim_lines, corrs):
            assoc_markers.append(
                {
                    "phenotype": p_names[pix],
                    "rsID": rs_ids[bim_line],
                    "bim_line_ix": bim_line,
                    "chr": chrs[bim_line],
                    "bp": on_chr_positions[bim_line],
                    "corr": c,
                }
            )

    return pd.DataFrame(assoc_markers)


# def marker_pheno_associations(
#     blockfile: str,
#     outdir: str,
#     pheno_path: str,
#     bim_path: str,
#     depth=1,
# ):
#     bim_df = pd.read_csv(bim_path, sep="\t", header=None)

#     rs_ids = bim_df[1].values
#     chrs = bim_df[0].values
#     on_chr_positions = bim_df[3].values

#     p_names = get_pheno_codes(pheno_path)
#     num_phen = len(p_names)
#     assoc_markers = []

#     anc_sets = global_ancestor_sets(
#         blockfile, outdir, reduced_indices=False, depth=depth
#     )

#     for pix in np.arange(0, num_phen) + BASE_INDEX:
#         for bim_line in anc_sets[pix]:
#             try:
#                 assoc_markers.append(
#                     {
#                         "phenotype": p_names[pix - BASE_INDEX],
#                         "rsID": rs_ids[bim_line],
#                         "bim_line_ix": bim_line,
#                         "chr": chrs[bim_line],
#                         "bp": on_chr_positions[bim_line],
#                     }
#                 )
#             except IndexError:
#                 print("pix: ", pix, "bim_line: ", bim_line)

#     return pd.DataFrame(assoc_markers)


def pag_edge_types(pag_path: str, pheno_path: str) -> dict[tuple[int, int], int]:
    all_edges = {}
    p_names = get_pheno_codes(pheno_path)
    pag = mmread(pag_path).tocsr()
    plr = pag.tolil().rows
    for j, r in enumerate(plr):
        for i in r:
            e = (pag[i, j], pag[j, i])
            if e not in all_edges:
                all_edges[e] = 0
            all_edges[e] += 1
    return all_edges


def pag_x_to_y_edge_types(pag_path: str, pheno_path: str) -> dict[tuple[int, int], int]:
    x_to_y = {}
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    pag = mmread(pag_path).tocsr()
    plr = pag.tolil().rows

    for j in range(num_phen):
        for i in plr[j]:
            if i >= num_phen:
                e = (pag[i, j], pag[j, i])
                if e not in x_to_y:
                    x_to_y[e] = 0
                x_to_y[e] += 1
    return x_to_y


def combine_all_pheno_and_plot():
    outdir = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/bdpc_d1_l6_a1e8/"
    blockfile = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.blocks"
    p_names = get_pheno_codes(
        "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/input.phen"
    )
    num_phen = len(p_names)

    gepm = global_epm(blockfile, outdir)
    geps = global_eps(blockfile, outdir)
    gr = merge_block_outputs(blockfile, outdir)
    # gr.write_mm(outdir + "all_merged")

    geps_global_ix = {}
    for k, v in geps.items():
        if k[0] >= k[1]:
            nk = f"{p_names[k[0] - 1]}_{p_names[k[1] - 1]}"
            geps_global_ix[nk] = [int(gr.gmi[ix]) for ix in v]

    # with open("./../../../../human_data/causality/parent_set_selection/bdpc_d1_l6_a1e8/bim_indices_of_pleiotropy_matrix_set_elements.json", 'w') as fout:
    #     json.dump(geps_global_ix, fout)

    np_gepm = np.zeros(shape=(num_phen, num_phen))
    pd_gepm = []
    for (i, j), c in gepm.items():
        np_gepm[i - 1, j - 1] = c
        pd_gepm.append({"y1": p_names[i - 1], "y2": p_names[j - 1], "count": c})
    pd_gepm = pd.DataFrame(pd_gepm)

    pd_pcm = []
    for i in range(1, num_phen + 1):
        for j in range(1, num_phen + 1):
            pd_pcm.append(
                {"y1": p_names[i - 1], "y2": p_names[j - 1], "v": gr.scm[(i, j)]}
            )
    pd_pcm = pd.DataFrame(pd_pcm)

    plt.figure(figsize=(20, 15))
    g = pd_pcm.pivot("y1", "y2", "v")
    mask = np.triu(np.ones_like(g, dtype=bool))
    np.fill_diagonal(mask, False)
    sns.heatmap(
        g,
        xticklabels=1,
        yticklabels=1,
        cmap="RdBu",
        annot=True,
        square=True,
        mask=mask,
        vmax=1,
        vmin=-1,
    )

    pd_pam = []
    for i in range(1, num_phen + 1):
        for j in range(1, num_phen + 1):
            if (i, j) in gr.sam:
                pd_pam.append({"y1": p_names[i - 1], "y2": p_names[j - 1], "v": 1})
            else:
                pd_pam.append({"y1": p_names[i - 1], "y2": p_names[j - 1], "v": 0})
    pd_pam = pd.DataFrame(pd_pam)

    plt.figure(figsize=(14, 14))
    g = pd_pam.pivot("y1", "y2", "v")
    mask = np.triu(np.ones_like(g, dtype=bool))
    np.fill_diagonal(mask, False)
    sns.heatmap(
        g,
        xticklabels=1,
        yticklabels=1,
        cbar=False,
        square=True,
        mask=mask,
        cmap="Blues",
        vmax=2,
    )

    num_phen = len(p_names)
    np_gepm = np.zeros(shape=(num_phen, num_phen))
    pd_gepm = []
    for (i, j), c in {k: len(v) for k, v in geps.items()}.items():
        np_gepm[i - 1, j - 1] = c
        pd_gepm.append({"y1": p_names[i - 1], "y2": p_names[j - 1], "count": c})
    pd_gepm = pd.DataFrame(pd_gepm)

    plt.figure(figsize=(20, 15))
    g = pd_gepm.pivot("y1", "y2", "count")
    mask = np.triu(np.ones_like(g, dtype=bool))
    np.fill_diagonal(mask, False)
    sns.heatmap(
        g,
        xticklabels=1,
        yticklabels=1,
        norm=matplotlib.colors.LogNorm(),
        cmap="Blues",
        annot=True,
        square=True,
        mask=mask,
        cbar_kws={"label": "# parent markers"},
    )

    pag_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/bdpc_d1_l6_a1e8/all_merged_pag.mtx"

    pag_coo = mmread(pag_path)

    pag = mmread(pag_path).tocsr()

    pd_pag = []

    for i in range(num_phen):
        for j in range(num_phen):
            v = 0
            if pag[i, j] == 2:
                v = 1
            elif pag[i, j] == 1:
                v = -1
            pd_pag.append({"y1": p_names[i], "y2": p_names[j], "edge": v})

    pd_pag = pd.DataFrame(pd_pag)

    plt.figure(figsize=(20, 15))
    g = pd_pag.pivot("y1", "y2", "edge")
    sns.heatmap(
        g,
        xticklabels=1,
        yticklabels=1,
        cmap="Blues",
        annot=False,
        square=True,
        cbar=False,
    )

    num_phen = len(p_names)
    geps = pag_exclusive_pleiotropy_sets(pag, num_phen)

    np_gepm = np.zeros(shape=(num_phen, num_phen))
    pd_gepm = []
    for (i, j), c in {k: len(v) for k, v in geps.items()}.items():
        np_gepm[i, j] = c
        pd_gepm.append({"y1": p_names[i], "y2": p_names[j], "count": c})
    pd_gepm = pd.DataFrame(pd_gepm)

    plt.figure(figsize=(20, 15))
    g = pd_gepm.pivot("y1", "y2", "count")
    mask = np.triu(np.ones_like(g, dtype=bool))
    np.fill_diagonal(mask, False)
    sns.heatmap(
        g,
        xticklabels=1,
        yticklabels=1,
        norm=matplotlib.colors.LogNorm(),
        cmap="Blues",
        annot=True,
        square=True,
        mask=mask,
        cbar_kws={"label": "# parent markers"},
    )


@dataclass
class Performance:
    mse: float
    bias: float
    var: float
    mse_tp: float
    bias_tp: float
    var_tp: float
    fdr: float
    tpr: float


@dataclass
class MR:
    method: str
    exposure: str
    outcome: str
    p: str
    estimate: str


mr_cn = [
    MR("ivw", "Exposure", "Outcome", "p", "est"),
    MR("egger", "Exposure", "Outcome", "p", "est"),
    MR(
        "mrpresso",
        "V1",
        "V2",
        "P-value",
        "Causal Estimate",
    ),
    MR("cause", "V1", "V2", "V4", "gamma"),
]


def make_link_type_dict(adj: np.array):
    n = adj.shape[0]
    link_types = {}
    for i in range(n):
        for j in range(i + 1, n):
            if (adj[i, j] != 0) and (adj[j, i] != 0):
                link_types[(i, j)] = (2, 2)
            elif (adj[i, j] != 0) and (adj[j, i] == 0):
                link_types[(i, j)] = (2, 3)
            elif (adj[i, j] == 0) and (adj[j, i] != 0):
                link_types[(i, j)] = (3, 2)
    return link_types


def make_adj_symmetric(adj: np.array):
    n = adj.shape[0]
    sym_adj = np.zeros_like(adj)
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j] != 0 or adj[j, i] != 0:
                sym_adj[i, j] = 1
                sym_adj[j, i] = 1
    return sym_adj


def pag_to_dag_directed(pag: np.array):
    n = pag.shape[0]
    dag = np.zeros_like(pag)
    for i in range(n):
        for j in range(n):
            if pag[i, j] == 2 and pag[j, i] == 3:
                dag[i, j] = 1
            elif pag[i, j] == 2 and pag[j, i] == 2:
                dag[i, j] = 1
                dag[j, i] = 1
    return dag


def pag_to_dag_possibly_directed(pag: np.array):
    n = pag.shape[0]
    dag = np.zeros_like(pag)
    for i in range(n):
        for j in range(n):
            if (pag[i, j] == 2 and pag[j, i] == 3) or (
                pag[i, j] == 2 and pag[j, i] == 1
            ):
                dag[i, j] = 1
            elif pag[i, j] == 2 and pag[j, i] == 2:
                dag[i, j] = 1
                dag[j, i] = 1
    return dag


def path_in_sem(adj: np.array):
    num_var = adj.shape[0]
    causal_paths = np.zeros(shape=(num_var, num_var))
    for start_node in range(num_var):
        q = queue.Queue()
        q.put(start_node)
        descendants = set()
        while not q.empty():
            current_node = q.get()
            # assuming that links are always
            # ordered towards variables with larger index
            for neighbor in range(current_node + 1, num_var):
                link = adj[current_node, neighbor] != 0
                if neighbor not in descendants and link:
                    descendants.add(neighbor)
                    q.put(neighbor)
        for j in descendants:
            causal_paths[start_node, j] = 1
    return causal_paths


@dataclass
class CiGwasRelativeOrientationPerformance:
    mr_pos_tpr: float
    mr_neg_tdr: float


def compare_ci_gwas_orientation_performance_to_mr(
    true_directed: np.array,
    true_bidirected: np.array,
    mr_links: np.array,
    cig_pag: np.array,
) -> CiGwasRelativeOrientationPerformance:
    """Compute fraction of shared correctly directed links between ci-gwas and mr

    Args:
        true_directed (np.array): true DAG
        true_bidirected (np.array): edges that should be bidirected in the true MAG
        mr_links (np.array): directed adjacencies from MR run
        cig_pag (np.array): PAG from ci-gwas run
    """
    sum_mr = 0
    sum_shared = 0
    sum_cig_only = 0
    sum_cig_only_correct = 0
    num_phen = true_directed.shape[0]
    for i in range(num_phen):
        for j in range(i, num_phen):
            mr_edge = (mr_links[i, j], mr_links[j, i])
            cig_edge = (cig_pag[i, j], cig_pag[j, i])
            mr_true = (mr_edge == (True, False) and true_directed[i, j]) or (
                mr_edge == (True, True) and true_bidirected[i, j]
            )
            cig_true = (cig_edge in [(2, 3), (2, 1)] and true_directed[i, j]) or (
                cig_edge == (2, 2) and true_bidirected[i, j]
            )
            if mr_true:
                sum_mr += 1
                if cig_true:
                    sum_shared += 1
            if cig_edge != (0, 0) and not mr_true:
                sum_cig_only += 1
                if cig_true:
                    sum_cig_only_correct += 1

    return CiGwasRelativeOrientationPerformance(
        sum_shared / sum_mr if sum_mr > 0 else np.nan,
        sum_cig_only_correct / sum_cig_only if sum_cig_only > 0 else np.nan,
    )


@dataclass
class OrientationPerformance:
    directed: float
    bidirected: float


def calculate_pxp_orientation_performance_mr(
    true_directed: np.array,
    true_bidirected: np.array,
    inferred_links: np.array,
):
    causal_paths = path_in_sem(true_directed)
    num_phen = true_directed.shape[0]
    inferred_uni = 0
    inferred_bi = 0
    correct_uni = 0
    correct_bi = 0
    for i in range(num_phen):
        for j in range(i, num_phen):
            inferred_edge = (inferred_links[i, j], inferred_links[j, i])
            if inferred_edge == (True, True):
                inferred_bi += 1
                correct_bi += true_bidirected[i, j]
            elif inferred_edge != (False, False):
                inferred_uni += 1
                if inferred_edge == (True, False) and causal_paths[i, j]:
                    correct_uni += 1
    return OrientationPerformance(
        correct_uni / inferred_uni if inferred_uni > 0 else np.nan,
        correct_bi / inferred_bi if inferred_bi > 0 else np.nan,
    )


def calculate_pxp_orientation_performance_ci_gwas(
    true_directed: np.array,
    true_bidirected: np.array,
    inferred_pag: np.array,
) -> OrientationPerformance:
    num_phen = true_directed.shape[0]
    inferred_uni = 0
    inferred_bi = 0
    correct_uni = 0
    correct_bi = 0
    for i in range(num_phen):
        for j in range(i, num_phen):
            inferred_edge = (inferred_pag[i, j], inferred_pag[j, i])
            if inferred_edge == (2, 2):
                inferred_bi += 1
                correct_bi += true_bidirected[i, j]
            elif inferred_edge != (0, 0):
                inferred_uni += 1
                if inferred_edge in [(2, 3), (2, 1)] and true_directed[i, j]:
                    correct_uni += 1
    return OrientationPerformance(
        correct_uni / inferred_uni if inferred_uni > 0 else np.nan,
        correct_bi / inferred_bi if inferred_bi > 0 else np.nan,
    )


def calulate_performance_metrics(
    true_adj: np.array, true_eff: np.array, est_adj: np.array, est_eff: np.array
):
    sym_est_adj = make_adj_symmetric(est_adj)
    true_adj_masked = np.ma.array(true_adj, mask=np.tri(true_adj.shape[0], k=0))
    est_adj_masked = np.ma.array(est_adj, mask=np.tri(true_adj.shape[0], k=0))
    sym_est_adj_masked = np.ma.array(sym_est_adj, mask=np.tri(true_adj.shape[0], k=0))
    est_eff_masked = np.ma.array(est_eff, mask=np.tri(true_adj.shape[0], k=0))
    true_eff_masked = np.ma.array(true_eff, mask=np.tri(true_adj.shape[0], k=0))

    # no need to subset here. for MR, we do it before
    sig_eff = est_eff
    tp_mat = (true_adj_masked != 0) & (est_adj_masked != 0)

    p = np.sum(true_adj_masked != 0)
    # f = np.sum(true_adj_masked == 0)
    tp = np.sum((true_adj_masked != 0) & (sym_est_adj_masked != 0))
    fp = np.sum((true_adj_masked == 0) & (sym_est_adj_masked != 0))
    # fn = np.sum((true_adj_masked != 0) & (sym_est_adj_masked == 0))
    mse = np.sum((true_eff - sig_eff) ** 2)
    bias = np.sum(true_eff - sig_eff)
    var = np.var(true_eff - sig_eff)
    mse_tp = np.sum((est_eff_masked[tp_mat] - true_eff_masked[tp_mat]) ** 2)
    bias_tp = np.sum(est_eff_masked[tp_mat] - true_eff_masked[tp_mat])
    var_tp = np.var(est_eff_masked[tp_mat] - true_eff_masked[tp_mat])
    if (fp + tp) == 0:
        fdr = np.nan
    else:
        fdr = fp / (fp + tp)
    tpr = tp / p
    if np.ma.is_masked(mse_tp):
        mse_tp = np.nan
    if np.ma.is_masked(var_tp):
        var_tp = np.nan
    if np.ma.is_masked(bias_tp):
        bias_tp = np.nan
    return Performance(mse, bias, var, mse_tp, bias_tp, var_tp, fdr, tpr)


def load_simulation_results() -> pd.DataFrame:
    pdir = f"/nfs/scistore17/robingrp/human_data/causality/bias_as_fn_of_alpha/sim_small_effects/"

    d = 1
    l = 6
    n_arr = [2000, 4000, 8000, 16000]
    m_arr = [200, 400, 800, 1600]
    # m_arr = [200, 400, 1600]
    e_arr = list(range(1, 9))
    rep_arr = list(range(1, 21))
    num_phen = 10

    rows = []

    for n, m, rep in itertools.product(n_arr, m_arr, rep_arr):
        # load true dag
        dag_path = pdir + f"./true_adj_mat_n{n}_SNP_{m}_it_{rep}.mtx"
        dag = mmread(dag_path).tocsr()
        pdag = dag[-(num_phen):, -(num_phen):].toarray()

        # load true causal effects
        eff_path = pdir + f"./True_causaleffect_n{n}_SNP_{m}_it_{rep}.mtx"
        eff = mmread(eff_path).tocsr()
        peff = eff[-(num_phen):, -(num_phen):].toarray()
        peff = np.triu(peff, k=1)

        # mr standalone
        for mr in mr_cn:
            mr_results = pd.read_csv(
                pdir + f"mr_res_{mr.method}_n{n}_SNP_{m}_it_{rep}.csv"
            )
            mr_results["i"] = mr_results[mr.exposure].apply(
                lambda x: int(x.split("Y")[1]) - 1
            )
            mr_results["j"] = mr_results[mr.outcome].apply(
                lambda x: int(x.split("Y")[1]) - 1
            )
            pvals = np.ones(shape=(num_phen, num_phen))
            pvals[mr_results["i"], mr_results["j"]] = mr_results[mr.p]
            effects = np.zeros(shape=(num_phen, num_phen))
            effects[mr_results["i"], mr_results["j"]] = mr_results[mr.estimate]
            for e in e_arr:
                adj = pvals < 10 ** (-e)
                eff_copy = np.copy(effects)
                eff_copy[~adj] = 0
                perf = calulate_performance_metrics(pdag, peff, adj, eff_copy)
                rows.append(
                    {
                        # "edge orientation": perf.correct_orientation,
                        "mse": perf.mse,
                        "var": perf.var,
                        "bias": perf.bias,
                        "mse_tp": perf.mse_tp,
                        "var_tp": perf.var_tp,
                        "bias_tp": perf.bias_tp,
                        # "fdr": perf.fdr,
                        "tpr": perf.tpr,
                        "n": n,
                        "m": m,
                        "rep": rep,
                        "alpha": 10 ** (-e),
                        "method": mr.method,
                    }
                )

        for e in e_arr:
            indir = pdir + f"simpc_d{d}_l{l}_e{e}_i{rep}_n{n}_m{m}/"

            mdim_path = indir + "skeleton.mdim"
            with open(mdim_path, "r") as fin:
                num_var, num_phen, max_level = [
                    int(elem) for elem in fin.readline().split()
                ]

            # load skeleton
            adj = np.fromfile(indir + "skeleton.adj", dtype=np.int32).reshape(
                num_var, num_var
            )
            padj = adj[-(num_phen):, -(num_phen):]

            # load pag
            pag_path = indir + "estimated_pag_mpu.mtx"
            est_pag = mmread(pag_path).tocsr()[-(num_phen):, -(num_phen):].toarray()
            est_dag = pag_to_dag_directed(est_pag)

            # load var indices
            # var_ixs = np.fromfile(indir + "skeleton.ixs", dtype=np.int32)

            # load aces
            pace = np.zeros(shape=(num_phen, num_phen))
            missing = []
            for i in range(1, num_phen + 1):
                for j in range(1, num_phen + 1):
                    if i == j:
                        continue
                    file = indir + f"estimated_causaleffect_i{i}_j{j}_mpu.csv"
                    try:
                        with open(file) as fin:
                            sym = fin.readline().strip()
                            if sym == "NaN":
                                pace[i - 1, j - 1] = 0.0
                            else:
                                pace[i - 1, j - 1] = eval(sym)
                    except FileNotFoundError as e:
                        missing.append((i, j))

            ixs_path = indir + "skeleton.ixs"
            var_ixs = np.fromfile(ixs_path, dtype=np.int32)

            # load true dag
            dag_path = pdir + f"./true_adj_mat_n{n}_SNP_{m}_it_{rep}.mtx"
            full_dag = mmread(dag_path).tocsr()

            # load pag
            full_pag = mmread(pag_path).tocsr()

            dag_mxp = full_dag.toarray()[:m, -num_phen:]
            pag_mxp_sub = full_pag[:, -num_phen:].toarray()
            pag_mxp = np.zeros_like(dag_mxp)
            for sub_ix, glob_ix in enumerate(var_ixs):
                if glob_ix < m:
                    pag_mxp[glob_ix, :] = pag_mxp_sub[sub_ix, :]

            mxp_p = np.sum(dag_mxp != 0)
            mxp_tp = np.sum((dag_mxp != 0) & (pag_mxp != 0))
            mxp_fp = np.sum((dag_mxp == 0) & (pag_mxp != 0))
            mxp_fdr = mxp_fp / (mxp_fp + mxp_tp)
            mxp_tpr = mxp_tp / mxp_p

            perf = calulate_performance_metrics(pdag, peff, padj, pace)

            # perf = calulate_performance_metrics(pdag, peff, est_dag, pace)
            rows.append(
                {
                    "marker-trait fdr": mxp_fdr,
                    "marker-trait tpr": mxp_tpr,
                    # "edge orientation": perf.correct_orientation,
                    "mse": perf.mse,
                    "var": perf.var,
                    "bias": perf.bias,
                    "mse_tp": perf.mse_tp,
                    "var_tp": perf.var_tp,
                    "bias_tp": perf.bias_tp,
                    "fdr": perf.fdr,
                    "tpr": perf.tpr,
                    "n": n,
                    "m": m,
                    "rep": rep,
                    "alpha": 10 ** (-e),
                    "method": "ci-gwas",
                }
            )

    return pd.DataFrame(rows)


@dataclass
class SimulationTruth:
    dag_mxp: pd.DataFrame
    dag_pxp: np.array
    bidirected: np.array
    effects: np.array


def load_real_data_simulation_truth(
    rep: int,
    num_p=10,
    num_m=600,
    wdir="/nfs/scistore17/robingrp/mrobinso/cigwas/ukb/sim/real_marker_sim/sim/",
):
    simdir = wdir + f"sim_{rep}/"
    true_snp_ids = []
    with open(simdir + "CV_snp_ids.txt", "r") as fin:
        for line in fin:
            true_snp_ids.append(line.strip())
    true_eff_file = glob(simdir + "true_causaleffects*")[0]
    true_eff_lines = []
    true_eff_lines_all = []
    with open(true_eff_file, "r") as fin:
        for line in fin:
            # we skip the latent vars here
            true_eff_lines.append(
                {k + 1: v for (k, v) in enumerate(line.strip().split()[2:])}
            )
            true_eff_lines_all.append(
                {k: v for (k, v) in enumerate(line.strip().split())}
            )
    true_dag_mxp = pd.DataFrame(true_eff_lines[:num_m], index=true_snp_ids, dtype=float)
    true_dag_pxp = np.triu(
        np.array(
            pd.DataFrame(
                true_eff_lines[-num_p:],
                index=list(range(1, num_p + 1)),
                dtype=float,
            )
        ),
        1,
    )
    true_dag_lpxlp = np.triu(
        np.array(pd.DataFrame(true_eff_lines_all[num_m:], dtype=float)), 1
    )
    # mark true 2-2 edges
    m = true_dag_lpxlp != 0
    true_bidirected = np.zeros(shape=(num_p, num_p))
    for i in range(2, num_p + 2):
        for j in range(i + 1, num_p + 2):
            if not m[i, j] and ((m[0, i] and m[0, j]) or (m[1, i] and m[1, j])):
                true_bidirected[i - 2, j - 2] = 1

    m = true_dag_pxp
    m = np.linalg.inv(np.eye(m.shape[0]) - m.T)
    true_eff = np.triu(m @ m.T, 1)

    return SimulationTruth(true_dag_mxp, true_dag_pxp, true_bidirected, true_eff)


def load_real_data_simulation_adj_performance(
    e_arr=[2, 4, 6, 8], rep_arr=list(range(1, 41)), num_p=10, max_level=3
) -> pd.DataFrame:
    """Load simulation results for ci-gwas and mr methods, calculate fdr, tpr for adjacencies."""
    rows = []
    wdir = "/nfs/scistore17/robingrp/mrobinso/cigwas/ukb/sim/real_marker_sim/"
    blockfile = wdir + "ukb22828_UKB_EST_v3_ldp08_estonia_intersect_m11000_chr1.blocks"
    bim_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_a1_forced.bim"
    rs_ids = pd.read_csv(bim_path, sep="\t", header=None)[1].values
    for rep in rep_arr:
        truth = load_real_data_simulation_truth(rep)
        for alpha_e in e_arr:
            # -------------------- CI-GWAS -----------------------------
            est_dir = wdir + f"cusk/sim{rep}_e{alpha_e}_l{max_level}_2stepsk/"
            adj_path = est_dir + "all_merged_sam.mtx"
            try:
                est_adj = mmread(adj_path).tocsr()
            except FileNotFoundError as err:
                print(err)
                continue

            gr = merge_block_outputs(blockfile, est_dir)
            glob_ixs = np.array(sorted(list(gr.gmi.values())))
            est_adj_mxp = pd.DataFrame(
                est_adj[num_p:, :num_p].toarray(),
                columns=list(range(1, num_p + 1)),
                index=rs_ids[glob_ixs],
                dtype=int,
            )
            est_adj_pxp = est_adj[:num_p, :num_p].toarray()
            est_adj_pxp_triu = np.triu(est_adj_pxp, 1)

            p = (truth.dag_mxp != 0).sum().sum()
            td_mxp_matched = truth.dag_mxp.reindex_like(est_adj_mxp)
            tp = ((est_adj_mxp != 0) & (td_mxp_matched != 0)).sum().sum()
            fp = ((est_adj_mxp != 0) & (td_mxp_matched == 0)).sum().sum()
            mxp_tpr = tp / p
            mxp_fdr = fp / (tp + fp)

            m = (truth.bidirected != 0) | (truth.dag_pxp != 0)
            p = np.sum(m)
            tp = np.sum(m & (est_adj_pxp_triu != 0))
            fp = np.sum(~m & (est_adj_pxp_triu != 0))
            pxp_tpr = tp / p
            pxp_fdr = fp / (tp + fp)

            rows.append(
                {
                    "x -> y fdr": mxp_fdr,
                    "x -> y tpr": mxp_tpr,
                    "y -> y fdr": pxp_fdr,
                    "y -> y tpr": pxp_tpr,
                    "rep": rep,
                    "alpha": 10 ** (-alpha_e),
                    "method": "ci-gwas",
                }
            )

            # -------------------- MR -----------------------------
            for mr_str in ["cause", "presso", "mvpresso", "ivw", "mvivw"]:
                mr_p = np.ones((num_p, num_p))
                mr_est = np.zeros((num_p, num_p))
                for outcome in range(0, num_p):
                    try:
                        fpaths = glob(
                            est_dir
                            + f"mr_{mr_str}_alpha*_sim{rep}_outcome_y{outcome + 1}_seed1000"
                        )
                        if len(fpaths) == 0:
                            continue
                        fpath = fpaths[0]
                        mr_res_df = pd.read_csv(fpath)
                        exposures = [int(s[1:]) - 1 for s in mr_res_df["Exposure"]]
                        mr_p[exposures, outcome] = mr_res_df["p"].values
                        mr_est[exposures, outcome] = mr_res_df["est"].values
                    except FileNotFoundError:
                        continue
                mr_links = mr_p <= (0.05 / ((num_p - 1) * 2))
                mr_adj = make_adj_symmetric(mr_links)

                m = (truth.bidirected != 0) | (truth.dag_pxp != 0)
                p = np.sum(m)
                tp = np.sum(m & (mr_adj != 0))
                pxp_tpr = tp / p

                rows.append(
                    {
                        "y -> y tpr": pxp_tpr,
                        "rep": rep,
                        "alpha": 10 ** (-alpha_e),
                        "method": mr_str,
                    }
                )
    return pd.DataFrame(rows)


def plot_real_data_simulation_adj_performance(
    e_arr=[2, 4, 6, 8], rep_arr=list(range(1, 41)), num_p=10, max_level=3, fig_path=None
):
    bar_xlabel_rotation = 45

    df = load_real_data_simulation_adj_performance(e_arr, rep_arr, num_p, max_level)
    gr = df.groupby(["method", "alpha"])
    means = gr.mean()
    stds = gr.std()

    alphas = 10.0 ** -np.array(e_arr[::-1])
    methods = ["ci-gwas", "cause", "presso", "mvpresso", "ivw", "mvivw"]
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    def plot_ci_gwas_bars(x_vals, metric, means, stds, ax, title, axhline=False):
        ax.set_title(title, **title_kw)
        x_sorted = sorted(list(x_vals))
        x = np.arange(len(x_vals))  # the label locations
        width = 0.5  # the width of the bars
        multiplier = 0

        handles = []
        for method in ["ci-gwas"]:
            mu = means.loc[method][metric]
            sig = stds.loc[method][metric]
            loc_x = np.array(
                [x_sorted.index(e) for e in means.loc[method].index.values]
            )
            offset = width * multiplier
            try:
                bars = ax.bar(
                    loc_x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "e)" or title == "d)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=bar_xlabel_rotation,
            ha="right",
            rotation_mode="anchor",
        )
        if axhline:
            ax.axhline(0.05, linestyle="dotted", color="gray")
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    def plot_bars(x_vals, metric, means, stds, ax, title, methods):
        ax.set_title(title, **title_kw)
        x_sorted = sorted(list(x_vals))
        x = np.arange(len(x_vals))  # the label locations
        width = 0.10  # the width of the bars
        multiplier = 0

        handles = []
        for method in methods:
            mu = means.loc[method][metric]
            sig = stds.loc[method][metric]
            loc_x = np.array(
                [x_sorted.index(e) for e in means.loc[method].index.values]
            )
            offset = width * multiplier
            try:
                bars = ax.bar(
                    loc_x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "g)" or title == "h)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x + width, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=bar_xlabel_rotation,
            ha="right",
            rotation_mode="anchor",
        )
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    fig = plt.figure(
        # layout="tight",
        figsize=(9, 6)
    )
    ax_dict = fig.subplot_mosaic(
        """
        XiX
        abc
        ddd
        """,
        empty_sentinel="X",
        sharex=False,
        # set the height ratios between the rows
        height_ratios=[0.1, 0.4, 0.5],
    )

    ax_dict["d"].set_xlabel(r"$\alpha$")

    plot_ci_gwas_bars(
        alphas, "x -> y fdr", means, stds, ax_dict["a"], "a)", axhline=True
    )
    ax_dict["a"].set_ylabel(r"$x - y \ FDR$")
    plot_ci_gwas_bars(alphas, "x -> y tpr", means, stds, ax_dict["b"], "b)")
    ax_dict["b"].set_ylabel(r"$x - y \ TPR$")
    plot_ci_gwas_bars(
        alphas, "y -> y fdr", means, stds, ax_dict["c"], "c)", axhline=True
    )
    ax_dict["c"].set_ylabel(r"$y - y \ FDR$")
    h = plot_bars(
        alphas, "y -> y tpr", means, stds, ax_dict["d"], "d)", methods=methods
    )
    ax_dict["d"].set_ylabel(r"$y - y \ TPR$")
    ax_dict["i"].legend(
        handles=h,
        labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "MR-PRESSO (MV)", "IVW", "IVW (MV)"],
        # labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "IVW Regression"],
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=3,
        # title="method",
    )
    ax_dict["i"].axis("off")
    # _ = [ax_dict[i].tick_params(labelbottom=False) for i in "de"]
    # _ = [ax_dict[i].sharex(ax_dict['f']) for i in 'de']
    fig.subplots_adjust(wspace=0.7, hspace=1.3)
    if fig_path is not None:
        plt.savefig(fig_path, bbox_inches="tight")


def load_real_data_simulation_results(
    e_arr=[3, 4, 6, 8],
    rep_arr=list(range(1, 41)),
    num_p=10,
) -> pd.DataFrame:
    alpha_str_d = {
        1: "0.1",
        2: "0.01",
        3: "0.001",
        4: "1e-04",
        5: "1e-05",
        6: "1e-06",
        7: "1e-07",
        8: "1e-08",
    }
    num_m = 600
    alpha_str_arr = [alpha_str_d[e] for e in e_arr]
    rows = []
    wdir = "/nfs/scistore17/robingrp/mrobinso/cigwas/ukb/sim/real_marker_sim/"
    blockfile = wdir + "ukb22828_UKB_EST_v3_ldp08_estonia_intersect_m11000_chr1.blocks"
    bim_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_a1_forced.bim"
    rs_ids = pd.read_csv(bim_path, sep="\t", header=None)[1].values
    for rep in rep_arr:
        simdir = f"/nfs/scistore17/robingrp/mrobinso/cigwas/ukb/sim/real_marker_sim/sim/sim_{rep}/"
        true_snp_ids = []
        with open(simdir + "CV_snp_ids.txt", "r") as fin:
            for line in fin:
                true_snp_ids.append(line.strip())
        true_eff_file = glob(simdir + "true_causaleffects*")[0]
        true_eff_lines = []
        true_eff_lines_all = []
        with open(true_eff_file, "r") as fin:
            for line in fin:
                # we skip the latent vars here
                true_eff_lines.append(
                    {k + 1: v for (k, v) in enumerate(line.strip().split()[2:])}
                )
                true_eff_lines_all.append(
                    {k: v for (k, v) in enumerate(line.strip().split())}
                )
        true_dag_mxp = pd.DataFrame(
            true_eff_lines[:num_m], index=true_snp_ids, dtype=float
        )
        true_dag_pxp = np.triu(
            np.array(
                pd.DataFrame(
                    true_eff_lines[-num_p:],
                    index=list(range(1, num_p + 1)),
                    dtype=float,
                )
            ),
            1,
        )
        true_dag_lpxlp = np.triu(
            np.array(pd.DataFrame(true_eff_lines_all[num_m:], dtype=float)), 1
        )
        # mark true 2-2 edges
        m = true_dag_lpxlp != 0
        true_bidirected = np.zeros(shape=(num_p, num_p))
        for i in range(2, num_p + 2):
            for j in range(i + 1, num_p + 2):
                if not m[i, j] and ((m[0, i] and m[0, j]) or (m[1, i] and m[1, j])):
                    true_bidirected[i - 2, j - 2] = 1

        for alpha_e in e_arr:
            # -------------------- CI-GWAS -----------------------------
            est_dir = f"/nfs/scistore17/robingrp/mrobinso/cigwas/ukb/sim/real_marker_sim/cusk/sim{rep}_e{alpha_e}/"
            pag_path = est_dir + "max_sep_min_pc_estimated_pag_cusk2.mtx"
            try:
                full_pag = mmread(pag_path).tocsr()
            except:
                continue
            outdir = wdir + f"cusk/sim{rep}_e{alpha_e}/"

            gr = merge_block_outputs(blockfile, outdir)
            glob_ixs = np.array(sorted(list(gr.gmi.values())))
            est_dag_mxp = pd.DataFrame(
                full_pag[num_p:, :num_p].toarray(),
                columns=list(range(1, num_p + 1)),
                index=rs_ids[glob_ixs],
                dtype=int,
            )
            est_pag_pxp = full_pag[:num_p, :num_p].toarray()
            est_pag_pxp_triu = np.triu(est_pag_pxp, 1)

            p = (true_dag_mxp != 0).sum().sum()
            td_mxp_matched = true_dag_mxp.reindex_like(est_dag_mxp)
            tp = ((est_dag_mxp != 0) & (td_mxp_matched != 0)).sum().sum()
            fp = ((est_dag_mxp != 0) & (td_mxp_matched == 0)).sum().sum()
            mxp_tpr = tp / p
            mxp_fdr = fp / (tp + fp)

            m = (true_bidirected != 0) | (true_dag_pxp != 0)
            p = np.sum(m)
            tp = np.sum(m & (est_pag_pxp_triu != 0))
            fp = np.sum(~m & (est_pag_pxp_triu != 0))
            pxp_tpr = tp / p
            pxp_fdr = fp / (tp + fp)

            orientation_perf = calculate_pxp_orientation_performance_ci_gwas(
                true_dag_pxp != 0, true_bidirected != 0, est_pag_pxp
            )

            # load aces
            pace = np.zeros(shape=(num_p, num_p))
            missing = []
            for i in range(1, num_p + 1):
                for j in range(1, num_p + 1):
                    if i == j:
                        continue
                    file = est_dir + f"ACE_i{i}_j{j}.csv"
                    try:
                        with open(file) as fin:
                            next(fin)
                            sym = fin.readline().strip()
                            if sym == "NA":
                                pace[i - 1, j - 1] = 0.0
                            else:
                                pace[i - 1, j - 1] = eval(sym)
                    except FileNotFoundError as err:
                        # print(err)
                        missing.append((i, j))

            m = true_dag_pxp
            m = np.linalg.inv(np.eye(m.shape[0]) - m.T)
            true_eff = np.triu(m @ m.T, 1)

            mse = np.mean((pace - true_eff) ** 2)

            tp_m = (true_eff != 0) & (pace != 0)
            correct_sign = np.mean(np.sign(true_eff[tp_m]) == np.sign(pace[tp_m]))

            rows.append(
                {
                    "x -> y fdr": mxp_fdr,
                    "x -> y tpr": mxp_tpr,
                    "-> orientation": orientation_perf.directed,
                    "<-> orientation": orientation_perf.bidirected,
                    "mse": mse,
                    "sign": correct_sign,
                    "y -> y fdr": pxp_fdr,
                    "y -> y tpr": pxp_tpr,
                    "rep": rep,
                    "alpha": 10 ** (-alpha_e),
                    "method": "ci-gwas",
                }
            )

            # -------------------- MR -----------------------------
            for mr_str in ["cause", "presso", "mvpresso", "ivw", "mvivw"]:
                mr_p = np.ones((num_p, num_p))
                mr_est = np.zeros((num_p, num_p))
                for outcome in range(0, num_p):
                    try:
                        fpaths = glob(
                            est_dir
                            + f"mr_{mr_str}_alpha*_sim{rep}_outcome_y{outcome + 1}_seed1000"
                        )
                        if len(fpaths) == 0:
                            continue
                        fpath = fpaths[0]
                        mr_res_df = pd.read_csv(fpath)
                        exposures = [int(s[1:]) - 1 for s in mr_res_df["Exposure"]]
                        mr_p[exposures, outcome] = mr_res_df["p"].values
                        mr_est[exposures, outcome] = mr_res_df["est"].values
                    except FileNotFoundError:
                        continue
                mr_links = mr_p <= (10.0**-alpha_e)
                mr_adj = make_adj_symmetric(mr_links)
                # are we doing this?
                mr_est[~mr_links] = 0.0

                m = (true_bidirected != 0) | (true_dag_pxp != 0)
                p = np.sum(m)
                tp = np.sum(m & (mr_adj != 0))
                pxp_tpr = tp / p

                mse = np.mean((mr_est - true_eff) ** 2)
                tp_m = (true_eff != 0) & (mr_est != 0)
                if np.sum(tp_m) == 0:
                    correct_sign = np.nan
                else:
                    correct_sign = np.mean(
                        np.sign(true_eff[tp_m]) == np.sign(mr_est[tp_m])
                    )

                orientation_perf = calculate_pxp_orientation_performance_mr(
                    true_dag_pxp != 0, true_bidirected != 0, mr_links
                )

                orientation_perf_rel_to_mr = (
                    compare_ci_gwas_orientation_performance_to_mr(
                        true_dag_pxp, true_bidirected, mr_links, est_pag_pxp
                    )
                )

                rows.append(
                    {
                        "mr_pos_tpr": orientation_perf_rel_to_mr.mr_pos_tpr,
                        "mr_neg_tdr": orientation_perf_rel_to_mr.mr_neg_tdr,
                        "-> orientation": orientation_perf.directed,
                        "<-> orientation": orientation_perf.bidirected,
                        "mse": mse,
                        "sign": correct_sign,
                        "y -> y tpr": pxp_tpr,
                        "rep": rep,
                        "alpha": 10 ** (-alpha_e),
                        "method": mr_str,
                    }
                )

    return pd.DataFrame(rows)


def load_n16k_m1600_simulation_results(
    mr=True, e_arr=list(range(2, 9))
) -> pd.DataFrame:
    data_dir = (
        "/nfs/scistore17/robingrp/nmachnik/mnt_home/projects/ci-gwas/simulations/data/"
    )
    cusk_dir = (
        "/nfs/scistore17/robingrp/nmachnik/mnt_home/projects/ci-gwas/simulations/cusk/"
    )
    mr_git_thr_pdir = "/nfs/scistore17/robingrp/nmachnik/mnt_home/projects/ci-gwas/simulations/mr_results/"

    d = 1
    l = 14
    n_arr = [16000]
    m_arr = [1600]
    # m_arr = [200, 400, 1600]
    alpha_str_d = {
        1: "0.1",
        2: "0.01",
        3: "0.001",
        4: "1e-04",
        5: "1e-05",
        6: "1e-06",
        7: "1e-07",
        8: "1e-08",
    }
    alpha_str_arr = [alpha_str_d[e] for e in e_arr]
    rep_arr = list(range(1, 21))
    num_phen = 10

    rows = []

    for n, m, rep in itertools.product(n_arr, m_arr, rep_arr):
        # load true dag
        dag_path = data_dir + f"./true_adj_mat_n{n}_SNP_{m}_it_{rep}.mtx"
        dag = mmread(dag_path).tocsr()
        pdag = dag[-(num_phen):, -(num_phen):].toarray()

        # load true causal effects
        eff_path = data_dir + f"./true_trait_causaleffects_n{n}_SNP_{m}_it_{rep}.mtx"
        eff = mmread(eff_path).tocsr()
        peff = eff[-(num_phen):, -(num_phen):].toarray()
        peff = np.triu(peff, k=1)
        peff[path_in_sem(pdag) == 0] = 0

        # load true dag
        dag_mxp = dag.toarray()[:m, -num_phen:]

        for a_str in alpha_str_arr:
            rows.append(
                {
                    "mse": np.sum(peff**2),
                    "n": n,
                    "m": m,
                    "rep": rep,
                    "alpha": float(a_str),
                    "method": "0-const",
                }
            )

        if mr:
            for mr in mr_cn:
                for a_str in alpha_str_arr:
                    try:
                        mr_file = (
                            mr_git_thr_pdir
                            + f"mr_{mr.method}_n{n}_SNP_{m}_alpha_{a_str}_it_{rep}.csv"
                        )
                        mr_results = pd.read_csv(mr_file)
                        try:
                            mr_results["i"] = mr_results[mr.exposure].apply(
                                lambda x: int(x.split("Y")[1]) - 1
                            )
                            mr_results["j"] = mr_results[mr.outcome].apply(
                                lambda x: int(x.split("Y")[1]) - 1
                            )
                        except KeyError as error:
                            print(f"failed in mr file: {mr_file}")
                            print(error)
                            continue

                        # git_path = (
                        #     mr_git_thr_pdir
                        #     + f"git_n{n}_SNP_{m}_alpha_{a_str}_it_{rep}.csv"
                        # )
                        # git_mxp = pd.read_csv(git_path).to_numpy()

                        # mxp_p = np.sum(dag_mxp != 0)
                        # mxp_tp = np.sum((dag_mxp != 0) & (git_mxp))
                        # mxp_fp = np.sum((dag_mxp == 0) & (git_mxp != 0))
                        # mxp_fdr = mxp_fp / (mxp_fp + mxp_tp)
                        # mxp_tpr = mxp_tp / mxp_p

                        pvals = np.ones(shape=(num_phen, num_phen))
                        try:
                            pvals[mr_results["i"], mr_results["j"]] = mr_results[mr.p]
                        except KeyError as err:
                            print(err)
                            print(mr_results)
                            print(mr_file)
                            exit()
                        effects = np.zeros(shape=(num_phen, num_phen))
                        effects[mr_results["i"], mr_results["j"]] = mr_results[
                            mr.estimate
                        ]
                        adj = pvals < float(a_str)
                        effects[~adj] = 0
                        perf = calulate_performance_metrics(pdag, peff, adj, effects)
                        rows.append(
                            {
                                # "marker-trait fdr": mxp_fdr,
                                # "marker-trait tpr": mxp_tpr,
                                "edge orientation": calculate_pxp_orientation_performance_mr(
                                    pdag != 0, adj
                                ),
                                "mse": perf.mse,
                                "var": perf.var,
                                "bias": perf.bias,
                                "mse_tp": perf.mse_tp,
                                "var_tp": perf.var_tp,
                                "bias_tp": perf.bias_tp,
                                # "fdr": np.nan,
                                "y -> y tpr": perf.tpr,
                                "n": n,
                                "m": m,
                                "rep": rep,
                                "alpha": float(a_str),
                                "method": f"{mr.method}",
                            }
                        )
                    except FileNotFoundError as error:
                        print(error)
                        rows.append(
                            {
                                "edge orientation": np.nan,
                                "mse": np.nan,
                                "var": np.nan,
                                "bias": np.nan,
                                "mse_tp": np.nan,
                                "var_tp": np.nan,
                                "bias_tp": np.nan,
                                "fdr": np.nan,
                                "tpr": np.nan,
                                "n": n,
                                "m": m,
                                "rep": rep,
                                "alpha": float(a_str),
                                "method": f"{mr.method}",
                            }
                        )

        for e in e_arr:
            indir = cusk_dir + f"simpc_d{d}_l{l}_e{e}_i{rep}_n{n}_m{m}/"

            mdim_path = indir + "skeleton.mdim"
            with open(mdim_path, "r") as fin:
                num_var, num_phen, max_level = [
                    int(elem) for elem in fin.readline().split()
                ]

            # load skeleton
            adj = np.fromfile(indir + "skeleton.adj", dtype=np.int32).reshape(
                num_var, num_var
            )
            padj = adj[-(num_phen):, -(num_phen):]

            # load pag
            pag_path = indir + "max_sep_min_pc_estimated_pag_cusk2.mtx"
            try:
                est_pag = mmread(pag_path).tocsr()[-(num_phen):, -(num_phen):].toarray()
            except FileNotFoundError as err:
                print(err)
                continue

            est_dag = pag_to_dag_directed(est_pag)

            # load var indices
            # var_ixs = np.fromfile(indir + "skeleton.ixs", dtype=np.int32)

            # load aces
            pace = np.zeros(shape=(num_phen, num_phen))
            missing = []
            for i in range(1, num_phen + 1):
                for j in range(1, num_phen + 1):
                    if i == j:
                        continue
                    file = indir + f"estimated_causaleffect_i{i}_j{j}_cusk2.csv"
                    try:
                        with open(file) as fin:
                            sym = fin.readline().strip()
                            if sym == "NaN":
                                pace[i - 1, j - 1] = 0.0
                            else:
                                pace[i - 1, j - 1] = eval(sym)
                    except FileNotFoundError as err:
                        missing.append((i, j))

            ixs_path = indir + "skeleton.ixs"
            var_ixs = np.fromfile(ixs_path, dtype=np.int32)

            # load true dag
            dag_path = data_dir + f"./true_adj_mat_n{n}_SNP_{m}_it_{rep}.mtx"
            full_dag = mmread(dag_path).tocsr()

            # load pag
            full_pag = mmread(pag_path).tocsr()

            dag_mxp = full_dag.toarray()[:m, -num_phen:]
            pag_mxp_sub = full_pag[:, -num_phen:].toarray()
            pag_mxp = np.zeros_like(dag_mxp)
            for sub_ix, glob_ix in enumerate(var_ixs):
                if glob_ix < m:
                    pag_mxp[glob_ix, :] = pag_mxp_sub[sub_ix, :]

            mxp_p = np.sum(dag_mxp != 0)
            mxp_tp = np.sum((dag_mxp != 0) & (pag_mxp != 0))
            mxp_fp = np.sum((dag_mxp == 0) & (pag_mxp != 0))
            mxp_fdr = mxp_fp / (mxp_fp + mxp_tp)
            mxp_tpr = mxp_tp / mxp_p

            perf = calulate_performance_metrics(pdag, peff, padj, pace)
            # perf = calulate_performance_metrics(pdag, peff, est_dag, pace)
            rows.append(
                {
                    "x -> y fdr": mxp_fdr,
                    "x -> y tpr": mxp_tpr,
                    "edge orientation": calculate_pxp_orientation_performance_ci_gwas(
                        pdag != 0, est_pag
                    ),
                    "mse": perf.mse,
                    "var": perf.var,
                    "bias": perf.bias,
                    "mse_tp": perf.mse_tp,
                    "var_tp": perf.var_tp,
                    "bias_tp": perf.bias_tp,
                    "y -> y fdr": perf.fdr,
                    "y -> y tpr": perf.tpr,
                    "n": n,
                    "m": m,
                    "rep": rep,
                    "alpha": 10 ** (-e),
                    "method": "ci-gwas",
                }
            )

    return pd.DataFrame(rows)


def plot_simulation_results():
    d = 1
    l = 6
    n_arr = [2000, 4000, 8000, 16000]
    m_arr = [200, 400, 800, 1600]
    # m_arr = [200, 400, 1600]
    e_arr = list(range(1, 9))
    rep_arr = list(range(1, 21))
    num_phen = 10

    df = load_n16k_m1600_simulation_results()
    dfs = df.loc[(df["n"] == 16000) & (df["m"] == 1600)]
    gr = dfs.groupby(["method", "alpha"])
    means = gr.mean()
    stds = gr.std()

    alphas = 10.0 ** -np.array(e_arr[::-1])
    # methods = ['ci-gwas', 'cause', 'mrpresso', 'egger', 'ivw']
    methods = ["ci-gwas", "cause", "mrpresso", "ivw"]
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    def plot_bars(x_vals, metric, means, stds, ax, title):
        ax.set_title(title, **title_kw)
        x = np.arange(len(x_vals))  # the label locations
        width = 0.15  # the width of the bars
        multiplier = 0

        handles = []
        for method in methods:
            mu = means.loc[method][metric]
            sig = stds.loc[method][metric]
            offset = width * multiplier
            try:
                bars = ax.bar(
                    x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "e)" or title == "d)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x + width, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=XLABEL_ROTATION,
            ha="right",
            rotation_mode="anchor",
        )
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    fig = plt.figure(layout="tight", figsize=(10, 8))
    ax_dict = fig.subplot_mosaic(
        """
        fa
        bc
        de
        """,
        empty_sentinel="X",
        sharex=True,
    )

    plot_bars(alphas, "fdr", means, stds, ax_dict["a"], "a)")
    plot_bars(alphas, "tpr", means, stds, ax_dict["b"], "b)")
    plot_bars(alphas, "mse", means, stds, ax_dict["c"], "c)")
    plot_bars(alphas, "bias", means, stds, ax_dict["d"], "d)")
    h = plot_bars(alphas, "var", means, stds, ax_dict["e"], "e)")
    ax_dict["f"].legend(
        handles=h,
        labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "IVW Regression"],
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=2,
        # title="method",
    )
    ax_dict["f"].axis("off")


def plot_real_data_simulation_results(
    rep_arr=list(range(1, 41)),
    e_arr=[3, 4, 6, 8],
    fig1_path=None,
    fig2_path=None,
):
    bar_xlabel_rotation = 45

    df = load_real_data_simulation_results(e_arr=e_arr, rep_arr=rep_arr)
    gr = df.groupby(["method", "alpha"])
    means = gr.mean()
    stds = gr.std()

    alphas = 10.0 ** -np.array(e_arr[::-1])
    # methods = ['ci-gwas', 'cause', 'mrpresso', 'egger', 'ivw']
    # methods = ["ci-gwas", "cause", "mrpresso", "ivw"]
    methods = ["ci-gwas", "cause", "presso", "mvpresso", "ivw", "mvivw"]
    # methods = ["ci-gwas"]
    mse_methods = ["ci-gwas", "cause", "presso", "mvpresso", "ivw", "mvivw"]
    rel_mr_methods = ["ci-gwas", "cause", "presso", "mvpresso", "ivw", "mvivw"]
    # mse_methods = ["ci-gwas"]
    # mse_methods = ["ci-gwas", "cause", "mrpresso", "ivw", "0-const"]
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    def plot_ci_gwas_bars(x_vals, metric, means, stds, ax, title, axhline=False):
        ax.set_title(title, **title_kw)
        x_sorted = sorted(list(x_vals))
        x = np.arange(len(x_vals))  # the label locations
        width = 0.5  # the width of the bars
        multiplier = 0

        handles = []
        for method in ["ci-gwas"]:
            mu = means.loc[method][metric]
            sig = stds.loc[method][metric]
            loc_x = np.array(
                [x_sorted.index(e) for e in means.loc[method].index.values]
            )
            offset = width * multiplier
            try:
                bars = ax.bar(
                    loc_x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "e)" or title == "d)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=bar_xlabel_rotation,
            ha="right",
            rotation_mode="anchor",
        )
        if axhline:
            ax.axhline(0.05, linestyle="dotted", color="gray")
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    def plot_bars(x_vals, metric, means, stds, ax, title, methods):
        ax.set_title(title, **title_kw)
        x_sorted = sorted(list(x_vals))
        x = np.arange(len(x_vals))  # the label locations
        width = 0.10  # the width of the bars
        multiplier = 0

        handles = []
        for method in methods:
            mu = means.loc[method][metric]
            sig = stds.loc[method][metric]
            loc_x = np.array(
                [x_sorted.index(e) for e in means.loc[method].index.values]
            )
            offset = width * multiplier
            try:
                bars = ax.bar(
                    loc_x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "g)" or title == "h)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x + width, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=bar_xlabel_rotation,
            ha="right",
            rotation_mode="anchor",
        )
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    # ---------- adjacency performance

    fig = plt.figure(
        # layout="tight",
        figsize=(9, 6)
    )
    ax_dict = fig.subplot_mosaic(
        """
        XiX
        abc
        ddd
        """,
        empty_sentinel="X",
        sharex=False,
        # set the height ratios between the rows
        height_ratios=[0.1, 0.4, 0.5],
    )

    ax_dict["d"].set_xlabel(r"$\alpha$")

    plot_ci_gwas_bars(
        alphas, "x -> y fdr", means, stds, ax_dict["a"], "a)", axhline=True
    )
    ax_dict["a"].set_ylabel(r"$x - y \ FDR$")
    plot_ci_gwas_bars(alphas, "x -> y tpr", means, stds, ax_dict["b"], "b)")
    ax_dict["b"].set_ylabel(r"$x - y \ TPR$")
    plot_ci_gwas_bars(
        alphas, "y -> y fdr", means, stds, ax_dict["c"], "c)", axhline=True
    )
    ax_dict["c"].set_ylabel(r"$y - y \ FDR$")
    h = plot_bars(
        alphas, "y -> y tpr", means, stds, ax_dict["d"], "d)", methods=methods
    )
    ax_dict["d"].set_ylabel(r"$y - y \ TPR$")
    # plot_bars(alphas, "edge orientation", means, stds, ax_dict["e"], "e)", methods=methods)
    # h = plot_bars(alphas, "mse", means, stds, ax_dict["f"], "f)", methods=mse_methods)
    # ax_dict["f"].set_ylabel("MSE")
    # plot_bars(alphas, "bias", means, stds, ax_dict["g"], "g)")
    # h = plot_bars(alphas, "var", means, stds, ax_dict["h"], "h)")
    ax_dict["i"].legend(
        handles=h,
        labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "MR-PRESSO (MV)", "IVW", "IVW (MV)"],
        # labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "IVW Regression"],
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=3,
        # title="method",
    )
    ax_dict["i"].axis("off")
    # _ = [ax_dict[i].tick_params(labelbottom=False) for i in "de"]
    # _ = [ax_dict[i].sharex(ax_dict['f']) for i in 'de']
    fig.subplots_adjust(wspace=0.7, hspace=1.3)
    if fig1_path is not None:
        plt.savefig(fig1_path, bbox_inches="tight")

    # ---------- orientation and ace performance

    fig = plt.figure(
        # layout="tight",
        figsize=(9, 6)
    )
    ax_dict = fig.subplot_mosaic(
        """
        XiX
        aaa
        bbb
        ccc
        """,
        empty_sentinel="X",
        sharex=False,
        # set the height ratios between the rows
        height_ratios=[0.1, 0.3, 0.3, 0.3],
    )

    ax_dict["c"].set_xlabel(r"$\alpha$")

    plot_bars(
        alphas, "mr_pos_tpr", means, stds, ax_dict["a"], "a)", methods=rel_mr_methods
    )
    ax_dict["a"].set_ylabel(r"$\rightarrow TPR_{MR+}$")
    plot_bars(
        alphas,
        "-> orientation",
        means,
        stds,
        ax_dict["b"],
        "b)",
        methods=rel_mr_methods,
    )
    ax_dict["b"].set_ylabel(r"$\rightarrow TDR$")
    h = plot_bars(alphas, "mse", means, stds, ax_dict["c"], "c)", methods=mse_methods)
    ax_dict["c"].set_ylabel("MSE")

    ax_dict["i"].legend(
        handles=h,
        labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "MR-PRESSO (MV)", "IVW", "IVW (MV)"],
        # labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "IVW Regression"],
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=3,
        # title="method",
    )
    ax_dict["i"].axis("off")
    _ = [ax_dict[i].tick_params(labelbottom=False) for i in "ab"]
    _ = [ax_dict[i].sharex(ax_dict["c"]) for i in "ab"]
    fig.subplots_adjust(wspace=0.7, hspace=1.3)
    if fig2_path is not None:
        plt.savefig(fig2_path, bbox_inches="tight")


def plot_simulation_results_figure_2():
    bar_xlabel_rotation = 45
    d = 1
    l = 6
    n_arr = [2000, 4000, 8000, 16000]
    m_arr = [200, 400, 800, 1600]
    # m_arr = [200, 400, 1600]
    # e_arr = list(range(2, 9))
    e_arr = [2, 4, 6, 8]
    rep_arr = list(range(1, 21))
    num_phen = 10

    df = load_n16k_m1600_simulation_results(mr=True, e_arr=e_arr)
    dfs = df.loc[(df["n"] == 16000) & (df["m"] == 1600)]
    gr = dfs.groupby(["method", "alpha"])
    means = gr.mean()
    stds = gr.std()

    alphas = 10.0 ** -np.array(e_arr[::-1])
    # methods = ['ci-gwas', 'cause', 'mrpresso', 'egger', 'ivw']
    methods = ["ci-gwas", "cause", "mrpresso", "ivw"]
    mse_methods = ["ci-gwas", "cause", "mrpresso", "ivw", "0-const"]
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    def plot_ci_gwas_bars(x_vals, metric, means, stds, ax, title, axhline=False):
        ax.set_title(title, **title_kw)
        x = np.arange(len(x_vals))  # the label locations
        width = 0.5  # the width of the bars
        multiplier = 0

        handles = []
        for method in ["ci-gwas"]:
            mu = means.loc[method][metric]
            sig = stds.loc[method][metric]
            offset = width * multiplier
            try:
                bars = ax.bar(
                    x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "e)" or title == "d)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=bar_xlabel_rotation,
            ha="right",
            rotation_mode="anchor",
        )
        if axhline:
            ax.axhline(0.05, linestyle="dotted", color="gray")
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    def plot_bars(x_vals, metric, means, stds, ax, title, methods):
        ax.set_title(title, **title_kw)
        x = np.arange(len(x_vals))  # the label locations
        width = 0.14  # the width of the bars
        multiplier = 0

        handles = []
        for method in methods:
            mu = means.loc[method][metric]
            sig = stds.loc[method][metric]
            offset = width * multiplier
            try:
                bars = ax.bar(
                    x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "g)" or title == "h)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x + width, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=bar_xlabel_rotation,
            ha="right",
            rotation_mode="anchor",
        )
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    fig = plt.figure(
        # layout="tight",
        figsize=(9, 10)
    )
    ax_dict = fig.subplot_mosaic(
        """
        XiX
        abc
        ddd
        eee
        fff
        """,
        empty_sentinel="X",
        sharex=False,
        # set the height ratios between the rows
        height_ratios=[0.1, 0.4, 0.5, 0.5, 0.5],
    )

    ax_dict["f"].set_xlabel(r"$\alpha$")

    plot_ci_gwas_bars(
        alphas, "x -> y fdr", means, stds, ax_dict["a"], "a)", axhline=True
    )
    ax_dict["a"].set_ylabel(r"x$\rightarrow$y FDR")
    plot_ci_gwas_bars(alphas, "x -> y tpr", means, stds, ax_dict["b"], "b)")
    ax_dict["b"].set_ylabel(r"x$\rightarrow$y TPR")
    plot_ci_gwas_bars(
        alphas, "y -> y fdr", means, stds, ax_dict["c"], "c)", axhline=True
    )
    ax_dict["c"].set_ylabel(r"y$\rightarrow$y FDR")
    plot_bars(alphas, "y -> y tpr", means, stds, ax_dict["d"], "d)", methods=methods)
    ax_dict["d"].set_ylabel(r"y$\rightarrow$y TPR")
    plot_bars(
        alphas, "edge orientation", means, stds, ax_dict["e"], "e)", methods=methods
    )
    h = plot_bars(alphas, "mse", means, stds, ax_dict["f"], "f)", methods=mse_methods)
    ax_dict["f"].set_ylabel("MSE")
    # plot_bars(alphas, "bias", means, stds, ax_dict["g"], "g)")
    # h = plot_bars(alphas, "var", means, stds, ax_dict["h"], "h)")
    ax_dict["i"].legend(
        handles=h,
        labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "IVW Regression", "0-const"],
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=5,
        # title="method",
    )
    ax_dict["i"].axis("off")
    _ = [ax_dict[i].tick_params(labelbottom=False) for i in "de"]
    _ = [ax_dict[i].sharex(ax_dict["f"]) for i in "de"]
    fig.subplots_adjust(wspace=0.4, hspace=1.1)


def plot_simulation_results_sup_tp_with_orientation():
    d = 1
    l = 6
    n_arr = [2000, 4000, 8000, 16000]
    m_arr = [200, 400, 800, 1600]
    # m_arr = [200, 400, 1600]
    e_arr = list(range(1, 9))
    rep_arr = list(range(1, 21))
    num_phen = 10

    df = load_n16k_m1600_simulation_results()
    dfs = df.loc[(df["n"] == 16000) & (df["m"] == 1600)]
    gr = dfs.groupby(["method", "alpha"])
    means = gr.mean()
    stds = gr.std()

    alphas = 10.0 ** -np.array(e_arr[::-1])
    # methods = ['ci-gwas', 'cause', 'mrpresso', 'egger', 'ivw']
    methods = ["ci-gwas", "cause", "mrpresso", "ivw"]
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    def plot_bars(x_vals, metric, means, stds, ax, title):
        ax.set_title(title, **title_kw)
        x = np.arange(len(x_vals))  # the label locations
        width = 0.15  # the width of the bars
        multiplier = 0

        handles = []
        for method in methods:
            mu = means.loc[method][metric]
            sig = stds.loc[method][metric]
            offset = width * multiplier
            try:
                bars = ax.bar(
                    x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "e)" or title == "f)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x + width, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=XLABEL_ROTATION,
            ha="right",
            rotation_mode="anchor",
        )
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    fig = plt.figure(layout="tight", figsize=(10, 8))
    ax_dict = fig.subplot_mosaic(
        """
        gg
        ab
        cd
        ef
        """,
        empty_sentinel="X",
        sharex=True,
        # set the height ratios between the rows
        height_ratios=[0.2, 1, 1, 1],
    )

    plot_bars(alphas, "fdr", means, stds, ax_dict["a"], "a)")
    plot_bars(alphas, "tpr", means, stds, ax_dict["b"], "b)")
    plot_bars(alphas, "edge orientation", means, stds, ax_dict["c"], "c)")
    plot_bars(alphas, "mse_tp", means, stds, ax_dict["d"], "d)")
    plot_bars(alphas, "bias_tp", means, stds, ax_dict["e"], "e)")
    h = plot_bars(alphas, "var_tp", means, stds, ax_dict["f"], "f)")
    ax_dict["g"].legend(
        handles=h,
        labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "IVW Regression"],
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=4,
        # title="method",
    )
    ax_dict["g"].axis("off")


def plot_simulation_results_sup_tp():
    d = 1
    l = 6
    n_arr = [2000, 4000, 8000, 16000]
    m_arr = [200, 400, 800, 1600]
    # m_arr = [200, 400, 1600]
    e_arr = list(range(1, 9))
    rep_arr = list(range(1, 21))
    num_phen = 10

    df = load_n16k_m1600_simulation_results(mr_git_thr=True)
    dfs = df.loc[(df["n"] == 16000) & (df["m"] == 1600)]
    gr = dfs.groupby(["method", "alpha"])
    means = gr.mean()
    stds = gr.std()

    alphas = 10.0 ** -np.array(e_arr[::-1])
    # methods = ['ci-gwas', 'cause', 'mrpresso', 'egger', 'ivw']
    methods = ["ci-gwas", "cause", "mrpresso", "ivw"]
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    def plot_bars(x_vals, metric, means, stds, ax, title):
        ax.set_title(title, **title_kw)
        x = np.arange(len(x_vals))  # the label locations
        width = 0.15  # the width of the bars
        multiplier = 0

        handles = []
        for method in methods:
            mu = means.loc[method][metric]
            sig = stds.loc[method][metric]
            offset = width * multiplier
            try:
                bars = ax.bar(
                    x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "b)" or title == "c)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x + width, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=XLABEL_ROTATION,
            ha="right",
            rotation_mode="anchor",
        )
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    fig = plt.figure(layout="tight", figsize=(10, 6))
    ax_dict = fig.subplot_mosaic(
        """
        fa
        bc
        """,
        empty_sentinel="X",
        sharex=True,
    )

    plot_bars(alphas, "mse_tp", means, stds, ax_dict["a"], "a)")
    plot_bars(alphas, "bias_tp", means, stds, ax_dict["b"], "b)")
    h = plot_bars(alphas, "var_tp", means, stds, ax_dict["c"], "c)")
    ax_dict["f"].legend(
        handles=h,
        labels=["CI-GWAS", "CAUSE", "MR-PRESSO", "IVW Regression"],
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=2,
        # title="method",
    )
    ax_dict["f"].axis("off")


def plot_simulation_results_sup(alpha=10 ** (-8)):
    d = 1
    l = 6
    n_arr = [2000, 4000, 8000, 16000]
    m_arr = [200, 400, 800, 1600]
    # m_arr = [200, 400, 1600]
    e_arr = list(range(1, 9))
    rep_arr = list(range(1, 21))
    num_phen = 10

    df = load_simulation_results()
    gr = df.groupby(["method", "alpha", "n", "m"])
    means = gr.mean()
    stds = gr.std()

    # methods = ['ci-gwas', 'cause', 'mrpresso', 'egger', 'ivw']
    methods = ["ci-gwas", "cause", "mrpresso", "ivw"]
    # methods = ['ci-gwas']
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    def plot_ci_gwas_bars(x_vals, metric, means, stds, ax, title, axhline=False):
        ax.set_title(title, **title_kw)
        x = np.arange(len(x_vals))  # the label locations
        width = 0.5  # the width of the bars
        multiplier = 0

        handles = []
        for method in ["ci-gwas"]:
            mu = means.loc[method, alpha][metric]
            sig = stds.loc[method, alpha][metric]
            offset = width * multiplier
            try:
                bars = ax.bar(
                    x + offset, mu, width, yerr=sig, capsize=1.5, label=method
                )
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == "e)" or title == "d)":
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x, x_vals)
        plt.setp(
            ax.get_xticklabels(),
            rotation=XLABEL_ROTATION,
            ha="right",
            rotation_mode="anchor",
        )
        if axhline:
            ax.axhline(0.05, linestyle="dashed", color="k")
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    def plot_bars(x_tup, metric, means, stds, ax, title):
        ax.set_title(title, **title_kw)
        x = np.arange(len(x_tup))  # the label locations
        width = 0.15  # the width of the bars
        multiplier = 0

        handles = []
        for method in methods:
            mu = means.loc[method, alpha][metric]
            sig = stds.loc[method, alpha][metric]
            offset = width * multiplier
            try:
                bars = ax.bar(x + offset, mu, width, yerr=sig, capsize=1, label=method)
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel(metric)
        if title == "g)" or title == "h)":
            ax.set_xlabel("(n, m)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # ax.set_title('Penguin attributes by species')
        ax.set_xticks(x + width, x_tup)
        plt.setp(
            ax.get_xticklabels(),
            rotation=XLABEL_ROTATION,
            ha="right",
            rotation_mode="anchor",
        )
        # ax.legend(loc='upper left', ncols=3)
        return handles

    x_tup = list(itertools.product(n_arr, m_arr))

    fig = plt.figure(layout="tight", figsize=(15, 8))
    ax_dict = fig.subplot_mosaic(
        """
        ii
        ab
        cd
        ef
        gh
        """,
        empty_sentinel="X",
        sharex=True,
        # set the height ratios between the rows
        height_ratios=[0.2, 0.6, 1, 1, 1],
    )

    plot_ci_gwas_bars(
        x_tup, "marker-trait fdr", means, stds, ax_dict["a"], "a)", axhline=True
    )
    plot_ci_gwas_bars(x_tup, "marker-trait tpr", means, stds, ax_dict["b"], "b)")
    plot_ci_gwas_bars(x_tup, "fdr", means, stds, ax_dict["c"], "c)", axhline=True)
    plot_bars(x_tup, "tpr", means, stds, ax_dict["d"], "d)")
    # plot_bars(x_tup, "edge orientation", means, stds, ax_dict["e"], "e)")
    plot_bars(x_tup, "mse", means, stds, ax_dict["f"], "e)")
    plot_bars(x_tup, "bias", means, stds, ax_dict["g"], "f)")
    h = plot_bars(x_tup, "var", means, stds, ax_dict["h"], "g)")
    ax_dict["i"].legend(
        handles=h,
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=4,
        title="method",
    )
    ax_dict["i"].axis("off")
    fig.suptitle(rf"$\alpha={{{alpha}}}$")


def plot_block_size_experiment_results():
    num_phen = 17
    block_sizes = np.arange(1, 13) * 10**3
    basepath = "/nfs/scistore17/robingrp/human_data/causality/block_size_effect/"
    bim_path = basepath + "ukb22828_UKB_EST_v3_ldp08.bim"
    bim = pd.read_csv(bim_path, sep="\t", header=None)

    bps = {}
    block_boundaries = {}

    for bs in block_sizes:
        blockfile = basepath + f"ukb22828_UKB_EST_v3_ldp08_m{bs}_chr1.blocks"
        outdir = basepath + f"bdpc_d1_l6_a1e4_m{bs}_chr1/"
        gr = merge_block_outputs(blockfile, outdir)
        print(bs, gr.gmi)
        bps[bs] = sorted(bim.loc[list(gr.gmi.values())][3].values)
        block_df = pd.read_csv(blockfile, sep="\t", names=["chr", "first", "last"])
        block_boundaries[bs] = np.array(sorted(bim.loc[block_df["first"]][3].values))

    num_markers_selected = [len(bps[bs]) for bs in block_sizes]

    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="tight", figsize=(10, 4))

    ax_dict = fig.subplot_mosaic(
        """
        aaab
        """,
        empty_sentinel="X",
        # set the height ratios between the rows
        # height_ratios=[1, 1, 0.3],
    )

    ax = ax_dict["a"]
    title = "a)"
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for k, v in bps.items():
        x = np.array(v) / 10**6
        y = np.ones_like(x) * k
        ax.plot(x, y, "v", color="k", alpha=0.5)
        ax.plot(
            block_boundaries[k] / 10**6,
            np.ones_like(block_boundaries[k]) * k,
            "|",
            color="r",
            alpha=0.5,
        )
    ax.set_ylabel("max block size")
    ax.set_xlabel("position on chr1 [Mbp]")
    ax.set_title(title, **title_kw)
    ax.set_yticks(block_sizes)

    ax = ax_dict["b"]
    title = "b)"
    ax.plot(block_sizes, num_markers_selected, "o-", color="k")
    ax.set_title(title, **title_kw)
    ax.set_ylabel("# selected markers")
    ax.set_xlabel("max block size")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_alpha_effect_on_davs_results():
    d = 1
    l = 6
    es = [2, 3, 4, 5, 6, 7, 8]
    nrows = len(es) - 1
    pheno_path = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/input.phen"

    ace = {}

    for i, e in enumerate(es):
        bdpc_ace_path_sk = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/bdpc_d{d}_l{l}_a1e{e}/all_merged_ACE_sk_mk3.mtx"

        ace[e] = load_ace(bdpc_ace_path_sk, pheno_path).flatten()

    gmin = np.min([np.min(a) for a in ace.values()])
    gmax = np.max([np.max(a) for a in ace.values()])
    gmax = max(abs(gmin), abs(gmax)) + 0.05

    diag = np.linspace(-gmax, gmax, 10)

    fig, axes = plt.subplots(nrows, nrows, sharex=True, sharey=True, figsize=(8, 8))

    for i in range(nrows):
        for j in range(nrows):
            if i < j:
                axes[i][j].axis("off")
                continue
            ey, ex = es[i + 1], es[j]
            ax = axes[i][j]
            ax.plot(diag, diag, ":", color="gray")
            x, y = ace[ex], ace[ey]
            ax.plot(x, y, "k.")
            ax.set_xlim(-gmax, gmax)
            ax.set_ylim(-gmax, gmax)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            rho = np.corrcoef(x, y)
            ax.set_title(rf"$\rho={np.round(rho[0, 1], 2)}$")
            if j == 0:
                axes[i][j].set_ylabel(rf"$\alpha=10^{{{ey}}}$")
            if i == nrows - 1:
                axes[i][j].set_xlabel(rf"$\alpha=10^{{{ex}}}$")

    # fig.suptitle(r"$ACE; k\leq3; SK$")
    plt.tight_layout()


def plot_compare_max_k_effect_on_ace():
    d = 1
    l = 6
    es = [5, 6, 7, 8]
    nrows = len(es)
    pheno_path = (
        f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/input.phen"
    )

    ace_mk = {}
    ace = {}

    for i, e in enumerate(es):
        bdpc_ace_path_sk_mk = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/bdpc_d{d}_l{l}_a1e{e}/all_merged_ACE_sk_mk3.mtx"
        bdpc_ace_path_sk = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/bdpc_d{d}_l{l}_a1e{e}/all_merged_ACE_sk.mtx"
        ace_mk[e] = load_ace(bdpc_ace_path_sk_mk, pheno_path).flatten()
        ace[e] = load_ace(bdpc_ace_path_sk, pheno_path).flatten()

    gmin = np.min([np.min(a) for a in ace.values()])
    gmax = np.max([np.max(a) for a in ace.values()])
    gmax = max(abs(gmin), abs(gmax)) + 0.05

    diag = np.linspace(-gmax, gmax, 10)

    fig, axes = plt.subplots(1, nrows, sharex=True, sharey=True, figsize=(10, 3))

    for i, e in enumerate(es):
        ax = axes[i]
        ax.plot(diag, diag, ":", color="gray")
        x, y = ace[e], ace_mk[e]
        ax.plot(x, y, "k.")
        ax.set_xlim(-gmax, gmax)
        ax.set_ylim(-gmax, gmax)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        rho = np.corrcoef(x, y)
        ax.set_title(rf"$\rho={np.round(rho[0, 1], 2)}; \alpha=10^{{-{e}}}$")
        ax.set_xlabel("$ACE$")
    axes[0].set_ylabel("$ACE_{k \leq 3}$")
    plt.tight_layout()


def plot_full_ukb_results_figure_3(
    pag_path: str,
    ace_path: str,
    pheno_path: str,
    ace_norm=None,
    p_thr=0.05,
    max_path_len=np.inf,
):
    e = 4
    d = 10000
    # d = 1
    l = 6
    outdir = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/bdpc_d{d}_l{l}_a1e{e}/"
    blockfile = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.blocks"
    pheno_path = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/input.phen"

    z2 = get_skeleton_pleiotropy_mat(outdir, blockfile, pheno_path, max_depth=2)
    # zinf = get_skeleton_pleiotropy_mat(outdir, blockfile, pheno_path, mat_type="union")

    fig = plt.figure(
        # layout="tight",
        figsize=(17, 8.5)
    )
    ax_dict = fig.subplot_mosaic(
        """
        abe
        cde
        """,
        empty_sentinel="X",
        width_ratios=[1, 1, 0.45],
    )

    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    plot_ace(
        ace_path,
        pheno_path,
        ax=ax_dict["d"],
        title="d)",
        cbar_kw=cbar_kw,
        title_kw=title_kw,
        cmap="RdBu",
        norm=ace_norm,
        cbar=False,
    )

    plot_pag(
        pag_path,
        pheno_path,
        ax=ax_dict["c"],
        title="c)",
        title_kw=title_kw,
        edge_encoding=two_common_edge_types,
        cbar_kw=cbar_kw,
    )

    # plot_cause_ci_gwas_ace_comparison(
    #     ace_path, pheno_path, p_thr=p_thr, ax=ax_dict["f"], title="f)", title_kw=title_kw
    # )
    # ax_dict['f'].set_ylim(-0.05, 0.25)

    plot_non_pleio_barplot(
        pag_path, pheno_path, ax=ax_dict["b"], title="b)", title_kw=title_kw
    )

    # norm = mpl.colors.LogNorm(vmin=1, vmax=np.max(z2))
    cmap = "PuBuGn"

    ax_dict["a"].text(0, -2.5, "a)", size=20, verticalalignment="top")
    plot_skeleton_pleiotropy_mat_z(
        z2,
        pheno_path,
        ax=ax_dict["a"],
        title="depth=2",
        # norm=norm,
        cmap=cmap,
        cbar=True,
        # aspect="auto",
    )
    # ax_dict["b"].text(0, -2.5, "b)", size=20, verticalalignment="top")
    # plot_skeleton_pleiotropy_mat_z(
    #     zinf,
    #     pheno_path,
    #     ax=ax_dict["b"],
    #     title="depth=inf",
    #     norm=norm,
    #     cmap=cmap,
    #     cbar=True,
    #     # aspect="auto",
    # )

    plot_ci_gwas_cause_ace_comparison_tri(
        pag_path,
        ace_path,
        pheno_path,
        p_thr,
        ax=ax_dict["e"],
        title="e)",
        title_kw=title_kw,
        cbar_kw=cbar_kw,
        max_path_len=max_path_len,
    )

    _ = [ax_dict[k].set_box_aspect(1) for k in "acd"]
    ax_dict["b"].set_box_aspect(1)
    fig.subplots_adjust(wspace=0.1, hspace=0.4)


def plot_ace_results_comp_cause_production(
    pag_path: str,
    ace_path: str,
    pheno_path: str,
    ace_norm=None,
    p_thr=0.05,
    max_path_len=np.inf,
):
    """
    AT = Asthma
    BMI = Body mass index
    BP = self-report high blood pressure
    CAD = cardiovascular disease
    CHOL = cholesterol
    DBP = diastolic blood pressure
    GLU = blood glucose
    HDL = high density lipoprotein
    HbA1c = hemoglobin A1C average blood sugar
    LDL = low-density lipoprotein
    SBP = systolic blood pressure
    ST = stroke
    T2D = type-2 diabetes
    TRIG = blood triglyceride levels
    WHR = waist-hip ratio
    """

    fig = plt.figure(layout="tight", figsize=(12, 14))
    ax_dict = fig.subplot_mosaic(
        """
        ab
        cd
        ee
        """,
        empty_sentinel="X",
        height_ratios=[1, 1, 0.7],
    )

    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    plot_ace(
        ace_path,
        pheno_path,
        ax=ax_dict["d"],
        title="d)",
        cbar_kw=cbar_kw,
        title_kw=title_kw,
        cmap="Spectral",
        norm=ace_norm,
    )

    plot_pleiotropy_mat(
        pag_path,
        pheno_path,
        ax=ax_dict["a"],
        title="a)",
        cbar_kw=cbar_kw,
        title_kw=title_kw,
        cmap="BuPu",
    )

    plot_pag(
        pag_path,
        pheno_path,
        ax=ax_dict["c"],
        title="c)",
        title_kw=title_kw,
        edge_encoding=two_common_edge_types,
        cbar_kw=cbar_kw,
    )

    plot_non_pleio_barplot(
        pag_path, pheno_path, ax=ax_dict["b"], title="b)", title_kw=title_kw
    )

    plot_direct_link_cause_comparison_wide(
        pag_path,
        pheno_path,
        p_thr=p_thr,
        title="e)",
        title_kw=title_kw,
        ax=ax_dict["e"],
        max_path_len=max_path_len,
    )


def cause_ci_gwas_ace_comparison_df(
    ace_path, pheno_path, p_thr=0.0001, title=None, title_kw={}, ax=None
):
    reg_cfg = {
        "skip_na": True,
        "rf2rf": False,
        "d2d": False,
        "rf2d": True,
        "d2rf": False,
    }

    pnames = get_pheno_codes(pheno_path)
    cause_ys = set(cause_gamma["y1"].values)
    cause_ys.update(set(cause_gamma["y2"].values))
    reg_pnames = cause_ys.intersection(set(pnames))

    ace = load_ace(ace_path, pheno_path)

    rows = []
    significant_rows = []
    for i, pi in enumerate(pnames):
        if pi not in reg_pnames:
            continue
        for j, pj in enumerate(pnames):
            if pj not in reg_pnames:
                continue
            if pi in risk_factors and pj in risk_factors and not reg_cfg["rf2rf"]:
                continue
            if pi in risk_factors and pj in diseases and not reg_cfg["rf2d"]:
                continue
            if pi in diseases and pj in diseases and not reg_cfg["d2d"]:
                continue
            if pi in diseases and pj in risk_factors and not reg_cfg["d2rf"]:
                continue
            cg = cause_gamma[
                (cause_gamma["y1"] == pi) & (cause_gamma["y2"] == pj)
            ].gamma.values
            cp = cause_gamma[
                (cause_gamma["y1"] == pi) & (cause_gamma["y2"] == pj)
            ].p_value.values
            if len(cg) > 1:
                raise ValueError("too many gammas")
            gamma = 0 if len(cg) == 0 else cg[0]
            p_val = 1 if len(cg) == 0 else cp[0]
            if reg_cfg["skip_na"] and (p_val >= p_thr and ace[i, j] == 0):
                continue
            rows.append({"y1": pi, "y2": pj, "gamma": gamma, "ace": ace[i, j]})
            if p_val <= p_thr:
                significant_rows.append(
                    {"y1": pi, "y2": pj, "gamma": gamma, "ace": ace[i, j]}
                )

    return pd.DataFrame(rows)


def plot_cause_ci_gwas_ace_comparison(
    ace_path, pheno_path, p_thr=0.0001, title=None, title_kw={}, ax=None
):
    reg_cfg = {
        "skip_na": True,
        "rf2rf": False,
        "d2d": False,
        "rf2d": True,
        "d2rf": False,
    }

    if ax is None:
        plt.figure(figsize=(10, 10))
        ax = plt.gca()

    pnames = get_pheno_codes(pheno_path)
    cause_ys = set(cause_gamma["y1"].values)
    cause_ys.update(set(cause_gamma["y2"].values))
    reg_pnames = cause_ys.intersection(set(pnames))

    ace = load_ace(ace_path, pheno_path)

    rows = []
    significant_rows = []
    for i, pi in enumerate(pnames):
        if pi not in reg_pnames:
            continue
        for j, pj in enumerate(pnames):
            if pj not in reg_pnames:
                continue
            if pi in risk_factors and pj in risk_factors and not reg_cfg["rf2rf"]:
                continue
            if pi in risk_factors and pj in diseases and not reg_cfg["rf2d"]:
                continue
            if pi in diseases and pj in diseases and not reg_cfg["d2d"]:
                continue
            if pi in diseases and pj in risk_factors and not reg_cfg["d2rf"]:
                continue
            cg = cause_gamma[
                (cause_gamma["y1"] == pi) & (cause_gamma["y2"] == pj)
            ].gamma.values
            cp = cause_gamma[
                (cause_gamma["y1"] == pi) & (cause_gamma["y2"] == pj)
            ].p_value.values
            if len(cg) > 1:
                raise ValueError("too many gammas")
            gamma = 0 if len(cg) == 0 else cg[0]
            p_val = 1 if len(cg) == 0 else cp[0]
            if reg_cfg["skip_na"] and (p_val >= p_thr and ace[i, j] == 0):
                continue
            rows.append({"y1": pi, "y2": pj, "gamma": gamma, "ace": ace[i, j]})
            if p_val <= p_thr:
                significant_rows.append(
                    {"y1": pi, "y2": pj, "gamma": gamma, "ace": ace[i, j]}
                )

    data = pd.DataFrame(rows)
    mr = data["gamma"].values
    ace_flat = data["ace"].values
    exposures = data["y1"].values
    outcomes = data["y2"].values

    sig_data = pd.DataFrame(significant_rows)
    sig_mr = sig_data["gamma"].values
    sig_ace_flat = sig_data["ace"].values

    ax.scatter(
        mr, ace_flat, color="#d8dcd6", edgecolors="#5729ce", s=80, alpha=0.5, zorder=10
    )
    ax.scatter(
        sig_mr,
        sig_ace_flat,
        color="#f10c45",
        edgecolors="#5729ce",
        s=80,
        alpha=0.5,
        zorder=11,
    )

    texts = []
    for i in range(len(mr)):
        txt = rf"{exposures[i]}$\rightarrow${outcomes[i]}"
        texts.append(ax.annotate(txt, (mr[i], ace_flat[i])))

    adjust_text(texts, arrowprops=dict(arrowstyle="->", color="black"))

    ax.grid(linestyle=":")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel("CI-GWAS ACE")
    ax.set_xlabel(r"CAUSE $\gamma$ (Morrison et al. 2020)")
    ax.set_title(title, **title_kw)
    ax.set_box_aspect(1)


def plot_non_pleio_barplot(
    pag_path: str, pheno_path: str, ax=None, title=None, title_kw=None
):
    if ax is None:
        ax = plt.gca()
    ps = pag_exclusive_pleiotropy_sets(pag_path, pheno_path, is_possible_child)
    p_names = get_pheno_codes(pheno_path)
    bd = [len(ps[i, i]) for i in range(len(p_names))]
    ax.bar(range(len(bd)), bd, color="gray")
    ax.set_ylabel("# non-pleiotropic parent markers")
    ax.set_xticks(np.arange(len(bd)), labels=p_names)
    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(linestyle=":", axis="y")
    ax.set_title(title, **title_kw)
    ax.set_box_aspect(1)


def plot_all_alpha_production_pags_and_ace(
    basepath: str, pheno_path: str, l=6, d=1, edge_encoding=two_common_edge_types
):
    fig = plt.figure(layout="tight", figsize=(12, 14))
    ax_dict = fig.subplot_mosaic(
        """
        ac
        bd
        op
        """,
        empty_sentinel="X",
        height_ratios=[1, 1, 0.02],
        # sharex=True,
        # sharey=True,
    )

    ace_norm = mpl.colors.SymLogNorm(vmin=-2.0, vmax=2.0, linthresh=0.01)
    # cbar_kw = {"fraction": 0.046, "pad": 0.04, "shrink": 0.5}
    cbar_kw = {"shrink": 0.4}
    # title_kw = {"loc": "left", "pad": 15, "size": 20}
    title_kw = {}

    pag_panels = ["a", "b"]
    ace_panels = ["c", "d"]
    alphas = [4, 8]

    for pos, (ace_panel, a) in enumerate(zip(ace_panels, alphas)):
        try:
            ace_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_ACE_mpu.mtx"
            # ace_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_ACE_mpu_q_pheno.mtx"
            pag_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_pag_mpu.mtx"

            ima = plot_ace(
                ace_path,
                pheno_path,
                ax=ax_dict[ace_panel],
                title=rf"$\alpha=10^{{-{a}}}$",
                cbar_kw=cbar_kw,
                title_kw=title_kw,
                cmap="PuOr",
                norm=ace_norm,
                cbar=False,
            )

        except FileNotFoundError as e:
            print(e)

    ax = ax_dict["p"]
    cbar = fig.colorbar(ima, cax=ax_dict["p"], **cbar_kw, location="bottom")
    cbar.ax.set_xlabel(r"$ACE \: (y_1 \rightarrow y_2)$")

    for pos, (pag_panel, a) in enumerate(zip(pag_panels, alphas)):
        try:
            ace_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_ACE_mpu.mtx"
            # ace_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_ACE_mpu_q_pheno.mtx"
            pag_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_pag_mpu.mtx"

            imp = plot_pag(
                pag_path,
                pheno_path,
                ax=ax_dict[pag_panel],
                title=rf"$\alpha=10^{{-{a}}}$",
                title_kw=title_kw,
                edge_encoding=edge_encoding,
                cbar_kw=cbar_kw,
                cbar=False,
            )

        except FileNotFoundError as e:
            print(e)

    ne = len(edge_encoding.int_rep)

    norm = mpl.colors.BoundaryNorm(np.linspace(0, ne, ne + 1), ne)
    fmt = mpl.ticker.FuncFormatter(lambda x, pos: edge_encoding.str_rep[norm(x)])

    cbar_kw["ticks"] = np.arange(ne) + 0.5
    cbar_kw["format"] = fmt
    ax = ax_dict["o"]
    cbar = ax.figure.colorbar(imp, cax=ax_dict["o"], location="bottom", **cbar_kw)
    cbar.ax.tick_params(rotation=45)


def plot_all_alpha_production_pags(basepath: str, pheno_path: str, l=6, d=1):
    fig = plt.figure(layout="tight", figsize=(17, 10))
    ax_dict = fig.subplot_mosaic(
        """
        abc
        def
        """,
        empty_sentinel="X",
        sharex=True,
        sharey=True,
    )

    ace_norm = mpl.colors.SymLogNorm(vmin=-2.0, vmax=2.0, linthresh=0.01)
    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    # title_kw = {"loc": "left", "pad": 15, "size": 20}
    title_kw = {}

    # alphas = [2, 3, 4, 5, 6, 7, 8]
    # panels = ["a", "b", "c", "d", "e", "f", "g", "h"]

    alphas = [2, 3, 4, 5, 7, 8]
    panels = ["a", "b", "c", "d", "e", "f"]

    for panel, a in zip(panels, alphas):
        try:
            ace_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_ACE_sk_mk3.mtx"
            pag_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_pag_sk_mk3.mtx"

            im = plot_pag(
                pag_path,
                pheno_path,
                ax=ax_dict[panel],
                title=rf"$\alpha=10^{{-{a}}}$",
                title_kw=title_kw,
                edge_encoding=five_common_edge_types,
                cbar_kw=cbar_kw,
            )

        except FileNotFoundError as e:
            print(e)

    # Create colorbar
    # cbarlabel = r"$ACE \: (y_1 \rightarrow y_2)$"
    # ax = ax_dict["f"]
    # cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")


def plot_all_alpha_production_aces(basepath: str, pheno_path: str, l=6, d=1):
    fig = plt.figure(layout="tight", figsize=(11, 10))
    ax_dict = fig.subplot_mosaic(
        """
        abc
        deX
        """,
        empty_sentinel="X",
        sharex=True,
        sharey=True,
    )

    ace_norm = mpl.colors.SymLogNorm(vmin=-2.0, vmax=2.0, linthresh=0.01)
    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    # title_kw = {"loc": "left", "pad": 15, "size": 20}
    title_kw = {}

    # alphas = [2, 3, 4, 5, 6, 7, 8]
    # panels = ["a", "b", "c", "d", "e", "f", "g", "h"]

    alphas = [3, 4, 5, 7, 8]
    panels = [
        "a",
        "b",
        "c",
        "d",
        "e",
    ]

    for panel, a in zip(panels, alphas):
        try:
            ace_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_ACE_sk_mk3.mtx"
            pag_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_pag_sk_mk3.mtx"

            im = plot_ace(
                ace_path,
                pheno_path,
                ax=ax_dict[panel],
                title=rf"$\alpha=10^{{-{a}}}$",
                cbar_kw=cbar_kw,
                title_kw=title_kw,
                cmap="PuOr",
                norm=ace_norm,
                cbar=False,
            )
        except FileNotFoundError as e:
            print(e)

    # Create colorbar
    # cbarlabel = r"$ACE \: (y_1 \rightarrow y_2)$"
    # ax = ax_dict["f"]
    # cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")


def plot_pag_and_ace(pag_path: str, ace_path: str, pheno_path: str):
    fig = plt.figure(layout="tight", figsize=(11, 10))
    ax_dict = fig.subplot_mosaic(
        """
        ab
        """,
        empty_sentinel="X",
        sharey=False,
    )

    ace_norm = mpl.colors.SymLogNorm(vmin=-2.0, vmax=2.0, linthresh=0.01)
    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    plot_ace(
        ace_path,
        pheno_path,
        ax=ax_dict["b"],
        title="b)",
        cbar_kw=cbar_kw,
        title_kw=title_kw,
        cmap="PuOr",
        norm=ace_norm,
        # ylabel=None,
    )

    plot_pag(
        pag_path,
        pheno_path,
        ax=ax_dict["a"],
        title="a)",
        title_kw=title_kw,
        edge_encoding=two_common_edge_types,
        cbar_kw=cbar_kw,
    )


def merge_parallel_davs_csvs(
    indir: str, infile_suffix: str, outfile: str, num_phen: int
):
    res = np.zeros(shape=(num_phen, num_phen))
    missing = []
    for i in range(1, num_phen + 1):
        for j in range(1, num_phen + 1):
            if i == j:
                continue
            file = indir + infile_suffix + f"_i{i}_j{j}.csv"
            # file = indir + f"all_merged_ACE_sk_mk3_i{i}_j{j}.csv"
            try:
                with open(file) as fin:
                    next(fin)
                    sym = fin.readline().strip()
                    if sym == "NA":
                        res[i - 1, j - 1] = np.nan
                    else:
                        res[i - 1, j - 1] = eval(sym)
            except FileNotFoundError as e:
                missing.append((i, j))

    res[np.isnan(res)] = 0

    scipy.io.mmwrite(outfile, scipy.sparse.coo_matrix(res))

    return missing


def plot_inf_depth_pleio(e=4):
    d = 10000
    l = 6
    outdir = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/bdpc_d{d}_l{l}_a1e{e}/"
    blockfile = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.blocks"
    pheno_path = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/input.phen"

    z1 = get_skeleton_pleiotropy_mat(outdir, blockfile, pheno_path, max_depth=1)
    z2 = get_skeleton_pleiotropy_mat(outdir, blockfile, pheno_path, max_depth=2)
    zinf = get_skeleton_pleiotropy_mat(outdir, blockfile, pheno_path, mat_type="union")

    fig = plt.figure(layout="tight", figsize=(15, 10))
    ax_dict = fig.subplot_mosaic(
        """
        abc
        dXX
        """,
        empty_sentinel="X",
    )

    norm = mpl.colors.LogNorm(vmin=1, vmax=np.max(zinf))
    cmap = "PuBuGn"

    plot_skeleton_pleiotropy_mat_z(
        z1,
        pheno_path,
        ax=ax_dict["a"],
        title="depth=1",
        norm=norm,
        cmap=cmap,
        cbar=False,
        aspect="auto",
    )
    ax_dict["a"].text(0, -2.5, "a)", size=20, verticalalignment="top")
    plot_skeleton_pleiotropy_mat_z(
        z2,
        pheno_path,
        ax=ax_dict["b"],
        title="depth=2",
        norm=norm,
        cmap=cmap,
        cbar=False,
        aspect="auto",
    )
    ax_dict["b"].text(0, -2.5, "b)", size=20, verticalalignment="top")
    plot_skeleton_pleiotropy_mat_z(
        zinf,
        pheno_path,
        ax=ax_dict["c"],
        title="depth=inf",
        norm=norm,
        cmap=cmap,
        cbar=True,
        aspect="auto",
    )
    ax_dict["c"].text(0, -2.5, "c)", size=20, verticalalignment="top")

    x = np.diag(z1)
    y = np.diag(zinf)

    ax = ax_dict["d"]

    # X = x.reshape(-1, 1)
    # linear_regressor = LinearRegression()
    # linear_regressor.fit(X, y)
    # linear_regressor.fit(X, y)
    # x_pred = np.linspace(np.min(x), np.max(x), 100).reshape(-1, 1)
    # y_pred = linear_regressor.predict(x_pred)
    # ax.plot(x_pred, y_pred, color="#696969", lw=2, linestyle="dashed")
    # beta = np.round(linear_regressor.coef_[0], 2)
    # mu = np.round(linear_regressor.intercept_, 2)
    # ax.set_title("c)", **title_kw)
    # ax.text(100, 0.2 * 10 ** 6, rf"$y={mu} + {beta}x$", fontsize=13)
    # ax.text(100, 0.1 * 10 ** 6, rf"$\rho={np.round(np.corrcoef(x, y)[0, 1], 4)}$", fontsize=13)
    ax.plot(x, y, "o")
    ax.set_ylabel("# shared ancestral markers")
    ax.set_xlabel("# parent markers")
    title_kw = {"loc": "left", "pad": 15, "size": 20}
    ax_dict["d"].set_title("d)", **title_kw)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_est_ukb_corr():
    basepath = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/bdpc_d1_l6_a1e2/pheno_sk"
    num_m = 0
    num_p = 9
    marker_offset = 0

    fig = plt.figure(layout="tight", figsize=(11, 5))
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    ax_dict = fig.subplot_mosaic(
        """
        ab
        """,
        empty_sentinel="X",
    )

    estonia_pnames = ["AT", "BMI", "CAD", "DBP", "HT", "SBP", "SMK", "ST", "T2D"]
    estonia_corr_path = (
        "/nfs/scistore17/robingrp/human_data/cigwas_estonia/trait-trait.txt"
    )

    p_corrs_ukb = load_corr_sparse(basepath, num_m, num_p, marker_offset)
    ukb_pnames = ["AT", "BMI", "CAD", "DBP", "HT", "SBP", "SMK", "ST", "T2D"]

    ukb_mat = np.zeros(shape=(9, 9))

    for (i, j), c in p_corrs_ukb.items():
        try:
            ei = estonia_pnames.index(ukb_pnames[i - 1])
            ej = estonia_pnames.index(ukb_pnames[j - 1])
            ukb_mat[ei, ej] = c
        except ValueError:
            continue

    est_mat = []
    with open(estonia_corr_path, "r") as fin:
        next(fin)
        for line in fin:
            est_mat.append([eval(e) for e in line.split()[1:]])
    est_mat = np.array(est_mat)
    for i in range(9):
        for j in range(i + 1, 9):
            est_mat[j, i] = est_mat[i, j]

    x = np.triu(ukb_mat, k=1).flatten()
    y = np.triu(est_mat, k=1).flatten()

    ax = ax_dict["a"]
    X = x.reshape(-1, 1)

    linear_regressor = LinearRegression()
    linear_regressor.fit(X, y)
    linear_regressor.fit(X, y)
    x_pred = np.linspace(np.min(x), np.max(x), 100).reshape(-1, 1)
    y_pred = linear_regressor.predict(x_pred)

    ax.plot(x_pred, y_pred, color="#696969", lw=2, linestyle="dashed")
    beta = np.round(linear_regressor.coef_[0], 2)
    mu = np.round(linear_regressor.intercept_, 2)
    ax.set_title("a)", **title_kw)
    ax.text(0, 0.55, rf"$y={mu} + {beta}x$", fontsize=13)
    ax.text(0, 0.5, rf"$\rho={np.round(np.corrcoef(x, y)[0, 1], 4)}$", fontsize=13)
    ax.plot(x, y, "o")
    ax.set_ylabel("trait-trait corr EST")
    ax.set_xlabel("trait-trait corr UKB")
    title_kw = {"loc": "left", "pad": 15, "size": 20}
    # ax_dict['d'].set_title("d)", **title_kw)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    z = est_mat - ukb_mat
    mask = ~np.tri(z.shape[0], k=-1, dtype=bool)
    z = np.ma.array(z, mask=mask)  # mask out the lower triangle
    cmap = plt.get_cmap("bwr")
    cmap.set_bad("w")  # default value is 'k'
    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    ax_dict["b"]

    heatmap(
        z,
        estonia_pnames,
        estonia_pnames,
        ax_dict["b"],
        cbarlabel=r"$\Delta corr$",
        cmap=cmap,
        cbar_kw=cbar_kw,
        vmin=-0.12,
        vmax=0.12,
        title_kw=title_kw,
        title="b)",
    )


def plot_ukb_age_sex_causal_path_aces():
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="tight", figsize=(7, 10))
    ax_dict = fig.subplot_mosaic(
        """
        abc
        def
        XxX
        """,
        sharex=True,
        height_ratios=[1, 1, 0.05],
        empty_sentinel="X",
    )

    d = 1
    l = 6
    e = 4

    subsets = [
        "f_first_q",
        "f_second_q",
        "f_third_q",
        "m_first_q",
        "m_second_q",
        "m_third_q",
    ]

    subset_print = {
        "m_third_q": rf"$m, T_{3}$",
        "m_second_q": rf"$m, T_{2}$",
        "m_first_q": rf"$m, T_{1}$",
        "f_third_q": rf"$f, T_{3}$",
        "f_second_q": rf"$f, T_{2}$",
        "f_first_q": rf"$f, T_{1}$",
    }

    pheno_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/age_sex/f_first_q.phen"

    for subset, ax_name in zip(subsets, ["a", "b", "c", "d", "e", "f"]):
        pset = "age_sex"
        outdir = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/{pset}/bdpc_{subset}_d{d}_l{l}_a1e{e}/"
        pag_path = outdir + "all_merged_pag_mpu.mtx"
        ace_path = outdir + "all_merged_ACE_mpu.mtx"
        # im = plot_ace_rf_to_d(ace_path, pag_path, pheno_path, title=subset_print[subset], ax=ax_dict[ax_name], cbar=False, norm=mpl.colors.SymLogNorm(vmin=-2.0, vmax=2.0, linthresh=0.01), )
        im = plot_ace_rf_to_d(
            ace_path,
            pag_path,
            pheno_path,
            title=subset_print[subset],
            ax=ax_dict[ax_name],
            cbar=False,
            cmap="Greens",
            vmin=0,
            vmax=0.2,
        )

    # cbar_kw = {"fraction": 0.046, "pad": 0.04, "shrink": 0.5}

    cbar_kw = {"fraction": 2, "pad": 0}

    ax = ax_dict["x"]
    cbar = ax.figure.colorbar(im, ax=ax, orientation="horizontal", **cbar_kw)
    cbar.ax.set_xlabel(r"$ACE \: (y_1 \rightarrow y_2)$")
    ax.axis("off")


def plot_ukb_age_sex_marker_positions():
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="tight", figsize=(17, 6))
    ax_dict = fig.subplot_mosaic(
        """
        yyyyyy
        """,
        sharey=False,
        width_ratios=[1, 1, 1, 1, 1, 1],
        empty_sentinel="X",
    )

    ax = ax_dict["y"]

    d = 1
    l = 6
    e = 4

    row_height = 10

    subsets = [
        "m_third_q",
        "m_second_q",
        "m_first_q",
        "f_third_q",
        "f_second_q",
        "f_first_q",
    ]

    subset_print = {
        "m_third_q": rf"$m, T_{3}$",
        "m_second_q": rf"$m, T_{2}$",
        "m_first_q": rf"$m, T_{1}$",
        "f_third_q": rf"$f, T_{3}$",
        "f_second_q": rf"$f, T_{2}$",
        "f_first_q": rf"$f, T_{1}$",
    }

    pheno_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/age_sex/f_first_q.phen"
    bim_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.bim"

    traits_with_parents = set()
    for row, subset in enumerate(subsets):
        pset = "age_sex"
        outdir = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/{pset}/bdpc_{subset}_d{d}_l{l}_a1e{e}/"
        traits_with_parents.update(
            set(
                pd.read_csv(
                    outdir + "marker_phen_assoc.csv", sep="\t"
                ).phenotype.unique()
            )
        )
    traits_with_parents = list(traits_with_parents)

    trait_slot_height = row_height / (len(traits_with_parents) + 1)

    colors = []
    cmap = plt.get_cmap("tab20", len(traits_with_parents))
    for i in range(cmap.N):
        rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for row, subset in enumerate(subsets):
        pset = "age_sex"
        outdir = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/{pset}/bdpc_{subset}_d{d}_l{l}_a1e{e}/"
        df = pd.read_csv(outdir + "marker_phen_assoc.csv", sep="\t")
        for trait_ix, trait in enumerate(traits_with_parents):
            sub_df = df[df["phenotype"] == trait]
            x = (
                sub_df.bp.values
                + np.array([global_chr_starts[c] for c in sub_df.chr.values])
            ) / 10**6
            y = np.ones_like(x) * (row * row_height) + trait_slot_height * (
                trait_ix + 1
            )
            ax.plot(x, y, "o", color=colors[trait_ix], alpha=0.8, label=trait)

    handles = []
    for color, trait in zip(colors, traits_with_parents):
        handles.append(mpatches.Patch(color=color, label=trait))

    ax.set_yticks(
        np.array(range(len(subsets))) * row_height + row_height * 0.5,
        [subset_print[s] for s in subsets],
    )
    ax.set_xticks(
        np.array([global_chr_starts[c] + 0.5 * chr_lengths[c] for c in range(1, 23)])
        / 10**6,
        [f"Chr {c}" for c in range(1, 23)],
    )

    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    ax.legend(
        handles=handles,
        loc="upper center",
        fancybox=True,
        shadow=False,
        ncol=len(traits_with_parents) / 2,
        bbox_to_anchor=(0.5, 1.2),
    )

    for x in np.array(list(global_chr_starts.values())) / 10**6:
        ax.axvline(x, color="gray", linestyle=":")

    for row in range(len(subsets)):
        ax.axhline(row * row_height, color="gray", linestyle=":")


def plot_age_sex_composite_figure():
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="tight", figsize=(17, 12))
    ax_dict = fig.subplot_mosaic(
        """
        yyyyyyy
        XXXXXXX
        abcdefx
        """,
        sharey=False,
        width_ratios=[1, 1, 1, 1, 1, 1, 0.05],
        height_ratios=[1, 0.01, 1],
        empty_sentinel="X",
    )

    d = 1
    l = 6
    e = 4

    subsets = [
        "f_first_q",
        "f_second_q",
        "f_third_q",
        "m_first_q",
        "m_second_q",
        "m_third_q",
    ]

    subset_print = {
        "m_third_q": rf"$m, T_{3}$",
        "m_second_q": rf"$m, T_{2}$",
        "m_first_q": rf"$m, T_{1}$",
        "f_third_q": rf"$f, T_{3}$",
        "f_second_q": rf"$f, T_{2}$",
        "f_first_q": rf"$f, T_{1}$",
    }

    pheno_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/age_sex/f_first_q.phen"

    for subset, ax_name in zip(subsets, ["a", "b", "c", "d", "e", "f"]):
        pset = "age_sex"
        outdir = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/{pset}/bdpc_{subset}_d{d}_l{l}_a1e{e}/"
        pag_path = outdir + "all_merged_pag_mk3_ns_lut_correct_n.mtx"
        ace_path = outdir + "all_merged_ACE_ns_lut_no_direction_forced_correct_n.mtx"
        # im = plot_ace_rf_to_d(ace_path, pag_path, pheno_path, title=subset_print[subset], ax=ax_dict[ax_name], cbar=False, norm=mpl.colors.SymLogNorm(vmin=-2.0, vmax=2.0, linthresh=0.01), )
        im = plot_ace_rf_to_d(
            ace_path,
            pag_path,
            pheno_path,
            title=subset_print[subset],
            ax=ax_dict[ax_name],
            cbar=False,
            cmap="Greens",
            vmin=0,
            vmax=0.2,
        )

    # cbar_kw = {"fraction": 0.046, "pad": 0.04, "shrink": 0.5}

    cbar_kw = {"fraction": 2, "pad": 0}

    ax = ax_dict["x"]
    cbar = ax.figure.colorbar(im, ax=ax, orientation="vertical", **cbar_kw)
    cbar.ax.set_ylabel(r"$ACE \: (y_1 \rightarrow y_2)$", rotation=-90, va="bottom")
    ax.axis("off")

    ax = ax_dict["y"]

    d = 1
    l = 6
    e = 4

    row_height = 10

    subsets = [
        "m_third_q",
        "m_second_q",
        "m_first_q",
        "f_third_q",
        "f_second_q",
        "f_first_q",
    ]

    pheno_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/age_sex/f_first_q.phen"
    bim_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.bim"

    traits_with_parents = set()
    for row, subset in enumerate(subsets):
        pset = "age_sex"
        outdir = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/{pset}/bdpc_{subset}_d{d}_l{l}_a1e{e}/"
        traits_with_parents.update(
            set(
                pd.read_csv(
                    outdir + "marker_phen_assoc.csv", sep="\t"
                ).phenotype.unique()
            )
        )
    traits_with_parents = list(traits_with_parents)

    trait_slot_height = row_height / (len(traits_with_parents) + 1)

    colors = []
    cmap = plt.get_cmap("tab20", len(traits_with_parents))
    for i in range(cmap.N):
        rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for row, subset in enumerate(subsets):
        pset = "age_sex"
        outdir = f"/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/{pset}/bdpc_{subset}_d{d}_l{l}_a1e{e}/"
        df = pd.read_csv(outdir + "marker_phen_assoc.csv", sep="\t")
        for trait_ix, trait in enumerate(traits_with_parents):
            sub_df = df[df["phenotype"] == trait]
            x = (
                sub_df.bp.values
                + np.array([global_chr_starts[c] for c in sub_df.chr.values])
            ) / 10**6
            y = np.ones_like(x) * (row * row_height) + trait_slot_height * (
                trait_ix + 1
            )
            ax.plot(x, y, "o", color=colors[trait_ix], alpha=0.8, label=trait)

    handles = []
    for color, trait in zip(colors, traits_with_parents):
        handles.append(mpatches.Patch(color=color, label=trait))

    ax.set_yticks(
        np.array(range(len(subsets))) * row_height + row_height * 0.5,
        [subset_print[s] for s in subsets],
    )
    ax.set_xticks(
        np.array([global_chr_starts[c] + 0.5 * chr_lengths[c] for c in range(1, 23)])
        / 10**6,
        [f"Chr {c}" for c in range(1, 23)],
    )

    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    ax.legend(
        handles=handles,
        loc="upper center",
        fancybox=True,
        shadow=False,
        ncol=len(traits_with_parents) / 2,
        bbox_to_anchor=(0.5, 1.2),
    )

    for x in np.array(list(global_chr_starts.values())) / 10**6:
        ax.axvline(x, color="gray", linestyle=":")

    for row in range(len(subsets)):
        ax.axhline(row * row_height, color="gray", linestyle=":")

    ax_dict["y"].set_title("a)", **title_kw)
    ax_dict["a"].text(-2, -2.5, "b)", size=20, verticalalignment="top")


def plot_est_ukb_marker_positions():
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="tight", figsize=(17, 12))
    ax_dict = fig.subplot_mosaic(
        """
        yyyyyyy
        """,
        sharey=False,
        width_ratios=[1, 1, 1, 1, 1, 1, 0.05],
        empty_sentinel="X",
    )

    d = 1
    l = 6
    e = 4

    subsets = [
        "f_first_q",
        "f_second_q",
        "f_third_q",
        "m_first_q",
        "m_second_q",
        "m_third_q",
    ]

    subset_print = {
        "m_third_q": rf"$m, T_{3}$",
        "m_second_q": rf"$m, T_{2}$",
        "m_first_q": rf"$m, T_{1}$",
        "f_third_q": rf"$f, T_{3}$",
        "f_second_q": rf"$f, T_{2}$",
        "f_first_q": rf"$f, T_{1}$",
    }

    p_names = ["AT", "CAD", "HT", "SMK", "ST", "T2D"]

    ax = ax_dict["y"]

    d = 1
    l = 6
    e = 4

    row_height = 10

    blockfile = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_m11000.blocks"
    bim_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_a1_forced.bim"
    common_out_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/age_sex_split/"

    ukb_assoc = {
        subset: marker_pheno_associations_with_pnames(
            blockfile,
            common_out_path + f"ukb/bdpc_d{d}_l{l}_{subset}_wo_bp_bmi_a1e{e}/",
            p_names,
            bim_path,
        )
        for subset in subsets
    }

    est_assoc = {
        subset: marker_pheno_associations_with_pnames(
            blockfile,
            common_out_path + f"est/bdpc_d{d}_l{l}_{subset}_wo_bp_bmi_a1e{e}/",
            p_names,
            bim_path,
        )
        for subset in subsets
    }

    traits_with_parents = set()

    for assoc in ukb_assoc.values():
        traits_with_parents.update(set(assoc.phenotype.unique()))
    for assoc in est_assoc.values():
        traits_with_parents.update(set(assoc.phenotype.unique()))
    traits_with_parents = list(traits_with_parents)

    trait_slot_height = row_height / (len(traits_with_parents) + 1)

    colors = []
    cmap = plt.get_cmap("tab20", len(traits_with_parents))
    for i in range(cmap.N):
        rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for set_id, assoc_set in zip(["UKB", "EST"], [ukb_assoc, est_assoc]):
        for subset_id, subset in enumerate(subsets):
            row = subset_id * 2 + int(set_id == "EST")
            df = assoc_set[subset]
            for trait_ix, trait in enumerate(traits_with_parents):
                sub_df = df[df["phenotype"] == trait]
                x = (
                    sub_df.bp.values
                    + np.array([global_chr_starts[c] for c in sub_df.chr.values])
                ) / 10**6
                y = np.ones_like(x) * (row * row_height) + trait_slot_height * (
                    trait_ix + 1
                )
                ax.plot(x, y, "o", color=colors[trait_ix], alpha=0.8, label=trait)

    handles = []
    for color, trait in zip(colors, traits_with_parents):
        handles.append(mpatches.Patch(color=color, label=trait))

    ylabels = []
    for subset in subsets:
        ylabels.append(subset_print[subset] + ", UKB")
        ylabels.append(subset_print[subset] + ", EST")

    ax.set_yticks(
        np.array(range(len(ylabels))) * row_height + row_height * 0.5, ylabels
    )
    ax.set_xticks(
        np.array([global_chr_starts[c] + 0.5 * chr_lengths[c] for c in range(1, 23)])
        / 10**6,
        [f"Chr {c}" for c in range(1, 23)],
    )

    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    ax.legend(
        handles=handles,
        loc="upper center",
        fancybox=True,
        shadow=False,
        ncol=len(traits_with_parents) / 2,
        bbox_to_anchor=(0.5, 1.1),
    )

    for x in np.array(list(global_chr_starts.values())) / 10**6:
        ax.axvline(x, color="gray", linestyle=":")

    for row in range(len(ylabels)):
        ax.axhline(row * row_height, color="gray", linestyle=":")

    # ax_dict['y'].set_title("a)", **title_kw)


def plot_ukb_age_sex_marker_positions_ss_vs_no_ss():
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="tight", figsize=(17, 12))
    ax_dict = fig.subplot_mosaic(
        """
        yyyyyyy
        """,
        sharey=False,
        width_ratios=[1, 1, 1, 1, 1, 1, 0.05],
        empty_sentinel="X",
    )

    d = 1
    l = 6
    e = 4

    subsets = [
        "f_first_q",
        "f_second_q",
        "f_third_q",
        "m_first_q",
        "m_second_q",
        "m_third_q",
    ]

    subset_print = {
        "m_third_q": rf"$m, T_{3}$",
        "m_second_q": rf"$m, T_{2}$",
        "m_first_q": rf"$m, T_{1}$",
        "f_third_q": rf"$f, T_{3}$",
        "f_second_q": rf"$f, T_{2}$",
        "f_first_q": rf"$f, T_{1}$",
    }

    p_names = ["AT", "BMI", "CAD", "HT", "SMK", "ST", "T2D"]

    ax = ax_dict["y"]

    d = 1
    l = 6
    e = 4

    row_height = 10

    blockfile = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_m11000.blocks"
    blockfile2 = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.blocks"
    bim_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_a1_forced.bim"
    common_out_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/age_sex_split/"
    cusk_out_path = (
        "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/age_sex/"
    )

    cusk_ss_assoc = {
        subset: marker_pheno_associations_with_pnames(
            blockfile,
            common_out_path + f"ukb/bdpc_d{d}_l{l}_{subset}_wo_bp_a1e{e}/",
            p_names,
            bim_path,
        )
        for subset in subsets
    }

    pheno_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/age_sex/f_first_q.phen"
    cusk_assoc = {
        subset: marker_pheno_associations(
            blockfile2,
            cusk_out_path + f"bdpc_{subset}_d{d}_l{l}_a1e{e}/",
            pheno_path,
            bim_path,
        )
        for subset in subsets
    }

    # traits_with_parents = set()

    # for assoc in cusk_ss_assoc.values():
    #     traits_with_parents.update(set(assoc.phenotype.unique()))
    # for assoc in cusk_assoc.values():
    #     traits_with_parents.update(set(assoc.phenotype.unique()))
    # traits_with_parents = list(traits_with_parents)

    trait_slot_height = row_height / (len(p_names) + 1)

    colors = []
    cmap = plt.get_cmap("tab20", len(p_names))
    for i in range(cmap.N):
        rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for set_id, assoc_set in zip(["cusk-ss", "cusk"], [cusk_ss_assoc, cusk_assoc]):
        for subset_id, subset in enumerate(subsets):
            row = subset_id * 2 + int(set_id == "cusk")
            df = assoc_set[subset]
            for trait_ix, trait in enumerate(p_names):
                sub_df = df[df["phenotype"] == trait]
                x = (
                    sub_df.bp.values
                    + np.array([global_chr_starts[c] for c in sub_df.chr.values])
                ) / 10**6
                y = np.ones_like(x) * (row * row_height) + trait_slot_height * (
                    trait_ix + 1
                )
                ax.plot(x, y, "o", color=colors[trait_ix], alpha=0.8, label=trait)

    handles = []
    for color, trait in zip(colors, p_names):
        handles.append(mpatches.Patch(color=color, label=trait))

    ylabels = []
    for subset in subsets:
        ylabels.append(subset_print[subset] + ", cusk-ss")
        ylabels.append(subset_print[subset] + ", cusk")

    ax.set_yticks(
        np.array(range(len(ylabels))) * row_height + row_height * 0.5, ylabels
    )
    ax.set_xticks(
        np.array([global_chr_starts[c] + 0.5 * chr_lengths[c] for c in range(1, 23)])
        / 10**6,
        [f"Chr {c}" for c in range(1, 23)],
    )

    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    ax.legend(
        handles=handles,
        loc="upper center",
        fancybox=True,
        shadow=False,
        ncol=len(p_names) / 2,
        bbox_to_anchor=(0.5, 1.2),
    )

    for x in np.array(list(global_chr_starts.values())) / 10**6:
        ax.axvline(x, color="gray", linestyle=":")

    for row in range(len(ylabels)):
        ax.axhline(row * row_height, color="gray", linestyle=":")

    ax_dict["y"].set_title("a)", **title_kw)


def plot_est_ukb_full_db_marker_positions():
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="tight", figsize=(17, 12))
    ax_dict = fig.subplot_mosaic(
        """
        yyyyyyy
        """,
        sharey=False,
        width_ratios=[1, 1, 1, 1, 1, 1, 0.05],
        empty_sentinel="X",
    )

    d = 1
    l = 6
    e = 4

    p_names = ["AT", "CAD", "HT", "SMK", "ST", "T2D"]

    ax = ax_dict["y"]

    d = 1
    l = 6
    e = 4

    row_height = 10

    blockfile = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_m11000.blocks"
    bim_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_a1_forced.bim"
    common_out_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/"

    est_assoc = marker_pheno_associations_with_pnames(
        blockfile,
        common_out_path + f"est/bdpc_d{d}_l{l}_a1e{e}_wo_bp_bmi/",
        p_names,
        bim_path,
    )
    ukb_assoc = marker_pheno_associations_with_pnames(
        blockfile,
        common_out_path + f"ukb/bdpc_d{d}_l{l}_a1e{e}_wo_bp_bmi/",
        p_names,
        bim_path,
    )

    trait_slot_height = row_height / (len(p_names) + 1)

    colors = []
    cmap = plt.get_cmap("tab20", len(p_names))
    for i in range(cmap.N):
        rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for row, assoc_set in enumerate([est_assoc, ukb_assoc]):
        df = assoc_set
        for trait_ix, trait in enumerate(p_names):
            sub_df = df[df["phenotype"] == trait]
            x = (
                sub_df.bp.values
                + np.array([global_chr_starts[c] for c in sub_df.chr.values])
            ) / 10**6
            y = np.ones_like(x) * (row * row_height) + trait_slot_height * (
                trait_ix + 1
            )
            ax.plot(x, y, "o", color=colors[trait_ix], alpha=0.8, label=trait)

    handles = []
    for color, trait in zip(colors, p_names):
        handles.append(mpatches.Patch(color=color, label=trait))

    ylabels = ["EST", "UKB"]

    ax.set_yticks(
        np.array(range(len(ylabels))) * row_height + row_height * 0.5, ylabels
    )
    ax.set_xticks(
        np.array([global_chr_starts[c] + 0.5 * chr_lengths[c] for c in range(1, 23)])
        / 10**6,
        [f"Chr {c}" for c in range(1, 23)],
    )

    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    ax.legend(
        handles=handles,
        loc="upper center",
        fancybox=True,
        shadow=False,
        ncol=len(p_names) / 2,
        bbox_to_anchor=(0.5, 1.2),
    )

    for x in np.array(list(global_chr_starts.values())) / 10**6:
        ax.axvline(x, color="gray", linestyle=":")

    for row in range(len(ylabels)):
        ax.axhline(row * row_height, color="gray", linestyle=":")

    ax_dict["y"].set_title("a)", **title_kw)


def plot_ukb_full_db_marker_positions_ss_vs_data():
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="tight", figsize=(17, 8))
    ax_dict = fig.subplot_mosaic(
        """
        yyyyyyy
        """,
        sharey=False,
        width_ratios=[1, 1, 1, 1, 1, 1, 0.05],
        empty_sentinel="X",
    )

    d = 1
    l = 6
    e = 4

    p_names = ["AT", "CAD", "ST", "T2D", "HT", "SMK"]

    ax = ax_dict["y"]

    d = 1
    l = 6
    e = 4

    row_height = 10

    ss_blocks = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_m11000.blocks"
    data_blocks = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection//ukb22828_UKB_EST_v3_ldp08.blocks"
    ss_bim = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.bim"
    data_bim = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb22828_UKB_EST_v3_ldp08_estonia_intersect_a1_forced.bim"
    pss_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/"

    data_assoc = marker_pheno_associations_with_pnames(
        data_blocks,
        pss_path + f"production/cusk_d{d}_l{l}_a1e{e}_est_intersect_wo_bp_bmi/",
        p_names,
        data_bim,
    )
    ss_assoc = marker_pheno_associations_with_pnames(
        ss_blocks,
        pss_path + f"estonian_comparison/ukb/bdpc_d{d}_l{l}_a1e{e}_wo_bp_bmi/",
        p_names,
        ss_bim,
    )

    trait_slot_height = row_height / (len(p_names) + 1)

    colors = []
    cmap = plt.get_cmap("tab20", len(p_names))
    for i in range(cmap.N):
        rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for row, assoc_set in enumerate([data_assoc, ss_assoc]):
        df = assoc_set
        for trait_ix, trait in enumerate(p_names):
            sub_df = df[df["phenotype"] == trait]
            x = (
                sub_df.bp.values
                + np.array([global_chr_starts[c] for c in sub_df.chr.values])
            ) / 10**6
            y = np.ones_like(x) * (row * row_height) + trait_slot_height * (
                trait_ix + 1
            )
            ax.plot(x, y, "o", color=colors[trait_ix], alpha=0.8, label=trait)

    handles = []
    for color, trait in zip(colors, p_names):
        handles.append(mpatches.Patch(color=color, label=trait))

    ylabels = ["raw-data", "sum-stat"]

    ax.set_yticks(
        np.array(range(len(ylabels))) * row_height + row_height * 0.5, ylabels
    )
    ax.set_xticks(
        np.array([global_chr_starts[c] + 0.5 * chr_lengths[c] for c in range(1, 23)])
        / 10**6,
        [f"Chr {c}" for c in range(1, 23)],
    )

    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    ax.legend(
        handles=handles,
        loc="upper center",
        fancybox=True,
        shadow=False,
        ncol=len(p_names) / 2,
        bbox_to_anchor=(0.5, 1.2),
    )

    for x in np.array(list(global_chr_starts.values())) / 10**6:
        ax.axvline(x, color="gray", linestyle=":")

    for row in range(len(ylabels)):
        ax.axhline(row * row_height, color="gray", linestyle=":")

    ax_dict["y"].set_title("a)", **title_kw)


def plot_ukb_full_db_marker_positions_6_vs_17():
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="tight", figsize=(17, 8))
    ax_dict = fig.subplot_mosaic(
        """
        yyyyyyy
        """,
        sharey=False,
        width_ratios=[1, 1, 1, 1, 1, 1, 0.05],
        empty_sentinel="X",
    )

    d = 1
    l = 6
    e = 4

    # small_set_names = ["AT", "CAD", "HT", "SMK", "ST", "T2D"]
    # full_set_names = [
    #     "WHR",
    #     "BMI",
    #     "SBP",
    #     "DBP",
    #     "HbA1c",
    #     "GLU",
    #     "CHOL",
    #     "TRIG",
    #     "HDL",
    #     "LDL",
    #     "AT",
    #     "CAD",
    #     "ST",
    #     "T2D",
    #     "HT",
    #     "ALC",
    #     "SMK"]

    ax = ax_dict["y"]

    d = 1
    l = 6
    e = 4

    row_height = 10

    small_set_pheno = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/est_intersect_wo_bp_bmi.phen"
    full_set_pheno = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/input.phen"
    blocks = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection//ukb22828_UKB_EST_v3_ldp08.blocks"
    bim = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.bim"
    pss_path = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/"

    small_set_assoc = marker_pheno_associations(
        blocks,
        pss_path + f"production/cusk_d{d}_l{l}_a1e{e}_est_intersect_wo_bp_bmi/",
        small_set_pheno,
        bim,
    )
    full_set_assoc = marker_pheno_associations(
        blocks, pss_path + f"production/bdpc_d{d}_l{l}_a1e{e}/", full_set_pheno, bim
    )
    p_names = list(
        set(full_set_assoc.phenotype.unique()) | set(small_set_assoc.phenotype.unique())
    )

    trait_slot_height = row_height / (len(p_names) + 1)

    colors = []
    cmap = plt.get_cmap("tab20", len(p_names))
    for i in range(cmap.N):
        rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for row, assoc_set in enumerate([small_set_assoc, full_set_assoc]):
        df = assoc_set
        for trait_ix, trait in enumerate(p_names):
            sub_df = df[df["phenotype"] == trait]
            x = (
                sub_df.bp.values
                + np.array([global_chr_starts[c] for c in sub_df.chr.values])
            ) / 10**6
            y = np.ones_like(x) * (row * row_height) + trait_slot_height * (
                trait_ix + 1
            )
            ax.plot(x, y, "o", color=colors[trait_ix], alpha=0.8, label=trait)

    handles = []
    for color, trait in zip(colors, p_names):
        handles.append(mpatches.Patch(color=color, label=trait))

    ylabels = ["6 traits", "17 traits"]

    ax.set_yticks(
        np.array(range(len(ylabels))) * row_height + row_height * 0.5, ylabels
    )
    ax.set_xticks(
        np.array([global_chr_starts[c] + 0.5 * chr_lengths[c] for c in range(1, 23)])
        / 10**6,
        [f"Chr {c}" for c in range(1, 23)],
    )

    plt.setp(
        ax.get_xticklabels(),
        rotation=XLABEL_ROTATION,
        ha="right",
        rotation_mode="anchor",
    )

    ax.legend(
        handles=handles,
        loc="upper center",
        fancybox=True,
        shadow=False,
        ncol=len(p_names) / 2,
        bbox_to_anchor=(0.5, 1.2),
    )

    for x in np.array(list(global_chr_starts.values())) / 10**6:
        ax.axvline(x, color="gray", linestyle=":")

    for row in range(len(ylabels)):
        ax.axhline(row * row_height, color="gray", linestyle=":")

    ax_dict["y"].set_title("a)", **title_kw)


def average_corrs():
    wdir = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/production/bdpc_d1_l6_a1e4/"

    with open(wdir + "all_merged.mdim", "r") as fin:
        line = next(fin)
        fields = line.split()
        num_var = int(fields[0])
        num_phen = int(fields[1])
        num_markers = num_var - num_phen

    corr = bdpc.mmread(wdir + "all_merged_scm.mtx")
    adj = bdpc.mmread(wdir + "all_merged_sam.mtx")

    corr = corr.toarray()
    adj = adj.toarray()
    pxp_corr = corr[:num_phen, :num_phen]
    pxp_adj = adj[:num_phen, :num_phen]
    mxp_corr = corr[:num_phen, num_phen:]
    mxp_adj = adj[:num_phen, num_phen:]
    mxm_corr = corr[num_phen:, num_phen:]
    mxm_adj = adj[num_phen:, num_phen:]

    print(
        "pxp corr mean, var: ",
        np.mean(np.abs(pxp_corr[pxp_adj == 1])),
        np.var(np.abs(pxp_corr[pxp_adj == 1])),
    )
    print(
        "mxp corr mean, var: ",
        np.mean(np.abs(mxp_corr[mxp_adj == 1])),
        np.var(np.abs(mxp_corr[mxp_adj == 1])),
    )


def plot_compare_correlations_matrices():
    sumstat_corrs = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/ukb/bdpc_d1_l6_a1e4_wo_bp_bmi/3_0_3967.all_corrs"
    rawdata_corrs = "/nfs/scistore17/robingrp/human_data/causality/parent_set_selection/estonian_comparison/bdpc_d1_l6_a1e4/3_0_3967.all_corrs"

    ssm = np.fromfile(sumstat_corrs, dtype=np.float32).reshape(3974, 3974)
    rdm = np.fromfile(rawdata_corrs, dtype=np.float32).reshape(3974, 3974)

    num_markers = 3968
    ssm_mxm = ssm[:num_markers, :num_markers]
    rdm_mxm = rdm[:num_markers, :num_markers]

    num_phen = 6
    ssm_mxp = ssm[:num_markers, num_markers:]
    rdm_mxp = rdm[:num_markers, num_markers:]

    ssm_pxp = ssm[num_markers:, num_markers:]
    rdm_pxp = rdm[num_markers:, num_markers:]

    fig = plt.figure(layout="constrained", figsize=(15, 5))
    ax_dict = fig.subplot_mosaic(
        """
        abc
        """
    )

    ax = ax_dict["a"]
    diag = np.linspace(
        np.min([ssm_mxm.min(), rdm_mxm.min()]),
        np.max([ssm_mxm.max(), rdm_mxm.max()]),
        10,
    )
    ax.plot(diag, diag, ":", color="gray")
    ax.plot(ssm_mxm.flatten(), rdm_mxm.flatten(), ".", rasterized=True)
    ax.set_title("MxM")
    ax.set_ylabel("sumstat")
    ax.set_xlabel("rawdata")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_box_aspect(1)

    ax = ax_dict["b"]
    diag = np.linspace(
        np.min([ssm_mxp.min(), rdm_mxp.min()]),
        np.max([ssm_mxp.max(), rdm_mxp.max()]),
        10,
    )
    ax.plot(diag, diag, ":", color="gray")
    ax.plot(ssm_mxp.flatten(), rdm_mxp.flatten(), ".", rasterized=True)
    ax.set_title("MxP")
    ax.set_ylabel("sumstat")
    ax.set_xlabel("rawdata")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_box_aspect(1)

    ax = ax_dict["c"]
    diag = np.linspace(
        np.min([ssm_pxp.min(), rdm_pxp.min()]),
        np.max([ssm_pxp.max(), rdm_pxp.max()]),
        10,
    )
    ax.plot(diag, diag, ":", color="gray")
    ax.plot(ssm_pxp.flatten(), rdm_pxp.flatten(), ".", rasterized=True)
    ax.set_title("PxP")
    ax.set_ylabel("sumstat")
    ax.set_xlabel("rawdata")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_box_aspect(1)
