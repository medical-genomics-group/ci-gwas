import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from dataclasses import dataclass
import json
import pandas as pd
import seaborn as sns
from dataclasses import dataclass
import queue
from scipy.io import mmread
import json

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

rng = np.random.default_rng()

BASE_INDEX = 1


def get_pheno_codes(phen_path):
    with open(phen_path, 'r') as fin:
        header = fin.readline()
    return header.strip().split("\t")[2:]


def get_block_out_stems(blockpath):
    res = []
    with open(blockpath, 'r') as fin:
        for line in fin:
            line = line.strip()
            fields = line.split("\t")
            res.append(f"{fields[0]}_{fields[1]}_{fields[2]}")
    return res


def load_mdim(basepath: str):
    with open(basepath + '.mdim', 'r') as fin:
        fields = fin.readline().strip().split("\t")
    return [int(f) for f in fields]


def make_dm_ix_to_sm_ix(num_m: int, num_p: int,
                        marker_offset: int) -> np.array:
    ixs = np.arange(num_m + num_p)
    new_ixs = np.zeros_like(ixs)
    new_ixs[np.where(ixs < num_m)] = ixs[np.where(
        ixs < num_m)] + marker_offset + num_p + BASE_INDEX
    new_ixs[np.where(
        ixs >= num_m)] = ixs[np.where(ixs >= num_m)] - num_m + BASE_INDEX
    return new_ixs


def load_mat_sparse(basepath: str, num_m: int, num_p: int, marker_offset,
                    dtype, suffix):
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


def load_sepset_sparse(basepath: str, num_m: int, num_p: int, max_level: int,
                       marker_offset):
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
                    raise ValueError('SepSet(x, y) contains x or y')
                res[(dm2sm[i], dm2sm[j])] = new_s

    return res


def load_corr_sparse(basepath: str, num_m: int, num_p: int, marker_offset):
    return load_mat_sparse(basepath, num_m, num_p, marker_offset, np.float32,
                           ".corr")


def load_adj_sparse(basepath: str, num_m: int, num_p: int, marker_offset=0):
    return load_mat_sparse(basepath, num_m, num_p, marker_offset, np.int32,
                           ".adj")


def load_global_marker_indices(basepath: str,
                               num_m: int,
                               num_p: int,
                               selected_marker_offset=0,
                               global_marker_offset=0):
    global_marker_indices = {}
    rel_to_block = np.fromfile(basepath + ".ixs", dtype=np.int32)
    dm2sm = make_dm_ix_to_sm_ix(num_m, num_p, selected_marker_offset)
    for dm_ix in range(len(dm2sm)):
        sm_ix = dm2sm[dm_ix]
        if sm_ix > num_p:
            global_marker_indices[
                sm_ix] = rel_to_block[dm_ix] + global_marker_offset
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
    a.update(b)


def merge_block_outputs(blockfile: str, outdir: str):
    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    bo = BlockOutput(basepaths[0])
    marker_offset = bo.num_markers()
    global_marker_offset = bo.block_size()

    sam = bo.sam()
    scm = bo.scm()
    ssm = bo.ssm()
    gmi = bo.gmi()

    for path in basepaths[1:]:
        bo = BlockOutput(path, marker_offset, global_marker_offset)
        add_sam(sam, bo.sam(), bo.num_phen())
        add_scm(scm, bo.scm())
        add_ssm(ssm, bo.ssm())
        add_gmi(gmi, bo.gmi())
        marker_offset += bo.num_markers()
        global_marker_offset += bo.block_size()

    return GlobalBdpcResult(sam, scm, ssm, gmi, marker_offset + bo.num_phen(),
                            bo.num_phen(), bo.max_level())


def global_epm(blockfile: str, outdir: str):
    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    bo = BlockOutput(basepaths[0])
    marker_offset = bo.num_markers()

    epm = bo.exclusive_pleiotropy_mat()

    for path in basepaths[1:]:
        bo = BlockOutput(path, marker_offset)
        marker_offset += bo.num_markers()
        for k, v in bo.exclusive_pleiotropy_mat().items():
            if k in epm:
                epm[k] += v
            else:
                epm[k] = v
    return epm


def global_eps(blockfile: str, outdir: str):
    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    bo = BlockOutput(basepaths[0])
    marker_offset = bo.num_markers()

    eps = bo.exclusive_pleiotropy_sets()

    for path in basepaths[1:]:
        bo = BlockOutput(path, marker_offset)
        marker_offset += bo.num_markers()
        for k, v in bo.exclusive_pleiotropy_sets().items():
            if k in eps:
                eps[k].update(v)
            else:
                eps[k] = v
    return eps


def pleitropy_mat(pag, num_phen):
    """Compute upper bound of markers that could affect each phenotype or combination of phenotypes
    """
    res = {}
    plr = pag.tolil().rows
    phens = set(np.arange(num_phen))
    for pix in range(num_phen):
        visited = set()
        q = queue.Queue()
        q.put(pix)
        while not q.empty():
            v1 = q.get()
            for v2 in plr[v1]:
                if not v2 in phens and not v2 in visited and pag[v2, v1] == 2:
                    q.put(v2)
                    visited.add(v2)
        res[pix] = visited
    return res


def exclusive_pleiotropy_sets(pag, num_phen):
    pm = pleitropy_mat(pag, num_phen)
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


class BlockOutput:

    def __init__(self, basepath: str, marker_offset=0, global_marker_offset=0):
        self.basepath = basepath
        self.mdim = load_mdim(basepath)
        # number of selected markers in all previous blocks
        self.marker_offset = marker_offset
        # .bim row index of first marker in block definition
        self.global_marker_offset = global_marker_offset

    def block_size(self) -> int:
        first, last = self.basepath.split("_")[-2:]
        return int(last) - int(first) + 1

    def max_level(self) -> int:
        return self.mdim[2]

    def num_markers(self) -> int:
        return self.mdim[0] - self.mdim[1]

    def num_phen(self) -> int:
        return self.mdim[1]

    def has_markers(self) -> bool:
        return self.num_markers() > 0

    def sam(self):
        return load_adj_sparse(self.basepath, self.num_markers(),
                               self.num_phen(), self.marker_offset)

    def scm(self):
        return load_corr_sparse(self.basepath, self.num_markers(),
                                self.num_phen(), self.marker_offset)

    def ssm(self):
        return load_sepset_sparse(self.basepath, self.num_markers(),
                                  self.num_phen(), self.max_level(),
                                  self.marker_offset)

    def gmi(self):
        return load_global_marker_indices(self.basepath, self.num_markers(),
                                          self.num_phen(), self.marker_offset,
                                          self.global_marker_offset)

    def marker_indices(self):
        first = self.num_phen() + self.marker_offset
        last = first + self.num_markers()
        return np.arange(first, last) + BASE_INDEX

    def pheno_indices(self):
        return np.arange(0, self.num_phen()) + BASE_INDEX

    def pleitropy_mat(self):
        """Compute upper bound of markers that could affect each phenotype or combination of phenotypes
        """
        res = {}
        adj = self.sam()
        phens = set(self.pheno_indices())
        for pix in self.pheno_indices():
            visited = set()
            q = queue.Queue()
            q.put(pix)
            while not q.empty():
                v1 = q.get()
                for v2 in self.marker_indices():
                    if (v1, v2) in adj and not v2 in visited:
                        q.put(v2)
                        visited.add(v2)
            res[pix] = visited
        return res

    def exclusive_pleiotropy_mat(self):
        pm = self.pleitropy_mat()
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

    def exclusive_pleiotropy_sets(self):
        pm = self.pleitropy_mat()
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
            N = M = len(set([t[0] for t in self.sam.keys()]))
            fout.write("%%MatrixMarket matrix coordinate integer general\n")
            fout.write(f"{N}\t{M}\t{L}\n")
            for (t1, t2), v in self.sam.items():
                fout.write(f"{t1}\t{t2}\t{v}\n")

        with open(basepath + "_scm" + ".mtx", "w") as fout:
            L = len(self.scm)
            N = M = len(set([t[0] for t in self.scm.keys()]))
            fout.write("%%MatrixMarket matrix coordinate real general\n")
            fout.write(f"{N}\t{M}\t{L}\n")
            for (t1, t2), v in self.scm.items():
                fout.write(f"{t1}\t{t2}\t{v}\n")

        with open(basepath + ".ssm", "w") as fout:
            for (t1, t2), v in self.ssm.items():
                outline = " ".join([str(e) for e in [t1, t2] + v])
                fout.write(outline + "\n")

        with open(basepath + ".mdim", "w") as fout:
            fout.write(f"{self.num_var}\t{self.num_phen}\t{self.max_level}\n")


def combine_all_pheno_and_plot():
    outdir = "/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/bdpc_d1_l6_a1e8/"
    blockfile = "/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.blocks"
    p_names = get_pheno_codes(
        "/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/input.phen"
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
        pd_gepm.append({
            "y1": p_names[i - 1],
            "y2": p_names[j - 1],
            "count": c
        })
    pd_gepm = pd.DataFrame(pd_gepm)

    pd_pcm = []
    for i in range(1, num_phen + 1):
        for j in range(1, num_phen + 1):
            pd_pcm.append({
                "y1": p_names[i - 1],
                "y2": p_names[j - 1],
                "v": gr.scm[(i, j)]
            })
    pd_pcm = pd.DataFrame(pd_pcm)

    plt.figure(figsize=(20, 15))
    g = pd_pcm.pivot("y1", "y2", "v")
    mask = np.triu(np.ones_like(g, dtype=bool))
    np.fill_diagonal(mask, False)
    sns.heatmap(g,
                xticklabels=1,
                yticklabels=1,
                cmap="RdBu",
                annot=True,
                square=True,
                mask=mask,
                vmax=1,
                vmin=-1)

    pd_pam = []
    for i in range(1, num_phen + 1):
        for j in range(1, num_phen + 1):
            if (i, j) in gr.sam:
                pd_pam.append({
                    "y1": p_names[i - 1],
                    "y2": p_names[j - 1],
                    "v": 1
                })
            else:
                pd_pam.append({
                    "y1": p_names[i - 1],
                    "y2": p_names[j - 1],
                    "v": 0
                })
    pd_pam = pd.DataFrame(pd_pam)

    plt.figure(figsize=(14, 14))
    g = pd_pam.pivot("y1", "y2", "v")
    mask = np.triu(np.ones_like(g, dtype=bool))
    np.fill_diagonal(mask, False)
    sns.heatmap(g,
                xticklabels=1,
                yticklabels=1,
                cbar=False,
                square=True,
                mask=mask,
                cmap="Blues",
                vmax=2)

    num_phen = len(p_names)
    np_gepm = np.zeros(shape=(num_phen, num_phen))
    pd_gepm = []
    for (i, j), c in {k: len(v) for k, v in geps.items()}.items():
        np_gepm[i - 1, j - 1] = c
        pd_gepm.append({
            "y1": p_names[i - 1],
            "y2": p_names[j - 1],
            "count": c
        })
    pd_gepm = pd.DataFrame(pd_gepm)

    plt.figure(figsize=(20, 15))
    g = pd_gepm.pivot("y1", "y2", "count")
    mask = np.triu(np.ones_like(g, dtype=bool))
    np.fill_diagonal(mask, False)
    sns.heatmap(g,
                xticklabels=1,
                yticklabels=1,
                norm=matplotlib.colors.LogNorm(),
                cmap="Blues",
                annot=True,
                square=True,
                mask=mask,
                cbar_kws={'label': '# parent markers'})

    pag_path = "/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/bdpc_d1_l6_a1e8/all_merged_pag.mtx"

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
    sns.heatmap(g,
                xticklabels=1,
                yticklabels=1,
                cmap="Blues",
                annot=False,
                square=True,
                cbar=False)

    num_phen = len(p_names)
    geps = exclusive_pleiotropy_sets(pag, num_phen)

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
    sns.heatmap(g,
                xticklabels=1,
                yticklabels=1,
                norm=matplotlib.colors.LogNorm(),
                cmap="Blues",
                annot=True,
                square=True,
                mask=mask,
                cbar_kws={'label': '# parent markers'})
