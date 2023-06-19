import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from dataclasses import dataclass
import json
import pandas as pd
import seaborn as sns
import queue
from scipy.io import mmread
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
    new_ixs[np.where(ixs < num_m)] = ixs[np.where(ixs < num_m)] + marker_offset + num_p + BASE_INDEX
    new_ixs[np.where(ixs >= num_m)] = ixs[np.where(ixs >= num_m)] - num_m + BASE_INDEX
    return new_ixs


def load_mat_sparse(basepath: str, num_m: int, num_p: int, marker_offset, dtype, suffix):
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

def load_sepset_sparse(basepath: str, num_m: int, num_p: int, max_level: int, marker_offset):
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
                res[(dm2sm[i], dm2sm[j])] = new_s

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

    for (i, j) in zip(dm2sm[not_zero[0]], dm2sm[not_zero[1]]):
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
    basepath: str, num_m: int, num_p: int, selected_marker_offset=0, global_marker_offset=0
):
    global_marker_indices = {}
    rel_to_block = np.fromfile(basepath + ".ixs", dtype=np.int32)
    dm2sm = make_dm_ix_to_sm_ix(num_m, num_p, selected_marker_offset)
    for dm_ix in range(len(dm2sm)):
        sm_ix = dm2sm[dm_ix]
        if sm_ix >= num_p:
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
        try:
            bo = BlockOutput(path, marker_offset, global_marker_offset)
        except FileNotFoundError:
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
        bo = BlockOutput(path, marker_offset)
        marker_offset += bo.num_markers()
        for k, v in bo.pheno_ancestor_sets(depth).items():
            if not reduced_indices:
                v = {frozenset({gr.gmi[ix] for ix in s}) for s in v}
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
                    if not v2 in phens and not v2 in visited and neighbor_fn(pag, v1, v2):
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
        return load_skeleton(self.basepath, self.num_markers(), self.num_phen(), self.marker_offset)

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
            self.basepath, self.num_markers(), self.num_phen(), self.max_level(), self.marker_offset
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
                    anc_set = set()
                    anc_set.add(parent)
                    q = queue.Queue()
                    q.put(parent)
                    next_q = queue.Queue()
                    for _ in range(depth - 1):
                        while not q.empty():
                            v1 = q.get()
                            for v2 in self.marker_indices():
                                if (v1, v2) in adj and not v2 in anc_set:
                                    next_q.put(v2)
                                    anc_set.add(v2)
                    res[pix].add(frozenset(anc_set))
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
        res = {(i, i): set() for i in self.pheno_indices()}
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
    title_kw=dict(),
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
        cm.set_bad("white")
        kwargs["cmap"] = cm

    # ax.grid(which="major", color="gray", linestyle=":", linewidth=0.5)

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    if cbar:
        # Create colorbar
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    else:
        cbar = None

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=50, ha="right", rotation_mode="anchor")

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
    mat_type="exclusive"
):
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
    **kwargs
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
        **kwargs
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
    poss_parents = pag_exclusive_pleiotropy_sets(pag_path, pheno_path, neighbor_fn, depth)
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


def load_ace(ace_path: str, pheno_path: str) -> np.array:
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    ace = mmread(ace_path).tocsr()
    z = [[0.0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(num_phen):
            z[i][j] = ace[i, j]
    return np.array(z)


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

six_edge_types = EdgeEncoding(
    [
        r"$y_1 \; \; \; y_2$",
        r"$y_1$ <-o $y_2$",
        r"$y_1$ o-> $y_2$",
        r"$y_1$ <-> $y_2$",
        r"$y_1$ -> $y_2$",
        r"$y_1$ <- $y_2$",
        r"$y_1$ - $y_2$",
    ],
    {
        (0, 0): 0,
        (1, 2): 1,
        (2, 1): 2,
        (2, 2): 3,
        (2, 3): 4,
        (3, 2): 5,
        (3, 3): 6,
    },
    mpl.colors.ListedColormap(
        np.array(
            ["#ffffff", "#000421", "#46204c", "#99305b", "#dd584b", "#fb9d20", "#e5ed05"]  # white
        )
    ),
)

six_edge_types = EdgeEncoding(
    [
        r"$y_1 \; \; \; y_2$",
        r"$y_1$ <-o $y_2$",
        r"$y_1$ o-> $y_2$",
        r"$y_1$ <-> $y_2$",
        r"$y_1$ -> $y_2$",
        r"$y_1$ <- $y_2$",
        r"$y_1$ - $y_2$",
    ],
    {
        (0, 0): 0,
        (1, 2): 1,
        (2, 1): 2,
        (2, 2): 3,
        (2, 3): 4,
        (3, 2): 5,
        (3, 3): 6,
    },
    mpl.colors.ListedColormap(
        np.array(
            ["#ffffff", "#000421", "#46204c", "#99305b", "#dd584b", "#fb9d20", "#e5ed05"]  # white
        )
    ),
)

five_common_edge_types = EdgeEncoding(
    [
        r"$y_1 \; \; \; y_2$",
        r"$y_1$ <-o $y_2$",
        r"$y_1$ o-> $y_2$",
        r"$y_1$ <-> $y_2$",
        r"$y_1$ -> $y_2$",
        r"$y_1$ <- $y_2$",
    ],
    {
        (0, 0): 0,
        (1, 2): 1,
        (2, 1): 2,
        (2, 2): 3,
        (2, 3): 4,
        (3, 2): 5,
    },
    mpl.colors.ListedColormap(
        np.array(
            [
                "#ffffff",  # white
                "#000421",
                "#5b2453",
                "#bf4056",
                "#f88a2e",
                "#e5ed05",
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
                "#b2df8a",
                "#a6cee3",
                "#1f78b4",
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
    cbar=True
):

    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)

    pag = mmread(pag_path).tocsr()

    z = [[0 for _ in range(num_phen)] for _ in range(num_phen)]

    for i in range(num_phen):
        for j in range(i):
            v1 = pag[i, j]
            v2 = pag[j, i]
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
        p_names,
        p_names,
        cmap=edge_encoding.cmap,
        norm=norm,
        cbar_kw=cbar_kw,
        # cbarlabel="Edge Type",
        xlabel=r"$y_2$",
        ylabel=r"$y_1$",
        title=title,
        title_kw=title_kw,
        ax=ax,
        cbar=cbar
    )


def marker_pheno_associations(
    blockfile: str,
    outdir: str,
    pag_path: str,
    pheno_path: str,
    bim_path: str,
    neighbor_fn=is_possible_child,
    depth=1,
):
    gr = merge_block_outputs(blockfile, outdir)
    rs_ids = pd.read_csv(bim_path, sep="\t", header=None)[1].values
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    geps = pag_exclusive_pleiotropy_sets(pag_path, pheno_path, neighbor_fn, depth)
    assoc_markers = []

    for pix in range(num_phen):
        for v in geps[(pix, pix)]:
            # pag is 0-based indexed, everything else uses 1-based indexing
            bim_line = int(gr.gmi[v + BASE_INDEX])
            assoc_markers.append(
                {"phenotype": p_names[pix], "rsID": rs_ids[bim_line], "bim_line_ix": bim_line}
            )

    return pd.DataFrame(assoc_markers)


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
        pd_gepm.append({"y1": p_names[i - 1], "y2": p_names[j - 1], "count": c})
    pd_gepm = pd.DataFrame(pd_gepm)

    pd_pcm = []
    for i in range(1, num_phen + 1):
        for j in range(1, num_phen + 1):
            pd_pcm.append({"y1": p_names[i - 1], "y2": p_names[j - 1], "v": gr.scm[(i, j)]})
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
        g, xticklabels=1, yticklabels=1, cbar=False, square=True, mask=mask, cmap="Blues", vmax=2
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
    sns.heatmap(g, xticklabels=1, yticklabels=1, cmap="Blues", annot=False, square=True, cbar=False)

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

def load_simulation_results() -> pd.DataFrame:
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

    def calulate_performance_metrics(
        true_adj: np.array,
        true_eff: np.array,
        est_adj: np.array,
        est_eff: np.array
    ):
        true_adj_masked = np.ma.array(true_adj, mask=np.tri(true_adj.shape[0], k=0))
        est_adj_masked = np.ma.array(est_adj, mask=np.tri(true_adj.shape[0], k=0))
        est_eff_masked = np.ma.array(est_eff, mask=np.tri(true_adj.shape[0], k=0))
        true_eff_masked = np.ma.array(true_eff, mask=np.tri(true_adj.shape[0], k=0))
        sig_eff = np.zeros_like(est_eff)
        sig_eff[est_adj != 0] = est_eff[est_adj != 0]
        tp_mat = (true_adj_masked != 0) & (est_adj_masked != 0)
        p = np.sum(true_adj_masked != 0)
        f = np.sum(true_adj_masked == 0)
        tp = np.sum((true_adj_masked != 0) & (est_adj_masked != 0))
        fp = np.sum((true_adj_masked == 0) & (est_adj_masked != 0))
        fn = np.sum((true_adj_masked != 0) & (est_adj_masked == 0))
        mse = np.sum((true_eff - sig_eff) ** 2)
        bias = np.sum(true_eff - sig_eff)
        var = np.var(true_eff - sig_eff)
        mse_tp = np.sum((est_eff_masked[tp_mat] - true_eff_masked[tp_mat]) ** 2)
        bias_tp = np.sum(est_eff_masked[tp_mat] - true_eff_masked[tp_mat])
        var_tp = np.var(est_eff_masked[tp_mat] - true_eff_masked[tp_mat])
        fdr = fp / (fp + tp)
        tpr = tp / p
        if np.ma.is_masked(mse_tp):
            mse_tp = np.nan
        if np.ma.is_masked(var_tp):
            var_tp = np.nan
        if np.ma.is_masked(bias_tp):
            bias_tp = np.nan
        return Performance(mse, bias, var, mse_tp, bias_tp, var_tp, fdr, tpr)

    @dataclass
    class MR:
        method: str
        exposure: str
        outcome: str
        p: str
        estimate: str

    mr_cn = [
        MR("ivw", "ivw.Exposure", "ivw.Outcome", "ivw.p", "ivw.est"),
        MR("egger", "egger.Exposure", "egger.Outcome", "egger.p", "egger.est"),
        MR("mrpresso", "mrpresso.V1", "mrpresso.V2", "mrpresso.P.value", "mrpresso.Causal.Estimate"),
        MR("cause", "CAUSE.V1", "CAUSE.V2", "CAUSE.V4", "CAUSE.gamma"),
    ]

    pdir = f"/nfs/scistore13/robingrp/human_data/causality/bias_as_fn_of_alpha/sim_small_effects/"

    d = 1
    l = 6
    n_arr = [2000, 4000, 8000, 16000]
    m_arr = [200, 400, 800, 1600]
    # m_arr = [200, 400, 1600]
    e_arr = list(range(1, 9))
    rep_arr = list(range(1, 21))
    num_phen = 10

    rows = []

    for (n, m, rep) in itertools.product(n_arr, m_arr, rep_arr):

        # load true dag
        dag_path = pdir + f"./true_adj_mat_n{n}_SNP_{m}_it_{rep}.mtx"
        dag = mmread(dag_path).tocsr()
        pdag = dag[-(num_phen):,-(num_phen):].toarray()

        # load true causal effects
        eff_path = pdir + f"./True_causaleffect_n{n}_SNP_{m}_it_{rep}.mtx"
        eff = mmread(dag_path).tocsr()
        peff = eff[-(num_phen):,-(num_phen):].toarray()

        for mr in mr_cn:
            mr_results = pd.read_csv(pdir + f"mr_res_{mr.method}_n{n}_SNP_{m}_it_{rep}.csv")
            mr_results['i'] = mr_results[mr.exposure].apply(lambda x: int(x.split("Y")[1]) - 1)
            mr_results['j'] = mr_results[mr.outcome].apply(lambda x: int(x.split("Y")[1]) - 1)
            pvals = np.ones(shape=(num_phen, num_phen))
            pvals[mr_results['i'], mr_results['j']] = mr_results[mr.p]
            effects = np.zeros(shape=(num_phen, num_phen))
            effects[mr_results['i'], mr_results['j']] = mr_results[mr.estimate]
            for e in e_arr:
                adj = pvals < 10**(-e)
                perf = calulate_performance_metrics(pdag, peff, adj, effects)
                rows.append({
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
                    "alpha": 10**(-e),
                    "method": mr.method,
                })

        for e in e_arr:
            indir = pdir + f"simpc_d{d}_l{l}_e{e}_i{rep}_n{n}_m{m}/"

            mdim_path = indir + "skeleton.mdim"
            with open(mdim_path, 'r') as fin:
                num_var, num_phen, max_level = [int(elem) for elem in fin.readline().split()]

            # load skeleton
            adj = np.fromfile(indir + "skeleton.adj", dtype=np.int32).reshape(num_var, num_var)
            padj = adj[-(num_phen):,-(num_phen):]

            # load var indices
            # var_ixs = np.fromfile(indir + "skeleton.ixs", dtype=np.int32)

            # load aces
            pace = np.zeros(shape=(num_phen, num_phen))
            missing = []
            for i in range(1, num_phen + 1):
                for j in range(1, num_phen + 1):
                    if i == j:
                        continue
                    file = indir + f"estimated_causaleffect_i{i}_j{j}.csv"
                    try:
                        with open(file) as fin:
                            sym = fin.readline().strip()
                            if sym == "NaN":
                                pace[i - 1, j - 1] = 0.0
                            else:
                                pace[i - 1, j - 1] = eval(sym)
                    except FileNotFoundError as e:
                        missing.append((i, j))

            perf = calulate_performance_metrics(pdag, peff, padj, pace)
            rows.append({
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
                "alpha": 10**(-e),
                "method": "ci-gwas",
            })

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

    df = load_simulation_results()
    dfs = df.loc[(df['n'] == 16000) & (df['m'] == 1600)]
    gr = dfs.groupby(["method", "alpha"])
    means = gr.mean()
    stds = gr.std()

    alphas = 10.0 ** -np.array(e_arr[::-1])
    # methods = ['ci-gwas', 'cause', 'mrpresso', 'egger', 'ivw']
    methods = ['ci-gwas', 'cause', 'mrpresso', 'ivw']
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
                bars = ax.bar(x + offset, mu, width, yerr=sig, capsize=1.5, label=method)
                handles.append(bars)
            except StopIteration:
                pass
            multiplier += 1

        ax.set_ylabel(metric)
        if title == 'e)' or title == 'd)':
            ax.set_xlabel(r"$\alpha$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks(x + width, x_vals)
        plt.setp(ax.get_xticklabels(), rotation=50, ha="right", rotation_mode="anchor")
        # ax.legend(loc='upper left', ncols=3)
        # ax.set_yscale("symlog")
        return handles

    fig = plt.figure(layout="constrained", figsize=(10, 8))
    ax_dict = fig.subplot_mosaic(
        """
        fa
        bc
        de
        """,
        empty_sentinel="X",
        sharex=True
    )

    plot_bars(alphas, "fdr", means, stds, ax_dict['a'], "a)")
    plot_bars(alphas, "tpr", means, stds, ax_dict['b'], "b)")
    plot_bars(alphas, "mse", means, stds, ax_dict['c'], "c)")
    plot_bars(alphas, "bias", means, stds, ax_dict['d'], "d)")
    h = plot_bars(alphas, "var", means, stds, ax_dict['e'], "e)")
    ax_dict["f"].legend(
        handles=h,
        labels=['CI-GWAS', 'CAUSE', 'MR-PRESSO', 'IVW Regression'],
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=2,
        # title="method",
    )
    ax_dict["f"].axis("off")
    plt.savefig('/nfs/scistore13/robingrp/human_data/causality/figures/figure_2.eps')

def plot_simulation_results_sup():
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

    alpha = 10**(-8)
    # methods = ['ci-gwas', 'cause', 'mrpresso', 'egger', 'ivw']
    methods = ['ci-gwas', 'cause', 'mrpresso', 'ivw']
    # methods = ['ci-gwas']

    def plot_bars(x_tup, metric, means, stds, ax):
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
        ax.set_xlabel("(n, m)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # ax.set_title('Penguin attributes by species')
        ax.set_xticks(x + width, x_tup)
        plt.setp(ax.get_xticklabels(), rotation=50, ha="right", rotation_mode="anchor")
        # ax.legend(loc='upper left', ncols=3)
        return handles

    x_tup = list(itertools.product(n_arr, m_arr))


    fig = plt.figure(layout="constrained", figsize=(15, 15))
    ax_dict = fig.subplot_mosaic(
        """
        ab
        cd
        ef
        """,
        empty_sentinel="X",
        sharex=True
    )

    plot_bars(x_tup, "fdr", means, stds, ax_dict['a'])
    plot_bars(x_tup, "tpr", means, stds, ax_dict['b'])
    plot_bars(x_tup, "mse", means, stds, ax_dict['c'])
    plot_bars(x_tup, "bias", means, stds, ax_dict['d'])
    h = plot_bars(x_tup, "var", means, stds, ax_dict['e'])
    ax_dict["f"].legend(
        handles=h,
        loc="center",
        fancybox=False,
        shadow=False,
        ncol=2,
        title="method",
    )
    ax_dict["f"].axis("off")
    fig.suptitle(rf"$\alpha={{{alpha}}}$")

def plot_block_size_experiment_results():
    block_sizes = np.arange(1, 12) * 10**3
    basepath = "/nfs/scistore13/robingrp/human_data/causality/block_size_effect/"
    bim_path = basepath + "ukb22828_UKB_EST_v3_ldp08.bim"
    bim = pd.read_csv(bim_path, sep="\t", header=None)

    bps = {}

    for bs in block_sizes:
        blockfile = basepath + f"ukb22828_UKB_EST_v3_ldp08_m{bs}_chr1.blocks"
        outdir = basepath + f"bdpc_d1_l6_a1e2_m{bs}_chr1/"
        gr = bdpc.merge_block_outputs(blockfile, outdir)
        bps[bs] = sorted(bim.loc[list(gr.gmi.values())][3].values)

    num_markers_selected = [len(bps[bs]) for bs in block_sizes]

    title_kw = {"loc": "left", "pad": 15, "size": 20}

    fig = plt.figure(layout="constrained", figsize=(10, 4))

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
        ax.plot(x, y, "|", color="k")
    ax.set_ylabel("max block size")
    ax.set_xlabel("position of selected marker [Mbp]")
    ax.set_title(title, **title_kw)

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
    es = [3, 4, 5, 6, 7, 8]
    nrows = len(es) - 1
    pheno_path = f"/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/input.phen"

    ace = {}

    for i, e in enumerate(es):
        bdpc_ace_path_sk = f"/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/bdpc_d{d}_l{l}_a1e{e}/all_merged_ACE_sk_mk3.mtx"

        ace[e] = bdpc.load_ace(bdpc_ace_path_sk, pheno_path).flatten()

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

    fig.suptitle(r"$ACE; k\leq3; SK$")
    plt.tight_layout()


def plot_compare_max_k_effect_on_ace():
    d = 1
    l = 6
    es = [5, 6, 7, 8]
    nrows = len(es)
    pheno_path = f"/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/input.phen"

    ace_mk = {}
    ace = {}

    for i, e in enumerate(es):
        bdpc_ace_path_sk_mk = f"/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/bdpc_d{d}_l{l}_a1e{e}/all_merged_ACE_sk_mk3.mtx"
        bdpc_ace_path_sk = f"/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/bdpc_d{d}_l{l}_a1e{e}/all_merged_ACE_sk.mtx"
        ace_mk[e] = bdpc.load_ace(bdpc_ace_path_sk_mk, pheno_path).flatten()
        ace[e] = bdpc.load_ace(bdpc_ace_path_sk, pheno_path).flatten()

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


def plot_ace_results_comp_cause_production(
    pag_path: str, ace_path: str, pheno_path: str, ace_norm=None
):

    # supp table 3 first three cols

    cause_gamma = pd.DataFrame(
        [
            {"y1": "bw".upper(), "y2": "t2d".upper(), "gamma": -0.274768438468858},
            {"y1": "bw".upper(), "y2": "ST".upper(), "gamma": -0.11133390989017},
            {"y1": "bw".upper(), "y2": "AT".upper(), "gamma": 0.0569668051462395},
            {"y1": "bw".upper(), "y2": "cad".upper(), "gamma": -0.142129232069397},
            {"y1": "bmi".upper(), "y2": "t2d".upper(), "gamma": 0.753001009486836},
            {"y1": "bmi".upper(), "y2": "ST".upper(), "gamma": 0.0724579259674975},
            {"y1": "bmi".upper(), "y2": "AT".upper(), "gamma": 0.127570936884617},
            {"y1": "bmi".upper(), "y2": "cad".upper(), "gamma": 0.254304518656026},
            {"y1": "HT".upper(), "y2": "t2d".upper(), "gamma": 0.0133993365533514},
            {"y1": "HT".upper(), "y2": "ST".upper(), "gamma": -0.0158921479966716},
            {"y1": "HT".upper(), "y2": "AT".upper(), "gamma": -0.00266884948627959},
            {"y1": "HT".upper(), "y2": "cad".upper(), "gamma": -0.0646821638302025},
            {"y1": "hdl".upper(), "y2": "t2d".upper(), "gamma": -0.157545489184377},
            {"y1": "hdl".upper(), "y2": "ST".upper(), "gamma": -0.0346983958224145},
            {"y1": "hdl".upper(), "y2": "AT".upper(), "gamma": 0.00954960164699658},
            {"y1": "hdl".upper(), "y2": "cad".upper(), "gamma": -0.203750511894064},
            {"y1": "ldl".upper(), "y2": "t2d".upper(), "gamma": -0.125971830503692},
            {"y1": "ldl".upper(), "y2": "ST".upper(), "gamma": 0.0618108230339071},
            {"y1": "ldl".upper(), "y2": "AT".upper(), "gamma": -0.0178475072825962},
            {"y1": "ldl".upper(), "y2": "cad".upper(), "gamma": 0.363195509520469},
            {"y1": "TRIG".upper(), "y2": "t2d".upper(), "gamma": 0.181573783879022},
            {"y1": "TRIG".upper(), "y2": "ST".upper(), "gamma": 0.0081104586416291},
            {"y1": "TRIG".upper(), "y2": "AT".upper(), "gamma": -0.0917176821192589},
            {"y1": "TRIG".upper(), "y2": "cad".upper(), "gamma": 0.281762788578154},
            {"y1": "ALC".upper(), "y2": "t2d".upper(), "gamma": 0.0692860734762807},
            {"y1": "ALC".upper(), "y2": "ST".upper(), "gamma": 0.125010784103481},
            {"y1": "ALC".upper(), "y2": "AT".upper(), "gamma": -0.070909459789924},
            {"y1": "ALC".upper(), "y2": "cad".upper(), "gamma": 0.0369981351799616},
            {"y1": "SMK".upper(), "y2": "t2d".upper(), "gamma": 0.153669125333255},
            {"y1": "SMK".upper(), "y2": "ST".upper(), "gamma": 0.287713602378084},
            {"y1": "SMK".upper(), "y2": "AT".upper(), "gamma": 0.129522157700426},
            {"y1": "SMK".upper(), "y2": "cad".upper(), "gamma": 0.479817266067469},
            {"y1": "bfp".upper(), "y2": "t2d".upper(), "gamma": 0.0647735280855574},
            {"y1": "bfp".upper(), "y2": "ST".upper(), "gamma": 0.0044428460841386},
            {"y1": "bfp".upper(), "y2": "AT".upper(), "gamma": 0.0957086524513934},
            {"y1": "bfp".upper(), "y2": "cad".upper(), "gamma": 0.133460312029649},
            {"y1": "fg".upper(), "y2": "t2d".upper(), "gamma": 1.32033326715937},
            {"y1": "fg".upper(), "y2": "ST".upper(), "gamma": 0.00957326536858055},
            {"y1": "fg".upper(), "y2": "AT".upper(), "gamma": -0.200412032138845},
            {"y1": "fg".upper(), "y2": "cad".upper(), "gamma": 0.111579686635068},
            {"y1": "dbp".upper(), "y2": "t2d".upper(), "gamma": 0.0204781118822844},
            {"y1": "dbp".upper(), "y2": "ST".upper(), "gamma": 0.0318219908469929},
            {"y1": "dbp".upper(), "y2": "AT".upper(), "gamma": 0.00387409082195085},
            {"y1": "dbp".upper(), "y2": "cad".upper(), "gamma": 0.0371147380330953},
            {
                "y1": "bp".upper(),
                "y2": "t2d".upper(),
                "gamma": 0.0142748064016252,
            },  # changed from sbp,
            {
                "y1": "bp".upper(),
                "y2": "ST".upper(),
                "gamma": 0.0216220295177384,
            },  # changed from sbp,
            {
                "y1": "bp".upper(),
                "y2": "AT".upper(),
                "gamma": 0.00246673566316583,
            },  # changed from sbp,
            {"y1": "bp".upper(), "y2": "cad".upper(), "gamma": 0.024807681904989},
        ]
    )  # changed from sbp)

    pnames = get_pheno_codes(pheno_path)
    cause_ys = set(cause_gamma["y1"].values)
    cause_ys.update(set(cause_gamma["y2"].values))
    reg_pnames = cause_ys.intersection(set(pnames))

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

    diseases = set(["AT", "ST", "T2D", "CAD"])

    risk_factors = set(
        [
            "BMI",
            "HT",
            "BP",
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

    reg_cfg = {
        "skip_na": False,
        "rf2rf": False,
        "d2d": False,
        "rf2d": True,
        "d2rf": False,
    }

    fig = plt.figure(layout="constrained", figsize=(17, 10))
    ax_dict = fig.subplot_mosaic(
        """
        abc
        deX
        """,
        empty_sentinel="X",
    )

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
    )

    plot_pleiotropy_mat(
        pag_path,
        pheno_path,
        ax=ax_dict["d"],
        title="d)",
        cbar_kw=cbar_kw,
        title_kw=title_kw,
        cmap="BuPu",
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

    plot_non_pleio_barplot(pag_path, pheno_path, ax=ax_dict["e"], title="e)", title_kw=title_kw)

    ax = ax_dict["c"]
    ace = load_ace(ace_path, pheno_path)
    # ace_flat = ace.flatten()
    # mr = ace_flat + rng.normal(0, 0.10, size=len(ace_flat))

    rows = []
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
            cg = cause_gamma[(cause_gamma["y1"] == pi) & (cause_gamma["y2"] == pj)].gamma.values
            if len(cg) > 1:
                raise ValueError("too many gammas")
            gamma = 0 if len(cg) == 0 else cg[0]
            if reg_cfg["skip_na"] and (gamma == 0 or ace[i, j] == 0):
                continue
            rows.append({"y1": pi, "y2": pj, "gamma": gamma, "ace": ace[i, j]})

    data = pd.DataFrame(rows)
    mr = data["gamma"].values
    ace_flat = data["ace"].values

    ax.scatter(mr, ace_flat, color="#d8dcd6", edgecolors="#5729ce", s=80, alpha=0.5, zorder=10)
    ax.grid(linestyle=":")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel("CI-GWAS ACE")
    ax.set_xlabel(r"CAUSE $\gamma$ (Morrison et al. 2020)")

    X = mr.reshape(-1, 1)
    y = ace_flat

    linear_regressor = LinearRegression()
    linear_regressor.fit(X, y)
    linear_regressor.fit(X, y)
    x_pred = np.linspace(np.min(mr), np.max(mr), 100).reshape(-1, 1)
    y_pred = linear_regressor.predict(x_pred)

    # Plot regression line.
    # * Logarithmic transformation is reverted by using the exponential one.
    ax.plot(x_pred, y_pred, color="#696969", lw=2, linestyle="dashed")
    beta = np.round(linear_regressor.coef_[0], 2)
    mu = np.round(linear_regressor.intercept_, 2)
    ax.set_title("c)", **title_kw)
    ax.text(0.2, 0.12, rf"$y={mu} + {beta}x$", fontsize=13)
    # fig.align_labels()
    # fig.tight_layout()


def plot_non_pleio_barplot(pag_path: str, pheno_path: str, ax=None, title=None, title_kw=None):
    if ax is None:
        ax = plt.gca()
    ps = pag_exclusive_pleiotropy_sets(pag_path, pheno_path, is_possible_child)
    p_names = get_pheno_codes(pheno_path)
    bd = [len(ps[i, i]) for i in range(len(p_names))]
    ax.bar(range(len(bd)), bd, color="purple")
    ax.set_ylabel("# non-pleiotropic parent markers")
    ax.set_xticks(np.arange(len(bd)), labels=p_names)
    plt.setp(ax.get_xticklabels(), rotation=50, ha="right", rotation_mode="anchor")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(linestyle=":", axis="y")
    ax.set_title(title, **title_kw)


def plot_all_alpha_production_pags_and_ace(basepath: str, pheno_path: str, l=6, d=1):
    fig = plt.figure(layout="constrained", figsize=(10, 35))
    ax_dict = fig.subplot_mosaic(
        """
        ah
        bi
        cj
        dk
        el
        fm
        gn
        op
        """,
        empty_sentinel="X",
        height_ratios=[1, 1, 1, 1, 1, 1, 1, 0.05],
        
    )

    ace_norm = mpl.colors.SymLogNorm(vmin=-2.0, vmax=2.0, linthresh=0.01)
    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    # title_kw = {"loc": "left", "pad": 15, "size": 20}
    title_kw = {}

    pag_panels = ["a", "b", "c", "d", "e", "f", "g"]
    ace_panels = ["h", "i", "j", "k", "l", "m", "n"]
    alphas = [2, 3, 4, 5, 6, 7, 8]

    for pag_panel, ace_panel, a in zip(pag_panels, ace_panels, alphas):
        try:
            ace_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_ACE_sk_mk3.mtx"
            pag_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_pag_sk_mk3.mtx"

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

    fig.colorbar(ima, cax=ax_dict['p'], **cbar_kw, location='bottom')

    for pag_panel, ace_panel, a in zip(pag_panels, ace_panels, alphas):
        try:
            ace_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_ACE_sk_mk3.mtx"
            pag_path = basepath + f"/bdpc_d{d}_l{l}_a1e{a}/all_merged_pag_sk_mk3.mtx"

            imp = plot_pag(
                pag_path,
                pheno_path,
                ax=ax_dict[pag_panel],
                title=rf"$\alpha=10^{{-{a}}}$",
                title_kw=title_kw,
                edge_encoding=five_common_edge_types,
                cbar_kw=cbar_kw,
                cbar=False
            )

        except FileNotFoundError as e:
            print(e)

    fig.colorbar(imp, cax=ax_dict['o'], location='bottom')
    # cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")


def plot_all_alpha_production_pags(basepath: str, pheno_path: str, l=6, d=1):
    fig = plt.figure(layout="constrained", figsize=(17, 10))
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
    panels = [
        "a",
        "b",
        "c",
        "d",
        "e",
        "f"
    ]

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
    fig = plt.figure(layout="constrained", figsize=(11, 10))
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
    fig = plt.figure(layout="constrained", figsize=(11, 10))
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

def merge_parallel_davs_csvs(indir: str, outdir: str, num_phen: int):
    res = np.zeros(shape=(num_phen, num_phen))
    missing = []
    for i in range(1, num_phen + 1):
        for j in range(1, num_phen + 1):
            if i == j:
                continue
            file = indir + f"all_merged_ACE_sk_mk3_i{i}_j{j}.csv"
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

    scipy.io.mmwrite(outdir + "all_merged_ACE_sk_mk3.mtx", scipy.sparse.coo_matrix(res))


def plot_inf_depth_pleio():
    d = 10000
    l = 6
    e = 2
    outdir = f"/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/bdpc_d{d}_l{l}_a1e{e}/"
    blockfile = "/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/ukb22828_UKB_EST_v3_ldp08.blocks"
    pheno_path = f"/nfs/scistore13/robingrp/human_data/causality/parent_set_selection/production/input.phen"

    z1 = get_skeleton_pleiotropy_mat(outdir, blockfile, pheno_path, max_depth=1)
    z2 = get_skeleton_pleiotropy_mat(outdir, blockfile, pheno_path, max_depth=2)
    zinf = get_skeleton_pleiotropy_mat(outdir, blockfile, pheno_path, mat_type="union")

    fig = plt.figure(layout="constrained", figsize=(15, 10))
    ax_dict = fig.subplot_mosaic(
        """
        abc
        dXX
        """,
        empty_sentinel="X"
    )

    norm = mpl.colors.LogNorm(vmin=1, vmax=np.max(zinf))
    cmap="PuBuGn"

    plot_skeleton_pleiotropy_mat_z(z1, pheno_path, ax=ax_dict['a'], title='depth=1', norm=norm, cmap=cmap, cbar=False, aspect='auto')
    ax_dict['a'].text(0, -2.5, "a)", size=20, verticalalignment='top')
    plot_skeleton_pleiotropy_mat_z(z2, pheno_path, ax=ax_dict['b'], title='depth=2', norm=norm, cmap=cmap, cbar=False, aspect='auto')
    ax_dict['b'].text(0, -2.5, "b)", size=20, verticalalignment='top')
    plot_skeleton_pleiotropy_mat_z(zinf, pheno_path, ax=ax_dict['c'], title='depth=inf', norm=norm, cmap=cmap, cbar=True, aspect='auto')
    ax_dict['c'].text(0, -2.5, "c)", size=20, verticalalignment='top')

    x = np.diag(z1)
    y = np.diag(zinf)

    ax = ax_dict['d']

    X = x.reshape(-1, 1)

    linear_regressor = LinearRegression()
    linear_regressor.fit(X, y)
    linear_regressor.fit(X, y)
    x_pred = np.linspace(np.min(x), np.max(x), 100).reshape(-1, 1)
    y_pred = linear_regressor.predict(x_pred)

    ax.plot(x_pred, y_pred, color="#696969", lw=2, linestyle="dashed")
    beta = np.round(linear_regressor.coef_[0], 2)
    mu = np.round(linear_regressor.intercept_, 2)
    # ax.set_title("c)", **title_kw)
    ax.text(100, 0.2 * 10 ** 6, rf"$y={mu} + {beta}x$", fontsize=13)
    ax.text(100, 0.1 * 10 ** 6, rf"$\rho={np.round(np.corrcoef(x, y)[0, 1], 4)}$", fontsize=13)
    ax.plot(x, y, 'o')
    ax.set_ylabel('# shared ancestral markers')
    ax.set_xlabel('# parent markers')
    title_kw = {"loc": "left", "pad": 15, "size": 20}
    ax_dict['d'].set_title("d)", **title_kw)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.savefig('/nfs/scistore13/robingrp/human_data/causality/figures/sup_figure_inf_depth_pleio.eps')


def plot_est_ukb_corr():
    basepath = '/nfs/scistore13/robingrp/human_data/causality/phenopc/phenopc_l14_a1e2/pheno_sk'
    num_m = 0
    num_p = 17
    marker_offset = 0

    fig = plt.figure(layout="constrained", figsize=(11, 5))
    title_kw = {"loc": "left", "pad": 15, "size": 20}

    ax_dict = fig.subplot_mosaic(
        """
        ab
        """,
        empty_sentinel="X",
    )

    estonia_pnames = ["AT","BMI","CAD","DBP","HT","SBP","SMK","ST","T2D"]
    estonia_corr_path = "/nfs/scistore13/robingrp/human_data/cigwas_estonia/trait-trait.txt"

    p_corrs_ukb = load_corr_sparse(basepath, num_m, num_p, marker_offset)
    ukb_pnames = ["WHR","BMI","SBP","DBP","HbA1c","GLU","CHOL","TRIG","HDL","LDL","AT","CAD","ST","T2D","HT","ALC","SMK"]

    ukb_mat = np.zeros(shape=(9, 9))

    for ((i, j), c) in p_corrs_ukb.items():
        try:    
            ei = estonia_pnames.index(ukb_pnames[i - 1])
            ej = estonia_pnames.index(ukb_pnames[j - 1])
            ukb_mat[ei, ej] = c
        except ValueError:
            continue

    est_mat = []
    with open(estonia_corr_path, 'r') as fin:
        next(fin)
        for line in fin:
            est_mat.append([eval(e) for e in line.split()[1:]])
    est_mat = np.array(est_mat)
    for i in range(9):
        for j in range(i + 1, 9):
            est_mat[j, i] = est_mat[i, j]

    x = np.triu(ukb_mat, k=1).flatten()
    y = np.triu(est_mat, k=1).flatten()

    ax = ax_dict['a']
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
    ax.plot(x, y, 'o')
    ax.set_ylabel('trait-trait corr EST')
    ax.set_xlabel('trait-trait corr UKB')
    title_kw = {"loc": "left", "pad": 15, "size": 20}
    # ax_dict['d'].set_title("d)", **title_kw)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    z = est_mat - ukb_mat
    mask = ~np.tri(z.shape[0], k=-1, dtype=bool)
    z = np.ma.array(z, mask=mask)  # mask out the lower triangle
    cmap = plt.get_cmap("coolwarm")
    cmap.set_bad("w")  # default value is 'k'
    cbar_kw = {"fraction": 0.046, "pad": 0.04}
    ax_dict['b']

    heatmap(
        z,
        estonia_pnames,
        estonia_pnames,
        ax_dict['b'],
        cbarlabel=r'$\Delta corr$',
        cmap=cmap,
        cbar_kw=cbar_kw,
        vmin=-0.12,
        vmax=0.12,
        title_kw=title_kw,
        title='b)'
        )

    plt.savefig('/nfs/scistore13/robingrp/human_data/causality/figures/sup_figure_est_ukb_corr.eps')