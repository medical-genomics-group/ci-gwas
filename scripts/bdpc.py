import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from dataclasses import dataclass
import json
import pandas as pd
import seaborn as sns
import queue
from scipy.io import mmread

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
        if sm_ix >= num_p:
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

    return GlobalBdpcResult(sam, scm, ssm, gmi, marker_offset + bo.num_phen(),
                            bo.num_phen(), bo.max_level())


def global_epm(blockfile: str, outdir: str):
    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    bo = BlockOutput(basepaths[0])
    marker_offset = bo.num_markers()

    epm = bo.exclusive_pleiotropy_mat()

    for path in basepaths[1:]:
        try:
            bo = BlockOutput(path, marker_offset)
        except FileNotFoundError:
            continue
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


def global_ancestor_sets(blockfile: str,
                         outdir: str,
                         reduced_indices=False,
                         depth=1):
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
    """Compute upper bound of markers that could affect each phenotype or combination of phenotypes
    """
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
                    if not v2 in phens and not v2 in visited and neighbor_fn(
                            pag, v1, v2):
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


def pag_exclusive_pleiotropy_sets(pag_path: str,
                                  pheno_path: str,
                                  neighbor_fn,
                                  depth=1):
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

    def pheno_parents(self):
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
        """For each phenotype, find all markers that are direct parents to it.
        """
        res = {}
        adj = self.sam()
        for pix in self.pheno_indices():
            res[pix] = set()
            for parent_candidate in self.marker_indices():
                if (pix, parent_candidate) in adj:
                    res[pix].add(parent_candidate)
        return res

    def exclusive_pleiotropy_mat(self):
        pm = self.pheno_parents()
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
        pm = self.pheno_parents()
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


def heatmap(data,
            row_labels,
            col_labels,
            ax=None,
            cbar_kw=None,
            cbarlabel="",
            xlabel=None,
            ylabel=None,
            title=None,
            title_kw=dict(),
            **kwargs):
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

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(),
             rotation=45,
             ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    if title:
        ax.set_title(title, **title_kw)

    return im, cbar


def plot_pleiotropy_mat(pag_path: str,
                        pheno_path: str,
                        neighbor_fn=is_possible_child,
                        depth=1,
                        ax=None,
                        cbar_kw=None,
                        title=None,
                        title_kw=dict(),
                        cmap='BuPu'):
    poss_parents = pag_exclusive_pleiotropy_sets(pag_path, pheno_path,
                                                 neighbor_fn, depth)
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    z = [[0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(i + 1):
            z[i][j] = len(poss_parents[(i, j)])
    z = np.array(z)
    mask = ~np.tri(z.shape[0], k=0, dtype=bool)
    z = np.ma.array(z, mask=mask)  # mask out the lower triangle
    cmap = mpl.colormaps[cmap]  # jet doesn't have white color
    cmap.set_bad('w')  # default value is 'k'
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
        ax=None)


def load_ace(ace_path: str, pheno_path: str) -> np.array:
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    ace = mmread(ace_path).tocsr()
    z = [[0.0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(num_phen):
            z[i][j] = ace[i, j]
    return np.array(z)


def plot_ace(ace_path: str,
             pheno_path: str,
             title=None,
             title_kw=dict(),
             cmap="bwr",
             cbarlabel=r"$ACE \: (y_1 \rightarrow y_2)$",
             cbar_kw=None,
             ax=None):
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    ace = mmread(ace_path).tocsr()
    z = [[0.0 for _ in range(num_phen)] for _ in range(num_phen)]
    for i in range(num_phen):
        for j in range(num_phen):
            z[i][j] = ace[i, j]
    max_z = np.max(z)
    im, _ = heatmap(np.array(z),
                    p_names,
                    p_names,
                    cmap=cmap,
                    cbarlabel=cbarlabel,
                    cbar_kw=cbar_kw,
                    vmin=-max_z,
                    vmax=max_z,
                    xlabel=r"$y_2$",
                    ylabel=r"$y_1$",
                    title=title,
                    title_kw=title_kw,
                    ax=ax)


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
        np.array([
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
        ])))

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
        np.array([
            "#ffffff",  # white
            "#000421",
            "#46204c",
            "#99305b",
            "#dd584b",
            "#fb9d20",
            "#e5ed05"
        ])))

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
        np.array([
            "#ffffff",  # white
            "#000421",
            "#46204c",
            "#99305b",
            "#dd584b",
            "#fb9d20",
            "#e5ed05"
        ])))

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
        np.array([
            "#ffffff",  # white
            "#000421",
            "#5b2453",
            "#bf4056",
            "#f88a2e",
            "#e5ed05",
        ])))


def plot_pag(pag_path: str,
             pheno_path: str,
             title=None,
             title_kw=dict(),
             cbar_kw=None,
             edge_encoding=all_edge_types,
             ax=None):

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
    fmt = mpl.ticker.FuncFormatter(
        lambda x, pos: edge_encoding.str_rep[norm(x)])

    if cbar_kw is None:
        cbar_kw = {}
    cbar_kw["ticks"] = np.arange(ne) + 0.5
    cbar_kw["format"] = fmt

    im, _ = heatmap(np.array(z),
                    p_names,
                    p_names,
                    cmap=edge_encoding.cmap,
                    norm=norm,
                    cbar_kw=cbar_kw,
                    cbarlabel="Edge Type",
                    xlabel=r"$y_2$",
                    ylabel=r"$y_1$",
                    title=title,
                    title_kw=title_kw,
                    ax=ax)


def marker_pheno_associations(blockfile: str,
                              outdir: str,
                              pag_path: str,
                              pheno_path: str,
                              bim_path: str,
                              neighbor_fn=is_possible_child,
                              depth=1):
    gr = merge_block_outputs(blockfile, outdir)
    rs_ids = pd.read_csv(bim_path, sep='\t', header=None)[1].values
    p_names = get_pheno_codes(pheno_path)
    num_phen = len(p_names)
    pag = mmread(pag_path).tocsr()
    geps = pag_exclusive_pleiotropy_sets(pag, pheno_path, neighbor_fn, depth)
    assoc_markers = []

    for pix in range(num_phen):
        for v in geps[(pix, pix)]:
            # pag is 0-based indexed, everything else uses 1-based indexing
            bim_line = int(gr.gmi[v + BASE_INDEX])
            assoc_markers.append({
                "phenotype": p_names[pix],
                "rsID": rs_ids[bim_line],
                "bim_line_ix": bim_line
            })

    return pd.DataFrame(assoc_markers)


def pag_edge_types(pag_path: str,
                   pheno_path: str) -> dict[tuple[int, int], int]:
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


def pag_x_to_y_edge_types(pag_path: str,
                          pheno_path: str) -> dict[tuple[int, int], int]:
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
    sns.heatmap(g,
                xticklabels=1,
                yticklabels=1,
                norm=matplotlib.colors.LogNorm(),
                cmap="Blues",
                annot=True,
                square=True,
                mask=mask,
                cbar_kws={'label': '# parent markers'})
