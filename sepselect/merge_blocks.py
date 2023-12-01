from dataclasses import dataclass
import queue
import numpy as np

BASE_INDEX = 1


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


def load_corr_sparse(basepath: str, num_m: int, num_p: int, marker_offset):
    return load_mat_sparse(basepath, num_m, num_p, marker_offset, np.float32, ".corr")


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


def load_adj_sparse(basepath: str, num_m: int, num_p: int, marker_offset=0):
    return load_mat_sparse(basepath, num_m, num_p, marker_offset, np.int32, ".adj")


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


class BlockOutput:
    def __init__(self, basepath: str, marker_offset=0, global_marker_offset=0):
        self.basepath = basepath
        self.mdim = load_mdim(basepath)
        # number of selected markers in all previous blocks
        self.marker_offset = marker_offset
        # .bim row index of first marker in block definition
        self.global_marker_offset = global_marker_offset

    def __str__(self) -> str:
        chrom, first, last = self.basepath.split("/")[-1].split("_")
        return f"{chrom}:{first}-{last}"

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

        with open(basepath + ".mdim", "w") as fout:
            fout.write(f"{self.num_var}\t{self.num_phen}\t{self.max_level}\n")

        np.array(sorted(list(self.gmi.values())), dtype=np.int32).tofile(
            basepath + ".ixs"
        )


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


def block_size(basepath: str) -> int:
    first, last = basepath.split("_")[-2:]
    return int(last) - int(first) + 1


def merge_block_outputs(blockfile: str, outdir: str):
    basepaths = [outdir + s for s in get_block_out_stems(blockfile)]

    try:
        bo = BlockOutput(basepaths[0])
        marker_offset = bo.num_markers()
        global_marker_offset = bo.block_size()
        sam = bo.sam()
        scm = bo.scm()
        gmi = bo.gmi()
    except FileNotFoundError:
        path = basepaths[0]
        print(f"Missing: {path}")
        global_marker_offset = block_size(path)
        marker_offset = 0
        sam = {}
        scm = {}
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
        add_gmi(gmi, bo.gmi())
        marker_offset += bo.num_markers()
        global_marker_offset += bo.block_size()

    return GlobalBdpcResult(
        sam, scm, gmi, marker_offset + bo.num_phen(), bo.num_phen(), bo.max_level()
    )
