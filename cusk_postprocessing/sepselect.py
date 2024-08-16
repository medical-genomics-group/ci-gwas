import numpy as np
import scipy
from scipy.stats import norm
from scipy.io import mmread, mmwrite
import sys


def pcorr(corr: np.array, var: list[int]):
    try:
        # return np.linalg.inv(corr[np.ix_(var, var)])
        return np.linalg.inv(corr[np.ix_(var, var)])
    except np.linalg.LinAlgError as err:
        print(f"vars: {var}")
        print(corr[np.ix_(var, var)])
        m = np.linalg.pinv(corr[np.ix_(var, var)])
        print(f"pcorr pinv: {fisher_z(-(m[0, 1] / np.sqrt(m[0, 0] * m[1, 1])))}")
        print(err)
        sys.exit()


def fisher_z(v):
    return np.abs(0.5 * np.log(np.abs((1 + v) / (1 - v))))


def alpha_thr(alpha: int, n: int, l: int):
    return norm.ppf(1 - (alpha / 2)) / np.sqrt(n - l - 3)


def independent(partial_correlation, alpha, n, l):
    return partial_correlation < alpha_thr(alpha, n, l)


def dag_to_pag_notation(dag: np.array):
    pag = np.zeros_like(dag)
    pag[dag != 0] = 2
    pag[dag.T != 0] = 3
    return pag


class CuskResults:
    def __init__(self, stem: str):
        with open(f"{stem}.mdim", "r") as fin:
            self.num_var, self.num_phen, self.max_level = [
                int(e) for e in next(fin).split()
            ]

        self.num_m = self.num_var - self.num_phen
        self.sepset_arr = np.fromfile(f"{stem}.sep", dtype=np.int32).reshape(
            self.num_var, self.num_var, self.max_level
        )
        self.ixs = np.fromfile(f"{stem}.ixs", dtype=np.int32)
        self.marker_ixs = set(self.ixs[: self.num_m])
        self.adj = np.fromfile(f"{stem}.adj", dtype=np.int32).reshape(
            self.num_var, self.num_var
        )
        self.adj = self.adj.astype(bool)
        self.corr = np.fromfile(f"{stem}.corr", dtype=np.float32).reshape(
            self.num_var, self.num_var
        )
        self.selected_markers = [e for e in self.ixs if e < self.num_m]
        self.max_level_maximal_sepsets = None
        self.maximal_sepset_arr = None
        self.max_level_minimal_sepsets = None
        self.minimal_sepset_arr = None
        self.max_level_minimal_pcorr_sepsets = None
        self.minimal_pcorr_sepset_arr = None
        self.unshielded_triples = None
        self.ambiguous_triples = None
        self.rfci_relevant_unshielded_triples = None

    def get_rfci_relevant_unshielded_triples(self):
        """Only unshielded triples with two or more traits,
        where the middle variable is always a trait

        Without duplicates, and always coded x, y, z, s
        s.t. x < z
        """
        if self.rfci_relevant_unshielded_triples is None:
            self.rfci_relevant_unshielded_triples = []
            for (x, y, z) in self.get_unshielded_triples():
                if not self.is_marker(y) and x < z and sum(self.is_marker(e) for e in [x, y, z]) < 2:
                    self.rfci_relevant_unshielded_triples.append([x, y, z])
            self.rfci_relevant_unshielded_triples = np.array(self.rfci_relevant_unshielded_triples, dtype=np.int32)
        return self.rfci_relevant_unshielded_triples

    def get_rfci_relevant_unshielded_triples_outer_pairs(self):
        res = set()
        for t in self.get_rfci_relevant_unshielded_triples():
            res.add((t[0], t[2]))
            res.add((t[2], t[0]))
        return res

    def is_marker(self, variable_ix: int) -> bool:
        return variable_ix < self.num_m

    def mark_ambiguous_triples(self):
        if self.maximal_sepset_arr is None or self.minimal_pcorr_sepset_arr is None:
            raise RuntimeError(
                "Cannot mark ambiguous triples without minimal pcorr and maximal sepsets"
            )
        self.ambiguous_triples = []
        for a, b, c in self.get_unshielded_triples():
            if np.any(self.maximal_sepset_arr[a, c] == b) and np.all(
                self.minimal_pcorr_sepset_arr[a, c] != b
            ):
                self.ambiguous_triples.append([a, b, c])
        self.ambiguous_triples = np.array(self.ambiguous_triples, dtype=np.int32)

    def to_file(self, stem: str):
        with open(stem + ".mdim", "w") as fout:
            fout.write(
                f"{self.num_var}\t{self.num_phen}\t{self.max_level_maximal_sepsets}\t{self.ambiguous_triples.shape[0]}\t{self.get_rfci_relevant_unshielded_triples.shape[0]}\n"
            )
        self.ixs.tofile(f"{stem}.ixs")
        self.adj.astype(np.int32).tofile(f"{stem}.adj")
        self.corr.tofile(f"{stem}.corr")
        self.maximal_sepset_arr.flatten().tofile(f"{stem}.sep")
        self.ambiguous_triples.tofile(f"{stem}.atr")
        self.get_rfci_relevant_unshielded_triples().tofile(f"{stem}.ut")

    def sepset_sizes(self) -> list[int]:
        sepset_sizes = []
        for i in range(self.num_var):
            for j in range(self.num_var):
                sepset_sizes.append(np.sum(self.sepset_arr[i, j] != -1))
        return sepset_sizes

    def trait_neighbors(self, node_ix: int) -> np.array:
        neighbors = self.neighbors(node_ix)
        return neighbors[neighbors >= self.num_m]

    def neighbors(self, parent_ix: int):
        return np.where(self.adj[parent_ix, :])[0]

    def non_neighbors(self, node_ix: int):
        return np.where(~self.adj[node_ix, :])[0]

    def adjacent(self, a: int, b: int):
        return self.adj[a, b] or self.adj[b, a]

    def get_unshielded_triples(self):
        if self.unshielded_triples is None:
            self.unshielded_triples = set()
            for a in range(self.num_var):
                for b in self.neighbors(a):
                    # a into sepset(b, c) (common cause)
                    for c in self.neighbors(a):
                        if b != c and not self.adjacent(b, c):
                            self.unshielded_triples.add((b, a, c))
                    # b into sepset(a, c) (node on causal path)
                    for c in self.neighbors(b):
                        if c != a and not self.adjacent(a, c):
                            self.unshielded_triples.add((a, b, c))
        return self.unshielded_triples

    def get_unshielded_triples_outer_pairs(self):
        res = set()
        for t in self.get_unshielded_triples():
            res.add((t[0], t[2]))
        return res

    def partial_correlation(self, vals):
        m = pcorr(self.corr, vals)
        return fisher_z(-(m[0, 1] / np.sqrt(np.abs(m[0, 0] * m[1, 1]))))

    def find_min_pcorr_sepsets_incr(self, alpha: float, num_samples: int):
        """
        Find the sepsets that reduces the partial correlation the most,
        by growing the sepset one variabel at a time.
        """
        sepsets = {}

        for i, j in remaining_pairs:
            neighbors_i = list(self.neighbors(i))
            ref_pc = self.partial_correlation([i, j])
            candidate_sepset = []
            found_sepset = independent(ref_pc, alpha, num_samples, 0)
            remaining_neighbors = neighbors_i.copy()

            for sepset_size in range(len(neighbors_i)):
                loc_ref_pc = 1.0
                add_neighbor = None
                for neighbor in remaining_neighbors:
                    sepset = candidate_sepset.copy()
                    sepset.append(neighbor)
                    loc_pc = self.partial_correlation([i, j] + sepset)
                    if found_sepset and loc_pc <= ref_pc:
                        ref_pc = loc_pc
                        add_neighbor = neighbor
                    elif not found_sepset and loc_pc <= loc_ref_pc:
                        loc_ref_pc = loc_pc
                        add_neighbor = neighbor

                if add_neighbor is None:
                    break
                candidate_sepset.append(add_neighbor)
                remaining_neighbors.remove(add_neighbor)
                if not found_sepset:
                    found_sepset = independent(
                        loc_ref_pc, alpha, num_samples, sepset_size
                    )

            sepsets[(i, j)] = candidate_sepset

        self.max_level_minimal_pcorr_sepsets = max(len(v) for v in sepsets.values())
        self.minimal_pcorr_sepset_arr = np.full(
            (self.num_var, self.num_var, self.max_level_minimal_pcorr_sepsets),
            -1,
            dtype=np.int32,
        )
        for (i, j), v in sepsets.items():
            for k, e in enumerate(v):
                self.minimal_pcorr_sepset_arr[i, j, k] = e

    def find_minimal_sepsets(self, alpha: float, num_samples: int):
        """
        Find the smallest sepset that separates.
        """
        sepsets = {}

        for i, j in self.get_unshielded_triples_outer_pairs():
            neighbors_i = list(self.neighbors(i))
            ref_pc = self.partial_correlation([i, j])
            candidate_sepset = []
            found_sepset = independent(ref_pc, alpha, num_samples, 0)
            remaining_neighbors = neighbors_i.copy()

            for sepset_size in range(len(neighbors_i)):
                if found_sepset:
                    break
                loc_ref_pc = 1.0
                add_neighbor = None
                for neighbor in remaining_neighbors:
                    sepset = candidate_sepset.copy()
                    sepset.append(neighbor)
                    loc_pc = self.partial_correlation([i, j] + sepset)
                    if loc_pc <= loc_ref_pc:
                        loc_ref_pc = loc_pc
                        add_neighbor = neighbor

                if add_neighbor is None:
                    raise RuntimeError("Failed to add any variable to sepset")
                else:
                    candidate_sepset.append(add_neighbor)
                    remaining_neighbors.remove(add_neighbor)
                    found_sepset = independent(
                        loc_ref_pc, alpha, num_samples, sepset_size
                    )

            sepsets[(i, j)] = candidate_sepset

        self.max_level_minimal_sepsets = max(len(v) for v in sepsets.values())
        self.minimal_sepset_arr = np.full(
            (self.num_var, self.num_var, self.max_level_minimal_sepsets),
            -1,
            dtype=np.int32,
        )
        for (i, j), v in sepsets.items():
            for k, e in enumerate(v):
                self.minimal_sepset_arr[i, j, k] = e

    def find_maximal_and_min_pcorr_sepsets_incr(self, alpha: float, num_samples: int):
        max_sepsets = {}
        min_pcorr_sepsets = {}
        pairs = self.get_rfci_relevant_unshielded_triples_outer_pairs()
        for (count, (i, j)) in enumerate(pairs):
            if count % 1000 == 0:
                print(f"Processing outer pair #{count + 1} out of {len(remaining_pairs)}", flush=True)
            remaining_neighbors = set(self.trait_neighbors(i))
            max_sepset_size = len(remaining_neighbors)
            candidate_sepset = []
            found_sepset = independent(
                self.partial_correlation([i, j]), alpha, num_samples, 0
            )
            found_minimum = False
            last_ref_pc = np.inf

            for sepset_size in range(1, max_sepset_size + 1):
                add_neighbor = None
                ref_pc = np.inf
                sepset = candidate_sepset.copy()
                sepset.append(None)
                for neighbor in remaining_neighbors:
                    sepset[-1] = neighbor
                    loc_pc = self.partial_correlation([i, j] + sepset)
                    if loc_pc <= ref_pc:
                        ref_pc = loc_pc
                        add_neighbor = neighbor

                if ref_pc > last_ref_pc and found_sepset and not found_minimum:
                    found_minimum = True
                    min_pcorr_sepsets[(i, j)] = candidate_sepset

                loc_idp = independent(ref_pc, alpha, num_samples, sepset_size)

                if not loc_idp and found_sepset:
                    break
                if loc_idp:
                    found_sepset = True

                last_ref_pc = ref_pc
                candidate_sepset.append(add_neighbor)
                remaining_neighbors.remove(add_neighbor)

            max_sepsets[(i, j)] = candidate_sepset

        self.max_sepsets = max_sepsets
        self.max_level_maximal_sepsets = max(len(v) for v in max_sepsets.values())
        self.maximal_sepset_arr = np.full(
            (self.num_var, self.num_var, self.max_level_maximal_sepsets),
            -1,
            dtype=np.int32,
        )
        for (i, j), v in max_sepsets.items():
            for k, e in enumerate(v):
                self.maximal_sepset_arr[i, j, k] = e

        self.min_sepsets = min_sepsets
        self.max_level_minimal_pcorr_sepsets = max(
            len(v) for v in min_pcorr_sepsets.values()
        )
        self.minimal_pcorr_sepset_arr = np.full(
            (self.num_var, self.num_var, self.max_level_minimal_pcorr_sepsets),
            -1,
            dtype=np.int32,
        )
        for (i, j), v in min_pcorr_sepsets.items():
            for k, e in enumerate(v):
                self.minimal_pcorr_sepset_arr[i, j, k] = e

    def find_maximal_sepsets_incr(self, alpha: float, num_samples: int):
        """
        Attempt to find, for each pair, the largest sepset that separates,
        starting from the empty set, adding neighbors one by one.
        Terminate if a sepset was found and no new addition leads to separation.
        """
        sepsets = {}

        for i, j in self.get_unshielded_triples_outer_pairs():
            neighbors_i = list(self.neighbors(i))
            candidate_sepset = []
            found_sepset = independent(
                self.partial_correlation([i, j]), alpha, num_samples, 0
            )
            remaining_neighbors = neighbors_i.copy()

            for sepset_size in range(len(neighbors_i)):
                add_neighbor = None
                ref_pc = np.inf
                for neighbor in remaining_neighbors:
                    sepset = candidate_sepset.copy()
                    sepset.append(neighbor)
                    loc_pc = self.partial_correlation([i, j] + sepset)
                    if loc_pc <= ref_pc:
                        ref_pc = loc_pc
                        add_neighbor = neighbor

                loc_idp = independent(ref_pc, alpha, num_samples, sepset_size)

                if not loc_idp and found_sepset:
                    break
                elif loc_idp:
                    found_sepset = True

                candidate_sepset.append(add_neighbor)
                remaining_neighbors.remove(add_neighbor)

            sepsets[(i, j)] = candidate_sepset

        self.max_level_maximal_sepsets = max(len(v) for v in sepsets.values())
        self.maximal_sepset_arr = np.full(
            (self.num_var, self.num_var, self.max_level_maximal_sepsets),
            -1,
            dtype=np.int32,
        )
        for (i, j), v in sepsets.items():
            for k, e in enumerate(v):
                self.maximal_sepset_arr[i, j, k] = e

    def find_maximal_sepsets_decr(self, alpha: float, num_samples: int):
        """
        Attempt to find, for each pair, the largest sepsets that separates
        """
        sepsets = {}

        for i, j in self.get_unshielded_triples_outer_pairs():
            neighbors_i = list(self.neighbors(i))

            full_con_pc = self.partial_correlation([i, j] + list(neighbors_i))
            ref_pc = full_con_pc
            ref_pc_sepset = neighbors_i
            candidate_sepset = neighbors_i

            while len(candidate_sepset) > 0 and not independent(
                ref_pc, alpha, num_samples, len(candidate_sepset)
            ):
                remove_ix = None
                ref_pc = 1.0
                for neighbor_ix in range(len(candidate_sepset)):
                    sepset = candidate_sepset.copy()
                    sepset.pop(neighbor_ix)
                    loc_pc = self.partial_correlation([i, j] + sepset)
                    if loc_pc <= ref_pc:
                        ref_pc = loc_pc
                        remove_ix = neighbor_ix
                if remove_ix is not None:
                    candidate_sepset.pop(remove_ix)
                    ref_pc_sepset = candidate_sepset.copy()
                else:
                    print("i, j:", i, j)
                    print("candidate sepset:", candidate_sepset)
                    raise RuntimeError("No element removed in sepset reduction round")

            sepsets[(i, j)] = ref_pc_sepset

        self.max_level_maximal_sepsets = max(len(v) for v in sepsets.values())
        self.maximal_sepset_arr = np.full(
            (self.num_var, self.num_var, self.max_level_maximal_sepsets),
            -1,
            dtype=np.int32,
        )
        for (i, j), v in sepsets.items():
            for k, e in enumerate(v):
                self.maximal_sepset_arr[i, j, k] = e


class MergedCuskResults(CuskResults):
    def __init__(self, stem, orientation_prior_file=None):
        with open(f"{stem}.mdim", "r") as fin:
            self.num_var, self.num_phen, self.max_level = [
                int(e) for e in next(fin).split()
            ]
        self.num_m = self.num_var - self.num_phen
        self.sepset_arr = None
        self.ixs = np.fromfile(f"{stem}.ixs", dtype=np.int32)
        self.adj = mmread(f"{stem}_sam.mtx").toarray()
        self.adj = self.adj.astype(bool)
        self.corr = mmread(f"{stem}_scm.mtx").toarray()
        # make sure that diagonal is 1
        np.fill_diagonal(self.corr, 1.0)
        self.selected_markers = None
        self.max_sepsets = None
        self.max_level_maximal_sepsets = None
        self.maximal_sepset_arr = None
        self.min_sepsets = None
        self.max_level_minimal_sepsets = None
        self.minimal_sepset_arr = None
        self.max_level_minimal_pcorr_sepsets = None
        self.minimal_pcorr_sepset_arr = None
        self.unshielded_triples = None
        self.ambiguous_triples = None
        self.rfci_relevant_unshielded_triples = None
        self.rm_collinear_markers()
        self.orientation_prior = np.zeros_like(self.adj)
        # marker -> trait
        self.orientation_prior[self.num_phen:, :self.num_phen] = self.adj[self.num_phen:, :self.num_phen]
        self.pag = None
        if orientation_prior_file is not None:
            orientation_prior = np.fromfile(orientation_prior_file, dtype=np.int32)
            assert orientation_prior.shape[0] == self.self.num_phen ** 2, "orientation prior has to have n_trait * n_trait entries"
            self.orientation_prior[:self.num_phen, :self.num_phen] = self.orientation_prior.reshape(self.self.num_phen, self.self.num_phen)

    def rm_collinear_markers(self):
        n_rm = 0
        curr_i = self.num_phen
        while curr_i < self.num_var:
            if np.sum(self.corr[curr_i, :] == 1) > 1:
                self.corr = np.delete(self.corr, curr_i, 0)
                self.corr = np.delete(self.corr, curr_i, 1)
                self.adj = np.delete(self.adj, curr_i, 0)
                self.adj = np.delete(self.adj, curr_i, 1)
                self.ixs = np.delete(self.ixs, curr_i - self.num_phen)
                self.num_var -= 1
                n_rm += 1
            else:
                curr_i += 1
        print(f"Removed {n_rm} collinear markers")

    def orient_v_structures(self, alpha: float, num_samples: int):
        self.pag = np.zeros_like(self.adj)
        self.pag[self.adj] = 1
        if self.max_sepsets is None:
            self.find_maximal_and_min_pcorr_sepsets_incr(alpha, num_samples)
        for (x, y, z) in self.get_rfci_relevant_unshielded_triples():
            orient = y not in self.max_sepsets[(x, z)] and y not in self.max_sepsets[(z, x)]
            # x-y edge:
            if self.orientation_prior[x, y] == 1:
                self.pag[x, y] = 2
                self.pag[y, x] = 3
            elif self.orientation_prior[y, x] == 1:
                self.pag[y, x] = 2
                self.pag[x, y] = 3
            elif orient:
                self.pag[x, y] = 2
            # z-y edge:
            if self.orientation_prior[z, y] == 1:
                self.pag[z, y] = 2
                self.pag[y, z] = 3
            elif self.orientation_prior[y, z] == 1:
                self.pag[y, z] = 2
                self.pag[z, y] = 3
            elif orient:
                self.pag[z, y] = 2

    def trait_neighbors(self, node_ix: int) -> np.array:
        neighbors = self.neighbors(node_ix)
        return neighbors[neighbors < self.num_phen]

    def is_marker(self, variable_ix: int) -> bool:
        return variable_ix > self.num_phen

    def to_file(self, stem: str):
        num_ambiguous_triples = self.ambiguous_triples.shape[0]
        num_triples = self.get_rfci_relevant_unshielded_triples().shape[0]

        with open(stem + ".mdim", "w") as fout:
            fout.write(
                f"{self.num_var}\t{self.num_phen}\t{self.max_level_maximal_sepsets}\t{num_ambiguous_triples}\t{num_triples}\n"
            )

        mmwrite(f"{stem}_sam.mtx", scipy.sparse.coo_matrix(self.adj.astype(np.int32)))
        mmwrite(f"{stem}_scm.mtx", scipy.sparse.coo_matrix(self.corr))
        mmwrite(f"{stem}_spm.mtx", scipy.sparse.coo_matrix(self.pag))

        self.ambiguous_triples.tofile(f"{stem}.atr")
        self.get_rfci_relevant_unshielded_triples().tofile(f"{stem}.ut")
        self.max_sepset_to_file(stem)

    def max_sepset_to_file(self, stem):
        with open(f"{stem}.ssm", "w") as fout:
            for i in range(self.num_var):
                for j in range(self.num_var):
                    ij_sepset = self.maximal_sepset_arr[i, j]
                    row = np.array([i, j] + list(ij_sepset[ij_sepset != -1]))
                    if len(row) == 2:
                        continue
                    row += 1
                    fout.write(" ".join([str(e) for e in row]) + "\n")


def orient_v_structures_merged(
    cusk1_result_stem: str, alpha: float, num_samples: int, orientation_prior_file=None
) -> MergedCuskResults:
    cr = MergedCuskResults(cusk1_result_stem, orientation_prior_file=orientation_prior_file)
    print("Orienting v-structures")
    cr.orient_v_structures(alpha=alpha, num_samples=num_samples)
    cr.mark_ambiguous_triples()
    return cr


def sepselect_merged(
    cusk1_result_stem: str, alpha: float, num_samples: int
) -> MergedCuskResults:
    cr = MergedCuskResults(cusk1_result_stem)
    print("Starting sepselect")
    cr.find_maximal_and_min_pcorr_sepsets_incr(alpha, num_samples)
    cr.mark_ambiguous_triples()
    return cr


def sepselect(cusk1_result_stem: str, alpha: float, num_samples: int) -> CuskResults:
    cr = CuskResults(cusk1_result_stem)
    cr.find_maximal_sepsets_incr(alpha, num_samples)
    cr.find_min_pcorr_sepsets_incr(alpha, num_samples)
    cr.mark_ambiguous_triples()
    return cr
