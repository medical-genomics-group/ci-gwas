#pragma once

#include <mps/io.h>

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * @brief Make subset of indices of phenotypes and markers that are ancestral at depth <=
 * max_depth.
 *
 * @param G
 * @param num_var
 * @param num_markers
 * @param max_depth
 * @return std::unordered_set<int>
 */
std::unordered_set<int> subset_variables(
    const std::vector<int> &G, const size_t num_var, const size_t num_markers, const int max_depth
);

std::vector<int> set_to_vec(const std::unordered_set<int>);

void direct_x_to_y(std::vector<int> &G, const size_t num_var, const size_t num_markers);

struct ReducedGCS
{
    size_t num_var;
    size_t num_phen;
    size_t max_level;
    std::vector<int> new_to_old_indices;
    std::vector<int> G;
    std::vector<float> C;
    std::vector<int> S;

    size_t num_markers() { return num_var - num_phen; }

    void to_file(std::string base)
    {
        std::ofstream fout;
        fout.open(base + ".mdim", std::ios::out);
        fout << num_var << "\t" << num_phen << "\t" << max_level << std::endl;
        fout.close();
        write_ints_to_binary(new_to_old_indices.data(), new_to_old_indices.size(), base + ".ixs");
        write_ints_to_binary(G.data(), G.size(), base + ".adj");
        write_floats_to_binary(C.data(), C.size(), base + ".corr");
        write_ints_to_binary(S.data(), S.size(), base + ".sep");
    }
};

/**
 * @brief Removes everything but ancestral and phenotype nodes from adajcency matrix, corr matrix and
 * separations sets.
 *
 * @param G         n*n adjacency matrix
 * @param C         n*n correlation matrix
 * @param S         n*n*ML separation set matrix
 * @param P         set of indices of nodes to be retained
 * @param num_var   number of variables in G, C, S
 * @param max_level max size of a separation set
 */
ReducedGCS reduce_gcs(
    const std::vector<int> &G,
    const std::vector<float> &C,
    const std::vector<int> &S,
    const std::unordered_set<int> &P,
    const size_t num_var,
    const size_t num_phen,
    const size_t max_level
);

ReducedGCS reduce_gcs(
    const std::vector<int> &G,
    const std::vector<float> &C,
    const std::vector<int> &S,
    const std::unordered_set<int> &P,
    const size_t num_var,
    const size_t num_phen,
    const size_t max_level,
    const std::vector<int> &index_map
);

class VariableSubsetIndices
{
   private:
    std::vector<int> new_to_old;
    std::unordered_map<int, int> old_to_new;

   public:
    VariableSubsetIndices(const std::unordered_set<int> &p)
    {
        new_to_old = set_to_vec(p);
        for (int i = 0; i < new_to_old.size(); i++)
        {
            old_to_new[new_to_old[i]] = i;
        }
    }

    VariableSubsetIndices(const std::unordered_set<int> &p, const std::vector<int> &index_map)
    {
        std::vector<int> new_to_old_tmp = set_to_vec(p);
        new_to_old = std::vector<int>(new_to_old_tmp.size());
        for (int i = 0; i < new_to_old.size(); i++)
        {
            new_to_old[i] = index_map[new_to_old_tmp[i]];
            old_to_new[new_to_old[i]] = i;
        }
    }

    int get_old_ix(const int new_ix) const { return new_to_old[new_ix]; }

    int get_new_ix(const int old_ix) { return old_to_new[old_ix]; }

    std::vector<int> new_to_old_indices() { return new_to_old; }
};
