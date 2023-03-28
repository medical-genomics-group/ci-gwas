#pragma once

#include <mps/io.h>

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

std::unordered_set<int> parent_set(
    const std::vector<int> &G, const size_t num_var, const size_t num_markers, const int max_depth
);

std::vector<int> set_to_vec(const std::unordered_set<int>);

void direct_x_to_y(std::vector<int> &G, const size_t num_var, const size_t num_markers);

struct ReducedGCS
{
    size_t num_var;
    size_t num_phen;
    size_t max_level;
    std::vector<int> G;
    std::vector<float> C;
    std::vector<int> S;

    void to_file(std::string base)
    {
        std::ofstream fout;
        fout.open(base + ".mdim", std::ios::out);
        fout << num_var << "\t" << num_phen << "\t" << max_level << std::endl;
        fout.close();
        write_ints_to_binary(G.data(), G.size(), base + ".adj");
        write_floats_to_binary(C.data(), C.size(), base + ".corr");
        write_ints_to_binary(S.data(), S.size(), base + ".sep");
    }
};

/**
 * @brief Removes non-phenotype and non-ancestral nodes from adajcency matrix, corr matrix and
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

class ParentSetIndices
{
   private:
    std::vector<int> new_to_old;
    std::unordered_map<int, int> old_to_new;

   public:
    ParentSetIndices(const std::unordered_set<int> &p)
    {
        new_to_old = set_to_vec(p);
        for (int i = 0; i < new_to_old.size(); i++)
        {
            old_to_new[new_to_old[i]] = i;
        }
    }

    int get_old_ix(const int new_ix) const { return new_to_old[new_ix]; }

    int get_new_ix(const int old_ix) { return old_to_new[old_ix]; }
};
