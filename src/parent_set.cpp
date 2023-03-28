#include <mps/parent_set.h>

#include <algorithm>
#include <queue>
#include <unordered_set>
#include <vector>

std::unordered_set<int> parent_set(
    const std::vector<int> &G, const size_t num_var, const size_t num_markers, const int max_depth
)
{
    std::unordered_set<int> global_parent_set;

    // do bfs, starting from each of the markers.
    for (int start_ix = num_markers; start_ix < num_var; start_ix++)
    {
        std::unordered_set<int> curr_parent_set;
        std::queue<int> q;
        std::queue<int> next_q;
        for (int i = num_markers; i < num_var; ++i)
        {
            curr_parent_set.insert(i);
        }

        q.push(start_ix);

        for (int depth = 0; depth < max_depth; depth++)
        {
            next_q = {};
            while (!q.empty())
            {
                int node_ix = q.front();
                q.pop();
                for (int cix = 0; cix < num_markers; cix++)
                {
                    if ((G[node_ix * num_var + cix] == 1) && (!curr_parent_set.contains(cix)))
                    {
                        curr_parent_set.insert(cix);
                        next_q.push(cix);
                    }
                }
            }
            q = next_q;
        }

        for (auto e : curr_parent_set)
        {
            global_parent_set.insert(e);
        }
    }

    return global_parent_set;
}

std::vector<int> set_to_vec(const std::unordered_set<int> s)
{
    std::vector<int> v(s.begin(), s.end());
    std::sort(v.begin(), v.end());
    return v;
}

void direct_x_to_y(std::vector<int> &G, const size_t num_var, const size_t num_markers)
{
    std::unordered_set<int> phen_ids;
    for (int i = num_markers; i < num_var; ++i)
    {
        phen_ids.insert(i);
    }

    for (int sink = num_markers; sink < num_var; ++sink)
    {
        for (int source = 0; source < num_markers; source++)
        {
            if ((G[sink * num_var + source] == 1) && (!phen_ids.contains(source)))
            {
                // direct: source -> sink
                G[source * num_var + sink] = 2;
                G[sink * num_var + source] = 3;
            }
        }
    }
}

ReducedGCS reduce_gcs(
    const std::vector<int> &G,
    const std::vector<float> &C,
    const std::vector<int> &S,
    const std::unordered_set<int> &P,
    const size_t num_var,
    const size_t num_phen,
    const size_t max_level
)
{
    ReducedGCS res;
    res.num_var = P.size();
    res.num_phen = num_phen;
    res.max_level = max_level;
    ParentSetIndices pix(P);
    std::vector<int> ordered_ixs = set_to_vec(P);

    for (auto i : ordered_ixs)
    {
        for (auto j : ordered_ixs)
        {
            res.G.push_back(G[i * num_var + j]);
            res.C.push_back(C[i * num_var + j]);
            int s_ij_size = 0;

            for (int l = 0; l < max_level; l++)
            {
                int sep_elem = S[(i * num_var + j) * max_level + l];
                if ((sep_elem != -1) && (P.contains(sep_elem)))
                {
                    res.S.push_back(pix.get_new_ix(sep_elem));
                    ++s_ij_size;
                }
            }

            for (int l = s_ij_size; l < max_level; l++)
            {
                res.S.push_back(-1);
            }
        }
    }

    return res;
}