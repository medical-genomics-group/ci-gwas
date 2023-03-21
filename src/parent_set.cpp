#include <mps/parent_set.h>

#include <algorithm>
#include <queue>
#include <unordered_set>
#include <vector>

std::vector<int> parent_set(
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

    std::vector<int> res;
    res.reserve(global_parent_set.size());
    for (auto it = global_parent_set.begin(); it != global_parent_set.end();)
    {
        res.push_back(std::move(global_parent_set.extract(it++).value()));
    }

    std::sort(res.begin(), res.end());

    return res;
}