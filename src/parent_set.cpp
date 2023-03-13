#include <mps/parent_set.h>
#include <queue>
#include <set>
#include <vector>

std::set<int> parent_set(
    const std::vector<int> &G,
    const size_t num_var,
    const size_t num_markers,
    const int max_depth)
{
    std::set<int> global_parent_set;

    // do bfs, starting from each of the markers.
    for (int start_ix = num_markers; start_ix < num_var; start_ix++)
    {
        std::set<int> curr_parent_set;
        std::queue<int> q;
        std::queue<int> next_q;
        curr_parent_set.insert(start_ix);
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
                        q.push(cix);
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