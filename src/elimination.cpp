// Elimination orders and treewidth bounds. None of this code colors anything:
// it only chooses the tree decomposition that the projection DP then uses.
#include "detail.hpp"
#include <bit>
#include <set>
#include <tuple>
namespace td3 {
namespace {
// Symmetric adjacency bitsets of the (later filled) graph, without self-loops.
struct BitGraph {
    std::size_t n, words;
    std::vector<std::vector<std::uint64_t>> a;
    BitGraph(const Graph &g) : n(g.names.size()), words((n + 63) / 64), a(n, std::vector<std::uint64_t>(words)) {
        for (auto [u, v] : g.edges)
            if (u != v) {
                set(u, v);
                set(v, u);
            }
    }
    void set(int u, int v) {
        a[u][static_cast<std::size_t>(v) / 64] |= 1ULL << (v % 64);
    }
    void reset(int u, int v) {
        a[u][static_cast<std::size_t>(v) / 64] &= ~(1ULL << (v % 64));
    }
    bool has(int u, int v) const {
        return (a[u][static_cast<std::size_t>(v) / 64] >> (v % 64)) & 1ULL;
    }
    std::vector<int> neighbors(int v) const {
        std::vector<int> out;
        for (std::size_t w = 0; w < words; ++w) {
            auto bits = a[v][w];
            while (bits) {
                out.push_back(static_cast<int>(w * 64 + static_cast<std::size_t>(std::countr_zero(bits))));
                bits &= bits - 1;
            }
        }
        return out;
    }
    // Make the neighbourhood of v a clique, then delete v (one elimination step).
    // Calls touched(x) for every vertex whose fill or degree score may change.
    template <class F> std::vector<int> eliminate(int v, Budget &budget, F touched) {
        auto ns = neighbors(v);
        for (int u : ns)
            touched(u);
        for (std::size_t i = 0; i < ns.size(); ++i)
            for (std::size_t j = i + 1; j < ns.size(); ++j) {
                int u = ns[i], w = ns[j];
                if (has(u, w))
                    continue;
                budget.check();
                // A new fill edge changes scores only at its endpoints or their common
                // neighbours; marking both endpoint neighbourhoods is a safe superset.
                for (int x : neighbors(u))
                    touched(x);
                for (int x : neighbors(w))
                    touched(x);
                set(u, w);
                set(w, u);
            }
        for (int u : ns)
            reset(u, v);
        std::fill(a[v].begin(), a[v].end(), 0);
        return ns;
    }
};
std::uint64_t bitgraph_bytes(std::size_t n) {
    return checked_product(checked_product(n, (n + 63) / 64), 8) + checked_product(n, 192);
}
} // namespace
std::vector<int> greedy_order(const Graph &g, Budget &budget, Greedy rule, const std::vector<int> *tie_rank) {
    const std::size_t n = g.names.size();
    detail::Memory memory(budget);
    memory.add(bitgraph_bytes(n));
    BitGraph a(g);
    std::vector<bool> active(n, true);
    // Score: (fill-in, degree, tie rank, vertex) for min-fill; (degree, 0, rank, vertex)
    // for min-degree. The smallest key is eliminated next.
    using Key = std::tuple<std::uint64_t, std::size_t, int, int>;
    std::vector<Key> scores(n);
    std::set<Key> queue;
    auto score = [&](int v) {
        auto ns = a.neighbors(v);
        std::uint64_t fill = 0;
        if (rule == Greedy::MinFill)
            for (std::size_t i = 0; i < ns.size(); ++i) {
                budget.check();
                for (std::size_t j = i + 1; j < ns.size(); ++j)
                    fill += !a.has(ns[i], ns[j]);
            }
        int rank = tie_rank ? (*tie_rank)[v] : v;
        return rule == Greedy::MinFill ? Key{fill, ns.size(), rank, v} : Key{ns.size(), 0, rank, v};
    };
    for (std::size_t v = 0; v < n; ++v) {
        scores[v] = score(static_cast<int>(v));
        queue.insert(scores[v]);
    }
    std::vector<int> order;
    order.reserve(n);
    std::vector<int> dirty;
    std::vector<bool> marked(n, false);
    for (std::size_t step = 0; step < n; ++step) {
        budget.check();
        int v = std::get<3>(*queue.begin());
        queue.erase(queue.begin());
        a.eliminate(v, budget, [&](int x) {
            if (!marked[x]) {
                marked[x] = true;
                dirty.push_back(x);
            }
        });
        active[v] = false;
        order.push_back(v);
        for (int u : dirty) {
            marked[u] = false;
            if (active[u]) {
                queue.erase(scores[u]);
                scores[u] = score(u);
                queue.insert(scores[u]);
            }
        }
        dirty.clear();
    }
    return order;
}
std::vector<int> mcs_order(const Graph &g, Budget &budget) {
    // Maximum cardinality search visits a vertex with the most visited neighbours;
    // the reverse visit order is an elimination order.
    const std::size_t n = g.names.size();
    std::vector<int> weight(n, 0), visit;
    std::vector<bool> done(n, false);
    std::set<std::pair<int, int>> queue; // (-weight, vertex)
    for (std::size_t v = 0; v < n; ++v)
        queue.insert({0, static_cast<int>(v)});
    while (!queue.empty()) {
        budget.check();
        int v = queue.begin()->second;
        queue.erase(queue.begin());
        done[v] = true;
        visit.push_back(v);
        for (int u : g.adj[v])
            if (u != v && !done[u]) {
                queue.erase({-weight[u], u});
                queue.insert({-++weight[u], u});
            }
    }
    return {visit.rbegin(), visit.rend()};
}
std::vector<int> exact_order(const Graph &g, Budget &budget, int cap) {
    // Bodlaender, Fomin, Koster, Kratsch and Thilikos (2006): for a set S of vertices
    // eliminated first, TW(S) = min over v in S of max(TW(S - v), |Q(S - v, v)|), where
    // Q(T, v) is the set of vertices outside T + v reachable from v through T. These
    // are exactly v's neighbours in the filled graph when it is eliminated after T.
    const std::size_t n = g.names.size();
    if (n > 26)
        throw ResourceError("exact treewidth search supports at most 26 vertices per component");
    if (n == 0 || cap <= 0)
        return {};
    std::vector<std::uint32_t> nb(n, 0);
    for (auto [u, v] : g.edges)
        if (u != v) {
            nb[u] |= 1U << v;
            nb[v] |= 1U << u;
        }
    const std::uint32_t full = n == 32 ? ~0U : (1U << n) - 1;
    const std::size_t states = static_cast<std::size_t>(full) + 1;
    detail::Memory memory(budget);
    memory.add(checked_product(states, 2));
    // tw[S] is min(cap, best width for eliminating S first); cap means "no width < cap".
    std::vector<std::uint8_t> tw(states), choice(states, 0xff);
    const auto limit = static_cast<std::uint8_t>(std::min(cap, 255));
    tw[0] = 0;
    for (std::uint32_t s = 1; s <= full && s != 0; ++s) {
        if ((s & 0x3fff) == 0)
            budget.check();
        std::uint8_t best = limit, chosen = 0xff;
        for (std::uint32_t rest = s; rest; rest &= rest - 1) {
            int v = std::countr_zero(rest);
            std::uint32_t t = s & ~(1U << v);
            if (tw[t] >= best)
                continue;
            // Breadth-first closure from v through T; collect the boundary outside S.
            std::uint32_t inside = 1U << v, frontier = inside, out = 0;
            while (frontier) {
                std::uint32_t reach = 0;
                for (std::uint32_t f = frontier; f; f &= f - 1)
                    reach |= nb[std::countr_zero(f)];
                out |= reach & ~s;
                frontier = reach & t & ~inside;
                inside |= frontier;
            }
            auto value = static_cast<std::uint8_t>(std::max<int>(tw[t], std::popcount(out)));
            if (value < best) {
                best = value;
                chosen = static_cast<std::uint8_t>(v);
            }
        }
        tw[s] = best;
        choice[s] = chosen;
        if (s == full)
            break;
    }
    if (tw[full] >= limit)
        return {};
    std::vector<int> reversed;
    for (std::uint32_t s = full; s; s &= ~(1U << choice[s]))
        reversed.push_back(choice[s]);
    return {reversed.rbegin(), reversed.rend()};
}
Decomposition from_elimination_order(const Graph &g, const std::vector<int> &order, Budget &budget) {
    const std::size_t n = g.names.size();
    if (order.size() != n)
        throw std::logic_error("elimination order is not a permutation");
    std::vector<int> position(n, -1);
    for (std::size_t i = 0; i < n; ++i) {
        if (order[i] < 0 || static_cast<std::size_t>(order[i]) >= n || position[order[i]] >= 0)
            throw std::logic_error("elimination order is not a permutation");
        position[order[i]] = static_cast<int>(i);
    }
    detail::Memory memory(budget);
    memory.add(bitgraph_bytes(n));
    BitGraph a(g);
    Decomposition d;
    d.bags.reserve(n);
    d.names.reserve(n);
    for (std::size_t step = 0; step < n; ++step) {
        budget.check();
        int v = order[step];
        auto bag = a.eliminate(v, budget, [](int) {});
        memory.add(checked_product(bag.size() + 1, 4) + 64);
        bag.push_back(v);
        std::sort(bag.begin(), bag.end());
        d.bags.push_back(std::move(bag));
        d.names.push_back("B" + std::to_string(step));
    }
    // bag(v) - v is a clique of later vertices, all present in the bag of the earliest
    // of them; attach there. Bags without later vertices are component roots.
    std::vector<int> component_roots;
    for (std::size_t t = 0; t < n; ++t) {
        int v = order[t], p = -1;
        for (int u : d.bags[t])
            if (u != v && (p < 0 || position[u] < position[p]))
                p = u;
        if (p >= 0)
            d.edges.push_back({static_cast<int>(t), position[p]});
        else
            component_roots.push_back(static_cast<int>(t));
    }
    for (std::size_t i = 1; i < component_roots.size(); ++i)
        d.edges.push_back({component_roots[i - 1], component_roots[i]});
    return d;
}
int treewidth_lower_bound(const Graph &g, Budget &budget) {
    // Minor-min-width: treewidth never increases under taking minors, and every graph
    // has treewidth at least its minimum degree. Repeatedly record the minimum degree,
    // then contract a minimum-degree vertex into its minimum-degree neighbour.
    const std::size_t n = g.names.size();
    if (n == 0)
        return -1;
    detail::Memory memory(budget);
    memory.add(checked_product(n, 128) + checked_product(g.edges.size(), 96));
    std::vector<std::set<int>> adj(n);
    for (auto [u, v] : g.edges)
        if (u != v) {
            adj[u].insert(v);
            adj[v].insert(u);
        }
    std::set<std::pair<std::size_t, int>> queue;
    std::vector<std::size_t> key(n); // degree under which each vertex is queued
    for (std::size_t v = 0; v < n; ++v) {
        key[v] = adj[v].size();
        queue.insert({key[v], static_cast<int>(v)});
    }
    int bound = 0;
    auto requeue = [&](int x) {
        queue.erase({key[x], x});
        key[x] = adj[x].size();
        queue.insert({key[x], x});
    };
    while (queue.size() > 1) {
        budget.check();
        auto [degree, v] = *queue.begin();
        queue.erase(queue.begin());
        bound = std::max(bound, static_cast<int>(degree));
        if (degree == 0)
            continue;
        int u = -1;
        for (int w : adj[v])
            if (u < 0 || adj[w].size() < adj[u].size())
                u = w;
        auto neighbours = std::move(adj[v]);
        adj[v].clear();
        for (int w : neighbours) {
            adj[w].erase(v);
            if (w != u) {
                adj[w].insert(u);
                adj[u].insert(w);
            }
        }
        for (int w : neighbours)
            requeue(w);
    }
    return bound;
}
} // namespace td3
