#include "detail.hpp"
#include <cmath>
#include <memory>
#include <random>
#include <set>
#include <tuple>
namespace td3 {
int Decomposition::width() const {
    int w = -1;
    for (auto &b : bags)
        w = std::max(w, static_cast<int>(b.size()) - 1);
    return w;
}
Decomposition read_decomposition(std::istream &in, const Graph &g) {
    Decomposition d;
    std::unordered_map<std::string, int> ids, vertices;
    for (std::size_t i = 0; i < g.names.size(); ++i)
        vertices[g.names[i]] = static_cast<int>(i);
    std::map<int, std::vector<int>> declarations;
    std::set<std::pair<int, int>> edges;
    detail::rows(in, [&](const detail::Row &r) {
        int a = detail::intern(r.left, ids, d.names);
        if (r.right.empty()) {
            if (!r.third)
                throw InputError("line " + std::to_string(r.line) +
                                 ": bag requires a third column (possibly empty)");
            std::vector<int> bag;
            if (!r.label.empty()) {
                std::size_t start = 0;
                for (;;) {
                    auto end = r.label.find(';', start);
                    auto name = r.label.substr(start, end == std::string::npos ? end : end - start);
                    auto it = vertices.find(name);
                    if (it == vertices.end())
                        throw InputError("line " + std::to_string(r.line) +
                                         ": unknown bag vertex " + name);
                    bag.push_back(it->second);
                    if (end == std::string::npos)
                        break;
                    start = end + 1;
                }
            }
            std::sort(bag.begin(), bag.end());
            if (std::adjacent_find(bag.begin(), bag.end()) != bag.end())
                throw InputError("duplicate vertex in bag");
            if (!declarations.emplace(a, bag).second)
                throw InputError("duplicate bag declaration");
        } else {
            if (!r.label.empty())
                throw InputError("tree edges must be unlabelled");
            int b = detail::intern(r.right, ids, d.names);
            if (a > b)
                std::swap(a, b);
            if (!edges.emplace(a, b).second)
                throw InputError("duplicate decomposition edge");
        }
    });
    if (declarations.size() != d.names.size())
        throw InputError("every tree node requires a bag declaration");
    d.bags.resize(d.names.size());
    for (auto &[i, bag] : declarations)
        d.bags[i] = std::move(bag);
    d.edges.assign(edges.begin(), edges.end());
    return d;
}
void write_decomposition(const Decomposition &d, const Graph &g, std::ostream &out) {
    for (auto [a, b] : d.edges)
        out << d.names[a] << ',' << d.names[b] << ",\n";
    for (std::size_t i = 0; i < d.bags.size(); ++i) {
        out << d.names[i] << ",,";
        bool first = true;
        for (int v : d.bags[i]) {
            if (!first)
                out << ';';
            first = false;
            out << g.names[v];
        }
        out << '\n';
    }
}
void validate_decomposition(const Graph &g, const Decomposition &d, Budget *budget) {
    auto tick = [&] {
        if (budget)
            budget->check();
    };
    tick();
    const auto n = g.names.size(), b = d.bags.size();
    std::unique_ptr<Reservation> workspace;
    if (budget) {
        std::uint64_t bytes = checked_product(n, 64) + checked_product(b, 192);
        for (auto &bag : d.bags)
            bytes += checked_product(bag.size(), 16);
        workspace = std::make_unique<Reservation>(*budget, bytes);
    }
    if (d.names.size() != b)
        throw InputError("bag name count mismatch");
    if (b == 0) {
        if (n || !d.edges.empty())
            throw InputError("empty decomposition does not cover graph");
        return;
    }
    if (d.edges.size() != b - 1)
        throw InputError("decomposition must have exactly bags-1 tree edges");
    std::vector<int> roots(b);
    std::iota(roots.begin(), roots.end(), 0);
    auto find = [&](int x) {
        while (roots[x] != x) {
            roots[x] = roots[roots[x]];
            x = roots[x];
        }
        return x;
    };
    std::vector<std::vector<int>> bags = d.bags, occurrences(n);
    std::set<std::string> names;
    for (std::size_t t = 0; t < b; ++t) {
        tick();
        if (!names.insert(d.names[t]).second)
            throw InputError("duplicate bag name");
        std::sort(bags[t].begin(), bags[t].end());
        if (std::adjacent_find(bags[t].begin(), bags[t].end()) != bags[t].end())
            throw InputError("duplicate bag vertex");
        for (int v : bags[t]) {
            if (v < 0 || static_cast<std::size_t>(v) >= n)
                throw InputError("bag vertex outside graph");
            occurrences[v].push_back(static_cast<int>(t));
        }
    }
    std::vector<std::size_t> shared_edges(n, 0);
    for (auto [u, v] : d.edges) {
        tick();
        if (u < 0 || v < 0 || static_cast<std::size_t>(u) >= b || static_cast<std::size_t>(v) >= b)
            throw InputError("tree endpoint outside decomposition");
        int ru = find(u), rv = find(v);
        if (ru == rv)
            throw InputError("decomposition contains a cycle");
        roots[ru] = rv;
        auto i = bags[u].begin(), j = bags[v].begin();
        while (i != bags[u].end() && j != bags[v].end()) {
            if (*i < *j)
                ++i;
            else if (*j < *i)
                ++j;
            else {
                ++shared_edges[*i];
                ++i;
                ++j;
            }
        }
    }
    for (std::size_t v = 0; v < n; ++v) {
        if (occurrences[v].empty())
            throw InputError("vertex not covered: " + g.names[v]);
        if (shared_edges[v] + 1 != occurrences[v].size())
            throw InputError("running intersection fails for vertex " + g.names[v]);
    }
    for (auto [u, v] : g.edges) {
        tick();
        bool covered = false;
        const auto &a = occurrences[u];
        const auto &c = occurrences[v];
        std::size_t i = 0, j = 0;
        while (i < a.size() && j < c.size()) {
            if (a[i] == c[j]) {
                covered = true;
                break;
            }
            if (a[i] < c[j])
                ++i;
            else
                ++j;
        }
        if (!covered)
            throw InputError("edge not covered: " + g.names[u] + "," + g.names[v]);
    }
}
Decomposition min_fill(const Graph &g, Budget &budget) {
    return from_elimination_order(g, greedy_order(g, budget, Greedy::MinFill), budget);
}
double dp_cost(const Decomposition &d) {
    double cost = 0;
    for (auto &bag : d.bags)
        cost += std::pow(3.0, static_cast<double>(bag.size()));
    return cost;
}
std::size_t contract_contained_bags(Decomposition &d) {
    // Contracting a tree edge s-t and labelling the merged node with the union of both
    // bags always yields a tree decomposition. When B_s is a subset of B_t the union is
    // B_t, so the width cannot grow. Each contraction removes one bag, so this ends.
    const std::size_t b = d.bags.size();
    for (auto &bag : d.bags)
        std::sort(bag.begin(), bag.end());
    std::vector<std::set<int>> adj(b);
    for (auto [u, v] : d.edges) {
        adj[u].insert(v);
        adj[v].insert(u);
    }
    std::vector<bool> alive(b, true);
    std::vector<std::pair<int, int>> work(d.edges.begin(), d.edges.end());
    auto contains = [&](int big, int small) {
        return std::includes(d.bags[big].begin(), d.bags[big].end(), d.bags[small].begin(),
                             d.bags[small].end());
    };
    std::size_t removed = 0;
    while (!work.empty()) {
        int s = work.back().first, t = work.back().second;
        work.pop_back();
        if (!alive[s] || !alive[t] || !adj[s].contains(t))
            continue;
        if (!contains(t, s)) {
            if (!contains(s, t))
                continue;
            std::swap(s, t);
        }
        // Merge s into t: every other neighbour r of s becomes a neighbour of t.
        for (int r : adj[s])
            if (r != t) {
                adj[r].erase(s);
                adj[r].insert(t);
                adj[t].insert(r);
                work.push_back({r, t});
            }
        adj[t].erase(s);
        adj[s].clear();
        alive[s] = false;
        ++removed;
    }
    if (!removed)
        return 0;
    std::vector<int> index(b, -1);
    Decomposition out;
    for (std::size_t t = 0; t < b; ++t)
        if (alive[t]) {
            index[t] = static_cast<int>(out.bags.size());
            out.bags.push_back(std::move(d.bags[t]));
            out.names.push_back(std::move(d.names[t]));
        }
    for (std::size_t t = 0; t < b; ++t)
        for (int u : adj[t])
            if (static_cast<int>(t) < u)
                out.edges.push_back({index[t], index[u]});
    d = std::move(out);
    return removed;
}
void arrange_tree(Decomposition &d) {
    const std::size_t b = d.bags.size();
    if (b == 0)
        return;
    for (auto &bag : d.bags)
        std::sort(bag.begin(), bag.end());
    // Edges with an empty separator carry no information: removing them leaves a forest
    // of valid decompositions of disjoint vertex sets, which may be re-linked freely.
    std::vector<std::vector<int>> adj(b);
    for (auto [u, v] : d.edges) {
        std::vector<int> shared;
        std::set_intersection(d.bags[u].begin(), d.bags[u].end(), d.bags[v].begin(),
                              d.bags[v].end(), std::back_inserter(shared));
        if (!shared.empty()) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
    }
    std::vector<int> piece(b, -1);
    std::vector<std::vector<int>> pieces;
    for (std::size_t first = 0; first < b; ++first) {
        if (piece[first] >= 0)
            continue;
        const int id = static_cast<int>(pieces.size());
        std::vector<int> nodes{static_cast<int>(first)};
        piece[first] = id;
        for (std::size_t i = 0; i < nodes.size(); ++i)
            for (int u : adj[nodes[i]])
                if (piece[u] < 0) {
                    piece[u] = id;
                    nodes.push_back(u);
                }
        pieces.push_back(std::move(nodes));
    }
    // A centre (minimum eccentricity) is found by peeling leaves layer by layer.
    std::vector<int> degree(b), centers;
    std::vector<bool> peeled(b, false);
    for (auto &nodes : pieces) {
        std::vector<int> layer;
        for (int x : nodes) {
            degree[x] = static_cast<int>(adj[x].size());
            if (degree[x] <= 1)
                layer.push_back(x);
        }
        std::size_t remaining = nodes.size();
        while (remaining > 2) {
            std::vector<int> next;
            for (int x : layer) {
                peeled[x] = true;
                --remaining;
                for (int y : adj[x])
                    if (!peeled[y] && --degree[y] == 1)
                        next.push_back(y);
            }
            layer = std::move(next);
        }
        int center = -1;
        for (int x : nodes)
            if (!peeled[x] && (center < 0 || x < center))
                center = x;
        centers.push_back(center);
    }
    std::size_t main = 0;
    for (std::size_t p = 1; p < pieces.size(); ++p)
        if (pieces[p].size() > pieces[main].size())
            main = p;
    const int root = centers[main];
    for (std::size_t p = 0; p < pieces.size(); ++p)
        if (p != main) {
            adj[root].push_back(centers[p]);
            adj[centers[p]].push_back(root);
        }
    for (auto &a : adj)
        std::sort(a.begin(), a.end());
    // Breadth-first renumbering: B0 is the root, and every parent precedes its children.
    std::vector<int> order{root}, index(b, -1), parent(b, -1);
    index[root] = 0;
    for (std::size_t i = 0; i < order.size(); ++i)
        for (int u : adj[order[i]])
            if (index[u] < 0) {
                index[u] = static_cast<int>(order.size());
                parent[u] = order[i];
                order.push_back(u);
            }
    Decomposition out;
    for (std::size_t i = 0; i < order.size(); ++i) {
        out.bags.push_back(std::move(d.bags[order[i]]));
        out.names.push_back("B" + std::to_string(i));
        if (i)
            out.edges.push_back({index[parent[order[i]]], static_cast<int>(i)});
    }
    d = std::move(out);
}
namespace {
std::vector<std::vector<int>> components(const Graph &g) {
    std::vector<std::vector<int>> out;
    std::vector<bool> seen(g.names.size(), false);
    for (std::size_t first = 0; first < g.names.size(); ++first)
        if (!seen[first]) {
            std::vector<int> c{static_cast<int>(first)};
            seen[first] = true;
            for (std::size_t i = 0; i < c.size(); ++i)
                for (int u : g.adj[c[i]])
                    if (!seen[u]) {
                        seen[u] = true;
                        c.push_back(u);
                    }
            std::sort(c.begin(), c.end());
            out.push_back(std::move(c));
        }
    return out;
}
Graph induced(const Graph &g, const std::vector<int> &vertices) {
    Graph sub;
    std::unordered_map<int, int> mapped;
    for (int v : vertices) {
        mapped[v] = static_cast<int>(sub.names.size());
        sub.names.push_back(g.names[v]);
        sub.labels.emplace_back();
    }
    sub.adj.resize(vertices.size());
    for (int v : vertices)
        for (int u : g.adj[v]) {
            auto it = mapped.find(u);
            if (it == mapped.end())
                continue;
            sub.adj[mapped[v]].push_back(it->second);
            if (v < u) {
                sub.edges.push_back({mapped[v], it->second});
                sub.edge_labels.emplace_back();
            } else if (v == u)
                sub.self_loop = true;
        }
    return sub;
}
std::vector<int> shuffled_ranks(std::size_t n, unsigned seed) {
    // Fisher-Yates driven by mt19937, whose output (unlike std::shuffle) is identical on
    // every standard library, so seeded runs are reproducible across compilers.
    std::vector<int> rank(n);
    std::iota(rank.begin(), rank.end(), 0);
    std::mt19937 rng(seed);
    for (std::size_t i = n; i > 1; --i)
        std::swap(rank[i - 1], rank[rng() % i]);
    return rank;
}
struct Built {
    Decomposition d;
    std::vector<int> order;
    TdCandidate candidate;
    int raw_width = -1;
};
} // namespace
Decomposition build_decomposition(const Graph &g, const Options &options, Budget &budget,
                                  TdInfo &info) {
    const std::size_t n = g.names.size();
    info = TdInfo{};
    auto parts = components(g);
    info.components = parts.size();
    info.lower_bound = treewidth_lower_bound(g, budget);
    auto evaluate = [&](const std::string &name, Clock::time_point start, std::vector<int> order) {
        Built built;
        built.d = from_elimination_order(g, order, budget);
        validate_decomposition(g, built.d, &budget);
        built.raw_width = built.d.width();
        built.candidate.raw_bags = built.d.bags.size();
        if (options.simplify) {
            contract_contained_bags(built.d);
            arrange_tree(built.d);
            validate_decomposition(g, built.d, &budget);
            if (built.d.width() > built.raw_width)
                throw std::logic_error("simplification increased the width");
        }
        built.order = std::move(order);
        auto &c = built.candidate;
        c.name = name;
        c.width = built.d.width();
        c.bags = built.d.bags.size();
        c.cost = dp_cost(built.d);
        c.ms = detail::elapsed(start);
        return built;
    };
    auto better = [](const TdCandidate &a, const TdCandidate &b) {
        return std::tie(a.width, a.cost, a.bags) < std::tie(b.width, b.cost, b.bags);
    };
    Built best;
    std::size_t best_index = 0;
    auto consider = [&](Built built) {
        info.candidates.push_back(built.candidate);
        if (info.candidates.size() == 1 || better(built.candidate, best.candidate)) {
            best = std::move(built);
            best_index = info.candidates.size() - 1;
        }
    };
    auto run = [&](const std::string &name, auto make_order) {
        auto start = Clock::now();
        consider(evaluate(name, start, make_order()));
    };
    run("min-fill", [&] { return greedy_order(g, budget, Greedy::MinFill); });
    if (options.td_method != TdMethod::MinFill) {
        run("min-degree", [&] { return greedy_order(g, budget, Greedy::MinDegree); });
        run("max-cardinality", [&] { return mcs_order(g, budget); });
        if (n <= 2000)
            for (int seed = 1; seed <= options.heuristic_seeds; ++seed) {
                auto rank = shuffled_ranks(n, static_cast<unsigned>(seed));
                run("min-fill/seed" + std::to_string(seed),
                    [&] { return greedy_order(g, budget, Greedy::MinFill, &rank); });
                run("min-degree/seed" + std::to_string(seed),
                    [&] { return greedy_order(g, budget, Greedy::MinDegree, &rank); });
            }
    }
    const bool exact_wanted = (options.td_method == TdMethod::Auto && options.exact_limit > 0) ||
                              options.td_method == TdMethod::Exact;
    if (info.lower_bound == best.candidate.width) {
        info.proof = "lower_bound_match";
    } else if (exact_wanted) {
        // Exact search, one connected component at a time: treewidth is the maximum over
        // components, and elimination orders of different components can be concatenated.
        auto start = Clock::now();
        std::vector<int> component(n, -1);
        for (std::size_t c = 0; c < parts.size(); ++c)
            for (int v : parts[c])
                component[v] = static_cast<int>(c);
        std::vector<std::vector<int>> sequences(parts.size());
        for (int v : best.order)
            sequences[component[v]].push_back(v);
        std::vector<int> part_width(parts.size(), -1);
        {
            auto raw = from_elimination_order(g, best.order, budget);
            for (std::size_t t = 0; t < raw.bags.size(); ++t) {
                int c = component[best.order[t]];
                part_width[c] = std::max(part_width[c], static_cast<int>(raw.bags[t].size()) - 1);
            }
        }
        bool improved = false, all_exact = true;
        int exact_max = -1, heuristic_max = -1;
        std::vector<int> order;
        for (std::size_t c = 0; c < parts.size(); ++c) {
            const auto &part = parts[c];
            if (part.size() > static_cast<std::size_t>(std::max(options.exact_limit, 0))) {
                if (options.td_method == TdMethod::Exact)
                    throw ResourceError("exact decomposition: a component has " +
                                        std::to_string(part.size()) +
                                        " vertices, above the exact limit " +
                                        std::to_string(options.exact_limit));
                all_exact = false;
                heuristic_max = std::max(heuristic_max, part_width[c]);
            } else {
                ++info.exact_components;
                auto sub = induced(g, part);
                auto exact = exact_order(sub, budget, part_width[c]);
                if (!exact.empty()) {
                    improved = true;
                    sequences[c].clear();
                    for (int v : exact)
                        sequences[c].push_back(part[v]);
                    part_width[c] = from_elimination_order(sub, exact, budget).width();
                }
                exact_max = std::max(exact_max, part_width[c]);
            }
            order.insert(order.end(), sequences[c].begin(), sequences[c].end());
        }
        if (improved)
            consider(evaluate("exact", start, std::move(order)));
        else {
            // The exact search confirmed the heuristic width; record what it cost.
            TdCandidate confirm;
            confirm.name = "exact search (confirms heuristic)";
            confirm.width = best.candidate.width;
            confirm.ms = detail::elapsed(start);
            info.candidates.push_back(confirm);
        }
        if (all_exact || exact_max >= heuristic_max)
            info.proof = "exact_search";
    }
    info.candidates[best_index].selected = true;
    info.method = best.candidate.name;
    info.width_status = info.proof.empty() ? "upper_bound" : "optimal";
    info.raw_width = best.raw_width;
    info.raw_bags = best.candidate.raw_bags;
    info.merged_bags = best.candidate.raw_bags - best.candidate.bags;
    return std::move(best.d);
}
} // namespace td3
