#include "td3/td3.hpp"
#include <algorithm>
#include <bit>
#include <numeric>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
using namespace td3;
static int assertions = 0;
void require(bool b, const std::string &message) {
    ++assertions;
    if (!b)
        throw std::runtime_error(message);
}
template <class F> void bad(F f) {
    bool caught = false;
    try {
        f();
    } catch (const InputError &) {
        caught = true;
    }
    require(caught, "expected invalid input");
}
Graph parse(const std::string &s) {
    std::istringstream in(s);
    return read_graph(in);
}
Graph make_graph(int n, const std::vector<std::pair<int, int>> &edges) {
    std::ostringstream s;
    for (int i = 0; i < n; ++i)
        s << i << ",\n";
    for (auto [a, b] : edges)
        s << a << "," << b << "\n";
    return parse(s.str());
}
bool oracle(const Graph &g) {
    if (g.self_loop)
        return false;
    std::vector<int> c(g.names.size(), -1);
    std::function<bool(std::size_t)> dfs = [&](std::size_t v) {
        if (v == c.size())
            return true;
        for (int color = 0; color < 3; ++color) {
            bool ok = true;
            for (int u : g.adj[v])
                if (c[u] == color)
                    ok = false;
            if (ok) {
                c[v] = color;
                if (dfs(v + 1))
                    return true;
                c[v] = -1;
            }
        }
        return false;
    };
    return dfs(0);
}
void check(const Graph &g, bool expected, const Options &o = {}, const Decomposition *d = nullptr) {
    auto r = solve(g, o, d);
    require(r.status == (expected ? Status::Colorable : Status::NotColorable),
            "wrong solve status=" + std::to_string(static_cast<int>(r.status)) +
                " n=" + std::to_string(g.names.size()) + " m=" + std::to_string(g.edges.size()) +
                " feasible=" + std::to_string(r.stats.feasible) + " " + r.message);
    if (expected)
        require(verify_coloring(g, r.coloring), "invalid witness");
}
void parser_tests() {
    auto g = parse("a,,name\r\nb,,other\r\nisolated,,third\r\na,b,edge\r\nb,a,edge\r\n");
    require(g.names.size() == 3 && g.edges.size() == 1, "isolates or duplicate edges lost");
    require(g.labels[0] == "name" && g.edge_labels[0] == "edge", "labels not preserved");
    std::ostringstream out;
    write_graph(g, out);
    auto again = parse(out.str());
    require(again.names == g.names && again.edges == g.edges, "CSV round trip");
    for (auto s : {"a>b\n", "a\n", ",b\n", "a;b,c\n", "a b,c\n", "a,,x\nb,\n", "a,b,x\nb,c\n",
                   "a,b,x\nb,a,y\n"})
        bad([&] { parse(s); });
    std::istringstream legacy("4\n1 2 3 4\n1 2\n2 3\n");
    auto old = read_legacy(legacy);
    require(old.names.size() == 4 && old.adj[3].empty(), "legacy isolate dropped");
    for (auto s : {"3\n1 2\n1 2\n", "2\n1 1\n1 3\n", "-2\n"})
        bad([&] {
            std::istringstream in(s);
            read_legacy(in);
        });
}
void td_tests() {
    auto g = parse("a,b\nb,c\nx,\n");
    std::istringstream in("left,right,\nright,iso,\nleft,,a;b\nright,,b;c\niso,,x\n");
    auto d = read_decomposition(in, g);
    validate_decomposition(g, d);
    check(g, true, {}, &d);
    std::ostringstream out;
    write_decomposition(d, g, out);
    std::istringstream input(out.str());
    validate_decomposition(g, read_decomposition(input, g));
    auto wrong = d;
    wrong.bags[0] = {0};
    bad([&] { validate_decomposition(g, wrong); });
    wrong = d;
    wrong.edges.push_back({0, 2});
    bad([&] { validate_decomposition(g, wrong); });
    wrong = d;
    wrong.bags[2] = {3, 3};
    bad([&] { validate_decomposition(g, wrong); });
    auto h = make_graph(3, {});
    Decomposition split{{"a", "b", "c"}, {{0, 1}, {1, 2}, {0, 2}}, {{0, 1}, {1, 2}}};
    bad([&] { validate_decomposition(h, split); });
    Decomposition branching{
        {"r", "a", "b", "c"}, {{1, 0}, {2, 1}, {3, 1}, {4}}, {{0, 1}, {0, 2}, {0, 3}}};
    auto star = make_graph(5, {{0, 1}, {1, 2}, {1, 3}});
    check(star, true, {}, &branching);
    std::istringstream unknown("B,,missing\n");
    bad([&] { read_decomposition(unknown, g); });
}
void solver_tests() {
    check(parse(""), true);
    check(parse("a,\nb,\n"), true);
    check(parse("a,a\n"), false);
    auto triangle = make_graph(3, {{0, 1}, {1, 2}, {0, 2}});
    check(triangle, true);
    auto k4 = make_graph(4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}});
    check(k4, false);
    for (auto mode : {TableMode::Dense, TableMode::Sparse, TableMode::Auto}) {
        Options o;
        o.preprocess = false;
        o.table_mode = mode;
        o.samples = true;
        check(triangle, true, o);
        check(k4, false, o);
        for (int i = 0; i < 14; ++i) {
            std::ostringstream path;
            path << TD3_SOURCE_DIR << "/graphs/graph" << (i < 10 ? "0" : "") << i << ".csv";
            std::ifstream in(path.str());
            require(bool(in), "missing migrated fixture");
            auto g = read_graph(in);
            check(g, i != 7 && i != 8 && i != 12, o);
            check(g, i != 7 && i != 8 && i != 12);
        }
    }
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < 5; ++i) {
        int j = (i + 1) % 5;
        edges.push_back({i, j});
        edges.push_back({i + 5, j});
        edges.push_back({j + 5, i});
        edges.push_back({10, i + 5});
    }
    check(make_graph(11, edges), false);
    Options low;
    low.memory_bytes = 1;
    require(solve(k4, low).status == Status::ResourceLimit, "memory limit called UNSAT");
    low = Options{};
    low.time_limit_seconds = 1e-12;
    require(solve(k4, low).status == Status::ResourceLimit, "timeout called UNSAT");
    bool overflow = false;
    try {
        checked_power3(41);
    } catch (const ResourceError &) {
        overflow = true;
    }
    require(overflow, "power overflow not checked");
    std::mt19937 rng(91234);
    for (int trial = 0; trial < 200; ++trial) {
        int n = 7 + int(rng() % 4);
        edges.clear();
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (rng() % 100 < 35)
                    edges.push_back({i, j});
        auto g = make_graph(n, edges);
        Options o;
        o.preprocess = false;
        check(g, oracle(g), o);
        check(g, oracle(g));
    }
}
void report_tests() {
    auto g = parse("safe,,</script><script>alert(1)</script>\n");
    Options o;
    o.preprocess = false;
    o.samples = true;
    auto r = solve(g, o);
    std::ostringstream out;
    write_report(g, r, out);
    require(out.str().find("type=\"application/json\"") != std::string::npos,
            "report must embed machine-readable results");
    require(out.str().find("</script><script>alert(1)</script>") == std::string::npos,
            "untrusted label escaped report data");
    require(result_json(g, r).find("\\u003c/script\\u003e") != std::string::npos,
            "JSON must be safe inside script element");
    std::vector<std::pair<int, int>> edges;
    for (int component = 0; component < 1000; ++component)
        for (int u = 0; u < 3; ++u)
            for (int v = 3; v < 6; ++v)
                edges.push_back({component * 6 + u, component * 6 + v});
    auto many = make_graph(6000, edges);
    o.preprocess = true;
    auto large = solve(many, o);
    require(large.status == Status::Colorable, "disconnected residual cores");
    std::size_t samples = 0;
    for (auto &b : large.bags)
        samples += b.samples.size();
    require(samples <= 10000, "sample cap must apply across all components");
}
void large_star_test() {
    std::vector<std::pair<int, int>> edges;
    Decomposition d;
    for (int leaf = 0; leaf < 1024; ++leaf) {
        edges.push_back({leaf, 1024});
        d.names.push_back("B" + std::to_string(leaf));
        d.bags.push_back({leaf, 1024});
        if (leaf)
            d.edges.push_back({leaf - 1, leaf});
    }
    auto star = make_graph(1025, edges);
    Options o;
    o.preprocess = false;
    check(star, true, o, &d);
}
// Independent treewidth oracle: the best width over all n! elimination orders.
int brute_treewidth(const Graph &g) {
    const int n = static_cast<int>(g.names.size());
    if (n == 0)
        return -1;
    std::vector<unsigned> nb(n, 0);
    for (auto [u, v] : g.edges)
        if (u != v) {
            nb[u] |= 1U << v;
            nb[v] |= 1U << u;
        }
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    int best = n;
    do {
        auto a = nb;
        unsigned alive = (1U << n) - 1;
        int width = 0;
        for (int v : order) {
            unsigned ns = a[v] & alive;
            width = std::max(width, std::popcount(ns));
            if (width >= best)
                break;
            for (unsigned r = ns; r; r &= r - 1)
                a[std::countr_zero(r)] |= ns & ~(1U << std::countr_zero(r));
            alive &= ~(1U << v);
        }
        best = std::min(best, width);
    } while (std::next_permutation(order.begin(), order.end()));
    return best;
}
bool contained_neighbours(const Decomposition &d) {
    for (auto [a, b] : d.edges) {
        auto x = d.bags[a], y = d.bags[b];
        std::sort(x.begin(), x.end());
        std::sort(y.begin(), y.end());
        if (std::includes(x.begin(), x.end(), y.begin(), y.end()) ||
            std::includes(y.begin(), y.end(), x.begin(), x.end()))
            return true;
    }
    return false;
}
Graph random_graph(std::mt19937 &rng, int n, int percent) {
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (static_cast<int>(rng() % 100) < percent)
                edges.push_back({i, j});
    return make_graph(n, edges);
}
Graph grid(int rows, int cols) {
    std::vector<std::pair<int, int>> edges;
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c) {
            if (c + 1 < cols)
                edges.push_back({r * cols + c, r * cols + c + 1});
            if (r + 1 < rows)
                edges.push_back({r * cols + c, (r + 1) * cols + c});
        }
    return make_graph(rows * cols, edges);
}
void decomposition_quality_tests() {
    Options big;
    Budget budget(big);
    std::mt19937 rng(20261117);
    // Elimination-order decompositions, contraction and re-rooting keep validity and width.
    for (int trial = 0; trial < 400; ++trial) {
        int n = 1 + static_cast<int>(rng() % 14);
        auto g = random_graph(rng, n, static_cast<int>(rng() % 70));
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        for (int i = n; i > 1; --i)
            std::swap(order[i - 1], order[rng() % i]);
        auto d = from_elimination_order(g, order, budget);
        validate_decomposition(g, d);
        require(static_cast<int>(d.bags.size()) == n, "one bag per eliminated vertex");
        auto original = d;
        int width = d.width();
        auto removed = contract_contained_bags(d);
        validate_decomposition(g, d);
        require(d.width() == width, "contraction changed the width");
        require(removed + d.bags.size() == original.bags.size(), "contraction count");
        require(!contained_neighbours(d), "contraction left a contained neighbour");
        require(static_cast<int>(d.bags.size()) <= std::max(n, 1), "reduced TD exceeds n bags");
        for (auto &bag : d.bags) // every surviving bag is one of the original bags
            require(std::find(original.bags.begin(), original.bags.end(), bag) !=
                        original.bags.end(),
                    "contraction invented a bag");
        arrange_tree(d);
        validate_decomposition(g, d);
        require(d.width() == width && !contained_neighbours(d), "arrangement changed bags");
        for (std::size_t i = 0; i < d.names.size(); ++i)
            require(d.names[i] == "B" + std::to_string(i), "breadth-first names");
        for (auto [p, c] : d.edges)
            require(p < c, "parents must precede children");
    }
    // Exact search and the lower bound against the brute-force oracle.
    for (int trial = 0; trial < 300; ++trial) {
        int n = static_cast<int>(rng() % 9);
        auto g = random_graph(rng, n, 15 + static_cast<int>(rng() % 70));
        int tw = brute_treewidth(g);
        require(treewidth_lower_bound(g, budget) <= tw, "lower bound exceeds treewidth");
        if (n) {
            auto order = exact_order(g, budget, n);
            require(!order.empty(), "exact search found nothing below n");
            require(from_elimination_order(g, order, budget).width() == tw, "exact width");
            require(exact_order(g, budget, tw).empty(), "exact search beat the optimum");
        }
        for (auto method : {TdMethod::MinFill, TdMethod::Best, TdMethod::Auto, TdMethod::Exact}) {
            Options o;
            o.td_method = method;
            TdInfo info;
            auto d = build_decomposition(g, o, budget, info);
            validate_decomposition(g, d);
            require(d.width() >= tw, "heuristic below treewidth");
            if (method == TdMethod::Auto || method == TdMethod::Exact)
                require(d.width() == tw && info.width_status == "optimal", "auto not optimal");
            if (info.width_status == "optimal")
                require(d.width() == tw, "claimed optimality is false");
            require(!contained_neighbours(d), "generated TD not simplified");
        }
    }
    // Known treewidths.
    auto cycle = make_graph(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}});
    auto k33 = make_graph(6, {{0, 3}, {0, 4}, {0, 5}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5}});
    std::vector<std::pair<int, int>> k5edges;
    for (int i = 0; i < 5; ++i)
        for (int j = i + 1; j < 5; ++j)
            k5edges.push_back({i, j});
    std::vector<std::pair<Graph, int>> known{{parse(""), -1},      {parse("a,\nb,\n"), 0},
                                             {cycle, 2},           {k33, 3},
                                             {grid(3, 3), 3},      {grid(4, 4), 4},
                                             {make_graph(5, k5edges), 4}};
    for (auto [name, tw] : std::vector<std::pair<std::string, int>>{
             {"graph07", 3}, {"graph11", 4}, {"graph13", 4}, {"groetzsch", 5}}) {
        std::ifstream in(std::string(TD3_SOURCE_DIR) + "/graphs/" + name + ".csv");
        known.push_back({read_graph(in), tw});
    }
    for (auto &[g, tw] : known) {
        Options o;
        o.preprocess = false;
        auto r = solve(g, o);
        require(r.decomposition.width() == tw && r.td.width_status == "optimal",
                "known treewidth " + std::to_string(tw) + " got " +
                    std::to_string(r.decomposition.width()));
    }
    // Method ordering on larger graphs: more search never gives a wider decomposition.
    for (int trial = 0; trial < 20; ++trial) {
        auto g = random_graph(rng, 16 + static_cast<int>(rng() % 5), 20);
        int widths[3];
        TdMethod methods[3] = {TdMethod::MinFill, TdMethod::Best, TdMethod::Auto};
        for (int i = 0; i < 3; ++i) {
            Options o;
            o.td_method = methods[i];
            TdInfo info;
            widths[i] = build_decomposition(g, o, budget, info).width();
            if (i == 2)
                require(info.width_status == "optimal", "auto must prove small components");
        }
        require(widths[1] <= widths[0] && widths[2] <= widths[1], "method ordering");
    }
    // The original single min-fill decomposition is reproducible exactly.
    {
        std::ifstream in(std::string(TD3_SOURCE_DIR) + "/graphs/graph20_3colorable.csv");
        auto g = read_graph(in);
        Options o;
        o.td_method = TdMethod::MinFill;
        o.simplify = false;
        TdInfo info;
        auto d = build_decomposition(g, o, budget, info);
        auto old = min_fill(g, budget);
        require(d.bags == old.bags && d.edges == old.edges, "raw min-fill changed");
        require(info.merged_bags == 0 && d.bags.size() == g.names.size(), "raw means unmerged");
    }
    // Disconnected inputs: isolated vertices hang below the largest piece's centre.
    {
        auto g = make_graph(9, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {5, 6}});
        Options o;
        o.preprocess = false;
        auto r = solve(g, o);
        require(r.status == Status::Colorable && r.td.components == 4, "disconnected solve");
        require(r.decomposition.width() == 2 && r.td.width_status == "optimal", "disconnected tw");
    }
    // Exact mode refuses a too-large component unless the lower bound already proves it.
    {
        auto g = grid(6, 6);
        Options o;
        o.preprocess = false;
        o.td_method = TdMethod::Exact;
        o.exact_limit = 10;
        auto r = solve(g, o);
        require(r.status == Status::ResourceLimit || r.td.proof == "lower_bound_match",
                "exact mode must not guess");
        o.td_method = TdMethod::Auto;
        r = solve(g, o);
        require(r.status == Status::Colorable && r.td.exact_components == 0, "auto fallback");
        require(r.td.width_status != "optimal" || r.td.lower_bound == r.decomposition.width(),
                "unproven grid width labelled optimal");
    }
    // Supplied decompositions stay untouched unless --simplify-td is requested.
    {
        auto g = make_graph(3, {{0, 1}, {1, 2}});
        Decomposition d{{"x", "y", "z"}, {{0, 1}, {0}, {1, 2}}, {{0, 1}, {0, 2}}};
        Options o;
        auto r = solve(g, o, &d);
        require(r.decomposition.bags.size() == 3 && r.td.method == "supplied", "supplied kept");
        o.simplify_supplied = true;
        r = solve(g, o, &d);
        require(r.decomposition.bags.size() == 2 && r.td.merged_bags == 1 &&
                    r.decomposition.names[0] == "x",
                "supplied simplification");
        require(r.status == Status::Colorable && verify_coloring(g, r.coloring), "simplified td");
    }
}
void exhaustive() {
    std::uint64_t count = 0;
    for (int n = 0; n <= 6; ++n) {
        std::vector<std::pair<int, int>> all;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                all.push_back({i, j});
        for (std::uint64_t mask = 0; mask < (1ULL << all.size()); ++mask) {
            std::vector<std::pair<int, int>> e;
            for (std::size_t i = 0; i < all.size(); ++i)
                if (mask & (1ULL << i))
                    e.push_back(all[i]);
            auto g = make_graph(n, e);
            const bool colorable = oracle(g);
            const int tw = brute_treewidth(g);
            Options o;
            o.preprocess = false;
            check(g, colorable, o); // default: auto (heuristics + exact search)
            check(g, colorable);    // with degree-two peeling
            auto r = solve(g, o);
            require(r.decomposition.width() == tw && r.td.width_status == "optimal",
                    "auto decomposition width is not the treewidth");
            require(!contained_neighbours(r.decomposition), "contained bags remain");
            for (auto method : {TdMethod::MinFill, TdMethod::Best}) {
                o.td_method = method;
                for (bool simplify : {true, false}) {
                    o.simplify = simplify;
                    auto h = solve(g, o);
                    require(h.decomposition.width() >= tw, "heuristic below treewidth");
                    require(h.td.width_status != "optimal" || h.decomposition.width() == tw,
                            "false optimality claim");
                    check(g, colorable, o);
                }
            }
            ++count;
        }
    }
    std::cout << "Exhaustive graphs checked: " << count << "\n";
}
int main(int argc, char **) {
    try {
        if (argc > 1)
            exhaustive();
        else {
            parser_tests();
            td_tests();
            solver_tests();
            report_tests();
            large_star_test();
            decomposition_quality_tests();
        }
        std::cout << "PASS " << assertions << " assertions\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
