#include "td3/td3.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <random>
#include <iostream>
#include <sstream>
using namespace td3;
struct Baseline {
    bool colorable = false;
    std::uint64_t comparisons = 0;
    double ms = 0;
};
// Benchmark-only reproduction of nested parent/child table scans. No use in the solver.
Baseline naive(const Graph &g, const Decomposition &d) {
    auto start = Clock::now();
    Baseline result;
    std::vector<std::vector<std::vector<int>>> tables(d.bags.size());
    for (int t = static_cast<int>(d.bags.size()) - 1; t >= 0; --t) {
        auto &bag = d.bags[t];
        std::vector<std::pair<int, int>> shared, local_edges;
        if (t + 1 < static_cast<int>(d.bags.size()))
            for (std::size_t i = 0; i < bag.size(); ++i)
                for (std::size_t j = 0; j < d.bags[t + 1].size(); ++j)
                    if (bag[i] == d.bags[t + 1][j])
                        shared.push_back({static_cast<int>(i), static_cast<int>(j)});
        for (std::size_t i = 0; i < bag.size(); ++i)
            for (std::size_t j = i + 1; j < bag.size(); ++j)
                if (std::binary_search(g.adj[bag[i]].begin(), g.adj[bag[i]].end(), bag[j]))
                    local_edges.push_back({static_cast<int>(i), static_cast<int>(j)});
        for (Code code = 0; code < checked_power3(bag.size()); ++code) {
            Code x = code;
            std::vector<int> colors(bag.size());
            for (auto &c : colors) {
                c = static_cast<int>(x % 3);
                x /= 3;
            }
            bool proper = true;
            for (auto [i, j] : local_edges)
                if (colors[i] == colors[j]) {
                    proper = false;
                    break;
                }
            if (!proper)
                continue;
            bool extendable = t + 1 == static_cast<int>(d.bags.size());
            if (!extendable)
                for (const auto &child : tables[t + 1]) {
                    ++result.comparisons;
                    bool match = true;
                    for (auto [i, j] : shared)
                        if (colors[i] != child[j]) {
                            match = false;
                            break;
                        }
                    if (match) {
                        extendable = true;
                        break;
                    }
                }
            if (extendable)
                tables[t].push_back(std::move(colors));
        }
    }
    result.colorable = tables.empty() || !tables[0].empty();
    result.ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    return result;
}
void measure(const std::string &family, int width, int bags, int repeats) {
    int n = bags + width;
    std::ostringstream input;
    for (int v = 0; v < n; ++v)
        input << v << ",\n";
    if (family == "banded_3colorable")
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n && v <= u + width; ++v)
                if (u % 3 != v % 3 && (v == u + 1 || ((u * 31 + v * 17) % 5 == 0)))
                    input << u << ',' << v << '\n';
    std::istringstream in(input.str());
    auto g = read_graph(in);
    Decomposition d;
    for (int t = 0; t < bags; ++t) {
        d.names.push_back("B" + std::to_string(t));
        std::vector<int> bag;
        for (int v = t; v <= t + width; ++v)
            bag.push_back(v);
        d.bags.push_back(bag);
        if (t)
            d.edges.push_back({t - 1, t});
    }
    validate_decomposition(g, d);
    std::vector<double> old_times, new_times, dp_times, preparation_times;
    Baseline old;
    Result fresh;
    Options options;
    options.preprocess = false;
    for (int i = 0; i < repeats; ++i) {
        old = naive(g, d);
        fresh = solve(g, options, &d);
        if (!old.colorable || fresh.status != Status::Colorable ||
            !verify_coloring(g, fresh.coloring))
            throw std::runtime_error("benchmark verdict mismatch");
        old_times.push_back(old.ms);
        new_times.push_back(fresh.stats.preparation_ms + fresh.stats.dp_ms);
        dp_times.push_back(fresh.stats.dp_ms);
        preparation_times.push_back(fresh.stats.preparation_ms);
    }
    std::sort(old_times.begin(), old_times.end());
    std::sort(new_times.begin(), new_times.end());
    std::sort(dp_times.begin(), dp_times.end());
    std::sort(preparation_times.begin(), preparation_times.end());
    if (family == "table_stress" && old.comparisons <= fresh.stats.lookups)
        throw std::runtime_error("projection did not reduce table operations");
    std::cout << family << ',' << width << ',' << bags << ',' << n << ',' << g.edges.size() << ','
              << old_times[repeats / 2] << ',' << new_times[repeats / 2] << ','
              << preparation_times[repeats / 2] << ',' << dp_times[repeats / 2] << ','
              << old.comparisons << ',' << fresh.stats.lookups << ',' << fresh.stats.projections
              << ',' << fresh.stats.peak_accounted_bytes << '\n';
}
// Decomposition quality: the original single min-fill order without simplification
// ("before") against the new construction methods, on identical graphs.
Graph synthetic(const std::string &kind, unsigned seed) {
    std::ostringstream s;
    std::mt19937 rng(seed);
    int n = 0;
    std::vector<std::pair<int, int>> edges;
    if (kind == "grid8x8") {
        n = 64;
        for (int r = 0; r < 8; ++r)
            for (int c = 0; c < 8; ++c) {
                if (c < 7)
                    edges.push_back({r * 8 + c, r * 8 + c + 1});
                if (r < 7)
                    edges.push_back({r * 8 + c, r * 8 + c + 8});
            }
    } else if (kind == "planted60") { // random 3-colorable graph, average degree about 4
        n = 60;
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v)
                if (u % 3 != v % 3 && rng() % 1000 < 100)
                    edges.push_back({u, v});
    } else if (kind == "banded400") { // known path decomposition of width 4
        n = 400;
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n && v <= u + 4; ++v)
                if (u % 3 != v % 3 && (v == u + 1 || rng() % 3 == 0))
                    edges.push_back({u, v});
    }
    for (int v = 0; v < n; ++v)
        s << v << ",\n";
    for (auto [u, v] : edges)
        s << u << ',' << v << '\n';
    std::istringstream in(s.str());
    return read_graph(in);
}
void decompositions(const std::string &source, int repeats) {
    std::cout << "graph,vertices,edges,configuration,width,width_status,lower_bound,bags,"
                 "dp_state_bound,construction_median_ms,dp_median_ms,verdict\n";
    std::vector<std::pair<std::string, Graph>> graphs;
    for (auto name : {"graph00", "graph01", "graph02", "graph03", "graph04", "graph05", "graph06",
                      "graph07", "graph08", "graph09", "graph10", "graph11", "graph12", "graph13",
                      "graph15_3colorable", "graph20_3colorable", "groetzsch"}) {
        std::ifstream in(source + "/graphs/" + name + ".csv");
        if (!in)
            throw std::runtime_error(std::string("missing fixture ") + name);
        graphs.push_back({name, read_graph(in)});
    }
    for (auto kind : {"grid8x8", "planted60", "banded400"})
        graphs.push_back({kind, synthetic(kind, 2026)});
    struct Config {
        const char *name;
        TdMethod method;
        bool simplify;
    };
    const Config configs[] = {{"before:minfill-raw", TdMethod::MinFill, false},
                              {"minfill+simplify", TdMethod::MinFill, true},
                              {"best", TdMethod::Best, true},
                              {"auto", TdMethod::Auto, true}};
    for (auto &[name, g] : graphs) {
        std::string verdict;
        for (const auto &config : configs) {
            Options options;
            options.preprocess = false;
            options.td_method = config.method;
            options.simplify = config.simplify;
            std::vector<double> construction, dp;
            Result r;
            for (int i = 0; i < repeats; ++i) {
                r = solve(g, options);
                if (r.status != Status::Colorable && r.status != Status::NotColorable)
                    throw std::runtime_error(name + ": " + r.message);
                if (r.status == Status::Colorable && !verify_coloring(g, r.coloring))
                    throw std::runtime_error(name + ": invalid coloring");
                construction.push_back(r.stats.decomposition_ms);
                dp.push_back(r.stats.preparation_ms + r.stats.dp_ms + r.stats.reconstruction_ms);
            }
            auto status = status_name(r.status);
            if (verdict.empty())
                verdict = status;
            else if (verdict != status)
                throw std::runtime_error(name + ": verdict differs between configurations");
            std::sort(construction.begin(), construction.end());
            std::sort(dp.begin(), dp.end());
            std::cout << name << ',' << g.names.size() << ',' << g.edges.size() << ','
                      << config.name << ',' << r.decomposition.width() << ','
                      << r.td.width_status << ',' << r.td.lower_bound << ','
                      << r.decomposition.bags.size() << ',' << dp_cost(r.decomposition) << ','
                      << construction[repeats / 2] << ',' << dp[repeats / 2] << ',' << status
                      << '\n';
        }
    }
}
int main(int argc, char **argv) {
    try {
        if (argc > 1 && std::string(argv[1]) == "--decompositions") {
            decompositions(argc > 2 ? argv[2] : TD3_SOURCE_DIR, 5);
            return 0;
        }
        bool quick = argc > 1 && std::string(argv[1]) == "--quick";
        std::cout << "family,width,bags,vertices,edges,naive_median_ms,projection_total_median_ms,"
                     "preparation_median_ms,projection_dp_median_ms,naive_"
                     "assignment_comparisons,projection_lookups,projected_assignments,peak_"
                     "accounted_bytes\n";
        if (quick) {
            measure("table_stress", 3, 5, 1);
            return 0;
        }
        for (int width = 2; width <= 7; ++width) {
            measure("table_stress", width, 8, 3);
            measure("banded_3colorable", width, 8, 3);
        }
        for (int bags : {16, 64, 256})
            measure("banded_3colorable", 4, bags, 3);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
