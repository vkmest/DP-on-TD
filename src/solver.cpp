#include "detail.hpp"
#include <functional>
#include <limits>
#include <queue>
namespace td3 {
void Budget::check() const {
    if (options_.time_limit_seconds > 0 &&
        std::chrono::duration<double>(Clock::now() - start_).count() >= options_.time_limit_seconds)
        throw ResourceError("time limit reached; colorability is unknown");
}
void Budget::acquire(std::uint64_t bytes) {
    check();
    if (bytes > options_.memory_bytes - live_)
        throw ResourceError("accounted working-memory limit reached; colorability is unknown");
    live_ += bytes;
    peak_ = std::max(peak_, live_);
}
void Budget::release(std::uint64_t bytes) {
    live_ -= bytes;
}
std::uint64_t checked_product(std::uint64_t a, std::uint64_t b) {
    if (b && a > std::numeric_limits<std::uint64_t>::max() / b)
        throw ResourceError("size arithmetic overflow");
    return a * b;
}
std::uint64_t checked_power3(std::size_t exponent) {
    std::uint64_t x = 1;
    for (std::size_t i = 0; i < exponent; ++i)
        x = checked_product(x, 3);
    return x;
}
bool verify_coloring(const Graph &g, const std::vector<int> &c) {
    if (c.size() != g.names.size())
        return false;
    for (int color : c)
        if (color < 0 || color > 2)
            return false;
    for (auto [u, v] : g.edges)
        if (c[u] == c[v])
            return false;
    return true;
}
namespace {
constexpr Code absent = std::numeric_limits<Code>::max();
struct Table {
    std::vector<Code> dense;
    std::vector<std::pair<Code, Code>> sparse;
    std::uint64_t count = 0;
    Code get(Code key) const {
        if (!dense.empty())
            return dense[static_cast<std::size_t>(key)];
        auto it = std::lower_bound(sparse.begin(), sparse.end(), key,
                                   [](const auto &entry, Code k) { return entry.first < k; });
        return it != sparse.end() && it->first == key ? it->second : absent;
    }
};
Code project(const std::vector<int> &colors, const std::vector<int> &positions) {
    Code key = 0, power = 1;
    for (int i : positions) {
        key += static_cast<Code>(colors[i]) * power;
        power *= 3;
    }
    return key;
}
std::vector<int> positions(const std::vector<int> &bag, const std::vector<int> &separator) {
    std::vector<int> p;
    for (int v : separator)
        p.push_back(static_cast<int>(std::lower_bound(bag.begin(), bag.end(), v) - bag.begin()));
    return p;
}
bool run_dp(const Graph &g, Decomposition &d, const Options &options, Budget &budget,
            Result &result, std::size_t &remaining_samples) {
    auto preparation_start = Clock::now();
    const auto b = d.bags.size();
    if (b == 0) {
        result.coloring.clear();
        return true;
    }
    detail::Memory memory(budget);
    std::uint64_t metadata = checked_product(b, 512);
    for (auto &bag : d.bags) {
        std::sort(bag.begin(), bag.end());
        auto k = bag.size();
        metadata += checked_product(checked_product(k + 1, k + 1), 32);
        checked_power3(k);
    }
    memory.add(metadata);
    if (options.samples)
        memory.add(checked_product(std::min<std::size_t>(remaining_samples, b * 100), 32));
    auto tree = detail::tree_adj(d);
    std::vector<int> parent(b, -1), order{0};
    std::vector<std::vector<int>> children(b), out_positions(b), parent_positions(b);
    for (std::size_t i = 0; i < order.size(); ++i) {
        int t = order[i];
        for (int u : tree[t])
            if (u != parent[t]) {
                parent[u] = t;
                children[t].push_back(u);
                order.push_back(u);
            }
    }
    std::vector<Table> tables(b);
    result.bags.resize(b);
    for (std::size_t t = 0; t < b; ++t) {
        auto &bs = result.bags[t];
        bs.parent = parent[t];
        if (parent[t] >= 0) {
            auto &p = d.bags[parent[t]];
            auto &bag = d.bags[t];
            std::set_intersection(bag.begin(), bag.end(), p.begin(), p.end(),
                                  std::back_inserter(bs.separator));
            out_positions[t] = positions(bag, bs.separator);
            parent_positions[t] = positions(p, bs.separator);
        }
    }
    // Global adjacency lookup is a separate preparation cost, not part of the
    // single-exponential recurrence. Once prepared, each constraint is local.
    std::vector<std::vector<std::vector<int>>> local_edges(b);
    for (std::size_t t = 0; t < b; ++t) {
        budget.check();
        const auto &bag = d.bags[t];
        local_edges[t].resize(bag.size());
        for (std::size_t i = 0; i < bag.size(); ++i)
            for (std::size_t j = 0; j < i; ++j)
                if (std::binary_search(g.adj[bag[i]].begin(), g.adj[bag[i]].end(), bag[j]))
                    local_edges[t][i].push_back(static_cast<int>(j));
    }
    result.stats.preparation_ms += detail::elapsed(preparation_start);
    auto start = Clock::now();
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        budget.check();
        int t = *it;
        auto &bag = d.bags[t];
        auto &table = tables[t];
        auto &bs = result.bags[t];
        bool impossible = false;
        for (int u : children[t])
            if (tables[u].count == 0)
                impossible = true;
        if (impossible)
            continue;
        const auto domain = checked_power3(bs.separator.size());
        const bool dense = options.table_mode == TableMode::Dense ||
                           (options.table_mode == TableMode::Auto && domain <= 81 &&
                            domain <= budget.available() / 8);
        if (dense) {
            memory.add(checked_product(domain, 8));
            table.dense.assign(static_cast<std::size_t>(domain), absent);
        }
        detail::Memory temporary(budget);
        std::map<Code, Code> sparse;
        const auto &earlier = local_edges[t];
        std::vector<int> colors(bag.size());
        std::uint64_t ticks = 0;
        std::function<void(std::size_t, Code, Code)> enumerate = [&](std::size_t index, Code code,
                                                                     Code power) {
            if ((++ticks & 1023) == 0)
                budget.check();
            if (index < bag.size()) {
                for (int c = 0; c < 3; ++c) {
                    bool ok = true;
                    for (int j : earlier[index])
                        if (colors[j] == c) {
                            ok = false;
                            break;
                        }
                    if (ok) {
                        colors[index] = c;
                        enumerate(index + 1, code + static_cast<Code>(c) * power, power * 3);
                    }
                }
                return;
            }
            ++bs.candidates;
            ++result.stats.candidates;
            for (int u : children[t]) {
                ++result.stats.lookups;
                if (tables[u].get(project(colors, parent_positions[u])) == absent)
                    return;
            }
            ++bs.feasible;
            ++result.stats.feasible;
            ++result.stats.projections;
            Code key = project(colors, out_positions[t]);
            if (!table.dense.empty()) {
                auto &entry = table.dense[static_cast<std::size_t>(key)];
                if (entry == absent) {
                    entry = code;
                    ++table.count;
                }
            } else if (!sparse.contains(key)) {
                temporary.add(64);
                sparse.emplace(key, code);
            }
        };
        enumerate(0, 0, 1);
        if (table.dense.empty()) {
            table.count = sparse.size();
            // Conversion is chosen from actual occupancy, not a guessed density.
            if (options.table_mode == TableMode::Auto && domain <= sparse.size() * 2ULL &&
                domain <= budget.available() / 8) {
                memory.add(checked_product(domain, 8));
                table.dense.assign(static_cast<std::size_t>(domain), absent);
                for (auto [key, witness] : sparse)
                    table.dense[static_cast<std::size_t>(key)] = witness;
            } else {
                memory.add(checked_product(sparse.size(), sizeof(std::pair<Code, Code>)));
                table.sparse.assign(sparse.begin(), sparse.end());
            }
        }
        bs.projected = table.count;
        bs.dense = !table.dense.empty();
        if (options.samples && remaining_samples) {
            auto keep = [&](Code key, Code code) {
                if (bs.samples.size() < 100 && remaining_samples) {
                    bs.samples.push_back({key, code});
                    --remaining_samples;
                }
            };
            if (bs.dense) {
                for (std::size_t i = 0;
                     i < table.dense.size() && bs.samples.size() < 100 && remaining_samples; ++i) {
                    if ((i & 1023) == 0)
                        budget.check();
                    if (table.dense[i] != absent)
                        keep(i, table.dense[i]);
                }
            } else
                for (auto [key, code] : table.sparse) {
                    keep(key, code);
                    if (bs.samples.size() >= 100 || !remaining_samples)
                        break;
                }
        }
    }
    result.stats.dp_ms += detail::elapsed(start);
    if (tables[0].count == 0)
        return false;
    start = Clock::now();
    result.coloring.assign(g.names.size(), -1);
    std::vector<Code> chosen(b, absent);
    chosen[0] = tables[0].get(0);
    for (int t : order) {
        budget.check();
        Code code = chosen[t];
        std::vector<int> colors(d.bags[t].size());
        for (std::size_t i = 0; i < colors.size(); ++i) {
            colors[i] = static_cast<int>(code % 3);
            code /= 3;
            int v = d.bags[t][i];
            if (result.coloring[v] >= 0 && result.coloring[v] != colors[i])
                throw std::logic_error("inconsistent reconstruction");
            result.coloring[v] = colors[i];
        }
        for (int u : children[t]) {
            chosen[u] = tables[u].get(project(colors, parent_positions[u]));
            if (chosen[u] == absent)
                throw std::logic_error("missing reconstruction witness");
        }
    }
    result.stats.reconstruction_ms += detail::elapsed(start);
    return true;
}
void add_stats(Stats &a, const Stats &b) {
    a.preparation_ms += b.preparation_ms;
    a.decomposition_ms += b.decomposition_ms;
    a.validation_ms += b.validation_ms;
    a.dp_ms += b.dp_ms;
    a.reconstruction_ms += b.reconstruction_ms;
    a.candidates += b.candidates;
    a.feasible += b.feasible;
    a.projections += b.projections;
    a.lookups += b.lookups;
}
// Combines construction records of residual components. The joint width is the
// maximum component width; it is optimal when every component width is optimal.
void merge_td(TdInfo &total, const TdInfo &part, int part_width, int current_width) {
    const bool first = total.components == 0;
    if (first || part_width > current_width) {
        total.method = part.method;
        total.candidates = part.candidates;
    }
    if (first || part.width_status != "optimal") {
        total.width_status = part.width_status;
        total.proof = part.proof;
    } else if (total.width_status == "optimal" && part.proof == "exact_search")
        total.proof = "exact_search";
    total.components += part.components;
    total.exact_components += part.exact_components;
    total.raw_bags += part.raw_bags;
    total.merged_bags += part.merged_bags;
    total.raw_width = std::max(total.raw_width, part.raw_width);
    total.lower_bound = std::max(total.lower_bound, part.lower_bound);
}
std::uint64_t decomposition_bytes(const Decomposition &d) {
    std::uint64_t size = checked_product(d.bags.size(), 192);
    for (auto &bag : d.bags)
        size += checked_product(bag.size(), 16);
    return size;
}
} // namespace
Result solve(const Graph &g, const Options &options, const Decomposition *supplied) {
    Result result;
    result.supplied_decomposition = supplied != nullptr;
    Budget budget(options);
    std::size_t remaining_samples = 10000;
    try {
        detail::Memory persistent(budget);
        persistent.add(checked_product(g.names.size(), 192) + checked_product(g.edges.size(), 48));
        budget.check();
        if (supplied) {
            persistent.add(decomposition_bytes(*supplied));
            result.decomposition = *supplied;
            auto start = Clock::now();
            validate_decomposition(g, result.decomposition, &budget);
            auto &td = result.td;
            td.method = "supplied";
            td.raw_width = result.decomposition.width();
            td.raw_bags = result.decomposition.bags.size();
            if (options.simplify_supplied) {
                td.merged_bags = contract_contained_bags(result.decomposition);
                validate_decomposition(g, result.decomposition, &budget);
            }
            td.lower_bound = treewidth_lower_bound(g, budget);
            td.components = 1;
            if (td.lower_bound == result.decomposition.width()) {
                td.width_status = "optimal";
                td.proof = "lower_bound_match";
            } else
                td.width_status = "upper_bound";
            result.stats.validation_ms = detail::elapsed(start);
        }
        auto generate = [&] {
            auto start = Clock::now();
            result.decomposition = build_decomposition(g, options, budget, result.td);
            persistent.add(decomposition_bytes(result.decomposition));
            result.stats.decomposition_ms = detail::elapsed(start);
            start = Clock::now();
            validate_decomposition(g, result.decomposition, &budget);
            result.stats.validation_ms = detail::elapsed(start);
        };
        if (g.self_loop) {
            if (!supplied && !options.preprocess)
                generate();
            result.status = Status::NotColorable;
            result.message = "A self-loop forbids a proper coloring.";
        } else if (supplied || !options.preprocess) {
            if (!supplied)
                generate();
            result.status =
                run_dp(g, result.decomposition, options, budget, result, remaining_samples)
                    ? Status::Colorable
                    : Status::NotColorable;
        } else {
            result.used_preprocessing = true;
            std::vector<bool> alive(g.names.size(), true);
            std::vector<int> degree(g.names.size()), peeled;
            std::queue<int> queue;
            for (std::size_t v = 0; v < g.names.size(); ++v) {
                degree[v] = static_cast<int>(g.adj[v].size());
                if (degree[v] <= 2)
                    queue.push(static_cast<int>(v));
            }
            while (!queue.empty()) {
                budget.check();
                int v = queue.front();
                queue.pop();
                if (!alive[v])
                    continue;
                alive[v] = false;
                peeled.push_back(v);
                for (int u : g.adj[v])
                    if (alive[u] && --degree[u] == 2)
                        queue.push(u);
            }
            result.stats.peeled_vertices = peeled.size();
            result.coloring.assign(g.names.size(), -1);
            result.status = Status::Colorable;
            std::vector<bool> seen(g.names.size(), false);
            std::vector<int> mapped(g.names.size(), -1);
            for (std::size_t first = 0; first < g.names.size(); ++first)
                if (alive[first] && !seen[first]) {
                    budget.check();
                    std::vector<int> component{static_cast<int>(first)};
                    seen[first] = true;
                    for (std::size_t i = 0; i < component.size(); ++i)
                        for (int u : g.adj[component[i]])
                            if (alive[u] && !seen[u]) {
                                seen[u] = true;
                                component.push_back(u);
                            }
                    std::sort(component.begin(), component.end());
                    Graph sub;
                    for (int v : component) {
                        mapped[v] = static_cast<int>(sub.names.size());
                        sub.names.push_back(g.names[v]);
                        sub.labels.emplace_back();
                    }
                    sub.adj.resize(component.size());
                    for (int v : component)
                        for (int u : g.adj[v])
                            if (alive[u]) {
                                sub.adj[mapped[v]].push_back(mapped[u]);
                                if (v < u) {
                                    sub.edges.push_back({mapped[v], mapped[u]});
                                    sub.edge_labels.emplace_back();
                                }
                            }
                    Result part;
                    auto start = Clock::now();
                    auto d = build_decomposition(sub, options, budget, part.td);
                    part.stats.decomposition_ms = detail::elapsed(start);
                    merge_td(result.td, part.td, d.width(), result.decomposition.width());
                    start = Clock::now();
                    validate_decomposition(sub, d, &budget);
                    part.stats.validation_ms = detail::elapsed(start);
                    detail::Memory component_memory(budget);
                    component_memory.add(decomposition_bytes(d) +
                                         checked_product(component.size(), 192) +
                                         checked_product(sub.edges.size(), 48));
                    bool ok = run_dp(sub, d, options, budget, part, remaining_samples);
                    add_stats(result.stats, part.stats);
                    ++result.stats.solved_components;
                    int offset = static_cast<int>(result.decomposition.bags.size());
                    for (std::size_t t = 0; t < d.bags.size(); ++t) {
                        persistent.add(
                            192 + checked_product(d.bags[t].size(), 16) +
                            checked_product(part.bags[t].samples.capacity(), sizeof(Sample)));
                        for (int &v : d.bags[t])
                            v = component[v];
                        for (int &v : part.bags[t].separator)
                            v = component[v];
                        // Each component root stays a DP root; the exported tree connects them by
                        // empty intersections.
                        if (part.bags[t].parent >= 0)
                            part.bags[t].parent += offset;
                        result.decomposition.names.push_back(
                            "B" + std::to_string(offset + static_cast<int>(t)));
                        result.decomposition.bags.push_back(std::move(d.bags[t]));
                        result.bags.push_back(std::move(part.bags[t]));
                    }
                    for (auto [u, v] : d.edges)
                        result.decomposition.edges.push_back({u + offset, v + offset});
                    if (offset && !d.bags.empty())
                        result.decomposition.edges.push_back({0, offset});
                    if (!ok) {
                        result.status = Status::NotColorable;
                        break;
                    }
                    for (std::size_t i = 0; i < component.size(); ++i)
                        result.coloring[component[i]] = part.coloring[i];
                }
            if (result.td.components && result.td.width_status != "optimal" &&
                result.td.lower_bound == result.decomposition.width()) {
                result.td.width_status = "optimal";
                result.td.proof = "lower_bound_match";
            }
            if (result.status == Status::Colorable) {
                auto start = Clock::now();
                for (auto it = peeled.rbegin(); it != peeled.rend(); ++it) {
                    budget.check();
                    bool used[3] = {false, false, false};
                    for (int u : g.adj[*it])
                        if (result.coloring[u] >= 0)
                            used[result.coloring[u]] = true;
                    int c = 0;
                    while (c < 3 && used[c])
                        ++c;
                    if (c == 3)
                        throw std::logic_error("peeling reconstruction failed");
                    result.coloring[*it] = c;
                }
                result.stats.reconstruction_ms += detail::elapsed(start);
            }
        }
        if (result.status == Status::Colorable && !verify_coloring(g, result.coloring))
            throw std::logic_error("final coloring verification failed");
        if (result.status != Status::Colorable)
            result.coloring.clear();
    } catch (const InputError &e) {
        result.status = Status::InvalidInput;
        result.message = e.what();
        result.coloring.clear();
    } catch (const ResourceError &e) {
        result.status = Status::ResourceLimit;
        result.message = e.what();
        result.coloring.clear();
    } catch (const std::bad_alloc &) {
        result.status = Status::ResourceLimit;
        result.message = "allocation failed; colorability is unknown";
        result.coloring.clear();
    } catch (const std::length_error &) {
        result.status = Status::ResourceLimit;
        result.message = "container size limit exceeded";
        result.coloring.clear();
    }
    result.stats.peak_accounted_bytes = budget.peak();
    return result;
}
} // namespace td3
