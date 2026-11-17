#include "detail.hpp"
namespace td3 {
Graph read_graph(std::istream &in) {
    Graph g;
    std::unordered_map<std::string, int> ids;
    std::map<std::pair<int, int>, std::string> edges;
    std::map<int, std::string> labels;
    detail::rows(in, [&](const detail::Row &r) {
        int u = detail::intern(r.left, ids, g.names);
        if (r.right.empty()) {
            if (!r.label.empty()) {
                auto [it, added] = labels.emplace(u, r.label);
                if (!added && it->second != r.label)
                    throw InputError("line " + std::to_string(r.line) +
                                     ": conflicting vertex labels");
            }
        } else {
            int v = detail::intern(r.right, ids, g.names);
            if (u > v)
                std::swap(u, v);
            auto [it, added] = edges.emplace(std::pair{u, v}, r.label);
            if (!added && it->second != r.label)
                throw InputError("line " + std::to_string(r.line) +
                                 ": conflicting duplicate edge labels");
        }
    });
    if (!labels.empty() && labels.size() != g.names.size())
        throw InputError("vertex labels must be supplied for all vertices or none");
    g.labels.resize(g.names.size());
    g.adj.resize(g.names.size());
    for (auto &[v, label] : labels)
        g.labels[v] = label;
    bool labelled = false, unlabelled = false;
    for (auto &[e, label] : edges) {
        labelled |= !label.empty();
        unlabelled |= label.empty();
        auto [u, v] = e;
        g.edges.push_back(e);
        g.edge_labels.push_back(label);
        if (u == v)
            g.self_loop = true;
        g.adj[u].push_back(v);
        if (u != v)
            g.adj[v].push_back(u);
    }
    if (labelled && unlabelled)
        throw InputError("edge labels must be supplied for all edges or none");
    for (auto &a : g.adj)
        std::sort(a.begin(), a.end());
    return g;
}
void write_graph(const Graph &g, std::ostream &out) {
    for (std::size_t i = 0; i < g.names.size(); ++i) {
        out << g.names[i] << ',';
        if (i < g.labels.size() && !g.labels[i].empty())
            out << ',' << g.labels[i];
        out << '\n';
    }
    for (std::size_t i = 0; i < g.edges.size(); ++i) {
        auto [u, v] = g.edges[i];
        out << g.names[u] << ',' << g.names[v];
        if (i < g.edge_labels.size() && !g.edge_labels[i].empty())
            out << ',' << g.edge_labels[i];
        out << '\n';
    }
}
Graph read_legacy(std::istream &in) {
    std::string line;
    long long n = -1;
    char extra;
    if (!std::getline(in, line))
        throw InputError("missing vertex count");
    std::istringstream count(line);
    if (!(count >> n) || n < 0 || n > std::numeric_limits<int>::max() || (count >> extra))
        throw InputError("invalid vertex count");
    if (!std::getline(in, line))
        throw InputError("missing weight line");
    std::istringstream weights(line);
    long long weight;
    for (long long i = 0; i < n; ++i)
        if (!(weights >> weight))
            throw InputError("weight count must equal vertex count");
    if (weights >> extra)
        throw InputError("extra weight fields");
    Graph g;
    g.names.reserve(static_cast<std::size_t>(n));
    g.labels.resize(static_cast<std::size_t>(n));
    g.adj.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        g.names.push_back(std::to_string(i + 1));
    std::map<std::pair<int, int>, bool> seen;
    std::size_t lineno = 2;
    while (std::getline(in, line)) {
        ++lineno;
        if (line.find_first_not_of(" \t\r") == std::string::npos)
            continue;
        long long u, v;
        std::istringstream row(line);
        if (!(row >> u >> v) || (row >> extra) || u < 1 || u > n || v < 1 || v > n)
            throw InputError("line " + std::to_string(lineno) + ": invalid legacy edge");
        if (u > v)
            std::swap(u, v);
        seen[{static_cast<int>(u - 1), static_cast<int>(v - 1)}] = true;
    }
    for (auto &[e, _] : seen) {
        auto [u, v] = e;
        g.edges.push_back(e);
        g.edge_labels.emplace_back();
        g.adj[u].push_back(v);
        if (u != v)
            g.adj[v].push_back(u);
        else
            g.self_loop = true;
    }
    for (auto &a : g.adj)
        std::sort(a.begin(), a.end());
    return g;
}
} // namespace td3
