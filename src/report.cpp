#include "report_template.hpp"
#include "td3/td3.hpp"
#include <iomanip>
#include <sstream>
namespace td3 {
std::string status_name(Status s) {
    switch (s) {
    case Status::Colorable:
        return "COLORABLE";
    case Status::NotColorable:
        return "NOT_COLORABLE";
    case Status::InvalidInput:
        return "INVALID_INPUT";
    case Status::ResourceLimit:
        return "RESOURCE_LIMIT";
    }
    return "INVALID_INPUT";
}
namespace {
std::string quote(const std::string &s) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : s) {
        switch (c) {
        case '"':
            out << "\\\"";
            break;
        case '\\':
            out << "\\\\";
            break;
        case '<':
            out << "\\u003c";
            break;
        case '>':
            out << "\\u003e";
            break;
        case '&':
            out << "\\u0026";
            break;
        default:
            if (c < 32)
                out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c)
                    << std::dec;
            else
                out << static_cast<char>(c);
        }
    }
    out << '"';
    return out.str();
}
template <class T> void numbers(std::ostream &out, const std::vector<T> &values) {
    out << '[';
    bool first = true;
    for (auto v : values) {
        if (!first)
            out << ',';
        first = false;
        out << v;
    }
    out << ']';
}
void colors(std::ostream &out, Code code, std::size_t n) {
    out << '[';
    for (std::size_t i = 0; i < n; ++i) {
        if (i)
            out << ',';
        out << code % 3;
        code /= 3;
    }
    out << ']';
}
} // namespace
std::string result_json(const Graph &g, const Result &r) {
    std::ostringstream out;
    out << std::setprecision(12) << "{\"version\":1,\"status\":" << quote(status_name(r.status))
        << ",\"message\":" << quote(r.message)
        << ",\"preprocessed\":" << (r.used_preprocessing ? "true" : "false")
        << ",\"supplied_decomposition\":" << (r.supplied_decomposition ? "true" : "false")
        << ",\"graph\":{\"vertices\":[";
    for (std::size_t i = 0; i < g.names.size(); ++i) {
        if (i)
            out << ',';
        out << "{\"id\":" << i << ",\"name\":" << quote(g.names[i])
            << ",\"label\":" << quote(i < g.labels.size() ? g.labels[i] : "") << '}';
    }
    out << "],\"edges\":[";
    for (std::size_t i = 0; i < g.edges.size(); ++i) {
        if (i)
            out << ',';
        auto [u, v] = g.edges[i];
        out << '[' << u << ',' << v << ']';
    }
    out << "],\"edge_labels\":[";
    for (std::size_t i = 0; i < g.edges.size(); ++i) {
        if (i)
            out << ',';
        out << quote(i < g.edge_labels.size() ? g.edge_labels[i] : "");
    }
    out << "]},\"coloring\":";
    if (r.status == Status::Colorable)
        numbers(out, r.coloring);
    else
        out << "null";
    const auto &s = r.stats;
    out << ",\"stats\":{\"decomposition_ms\":" << s.decomposition_ms
        << ",\"preparation_ms\":" << s.preparation_ms << ",\"validation_ms\":" << s.validation_ms
        << ",\"dp_ms\":" << s.dp_ms << ",\"reconstruction_ms\":" << s.reconstruction_ms
        << ",\"candidates\":" << s.candidates << ",\"feasible\":" << s.feasible
        << ",\"projections\":" << s.projections << ",\"lookups\":" << s.lookups
        << ",\"peak_accounted_bytes\":" << s.peak_accounted_bytes
        << ",\"peeled_vertices\":" << s.peeled_vertices
        << ",\"solved_components\":" << s.solved_components
        << "},\"decomposition\":{\"width\":" << r.decomposition.width();
    const auto &td = r.td;
    out << ",\"method\":" << quote(td.method) << ",\"width_status\":" << quote(td.width_status)
        << ",\"proof\":" << quote(td.proof) << ",\"lower_bound\":" << td.lower_bound
        << ",\"raw_width\":" << td.raw_width << ",\"raw_bags\":" << td.raw_bags
        << ",\"merged_bags\":" << td.merged_bags << ",\"components\":" << td.components
        << ",\"exact_components\":" << td.exact_components << ",\"candidates\":[";
    for (std::size_t i = 0; i < td.candidates.size(); ++i) {
        const auto &c = td.candidates[i];
        if (i)
            out << ',';
        out << "{\"name\":" << quote(c.name) << ",\"width\":" << c.width
            << ",\"raw_bags\":" << c.raw_bags << ",\"bags\":" << c.bags << ",\"cost\":" << c.cost
            << ",\"ms\":" << c.ms << ",\"selected\":" << (c.selected ? "true" : "false") << '}';
    }
    out << "],\"bags\":[";
    for (std::size_t t = 0; t < r.decomposition.bags.size(); ++t) {
        if (t)
            out << ',';
        BagStats empty;
        const auto &bs = t < r.bags.size() ? r.bags[t] : empty;
        out << "{\"name\":" << quote(r.decomposition.names[t]) << ",\"vertices\":";
        numbers(out, r.decomposition.bags[t]);
        out << ",\"parent\":" << bs.parent << ",\"separator\":";
        numbers(out, bs.separator);
        out << ",\"candidates\":" << bs.candidates << ",\"feasible\":" << bs.feasible
            << ",\"projected\":" << bs.projected
            << ",\"storage\":" << quote(bs.dense ? "dense" : "sparse") << ",\"samples\":[";
        for (std::size_t i = 0; i < bs.samples.size(); ++i) {
            if (i)
                out << ',';
            const auto &sample = bs.samples[i];
            out << "{\"separator_colors\":";
            colors(out, sample.separator, bs.separator.size());
            out << ",\"bag_colors\":";
            colors(out, sample.witness, r.decomposition.bags[t].size());
            out << '}';
        }
        out << "]}";
    }
    out << "],\"edges\":[";
    for (std::size_t i = 0; i < r.decomposition.edges.size(); ++i) {
        if (i)
            out << ',';
        auto [u, v] = r.decomposition.edges[i];
        out << '[' << u << ',' << v << ']';
    }
    out << "]}}";
    return out.str();
}
void write_report(const Graph &g, const Result &r, std::ostream &out) {
    std::string_view html(report_template), marker = "__TD3_DATA__";
    auto at = html.find(marker);
    out << html.substr(0, at) << result_json(g, r) << html.substr(at + marker.size());
}
} // namespace td3
