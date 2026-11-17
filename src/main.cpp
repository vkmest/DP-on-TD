#include "td3/td3.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
using namespace td3;
namespace fs = std::filesystem;
namespace {
constexpr auto help = R"(td3color - exact 3-colorability by tree-decomposition projection DP

td3color solve GRAPH.csv [--td TD.csv] [--report REPORT.html] [--json RESULT.json]
    [--write-td TD.csv] [--memory-mib 1024] [--time-limit SECONDS]
    [--preprocess | --no-preprocess] [--table auto|dense|sparse]
    [--td-method auto|minfill|best|exact] [--exact-limit 20] [--seeds 4]
    [--no-simplify] [--simplify-td]
td3color validate-td GRAPH.csv TD.csv
td3color convert-legacy OLD.txt NEW.csv

Reports and --write-td default to no preprocessing for teaching and interoperability.
--td always preserves the supplied decomposition. --write-td requires no preprocessing.
Decomposition methods (when no --td is given):
  minfill  one deterministic minimum-fill elimination order
  best     min-fill, min-degree, max-cardinality search and seeded tie-breaks;
           keep the smallest width, then the fewest DP states
  auto     best, then exact minimum-width search on components of at most
           --exact-limit vertices (default; 0 disables, at most 26)
  exact    like auto, but fail with RESOURCE_LIMIT if a component is too large
Generated decompositions are simplified (contained bags merged, rooted at a centre)
unless --no-simplify is given. --simplify-td also merges contained bags of a --td.
Exit codes: 0 colorable/valid, 1 not colorable, 2 invalid input or I/O, 3 resource limit.
The memory limit accounts major solver buffers, not total process RSS or input parsing.
)";
std::ifstream input(const std::string &path) {
    std::ifstream in(path);
    if (!in)
        throw InputError("cannot open input: " + path);
    return in;
}
template <class F> void output(const std::string &path, F write) {
    if (path.empty())
        return;
    fs::path p(path);
    if (p.has_parent_path())
        fs::create_directories(p.parent_path());
    std::ofstream out(p, std::ios::binary);
    if (!out)
        throw InputError("cannot create output: " + path);
    write(out);
    out.close();
    if (!out)
        throw InputError("cannot finish output: " + path);
}
std::string canonical(const std::string &p) {
    auto s = fs::weakly_canonical(fs::absolute(p)).generic_string();
#ifdef _WIN32
    for (char &c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
#endif
    return s;
}
void distinct(const std::vector<std::string> &inputs, const std::vector<std::string> &outputs) {
    std::set<std::string> paths;
    std::vector<std::string> existing;
    for (auto &p : inputs)
        if (!p.empty()) {
            paths.insert(canonical(p));
            if (fs::exists(p))
                existing.push_back(p);
        }
    for (auto &p : outputs) {
        if (p.empty())
            continue;
        if (!paths.insert(canonical(p)).second)
            throw InputError("output paths must be distinct and must not overwrite inputs");
        if (fs::exists(p)) {
            for (const auto &prior : existing)
                if (fs::equivalent(p, prior))
                    throw InputError("output aliases an input or another output (hard link)");
            existing.push_back(p);
        }
    }
}
std::uint64_t integer(const std::string &s) {
    if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
        throw InputError("expected nonnegative integer: " + s);
    try {
        return std::stoull(s);
    } catch (...) {
        throw InputError("integer out of range");
    }
}
} // namespace
int main(int argc, char **argv) {
    Graph g;
    Result result;
    std::string json, report, td_path, td_output;
    try {
        if (argc == 1 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
            std::cout << help;
            return 0;
        }
        std::string command = argv[1];
        if (command == "convert-legacy") {
            if (argc != 4)
                throw InputError("convert-legacy requires input and output paths");
            distinct({argv[2]}, {argv[3]});
            auto in = input(argv[2]);
            auto old = read_legacy(in);
            output(argv[3], [&](auto &out) { write_graph(old, out); });
            std::cout << "Converted " << old.names.size() << " vertices and " << old.edges.size()
                      << " edges.\n";
            return 0;
        }
        if (command == "validate-td") {
            if (argc != 4)
                throw InputError("validate-td requires graph and decomposition paths");
            auto in = input(argv[2]);
            g = read_graph(in);
            auto tin = input(argv[3]);
            auto td = read_decomposition(tin, g);
            validate_decomposition(g, td);
            std::cout << "VALID width=" << td.width() << " bags=" << td.bags.size() << '\n';
            return 0;
        }
        if (command != "solve" || argc < 3)
            throw InputError("expected solve, validate-td or convert-legacy; use --help");
        Options options;
        bool explicit_preprocess = false;
        std::set<std::string> flags;
        for (int i = 3; i < argc; ++i) {
            std::string arg = argv[i];
            if (!flags.insert(arg).second)
                throw InputError("duplicate option: " + arg);
            auto value = [&]() {
                if (++i >= argc)
                    throw InputError("missing value for " + arg);
                return std::string(argv[i]);
            };
            if (arg == "--td")
                td_path = value();
            else if (arg == "--json")
                json = value();
            else if (arg == "--report")
                report = value();
            else if (arg == "--write-td")
                td_output = value();
            else if (arg == "--memory-mib")
                options.memory_bytes = checked_product(integer(value()), 1024 * 1024);
            else if (arg == "--time-limit") {
                auto s = value();
                std::size_t used = 0;
                try {
                    options.time_limit_seconds = std::stod(s, &used);
                } catch (...) {
                    throw InputError("invalid time limit");
                }
                if (used != s.size() || !std::isfinite(options.time_limit_seconds) ||
                    options.time_limit_seconds < 0)
                    throw InputError("time limit must be finite and nonnegative");
            } else if (arg == "--preprocess" || arg == "--no-preprocess") {
                if (explicit_preprocess)
                    throw InputError("conflicting preprocessing options");
                explicit_preprocess = true;
                options.preprocess = arg == "--preprocess";
            } else if (arg == "--table") {
                auto s = value();
                if (s == "auto")
                    options.table_mode = TableMode::Auto;
                else if (s == "dense")
                    options.table_mode = TableMode::Dense;
                else if (s == "sparse")
                    options.table_mode = TableMode::Sparse;
                else
                    throw InputError("table must be auto, dense or sparse");
            } else if (arg == "--td-method") {
                auto s = value();
                if (s == "auto")
                    options.td_method = TdMethod::Auto;
                else if (s == "minfill")
                    options.td_method = TdMethod::MinFill;
                else if (s == "best")
                    options.td_method = TdMethod::Best;
                else if (s == "exact")
                    options.td_method = TdMethod::Exact;
                else
                    throw InputError("td-method must be auto, minfill, best or exact");
            } else if (arg == "--exact-limit") {
                auto limit = integer(value());
                if (limit > 26)
                    throw InputError("exact-limit must be at most 26");
                options.exact_limit = static_cast<int>(limit);
            } else if (arg == "--seeds") {
                auto seeds = integer(value());
                if (seeds > 1000)
                    throw InputError("seeds must be at most 1000");
                options.heuristic_seeds = static_cast<int>(seeds);
            } else if (arg == "--no-simplify")
                options.simplify = false;
            else if (arg == "--simplify-td")
                options.simplify_supplied = true;
            else
                throw InputError("unknown option: " + arg);
        }
        if (!explicit_preprocess && (!report.empty() || !td_output.empty()))
            options.preprocess = false;
        if (!td_path.empty())
            options.preprocess = false;
        if (!td_output.empty() && options.preprocess)
            throw InputError("--write-td requires --no-preprocess");
        distinct({argv[2], td_path}, {json, report, td_output});
        auto in = input(argv[2]);
        g = read_graph(in);
        Decomposition td;
        if (!td_path.empty()) {
            auto tin = input(td_path);
            td = read_decomposition(tin, g);
        }
        options.samples = !report.empty();
        result = solve(g, options, td_path.empty() ? nullptr : &td);
        output(json, [&](auto &out) { out << result_json(g, result) << '\n'; });
        output(report, [&](auto &out) { write_report(g, result, out); });
        if (!td_output.empty() &&
            (result.status == Status::Colorable || result.status == Status::NotColorable)) {
            output(td_output,
                   [&](auto &out) { write_decomposition(result.decomposition, g, out); });
        }
        std::cout << status_name(result.status) << '\n';
        if (result.status == Status::Colorable)
            for (int c = 0; c < 3; ++c) {
                bool first = true;
                for (std::size_t v = 0; v < g.names.size(); ++v)
                    if (result.coloring[v] == c) {
                        if (!first)
                            std::cout << ' ';
                        std::cout << g.names[v];
                        first = false;
                    }
                std::cout << '\n';
            }
        std::cerr << std::fixed << std::setprecision(3) << "width=" << result.decomposition.width()
                  << " width_status=" << result.td.width_status
                  << " lower_bound=" << result.td.lower_bound << " td_method=" << result.td.method
                  << " bags=" << result.decomposition.bags.size() << " dp_ms=" << result.stats.dp_ms
                  << " decomposition_ms=" << result.stats.decomposition_ms
                  << " preparation_ms=" << result.stats.preparation_ms
                  << " lookups=" << result.stats.lookups
                  << " accounted_peak_bytes=" << result.stats.peak_accounted_bytes << '\n';
        if (!result.message.empty())
            std::cerr << result.message << '\n';
        return static_cast<int>(result.status);
    } catch (const ResourceError &e) {
        std::cerr << "RESOURCE_LIMIT: " << e.what() << '\n';
        return 3;
    } catch (const std::bad_alloc &) {
        std::cerr << "RESOURCE_LIMIT: allocation failed\n";
        return 3;
    } catch (const std::exception &e) {
        std::cerr << "INVALID_INPUT: " << e.what() << '\n';
        return 2;
    }
}
