#pragma once
#include <chrono>
#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace td3 {
using Code = std::uint64_t;
using Clock = std::chrono::steady_clock;
struct InputError : std::runtime_error {
    using runtime_error::runtime_error;
};
struct ResourceError : std::runtime_error {
    using runtime_error::runtime_error;
};
struct Graph {
    std::vector<std::string> names, labels;
    std::vector<std::vector<int>> adj;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::string> edge_labels;
    bool self_loop = false;
};
struct Decomposition {
    std::vector<std::string> names;
    std::vector<std::vector<int>> bags;
    std::vector<std::pair<int, int>> edges;
    int width() const;
};
enum class TableMode { Auto, Dense, Sparse };
// How a decomposition is constructed when none is supplied.
//   MinFill: the single deterministic minimum-fill elimination order (original behaviour).
//   Best:    several elimination heuristics; keep the lowest width, then the lowest DP cost.
//   Auto:    Best, then exact minimum-width search for components within exact_limit.
//   Exact:   like Auto, but every component must be solved exactly (else RESOURCE_LIMIT).
enum class TdMethod { Auto, MinFill, Best, Exact };
struct Options {
    std::uint64_t memory_bytes = 1024ULL * 1024 * 1024;
    double time_limit_seconds = 0;
    bool preprocess = true;
    bool samples = false;
    TableMode table_mode = TableMode::Auto;
    TdMethod td_method = TdMethod::Auto;
    int exact_limit = 20;          // max vertices per connected component for exact search
    int heuristic_seeds = 4;       // extra randomized tie-break runs per greedy heuristic
    bool simplify = true;          // contract contained bags and re-root generated TDs
    bool simplify_supplied = false; // also contract a supplied TD (names of survivors kept)
};
// One decomposition considered during construction.
struct TdCandidate {
    std::string name;
    int width = -1;
    std::size_t raw_bags = 0, bags = 0;
    double cost = 0; // sum over bags of 3^|bag| (upper bound on enumerated DP states)
    double ms = 0;
    bool selected = false;
};
// How the decomposition was obtained, and what is known about its width.
struct TdInfo {
    std::string method = "none";             // e.g. "min-fill", "exact", "supplied"
    std::string width_status = "none";       // "optimal" or "upper_bound" or "none"
    std::string proof;                       // "exact_search", "lower_bound_match" or ""
    int lower_bound = -1;                    // proven lower bound on treewidth
    int raw_width = -1;                      // width before simplification
    std::size_t raw_bags = 0;                // bag count before simplification
    std::size_t merged_bags = 0;             // bags removed by contraction
    std::size_t exact_components = 0, components = 0;
    std::vector<TdCandidate> candidates;
};
struct Stats {
    double preparation_ms = 0;
    double decomposition_ms = 0, validation_ms = 0, dp_ms = 0, reconstruction_ms = 0;
    std::uint64_t candidates = 0, feasible = 0, projections = 0, lookups = 0,
                  peak_accounted_bytes = 0;
    std::size_t peeled_vertices = 0, solved_components = 0;
};
struct Sample {
    Code separator = 0, witness = 0;
};
struct BagStats {
    int parent = -1;
    std::vector<int> separator;
    std::uint64_t candidates = 0, feasible = 0, projected = 0;
    bool dense = false;
    std::vector<Sample> samples;
};
enum class Status { Colorable, NotColorable, InvalidInput, ResourceLimit };
struct Result {
    Status status = Status::InvalidInput;
    std::string message;
    std::vector<int> coloring;
    Decomposition decomposition;
    TdInfo td;
    std::vector<BagStats> bags;
    Stats stats;
    bool supplied_decomposition = false, used_preprocessing = false;
};
class Budget {
    Options options_;
    Clock::time_point start_ = Clock::now();
    std::uint64_t live_ = 0, peak_ = 0;

  public:
    explicit Budget(const Options &options) : options_(options) {}
    void check() const;
    void acquire(std::uint64_t bytes);
    void release(std::uint64_t bytes);
    std::uint64_t peak() const {
        return peak_;
    }
    std::uint64_t available() const {
        return options_.memory_bytes - live_;
    }
};
struct Reservation {
    Budget *budget;
    std::uint64_t bytes;
    Reservation(Budget &b, std::uint64_t n) : budget(&b), bytes(n) {
        b.acquire(n);
    }
    ~Reservation() {
        budget->release(bytes);
    }
    Reservation(const Reservation &) = delete;
    Reservation &operator=(const Reservation &) = delete;
};
std::uint64_t checked_power3(std::size_t exponent);
std::uint64_t checked_product(std::uint64_t a, std::uint64_t b);
Graph read_graph(std::istream &input);
void write_graph(const Graph &graph, std::ostream &output);
Graph read_legacy(std::istream &input);
Decomposition read_decomposition(std::istream &input, const Graph &graph);
void write_decomposition(const Decomposition &, const Graph &, std::ostream &);
void validate_decomposition(const Graph &, const Decomposition &, Budget *budget = nullptr);
Decomposition min_fill(const Graph &, Budget &);
// Elimination orders: order[i] is the i-th eliminated vertex.
enum class Greedy { MinFill, MinDegree };
std::vector<int> greedy_order(const Graph &, Budget &, Greedy rule,
                              const std::vector<int> *tie_rank = nullptr);
std::vector<int> mcs_order(const Graph &, Budget &);
// Exact minimum-width elimination order of one graph by O*(2^n) subset DP.
// Searches only widths below `cap`; returns an empty vector if none exists.
std::vector<int> exact_order(const Graph &, Budget &, int cap);
// Standard elimination-order decomposition: bag(v) = v plus its later filled neighbours.
Decomposition from_elimination_order(const Graph &, const std::vector<int> &order, Budget &);
// Contract every tree edge whose one bag contains the other. Width never increases.
// Returns the number of removed bags. Names of surviving bags are kept.
std::size_t contract_contained_bags(Decomposition &);
// Re-link independent pieces under the largest piece's centre, root at that centre,
// and renumber bags B0 (root), B1, ... in breadth-first order.
void arrange_tree(Decomposition &);
// Minor-min-width (MMD+) lower bound on treewidth; -1 for the empty graph.
int treewidth_lower_bound(const Graph &, Budget &);
double dp_cost(const Decomposition &);
// Builds (and validates) a decomposition for a graph that has no self-loop.
Decomposition build_decomposition(const Graph &, const Options &, Budget &, TdInfo &);
Result solve(const Graph &, const Options &options = {}, const Decomposition *supplied = nullptr);
bool verify_coloring(const Graph &, const std::vector<int> &);
std::string status_name(Status);
std::string result_json(const Graph &, const Result &);
void write_report(const Graph &, const Result &, std::ostream &);
} // namespace td3
