#pragma once
#include "td3/td3.hpp"
#include <algorithm>
#include <cctype>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <unordered_map>
namespace td3::detail {
struct Row {
    std::string left, right, label;
    bool third = false;
    std::size_t line = 0;
};
inline void name_ok(const std::string &s, std::size_t line) {
    if (s.empty())
        throw InputError("line " + std::to_string(line) + ": empty vertex name");
    for (unsigned char c : s)
        if (std::isspace(c) || c == ',' || c == '>' || c == ';' || c == '"' || c < 32)
            throw InputError("line " + std::to_string(line) +
                             ": vertex name contains whitespace or a reserved delimiter");
}
template <class F> void rows(std::istream &in, F consume) {
    std::string s;
    std::size_t line = 0;
    while (std::getline(in, s)) {
        ++line;
        if (!s.empty() && s.back() == '\r')
            s.pop_back();
        if (line == 1 && s.starts_with("\xef\xbb\xbf"))
            s.erase(0, 3);
        if (s.empty())
            continue;
        auto first = s.find_first_of(",>");
        if (first == std::string::npos)
            throw InputError("line " + std::to_string(line) + ": expected comma");
        if (s[first] == '>')
            throw InputError("line " + std::to_string(line) + ": directed graphs are unsupported");
        auto second = s.find(',', first + 1);
        Row r{s.substr(0, first),
              s.substr(first + 1, second == std::string::npos ? second : second - first - 1), "",
              second != std::string::npos, line};
        if (r.third)
            r.label = s.substr(second + 1);
        name_ok(r.left, line);
        if (!r.right.empty())
            name_ok(r.right, line);
        consume(r);
    }
    if (in.bad())
        throw InputError("failed while reading input");
}
inline int intern(const std::string &name, std::unordered_map<std::string, int> &ids,
                  std::vector<std::string> &names) {
    auto it = ids.find(name);
    if (it != ids.end())
        return it->second;
    if (names.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw ResourceError("too many vertices");
    int id = static_cast<int>(names.size());
    names.push_back(name);
    ids.emplace(name, id);
    return id;
}
struct Memory {
    Budget &budget;
    std::uint64_t bytes = 0;
    explicit Memory(Budget &b) : budget(b) {}
    void add(std::uint64_t n) {
        budget.acquire(n);
        bytes += n;
    }
    void remove(std::uint64_t n) {
        budget.release(n);
        bytes -= n;
    }
    ~Memory() {
        budget.release(bytes);
    }
};
inline double elapsed(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
inline std::vector<std::vector<int>> tree_adj(const Decomposition &d) {
    std::vector<std::vector<int>> a(d.bags.size());
    for (auto [u, v] : d.edges) {
        a[u].push_back(v);
        a[v].push_back(u);
    }
    return a;
}
} // namespace td3::detail
