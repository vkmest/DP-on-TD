# DP on tree decompositions: 3-colorability

An exact 3-colorability solver in C++20. It runs dynamic programming over a tree
decomposition and reconstructs a coloring when one exists. For each child bag,
the feasible colorings are projected onto the vertices shared with the parent,
and one full coloring is kept as a witness for each projection. The DP takes
**O(b·(k+1)²·3^(k+1))** time for b bags of width k.
Details and proofs: [docs/algorithm.md](docs/algorithm.md).

## Build

You need CMake 3.20+ and a C++20 compiler (GCC 10+, Clang 12+ or Visual Studio 2022).

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release        # optional: run the tests
```

The executable is `build/td3color`. With Visual Studio it is `build\Release\td3color.exe`.

## Use

```sh
td3color solve graphs/graph13.csv                        # verdict + coloring
td3color solve graphs/graph13.csv --report out.html      # interactive offline report
td3color solve graphs/graph13.csv --json out.json        # machine-readable result
td3color solve graphs/graph13.csv --write-td td.csv      # export the decomposition
td3color solve graphs/graph13.csv --td td.csv            # use your own decomposition
td3color validate-td graphs/graph13.csv td.csv           # check a decomposition
```

`solve` prints `COLORABLE` followed by three lines, one per color class, or
`NOT_COLORABLE`. Width, timing and the status of the width are printed to stderr.
Exit codes: `0` colorable or valid, `1` not colorable, `2` invalid input,
`3` resource limit. A resource limit means the answer is unknown.

| Option | Meaning |
|---|---|
| `--td-method auto\|best\|minfill\|exact` | How the decomposition is built. The default `auto` compares several heuristics, then runs an exact minimum-width search on components of up to 20 vertices. |
| `--exact-limit N` | Largest component (in vertices) for exact search. Default 20, max 26, 0 disables it. |
| `--no-simplify` | Keep the raw one-bag-per-vertex decomposition instead of merging contained bags. |
| `--simplify-td` | Also merge contained bags of a decomposition given with `--td`. |
| `--preprocess` / `--no-preprocess` | Remove vertices of degree ≤ 2 before the DP. On by default, off for `--report` and `--write-td`. |
| `--time-limit S`, `--memory-mib M` | Resource limits. Defaults: no time limit, 1024 MiB. |
| `--table auto\|dense\|sparse` | How the DP tables are stored. |

Every run reports whether the decomposition width is **proven optimal** (by exact
search or a matching lower bound) or only a **heuristic upper bound**.

## Input format

Graphs use the Glasgow CSV format: an optional `vertex,` line for each vertex
(this keeps isolated vertices) and one `u,v` line for each edge. There is no
header line.

```csv
a,
b,
c,
a,b
b,c
```

A decomposition lists its tree edges, then gives each bag's vertices separated by `;`:

```csv
B0,B1,
B0,,a;b
B1,,b;c
```

## Examples

`graphs/` contains the test graphs:

- 3-colorable: `graph00`–`graph06`, `graph09`–`graph11`, `graph13` (Wagner),
  `graph15_3colorable`, `graph20_3colorable` and `example_2tree`. The last one has
  a width-2 decomposition that is easy to draw: `example_2tree.td.csv`.
- Not 3-colorable: `graph07` (K4), `graph08`, `graph12` and `groetzsch`.
