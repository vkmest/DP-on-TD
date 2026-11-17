# Correctness and complexity

Let (T, {B_t}) be a validated tree decomposition, rooted arbitrarily. Let V_t be
the union of bags in t's subtree. F_t contains exactly the colorings of B_t that
extend to a proper coloring of G[V_t]. Color identities 0, 1 and 2 are global.

For child u of t, S_u = B_u ∩ B_t. The outgoing table stores

    P_u[s] = one a in F_u with a restricted to S_u equal to s,

or an absent entry. Keys and full-bag witnesses are base-3 codes over sorted
internal vertex IDs. The sorted order is also exported in each report bag.

The recurrence is:

    c in F_t iff c properly colors G[B_t]
                   and P_u[c restricted to S_u] exists for every child u.

At a leaf, local propriety is sufficient. At an internal node, the running-
intersection property ensures that the child subtree can interact with the rest
only through its separator. Distinct child interiors share neither a vertex nor
an edge connecting them outside the parent bag: either would contradict running
intersection and edge coverage. Hence agreeing extensions can be combined. The
reverse direction follows by restricting any proper subtree coloring. Induction
proves the recurrence. A feasible root coloring therefore exists exactly when the
whole graph is 3-colorable.

One witness per separator coloring suffices: outside the child subtree, only the
separator can matter. Reconstruction retrieves a child witness using the chosen
parent colors, then repeats. No all-children pointer list is needed. The finished
coloring is independently checked against every original vertex and edge.

## Bounds for the DP phase

Let b be the bag count and k+1 the maximum bag size. For each bag, enumerate at
most 3^(k+1) assignments; local constraint checks and key construction cost at most
O((k+1)²) per candidate. Each tree edge causes one projection of its child and
lookups from its parent, not a Cartesian product of assignment tables. Summing
child degrees gives b-1, not b². Sparse map construction/sorting/search introduces
at most O(log(3^(k+1))) = O(k+1) per operation, within the conservative bound:

    Time: O(b (k+1)² 3^(k+1)).

Retained messages occupy O(sum_{t != root} 3^|S_t|) words; sparse construction and
its compact replacement coexist temporarily for one bag. Metadata is O(b(k+1)²).
Since the polynomial bag factor is bounded by a constant times 3^(k+1), overall:

    Space: O(n+m+b 3^(k+1)) words.
    Reconstruction: O(b(k+1)).

These recurrence bounds assume bag-local edge constraints and projection positions
have been prepared. The implementation reports that work separately as
`preparation_ms`: constructing local constraints by binary search in sorted graph
adjacency lists takes O(b(k+1)² log(n+1)) time. Including this preparation, the
bound is O(b(k+1)² (log(n+1) + 3^(k+1))) plus graph/decomposition I/O and validation.
This distinction matters even at width one, for example in a large star.

These are bounds in the width **used by the solver**. When every connected
component has at most 20 vertices (the default exact limit), that width is the
treewidth, proven by exact search. Otherwise it is a heuristic upper bound, unless
it equals the lower bound described below. Every run reports which case applies.
The base-3 dependence matches the conditional lower bound for 3-colorability
parameterized by treewidth, assuming SETH; this does not assert optimal polynomial
factors, optimal memory, or fastest runtime.
Source: [Lokshtanov, Marx, Saurabh](https://arxiv.org/abs/1007.5450).

A tighter count also holds but is not needed: enumerating one bag visits at most
1.5 · 3^(k+1) partial assignments, each checking at most k earlier neighbours, and
each child lookup costs O(k+1). Summed over all bags this gives O(b (k+1) 3^(k+1));
the stated (k+1)² bound is the conservative version.

## Construction of the decomposition

The DP is correct for **any** valid tree decomposition. Construction only affects
the running time, through the width k and the bag count b. It is a separate phase
with its own timing (`decomposition_ms`), and it never adds colouring constraints:
fill edges live in a private copy of the adjacency and are discarded.

### From an elimination order to a decomposition

An elimination order lists the vertices v1, ..., vn. Eliminating v makes its
remaining neighbours a clique (the new edges are *fill edges*) and deletes v. The
bag of v is v together with its neighbours at elimination time. Those neighbours
form a clique of later vertices, so all of them lie in the bag of the earliest of
them; attach v's bag there. Bags without later neighbours are roots of separate
components and are linked by edges with an empty intersection. This classical
construction is a tree decomposition of the original graph whose width is the
largest number of later neighbours of any vertex.

### Candidate orders and selection

| `--td-method` | what runs |
|---|---|
| `minfill` | minimum fill-in only (fewest missing edges among neighbours; ties by degree, then ID) — the original behaviour |
| `best` | minimum fill-in, minimum degree, maximum cardinality search, and both greedy rules again with 4 seeded random tie-breaks (graphs up to 2,000 vertices) |
| `auto` (default) | `best`, then exact search on components of at most 20 vertices |
| `exact` | as `auto`, but a component above the limit gives RESOURCE_LIMIT instead of a heuristic |

Greedy scores are cached in an ordered set; after an elimination only vertices whose
fill-in or degree can change are rescored. Seeded tie-breaks use a Fisher-Yates
shuffle driven by `mt19937`, so results are identical on every compiler. Each
candidate is validated, simplified, validated again, and ranked by

1. **width** (the exponent of the DP),
2. **Σ 3^|bag|**, an upper bound on the number of colourings the DP enumerates,
3. bag count.

Width always comes first. One bag containing every vertex would minimise the bag
count but maximise the DP cost, so bag count is only the final tie-break.

### Simplification: merging contained bags

Elimination creates one bag per vertex, and many of them are subsets of a
neighbouring bag. Such bags add no information but still cost up to 3^|bag| DP work.
Simplification repeats the following until no case remains:

> if tree neighbours s and t satisfy B_s ⊆ B_t, delete s and connect every other
> neighbour of s to t.

**Validity.** This contracts the tree edge s–t and labels the merged node with
B_s ∪ B_t = B_t. Contracting an edge keeps a tree a tree. Every vertex and edge
covered by B_s is still covered by B_t. For running intersection, take a vertex x
and the connected subtree T_x of bags containing x. If x ∉ B_s, T_x does not use s
and is unchanged. If x ∈ B_s then also x ∈ B_t, so the edge s–t lies inside T_x, and
contracting an edge inside a connected subtree leaves it connected.

**Width.** Every surviving bag is an original, unchanged bag, so the width cannot
grow. It cannot shrink either: a removed bag B_s is contained in a surviving B_t
with |B_t| ≥ |B_s|, so some bag of maximum size always survives.

**Size.** Afterwards no bag is contained in a neighbour, which forces at most n bags
(n ≥ 1). Root the tree anywhere. Every non-root bag has a vertex x missing from its
parent, so that bag is the top node of T_x. Distinct bags are tops of distinct
subtrees, and the root is the top of the subtrees of its own vertices.

Each merge removes one bag, so there are fewer than b merges; each re-checks the
affected tree edges with `std::includes` in O(k) time. Conservatively this is
O(b² (k + log b)), far below the construction cost in practice.

### Rooting and numbering

Re-rooting never changes the DP cost: a tree edge has the same separator whichever
end is the parent, so every projected table has the same key set. Rooting only
changes how the tree reads. After simplification:

1. edges with an empty separator are removed; they only join independent pieces
   (for example connected components), and removing them leaves valid
   decompositions of disjoint vertex sets that may be re-linked arbitrarily;
2. each piece's centre (minimum eccentricity) is found by peeling leaves layer by
   layer, and every other piece is attached to the centre of the largest piece;
3. bags are renumbered breadth-first from that root: B0 is the root and every parent
   has a smaller number than its children.

The report draws the tree with a tidy layout: leaves take consecutive columns in
depth-first order, each parent is centred over its children, and each tree edge is
labelled with its separator size |B_parent ∩ B_child|.

### Exact minimum width

For a set S of vertices eliminated first and a vertex v ∉ S eliminated next, v's
neighbours in the filled graph are exactly

    Q(S, v) = the vertices outside S ∪ {v} reachable from v by a path
              whose interior lies in S.

Hence the best width of any order that eliminates S first satisfies

    TW(∅) = 0,   TW(S) = min over v ∈ S of max( TW(S − v), |Q(S − v, v)| ),

and a connected component C has treewidth TW(C) (Bodlaender, Fomin, Koster,
Kratsch and Thilikos, *On exact algorithms for treewidth*, ESA 2006). The
implementation stores TW for all 2^|C| subsets as bytes, computes Q by a bitmask
breadth-first search, and remembers the minimising v so that an optimal order can
be read back. It is given the best heuristic width h and only searches for widths
below h; if there are none, h is proven optimal for that component.

* Time O(2^c c²) word operations for a component of c vertices (each closure visits
  every vertex at most once); memory 2 · 2^c bytes, counted in the memory budget.
* Default limit 20 vertices per component (2 MiB; `graph20_3colorable` takes about
  10 ms). `--exact-limit` accepts at most 26; `--exact-limit 0` disables the search.
* Treewidth is the maximum over components and orders of different components can
  be concatenated, so components are searched independently. If a component exceeds
  the limit, `auto` keeps its heuristic order. The result is still proven optimal
  when an exactly solved component already has width ≥ every heuristic component.

### Lower bound and honest reporting

The *minor-min-width* (MMD+) lower bound repeatedly records the minimum degree and
contracts a minimum-degree vertex into its minimum-degree neighbour. Every graph has
treewidth at least its minimum degree, and contraction produces a minor, whose
treewidth is no larger, so the largest recorded degree is a lower bound. It costs
O(m log² n) and is computed on every run.

Each result records `width_status`:

* `optimal`, `proof = exact_search`: exact search covered every component, or a
  component that attains the width;
* `optimal`, `proof = lower_bound_match`: the width equals the MMD+ lower bound;
* `upper_bound`: neither proof was found, so `lower_bound ≤ treewidth ≤ width`.

A supplied decomposition is used as given (or merged with `--simplify-td`) and is
reported as `optimal` only when it matches the lower bound. With peeling enabled,
these statements refer to the residual core, not the input graph.

### Construction costs

One greedy run is conservatively O(n⁴) time with O(n²/64) words of bitsets; cached
scores make typical runs much faster. `best` performs 3 + 2·4 runs on graphs up to
2,000 vertices and 3 runs above that. Maximum cardinality search is O(m log n).
Every candidate is validated in O(n + m + b(k+1) log b). All of this is separate from
the DP bounds above. Supplied decompositions bypass construction.

## Validation

Validation checks a connected acyclic tree, coverage, and running intersection.
For each graph vertex, its occurrence bags induce a forest in T. They form a
connected subtree iff the number of tree edges whose endpoints both contain it
equals the number of occurrences minus one. Edge coverage uses intersections of
the sorted occurrence lists. Its cost is reported independently of DP.

## Optional reductions

After checking self-loops, remove vertices of degree at most two repeatedly.
Color the residual components with the same projection DP. Restore removed
vertices in reverse order: at most two neighbors are already colored, leaving at
least one of three colors available. Therefore each removal preserves whether a
3-coloring exists, and reconstruction is exact. Residual decompositions need not
cover removed vertices; they are explicitly labelled as residual in reports.
Supplied TD and TD export modes preserve the original graph instead.

The empty graph has the empty coloring and width -1. A self-loop is an immediate
negative certificate. Timeout, allocation failure and encoding overflow are
unknown outcomes, never negative certificates.
