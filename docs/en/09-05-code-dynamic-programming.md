# Dynamic Programming in Python

## Overview

This chapter turns the dynamic programming (DP) models of
{doc}`05-01-dp-fundamentals`, {doc}`05-02-dp-equipment-replacement`,
{doc}`05-03-dp-resource-allocation` and {doc}`05-04-dp-shortest-paths` into Python.
A DP model consists of stages, states, decisions, a transition function, a
stage cost (or return) and boundary conditions ({prf:ref}`def-05-01-dp-model`),
and once these are written down the recursion is always the same. The
chapter therefore writes one small **stage-wise solver** for the backward
recursion {eq}`eq-05-01-recurrence` and one for the forward recursion, both
filling the tables $V$ and $P$ and counting the additions and comparisons
they perform, and a second style, **memoized recursion** with
`functools.lru_cache`, in which the recurrence is written as a recursive
function and Python remembers the values already computed. Both are applied
to the staged diamond network with and without a turn penalty, to equipment
replacement in its two formulations, to one- and two-resource allocation,
the knapsack problem and capital budgeting, and to shortest paths in
acyclic networks with a topological order. Every DP result is checked by
brute-force enumeration, and the measured operation counts are compared with
the formulas of the theory chapters. The standard library is enough;
NetworkX ({doc}`09-01-code-setup-solvers`) is used once, as a check. The
code is not run when the book is
built: copy the blocks, in order, into a script or a notebook and compare
the printed results with the `# Expected:` comments.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- encode a DP model as Python functions for states, decisions, transitions,
  costs and boundary values, and solve it with a generic backward or
  forward stage-wise solver;
- write the same recurrence as a memoized recursive function with
  `functools.lru_cache`;
- recover one or all optimal decision sequences from the policy table;
- check a DP solution by brute-force enumeration and measure its
  operation count against the predicted one; and
- compute a topological order and shortest paths in an acyclic network.
:::

(sec-09-05-setup)=
## Setup

```python
import itertools
import math
import random
from collections import deque
from functools import lru_cache

import networkx as nx

INF = math.inf
```

(sec-09-05-solver)=
## A stage-wise DP solver

The backward recursion {eq}`eq-05-01-recurrence` evaluates, for every stage
$k$ from the last to the first and every state $y\in\mathcal Y_k$,

$$
V_k(y)=\min_{d\in\mathcal D_k(y)}\bigl\{c_k(y,d)+V_{k+1}\bigl(t_k(y,d)\bigr)\bigr\},
$$

starting from the boundary values $V_{N}(y)=g(y)$. The solver below takes
the elements of the model as Python functions, so that one piece of code
serves every example of this chapter:

- `states(k)`: the states $\mathcal Y_k$ of stage $k$;
- `decisions(k, y)`: the allowed decisions $\mathcal D_k(y)$;
- `transition(k, y, d)`: the next state $t_k(y,d)$;
- `cost(k, y, d)`: the stage cost $c_k(y,d)$ (or the stage return when
  maximizing);
- `terminal(y)`: the boundary value $g(y)$ at the stage after the last.

The policy table keeps **all** optimal decisions of each state, because ties
produce alternative optima, which the examples ask for. Every evaluation
$c_k(y,d)+V_{k+1}(\cdot)$ counts as one addition, and choosing the best of
$m$ options as $m-1$ comparisons, the conventions used in the effort
calculations of {doc}`05-01-dp-fundamentals`.

```python
def solve_backward(
    first, last, states, decisions, transition, cost, terminal, maximize=False
):
    """Backward DP recursion over the stages first, ..., last.

    Parameters
    ----------
    first, last : int
        First and last stage with a decision; stage ``last + 1`` holds the
        boundary values.
    states : callable
        ``states(k)`` -> iterable of the states of stage k.
    decisions : callable
        ``decisions(k, y)`` -> iterable of the decisions allowed in y.
    transition : callable
        ``transition(k, y, d)`` -> state of stage k + 1.
    cost : callable
        ``cost(k, y, d)`` -> stage cost (stage return if maximize).
    terminal : callable
        ``terminal(y)`` -> boundary value V_{last+1}(y).
    maximize : bool
        Maximize returns instead of minimizing costs.

    Returns
    -------
    V : dict
        ``V[k][y]``, the optimal value function.
    P : dict
        ``P[k][y]``, the list of all optimal decisions.
    ops : dict
        Numbers of additions and comparisons performed.
    """
    if last < first:
        raise ValueError("need first <= last")
    best = max if maximize else min
    V = {last + 1: {y: terminal(y) for y in states(last + 1)}}
    P = {}
    ops = {"additions": 0, "comparisons": 0}
    for k in range(last, first - 1, -1):
        V[k], P[k] = {}, {}
        for y in states(k):
            options = {}
            for d in decisions(k, y):
                options[d] = cost(k, y, d) + V[k + 1][transition(k, y, d)]
                ops["additions"] += 1
            if not options:
                raise ValueError(f"state {y} of stage {k} has no decision")
            ops["comparisons"] += len(options) - 1
            V[k][y] = best(options.values())
            P[k][y] = [d for d, v in options.items() if v == V[k][y]]
    return V, P, ops


def optimal_plans(P, first, last, y, transition):
    """All optimal decision sequences from state y of stage `first`."""
    if first > last:
        return [[]]
    return [
        [d] + rest
        for d in P[first][y]
        for rest in optimal_plans(
            P, first + 1, last, transition(first, y, d), transition
        )
    ]
```

The forward recursion is written in the "reaching" form: instead of asking
where the best path *to* a state comes from, every state of stage $k$
passes its value on to the states it can reach. The result is the same
function $F_k(y)$ as in {eq}`eq-05-01-lattice-forward`, with a list of best
predecessors for each state.

```python
def solve_forward(first, last, y0, decisions, transition, cost):
    """Forward DP: F[k][y] = least cost from (first, y0) to state y of k.

    Returns F, the best predecessors ``pred[k][y] = [(y_prev, d), ...]``
    and the operation counts.
    """
    F = {first: {y0: 0}}
    pred = {first: {y0: []}}
    ops = {"additions": 0, "comparisons": 0}
    for k in range(first, last + 1):
        F[k + 1], pred[k + 1] = {}, {}
        for y, value in F[k].items():
            for d in decisions(k, y):
                y_next = transition(k, y, d)
                candidate = value + cost(k, y, d)
                ops["additions"] += 1
                if y_next in F[k + 1]:
                    ops["comparisons"] += 1
                if candidate < F[k + 1].get(y_next, INF):
                    F[k + 1][y_next] = candidate
                    pred[k + 1][y_next] = [(y, d)]
                elif candidate == F[k + 1].get(y_next):  # inf: not reached
                    pred[k + 1][y_next].append((y, d))
    return F, pred, ops
```

(sec-09-05-lattice)=
## The staged diamond network

### Data

The network of {prf:ref}`ex-05-01-diamond` is stored with the vertex names of
{numref}`fig-05-01-diamond-network` and their lattice coordinates $(k,y)$:
column $k$, level $y$. The decision "up" (U) leads from $(k,y)$ to
$(k+1,y+1)$, "down" (D) to $(k+1,y-1)$, and `arc_time` returns
$c^{\uparrow}(k,y)$ or $c^{\downarrow}(k,y)$, or $\infty$ when the arc does
not exist.

```python
coord = {
    "A": (0, 0), "C": (1, 1), "D": (1, -1), "E": (2, 2), "F": (2, 0),
    "G": (2, -2), "H": (3, 3), "I": (3, 1), "J": (3, -1), "K": (3, -3),
    "L": (4, 2), "M": (4, 0), "N": (4, -2), "P": (5, 1), "Q": (5, -1),
    "B": (6, 0),
}  # fmt: skip
name = {kv: v for v, kv in coord.items()}
# travel times in minutes, by (tail, head)
times = {
    ("A", "C"): 1, ("A", "D"): 0, ("C", "E"): 5, ("C", "F"): 4,
    ("D", "F"): 7, ("D", "G"): 3, ("E", "H"): 2, ("E", "I"): 1,
    ("F", "I"): 1, ("F", "J"): 2, ("G", "J"): 5, ("G", "K"): 4,
    ("H", "L"): 3, ("I", "L"): 3, ("I", "M"): 4, ("J", "M"): 2,
    ("J", "N"): 2, ("K", "N"): 2, ("L", "P"): 5, ("M", "P"): 2,
    ("M", "Q"): 8, ("N", "Q"): 4, ("P", "B"): 2, ("Q", "B"): 1,
}  # fmt: skip
STEP = {"U": 1, "D": -1}


def arc_time(k: int, y: int, d: str) -> float:
    """c_up(k, y) for d = 'U', c_down(k, y) for d = 'D'; inf if no arc."""
    if d not in STEP:
        raise ValueError(f"decision must be 'U' or 'D', not {d!r}")
    tail, head = name.get((k, y)), name.get((k + 1, y + STEP[d]))
    return times.get((tail, head), INF)


def levels(k: int) -> list:
    """States of stage k: the levels of the vertices in column k."""
    return sorted((y for kk, y in coord.values() if kk == k), reverse=True)


def lattice_moves(k: int, y: int) -> list:
    """Decisions with an arc out of vertex (k, y)."""
    return [d for d in "UD" if arc_time(k, y, d) < INF]


def lattice_step(k: int, y: int, d: str) -> int:
    """Transition t_k(y, d) = y + 1 or y - 1."""
    return y + STEP[d]
```

### Backward and forward recursion

With these functions the model of {ref}`sec-05-01-backward` is one call:
stages 0 to 5, boundary $V_6(0)=0$, answer cell $V_0(0)$.

```python
V, P, ops = solve_backward(
    0, 5, levels, lattice_moves, lattice_step, arc_time, terminal=lambda y: 0
)
print(V[0][0], ops)
# Expected: 13 {'additions': 24, 'comparisons': 9}
print(V[2][-2], P[2][-2], V[1][1], P[1][1])
# Expected: 11 ['U', 'D'] 12 ['D']
for plan in optimal_plans(P, 0, 5, 0, lattice_step):
    y, path = 0, ["A"]
    for k, d in enumerate(plan):
        y = lattice_step(k, y, d)
        path.append(name[(k + 1, y)])
    print("".join(plan), "-".join(path))
# Expected: UDDUUD A-C-F-J-M-P-B
```

The least time is 13 min on $A$–$C$–$F$–$J$–$M$–$P$–$B$, the tie at $G$
appears in the policy table as `['U', 'D']`, and the counts are those of
{prf:ref}`ex-05-01-backward`: 24 additions (one per arc) and 9 comparisons
(one per vertex with two outgoing arcs). The forward recursion gives
$F_6(0)=13$, and $F_k(y)+V_k(y)=13$ exactly on the optimal path.

```python
F, pred, ops_f = solve_forward(0, 5, 0, lattice_moves, lattice_step, arc_time)
print(F[6][0], ops_f)
# Expected: 13 {'additions': 24, 'comparisons': 9}
on_path = [v for v, (k, y) in coord.items() if F[k][y] + V[k][y] == 13]
print(on_path)
# Expected: ['A', 'C', 'F', 'J', 'M', 'P', 'B']
print(F[4][-2], pred[4][-2])
# Expected: 9 [(-1, 'D'), (-3, 'U')]
```

The tie at $N$ (best predecessor $J$ or $K$) is the one noted in
{numref}`tab-05-01-forward`.

### Memoized recursion

The same recurrence {eq}`eq-05-01-lattice-backward` can be written as a
recursive function. Without precautions the recursion would solve the
subproblem of a vertex such as $F$ once for every path that reaches it,
which is complete enumeration again. The decorator `lru_cache` stores
every value the first time it is computed and returns the stored value on
later calls, so each vertex is evaluated once: memoization turns the
recursion into DP. The order of evaluation is now chosen by Python (depth
first from the answer cell), and only the states reachable from the answer
cell are evaluated.

```python
@lru_cache(maxsize=None)
def V_memo(k: int, y: int) -> float:
    """Least time from vertex (k, y) to B, by memoized recursion."""
    if (k, y) == coord["B"]:
        return 0
    return min(
        arc_time(k, y, d) + V_memo(k + 1, lattice_step(k, y, d))
        for d in lattice_moves(k, y)
    )


print(V_memo(0, 0), V_memo.cache_info())
# Expected: 13 CacheInfo(hits=9, misses=16, maxsize=None, currsize=16)
```

The cache holds 16 values, one per vertex, and 9 calls were answered from
it instead of being recomputed.

### Brute-force check

All paths from $A$ to $B$ are the strings of three U and three D that stay
inside the diamond; there are $\binom63=20$ of them.

```python
def path_time(plan: str, turn_penalty: float = 0) -> float:
    """Travel time of a U/D string from A, plus a penalty per turn."""
    y, total = 0, 0
    for k, d in enumerate(plan):
        total += arc_time(k, y, d)
        y = lattice_step(k, y, d)
    turns = sum(a != b for a, b in zip(plan, plan[1:]))
    return total + turn_penalty * turns


plans = ["".join(p) for p in itertools.product("UD", repeat=6)]
plans = [p for p in plans if path_time(p) < INF and p.count("U") == 3]
ranked = sorted((path_time(p), p) for p in plans)
print(len(plans), ranked[:3])
# Expected: 20 [(13, 'UDDUUD'), (14, 'DDDUUU'), (14, 'DDUUUD')]
```

### State augmentation: the turn penalty

In {prf:ref}`ex-05-01-turn-penalty` every change of direction costs 3 min, and
the state must record the direction $\omega$ in which a vertex is left. A
state of stage $k$ is now a pair $(y,\omega)$ for which the arc exists;
the decision is the direction $d$ in which the *next* vertex is left, which
costs the arc time plus 3 if $d\ne\omega$, as in {eq}`eq-05-01-turn`. The
boundary is at stage 5: $V_5(y,\omega)$ is the time of the last arc.

```python
def turn_states(k: int) -> list:
    """Augmented states (y, omega) with an arc out of (k, y) in omega."""
    return [(y, w) for y in levels(k) for w in lattice_moves(k, y)]


def turn_moves(k: int, state: tuple) -> list:
    """Directions in which the next vertex can be left."""
    y, w = state
    return lattice_moves(k + 1, lattice_step(k, y, w))


def turn_step(k: int, state: tuple, d: str) -> tuple:
    """Next state: the next vertex, left in direction d."""
    y, w = state
    return (lattice_step(k, y, w), d)


def turn_cost(k: int, state: tuple, d: str) -> float:
    """Arc time out of (k, y) plus 3 min if the traveler turns next."""
    y, w = state
    return arc_time(k, y, w) + (3 if d != w else 0)


Vt, Pt, ops_t = solve_backward(
    0, 4, turn_states, turn_moves, turn_step, turn_cost,
    terminal=lambda state: arc_time(5, *state),
)  # fmt: skip
print(Vt[0], ops_t)
# Expected: {(0, 'U'): 20, (0, 'D'): 17} {'additions': 34, 'comparisons': 12}
print(Vt[3][(-1, "D")], Pt[3][(-1, "D")], Vt[4][(0, "U")], Pt[4][(0, "U")])
# Expected: 10 ['U'] 7 ['D']
start = min(Vt[0], key=Vt[0].get)
plan = start[1] + "".join(optimal_plans(Pt, 0, 4, start, turn_step)[0])
print(plan, path_time(plan), path_time(plan, turn_penalty=3))
# Expected: DDDUUU 14 17
print(sorted(path_time(p, turn_penalty=3) for p in plans)[:3])
# Expected: [17, 20, 20]
```

Leaving $A$ downward gives 17 min, against 20 for leaving upward; the
optimal path DDDUUU is $A$–$D$–$G$–$K$–$N$–$Q$–$B$ with 14 min of travel
and one turn. The entries $V_3(-1,\mathrm D)=10$ (U) and
$V_4(0,\mathrm U)=7$ (D) are those computed in the example, and the brute
force confirms the unique optimum, followed by two paths of 20 min. The
augmented model has one state per arc instead of one per vertex: 22 states
in stages 0 to 4 instead of 15 vertices, and 34 additions and 12
comparisons instead of 24 and 9.

### Operation counts: measured and predicted

For a diamond lattice with $N$ stages ($N$ even, $m=N/2$),
{eq}`eq-05-01-count-dp` predicts $2m(m+1)$ additions and $m^2$ comparisons for
DP, and {eq}`eq-05-01-count-enumeration` predicts $(N-1)\binom{N}{N/2}$
additions for complete enumeration. The function below builds a lattice of
any even size with random arc times and runs the same solver on it.

```python
def random_lattice(N: int, seed: int = 0):
    """Diamond lattice with N stages and random integer arc times."""
    if N % 2:
        raise ValueError("N must be even")
    rng = random.Random(seed)

    def inside(k, y):
        """True if level y is a vertex of column k."""
        return abs(y) <= min(k, N - k) and (k + y) % 2 == 0

    def states(k):
        """Levels of the vertices of column k."""
        return [y for y in range(-N, N + 1) if inside(k, y)]

    t = {
        (k, y, d): rng.randint(1, 9)
        for k in range(N)
        for y in states(k)
        for d in "UD"
        if inside(k + 1, y + STEP[d])
    }

    def moves(k, y):
        """Directions with an arc out of (k, y)."""
        return [d for d in "UD" if (k, y, d) in t]

    return states, moves, lambda k, y, d: t[(k, y, d)]


for N in [6, 10, 20, 30]:
    states, moves, times_N = random_lattice(N)
    _, _, ops_N = solve_backward(
        0, N - 1, states, moves, lattice_step, times_N, lambda y: 0
    )
    m = N // 2
    print(N, ops_N, 2 * m * (m + 1), m * m, (N - 1) * math.comb(N, m))
# Expected:
#   6 {'additions': 24, 'comparisons': 9} 24 9 100
#   10 {'additions': 60, 'comparisons': 25} 60 25 2268
#   20 {'additions': 220, 'comparisons': 100} 220 100 3510364
#   30 {'additions': 480, 'comparisons': 225} 480 225 4498408080
```

The measured counts equal the predictions of {numref}`tab-05-01-counts` for
every $N$, while the enumeration count grows from 100 to about
$4.5\times10^9$ additions.

(sec-09-05-replacement)=
## Equipment replacement

### Formulation with year and age

The data of {prf:ref}`ex-05-02-replacement` are in thousands of dollars, by
age $i$: operating cost $c(i)$, trade-in value $t(i)$ and salvage value
$s(i)$; a new machine costs $p=50$, the horizon is $N=5$ years and the
machine is $i_0=2$ years old.

```python
N, i0, price = 5, 2, 50
c_op = {0: 10, 1: 13, 2: 20, 3: 40, 4: 70, 5: 100, 6: 100}
trade_in = {1: 32, 2: 21, 3: 11, 4: 5, 5: 0, 6: 0}
salvage = {1: 25, 2: 17, 3: 8, 4: 0, 5: 0, 6: 0, 7: 0}


def ages(k: int) -> list:
    """Possible ages I_k = {1, ..., k-1} and i0 + k - 1 at the start of k."""
    return sorted(set(range(1, k)) | {i0 + k - 1})


def keep_or_buy(k: int, i: int) -> list:
    """Decisions in every state: buy a new machine (B) or keep (K)."""
    return ["B", "K"]


def age_next(k: int, i: int, d: str) -> int:
    """Age at the start of year k + 1: 1 after buying, i + 1 after keeping."""
    if d not in ("B", "K"):
        raise ValueError(f"decision must be 'B' or 'K', not {d!r}")
    return 1 if d == "B" else i + 1


def year_cost(k: int, i: int, d: str) -> float:
    """Cost of year k: p - t(i) + c(0) when buying, c(i) when keeping."""
    return price - trade_in[i] + c_op[0] if d == "B" else c_op[i]
```

The recurrence {eq}`eq-05-02-recurrence` with the boundary
{eq}`eq-05-02-boundary`, $V_6(i)=-s(i)$, is again one call.

```python
Vr, Pr, ops_r = solve_backward(
    1, N, ages, keep_or_buy, age_next, year_cost, lambda i: -salvage[i]
)
print(Vr[1][i0], ops_r)
# Expected: 115 {'additions': 30, 'comparisons': 15}
print(Vr[5][3], Pr[5][3], Vr[2][1], Pr[2][1])
# Expected: 24 ['B'] 76 ['B', 'K']
print(["".join(p) for p in optimal_plans(Pr, 1, N, i0, age_next)])
# Expected: ['BBKBK', 'BKBBK']
```

The minimum cost is \$115 000 with the two plans $\mathrm{BBKBK}$ and $\mathrm{BKBBK}$ of the
example; $V_5(3)=24$ (buy) and the tie $V_2(1)=76$ appear as computed in
{numref}`tab-05-02-two-arg`. The work is $N(N+1)=30$ additions and
$N(N+1)/2=15$ comparisons ({ref}`sec-05-02-effort`). The brute force evaluates
all $2^5=32$ plans.

```python
def plan_cost(plan: str) -> float:
    """Total cost of a keep/buy plan, net of the final salvage value."""
    i, total = i0, 0
    for k, d in enumerate(plan, start=1):
        total += year_cost(k, i, d)
        i = age_next(k, i, d)
    return total - salvage[i]


ranked = sorted(
    (plan_cost("".join(p)), "".join(p))
    for p in itertools.product("BK", repeat=N)
)
print(ranked[:3], plan_cost("KKKKK"))
# Expected: [(115, 'BBKBK'), (115, 'BKBBK'), (116, 'BKBKK')] 330
```

### Formulation with one argument

In the formulation {eq}`eq-05-02-one-arg` the value $U_k$ depends on the year
alone: the firm holds a one-year-old machine at the start of year $k$ and
decides how long to keep it. The recurrence reaches forward to any later
year, so it is natural to write it as a memoized recursive function.
Running sums of the operating costs are not even needed for this small
example.

```python
@lru_cache(maxsize=None)
def U(k: int) -> float:
    """Least cost of years k..N starting year k with a one-year-old machine."""
    if k == N + 1:
        return -salvage[1]
    keep_to_end = (
        sum(c_op[a] for a in range(1, N + 2 - k)) - salvage[N + 2 - k]
    )
    trade_at = [
        sum(c_op[a] for a in range(1, j - k + 1))
        + price - trade_in[j - k + 1] + c_op[0] + U(j + 1)
        for j in range(k, N + 1)
    ]  # fmt: skip
    return min([keep_to_end] + trade_at)


def start_value() -> float:
    """V_1(i0): keep the original machine j - 1 years, then trade it."""
    keep_to_end = sum(c_op[i0 + a] for a in range(N)) - salvage[i0 + N]
    trade_at = [
        sum(c_op[i0 + a] for a in range(j - 1))
        + price - trade_in[i0 + j - 1] + c_op[0] + U(j + 1)
        for j in range(1, N + 1)
    ]  # fmt: skip
    return min([keep_to_end] + trade_at)


print([U(k) for k in range(N + 1, 1, -1)], start_value())
# Expected: [-25, -4, 24, 48, 76] 115
print(all(U(k) == Vr[k][1] for k in range(2, N + 2)))
# Expected: True
```

The values $U_6,\dots,U_2=-25,-4,24,48,76$ and $V_1(2)=115$ are those of
{numref}`tab-05-02-one-arg`, and $U_k=V_k(1)$ for every $k$, as the definition
requires.

(sec-09-05-allocation)=
## Resource allocation

### One resource

In {prf:ref}`ex-05-03-one-resource` eight units are allocated among four
activities with the returns of {numref}`tab-05-03-returns`. The activities
are the stages, the units left $y$ the state, and the allocation $x\le y$
the decision. With the equality constraint the boundary is
$V_4(y)=r_4(y)$ ({eq}`eq-05-03-boundary`), so the decisions run over stages 1
to 3, and only the state $y=8$ is needed at stage 1.

```python
returns = [
    [0, 3, 7, 10, 12, 13, 14, 14, 14],
    [0, 1, 2, 4, 8, 13, 17, 19, 20],
    [0, 2, 4, 6, 8, 10, 12, 14, 16],
    [0, 1, 3, 6, 9, 12, 14, 16, 17],
]
X = 8
Va, Pa, ops_a = solve_backward(
    1, 3,
    states=lambda k: [X] if k == 1 else range(X + 1),
    decisions=lambda k, y: range(y + 1),
    transition=lambda k, y, x: y - x,
    cost=lambda k, y, x: returns[k - 1][x],
    terminal=lambda y: returns[3][y],
    maximize=True,
)  # fmt: skip
plan = optimal_plans(Pa, 1, 3, X, lambda k, y, x: y - x)
print(Va[1][X], plan, ops_a["additions"])
# Expected: 24 [[2, 6, 0]] 99
```

The plan gives $x_1=2$, $x_2=6$, $x_3=0$, and $x_4=8-2-6-0=0$: return 24,
unique. The 99 additions are the count of {ref}`sec-05-03-one-resource`: two
full stages of $(X+1)(X+2)/2=45$ options and 9 options at stage 1. The
brute force runs over all $\binom{11}{3}=165$ allocations of exactly eight
units.

```python
allocations = [
    x for x in itertools.product(range(X + 1), repeat=4) if sum(x) == X
]
values = [sum(returns[k][x[k]] for k in range(4)) for x in allocations]
best = max(values)
print(
    len(allocations),
    best,
    [x for x, v in zip(allocations, values) if v == best],
)
# Expected: 165 24 [(2, 6, 0, 0)]
```

### Two resources

With two resources the state is the pair $(y_1,y_2)$ of units left
({eq}`eq-05-03-two-resources`). The data of {prf:ref}`ex-05-03-two-resources` are
three $4\times4$ tables $r_k(x_1,x_2)$, rows $x_1$, columns $x_2$.

```python
r2 = [
    [
        [0, 1, 3, 6],
        [4, 5, 6, 7],
        [5, 6, 7, 8],
        [6, 7, 8, 9],
    ],
    [
        [0, 2, 4, 6],
        [1, 4, 6, 7],
        [4, 6, 8, 9],
        [6, 8, 10, 11],
    ],
    [
        [0, 3, 5, 8],
        [2, 5, 7, 9],
        [4, 7, 9, 11],
        [6, 9, 11, 13],
    ],
]
X1 = X2 = 3


def pairs_up_to(y: tuple) -> list:
    """All allocations (x1, x2) <= (y1, y2)."""
    return [(a, b) for a in range(y[0] + 1) for b in range(y[1] + 1)]


V2, P2, ops_2 = solve_backward(
    1, 2,
    states=lambda k: [(X1, X2)] if k == 1 else pairs_up_to((X1, X2)),
    decisions=lambda k, y: pairs_up_to(y),
    transition=lambda k, y, x: (y[0] - x[0], y[1] - x[1]),
    cost=lambda k, y, x: r2[k - 1][x[0]][x[1]],
    terminal=lambda y: r2[2][y[0]][y[1]],
    maximize=True,
)  # fmt: skip
print(V2[1][(X1, X2)], P2[1][(X1, X2)], P2[2][(2, 3)])
# Expected: 16 [(1, 0)] [(2, 0)]
print(ops_2["additions"] - 16, V2[2][(1, 1)], V2[2][(3, 1)])
# Expected: 100 5 9
```

Activity 1 gets $(1,0)$, activity 2 then gets $(2,0)$ from $(2,3)$, and
activity 3 the rest, $(0,3)$: return $4+4+8=16$. Stage 2 alone needed
$100$ additions, the count {eq}`eq-05-03-two-count`; with one resource of
three units it would have been 10. The values $V_2(1,1)=5$ and $V_2(3,1)=9$
agree with the example.

### The knapsack problem

The 0–1 knapsack problem {eq}`eq-05-03-knapsack` is short enough to write as a
memoized recursion. In {prf:ref}`ex-08-07-a8-q4` five housing applications
have need scores $h_k$ and need $a_k$ volunteers, and at most 23 volunteers
are available. `best_score(k, y)` is $V_k(y)$, the best total score of
applications $k,\dots,5$ with $y$ volunteers left.

```python
score = [78, 64, 68, 62, 85]
volunteers = [7, 4, 6, 5, 8]


@lru_cache(maxsize=None)
def best_score(k: int, y: int) -> int:
    """V_k(y) of the 0-1 knapsack, applications numbered 1, ..., 5."""
    if y < 0:
        raise ValueError("capacity must be nonnegative")
    if k > len(score):
        return 0
    reject = best_score(k + 1, y)
    if volunteers[k - 1] > y:
        return reject
    return max(reject, score[k - 1] + best_score(k + 1, y - volunteers[k - 1]))


def approved(y: int) -> list:
    """Trace the decisions of best_score forward from stage 1."""
    chosen = []
    for k in range(1, len(score) + 1):
        if best_score(k, y) != best_score(k + 1, y):
            chosen.append(k)
            y -= volunteers[k - 1]
    return chosen


print(best_score(1, 23), approved(23), best_score.cache_info().currsize)
# Expected: 279 [2, 3, 4, 5] 51
subsets = [
    S
    for n in range(6)
    for S in itertools.combinations(range(1, 6), n)
    if sum(volunteers[k - 1] for k in S) <= 23
]
print(max(sum(score[k - 1] for k in S) for S in subsets))
# Expected: 279
```

Applications 2, 3, 4 and 5 are approved for a total score of 279, using all
23 volunteers. The recursion evaluated only the states it reached from
$(1,23)$, far fewer than the $5\times24$ cells of a full table.

### Capital budgeting

In {prf:ref}`ex-05-03-capital-budgeting` four activities buy units of two
resources at \$3 and \$1 from a budget of \$6. The return tables
$r_k(x_1,x_2)$ are given for the affordable pairs only; a row of the list
below holds $r_k(x_1,0),r_k(x_1,1),\dots$.

```python
budget_returns = [
    [
        [0, 5, 9, 11, 12, 13, 14],
        [5, 9, 13, 18],
        [10],
    ],
    [
        [0, 3, 5, 6, 7, 8, 9],
        [2, 5, 9, 12],
        [7],
    ],
    [
        [0, 1, 3, 6, 10, 12, 13],
        [4, 6, 8, 11],
        [11],
    ],
    [
        [0, 1, 3, 6, 10, 13, 14],
        [7, 10, 14, 17],
        [13],
    ],
]
c1, c2, budget = 3, 1, 6


def affordable(k: int, b: int) -> list:
    """Pairs (x1, x2) with a return and cost c1*x1 + c2*x2 <= b."""
    table = budget_returns[k - 1]
    return [
        (x1, x2)
        for x1 in range(len(table))
        for x2 in range(len(table[x1]))
        if c1 * x1 + c2 * x2 <= b
    ]
```

**Formulation 2** ({eq}`eq-05-03-budget-2`) first computes the best return
$\hat r_k(\beta)$ of each activity from at most $\beta$ dollars
({eq}`eq-05-03-hat-r`), together with the pairs that reach it, and then
allocates money as a single resource with the boundary $V_5(b)=0$ (unspent
money allowed).

```python
def best_return(k: int, beta: int) -> tuple[int, list]:
    """r_hat_k(beta) and the pairs (x1, x2) that attain it."""
    options = {
        x: budget_returns[k - 1][x[0]][x[1]] for x in affordable(k, beta)
    }
    top = max(options.values())
    return top, [x for x, v in options.items() if v == top]


Vb, Pb, ops_b = solve_backward(
    1, 4,
    states=lambda k: [budget] if k == 1 else range(budget + 1),
    decisions=lambda k, b: range(b + 1),
    transition=lambda k, b, beta: b - beta,
    cost=lambda k, b, beta: best_return(k, beta)[0],
    terminal=lambda b: 0,
    maximize=True,
)  # fmt: skip
print(Vb[1][budget], [best_return(4, beta)[0] for beta in range(budget + 1)])
# Expected: 19 [0, 1, 3, 7, 10, 14, 17]
money_plans = optimal_plans(Pb, 1, 4, budget, lambda k, b, beta: b - beta)
plans = [
    pairs
    for betas in money_plans
    for pairs in itertools.product(
        *(best_return(k, beta)[1] for k, beta in enumerate(betas, start=1))
    )
]
for plan in plans:
    print(plan)
# Expected:
#   ((0, 1), (0, 0), (0, 0), (1, 2))
#   ((0, 2), (0, 0), (0, 0), (0, 4))
#   ((0, 2), (0, 0), (0, 0), (1, 1))
#   ((0, 2), (0, 0), (0, 4), (0, 0))
#   ((0, 2), (0, 1), (0, 0), (1, 0))
```

The maximum return is 19, and the five plans are those of the table in
{prf:ref}`ex-05-03-capital-budgeting`. **Formulation 1** ({eq}`eq-05-03-budget-1`)
keeps both purchases as the decision of a stage; it gives the same value
function, with every affordable pair evaluated for every amount of money
left.

```python
Vb1, Pb1, ops_b1 = solve_backward(
    1, 4,
    states=lambda k: range(budget + 1),
    decisions=affordable,
    transition=lambda k, b, x: b - c1 * x[0] - c2 * x[1],
    cost=lambda k, b, x: budget_returns[k - 1][x[0]][x[1]],
    terminal=lambda b: 0,
    maximize=True,
)  # fmt: skip
print(
    all(Vb1[k][b] == Vb[k][b] for k in range(2, 5) for b in range(budget + 1))
)
# Expected: True
print(Vb1[1][budget], ops_b1["additions"] // 4)
# Expected: 19 39
```

Per stage, formulation 1 performs 39 evaluations. The evaluation counts of
{numref}`tab-05-03-budget-effort` depend only on the prices and the budget,
so they can be computed for larger budgets without data: formulation 1
evaluates, for every $b\le B$, all pairs with $3x_1+x_2\le b$;
formulation 2 evaluates each pair with $3x_1+x_2\le B$ once and then the
$(B+1)(B+2)/2$ splits of the money.

```python
def n_pairs(b: int) -> int:
    """Number of pairs (x1, x2) >= 0 with c1*x1 + c2*x2 <= b."""
    return sum((b - c1 * x1) // c2 + 1 for x1 in range(b // c1 + 1))


for B in [6, 12, 30, 60, 120]:
    form1 = sum(n_pairs(b) for b in range(B + 1))
    print(B, form1, n_pairs(B), (B + 1) * (B + 2) // 2)
# Expected:
#   6 39 12 28
#   12 185 35 91
#   30 1991 176 496
#   60 13881 651 1891
#   120 103361 2501 7381
```

(sec-09-05-acyclic)=
## Shortest paths in acyclic networks

### Topological order

{prf:ref}`thm-05-04-topological-order` numbers the nodes of an acyclic
network so that every arc points forward: repeatedly take a node without
incoming arcs (a source node, {prf:ref}`lem-05-04-source-node`) and delete it
with its arcs. The function keeps the number of incoming arcs of every node
and a queue of the current source nodes; if the queue empties before all
nodes are numbered, the remaining nodes contain a cycle.

```python
def topological_order(nodes: list, arcs: list) -> list:
    """Topological order by repeated deletion of source nodes.

    Raises ValueError if the network has a cycle.
    """
    indegree = {v: 0 for v in nodes}
    succ = {v: [] for v in nodes}
    for i, j, _ in arcs:
        indegree[j] += 1
        succ[i].append(j)
    sources = deque(v for v in nodes if indegree[v] == 0)
    order = []
    while sources:
        i = sources.popleft()
        order.append(i)
        for j in succ[i]:
            indegree[j] -= 1
            if indegree[j] == 0:
                sources.append(j)
    if len(order) < len(nodes):
        raise ValueError("the network has a cycle")
    return order
```

### The one-pass recursion

In topological order the forward recursion {eq}`eq-05-04-acyclic-forward`
needs one pass: $V(j)=\min_i\{V(i)+c_{ij}\}$ over the arcs into $j$, with
$V$ of the origin equal to 0. The backward recursion
{eq}`eq-05-04-acyclic-backward` runs through the order in reverse. Both count
their additions and comparisons.

```python
def acyclic_shortest_paths(nodes: list, arcs: list, s):
    """Distances V and predecessors P from s in an acyclic network."""
    order = topological_order(nodes, arcs)
    into = {v: [] for v in nodes}
    for i, j, c in arcs:
        into[j].append((i, c))
    V = {v: INF for v in nodes}
    P = {v: None for v in nodes}
    V[s] = 0
    ops = {"additions": 0, "comparisons": 0}
    for j in order[order.index(s) + 1 :]:
        options = [(V[i] + c, i) for i, c in into[j] if V[i] < INF]
        ops["additions"] += len(options)
        ops["comparisons"] += max(len(options) - 1, 0)
        if options:
            V[j], P[j] = min(options, key=lambda o: o[0])
    return V, P, ops


def acyclic_to_target(nodes: list, arcs: list, t) -> dict:
    """Distances V_bar(i) from every node i to t (backward recursion)."""
    order = topological_order(nodes, arcs)
    out = {v: [] for v in nodes}
    for i, j, c in arcs:
        out[i].append((j, c))
    V_bar = {v: INF for v in nodes}
    V_bar[t] = 0
    for i in reversed(order):
        for j, c in out[i]:
            V_bar[i] = min(V_bar[i], c + V_bar[j])
    return V_bar
```

For the network of {prf:ref}`ex-05-04-acyclic` (lengths in km) the order is
$D,C,A,E,F,B$, the shortest path from $D$ to $B$ is $D$–$C$–$A$–$E$–$F$–$B$
of 8 km, and the recursion uses 9 additions and 4 comparisons. Backtracking
uses `trace` below, the same backward reading of predecessors as in
{doc}`09-04-code-networks`.

```python
acyclic_nodes = ["A", "B", "C", "D", "E", "F"]
acyclic_arcs = [
    ("D", "A", 4), ("D", "C", 2), ("C", "A", 1), ("C", "E", 6),
    ("A", "E", 2), ("A", "F", 5), ("E", "F", 2), ("E", "B", 6),
    ("F", "B", 1),
]  # fmt: skip


def trace(P: dict, j) -> list:
    """Path to j read backwards from the predecessors P."""
    path = [j]
    while P[path[-1]] is not None:
        path.append(P[path[-1]])
    return path[::-1]


print(topological_order(acyclic_nodes, acyclic_arcs))
# Expected: ['D', 'C', 'A', 'E', 'F', 'B']
Vd, Pd, ops_d = acyclic_shortest_paths(acyclic_nodes, acyclic_arcs, "D")
order = topological_order(acyclic_nodes, acyclic_arcs)
print([Vd[v] for v in order], ops_d)
# Expected: [0, 2, 3, 5, 7, 8] {'additions': 9, 'comparisons': 4}
print("-".join(trace(Pd, "B")))
# Expected: D-C-A-E-F-B
V_bar = acyclic_to_target(acyclic_nodes, acyclic_arcs, "B")
print([V_bar[v] for v in order], [Vd[v] + V_bar[v] for v in order])
# Expected: [8, 6, 5, 3, 1, 0] [8, 8, 8, 8, 8, 8]
```

Every node lies on the shortest $D$–$B$ path, so $V(j)+\bar V(j)=8$
everywhere, the check of the example. NetworkX finds the same order and
the same path.

```python
G_acyclic = nx.DiGraph()
G_acyclic.add_weighted_edges_from(acyclic_arcs)
print(
    list(nx.topological_sort(G_acyclic)),
    nx.shortest_path(G_acyclic, "D", "B", weight="weight"),
)
# Expected: ['D', 'C', 'A', 'E', 'F', 'B'] ['D', 'C', 'A', 'E', 'F', 'B']
```

### Equipment replacement as an acyclic shortest path

{ref}`sec-05-02-network` draws the replacement problem as a network with a
node for every pair (year, age), keep and buy arcs, and terminal arcs of
length $-s(i)$. The network is acyclic, so the one-pass recursion applies
despite the negative terminal arcs, and the distance to the terminal node
is the minimum cost of {ref}`sec-09-05-replacement`.

```python
rep_arcs = []
for k in range(1, N + 1):
    for i in ages(k):
        for d in "BK":
            rep_arcs.append(
                ((k, i), (k + 1, age_next(k, i, d)), year_cost(k, i, d))
            )
rep_arcs += [((N + 1, i), "end", -salvage[i]) for i in ages(N + 1)]
rep_nodes = sorted({a for arc in rep_arcs for a in arc[:2]}, key=str)
V_rep, P_rep, _ = acyclic_shortest_paths(rep_nodes, rep_arcs, (1, i0))
print(len(rep_nodes), len(rep_arcs), V_rep["end"], trace(P_rep, "end")[1:-1])
# Expected: 22 36 115 [(2, 1), (3, 2), (4, 1), (5, 1), (6, 2)]
```

The network has 21 nodes and 30 decision arcs plus the terminal node and
its six arcs; the shortest path costs 115 and is one of the two optimal
plans.

### Operation count of the one-pass recursion

In a network with an arc from every node to every higher-numbered node,
{eq}`eq-05-04-acyclic-count` predicts $n(n-1)/2$ additions and
$(n-1)(n-2)/2$ comparisons.

```python
for n in [5, 10, 50]:
    rng = random.Random(n)
    complete = [
        (i, j, rng.randint(1, 20))
        for i in range(1, n + 1)
        for j in range(i + 1, n + 1)
    ]
    _, _, ops_n = acyclic_shortest_paths(list(range(1, n + 1)), complete, 1)
    print(n, ops_n, n * (n - 1) // 2, (n - 1) * (n - 2) // 2)
# Expected:
#   5 {'additions': 10, 'comparisons': 6} 10 6
#   10 {'additions': 45, 'comparisons': 36} 45 36
#   50 {'additions': 1225, 'comparisons': 1176} 1225 1176
```

## Summary

- A DP model is five functions (states, decisions, transition, stage cost,
  boundary value); one backward solver for {eq}`eq-05-01-recurrence` and one
  forward solver produce $V$, $P$ (all optimal decisions) and the operation
  counts for every example.
- Memoized recursion with `lru_cache` evaluates each state once; it suits
  recurrences with long-range decisions ({eq}`eq-05-02-one-arg`) and problems
  where only a few states are reachable (the knapsack).
- The diamond network gives 13 min with 24 additions and 9 comparisons,
  the turn-penalty model 17 min; measured counts equal
  {eq}`eq-05-01-count-dp` for all lattice sizes.
- Equipment replacement gives \$115 000 with the plans $\mathrm{BBKBK}$ and $\mathrm{BKBBK}$ in
  both formulations, and as an acyclic shortest path.
- Resource allocation gives 24 (one resource), 16 (two resources, 100
  options at stage 2), 279 (knapsack) and 19 with five plans (capital
  budgeting); the evaluation counts of {numref}`tab-05-03-budget-effort` follow
  from counting affordable pairs.
- A topological order by deletion of source nodes gives the one-pass
  shortest-path recursion of {doc}`05-04-dp-shortest-paths`; every result is
  confirmed by brute force or NetworkX.

## References

- Bellman, R. (1957). *Dynamic programming*. Princeton University Press.
- Dreyfus, S. E., & Law, A. M. (1977). *The art and theory of dynamic
  programming*. Academic Press. (Chapters 1–4: staged networks, equipment
  replacement, resource allocation, shortest paths.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 11: dynamic programming.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 12: deterministic dynamic programming.)
