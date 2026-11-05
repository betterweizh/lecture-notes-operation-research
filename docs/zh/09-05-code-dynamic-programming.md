# 用 Python 实现动态规划

## 概述

本章把{doc}`05-01-dp-fundamentals`、{doc}`05-02-dp-equipment-replacement`、{doc}`05-03-dp-resource-allocation`和{doc}`05-04-dp-shortest-paths`中的动态规划（dynamic programming, DP）模型转写为 Python 代码。一个 DP 模型由阶段（stage）、状态（state）、决策（decision）、状态转移函数（transition function）、阶段费用（stage cost，或收益）和边界条件（boundary condition）组成（见{prf:ref}`def-05-01-dp-model`），这些要素一旦写出，递推过程总是相同的。因此本章为逆序递推（backward recursion）{eq}`eq-05-01-recurrence`编写一个小型的**逐阶段求解器**，再为顺序递推（forward recursion）编写一个，两者都填写表 $V$ 和 $P$，并统计所做的加法与比较次数；此外还采用第二种风格，即使用 `functools.lru_cache` 的**记忆化递归**（memoized recursion），把递推关系写成递归函数，由 Python 记住已经算出的值。这两种方法用于带转弯惩罚（turn penalty）和不带转弯惩罚的菱形阶段网络、两种形式的设备更新（equipment replacement）问题、单资源与双资源分配、背包问题（knapsack problem）和资本预算（capital budgeting），以及用拓扑序（topological order）求解无圈网络中的最短路。每个 DP 结果都用穷举法核对，实测的运算量（operation count）还与理论各章中的公式相比较。标准库就已足够；NetworkX（见{doc}`09-01-code-setup-solvers`）只用到一次，作为检验。构建本书时不运行这些代码：请把各代码块按顺序复制到脚本或笔记本中，并将打印结果与 `# Expected:` 注释相比较。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 把 DP 模型编码为表示状态、决策、状态转移、费用和边界值的 Python 函数，并用通用的逆序或顺序逐阶段求解器求解；
- 用 `functools.lru_cache` 把同一递推关系写成记忆化的递归函数；
- 从策略表中恢复一个或全部最优决策序列；
- 用穷举法检验 DP 的解，并把实测运算量与预测值相比较；
- 计算无圈网络的拓扑序及其中的最短路。
:::

(sec-09-05-setup)=
## 环境准备

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
## 逐阶段 DP 求解器

逆序递推{eq}`eq-05-01-recurrence`对每个阶段 $k$（从最后一个到第一个）和每个状态 $y\in\mathcal Y_k$ 计算

$$
V_k(y)=\min_{d\in\mathcal D_k(y)}\bigl\{c_k(y,d)+V_{k+1}\bigl(t_k(y,d)\bigr)\bigr\},
$$

其起点是边界值 $V_{N}(y)=g(y)$。下面的求解器以 Python 函数的形式接受模型的各个要素，这样同一段代码就能用于本章的每个例子：

- `states(k)`：阶段 $k$ 的状态集合 $\mathcal Y_k$；
- `decisions(k, y)`：允许的决策 $\mathcal D_k(y)$；
- `transition(k, y, d)`：下一个状态 $t_k(y,d)$；
- `cost(k, y, d)`：阶段费用 $c_k(y,d)$（最大化时为阶段收益）；
- `terminal(y)`：最后一个阶段之后的边界值 $g(y)$。

策略表为每个状态保存**全部**最优决策，因为并列会产生多重最优解（alternative optima），而例子中正需要这些解。每次计算 $c_k(y,d)+V_{k+1}(\cdot)$ 计为一次加法，从 $m$ 个选项中选出最好的计为 $m-1$ 次比较，这与{doc}`05-01-dp-fundamentals`中计算工作量时所用的约定相同。

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

顺序递推写成“推送”形式：不是追问到达某个状态的最好路从何而来，而是让阶段 $k$ 的每个状态把自己的值传递给它能到达的各个状态。结果是与{eq}`eq-05-01-lattice-forward`中相同的函数 $F_k(y)$，并为每个状态附带一个最好前驱的列表。

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
## 菱形阶段网络

### 数据

{prf:ref}`ex-05-01-diamond`中的网络以{numref}`fig-05-01-diamond-network`中的顶点名称及其格点坐标 $(k,y)$（第 $k$ 列、第 $y$ 层）存储。决策“上”（U）从 $(k,y)$ 通向 $(k+1,y+1)$，“下”（D）通向 $(k+1,y-1)$；`arc_time` 返回 $c^{\uparrow}(k,y)$ 或 $c^{\downarrow}(k,y)$，若弧不存在则返回 $\infty$。

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

### 逆序递推与顺序递推

有了这些函数，{ref}`sec-05-01-backward`中的模型只需一次调用：阶段 0 到 5，边界 $V_6(0)=0$，所求值为 $V_0(0)$。

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

最少时间为 13 分钟，路线为 $A$–$C$–$F$–$J$–$M$–$P$–$B$；$G$ 处的并列在策略表中表现为 `['U', 'D']`；运算次数与{prf:ref}`ex-05-01-backward`相同：24 次加法（每条弧一次）和 9 次比较（每个有两条出弧的顶点一次）。顺序递推给出 $F_6(0)=13$，并且恰在最优路上有 $F_k(y)+V_k(y)=13$。

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

$N$ 处的并列（最好前驱为 $J$ 或 $K$）就是{numref}`tab-05-01-forward`中指出的那个并列。

### 记忆化递归

同一个递推关系{eq}`eq-05-01-lattice-backward`可以写成递归函数。如果不加防范，递归会对到达某个顶点（例如 $F$）的每一条路都重新求解该顶点的子问题，这又回到了完全枚举。装饰器 `lru_cache` 在每个值第一次算出时将其存储，以后调用时直接返回所存的值，因此每个顶点只计算一次：记忆化把递归变成了 DP。现在计算顺序由 Python 决定（从所求值出发深度优先），并且只计算从所求值可以到达的状态。

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

缓存中保存了 16 个值，每个顶点一个；有 9 次调用直接从缓存得到答案，而不必重新计算。

### 穷举检验

从 $A$ 到 $B$ 的所有路，就是由三个 U 和三个 D 组成且始终位于菱形之内的字符串；共有 $\binom63=20$ 条。

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

### 状态扩充：转弯惩罚

在{prf:ref}`ex-05-01-turn-penalty`中，每改变一次方向要花费 3 分钟，状态必须记录离开某个顶点的方向 $\omega$。现在阶段 $k$ 的状态是使相应弧存在的二元组 $(y,\omega)$；决策是离开*下一个*顶点的方向 $d$，其费用为弧的时间，若 $d\ne\omega$ 则再加 3，如{eq}`eq-05-01-turn`所示。边界在阶段 5：$V_5(y,\omega)$ 是最后一条弧的时间。

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

从 $A$ 向下出发得到 17 分钟，向上出发则为 20 分钟；最优路 DDDUUU 为 $A$–$D$–$G$–$K$–$N$–$Q$–$B$，行驶时间 14 分钟，转弯一次。元素 $V_3(-1,\mathrm D)=10$（U）和 $V_4(0,\mathrm U)=7$（D）与该例中的计算相同，穷举也证实了最优解唯一，其次是两条 20 分钟的路。扩充后的模型每条弧一个状态，而不是每个顶点一个状态：阶段 0 到 4 共有 22 个状态，而不是 15 个顶点；需要 34 次加法和 12 次比较，而不是 24 次和 9 次。

### 运算量：实测与预测

对于有 $N$ 个阶段的菱形格点网络（$N$ 为偶数，$m=N/2$），{eq}`eq-05-01-count-dp`预测 DP 需要 $2m(m+1)$ 次加法和 $m^2$ 次比较，{eq}`eq-05-01-count-enumeration`预测完全枚举需要 $(N-1)\binom{N}{N/2}$ 次加法。下面的函数构建任意偶数规模、弧时间随机的格点网络，并在其上运行同一个求解器。

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

对每个 $N$，实测次数都等于{numref}`tab-05-01-counts`中的预测值，而枚举的加法次数从 100 增长到约 $4.5\times10^9$。

(sec-09-05-replacement)=
## 设备更新

### 以年份和役龄为状态的形式

{prf:ref}`ex-05-02-replacement`中的数据以千美元计，按役龄（age）$i$ 给出：运行费用 $c(i)$、旧设备折价 $t(i)$ 和残值 $s(i)$；新机器的价格为 $p=50$，计划期为 $N=5$ 年，机器当前役龄为 $i_0=2$ 年。

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

带边界条件{eq}`eq-05-02-boundary`即 $V_6(i)=-s(i)$ 的递推式{eq}`eq-05-02-recurrence`，同样只需一次调用。

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

最小费用为 11.5 万美元，对应该例中的两个方案 $\mathrm{BBKBK}$ 和 $\mathrm{BKBBK}$；$V_5(3)=24$（购买）以及并列 $V_2(1)=76$ 与{numref}`tab-05-02-two-arg`中的计算相同。工作量为 $N(N+1)=30$ 次加法和 $N(N+1)/2=15$ 次比较（见{ref}`sec-05-02-effort`）。穷举法计算全部 $2^5=32$ 个方案。

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

### 单变量形式

在形式{eq}`eq-05-02-one-arg`中，值 $U_k$ 只依赖于年份：企业在第 $k$ 年年初持有一台役龄为一年的机器，决定将其保留多久。该递推式可以向前跳到以后的任何一年，因此自然地写成记忆化的递归函数。对于这个小例子，甚至不需要运行费用的累加和。

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

值 $U_6,\dots,U_2=-25,-4,24,48,76$ 和 $V_1(2)=115$ 与{numref}`tab-05-02-one-arg`相同，并且按定义对每个 $k$ 都有 $U_k=V_k(1)$。

(sec-09-05-allocation)=
## 资源分配

### 单资源

在{prf:ref}`ex-05-03-one-resource`中，把八个单位的资源分配给四项活动，收益见{numref}`tab-05-03-returns`。活动是阶段，剩余单位数 $y$ 是状态，分配量 $x\le y$ 是决策。由于约束是等式，边界为 $V_4(y)=r_4(y)$（{eq}`eq-05-03-boundary`），因此决策只在阶段 1 到 3 中进行，并且阶段 1 只需要状态 $y=8$。

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

该方案给出 $x_1=2$、$x_2=6$、$x_3=0$ 和 $x_4=8-2-6-0=0$：收益为 24，且唯一。99 次加法与{ref}`sec-05-03-one-resource`中的计数一致：两个完整阶段各有 $(X+1)(X+2)/2=45$ 个选项，阶段 1 有 9 个选项。穷举法遍历恰好分配八个单位的全部 $\binom{11}{3}=165$ 种分配方案。

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

### 双资源

有两种资源时，状态是剩余单位数构成的二元组 $(y_1,y_2)$（{eq}`eq-05-03-two-resources`）。{prf:ref}`ex-05-03-two-resources`的数据是三张 $4\times4$ 的表 $r_k(x_1,x_2)$，行对应 $x_1$，列对应 $x_2$。

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

活动 1 得到 $(1,0)$，随后活动 2 从 $(2,3)$ 中得到 $(2,0)$，活动 3 得到其余的 $(0,3)$：收益为 $4+4+8=16$。仅阶段 2 就需要 $100$ 次加法，即计数{eq}`eq-05-03-two-count`；若只有一种三个单位的资源，则只需 10 次。值 $V_2(1,1)=5$ 和 $V_2(3,1)=9$ 与该例一致。

### 背包问题

0–1 背包问题（0–1 knapsack problem）{eq}`eq-05-03-knapsack`简短到可以写成记忆化递归。在{prf:ref}`ex-08-07-a8-q4`中，五份住房申请的需求评分为 $h_k$，各需要 $a_k$ 名志愿者，而可用的志愿者至多 23 名。`best_score(k, y)` 就是 $V_k(y)$，即剩余 $y$ 名志愿者时申请 $k,\dots,5$ 的最高总评分。

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

批准申请 2、3、4 和 5，总评分为 279，用完全部 23 名志愿者。递归只计算了从 $(1,23)$ 可以到达的状态，远少于完整表格的 $5\times24$ 个单元。

### 资本预算

在{prf:ref}`ex-05-03-capital-budgeting`中，四项活动从 6 美元的预算中购买两种资源，单价分别为 3 美元和 1 美元。收益表 $r_k(x_1,x_2)$ 只对买得起的组合给出；下面列表中的一行存放 $r_k(x_1,0),r_k(x_1,1),\dots$。

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

**形式 2**（{eq}`eq-05-03-budget-2`）先计算每项活动用至多 $\beta$ 美元所能得到的最大收益 $\hat r_k(\beta)$（{eq}`eq-05-03-hat-r`）以及达到该收益的组合，再把资金作为单一资源进行分配，边界为 $V_5(b)=0$（允许有未花完的资金）。

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

最大收益为 19，五个方案就是{prf:ref}`ex-05-03-capital-budgeting`表格中的方案。**形式 1**（{eq}`eq-05-03-budget-1`）把两种购买量都作为一个阶段的决策；它给出相同的值函数（value function），但要对每一个剩余资金额计算每一个买得起的组合。

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

形式 1 每个阶段进行 39 次计算。{numref}`tab-05-03-budget-effort`中的计算次数只依赖于价格和预算，因此对于更大的预算，无需数据就能算出：形式 1 对每个 $b\le B$ 计算满足 $3x_1+x_2\le b$ 的全部组合；形式 2 对满足 $3x_1+x_2\le B$ 的每个组合计算一次，然后计算资金的 $(B+1)(B+2)/2$ 种划分。

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
## 无圈网络中的最短路

### 拓扑序

{prf:ref}`thm-05-04-topological-order`给无圈网络（acyclic network）的节点编号，使每条弧都指向前方：反复取出一个没有入弧的节点（源节点，见{prf:ref}`lem-05-04-source-node`），并将它连同其弧一起删除。该函数记录每个节点的入弧数以及当前源节点构成的队列；若在所有节点编号之前队列已空，则剩下的节点中含有圈。

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

### 单遍递推

按拓扑序，顺序递推{eq}`eq-05-04-acyclic-forward`只需一遍：对进入 $j$ 的弧取 $V(j)=\min_i\{V(i)+c_{ij}\}$，起点的 $V$ 等于 0。逆序递推{eq}`eq-05-04-acyclic-backward`按相反的顺序遍历。两者都统计各自的加法和比较次数。

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

对于{prf:ref}`ex-05-04-acyclic`中的网络（长度以千米计），拓扑序为 $D,C,A,E,F,B$，从 $D$ 到 $B$ 的最短路为 $D$–$C$–$A$–$E$–$F$–$B$，长 8 千米，递推使用了 9 次加法和 4 次比较。回溯使用下面的 `trace`，与{doc}`09-04-code-networks`中一样，沿前驱节点反向读出。

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

每个节点都位于 $D$–$B$ 最短路上，因此处处有 $V(j)+\bar V(j)=8$，这正是该例中的检验。NetworkX 求得相同的顺序和相同的路。

```python
G_acyclic = nx.DiGraph()
G_acyclic.add_weighted_edges_from(acyclic_arcs)
print(
    list(nx.topological_sort(G_acyclic)),
    nx.shortest_path(G_acyclic, "D", "B", weight="weight"),
)
# Expected: ['D', 'C', 'A', 'E', 'F', 'B'] ['D', 'C', 'A', 'E', 'F', 'B']
```

### 作为无圈网络最短路的设备更新问题

{ref}`sec-05-02-network`把设备更新问题画成一个网络：每个（年份，役龄）二元组对应一个节点，有保留弧和购买弧，以及长度为 $-s(i)$ 的终端弧。该网络是无圈的，因此尽管终端弧的长度为负，单遍递推仍然适用，而到终端节点的距离就是{ref}`sec-09-05-replacement`中的最小费用。

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

该网络有 21 个节点和 30 条决策弧，另有终端节点及其六条弧；最短路的费用为 115，对应两个最优方案之一。

### 单遍递推的运算量

在每个节点都有弧通向每个编号更大的节点的网络中，{eq}`eq-05-04-acyclic-count`预测需要 $n(n-1)/2$ 次加法和 $(n-1)(n-2)/2$ 次比较。

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

## 小结

- DP 模型就是五个函数（状态、决策、状态转移、阶段费用、边界值）；一个求解{eq}`eq-05-01-recurrence`的逆序求解器和一个顺序求解器为每个例子给出 $V$、$P$（全部最优决策）和运算次数。
- 使用 `lru_cache` 的记忆化递归对每个状态只计算一次；它适用于决策跨度较大的递推关系（{eq}`eq-05-02-one-arg`），以及只有少数状态可达的问题（背包问题）。
- 菱形网络给出 13 分钟，需要 24 次加法和 9 次比较；转弯惩罚模型给出 17 分钟；对所有格点规模，实测次数都等于{eq}`eq-05-01-count-dp`。
- 设备更新问题在两种形式下以及作为无圈网络最短路时，都给出 11.5 万美元和方案 $\mathrm{BBKBK}$、$\mathrm{BKBBK}$。
- 资源分配给出 24（单资源）、16（双资源，阶段 2 有 100 个选项）、279（背包问题）和带五个方案的 19（资本预算）；{numref}`tab-05-03-budget-effort`中的计算次数可通过统计买得起的组合得到。
- 通过删除源节点得到的拓扑序给出{doc}`05-04-dp-shortest-paths`中的单遍最短路递推；每个结果都经穷举法或 NetworkX 证实。

## 参考文献

- Bellman, R. (1957). *Dynamic programming*. Princeton University Press.
- Dreyfus, S. E., & Law, A. M. (1977). *The art and theory of dynamic
  programming*. Academic Press. (Chapters 1–4: staged networks, equipment
  replacement, resource allocation, shortest paths.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 11: dynamic programming.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 12: deterministic dynamic programming.)
