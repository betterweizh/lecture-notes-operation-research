# 用 Python 实现网络算法

## 概述

本章用 Python 实现{doc}`04-01-networks-mst`、{doc}`04-02-shortest-paths`、{doc}`04-03-maximum-flow`、{doc}`04-04-minimum-cost-flow`和{doc}`05-04-dp-shortest-paths`中的网络算法。首先介绍存储网络的数据结构（弧表（arc list）、节点–弧关联矩阵（node–arc incidence matrix）、长度矩阵（length matrix）、邻接表（adjacency list））及其在 NetworkX 中的对应物。随后按照本书的步骤和记号逐一写出各个算法：带标号的 Prim 算法和使用并查集（union–find）结构的 Kruskal 算法；数组版本和堆版本的 Dijkstra 算法，包括前驱标号 $P(j)$ 和路的回溯；以中间节点为外层循环的 Floyd–Warshall 算法；带负圈（negative cycle）检测的 Bellman–Ford 方法；使用标号法（labeling procedure）和广度优先搜索（breadth-first search，即 Edmonds–Karp 规则）的 Ford–Fulkerson 方法，以及最小割（minimum cut）；最后是求解最小费用流的连续最短路算法（successive shortest paths）。每个实现都在理论各章的例题上运行，并与 NetworkX 的结果核对；对于三个网络线性规划问题（linear program, LP），还与使用关联矩阵的 `scipy.optimize.linprog` 核对，其用法与{doc}`09-02-code-lp-modeling`相同。构建本书时不运行这些代码：请把各代码块按顺序复制到脚本或笔记本中，并将打印结果与 `# Expected:` 注释相比较。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 把网络存储为弧表、关联矩阵、长度矩阵或邻接表，并在这些形式与 NetworkX 图之间相互转换；
- 实现 Prim 算法和 Kruskal 算法、Dijkstra 算法（数组版本和堆版本）、Floyd–Warshall 算法和 Bellman–Ford 算法，并由前驱标号回溯出路；
- 用标号法或广度优先搜索实现 Ford–Fulkerson 方法，由残余容量（residual capacity）恢复弧流量，并读出最小割；
- 用连续最短路算法求解带流量下界的最小费用流问题（minimum-cost flow problem）；
- 用 NetworkX 以及由节点–弧关联矩阵构建的线性规划问题检验每个结果。
:::

(sec-09-04-setup)=
## 环境准备

本章所有代码块共享同一个命名空间。第一个代码块导入后面用到的全部内容。节点就是本书例题中的整数（或名称），这样打印出的标号可以直接与理论各章的表格相比较。`INF` 代表 $+\infty$：表示不存在的弧、尚未设定的标号或无界的容量。

```python
import heapq
import itertools
import math
from collections import deque

import networkx as nx
import numpy as np
from scipy.optimize import linprog

np.set_printoptions(precision=3, suppress=True, linewidth=120)
INF = math.inf
```

(sec-09-04-representations)=
## 网络的表示

这里为{prf:ref}`ex-04-01-road-network`中的有向道路网络构建{doc}`04-01-networks-mst`中的四种数据结构。由三元组 $(i,j,c_{ij})$ 组成的**弧表**是本章所有函数的输入格式；其他结构都由它导出。

```python
# Directed road network: arcs (i, j, c_ij) with lengths in miles
road_nodes = [1, 2, 3, 4, 5]
road_arcs = [
    (1, 2, 0.5),
    (1, 4, 2.5),
    (1, 5, 5.0),
    (2, 4, 1.7),
    (2, 5, 1.2),
    (3, 1, 1.0),
    (3, 4, 2.0),
    (4, 5, 3.0),
]
```

**节点–弧关联矩阵**每个节点一行、每条弧一列；弧 $(i,j)$ 的列为 $\ve_i-\ve_j$，在弧尾处为 $+1$，在弧头处为 $-1$。它是下文每个网络流 LP 的约束矩阵。

```python
def incidence_matrix(nodes: list, arcs: list) -> np.ndarray:
    """Node-arc incidence matrix of a directed network.

    Parameters
    ----------
    nodes : list
        Node labels; row ``k`` belongs to ``nodes[k]``.
    arcs : list of tuple
        Arcs ``(i, j, ...)``; column ``a`` belongs to ``arcs[a]``. Only the
        first two entries (tail and head) are used.

    Returns
    -------
    numpy.ndarray
        Matrix with +1 in the row of the tail and -1 in the row of the head
        of every arc.
    """
    row = {node: k for k, node in enumerate(nodes)}
    A = np.zeros((len(nodes), len(arcs)))
    for a, (i, j, *_) in enumerate(arcs):
        if i not in row or j not in row:
            raise ValueError(f"arc ({i}, {j}) has an unknown end node")
        if i == j:
            raise ValueError(f"loop at node {i} is not allowed")
        A[row[i], a] = 1.0
        A[row[j], a] = -1.0
    return A


A_road = incidence_matrix(road_nodes, road_arcs)
print(A_road)
# Expected:
#   [[ 1.  1.  1.  0.  0. -1.  0.  0.]
#    [-1.  0.  0.  1.  1.  0.  0.  0.]
#    [ 0.  0.  0.  0.  0.  1.  1.  0.]
#    [ 0. -1.  0. -1.  0.  0. -1.  1.]
#    [ 0.  0. -1.  0. -1.  0.  0. -1.]]
print(A_road.sum(axis=0), np.linalg.matrix_rank(A_road))
# Expected: [0. 0. 0. 0. 0. 0. 0. 0.] 4
```

该矩阵与{numref}`tab-04-01-incidence`一致。每一列之和都为零，因此 $\vone^\T\mA=\vzero^\T$；由于网络是连通的，其秩为 $\abs{\mathcal N}-1=4$。

**长度矩阵** $\mD_0$ 在弧对应的位置上是 $c_{ij}$，对角线上为 0，其余位置为 $+\infty$；它是数组版本的 Dijkstra 算法和 Floyd–Warshall 算法的输入。**邻接表** $\mathcal A(i)$ 为每个节点存储从它出发的弧，适合每次扫描一个节点的弧的算法。对于无向网络，每条边按两个方向各存储一次。

```python
def length_matrix(nodes: list, arcs: list, directed: bool = True):
    """Length matrix D0: c_ij for arcs, 0 on the diagonal, inf otherwise.

    Nodes must be numbered 1, ..., n; row i - 1 belongs to node i.
    """
    n = len(nodes)
    if sorted(nodes) != list(range(1, n + 1)):
        raise ValueError("length_matrix needs nodes numbered 1, ..., n")
    D = np.full((n, n), INF)
    np.fill_diagonal(D, 0.0)
    for i, j, c in arcs:
        D[i - 1, j - 1] = c
        if not directed:
            D[j - 1, i - 1] = c
    return D


def adjacency_lists(arcs: list, directed: bool = True) -> dict:
    """Adjacency lists {i: [(j, c_ij), ...]}, sorted by head node."""
    succ = {}
    for i, j, c in arcs:
        succ.setdefault(i, []).append((j, c))
        succ.setdefault(j, [])
        if not directed:
            succ[j].append((i, c))
    return {i: sorted(succ[i]) for i in sorted(succ)}


print(length_matrix(road_nodes, road_arcs))
# Expected:
#   [[0.  0.5 inf 2.5 5. ]
#    [inf 0.  inf 1.7 1.2]
#    [1.  inf 0.  2.  inf]
#    [inf inf inf 0.  3. ]
#    [inf inf inf inf 0. ]]
print(
    {i: [j for j, _ in lst] for i, lst in adjacency_lists(road_arcs).items()}
)
# Expected: {1: [2, 4, 5], 2: [4, 5], 3: [1, 4], 4: [5], 5: []}
```

邻接表为 $\mathcal A(1)=\{2,4,5\}$，$\mathcal A(2)=\{4,5\}$，$\mathcal A(3)=\{1,4\}$，$\mathcal A(4)=\{5\}$，$\mathcal A(5)=\varnothing$，与{doc}`04-01-networks-mst`相同。NetworkX 的 `DiGraph` 是一个由邻接表构成的字典，弧的属性存放在 `"weight"` 或 `"capacity"` 之类的键下。有两个约定与本书不同，值得核对一次：`nx.incidence_matrix(..., oriented=True)` 在弧尾处放 $-1$、在弧头处放 $+1$，符号与本书相反，并且返回的是稀疏矩阵。

```python
G_road = nx.DiGraph()
G_road.add_weighted_edges_from(road_arcs)
B = nx.incidence_matrix(
    G_road,
    nodelist=road_nodes,
    edgelist=[a[:2] for a in road_arcs],
    oriented=True,
).toarray()
print(np.array_equal(B, -A_road))
# Expected: True
print(nx.to_numpy_array(G_road, nodelist=road_nodes, weight=None))
# Expected:
#   [[0. 1. 0. 1. 1.]
#    [0. 0. 0. 1. 1.]
#    [1. 0. 0. 1. 0.]
#    [0. 0. 0. 0. 1.]
#    [0. 0. 0. 0. 0.]]
print([G_road.out_degree(i) for i in road_nodes])
# Expected: [3, 2, 2, 1, 0]
```

0–1 邻接矩阵和出度 3、2、2、1、0 与{prf:ref}`ex-04-01-road-network`一致。

(sec-09-04-mst)=
## 最小生成树

{prf:ref}`ex-04-01-sewage`中的污水管网是无向网络；每条连接以 $(i,j,c_{ij})$（$i<j$）的形式存储一次，费用以百万美元计。

```python
sewage_nodes = [1, 2, 3, 4, 5, 6, 7]
sewage_edges = [
    (1, 2, 1),
    (1, 3, 5),
    (1, 4, 7),
    (1, 5, 9),
    (2, 3, 6),
    (2, 4, 4),
    (2, 5, 3),
    (3, 4, 5),
    (3, 6, 10),
    (3, 7, 8),
    (4, 5, 8),
    (4, 6, 3),
    (5, 7, 3),
    (6, 7, 4),
]
```

### 带标号的 Prim 算法

{prf:ref}`alg-04-01-prim`从起始节点 $r$ 出发生长一棵树。{prf:ref}`thm-04-01-prim`之后所述的 $O(\abs{\mathcal N}^2)$ 版本为每个尚未连接的节点 $j$ 保存标号 $[\gamma(j),P(j)]$：即从 $j$ 到已连接集合 $C_k$ 的最便宜边的费用，以及该边在 $C_k$ 中的另一端节点。下一个加入的是 $\gamma$ 最小的节点（并列时取编号较小者），并且只需检查新加入节点处的边。标号只在严格改进时才改变，因此在并列的边中保留最先找到的那条；按编号递增的顺序扫描节点，就重现了{prf:ref}`ex-04-01-sewage-prim`中的并列规则。

```python
def prim(nodes: list, edges: list, root) -> list:
    """Minimum spanning tree by Prim's algorithm with labels [gamma, P].

    Parameters
    ----------
    nodes : list
        Node labels (sortable).
    edges : list of tuple
        Undirected edges ``(i, j, c_ij)``.
    root
        Start node r.

    Returns
    -------
    list of tuple
        Tree edges ``(P(j), j, c)`` in the order in which they are added.
    """
    if root not in nodes:
        raise ValueError(f"start node {root} is not a node")
    succ = adjacency_lists(edges, directed=False)
    gamma = {j: INF for j in nodes}
    P = {j: None for j in nodes}
    connected = {root}
    tree = []
    i = root
    while len(connected) < len(nodes):
        # update the labels of the unconnected neighbors of the new node i
        for j, c in succ.get(i, []):
            if j not in connected and c < gamma[j]:
                gamma[j], P[j] = c, i
        candidates = [j for j in nodes if j not in connected]
        j_star = min(candidates, key=lambda j: (gamma[j], j))
        if gamma[j_star] == INF:
            raise ValueError("the network is not connected")
        tree.append((P[j_star], j_star, gamma[j_star]))
        connected.add(j_star)
        i = j_star
    return tree


tree_prim = prim(sewage_nodes, sewage_edges, root=7)
print(tree_prim)
# Expected: [(7, 5, 3), (5, 2, 3), (2, 1, 1), (2, 4, 4), (4, 6, 3), (1, 3, 5)]
print(sum(c for *_, c in tree_prim))
# Expected: 19
```

各连接按{numref}`tab-04-01-prim`中的顺序 5–7、2–5、1–2、2–4、4–6、1–3 加入，树的费用为 1 900 万美元。

### 使用并查集的 Kruskal 算法

{prf:ref}`alg-04-01-kruskal`需要判断一条边的两个端点是否位于当前已构建的森林（forest）的同一个连通分支中。**并查集**结构可以回答这个问题：每个连通分支是一棵由父指针构成的树，其根节点为该分支命名。`find` 沿指针走到根（并在途中把路径减半，使以后的查找更短）；`union` 把较小分支的根挂到较大分支的根之下。这是对{prf:ref}`thm-04-01-kruskal`之后所分析的重命名方案的改进；两者都使树保持较浅，而采用路径减半后，每次操作的代价几乎是常数。

```python
class UnionFind:
    """Disjoint sets of nodes with union by size and path halving."""

    def __init__(self, nodes: list) -> None:
        self.parent = {i: i for i in nodes}
        self.size = {i: 1 for i in nodes}

    def find(self, i):
        """Return the name (root node) of the component containing i."""
        if i not in self.parent:
            raise ValueError(f"unknown node {i}")
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, i, j) -> bool:
        """Merge the components of i and j; False if they are the same."""
        ri, rj = self.find(i), self.find(j)
        if ri == rj:
            return False
        if self.size[ri] < self.size[rj]:
            ri, rj = rj, ri
        self.parent[rj] = ri
        self.size[ri] += self.size[rj]
        return True


def kruskal(nodes: list, edges: list) -> tuple[list, list]:
    """Minimum spanning tree (forest) by Kruskal's algorithm.

    Edges are examined in order of cost, ties by end nodes. Returns the
    tree edges and a log of ``(edge, accepted)`` decisions up to the last
    accepted edge.
    """
    components = UnionFind(nodes)
    tree, log = [], []
    for i, j, c in sorted(edges, key=lambda e: (e[2], e[0], e[1])):
        if len(tree) == len(nodes) - 1:
            break
        accepted = components.union(i, j)
        log.append(((i, j, c), accepted))
        if accepted:
            tree.append((i, j, c))
    return tree, log


tree_kruskal, log = kruskal(sewage_nodes, sewage_edges)
for (i, j, c), accepted in log:
    print(f"{i}-{j} ({c}): {'accept' if accepted else 'reject'}")
# Expected:
#   1-2 (1): accept
#   2-5 (3): accept
#   4-6 (3): accept
#   5-7 (3): accept
#   2-4 (4): accept
#   6-7 (4): reject
#   1-3 (5): accept
print(sum(c for *_, c in tree_kruskal))
# Expected: 19
```

各项决定与{numref}`tab-04-01-kruskal`相同：6–7 被拒绝，因为它会形成圈 6–7–5–2–4–6；算法在加入 1–3 后停止，不再查看其余七条连接。

### 多重最优解与穷举检验

{prf:ref}`ex-04-01-sewage-kruskal`断言该网络有 1 020 棵生成树，其中恰有四棵的费用为 1 900 万美元。这两个断言都容易用计算机检验。矩阵树定理（matrix-tree theorem）指出，生成树的个数等于拉普拉斯矩阵（Laplacian matrix）$\mL=\diag(\deg)-(\text{邻接矩阵})$ 的任一代数余子式；穷举法则用并查集结构测试全部 $\binom{14}{6}=3003$ 个六条连接的集合（七个节点上的六条连接构成生成树，当且仅当其中没有一条形成圈）。

```python
G_sewage = nx.Graph()
G_sewage.add_weighted_edges_from(sewage_edges)
L = nx.laplacian_matrix(G_sewage, nodelist=sewage_nodes, weight=None)
print(round(np.linalg.det(L.toarray()[1:, 1:])))
# Expected: 1020


def is_spanning_tree(nodes: list, edges: list) -> bool:
    """True if the edges form a spanning tree of the nodes."""
    components = UnionFind(nodes)
    acyclic = all(components.union(i, j) for i, j, _ in edges)
    return acyclic and len(edges) == len(nodes) - 1


trees = [
    T
    for T in itertools.combinations(sewage_edges, len(sewage_nodes) - 1)
    if is_spanning_tree(sewage_nodes, T)
]
costs = [sum(c for *_, c in T) for T in trees]
optimal = [T for T, z in zip(trees, costs) if z == min(costs)]
print(len(trees), min(costs), len(optimal))
# Expected: 1020 19 4
for T in optimal:
    print([f"{i}-{j}" for i, j, _ in T])
# Expected:
#   ['1-2', '1-3', '2-4', '2-5', '4-6', '5-7']
#   ['1-2', '1-3', '2-5', '4-6', '5-7', '6-7']
#   ['1-2', '2-4', '2-5', '3-4', '4-6', '5-7']
#   ['1-2', '2-5', '3-4', '4-6', '5-7', '6-7']
```

这就是{prf:ref}`ex-04-01-sewage-kruskal`中的四种布局。NetworkX 在同一个函数背后实现了 Kruskal 算法、Prim 算法和 Borůvka 算法。三者都返回权重为 19 的树；返回四棵最优树中的哪一棵，取决于内部的并列处理方式。

```python
for method in ["kruskal", "prim", "boruvka"]:
    T = nx.minimum_spanning_tree(G_sewage, algorithm=method)
    print(method, T.size(weight="weight"), sorted(T.edges()))
# Expected:
#   kruskal 19.0 [(1, 2), (1, 3), (2, 4), (2, 5), (4, 6), (5, 7)]
#   prim 19.0 [(1, 2), (1, 3), (2, 4), (2, 5), (4, 6), (5, 7)]
#   boruvka 19.0 [(1, 2), (1, 3), (2, 4), (2, 5), (4, 6), (5, 7)]
```

(sec-09-04-dijkstra)=
## Dijkstra 算法

{prf:ref}`ex-04-02-route`中的网络有七个节点、十二条弧，长度以分钟计。

```python
route_nodes = [1, 2, 3, 4, 5, 6, 7]
route_arcs = [
    (1, 2, 2),
    (1, 4, 4),
    (2, 5, 5),
    (3, 1, 10),
    (3, 2, 11),
    (3, 6, 7),
    (3, 7, 5),
    (4, 3, 3),
    (4, 6, 1),
    (5, 3, 8),
    (5, 7, 6),
    (6, 7, 9),
]
D_route = length_matrix(route_nodes, route_arcs)
```

### 基于长度矩阵的数组版本

第一个实现在长度矩阵上按照{prf:ref}`alg-04-02-dijkstra`进行。每次迭代扫描最近获得永久标号的节点 $i$，改进尚无永久标号的节点的标号 $[V(j),P(j)]$，然后把最小的临时标号（temporary label）变为永久标号（permanent label）（并列时取节点编号较小者）。临时标号节点是 $V(j)$ 有限但尚未获得永久标号的节点。寻找最小标号时要扫描全部节点，每次迭代为 $O(\abs{\mathcal N})$，总计为 $O(\abs{\mathcal N}^2)$。该函数还返回各节点获得永久标号的顺序。

```python
def dijkstra_matrix(D: np.ndarray, s: int, t: int | None = None):
    """Dijkstra's algorithm on a length matrix (array version).

    Parameters
    ----------
    D : numpy.ndarray
        Length matrix; entry ``D[i-1, j-1]`` is c_ij (inf if no arc).
    s : int
        Origin node (nodes are numbered 1, ..., n).
    t : int, optional
        Destination; the algorithm stops when t becomes permanent.

    Returns
    -------
    V, P : dict
        Labels V(j) and predecessors P(j) (None for s and unreached nodes).
    order : list
        Nodes in the order in which their labels became permanent.
    """
    n = D.shape[0]
    if D.shape != (n, n) or not 1 <= s <= n:
        raise ValueError("D must be square and s a node 1, ..., n")
    if np.any(D < 0):
        raise ValueError("Dijkstra's algorithm needs lengths c_ij >= 0")
    V = {j: INF for j in range(1, n + 1)}
    P = {j: None for j in range(1, n + 1)}
    V[s] = 0.0
    permanent, order, i = {s}, [s], s
    while i != t and len(permanent) < n:
        # scan node i: try to improve the labels of non-permanent nodes
        for j in range(1, n + 1):
            if j not in permanent and V[i] + D[i - 1, j - 1] < V[j]:
                V[j], P[j] = V[i] + float(D[i - 1, j - 1]), i
        temporary = [j for j in V if j not in permanent and V[j] < INF]
        if not temporary:
            break
        i = min(temporary, key=lambda j: (V[j], j))
        permanent.add(i)
        order.append(i)
    return V, P, order


def trace_path(P: dict, s, j) -> list:
    """Path from s to j read backwards from predecessor labels P."""
    path = [j]
    while path[-1] != s:
        # P may also lack unreached nodes (dijkstra_heap).
        if P.get(path[-1]) is None:
            raise ValueError(f"node {j} is not reached from {s}")
        path.append(P[path[-1]])
    return path[::-1]


V, P, order = dijkstra_matrix(D_route, s=1, t=7)
print(order)
# Expected: [1, 2, 4, 6, 3, 5, 7]
print(V)
# Expected: {1: 0.0, 2: 2.0, 3: 7.0, 4: 4.0, 5: 7.0, 6: 5.0, 7: 12.0}
print(trace_path(P, 1, 7), V[7])
# Expected: [1, 4, 3, 7] 12.0
```

各节点按{numref}`tab-04-02-dijkstra`中的顺序 1、2、4、6、3、5、7 获得永久标号，长度为 12 分钟的路 $1\to4\to3\to7$ 由 $P(7)=3$、$P(3)=4$、$P(4)=1$ 反向读出。由于节点 7 是最远的节点，运行结束时所有标号都已是最终值。

### 基于邻接表的堆版本

对于稀疏网络，最小的临时标号保存在二叉堆（binary heap，Python 的 `heapq`）中。当某个标号得到改进时，把新的二元组 $(V(j),j)$ 压入堆中，旧的二元组留在堆里；当它被弹出时，由于其节点那时已获得永久标号，可以识别出来并跳过。这种“惰性删除”使代码简短，运行时间为 $O(\abs{\mathcal A}\log\abs{\mathcal N})$。并列规则自动满足：$V$ 相等的二元组按节点编号的顺序弹出。

```python
def dijkstra_heap(succ: dict, s):
    """Dijkstra's algorithm with a binary heap on adjacency lists.

    Parameters
    ----------
    succ : dict
        Adjacency lists ``{i: [(j, c_ij), ...]}``.
    s
        Origin node.

    Returns
    -------
    V, P : dict
        Distances and predecessors of the nodes reachable from s.
    """
    if s not in succ:
        raise ValueError(f"origin {s} is not a node")
    if any(c < 0 for arcs in succ.values() for _, c in arcs):
        raise ValueError("Dijkstra's algorithm needs lengths c_ij >= 0")
    V, P = {s: 0}, {s: None}
    heap = [(0, s)]
    permanent = set()
    while heap:
        v, i = heapq.heappop(heap)
        if i in permanent:
            continue
        permanent.add(i)
        for j, c in succ[i]:
            if j not in permanent and v + c < V.get(j, INF):
                V[j], P[j] = v + c, i
                heapq.heappush(heap, (V[j], j))
    return V, P


succ_route = adjacency_lists(route_arcs)
V_heap, P_heap = dijkstra_heap(succ_route, 1)
print(dict(sorted(V_heap.items())))
# Expected: {1: 0, 2: 2, 3: 7, 4: 4, 5: 7, 6: 5, 7: 12}
print(sorted((P_heap[j], j) for j in P_heap if P_heap[j] is not None))
# Expected: [(1, 2), (1, 4), (2, 5), (3, 7), (4, 3), (4, 6)]
```

弧 $(P(j),j)$ 构成{numref}`fig-04-02-dijkstra-tree`中的最短路树（shortest-path tree）。NetworkX 返回相同的距离；它的前驱字典列出最短路上的*全部*前驱节点，在这里每个节点只有一个。

```python
G_route = nx.DiGraph()
G_route.add_weighted_edges_from(route_arcs)
pred, dist = nx.dijkstra_predecessor_and_distance(G_route, 1)
print(dist == V_heap, pred[7], nx.dijkstra_path(G_route, 1, 7))
# Expected: True [3] [1, 4, 3, 7]
```

### 用约化长度进行最优性检验

{prf:ref}`thm-04-02-optimality`给出了一个与标号的求法无关的证书：若每个约化长度（reduced length）$\bar c_{ij}=c_{ij}+V(i)-V(j)$ 都非负，则这些标号就是距离。最短路树上的弧满足 $\bar c_{ij}=0$。

```python
reduced = {(i, j): c + V_heap[i] - V_heap[j] for i, j, c in route_arcs}
print(min(reduced.values()))
# Expected: 0
print([a for a, r in reduced.items() if r == 0])
# Expected: [(1, 2), (1, 4), (2, 5), (3, 7), (4, 3), (4, 6)]
```

### 最短路 LP

由 LP {eq}`eq-04-02-lp`也能得到同样的答案：从 $s=1$ 到 $t=7$ 输送一个单位的流，以关联矩阵为约束矩阵。HiGHS 在 `res.eqlin.marginals` 中返回各平衡方程的对偶值 $\pi_i$；令 $V(j)=\pi_s-\pi_j$，它们就是{eq}`eq-04-02-dual`中的节点势（node potential）。

```python
A_route = incidence_matrix(route_nodes, route_arcs)
cost_route = np.array([c for *_, c in route_arcs], dtype=float)
b_route = np.zeros(len(route_nodes))
b_route[0], b_route[-1] = 1.0, -1.0
res = linprog(cost_route, A_eq=A_route, b_eq=b_route, bounds=(0, None))
print(res.status, res.fun)
# Expected: 0 12.0
print([route_arcs[a][:2] for a in np.flatnonzero(res.x > 1e-9)])
# Expected: [(1, 4), (3, 7), (4, 3)]
pi = res.eqlin.marginals
print(pi[0] - pi[-1])
# Expected: 12.0
```

该流使用弧 $(1,4)$、$(4,3)$、$(3,7)$，对偶目标函数值 $\pi_1-\pi_7=V(7)$ 等于 12，这正是强对偶性所要求的。（其他节点势并不唯一：不在该路上的节点可以取任何与 $\pi_i-\pi_j\le c_{ij}$ 相容的值。）

### 一条负长度的弧

两个 Dijkstra 函数都拒绝负长度。{prf:ref}`ex-04-02-negative-arc`中的网络说明了原因。给所有长度加上一个常数并不能解决问题：所有长度都增加 2 后，Dijkstra 算法找到路 $1\to2$，但在原网络中 $1\to3\to2$ 更短。{ref}`sec-09-04-bellman-ford`中的 Bellman–Ford 方法给出正确的距离 $V(2)=1$。

```python
neg_arcs = [(1, 2, 2), (1, 3, 3), (3, 2, -2)]
try:
    dijkstra_heap(adjacency_lists(neg_arcs), 1)
except ValueError as err:
    print(err)
# Expected: Dijkstra's algorithm needs lengths c_ij >= 0
shifted = [(i, j, c + 2) for i, j, c in neg_arcs]
V_s, P_s = dijkstra_heap(adjacency_lists(shifted), 1)
print(trace_path(P_s, 1, 2))
# Expected: [1, 2]
```

(sec-09-04-floyd)=
## 所有点对：Floyd–Warshall 算法

{prf:ref}`alg-04-02-floyd`对每个中间节点 $k$ 把整个距离矩阵更新一次。递推式{eq}`eq-04-02-floyd`由完整的矩阵 $\mD_{k-1}$ 构造 $\mD_k$，因此关于 $k$ 的循环必须是**外层**循环。借助 NumPy，关于 $i$ 和 $j$ 的两层内循环变成一次数组运算：`D[:, k, None] + D[None, k, :]` 就是所有和 $d_{ik}+d_{kj}$ 构成的矩阵。第 $k$ 次迭代中第 $k$ 行和第 $k$ 列不变，因此原地更新是安全的。前驱矩阵（predecessor matrix）存储从 $i$ 到 $j$ 的当前已知最好路上 $j$ 之前的节点，用 0 表示“无前驱”（即本书中的短横线）。第 $k$ 次迭代后若出现负的对角元素，就说明存在负圈。

```python
def floyd_warshall(D0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """All-pairs shortest paths by the Floyd-Warshall algorithm.

    Parameters
    ----------
    D0 : numpy.ndarray
        Length matrix (0 on the diagonal, inf for missing arcs); nodes are
        numbered 1, ..., n.

    Returns
    -------
    D, P : numpy.ndarray
        Distances d_ij and predecessors p_ij (0 where there is none).

    Raises
    ------
    ValueError
        If a negative cycle is found (negative diagonal entry).
    """
    n = D0.shape[0]
    if D0.shape != (n, n) or np.any(np.diag(D0) != 0):
        raise ValueError("D0 must be square with a zero diagonal")
    D = D0.astype(float).copy()
    # p_ij = i for every arc (i, j) (also of length 0), 0 elsewhere
    arc = np.isfinite(D) & ~np.eye(n, dtype=bool)
    P = np.where(arc, np.arange(1, n + 1)[:, None], 0)
    for k in range(n):
        # all sums d_ik + d_kj at once; update where they are shorter
        through_k = D[:, k, None] + D[None, k, :]
        better = through_k < D
        D[better] = through_k[better]
        P[better] = np.broadcast_to(P[k, :], (n, n))[better]
        if np.any(np.diag(D) < 0):
            i = int(np.argmax(np.diag(D) < 0)) + 1
            raise ValueError(
                f"negative cycle through node {i} (iteration k={k + 1})"
            )
    return D, P
```

把该函数用于{prf:ref}`ex-04-02-floyd`中的快递网络，返回该例中的矩阵 $\mD_4$ 和 $\mP_4$。

```python
courier_nodes = [1, 2, 3, 4]
courier_arcs = [
    (1, 2, 5),
    (1, 4, -3),
    (2, 3, 3),
    (3, 1, 3),
    (3, 4, 7),
    (4, 3, 2),
]
D4, P4 = floyd_warshall(length_matrix(courier_nodes, courier_arcs))
print(D4)
# Expected:
#   [[ 0.  5. -1. -3.]
#    [ 6.  0.  3.  3.]
#    [ 3.  8.  0.  0.]
#    [ 5. 10.  2.  0.]]
print(P4)
# Expected:
#   [[0 1 4 1]
#    [3 0 2 1]
#    [3 1 0 1]
#    [3 1 4 0]]
```

与 Dijkstra 算法的标号完全一样，路可以在 $\mP$ 的第 $i$ 行中反向读出：第 4 行给出 $2\leftarrow p_{42}=1\leftarrow p_{41}=3\leftarrow p_{43}=4$。

```python
def fw_path(P: np.ndarray, i: int, j: int) -> list:
    """Path from i to j read backwards in row i of the predecessor matrix."""
    if P[i - 1, j - 1] == 0 and i != j:
        raise ValueError(f"no path from {i} to {j}")
    path = [j]
    while path[-1] != i:
        path.append(int(P[i - 1, path[-1] - 1]))
    return path[::-1]


print(fw_path(P4, 4, 2), fw_path(P4, 2, 4), fw_path(P4, 1, 3))
# Expected: [4, 3, 1, 2] [2, 3, 1, 4] [1, 4, 3]
```

这些路线是 $4\to3\to1\to2$（费用 10，即 1 000 美元）、$2\to3\to1\to4$（费用 3）和 $1\to4\to3$（净收益 100 美元）。当 $c_{14}=-6$ 时，圈 $1\to4\to3\to1$ 的费用为 $-1$，与该例相同，检测在第 3 次迭代时发现了它。

```python
courier_neg = [
    (1, 2, 5),
    (1, 4, -6),
    (2, 3, 3),
    (3, 1, 3),
    (3, 4, 7),
    (4, 3, 2),
]
try:
    floyd_warshall(length_matrix(courier_nodes, courier_neg))
except ValueError as err:
    print(err)
# Expected: negative cycle through node 4 (iteration k=3)
```

**为什么循环顺序很重要。** 一种诱人的写法把关于 $k$ 的循环放在最内层：对每个点对 $(i,j)$，把所有中间节点各尝试一次，$d_{ij}\leftarrow\min_k\{d_{ik}+d_{kj}\}$。这样计算 $d_{ij}$ 所用的元素 $d_{ik}$、$d_{kj}$ 可能还远未达到最终值，而该点对再也不会被重新访问。在{prf:ref}`ex-04-02-route`的网络上，这种写法漏掉了从节点 2 出发的两条路线。

```python
def floyd_wrong_order(D0: np.ndarray) -> np.ndarray:
    """Incorrect variant with k innermost (for comparison only)."""
    D = D0.astype(float).copy()
    n = D.shape[0]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                D[i, j] = min(D[i, j], D[i, k] + D[k, j])
    return D


D_right, P_right = floyd_warshall(D_route)
D_wrong = floyd_wrong_order(D_route)
print(np.argwhere(D_right != D_wrong) + 1)
# Expected:
#   [[2 1]
#    [2 4]]
print(D_right[1, 0], D_wrong[1, 0], fw_path(P_right, 2, 1))
# Expected: 23.0 inf [2, 5, 3, 1]
```

长度为 $5+8+10=23$ 的路线 $2\to5\to3\to1$ 是存在的，但当错误的循环处理点对 $(2,1)$ 时，它所需要的元素 $d_{23}$ 和 $d_{53}$ 仍为 $\infty$ 或过大。NetworkX 与正确的版本一致。`floyd_warshall_numpy` 按 `nodelist` 的顺序返回矩阵，`floyd_warshall_predecessor_and_distance` 则返回字典。

```python
G_courier = nx.DiGraph()
G_courier.add_weighted_edges_from(courier_arcs)
D_nx = nx.floyd_warshall_numpy(G_courier, nodelist=courier_nodes)
print(np.array_equal(D_nx, D4))
# Expected: True
pred_nx, dist_nx = nx.floyd_warshall_predecessor_and_distance(G_courier)
print(pred_nx[4][2], dist_nx[4][2])
# Expected: 1 10
print(
    np.array_equal(
        nx.floyd_warshall_numpy(G_route, nodelist=route_nodes), D_right
    )
)
# Expected: True
```

(sec-09-04-bellman-ford)=
## 负长度的弧：Bellman–Ford 算法

{doc}`05-04-dp-shortest-paths`中的过程 1 用递推式{eq}`eq-05-04-bellman-ford`计算 $V^{(k)}(j)$，即从起点到 $j$、至多含 $k$ 条弧的最短途径（walk）的长度。每次迭代只使用上一次迭代的值（下面代码中的副本 `new`），因此第 $k$ 次迭代确实对应“至多 $k$ 条弧”。由{prf:ref}`thm-05-04-successive`，这些值在第 $n$ 次迭代之前收敛，当且仅当不存在可达的负圈；若在第 $n$ 次迭代中这些值仍有变化，则从某个值发生变化的节点沿前驱节点回溯，就会进入一个负圈。

```python
def bellman_ford(nodes: list, arcs: list, s):
    """Shortest paths from s by successive approximation (Bellman-Ford).

    Parameters
    ----------
    nodes : list
        Node labels.
    arcs : list of tuple
        Arcs ``(i, j, c_ij)``; lengths may be negative.
    s
        Origin node.

    Returns
    -------
    V, P : dict
        Distances V(j) and predecessors P(j).
    history : list of dict
        The values V^(k) for k = 0, 1, ... up to convergence.

    Raises
    ------
    ValueError
        If a negative cycle can be reached from s; the message names it.
    """
    if s not in nodes:
        raise ValueError(f"origin {s} is not a node")
    V = {v: INF for v in nodes}
    P = {v: None for v in nodes}
    V[s] = 0
    history = [dict(V)]
    for _ in range(len(nodes)):
        new, newP = dict(V), dict(P)
        for i, j, c in arcs:
            if V[i] + c < new[j]:
                new[j], newP[j] = V[i] + c, i
        history.append(new)
        if new == V:
            return V, P, history
        changed = [v for v in nodes if new[v] != V[v]]
        V, P = new, newP
    # values still change in iteration n: walk back into the cycle
    v = changed[0]
    for _ in range(len(nodes)):
        v = P[v]
    cycle = [v]
    while P[cycle[-1]] != v:
        cycle.append(P[cycle[-1]])
    cycle = [v] + cycle[::-1]
    raise ValueError(f"negative cycle {cycle} reachable from {s}")
```

在{prf:ref}`ex-05-04-bellman-ford`的网络上（长度以百美元计，负值表示利润），`history` 的各行就是{numref}`tab-05-04-bellman-ford`的各行；第 4 次迭代重复第 3 次迭代的结果。

```python
bf_nodes = [1, 2, 3, 4]
bf_arcs = [(1, 2, 2), (1, 3, 5), (2, 3, 6), (3, 2, -4), (2, 4, 3), (3, 4, 1)]
V_bf, P_bf, history = bellman_ford(bf_nodes, bf_arcs, 1)
for k, values in enumerate(history):
    print(k, [values[j] for j in bf_nodes])
# Expected:
#   0 [0, inf, inf, inf]
#   1 [0, 2, 5, inf]
#   2 [0, 1, 5, 5]
#   3 [0, 1, 5, 4]
#   4 [0, 1, 5, 4]
print(P_bf, trace_path(P_bf, 1, 4))
# Expected: {1: None, 2: 3, 3: 1, 4: 2} [1, 3, 2, 4]
```

到节点 4 的最短路是 $1\to3\to2\to4$，长度为 $5-4+3=4$（400 美元），由于它含三条弧，直到第 3 次迭代才被找到。当 $c_{23}=3$ 时，圈 $2\to3\to2$ 的长度为 $-1$，这些值在第 $4=n$ 次迭代中仍有变化，前驱节点引向该圈：

```python
bf_cycle = [(1, 2, 2), (1, 3, 5), (2, 3, 3), (3, 2, -4), (2, 4, 3), (3, 4, 1)]
try:
    bellman_ford(bf_nodes, bf_cycle, 1)
except ValueError as err:
    print(err)
# Expected: negative cycle [2, 3, 2] reachable from 1
```

NetworkX 给出相同的距离；当存在可达的负圈时，它抛出 `NetworkXUnbounded`，并且可以返回这样一个圈。同一函数也解决了{prf:ref}`ex-04-02-negative-arc`中的负长度弧例子。

```python
G_bf = nx.DiGraph()
G_bf.add_weighted_edges_from(bf_arcs)
pred_bf, dist_bf = nx.bellman_ford_predecessor_and_distance(G_bf, 1)
print(dist_bf == V_bf)
# Expected: True
G_cycle = nx.DiGraph()
G_cycle.add_weighted_edges_from(bf_cycle)
try:
    nx.bellman_ford_predecessor_and_distance(G_cycle, 1)
except nx.NetworkXUnbounded as err:
    print("NetworkXUnbounded:", err)
# Expected: NetworkXUnbounded: Negative cycle detected.
print(nx.find_negative_cycle(G_cycle, 1))
# Expected: [3, 2, 3]
print(bellman_ford([1, 2, 3], neg_arcs, 1)[0])
# Expected: {1: 0, 2: 1, 3: 3}
```

(sec-09-04-max-flow)=
## 最大流

### 残余容量

流算法只在残余容量上操作。按照{doc}`04-03-maximum-flow`，每条弧 $(i,j)$ 都配有一条反向弧 $(j,i)$，若该反向弧本不存在，则其容量为 0；零流的残余容量{eq}`eq-04-03-residual`为 $r_{ij}=u_{ij}$。沿一条路输送 $\Delta$ 个单位，会使路上每条弧的 $r_{ij}$ 减少 $\Delta$、$r_{ji}$ 增加 $\Delta$。最后用{eq}`eq-04-03-recover`，即 $x_{ij}=\max\{0,u_{ij}-r_{ij}\}$，恢复弧流量。

```python
def residual_capacities(arcs: list) -> dict:
    """Residual capacities {(i, j): r_ij} of the zero flow.

    Every arc (i, j, u_ij) also gets the reverse pair (j, i), with
    capacity 0 if that arc does not exist.
    """
    r, seen = {}, set()
    for i, j, u in arcs:
        if u < 0:
            raise ValueError(f"negative capacity on arc ({i}, {j})")
        if (i, j) in seen:
            raise ValueError(f"parallel arcs ({i}, {j}): merge them into one")
        seen.add((i, j))
        r[(i, j)] = r.get((i, j), 0) + u
        r.setdefault((j, i), 0)
    return r


def recover_flows(arcs: list, r: dict) -> dict:
    """Arc flows x_ij = max(0, u_ij - r_ij) from residual capacities."""
    return {(i, j): max(0, u - r[(i, j)]) for i, j, u in arcs}
```

### 寻找增广路的两种方法

**标号法**（{prf:ref}`alg-04-03-labelling`）是一种深度优先搜索：从当前节点 $i$ 出发，给残余容量 $r_{ik}$ 最大的未标号、未删除节点 $k$ 标号（并列时取编号最小者）并移到该节点；遇到死路时删除该节点并退回到其前驱节点。**广度优先搜索**则按到达节点的先后顺序给节点标号，因此找到的第一条到 $t$ 的路所含的弧最少；这就是{doc}`04-03-maximum-flow`中讨论该方法的运行时间时提到的 Edmonds–Karp 规则。两个函数都返回路（或 `None`）以及从 $s$ 可到达的节点集合 $S$；当不存在路时，$S$ 就是一个最小割的源点一侧。

```python
def labeling_path(r: dict, s, t):
    """Augmenting path by the labeling procedure (largest residual first).

    Returns ``(path, S)``: the node list s..t or None, and the set of nodes
    labeled or deleted (all nodes reachable from s if path is None).
    """
    succ = {}
    for i, j in r:
        succ.setdefault(i, []).append(j)
    P = {s: None}
    deleted = set()
    i = s
    while True:
        J = [j for j in sorted(succ.get(i, [])) if r[(i, j)] > 0]
        J = [j for j in J if j not in P and j not in deleted]
        if J:
            k = max(J, key=lambda j: r[(i, j)])
            P[k] = i
            if k == t:
                return trace_path(P, s, t), set(P) | deleted
            i = k
        elif i == s:
            return None, set(P) | deleted
        else:
            deleted.add(i)
            i = P.pop(i)


def bfs_path(r: dict, s, t):
    """Augmenting path with the fewest arcs (breadth-first search)."""
    succ = {}
    for i, j in r:
        succ.setdefault(i, []).append(j)
    P = {s: None}
    queue = deque([s])
    while queue:
        i = queue.popleft()
        for j in sorted(succ.get(i, [])):
            if r[(i, j)] > 0 and j not in P:
                P[j] = i
                if j == t:
                    return trace_path(P, s, t), set(P)
                queue.append(j)
    return None, set(P)
```

`max` 返回多个最大值中的第一个，而候选节点已经排序，因此并列时取编号最小的节点，正如算法所要求的。

### Ford–Fulkerson 方法

{prf:ref}`alg-04-03-ford-fulkerson`反复执行“找到一条路，按其残余容量 $\Delta$ 增广”，直到不存在路为止。搜索函数是一个参数，因此同一段代码既可运行标号法版本，也可运行 Edmonds–Karp 版本。

```python
def ford_fulkerson(arcs: list, s, t, search=labeling_path):
    """Maximum flow and minimum cut by the Ford-Fulkerson method.

    Parameters
    ----------
    arcs : list of tuple
        Arcs ``(i, j, u_ij)`` with finite capacities.
    s, t
        Source and sink.
    search : callable
        ``search(r, s, t) -> (path, S)``, e.g. labeling_path or bfs_path.

    Returns
    -------
    v : number
        Maximum flow value.
    x : dict
        Arc flows.
    S : set
        Source side of a minimum cut.
    log : list of tuple
        The augmenting paths with their amounts Delta.
    """
    if s == t:
        raise ValueError("source and sink must differ")
    r = residual_capacities(arcs)
    v, log = 0, []
    while True:
        path, S = search(r, s, t)
        if path is None:
            return v, recover_flows(arcs, r), S, log
        delta = min(r[(i, j)] for i, j in zip(path, path[1:]))
        for i, j in zip(path, path[1:]):
            r[(i, j)] -= delta
            r[(j, i)] += delta
        v += delta
        log.append((path, delta))


def cut_capacity(arcs: list, S: set) -> float:
    """Capacity u[S, S-bar] of the arcs leaving the node set S."""
    return sum(u for i, j, u in arcs if i in S and j not in S)
```

在{prf:ref}`ex-04-03-diamond`的菱形网络上，标号法版本找到了该例中的三条路，包括第 3 次迭代中从节点 2 的回退，并得到割 $S=\{1\}$。

```python
diamond = [(1, 2, 4), (1, 3, 2), (3, 2, 3), (2, 4, 5), (3, 4, 1)]
v, x, S, log = ford_fulkerson(diamond, 1, 4)
print(log)
# Expected: [([1, 2, 4], 4), ([1, 3, 2, 4], 1), ([1, 3, 4], 1)]
print(v, x)
# Expected: 6 {(1, 2): 4, (1, 3): 2, (3, 2): 1, (2, 4): 5, (3, 4): 1}
print(S, cut_capacity(diamond, S))
# Expected: {1} 6
```

流量 $x_{12}=4$、$x_{13}=2$、$x_{32}=1$、$x_{24}=5$、$x_{34}=1$ 以及容量为 6 的最小割都与该例一致。广度优先搜索先取含两条弧的路 $1\to3\to4$，再取含三条弧的路 $1\to3\to2\to4$；流值和各弧流量都相同。

```python
v_bfs, x_bfs, S_bfs, log_bfs = ford_fulkerson(diamond, 1, 4, search=bfs_path)
print(log_bfs)
# Expected: [([1, 2, 4], 4), ([1, 3, 4], 1), ([1, 3, 2, 4], 1)]
print(v_bfs, x_bfs == x, S_bfs)
# Expected: 6 True {1}
```

### 输油管网

{prf:ref}`ex-04-03-pipeline`中的输油管网需要{ref}`sec-04-03-modelling`中的变换：每条双向管道变为两条方向相反的弧，并增加一个超级源点（super-source）0 和一个超级汇点（super-sink）9，分别与各炼油厂和各终端相连。这些弧的容量为无穷大；如本书所建议，任何大于其余弧总容量的数都可以代替 $\infty$，这样所有残余容量都保持有限。

```python
two_way = [(1, 4, 20), (3, 5, 15), (4, 5, 20), (4, 6, 10), (5, 6, 30)]
one_way = [
    (2, 4, 10),
    (2, 5, 20),
    (4, 7, 10),
    (5, 8, 30),
    (6, 7, 50),
    (6, 8, 20),
]
pipes = two_way + [(j, i, u) for i, j, u in two_way] + one_way
big = sum(u for *_, u in pipes) + 1
pipeline = pipes + [
    (0, 1, big),
    (0, 2, big),
    (0, 3, big),
    (7, 9, big),
    (8, 9, big),
]
print(len(pipeline), big)
# Expected: 21 331
v_pipe, x_pipe, S_pipe, log_pipe = ford_fulkerson(
    pipeline, 0, 9, search=bfs_path
)
print(v_pipe, len(log_pipe), sorted(S_pipe), cut_capacity(pipeline, S_pipe))
# Expected: 65 6 [0, 1, 2, 3] 65
```

最大流为 65 千桶/日，最小割是{numref}`fig-04-03-pipeline`中的炼油厂割 $S=\{0,1,2,3\}$，其容量为 $20+10+20+15=65$。

### 用 NetworkX 和 LP 检验

`nx.maximum_flow` 返回流值和弧流量的字典；`nx.minimum_cut` 返回最小割的容量及其两侧。没有 `"capacity"` 属性的弧容量为无穷大，因此在那里超级源点和超级汇点的弧不需要代替值。

```python
def flow_graph(arcs: list, uncapacitated: list = ()) -> nx.DiGraph:
    """DiGraph with arc capacities; arcs in `uncapacitated` get none."""
    G = nx.DiGraph()
    for i, j, u in arcs:
        G.add_edge(i, j, capacity=u)
    G.add_edges_from(uncapacitated)
    return G


G_diamond = flow_graph(diamond)
value, flow_dict = nx.maximum_flow(G_diamond, 1, 4)
cut_value, (S_nx, T_nx) = nx.minimum_cut(G_diamond, 1, 4)
print(value, flow_dict[3], cut_value, S_nx)
# Expected: 6 {2: 1, 4: 1} 6 {1, 2, 3}
G_pipe = flow_graph(pipes, [(0, 1), (0, 2), (0, 3), (7, 9), (8, 9)])
print(
    nx.maximum_flow_value(G_pipe, 0, 9),
    sorted(nx.minimum_cut(G_pipe, 0, 9)[1][0]),
)
# Expected: 65 [0, 1, 2, 3]
```

当存在多个最小割时（菱形网络有两个），NetworkX 可能返回不同的最小割；它们的容量都是 6。最后，LP {eq}`eq-04-03-maxflow-lp`为 $\max\{v:\mA\vx=v(\ve_s-\ve_t),\ \vzero\le\vx\le\vu\}$。以 $(\vx,v)$ 为变量，等式约束为 $[\mA\ \ -(\ve_s-\ve_t)]\,(\vx,v)=\vzero$，`linprog` 最小化 $-v$。

```python
def max_flow_lp(nodes: list, arcs: list, s, t):
    """Maximum flow value and arc flows from the LP with incidence matrix."""
    A = incidence_matrix(nodes, arcs)
    e = np.zeros(len(nodes))
    e[nodes.index(s)], e[nodes.index(t)] = 1.0, -1.0
    A_eq = np.column_stack([A, -e])
    cost = np.zeros(len(arcs) + 1)
    cost[-1] = -1.0
    bounds = [(0, u) for *_, u in arcs] + [(0, None)]
    res = linprog(cost, A_eq=A_eq, b_eq=np.zeros(len(nodes)), bounds=bounds)
    if res.status != 0:
        raise ValueError(res.message)
    return res.x[-1], res.x[:-1]


print(max_flow_lp([1, 2, 3, 4], diamond, 1, 4)[0])
# Expected: 6.0
print(max_flow_lp(list(range(10)), pipeline, 0, 9)[0])
# Expected: 65.0
```

(sec-09-04-min-cost-flow)=
## 最小费用流

### 连续最短路算法

{prf:ref}`alg-04-04-ssp`增加一个超级源点 $s$，从它向每个供应节点 $i$ 引一条容量为 $b_i$ 的弧 $(s,i)$，再增加一个超级汇点 $t$，从每个需求节点 $j$ 向它引一条容量为 $-b_j$ 的弧 $(j,t)$，这些弧的费用都为 0；然后反复沿残量网络（residual network）中一条最便宜的 $s$–$t$ 路输送尽可能多的流。残量网络中的弧 $(i,j)$ 在给弧 $(i,j)$ 增加流量时费用为 $c_{ij}$，在抵消弧 $(j,i)$ 上的流量时费用为 $-c_{ij}$。由于存在这些负费用，最短路用{ref}`sec-09-04-bellman-ford`中的 `bellman_ford` 求出。与算法本身一样，该函数假定流量下界为零且费用非负；流量下界须先行消去（见{ref}`sec-04-04-lower-bounds`）。

```python
def successive_shortest_paths(nodes: list, arcs: list, b: dict):
    """Minimum-cost flow by successive shortest paths.

    Parameters
    ----------
    nodes : list
        Node labels (not "s" or "t").
    arcs : list of tuple
        Arcs ``(i, j, c_ij, u_ij)`` with c_ij >= 0; u_ij may be INF.
    b : dict
        Net supplies b_i (positive supply, negative demand), sum zero.

    Returns
    -------
    z : number
        Minimum cost.
    x : dict
        Optimal arc flows.
    log : list of tuple
        (path, path cost, Delta, cumulative cost) for every augmentation.

    Raises
    ------
    ValueError
        For unbalanced data, negative costs, parallel or opposite arcs
        (split one of them with an extra node) or an infeasible problem.
    """
    if sum(b.values()) != 0:
        raise ValueError("net supplies must add up to zero")
    if any(c < 0 for _, _, c, _ in arcs):
        raise ValueError("successive shortest paths needs costs c_ij >= 0")
    # One residual pair per node pair: (i, j) and (j, i) must not both be
    # arcs, and no arc may be repeated.
    pairs = [(i, j) for i, j, _, _ in arcs]
    if len(set(pairs) | {(j, i) for i, j in pairs}) != 2 * len(pairs):
        raise ValueError("parallel or opposite arcs: split one with a node")
    r, cost = {}, {}

    def add_arc(i, j, c, u):
        """Residual pair: (i, j) with r = u, cost c; (j, i) with r = 0, -c."""
        r[(i, j)], cost[(i, j)] = u, c
        r[(j, i)], cost[(j, i)] = 0, -c

    for i, j, c, u in arcs:
        add_arc(i, j, c, u)
    for i in nodes:
        if b.get(i, 0) > 0:
            add_arc("s", i, 0, b[i])
        elif b.get(i, 0) < 0:
            add_arc(i, "t", 0, -b[i])
    z, log = 0, []
    while True:
        residual = [(i, j, cost[(i, j)]) for (i, j) in r if r[(i, j)] > 0]
        V, P, _ = bellman_ford(["s", "t", *nodes], residual, "s")
        if V["t"] == INF:
            break
        path = trace_path(P, "s", "t")
        delta = min(r[(i, j)] for i, j in zip(path, path[1:]))
        for i, j in zip(path, path[1:]):
            r[(i, j)] -= delta
            r[(j, i)] += delta
        z += delta * V["t"]
        log.append((path, V["t"], delta, z))
    if any(r[("s", i)] > 0 for i in nodes if b.get(i, 0) > 0):
        raise ValueError("the problem is infeasible")
    x = {(i, j): r[(j, i)] for i, j, _, _ in arcs}
    return z, x, log
```

弧 $(i,j)$ 上的流量由反向残余容量 $r_{ji}$ 读出；这就是{eq}`eq-04-03-recover`在不含方向相反弧对的网络中的形式。在{prf:ref}`ex-04-04-transportation`的运输问题上，出现了{numref}`tab-04-04-ssp`中的四次增广，包括使用后向弧（backward arc）$1\to B$ 的第四条路 $s\to A\to1\to B\to3\to t$。

```python
transport_nodes = ["A", "B", "1", "2", "3"]
transport_b = {"A": 40, "B": 50, "1": -20, "2": -25, "3": -45}
transport_cost = {
    ("A", "1"): 4,
    ("A", "2"): 6,
    ("A", "3"): 9,
    ("B", "1"): 2,
    ("B", "2"): 7,
    ("B", "3"): 5,
}
transport_arcs = [(i, j, c, INF) for (i, j), c in transport_cost.items()]
z_tr, x_tr, log_tr = successive_shortest_paths(
    transport_nodes, transport_arcs, transport_b
)
for path, path_cost, delta, cumulative in log_tr:
    print("-".join(path), path_cost, delta, cumulative)
# Expected:
#   s-B-1-t 2 20 40
#   s-B-3-t 5 30 190
#   s-A-2-t 6 25 340
#   s-A-1-B-3-t 7 15 445
print(z_tr, {f"{i}{j}": f for (i, j), f in x_tr.items() if f > 0})
# Expected: 445 {'A1': 15, 'A2': 25, 'B1': 5, 'B3': 45}
```

方案 $x_{A1}=15$、$x_{A2}=25$、$x_{B1}=5$、$x_{B3}=45$ 的费用为每周 445 美元。

### 带流量下界的供应链

{prf:ref}`ex-04-04-supply-chain`中的供应链有流量下界和上界。如{ref}`sec-04-04-feasibility`所示，由于弧 $(4,6)$ 的下界为 100，原始数据是不可行的；修正后的数据取 $0\le x_{46}\le120$。LP {eq}`eq-04-04-mcf`直接处理界：`linprog` 对每条弧接受一个 `(l, u)` 对。求解状态 2 表示“不可行”。

```python
chain_nodes = [1, 2, 3, 4, 5, 6]
chain_b = {1: 100, 2: 200, 3: 50, 4: -150, 5: -80, 6: -120}
# arcs (i, j, c_ij, l_ij, u_ij); u = INF means no upper bound
chain_given = [
    (1, 2, 3, 0, INF),
    (1, 3, 4, 0, INF),
    (1, 4, 1, 50, 80),
    (2, 3, 5, 0, INF),
    (2, 5, 6, 0, INF),
    (3, 4, 1, 70, 120),
    (3, 5, 2, 0, INF),
    (4, 6, 2, 100, 120),
    (5, 6, 4, 0, INF),
]
chain = [a if a[:2] != (4, 6) else (4, 6, 2, 0, 120) for a in chain_given]


def min_cost_flow_lp(nodes: list, arcs: list, b: dict):
    """Solve the minimum-cost flow LP with bounds; return linprog's result."""
    A = incidence_matrix(nodes, arcs)
    cost = [c for _, _, c, _, _ in arcs]
    bounds = [(l, None if u == INF else u) for *_, l, u in arcs]
    return linprog(cost, A_eq=A, b_eq=[b[i] for i in nodes], bounds=bounds)


print(min_cost_flow_lp(chain_nodes, chain_given, chain_b).status)
# Expected: 2
res_chain = min_cost_flow_lp(chain_nodes, chain, chain_b)
print(res_chain.fun, res_chain.x)
# Expected: 1810.0 [  0.  20.  80.  50. 150. 120.   0.  50.  70.]
pi = res_chain.eqlin.marginals - res_chain.eqlin.marginals[-1]
print(pi)
# Expected: [ 9. 10.  5.  2.  4.  0.]
print(np.array([c - pi[i - 1] + pi[j - 1] for i, j, c, _, _ in chain]))
# Expected: [ 4.  0. -6.  0.  0. -2.  1.  0.  0.]
```

最优费用为 1 810 美元，各流量见{numref}`tab-04-04-supply-chain-solution`。对偶值规范化为 $\pi_6=0$ 后就是节点势 $\pi=(9,10,5,2,4,0)$；检验数 $\bar c_{ij}=c_{ij}-\pi_i+\pi_j$ 在严格位于界内的五条弧上为 0，在处于下界 0 的两条弧上为正，在处于上界的两条弧上为 $-6$ 和 $-2$，这正是{prf:ref}`thm-04-04-reduced-cost`中的最优性条件。

连续最短路算法要求流量下界为零。代换{eq}`eq-04-04-shift`，即 $x_{ij}=l_{ij}+x'_{ij}$，把必须输送的流量移入净供应量{eq}`eq-04-04-shifted-supply`，并产生常数费用 $\sum c_{ij}l_{ij}$。

```python
def remove_lower_bounds(arcs: list, b: dict):
    """Shift x = l + x'; return arcs (i, j, c, u - l), supplies, constant."""
    b_new = dict(b)
    shifted, constant = [], 0
    for i, j, c, l, u in arcs:
        if l > u:
            raise ValueError(f"lower bound above upper bound on ({i}, {j})")
        b_new[i] -= l
        b_new[j] += l
        constant += c * l
        shifted.append((i, j, c, u - l))
    return shifted, b_new, constant


chain_shifted, chain_b_shifted, constant = remove_lower_bounds(chain, chain_b)
print(chain_b_shifted, constant)
# Expected: {1: 50, 2: 200, 3: -20, 4: -30, 5: -80, 6: -120} 120
z_chain, x_shifted, log_chain = successive_shortest_paths(
    chain_nodes, chain_shifted, chain_b_shifted
)
x_chain = {(i, j): l + x_shifted[(i, j)] for i, j, _, l, _ in chain}
print(len(log_chain), z_chain + constant)
# Expected: 5 1810
print(list(x_chain.values()))
# Expected: [0, 20, 80, 50, 150, 120, 0, 50, 70]
```

平移后的供应量 $b'_1=50$、$b'_3=-20$、$b'_4=-30$ 和常数 120 美元与{ref}`sec-04-04-lower-bounds`相同；把该常数加到平移后问题的费用上，又得到 1 810 美元，流量也相同。NetworkX 的 `min_cost_flow` 和 `network_simplex` 也要求流量下界为零。它们用 `demand` 描述节点，即 $-b_i$（需求节点为正），并在数据不可行时（例如弧 $(4,6)$ 上原来的下界）抛出 `NetworkXUnfeasible`。

```python
def mcf_graph(arcs: list, b: dict) -> nx.DiGraph:
    """DiGraph for networkx min-cost flow: demand = -b_i, weight, capacity."""
    G = nx.DiGraph()
    for i, bi in b.items():
        G.add_node(i, demand=-bi)
    for i, j, c, u in arcs:
        G.add_edge(i, j, weight=c, **({} if u == INF else {"capacity": u}))
    return G


cost_nx, flow_nx = nx.network_simplex(
    mcf_graph(chain_shifted, chain_b_shifted)
)
print(
    cost_nx + constant,
    nx.min_cost_flow_cost(mcf_graph(chain_shifted, chain_b_shifted))
    + constant,
)
# Expected: 1810 1810
given_shifted, given_b, _ = remove_lower_bounds(chain_given, chain_b)
try:
    nx.network_simplex(mcf_graph(given_shifted, given_b))
except nx.NetworkXUnfeasible as err:
    print("NetworkXUnfeasible:", err)
# Expected: NetworkXUnfeasible: no flow satisfies all node demands
```

## 小结

- 网络可以存储为弧表（输入格式）、节点–弧关联矩阵（所有流 LP 的约束矩阵，{numref}`tab-04-01-incidence`）、长度矩阵（稠密算法）或邻接表（稀疏算法）；NetworkX 的有向关联矩阵符号相反。
- 带标号 $[\gamma(j),P(j)]$ 的 Prim 算法和使用并查集的 Kruskal 算法都返回{prf:ref}`ex-04-01-sewage`中费用为 1 900 万美元的树；穷举计数证实共有 1 020 棵生成树，其中四棵最优。
- Dijkstra 算法（数组版本和堆版本）给出长度为 12 的路 $1\to4\to3\to7$（{prf:ref}`ex-04-02-dijkstra`）；约化长度{eq}`eq-04-02-reduced-length`和 LP {eq}`eq-04-02-lp`为其最优性提供了证书。
- Floyd–Warshall 算法必须把中间节点 $k$ 放在外层循环；向量化版本重现了{prf:ref}`ex-04-02-floyd`中的 $\mD_4$ 和 $\mP_4$ 并检测出负圈，而把 $k$ 放在最内层的变体会漏掉路线。
- Bellman–Ford 算法重现了{numref}`tab-05-04-bellman-ford`，并在第 $n$ 次迭代中值仍有变化时指出负圈。
- 使用标号法或广度优先搜索的 Ford–Fulkerson 方法给出{prf:ref}`ex-04-03-diamond`的最大流 6 和割 $S=\{1\}$，以及输油管网的最大流 65；流量用{eq}`eq-04-03-recover`恢复。
- 连续最短路算法重现了{numref}`tab-04-04-ssp`，并在消去流量下界后给出费用为 1 810 美元的供应链方案；LP 的对偶值给出节点势 $\pi=(9,10,5,2,4,0)$。

## 参考文献

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows:
  Theory, algorithms, and applications*. Prentice Hall. (Chapter 2: network
  representations; Chapters 4–5: shortest paths, Dijkstra, label-correcting
  methods; Chapters 6–7: maximum flows, Edmonds–Karp; Chapter 9: successive
  shortest paths; Chapter 13: minimum spanning trees and union–find.)
- Bellman, R. (1958). On a routing problem. *Quarterly of Applied
  Mathematics, 16*(1), 87–90. <https://doi.org/10.1090/qam/102435>
- Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
  *Numerische Mathematik, 1*(1), 269–271.
  <https://doi.org/10.1007/BF01386390>
- Edmonds, J., & Karp, R. M. (1972). Theoretical improvements in
  algorithmic efficiency for network flow problems. *Journal of the ACM,
  19*(2), 248–264. <https://doi.org/10.1145/321694.321699>
- Floyd, R. W. (1962). Algorithm 97: Shortest path. *Communications of the
  ACM, 5*(6), 345. <https://doi.org/10.1145/367766.368168>
- Ford, L. R., Jr. (1956). *Network flow theory* (Report P-923). RAND
  Corporation.
- Ford, L. R., & Fulkerson, D. R. (1956). Maximal flow through a network.
  *Canadian Journal of Mathematics, 8*, 399–404.
  <https://doi.org/10.4153/CJM-1956-045-5>
- Kruskal, J. B. (1956). On the shortest spanning subtree of a graph and the
  traveling salesman problem. *Proceedings of the American Mathematical
  Society, 7*(1), 48–50. <https://doi.org/10.1090/S0002-9939-1956-0078686-7>
- Prim, R. C. (1957). Shortest connection networks and some
  generalizations. *Bell System Technical Journal, 36*(6), 1389–1401.
  <https://doi.org/10.1002/j.1538-7305.1957.tb01515.x>
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 6: network models, Sections 6.2–6.4.)
- Warshall, S. (1962). A theorem on Boolean matrices. *Journal of the ACM,
  9*(1), 11–12.
