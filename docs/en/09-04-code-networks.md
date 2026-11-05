# Network Algorithms in Python

## Overview

This chapter implements the network algorithms of {doc}`04-01-networks-mst`,
{doc}`04-02-shortest-paths`, {doc}`04-03-maximum-flow`, {doc}`04-04-minimum-cost-flow`
and {doc}`05-04-dp-shortest-paths` in Python. It starts with the data
structures that store a network (arc list, node–arc incidence matrix,
length matrix, adjacency lists) and with their NetworkX counterparts. Each
algorithm is then written out from the book's steps and notation: Prim's
algorithm with labels and Kruskal's algorithm with a union–find structure;
Dijkstra's algorithm in an array version and a heap version, with the
predecessor labels $P(j)$ and path recovery; the Floyd–Warshall algorithm
with the intermediate node in the outer loop; the Bellman–Ford method with
negative-cycle detection; the Ford–Fulkerson method with the labeling
procedure and with breadth-first search (Edmonds–Karp), together with the
minimum cut; and the successive-shortest-path method for minimum-cost flows.
Every implementation is run on the worked examples of the theory chapters
and checked against NetworkX and, for the three network linear programs
(LPs), against
`scipy.optimize.linprog` with the incidence matrix, used as in
{doc}`09-02-code-lp-modeling`. The code is not run when
the book is built: copy the blocks, in order, into a script or a notebook
and compare the printed results with the `# Expected:` comments.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- store a network as an arc list, an incidence matrix, a length matrix or
  adjacency lists, and convert between these forms and a NetworkX graph;
- implement Prim's and Kruskal's algorithms, Dijkstra's algorithm (array
  and heap versions), Floyd–Warshall and Bellman–Ford, and recover paths
  from predecessor labels;
- implement the Ford–Fulkerson method with labeling or breadth-first
  search, recover arc flows from residual capacities and read off a minimum
  cut;
- solve a minimum-cost flow problem with lower bounds by successive
  shortest paths; and
- check every result with NetworkX and with the linear programs built from
  the node–arc incidence matrix.
:::

(sec-09-04-setup)=
## Setup

All blocks of the chapter share one namespace. The first block imports
everything that is used later. Nodes are the integers (or names) of the
book's examples, so that printed labels can be compared with the tables of
the theory chapters directly. `INF` stands for $+\infty$: a missing arc, a
label that is not set yet, or an unbounded capacity.

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
## Representing a network

The four data structures of {doc}`04-01-networks-mst` are built here for the
directed road network of {prf:ref}`ex-04-01-road-network`. The **arc list** of
triples $(i,j,c_{ij})$ is the input format of all functions in this
chapter; the other structures are derived from it.

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

The **node–arc incidence matrix** has one row per node and one column per
arc; the column of arc $(i,j)$ is $\ve_i-\ve_j$, with $+1$ at the tail and
$-1$ at the head. It is the constraint matrix of every network-flow LP
below.

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

The matrix agrees with {numref}`tab-04-01-incidence`. Every column adds up to
zero, so $\vone^\T\mA=\vzero^\T$, and the rank is
$\abs{\mathcal N}-1=4$ because the network is connected.

The **length matrix** $\mD_0$ holds $c_{ij}$ for the arcs, 0 on the diagonal
and $+\infty$ elsewhere; it is the input of the array version of Dijkstra's
algorithm and of Floyd–Warshall. **Adjacency lists** $\mathcal A(i)$ hold,
for every node, the arcs that leave it, and suit algorithms that scan the
arcs of one node at a time. For an undirected network every edge is stored
in both directions.

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

The adjacency lists are $\mathcal A(1)=\{2,4,5\}$, $\mathcal A(2)=\{4,5\}$,
$\mathcal A(3)=\{1,4\}$, $\mathcal A(4)=\{5\}$, $\mathcal A(5)=\varnothing$,
as in {doc}`04-01-networks-mst`. A NetworkX `DiGraph` is a dictionary of
adjacency lists with the arc attributes stored under keys such as
`"weight"` or `"capacity"`. Two conventions differ from the book and are
worth checking once: `nx.incidence_matrix(..., oriented=True)` puts $-1$ at
the tail and $+1$ at the head, the opposite sign, and returns a sparse
matrix.

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

The 0–1 adjacency matrix and the out-degrees 3, 2, 2, 1, 0 match
{prf:ref}`ex-04-01-road-network`.

(sec-09-04-mst)=
## Minimum spanning trees

The sewage network of {prf:ref}`ex-04-01-sewage` is undirected; each link is
stored once as $(i,j,c_{ij})$ with $i<j$ and costs in \$ million.

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

### Prim's algorithm with labels

{prf:ref}`alg-04-01-prim` grows one tree from a start node $r$. The
$O(\abs{\mathcal N}^2)$ version described after
{prf:ref}`thm-04-01-prim` keeps for every unconnected node $j$ the label
$[\gamma(j),P(j)]$: the cost of the cheapest edge from $j$ to the connected
set $C_k$ and the node of $C_k$ at its other end. The node with the smallest
$\gamma$ joins next (ties: the lower-numbered node), and only the edges at
the new node need to be examined. A label changes only on a strict
improvement, so among tied edges the one found first is kept; with the nodes
scanned in increasing order this reproduces the tie rule of
{prf:ref}`ex-04-01-sewage-prim`.

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

The links enter in the order 5–7, 2–5, 1–2, 2–4, 4–6, 1–3 of
{numref}`tab-04-01-prim`, and the tree costs \$19 million.

### Kruskal's algorithm with union–find

{prf:ref}`alg-04-01-kruskal` needs to know whether the two ends of an edge lie
in the same component of the forest built so far. A **union–find** structure
answers this: every component is a tree of parent pointers whose root names
the component. `find` follows the pointers to the root (and halves the path
on the way, so later searches are shorter); `union` hangs the root of the
smaller component below the root of the larger one. This is a refinement of
the renaming scheme analyzed after {prf:ref}`thm-04-01-kruskal`; both keep
the trees shallow, and with path halving the cost per operation is almost
constant.

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

The decisions are those of {numref}`tab-04-01-kruskal`: 6–7 is rejected
because it would close the cycle 6–7–5–2–4–6, and the algorithm stops after
1–3 without looking at the seven remaining links.

### Alternative optima and a brute-force check

{prf:ref}`ex-04-01-sewage-kruskal` claims that the network has 1 020 spanning
trees and exactly four of them cost \$19 million. Both claims are easy to
check by computer. The matrix-tree theorem says that the number of
spanning trees equals any cofactor of the Laplacian matrix
$\mL=\diag(\deg)-(\text{adjacency matrix})$; the brute force tests all
$\binom{14}{6}=3003$ sets of six links with the union–find structure (six
links on seven nodes form a spanning tree exactly when none of them closes a
cycle).

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

These are the four layouts of {prf:ref}`ex-04-01-sewage-kruskal`. NetworkX
implements Kruskal's, Prim's and Borůvka's algorithms behind one function.
All three return a tree of weight 19; which of the four optimal trees comes
back depends on the internal tie-breaking.

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
## Dijkstra's algorithm

The network of {prf:ref}`ex-04-02-route` has seven nodes and twelve arcs, with
lengths in minutes.

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

### Array version on the length matrix

The first implementation follows {prf:ref}`alg-04-02-dijkstra` on the length
matrix. Each iteration scans the most recent permanent node $i$, improving
the labels $[V(j),P(j)]$ of the nodes without a permanent label, and then
makes the smallest temporary label permanent (ties: smaller node number).
A temporary label is a node with a finite $V(j)$ that is not yet permanent.
The search for the smallest label scans all nodes, $O(\abs{\mathcal N})$ per
iteration and $O(\abs{\mathcal N}^2)$ in total. The function also returns
the order in which the nodes became permanent.

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

The nodes become permanent in the order 1, 2, 4, 6, 3, 5, 7 of
{numref}`tab-04-02-dijkstra`, and the path $1\to4\to3\to7$ of length 12
minutes is read backwards from $P(7)=3$, $P(3)=4$, $P(4)=1$. Because node 7
is the farthest node, all labels are final when the run stops.

### Heap version on adjacency lists

For sparse networks the smallest temporary label is kept in a binary heap
(Python's `heapq`). When a label improves, the new pair $(V(j),j)$ is
pushed and the old one is left in the heap; it is recognized and skipped
when it is popped, because its node is permanent by then. This "lazy
deletion" keeps the code short and the running time
$O(\abs{\mathcal A}\log\abs{\mathcal N})$. The tie rule is automatic:
pairs with equal $V$ are popped in order of the node number.

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

The arcs $(P(j),j)$ are the shortest-path tree of
{numref}`fig-04-02-dijkstra-tree`. NetworkX returns the same distances; its
predecessor dictionary lists *all* predecessors on shortest paths, which
here is a single one per node.

```python
G_route = nx.DiGraph()
G_route.add_weighted_edges_from(route_arcs)
pred, dist = nx.dijkstra_predecessor_and_distance(G_route, 1)
print(dist == V_heap, pred[7], nx.dijkstra_path(G_route, 1, 7))
# Expected: True [3] [1, 4, 3, 7]
```

### Optimality check with reduced lengths

{prf:ref}`thm-04-02-optimality` gives a certificate that does not depend on
how the labels were found: the labels are distances if every reduced length
$\bar c_{ij}=c_{ij}+V(i)-V(j)$ is nonnegative. The arcs of the
shortest-path tree have $\bar c_{ij}=0$.

```python
reduced = {(i, j): c + V_heap[i] - V_heap[j] for i, j, c in route_arcs}
print(min(reduced.values()))
# Expected: 0
print([a for a, r in reduced.items() if r == 0])
# Expected: [(1, 2), (1, 4), (2, 5), (3, 7), (4, 3), (4, 6)]
```

### The shortest-path LP

The same answer comes from the LP {eq}`eq-04-02-lp`: one unit of flow from
$s=1$ to $t=7$, with the incidence matrix as constraint matrix. HiGHS
returns the dual values $\pi_i$ of the balance equations in
`res.eqlin.marginals`; with $V(j)=\pi_s-\pi_j$ they are the node potentials
of {eq}`eq-04-02-dual`.

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

The flow uses the arcs $(1,4)$, $(4,3)$, $(3,7)$ and the dual objective
$\pi_1-\pi_7=V(7)$ equals 12, as strong duality requires. (The other
potentials are not unique: nodes that are not on the path may get any
values consistent with $\pi_i-\pi_j\le c_{ij}$.)

### A negative arc

Both Dijkstra functions refuse negative lengths. The network of
{prf:ref}`ex-04-02-negative-arc` shows why. Adding a constant to all lengths
is no cure: with every length raised by 2, Dijkstra's algorithm finds the
path $1\to2$, but in the original network $1\to3\to2$ is shorter. The
Bellman–Ford method of {ref}`sec-09-04-bellman-ford` gives the correct
distance $V(2)=1$.

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
## All pairs: Floyd–Warshall

{prf:ref}`alg-04-02-floyd` updates the whole distance matrix once for each
intermediate node $k$. The recurrence {eq}`eq-04-02-floyd` builds $\mD_k$ from
the complete matrix $\mD_{k-1}$, so the loop over $k$ must be the
**outer** loop. With NumPy the two inner loops over $i$ and $j$ become one
array operation: `D[:, k, None] + D[None, k, :]` is the matrix of all sums
$d_{ik}+d_{kj}$. Row and column $k$ do not change in iteration $k$, so
updating in place is safe. The predecessor matrix stores the node before
$j$ on the best known path from $i$ to $j$, with 0 for "no predecessor"
(the dash of the book). A negative diagonal entry after iteration $k$
reveals a negative cycle.

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

Applied to the courier network of {prf:ref}`ex-04-02-floyd`, the function
returns the matrices $\mD_4$ and $\mP_4$ of the example.

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

A path is read backwards in row $i$ of $\mP$, exactly like the labels of
Dijkstra's algorithm: row 4 gives
$2\leftarrow p_{42}=1\leftarrow p_{41}=3\leftarrow p_{43}=4$.

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

These are the routes $4\to3\to1\to2$ (cost 10, that is, \$1 000),
$2\to3\to1\to4$ (cost 3) and $1\to4\to3$ (a net gain of \$100). With
$c_{14}=-6$ the cycle $1\to4\to3\to1$ has cost $-1$, and the check finds it
in iteration 3, as in the example.

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

**Why the loop order matters.** A tempting version puts the loop over $k$
innermost: for each pair $(i,j)$ it tries all intermediate nodes once,
$d_{ij}\leftarrow\min_k\{d_{ik}+d_{kj}\}$. This computes $d_{ij}$ from
entries $d_{ik}$, $d_{kj}$ that may still be far from final, and the pair
is never revisited. On the network of {prf:ref}`ex-04-02-route` it misses two
routes out of node 2.

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

The route $2\to5\to3\to1$ of length $5+8+10=23$ exists, but when the
incorrect loop handles the pair $(2,1)$ the entries $d_{23}$ and $d_{53}$ it
would need are still $\infty$ or too large. NetworkX agrees with the correct
version. `floyd_warshall_numpy` returns the matrix in the order of
`nodelist`, and `floyd_warshall_predecessor_and_distance` returns
dictionaries.

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
## Negative arcs: Bellman–Ford

Procedure 1 of {doc}`05-04-dp-shortest-paths` computes $V^{(k)}(j)$, the
length of a shortest walk from the origin to $j$ with at most $k$ arcs, by
the recurrence {eq}`eq-05-04-bellman-ford`. Each iteration uses only the
values of the previous iteration (the copy `new` below), so iteration $k$
really is "at most $k$ arcs". By {prf:ref}`thm-05-04-successive` the values
converge by iteration $n$ if and only if no negative cycle can be reached;
if they still change in iteration $n$, following predecessors from a node
that changed leads into a negative cycle.

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

On the network of {prf:ref}`ex-05-04-bellman-ford` (lengths in \$100,
negative = profit) the rows of `history` are the rows of
{numref}`tab-05-04-bellman-ford`; iteration 4 repeats iteration 3.

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

The shortest path to node 4 is $1\to3\to2\to4$ with length $5-4+3=4$
(\$400), found only in iteration 3 because it has three arcs. With
$c_{23}=3$ the cycle $2\to3\to2$ has length $-1$, the values change in
iteration $4=n$, and the predecessors lead into the cycle:

```python
bf_cycle = [(1, 2, 2), (1, 3, 5), (2, 3, 3), (3, 2, -4), (2, 4, 3), (3, 4, 1)]
try:
    bellman_ford(bf_nodes, bf_cycle, 1)
except ValueError as err:
    print(err)
# Expected: negative cycle [2, 3, 2] reachable from 1
```

NetworkX gives the same distances, raises `NetworkXUnbounded` when a
negative cycle is reachable, and can return such a cycle. The same function
also settles the negative-arc example of {prf:ref}`ex-04-02-negative-arc`.

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
## Maximum flow

### Residual capacities

The flow algorithms work on residual capacities only. Following
{doc}`04-03-maximum-flow`, every arc $(i,j)$ gets a reverse pair $(j,i)$, with
capacity 0 if that arc does not exist, and the residual capacity
{eq}`eq-04-03-residual` of the zero flow is $r_{ij}=u_{ij}$. Sending $\Delta$
units along a path lowers $r_{ij}$ and raises $r_{ji}$ by $\Delta$ on every
arc of the path. At the end the arc flows are recovered with
{eq}`eq-04-03-recover`, $x_{ij}=\max\{0,u_{ij}-r_{ij}\}$.

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

### Two ways to find an augmenting path

The **labeling procedure** {prf:ref}`alg-04-03-labelling` is a depth-first
search: from the current node $i$ it labels the unlabeled, undeleted node
$k$ with the largest residual $r_{ik}$ (ties: smallest number) and moves
there; at a dead end it deletes the node and backs up to its predecessor.
**Breadth-first search** instead labels nodes in the order in which
they are reached, so the first path found to $t$ has the fewest arcs; this
is the Edmonds–Karp rule discussed with the running time of the method in
{doc}`04-03-maximum-flow`. Both functions return the path (or `None`) and the set $S$ of
nodes reached from $s$, which is the source side of a minimum cut when no
path exists.

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

`max` returns the first of several maxima, and the candidates are sorted,
so ties go to the smallest node number as the algorithm requires.

### The Ford–Fulkerson method

{prf:ref}`alg-04-03-ford-fulkerson` repeats "find a path, augment by its
residual capacity $\Delta$" until no path is left. The search function is a
parameter, so the same code runs the labeling version and the Edmonds–Karp
version.

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

On the diamond network of {prf:ref}`ex-04-03-diamond` the labeling version
finds the three paths of the example, including the backtracking from node
2 in iteration 3, and the cut $S=\{1\}$.

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

The flows $x_{12}=4$, $x_{13}=2$, $x_{32}=1$, $x_{24}=5$, $x_{34}=1$ and
the minimum cut of capacity 6 agree with the example. Breadth-first search
takes the two-arc path $1\to3\to4$ before the three-arc path
$1\to3\to2\to4$; the value and the flows are the same.

```python
v_bfs, x_bfs, S_bfs, log_bfs = ford_fulkerson(diamond, 1, 4, search=bfs_path)
print(log_bfs)
# Expected: [([1, 2, 4], 4), ([1, 3, 4], 1), ([1, 3, 2, 4], 1)]
print(v_bfs, x_bfs == x, S_bfs)
# Expected: 6 True {1}
```

### The pipeline network

The pipeline network of {prf:ref}`ex-04-03-pipeline` needs the transformations
of {ref}`sec-04-03-modelling`: each two-way pipeline becomes two opposite
arcs, and a super-source 0 and a super-sink 9 are joined to the refineries
and terminals. Their capacities are infinite; as the book suggests, any
number larger than the total capacity of the other arcs can stand in for
$\infty$, which keeps all residuals finite.

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

The maximum flow is 65 thousand barrels per day, and the minimum cut is the
refinery cut $S=\{0,1,2,3\}$ of {numref}`fig-04-03-pipeline`, of capacity
$20+10+20+15=65$.

### Checks with NetworkX and with the LP

`nx.maximum_flow` returns the value and a dictionary of arc flows;
`nx.minimum_cut` returns the capacity and the two sides of a minimum cut.
An arc without a `"capacity"` attribute has infinite capacity, so the
super-source and super-sink arcs need no stand-in value there.

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

NetworkX may return a different minimum cut when there are several (the
diamond network has two); every one has the capacity 6. Finally, the LP
{eq}`eq-04-03-maxflow-lp` reads $\max\{v:\mA\vx=v(\ve_s-\ve_t),\
\vzero\le\vx\le\vu\}$. With the variables $(\vx,v)$ the equality
constraints are $[\mA\ \ -(\ve_s-\ve_t)]\,(\vx,v)=\vzero$, and `linprog`
minimizes $-v$.

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
## Minimum-cost flow

### Successive shortest paths

{prf:ref}`alg-04-04-ssp` adds a super-source $s$ with an arc $(s,i)$ of
capacity $b_i$ to every supply node and a super-sink $t$ with an arc
$(j,t)$ of capacity $-b_j$ from every demand node, both of cost 0, and then
repeatedly sends as much flow as possible along a cheapest $s$–$t$ path of
the residual network. A residual arc $(i,j)$ costs $c_{ij}$ when it adds
flow to the arc $(i,j)$ and $-c_{ij}$ when it cancels flow on the arc
$(j,i)$. Because of these negative costs the shortest paths are found with
`bellman_ford` from {ref}`sec-09-04-bellman-ford`. The function assumes zero
lower bounds and nonnegative costs, as the algorithm does; lower bounds
are removed first ({ref}`sec-04-04-lower-bounds`).

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

The flow on arc $(i,j)$ is read from the reverse residual $r_{ji}$; this is
{eq}`eq-04-03-recover` for a network without opposite arcs. On the
transportation problem of {prf:ref}`ex-04-04-transportation` the four
augmentations of {numref}`tab-04-04-ssp` appear, including the fourth path
$s\to A\to1\to B\to3\to t$, which uses the backward arc $1\to B$.

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

The plan $x_{A1}=15$, $x_{A2}=25$, $x_{B1}=5$, $x_{B3}=45$ costs \$445 per
week.

### The supply chain with lower bounds

The supply chain of {prf:ref}`ex-04-04-supply-chain` has lower and upper
bounds. As {ref}`sec-04-04-feasibility` shows, the data as given are
infeasible because of the lower bound 100 on arc $(4,6)$; the repaired data
use $0\le x_{46}\le120$. The LP {eq}`eq-04-04-mcf` handles bounds directly:
`linprog` takes one `(l, u)` pair per arc. Status 2 means "infeasible".

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

The optimal cost is \$1 810 with the flows of
{numref}`tab-04-04-supply-chain-solution`. The duals, normalized to
$\pi_6=0$, are the potentials $\pi=(9,10,5,2,4,0)$, and the reduced costs
$\bar c_{ij}=c_{ij}-\pi_i+\pi_j$ are 0 on the five arcs strictly inside
their bounds, positive on the two arcs at the lower bound 0 and $-6$, $-2$
on the two arcs at their upper bounds, the optimality conditions of
{prf:ref}`thm-04-04-reduced-cost`.

Successive shortest paths needs zero lower bounds. The substitution
{eq}`eq-04-04-shift`, $x_{ij}=l_{ij}+x'_{ij}$, moves the compulsory flows into
the net supplies {eq}`eq-04-04-shifted-supply` and a constant cost
$\sum c_{ij}l_{ij}$.

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

The shifted supplies $b'_1=50$, $b'_3=-20$, $b'_4=-30$ and the constant
\$120 are those of {ref}`sec-04-04-lower-bounds`, and adding the constant to
the cost of the shifted problem gives \$1 810 again, with the same flows.
NetworkX's `min_cost_flow` and `network_simplex` also need zero lower
bounds. They describe a node by its `demand`, which is $-b_i$ (positive for
a demand node), and raise `NetworkXUnfeasible` for infeasible data, such as
the original lower bound on $(4,6)$.

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

## Summary

- A network is stored as an arc list (input format), a node–arc incidence
  matrix (the constraint matrix of all flow LPs, {numref}`tab-04-01-incidence`),
  a length matrix (dense algorithms) or adjacency lists (sparse
  algorithms); NetworkX's oriented incidence matrix has the opposite sign.
- Prim's algorithm with labels $[\gamma(j),P(j)]$ and Kruskal's algorithm
  with union–find both return the \$19 million tree of
  {prf:ref}`ex-04-01-sewage`; a brute-force count confirms 1 020 spanning trees
  and four optimal ones.
- Dijkstra's algorithm (array and heap versions) gives the path
  $1\to4\to3\to7$ of length 12 ({prf:ref}`ex-04-02-dijkstra`); the reduced
  lengths {eq}`eq-04-02-reduced-length` and the LP {eq}`eq-04-02-lp` certify it.
- Floyd–Warshall must have the intermediate node $k$ in the outer loop; the
  vectorized version reproduces $\mD_4$ and $\mP_4$ of
  {prf:ref}`ex-04-02-floyd` and detects the negative cycle, while the
  $k$-innermost variant misses routes.
- Bellman–Ford reproduces {numref}`tab-05-04-bellman-ford` and names the
  negative cycle when the values still change in iteration $n$.
- Ford–Fulkerson with the labeling procedure or with breadth-first search
  gives the maximum
  flow 6 and the cut $S=\{1\}$ of {prf:ref}`ex-04-03-diamond`, and 65 for the
  pipeline; flows are recovered with {eq}`eq-04-03-recover`.
- Successive shortest paths reproduces {numref}`tab-04-04-ssp` and, after
  removing lower bounds, the \$1 810 supply-chain plan; LP duals give the
  potentials $\pi=(9,10,5,2,4,0)$.

## References

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
