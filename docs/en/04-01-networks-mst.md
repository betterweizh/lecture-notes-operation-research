# Networks and Minimum Spanning Trees

## Overview

Many infrastructure systems are made of places and the links between them:
road junctions and road segments, pumping stations and pipelines, switches
and cables. A **network** records exactly this structure, and a large family
of optimization models (shortest paths, maximum flows, minimum-cost flows)
is built on it. This chapter fixes the vocabulary used in the network
chapters (nodes, arcs, paths, cycles, connectivity, trees), proves the basic
facts about trees (a tree on $\abs{\mathcal N}$ nodes has
$\abs{\mathcal N}-1$ arcs, and adding or deleting one arc behaves in a
predictable way), and shows how networks are stored. It then solves the first
network optimization problem: connect all nodes at the least total cost, the
**minimum spanning tree** (MST) problem. Two optimality results, the cut
property and the cycle property, explain why greedy choices are optimal here;
from them follow the correctness of **Prim's algorithm**, which grows one
tree from a start node, and of **Kruskal's algorithm**, which adds the
cheapest links one at a time. Both are run in full on a sewage network with
seven stations. The algorithms are combinatorial and need no optimization
theory; only {ref}`sec-04-01-representations` uses matrices and rank from
{doc}`01-01-linear-algebra`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- define directed and undirected networks, paths, chains, cycles,
  connectivity, trees and spanning trees, and recognize them in a drawing;
- prove the equivalent characterizations of a tree and use the exchange
  properties of spanning trees;
- store a network as an arc list, a node–arc incidence matrix, an adjacency
  matrix or adjacency lists, and choose between them;
- state and prove the cut and cycle properties and the optimality conditions
  of a minimum spanning tree;
- run Prim's and Kruskal's algorithms by hand, justify their correctness and
  state their running times; and
- decide whether a minimum spanning tree is unique and list the alternatives
  when it is not.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $G=(\mathcal N,\mathcal A)$ | network with node set $\mathcal N$ and arc (edge) set $\mathcal A$ |
| $(i,j)$ | arc with tail $i$ and head $j$; in an undirected network the edge joining $i$ and $j$, also written $i$–$j$ |
| $\abs{\mathcal N}$, $\abs{\mathcal A}$ | number of nodes, number of arcs |
| $c_{ij}$ | length or cost of arc $(i,j)$ |
| $\deg(i)$ | degree of node $i$ (number of edges at $i$) |
| $\mathcal T$, $c(\mathcal T)$ | a (spanning) tree, given by its edge set, and its total cost $\sum_{(i,j)\in\mathcal T}c_{ij}$ |
| $[S,\bar S]$ | cut: the edges with one end in $S\subset\mathcal N$ and the other in $\bar S=\mathcal N\setminus S$ |
| $\mathcal W$ | a cycle |
| $C_k$, $\bar C_k$ | connected and unconnected node sets in iteration $k$ of Prim's algorithm |
| $[\gamma(j),P(j)]$ | label of an unconnected node in Prim's algorithm: cost of its cheapest edge to $C_k$ and the node of $C_k$ at the other end |
| $\mA$ | node–arc incidence matrix (the constraint matrix of the network-flow models) |

(sec-04-01-graphs)=
## Graphs and directed graphs

:::{prf:definition} Network
:label: def-04-01-network

A **network** (or graph) $G=(\mathcal N,\mathcal A)$ consists of a finite set
$\mathcal N$ of **nodes** (vertices) and a set $\mathcal A$ of **arcs**
(links), each joining two distinct nodes.

- In a **directed network** (digraph) every arc is an ordered pair $(i,j)$
  with **tail** $i$ and **head** $j$; it can be traversed only from $i$ to
  $j$.
- In an **undirected network** every arc is an unordered pair, called an
  **edge**; the edge $(i,j)$ is the same as $(j,i)$ and can be traversed both
  ways.

Numbers attached to the arcs, such as a length, cost or capacity, are the
**attributes** of the arcs. The length or cost of arc $(i,j)$ is written
$c_{ij}$.
:::

Throughout the network chapters there is at most one arc from $i$ to $j$ and
no arc from a node to itself, so an arc is identified by its two end nodes.
Two nodes joined by an arc are **adjacent**, and the arc is **incident** to
both. The **degree** $\deg(i)$ of a node in an undirected network is the
number of edges incident to it. Every edge has two ends, so

$$
\sum_{i\in\mathcal N}\deg(i)=2\abs{\mathcal A}.
$$ (eq-04-01-handshake)

In a directed network the **out-degree** of $i$ counts the arcs with tail $i$
and the **in-degree** the arcs with head $i$; each sum over all nodes equals
$\abs{\mathcal A}$. An undirected network can always be turned into a directed
one by replacing each edge $(i,j)$ with the two arcs $(i,j)$ and $(j,i)$, both
with the attribute of the edge; the shortest-path and flow algorithms of
{doc}`04-02-shortest-paths` and {doc}`04-03-maximum-flow` work on this directed
version.

:::{prf:example} A small road network
:label: ex-04-01-road-network

Five road junctions are joined by eight road segments
({numref}`fig-04-01-road-network`). With the junctions numbered 1 to 5,

$$
\begin{aligned}
\mathcal N&=\{1,2,3,4,5\},\\
\mathcal A&=\{(3,1),(3,4),(1,2),(1,4),(1,5),(2,5),(2,4),(4,5)\},
\end{aligned}
$$

and the segment lengths $c_{ij}$ in miles are

| Arc $(i,j)$ | $(3,1)$ | $(3,4)$ | $(1,2)$ | $(1,4)$ | $(1,5)$ | $(2,5)$ | $(2,4)$ | $(4,5)$ |
|---|---|---|---|---|---|---|---|---|
| $c_{ij}$ | 1 | 2 | 0.5 | 2.5 | 5 | 1.2 | 1.7 | 3 |

Find the degrees of the nodes when all roads are two-way (panel (a)), and
the out- and in-degrees when every road is one-way in the direction in which
its arc is listed (panel (b)). Check both against the number of arcs.

*(Some drawings of this network orient the arc between nodes 2 and 4 from 4
to 2; here it is $(2,4)$ throughout, as in the arc set above.)*

```{figure} figures/04-01-road-network.svg
:name: fig-04-01-road-network
:width: 700px
:class: dark-light
:alt: Five nodes joined by eight links with lengths, drawn once without and once with arrow heads

The road network of {prf:ref}`ex-04-01-road-network`: (a) undirected, with
segment lengths in miles; (b) directed, every arc oriented as listed in
$\mathcal A$.
```
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

As an undirected network (all roads two-way, panel (a)) the degrees are
$\deg(1)=4$, $\deg(2)=3$, $\deg(3)=2$, $\deg(4)=4$ and $\deg(5)=3$; they add
up to $16=2\times 8$, as {eq}`eq-04-01-handshake` requires.

If every road is one-way in the direction in which its arc is listed
(panel (b)), the network is directed. The out-degrees of nodes 1 to 5 are
3, 2, 2, 1, 0 and the in-degrees 1, 1, 0, 3, 3; both lists add up to
$\abs{\mathcal A}=8$. Node 3 has no entering arc and node 5 no leaving arc.
:::

(sec-04-01-paths)=
## Paths, cycles and connectivity

:::{prf:definition} Paths, chains and cycles
:label: def-04-01-path

Let $G=(\mathcal N,\mathcal A)$ be a network.

- In an undirected network, a **path** from $i$ to $j$ is a sequence of
  distinct nodes $i=i_1,i_2,\dots,i_p=j$ such that every consecutive pair
  $i_k,i_{k+1}$ is joined by an edge. Its edges are distinct, and it has
  $p-1$ of them.
- In a directed network, a sequence with the same property in which arcs may
  be used in either direction is a **chain** (an undirected path). If every
  arc is used in its own direction, $(i_k,i_{k+1})\in\mathcal A$ for all $k$,
  it is a **directed path**.
- A **cycle** is a path that returns to its first node: distinct nodes
  $i_1,\dots,i_p$ with $p\ge3$, edges between consecutive nodes and an edge
  $(i_p,i_1)$. A **directed cycle** uses the arcs $(i_1,i_2),\dots,
  (i_{p-1},i_p),(i_p,i_1)$ in their own direction ($p\ge2$). A directed
  network without directed cycles is **acyclic**.

The **length** of a path, chain or cycle is the sum of the attributes of its
arcs.
:::

:::{prf:definition} Connectivity and sub-networks
:label: def-04-01-connected

An undirected network is **connected** if every two distinct nodes are
joined by at least one path. A directed network is (weakly) connected if its
undirected version is connected, and **strongly connected** if there is a
directed path from every node to every other node. A **sub-network** of $G$
is a network $G'=(\mathcal N',\mathcal A')$ with $\mathcal N'\subseteq\mathcal
N$, $\mathcal A'\subseteq\mathcal A$ and both ends of every arc of
$\mathcal A'$ in $\mathcal N'$; it is **spanning** if $\mathcal N'=\mathcal
N$. The maximal connected sub-networks of $G$ are its **components**.
:::

In the road network of {prf:ref}`ex-04-01-road-network`, read as an undirected
network ({numref}`fig-04-01-terminology`):

- 3–1–5 is a path from node 3 to node 5 of length $1+5=6$ miles (panel a);
- 1–2–4–1 is a cycle of length $0.5+1.7+2.5=4.7$ miles (panel b). The route
  3–1–2–4–1–5 joins nodes 3 and 5 but visits node 1 twice, so it is not a
  path; it contains the cycle 1–2–4–1, and cutting the cycle out leaves the
  path 3–1–5;
- the network is connected, but the spanning sub-network obtained by
  deleting the three edges at node 5 is not: node 5 forms a component of its
  own (panel c).

Now read the network as directed (panel (b) of
{numref}`fig-04-01-road-network`). The sequence 3, 1, 2, 5 is a directed path of
length $1+0.5+1.2=2.7$ miles (panel f), and there are four more directed
paths from 3 to 5: 3→1→2→4→5, 3→1→4→5, 3→1→5 and 3→4→5. The cycle 1–2–4–1 is
not a directed cycle, because closing it would need the arc $(4,1)$, and
$\mathcal A$ contains $(1,4)$ instead; the nodes, edges and arcs are the same,
only the directions differ. In fact the directed network is acyclic: in the
order 3, 1, 2, 4, 5 every arc leads from an earlier node to a later one, so
no directed path can return to its start. It is connected but not strongly
connected, because no arc enters node 3 and no arc leaves node 5.

```{figure} figures/04-01-terminology.svg
:name: fig-04-01-terminology
:width: 700px
:class: dark-light
:alt: Six copies of the five-node network, each highlighting a different sub-network

Terminology on the road network; the highlighted edges and nodes form the
object, dashed gray edges are the rest of the network. (a) The path 3–1–5.
(b) The cycle 1–2–4–1. (c) A spanning sub-network that is not connected
(node 5 is isolated). (d) A tree that does not span the network (node 5 is
missing). (e) A spanning tree. (f) The directed path 3→1→2→5 in the directed
version.
```

(sec-04-01-trees)=
## Trees

:::{prf:definition} Trees, forests and spanning trees
:label: def-04-01-tree

An undirected network is a **forest** if it has no cycle, and a **tree** if
it is a connected forest. A **spanning tree** of a network $G$ is a spanning
sub-network of $G$ that is a tree. A node of degree 1 in a tree is a
**leaf**. A tree is written $\mathcal T$ and identified with its set of
edges.
:::

In {numref}`fig-04-01-terminology`, panel (d) is a tree on the nodes 1 to 4 but
not a spanning tree of the road network, while panel (e) is a spanning tree.
Trees are the cheapest way to keep a set of nodes connected, and the next
results describe them completely.

:::{prf:lemma} Leaves
:label: lem-04-01-leaves

Every tree with at least two nodes has at least two leaves.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

A tree with two or more nodes is connected and therefore has an edge. Take a
path $i_1,\dots,i_p$ with the largest number of edges ($p\ge2$; it exists
because the network is finite). Node $i_1$ is adjacent to $i_2$. It has no
neighbor outside the path, since that neighbor would extend the path, and
no neighbor $i_k$ with $k\ge3$, since the edge $(i_1,i_k)$ and the path
$i_1,\dots,i_k$ would form a cycle. Hence $\deg(i_1)=1$, and in the same way
$\deg(i_p)=1$.
:::

:::{prf:theorem} Characterizations of a tree
:label: thm-04-01-tree

For an undirected network $G=(\mathcal N,\mathcal A)$ the following are
equivalent:

1. $G$ is a tree (connected and without cycles);
2. $G$ is connected and $\abs{\mathcal A}=\abs{\mathcal N}-1$;
3. $G$ has no cycle and $\abs{\mathcal A}=\abs{\mathcal N}-1$;
4. every two nodes of $G$ are joined by exactly one path.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

*A tree has $\abs{\mathcal N}-1$ edges.* By induction on
$\abs{\mathcal N}$. A tree with one node has no edge. If
$\abs{\mathcal N}\ge2$, delete a leaf and its edge ({prf:ref}`lem-04-01-leaves`).
What remains has no cycle, and it is connected, because a leaf cannot be an
interior node of a path (interior nodes have degree at least 2). By induction
it has $\abs{\mathcal N}-2$ edges, so the tree has $\abs{\mathcal N}-1$. This
proves 1 ⇒ 2 and 1 ⇒ 3. Applied to each component, it also shows that a
forest with $q$ components has $\abs{\mathcal N}-q$ edges.

*3 ⇒ 1.* A network without cycles is a forest; with $q$ components it has
$\abs{\mathcal N}-q$ edges. Then $\abs{\mathcal N}-q=\abs{\mathcal N}-1$ gives
$q=1$: the network is connected.

*2 ⇒ 1.* If $G$ had a cycle, deleting one edge of the cycle would leave it
connected (the rest of the cycle still joins the two ends of the deleted
edge). Repeating this until no cycle is left produces a tree with fewer than
$\abs{\mathcal N}-1$ edges, which contradicts 1 ⇒ 2.

*1 ⇒ 4.* A tree is connected, so there is at least one path. Suppose two
different paths $P$ and $Q$ join $i$ and $j$. Walking from $i$, let $u$ be the
last node before the two paths separate, and let $w$ be the first node after
$u$ on $P$ that also lies on $Q$ ($w$ exists because $j$ lies on both). The
part of $P$ from $u$ to $w$ and the part of $Q$ from $w$ back to $u$ share
only $u$ and $w$, and they are not the same single edge, so together they
form a cycle, a contradiction.

*4 ⇒ 1.* The network is connected. A cycle would give two paths between two
of its consecutive nodes (the edge between them and the rest of the cycle),
contradicting uniqueness.
:::

The theorem has two consequences that drive every spanning-tree algorithm.

:::{prf:corollary} Exchange properties of spanning trees
:label: cor-04-01-exchange

Let $\mathcal T$ be a spanning tree of a connected network $G$.

1. Adding an edge $(k,l)\notin\mathcal T$ to $\mathcal T$ creates exactly one
   cycle: the edge $(k,l)$ together with the unique tree path from $k$ to
   $l$. Deleting any edge of this cycle gives another spanning tree.
2. Deleting an edge $(k,l)\in\mathcal T$ splits $\mathcal T$ into exactly two
   trees, with node sets $S\ni k$ and $\bar S\ni l$. Adding any edge of $G$
   with one end in $S$ and the other in $\bar S$ gives another spanning tree.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

1. A cycle through the new edge consists of $(k,l)$ and a path from $l$ back
   to $k$ in $\mathcal T$, and there is only one such path (Theorem
   {prf:ref}`thm-04-01-tree`, part 4). Deleting an edge of this cycle leaves a
   connected network with $\abs{\mathcal N}-1$ edges, which is a tree by
   part 2.
2. After the deletion there are $\abs{\mathcal N}-2$ edges and no cycle, so
   the forest has $\abs{\mathcal N}-(\abs{\mathcal N}-2)=2$ components. An
   edge between them reconnects the network with $\abs{\mathcal N}-1$ edges.
:::

Every connected network has a spanning tree: delete edges of cycles, one at a
time, until no cycle is left (the argument of 2 ⇒ 1). It usually has very many.
By **Kirchhoff's matrix-tree theorem** the number of spanning trees equals any
cofactor of the Laplacian matrix (degrees on the diagonal, $-1$ for each pair
of adjacent nodes); for the road network this gives 40 spanning trees. For the
complete network on $\abs{\mathcal N}$ nodes (every pair joined) the count is
$\abs{\mathcal N}^{\abs{\mathcal N}-2}$ (Cayley's formula), already about
$2.6\times10^{23}$ for 20 nodes, so spanning-tree problems cannot be solved by
listing the trees.

(sec-04-01-representations)=
## Representing a network

A drawing is convenient for people; an algorithm needs one of the following
data structures.

**Arc list.** The arcs are listed as triples $(i,j,c_{ij})$. This needs
storage proportional to $\abs{\mathcal A}$, is the usual input format, and is
what Kruskal's algorithm sorts.

**Node–arc incidence matrix.** The matrix $\mA\in\R^{\abs{\mathcal N}\times
\abs{\mathcal A}}$ has one row per node and one column per arc. The column of
arc $(i,j)$ is $\ve_i-\ve_j$: it has $+1$ in row $i$ (the tail), $-1$ in row
$j$ (the head) and zeros elsewhere. {numref}`tab-04-01-incidence` shows it for
the directed road network. Row $i$ has as many entries $+1$ as the out-degree
of $i$ and as many $-1$ as its in-degree.

```{list-table} Node–arc incidence matrix of the directed road network
:name: tab-04-01-incidence
:header-rows: 1
:stub-columns: 1

* - Node
  - $(1,2)$
  - $(1,4)$
  - $(1,5)$
  - $(2,4)$
  - $(2,5)$
  - $(3,1)$
  - $(3,4)$
  - $(4,5)$
* - 1
  - 1
  - 1
  - 1
  - 0
  - 0
  - −1
  - 0
  - 0
* - 2
  - −1
  - 0
  - 0
  - 1
  - 1
  - 0
  - 0
  - 0
* - 3
  - 0
  - 0
  - 0
  - 0
  - 0
  - 1
  - 1
  - 0
* - 4
  - 0
  - −1
  - 0
  - −1
  - 0
  - 0
  - −1
  - 1
* - 5
  - 0
  - 0
  - −1
  - 0
  - −1
  - 0
  - 0
  - −1
```

Every column contains one $+1$ and one $-1$, so the rows add up to zero,
$\vone^\T\mA=\vzero^\T$, and $\rank\mA\le\abs{\mathcal N}-1$. For a connected
network the rank is exactly $\abs{\mathcal N}-1$ (here 4). The matrix is the
constraint matrix of the flow-conservation equations $\mA\vx=\vb$ of every
network-flow model ({doc}`04-02-shortest-paths`, {doc}`04-03-maximum-flow`,
{doc}`04-04-minimum-cost-flow`), and its special structure is the reason why
these linear programs have integer optimal solutions
({prf:ref}`thm-04-04-tu`, {prf:ref}`thm-04-04-integrality`). As a storage scheme it is wasteful: only
$2\abs{\mathcal A}$ of its $\abs{\mathcal N}\abs{\mathcal A}$ entries are
nonzero. For an undirected network both entries of a column are $+1$.

**Adjacency matrix.** The $\abs{\mathcal N}\times\abs{\mathcal N}$ matrix with
entry 1 in position $(i,j)$ when $(i,j)\in\mathcal A$ and 0 otherwise; it is
symmetric for an undirected network. Its weighted form holds $c_{ij}$ for the
arcs, 0 on the diagonal and $+\infty$ for missing arcs; this **length
matrix** is the starting point of the Floyd–Warshall algorithm
({prf:ref}`alg-04-02-floyd`). Storage is $\abs{\mathcal N}^2$ and testing whether
an arc exists takes constant time, which suits dense networks.

**Adjacency lists.** For each node $i$ the list $\mathcal A(i)=\{j:(i,j)\in
\mathcal A\}$ of nodes reached by an arc out of $i$, with the arc attributes.
For the directed road network,

$$
\begin{aligned}
\mathcal A(1)&=\{2,4,5\},\\
\mathcal A(2)&=\{4,5\},\\
\mathcal A(3)&=\{1,4\},\\
\mathcal A(4)&=\{5\},\\
\mathcal A(5)&=\varnothing .
\end{aligned}
$$

Storage is proportional to $\abs{\mathcal N}+\abs{\mathcal A}$, and scanning
the arcs out of $i$ takes time proportional to $\abs{\mathcal A(i)}$. Stored
in one array sorted by tail node with a pointer to the first arc of each node,
this is the **forward-star** representation. Prim's and Dijkstra's
algorithms, which repeatedly scan the arcs of one node, are implemented with
adjacency lists; {doc}`09-04-code-networks` shows the corresponding Python data
structures.

(sec-04-01-mst-problem)=
## The minimum spanning tree problem

:::{prf:definition} Minimum spanning tree
:label: def-04-01-mst

Let $G=(\mathcal N,\mathcal A)$ be a connected undirected network with edge
costs $c_{ij}$. The **minimum spanning tree problem** is to find a spanning
tree $\mathcal T$ of $G$ with the least total cost

$$
c(\mathcal T)=\sum_{(i,j)\in\mathcal T}c_{ij}.
$$ (eq-04-01-tree-cost)
:::

The task "connect all nodes as cheaply as possible" asks for a cheapest
connected spanning sub-network; with positive costs such a sub-network has no
cycle (an edge of a cycle could be removed without disconnecting anything), so
it is a spanning tree. Two further remarks will be useful. First, every
spanning tree has the same number $\abs{\mathcal N}-1$ of edges, so adding a
constant to every cost changes all tree costs by the same amount and leaves
the optimal trees unchanged; negative costs are therefore harmless. Second,
the algorithms below only compare costs, so the costs may be any real numbers.

Typical applications are the design of networks whose only requirement is
connectivity:

- **infrastructure layout**: pipelines from a reservoir or to a treatment
  plant, cable and fiber networks, power lines, rural roads, where the cost of
  a link is its construction cost;
- **clustering**: with costs measuring dissimilarity, deleting the $q-1$ most
  expensive edges of a minimum spanning tree leaves $q$ clusters
  (single-linkage clustering); {prf:ref}`ex-08-05-machine-cells` uses this to
  group machines into manufacturing cells;
- **approximation**: a minimum spanning tree is the starting point of
  heuristics for harder problems, such as the traveling salesman problem.

The problem is stated for undirected networks. Its directed counterpart, a
cheapest set of arcs that reaches every node from a root (a minimum-cost
arborescence), needs different algorithms and is not treated here.

:::{prf:example} A sewage collection network
:label: ex-04-01-sewage

Seven stations of a deep-tunnel sewage system must be linked by tunnels so
that wastewater from every station can reach the main treatment station,
node 7. Fourteen tunnels are technically possible
({numref}`fig-04-01-sewage-network`); pairs of stations without a link would be
too expensive to connect directly. The construction costs in \$ million are

| Link | 1–2 | 1–3 | 1–4 | 1–5 | 2–3 | 2–4 | 2–5 | 3–4 | 3–6 | 3–7 | 4–5 | 4–6 | 5–7 | 6–7 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Cost | 1 | 5 | 7 | 9 | 6 | 4 | 3 | 5 | 10 | 8 | 8 | 3 | 3 | 4 |

Which tunnels should be built so that all stations are connected at the least
total cost?

```{figure} figures/04-01-sewage-network.svg
:name: fig-04-01-sewage-network
:width: 560px
:class: dark-light
:alt: Seven stations joined by fourteen possible tunnel links with their costs; node 7 is marked

The sewage network of {prf:ref}`ex-04-01-sewage`. Labels are construction costs
in \$ million; node 7 is the main station.
```
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Any connected layout lets every station reach node 7, so the cheapest layout
is a minimum spanning tree. It has $7-1=6$ links, chosen from 14. Of the
$\binom{14}{6}=3003$ sets of six links, the matrix-tree theorem shows that
1 020 are spanning trees. One of them is

$$
\begin{gathered}
\mathcal T=\{1\text{–}2,\ 2\text{–}5,\ 5\text{–}7,\ 1\text{–}3,\ 3\text{–}6,\ 4\text{–}5\},\\
c(\mathcal T)=1+3+3+5+10+8=30,
\end{gathered}
$$

a feasible layout costing \$30 million. Prim's algorithm
({prf:ref}`ex-04-01-sewage-prim`) and Kruskal's algorithm
({prf:ref}`ex-04-01-sewage-kruskal`) find the optimum without listing the trees.
:::

(sec-04-01-optimality)=
## Optimality conditions: the cut and cycle properties

For a set $S$ of nodes with $\varnothing\ne S\ne\mathcal N$, the **cut**
$[S,\bar S]$ is the set of edges with one end in $S$ and the other in
$\bar S=\mathcal N\setminus S$. Every spanning tree contains at least one edge
of every cut, since it must connect $S$ to $\bar S$. Two simple exchange
arguments show which edges belong to minimum spanning trees.

:::{prf:theorem} Cut property
:label: thm-04-01-cut-property

Let $F$ be a set of edges contained in some minimum spanning tree, and let
$[S,\bar S]$ be a cut that contains no edge of $F$. If $(k,l)$ is an edge of
least cost in $[S,\bar S]$, then $F\cup\{(k,l)\}$ is contained in some minimum
spanning tree. If $(k,l)$ is the only edge of least cost in the cut, every
minimum spanning tree contains it.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\mathcal T^*$ be a minimum spanning tree with $F\subseteq\mathcal T^*$.
If $(k,l)\in\mathcal T^*$ there is nothing to prove. Otherwise, adding
$(k,l)$ to $\mathcal T^*$ closes a cycle ({prf:ref}`cor-04-01-exchange`, part 1)
that leaves $S$ through $(k,l)$ and must come back through another edge
$(p,q)$ of the cut; $(p,q)$ is a tree edge, and $(p,q)\notin F$ because no
edge of $F$ lies in the cut. The tree $\mathcal T'=\mathcal T^*\cup\{(k,l)\}
\setminus\{(p,q)\}$ contains $F\cup\{(k,l)\}$ and costs

$$
c(\mathcal T')=c(\mathcal T^*)+c_{kl}-c_{pq}\le c(\mathcal T^*),
$$

because $(k,l)$ is a cheapest edge of the cut. Hence $\mathcal T'$ is also a
minimum spanning tree. If $(k,l)$ is strictly cheaper than every other edge
of the cut, the inequality is strict for every minimum spanning tree
$\mathcal T^*$ that avoids $(k,l)$ (take $F=\varnothing$), which is
impossible; so every minimum spanning tree contains $(k,l)$.
:::

:::{prf:theorem} Cycle property
:label: thm-04-01-cycle-property

Let $\mathcal W$ be a cycle of $G$ and $(k,l)$ an edge of $\mathcal W$ that is
strictly more expensive than every other edge of $\mathcal W$. Then no minimum
spanning tree contains $(k,l)$. If $(k,l)$ is only a most expensive edge of
$\mathcal W$ (ties allowed), some
minimum spanning tree avoids it.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\mathcal T^*$ be a minimum spanning tree containing $(k,l)$. Deleting
$(k,l)$ splits it into two trees with node sets $S\ni k$ and $\bar S\ni l$
({prf:ref}`cor-04-01-exchange`, part 2). The rest of the cycle, $\mathcal W$ without
$(k,l)$, is a path from $k$ to $l$, so it contains an edge $(p,q)\ne(k,l)$
with one end in $S$ and one in $\bar S$. Then
$\mathcal T'=\mathcal T^*\setminus\{(k,l)\}\cup\{(p,q)\}$ is a spanning tree
with $c(\mathcal T')=c(\mathcal T^*)-c_{kl}+c_{pq}$. If $c_{pq}<c_{kl}$ this
contradicts the optimality of $\mathcal T^*$; if only $c_{pq}\le c_{kl}$,
$\mathcal T'$ is a minimum spanning tree without $(k,l)$.
:::

Together the two properties give a test that a proposed tree is optimal,
which needs only one cycle per non-tree edge.

:::{prf:theorem} Optimality conditions
:label: thm-04-01-mst-optimality

A spanning tree $\mathcal T^*$ is a minimum spanning tree if and only if, for
every edge $(k,l)\notin\mathcal T^*$, every edge $(i,j)$ on the tree path
between $k$ and $l$ satisfies $c_{ij}\le c_{kl}$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

*Necessity.* If $c_{ij}>c_{kl}$ for an edge $(i,j)$ on the tree path of a
non-tree edge $(k,l)$, then $\mathcal T^*\cup\{(k,l)\}\setminus\{(i,j)\}$ is
a spanning tree ({prf:ref}`cor-04-01-exchange`, part 1) that is cheaper than
$\mathcal T^*$.

*Sufficiency.* Let $\mathcal T^*$ satisfy the conditions, and among all
minimum spanning trees let $\mathcal T^o$ share the most edges with
$\mathcal T^*$. Suppose $\mathcal T^o\ne\mathcal T^*$ and take
$(i,j)\in\mathcal T^*\setminus\mathcal T^o$. Deleting $(i,j)$ from
$\mathcal T^*$ defines a cut $[S,\bar S]$ in which $(i,j)$ is the only edge
of $\mathcal T^*$. Adding $(i,j)$ to $\mathcal T^o$ closes a cycle that
crosses this cut at a second edge $(k,l)\in\mathcal T^o$, and
$(k,l)\notin\mathcal T^*$. The tree path in $\mathcal T^*$ between $k$ and $l$
must cross the cut, so it uses $(i,j)$, and the conditions give
$c_{ij}\le c_{kl}$. The necessity part applied to $\mathcal T^o$ (with the
non-tree edge $(i,j)$ and the edge $(k,l)$ on its cycle) gives
$c_{kl}\le c_{ij}$. Hence $c_{ij}=c_{kl}$, and $\mathcal T^o\cup\{(i,j)\}
\setminus\{(k,l)\}$ is a minimum spanning tree sharing one more edge with
$\mathcal T^*$, a contradiction. So $\mathcal T^*=\mathcal T^o$ is optimal.
:::

:::{prf:corollary} Uniqueness
:label: cor-04-01-unique

If every non-tree edge $(k,l)$ is strictly more expensive than every edge on
its tree path, $\mathcal T^*$ is the only minimum spanning tree; in particular
the minimum spanning tree is unique when all edge costs are different. If a
non-tree edge $(k,l)$ costs the same as an edge $(i,j)$ on its tree path,
$\mathcal T^*\cup\{(k,l)\}\setminus\{(i,j)\}$ is another minimum spanning
tree.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

The sufficiency part of the proof of {prf:ref}`thm-04-01-mst-optimality` shows that a minimum spanning tree
$\mathcal T^o\ne\mathcal T^*$ produces a non-tree edge $(k,l)$ with
$c_{kl}=c_{ij}$ for an edge $(i,j)$ on its tree path, which strict
inequalities exclude. The last statement is the exchange of
{prf:ref}`cor-04-01-exchange`, which does not change the cost.
:::

(sec-04-01-prim)=
## Prim's algorithm

The cut property suggests growing a single tree. Start from any node; at every
iteration look at the cut between the nodes already connected and the rest,
and add a cheapest edge of that cut. The resulting "connected set" method
was published by Prim (1957); it had been found earlier by Jarník (1930).

:::{prf:algorithm} Prim's algorithm
:label: alg-04-01-prim

**Input:** a connected undirected network $G=(\mathcal N,\mathcal A)$ with
costs $c_{ij}$, and a start node $r$.
**Output:** a minimum spanning tree $\mathcal T$.

1. **Initialize.** Set $C_1=\{r\}$, $\bar C_1=\mathcal N\setminus\{r\}$,
   $\mathcal T=\varnothing$ and $k=1$.
2. **Stopping test.** If $\bar C_k=\varnothing$, stop: $\mathcal T$ is a
   minimum spanning tree.
3. **Select.** Choose an edge $(i^*,j^*)$ with $i^*\in C_k$ and
   $j^*\in\bar C_k$ of least cost,
   $c_{i^*j^*}=\min\{c_{ij}:(i,j)\in[C_k,\bar C_k]\}$ (ties broken
   arbitrarily).
4. **Connect.** Set $\mathcal T\leftarrow\mathcal T\cup\{(i^*,j^*)\}$,
   $C_{k+1}=C_k\cup\{j^*\}$, $\bar C_{k+1}=\bar C_k\setminus\{j^*\}$,
   $k\leftarrow k+1$, and go to Step 2.
:::

No test for cycles is needed: the new edge has one end in the connected set
and the other, $j^*$, outside it, and a cycle through $j^*$ would need two tree
edges at $j^*$.
*(Statements of the algorithm that also require "adding $j^*$ must not create
a cycle" impose a condition that always holds.)* Each iteration moves one
node from $\bar C_k$ to $C_k$, so the algorithm stops after exactly
$\abs{\mathcal N}-1$ iterations with $\abs{\mathcal N}-1$ edges, as
{prf:ref}`thm-04-01-tree` requires.

:::{prf:theorem} Correctness of Prim's algorithm
:label: thm-04-01-prim

On a connected network, Prim's algorithm returns a minimum spanning tree,
whatever the start node and the tie-breaking.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

We show by induction that, at the start of every iteration, the current edge
set $\mathcal T$ is contained in some minimum spanning tree. This holds for
$\mathcal T=\varnothing$. In iteration $k$ every edge of $\mathcal T$ has both
ends in $C_k$, so no edge of $\mathcal T$ lies in the cut $[C_k,\bar C_k]$, and
the selected edge is a cheapest edge of that cut. The cut property
({prf:ref}`thm-04-01-cut-property` with $F=\mathcal T$, $S=C_k$) shows that
$\mathcal T\cup\{(i^*,j^*)\}$ is again contained in a minimum spanning tree.
The cut is never empty because $G$ is connected. At the end $\mathcal T$ is a
spanning tree (connected, $\abs{\mathcal N}-1$ edges) contained in a minimum
spanning tree, which also has $\abs{\mathcal N}-1$ edges; the two are equal.
:::

**Running time.** Searching the whole cut in every iteration costs
$O(\abs{\mathcal A})$ per iteration and $O(\abs{\mathcal N}\abs{\mathcal A})$
in total. It is better to keep, for every unconnected node $j$, a label
$[\gamma(j),P(j)]$: the cost $\gamma(j)$ of the cheapest edge from $j$ to the
connected set and the node $P(j)\in C_k$ at its other end ($\gamma(j)=\infty$
if there is none). Step 3 then picks the unconnected node with the smallest
$\gamma$, and after $j^*$ joins, only the edges at $j^*$ need checking: for
each neighbor $j\in\bar C_{k+1}$ with $c_{j^*j}<\gamma(j)$, set
$\gamma(j)=c_{j^*j}$ and $P(j)=j^*$. With the labels in an array, each
selection takes $O(\abs{\mathcal N})$ and the whole algorithm
$O(\abs{\mathcal N}^2)$, the best possible for dense networks. With the labels
in a binary heap the time is $O(\abs{\mathcal A}\log\abs{\mathcal N})$, and
with a Fibonacci heap $O(\abs{\mathcal A}+\abs{\mathcal N}\log\abs{\mathcal N})$
(Fredman & Tarjan, 1987). The label scheme is the same as in
Dijkstra's shortest-path algorithm ({prf:ref}`alg-04-02-dijkstra`); only the update
differs, $c_{j^*j}$ here against $V(j^*)+c_{j^*j}$ there.

:::{prf:example} Prim's algorithm on the sewage network
:label: ex-04-01-sewage-prim

Apply Prim's algorithm to the sewage network of {prf:ref}`ex-04-01-sewage`,
starting from the main station, node 7. When several links tie for the least
cost, take the one that brings in the lower-numbered station, and then the
one whose end in $C_k$ has the lower number.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

{numref}`tab-04-01-prim` lists the iterations. In iteration 1 the cut
$[\{7\},\{1,\dots,6\}]$ contains 5–7 (3), 6–7 (4) and 3–7 (8); the cheapest
is 5–7. In iteration 2 the cut of $C_2=\{5,7\}$ contains 1–5 (9), 2–5 (3),
4–5 (8), 3–7 (8) and 6–7 (4), and 2–5 enters. In iteration 3, with
$C_3=\{2,5,7\}$ ({numref}`fig-04-01-prim-steps`a), the cut contains 1–2 (1),
2–3 (6), 2–4 (4), 1–5 (9), 4–5 (8), 3–7 (8) and 6–7 (4), and 1–2 enters. In
iteration 4 two links of cost 4 tie, 2–4 and 6–7; the rule picks 2–4 (it
brings in station 4). In iteration 5 ({numref}`fig-04-01-prim-steps`b) the
cheapest link of the cut is 4–6 (3). In iteration 6 only station 3 is left;
its links to the connected set are 1–3 (5), 3–4 (5), 2–3 (6), 3–7 (8) and
3–6 (10), and the tie between 1–3 and 3–4 is broken in favor of 1–3.

```{list-table} Prim's algorithm on the sewage network, started at node 7 (costs in \$ million)
:name: tab-04-01-prim
:header-rows: 1

* - Iteration $k$
  - Connected set $C_k$
  - Cheapest links in $[C_k,\bar C_k]$
  - Link added
  - Cost so far
* - 1
  - $\{7\}$
  - 5–7 (3)
  - 5–7
  - 3
* - 2
  - $\{5,7\}$
  - 2–5 (3)
  - 2–5
  - 6
* - 3
  - $\{2,5,7\}$
  - 1–2 (1)
  - 1–2
  - 7
* - 4
  - $\{1,2,5,7\}$
  - 2–4 (4), 6–7 (4)
  - 2–4
  - 11
* - 5
  - $\{1,2,4,5,7\}$
  - 4–6 (3)
  - 4–6
  - 14
* - 6
  - $\{1,2,4,5,6,7\}$
  - 1–3 (5), 3–4 (5)
  - 1–3
  - 19
```

After six iterations $\bar C_7=\varnothing$. The tunnels to build are

$$
\begin{gathered}
\mathcal T^*=\{5\text{–}7,\ 2\text{–}5,\ 1\text{–}2,\ 2\text{–}4,\ 4\text{–}6,\ 1\text{–}3\},\\
c(\mathcal T^*)=3+3+1+4+3+5=19,
\end{gathered}
$$

a total of \$19 million ({numref}`fig-04-01-sewage-mst`), \$11 million less than
the feasible layout of {prf:ref}`ex-04-01-sewage`. The two ties already show that
the optimal layout need not be unique: choosing 6–7 in iteration 4 or 3–4 in
iteration 6 leads to other trees of the same cost (examined in
{prf:ref}`ex-04-01-sewage-kruskal`).

```{figure} figures/04-01-prim-steps.svg
:name: fig-04-01-prim-steps
:width: 560px
:class: dark-light
:alt: Two copies of the sewage network with the connected set filled, tree links thick and the links of the cut dashed

Prim's algorithm on the sewage network. Filled nodes form the connected set,
thick links the tree so far, dashed links the cut $[C_k,\bar C_k]$ with their
costs; gray links are not in the cut. (a) $C_3=\{2,5,7\}$: link 1–2 (cost 1)
enters next. (b) $C_5=\{1,2,4,5,7\}$: link 4–6 (cost 3) enters next.
```

```{figure} figures/04-01-sewage-mst.svg
:name: fig-04-01-sewage-mst
:width: 560px
:class: dark-light
:alt: The sewage network with six thick tree links and two dashed alternative links

The minimum spanning tree of the sewage network (thick, total cost \$19
million). The dashed links 6–7 and 3–4 tie with the tree links 2–4 and 1–3 on
their cycles; exchanging them gives the other three optimal trees.
```
:::

(sec-04-01-kruskal)=
## Kruskal's algorithm

Kruskal's algorithm (Kruskal, 1956) looks at the edges instead of the
nodes: it considers them from the cheapest to the most expensive and keeps an
edge unless it would close a cycle with the edges already kept. In between,
the kept edges form a forest whose components grow and merge.

:::{prf:algorithm} Kruskal's algorithm
:label: alg-04-01-kruskal

**Input:** an undirected network $G=(\mathcal N,\mathcal A)$ with costs
$c_{ij}$.
**Output:** a minimum spanning tree $\mathcal T$ (a minimum spanning forest if
$G$ is not connected).

1. **Sort.** List the edges in nondecreasing order of cost,
   $e_1,e_2,\dots,e_{\abs{\mathcal A}}$. Set $\mathcal T=\varnothing$ (every
   node is a component of its own) and $k=1$.
2. **Stopping test.** If $\abs{\mathcal T}=\abs{\mathcal N}-1$ or
   $k>\abs{\mathcal A}$, stop.
3. **Test the edge.** Let $e_k=(i,j)$. If $i$ and $j$ lie in different
   components of $(\mathcal N,\mathcal T)$, set
   $\mathcal T\leftarrow\mathcal T\cup\{e_k\}$ and merge the two components;
   otherwise reject $e_k$, which would close a cycle.
4. **Next edge.** Set $k\leftarrow k+1$ and go to Step 2.
:::

:::{prf:theorem} Correctness of Kruskal's algorithm
:label: thm-04-01-kruskal

On a connected network, Kruskal's algorithm returns a minimum spanning tree.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Accepted edges never close a cycle, so $\mathcal T$ stays a forest. Suppose
$e_k=(i,j)$ is accepted, and let $S$ be the component of $i$ just before. No
edge of $\mathcal T$ lies in the cut $[S,\bar S]$. Any edge of the cut cheaper
than $e_k$ was examined earlier; at that time its ends were in different
components too (components only merge), so it would have been accepted and
would now be an edge of $\mathcal T$ in the cut, which is impossible. Hence
$e_k$ is a cheapest edge of the cut, and by induction with the cut property
({prf:ref}`thm-04-01-cut-property`) $\mathcal T$ is always contained in a minimum
spanning tree. If the algorithm stopped with fewer than $\abs{\mathcal N}-1$
edges, some edge of $G$ would join two of the final components, and that edge
would have been accepted when examined. So it ends with a spanning tree.
:::

**Running time and union–find.** Sorting takes
$O(\abs{\mathcal A}\log\abs{\mathcal A})=O(\abs{\mathcal A}\log
\abs{\mathcal N})$ time. Step 3 needs to know whether two nodes are in the same
component, and to merge components. A simple **union–find** scheme stores for
every node the name of its component: the test compares two names, and a merge
renames the nodes of the smaller component. A node is renamed only when its
component at least doubles in size, so at most $\log_2\abs{\mathcal N}$ times,
and all renaming costs $O(\abs{\mathcal N}\log\abs{\mathcal N})$. Sorting
dominates, and Kruskal's algorithm runs in $O(\abs{\mathcal A}\log
\abs{\mathcal N})$ time. It is the natural choice for sparse networks given as
arc lists, while the array version of Prim's algorithm suits dense networks.
(With union–find trees and path compression the merging becomes almost linear;
see Ahuja et al., 1993, Section 13.4.)

:::{prf:example} Kruskal's algorithm and alternative optima for the sewage network
:label: ex-04-01-sewage-kruskal

Apply Kruskal's algorithm to the sewage network of {prf:ref}`ex-04-01-sewage`,
taking tied links in increasing order of their end nodes, and check whether
the optimal layout is unique.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The sorted list begins 1–2 (1), 2–5 (3), 4–6 (3), 5–7 (3), 2–4 (4), 6–7 (4),
1–3 (5), 3–4 (5). {numref}`tab-04-01-kruskal` shows the decisions. After 2–4 is
accepted, nodes 6 and 7 are already in the same component, so 6–7 would close
the cycle 6–7–5–2–4–6 and is rejected. Accepting 1–3 completes six links,
and the algorithm stops without examining the remaining seven links.

```{list-table} Kruskal's algorithm on the sewage network (costs in \$ million)
:name: tab-04-01-kruskal
:header-rows: 1

* - Link (cost)
  - Decision
  - Components after the step
* - 1–2 (1)
  - accept
  - $\{1,2\},\{3\},\{4\},\{5\},\{6\},\{7\}$
* - 2–5 (3)
  - accept
  - $\{1,2,5\},\{3\},\{4\},\{6\},\{7\}$
* - 4–6 (3)
  - accept
  - $\{1,2,5\},\{3\},\{4,6\},\{7\}$
* - 5–7 (3)
  - accept
  - $\{1,2,5,7\},\{3\},\{4,6\}$
* - 2–4 (4)
  - accept
  - $\{1,2,4,5,6,7\},\{3\}$
* - 6–7 (4)
  - reject (cycle)
  - unchanged
* - 1–3 (5)
  - accept
  - $\{1,2,3,4,5,6,7\}$
```

The accepted links are exactly the tree $\mathcal T^*$ found by Prim's
algorithm, with cost $1+3+3+3+4+5=19$, that is, \$19 million.

*Optimality and uniqueness.* {numref}`tab-04-01-sewage-check` applies
{prf:ref}`thm-04-01-mst-optimality`: for every link outside $\mathcal T^*$ it gives
the tree path between its ends and the most expensive link on that path.

```{list-table} Non-tree links of the sewage network and their tree paths
:name: tab-04-01-sewage-check
:header-rows: 1

* - Non-tree link (cost)
  - Tree path
  - Largest cost on the path
  - Comparison
* - 1–4 (7)
  - 1–2–4
  - 4
  - $7>4$
* - 1–5 (9)
  - 1–2–5
  - 3
  - $9>3$
* - 2–3 (6)
  - 2–1–3
  - 5
  - $6>5$
* - 3–4 (5)
  - 3–1–2–4
  - 5 (link 1–3)
  - tie
* - 3–6 (10)
  - 3–1–2–4–6
  - 5
  - $10>5$
* - 3–7 (8)
  - 3–1–2–5–7
  - 5
  - $8>5$
* - 4–5 (8)
  - 4–2–5
  - 4
  - $8>4$
* - 6–7 (4)
  - 6–4–2–5–7
  - 4 (link 2–4)
  - tie
```

No non-tree link is cheaper than a link on its tree path, so $\mathcal T^*$ is
optimal. There are two ties, so by {prf:ref}`cor-04-01-unique` the optimum is not
unique: 6–7 can replace 2–4, and 3–4 can replace 1–3. The two exchanges are
independent (after the first, the tree path of 3–4 is 3–1–2–5–7–6–4 and still
contains 1–3 as its most expensive link), so there are four optimal layouts,
each costing \$19 million:

$$
\begin{gathered}
\{1\text{–}2,\ 1\text{–}3,\ 2\text{–}4,\ 2\text{–}5,\ 4\text{–}6,\ 5\text{–}7\},\\
\{1\text{–}2,\ 1\text{–}3,\ 2\text{–}5,\ 4\text{–}6,\ 5\text{–}7,\ 6\text{–}7\},\\
\{1\text{–}2,\ 2\text{–}4,\ 2\text{–}5,\ 3\text{–}4,\ 4\text{–}6,\ 5\text{–}7\},\\
\{1\text{–}2,\ 2\text{–}5,\ 3\text{–}4,\ 4\text{–}6,\ 5\text{–}7,\ 6\text{–}7\}.
\end{gathered}
$$

A check of all 1 020 spanning trees confirms that these four, and no others,
cost \$19 million. Prim's algorithm reaches each of them with a suitable
tie-breaking, and Kruskal's algorithm does so when the tied links are listed in
another order.
:::

(sec-04-01-boruvka)=
## Borůvka's algorithm

The oldest minimum spanning tree method, published by Borůvka (1926) for the
design of an electricity network, works on all components at once. In each
round every component of the current forest (at the start, every node) selects
a cheapest edge leaving it, and all selected edges are added together; each
of them is justified by the cut property for its own component. Ties must be
broken by one fixed rule, for example by edge index, so that the selected edges
cannot form a cycle. Each round at least halves the number of components, so
there are at most $\lceil\log_2\abs{\mathcal N}\rceil$ rounds of
$O(\abs{\mathcal A})$ work, $O(\abs{\mathcal A}\log\abs{\mathcal N})$ in total.
On the sewage network, with ties broken towards the lower-numbered other end,
the first round selects 1–2, 1–3, 2–5, 4–6 and 5–7 and leaves the components
$\{1,2,3,5,7\}$ and $\{4,6\}$; the second round joins them with 2–4, giving the
\$19 million tree once more. Because every component acts independently, the
method is the basis of parallel and distributed MST algorithms.

## Summary and Common Mistakes

### Summary

- A network $G=(\mathcal N,\mathcal A)$ has nodes and arcs (directed) or edges
  (undirected); paths, chains, directed paths, cycles and connectivity are
  defined in {prf:ref}`def-04-01-path` and {prf:ref}`def-04-01-connected`, and the
  degrees satisfy {eq}`eq-04-01-handshake`.
- A tree is a connected network without cycles; equivalently it is connected
  with $\abs{\mathcal N}-1$ edges, acyclic with $\abs{\mathcal N}-1$ edges, or
  has a unique path between every two nodes ({prf:ref}`thm-04-01-tree`). Adding a
  non-tree edge creates one cycle and deleting a tree edge creates one cut
  ({prf:ref}`cor-04-01-exchange`).
- Networks are stored as arc lists, node–arc incidence matrices (rank
  $\abs{\mathcal N}-1$ when connected, {numref}`tab-04-01-incidence`), adjacency
  or length matrices, or adjacency lists.
- The minimum spanning tree problem {eq}`eq-04-01-tree-cost` asks for a cheapest
  spanning tree. A cheapest edge of any cut belongs to some optimal tree
  ({prf:ref}`thm-04-01-cut-property`); a strictly most expensive edge of a cycle
  belongs to none ({prf:ref}`thm-04-01-cycle-property`); a tree is optimal if and
  only if no non-tree edge is cheaper than an edge on its tree path
  ({prf:ref}`thm-04-01-mst-optimality`), and ties decide uniqueness
  ({prf:ref}`cor-04-01-unique`).
- Prim's algorithm ({prf:ref}`alg-04-01-prim`) grows one tree in
  $\abs{\mathcal N}-1$ iterations, in $O(\abs{\mathcal N}^2)$ time with labels
  or $O(\abs{\mathcal A}\log\abs{\mathcal N})$ with a heap; Kruskal's algorithm
  ({prf:ref}`alg-04-01-kruskal`) scans the sorted edges with union–find in
  $O(\abs{\mathcal A}\log\abs{\mathcal N})$; Borůvka's algorithm merges all
  components in parallel rounds.
- For the sewage network both algorithms give a \$19 million layout
  ({numref}`fig-04-01-sewage-mst`), one of four optimal trees among 1 020
  spanning trees.

### Common mistakes

- **Treating a minimum spanning tree as a set of shortest routes.** The tree
  minimizes the total cost of the links, not the distance from one node to the
  others. In the sewage network station 3 reaches node 7 along 3–1–2–5–7, at a
  cost of 12, although the direct link 3–7 costs 8. Shortest routes are the
  subject of {doc}`04-02-shortest-paths` (compare {prf:ref}`prop-04-02-tree`).
- **Looking only at the last node added.** Prim's algorithm selects the
  cheapest edge out of the whole connected set $C_k$, not out of the node added
  most recently; the latter "nearest neighbor" rule can miss cheap edges at
  earlier nodes.
- **Rejecting an edge in Kruskal's algorithm because one of its ends is already
  in the tree.** An edge is rejected only when both ends are in the same
  component. Early on the forest has many components, and an edge touching a
  connected node may still join two of them.
- **Claiming uniqueness after one run.** A single run returns one optimal
  tree. Whether it is the only one follows from
  {prf:ref}`cor-04-01-unique`: look for non-tree edges that tie with an edge on
  their tree path, as in {prf:ref}`ex-04-01-sewage-kruskal`.
- **Counting iterations or edges wrongly.** A spanning tree on
  $\abs{\mathcal N}$ nodes has $\abs{\mathcal N}-1$ edges, so Prim's algorithm
  runs $\abs{\mathcal N}-1$ iterations and Kruskal's algorithm can stop after
  $\abs{\mathcal N}-1$ acceptances.
- **Applying the algorithms to directed networks.** Prim's and Kruskal's
  algorithms ignore arc directions; a cheapest structure that reaches every
  node from a root along directed arcs is a different problem.

## References

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows: Theory,
  algorithms, and applications*. Prentice Hall. (Chapter 2: notation,
  definitions and network representations; Chapter 13: minimum spanning trees,
  optimality conditions, Kruskal's, Prim's and Sollin's algorithms.)
- Borůvka, O. (1926). O jistém problému minimálním [On a certain minimal
  problem]. *Práce Moravské přírodovědecké společnosti, 3*, 37–58.
- Fredman, M. L., & Tarjan, R. E. (1987). Fibonacci heaps and their uses in
  improved network optimization algorithms. *Journal of the ACM, 34*(3),
  596–615.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Sections 10.2 and 10.4: terminology of
  networks and the minimum spanning tree problem.)
- Jarník, V. (1930). O jistém problému minimálním [About a certain minimal
  problem]. *Práce Moravské přírodovědecké společnosti, 6*(4), 57–63.
- Kruskal, J. B. (1956). On the shortest spanning subtree of a graph and the
  traveling salesman problem. *Proceedings of the American Mathematical
  Society, 7*(1), 48–50. <https://doi.org/10.1090/S0002-9939-1956-0078686-7>
- Prim, R. C. (1957). Shortest connection networks and some generalizations.
  *Bell System Technical Journal, 36*(6), 1389–1401.
  <https://doi.org/10.1002/j.1538-7305.1957.tb01515.x>
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Sections 6.1–6.2: scope and definition of network models, minimal
  spanning tree algorithm.)
