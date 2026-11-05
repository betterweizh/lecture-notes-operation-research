# Network Optimization

This chapter works through the examples and problems of the four chapters
on network optimization: {doc}`04-01-networks-mst`,
{doc}`04-02-shortest-paths`, {doc}`04-03-maximum-flow` and
{doc}`04-04-minimum-cost-flow`.

Every item states its data in full (networks as arc lists, most of them also
drawn) and is followed by a folded solution. Nodes of flow
networks are named so that they do not clash with the book's symbols (no
node is called $u$, $v$ or $x$).

(sec-08-05-terminology)=
## Terminology and representation

:::{prf:example} Paths, cycles and trees in two directed networks
:label: ex-08-05-directed-terminology

The two directed networks of {numref}`fig-08-05-directed-networks` are

- (i) nodes $1,\dots,5$, arcs $(1,2)$, $(2,5)$, $(5,1)$, $(1,3)$, $(3,4)$,
  $(3,5)$, $(4,5)$, $(4,2)$;
- (ii) nodes $1,\dots,4$, arcs $(1,3)$, $(1,2)$, $(2,3)$, $(3,4)$, $(2,4)$.

For each network give (a) a directed path, (b) a directed cycle, (c) a tree
and (d) a spanning tree.

```{figure} figures/08-05-directed-networks.svg
:name: fig-08-05-directed-networks
:width: 520px
:class: dark-light
:alt: Two small directed networks, a five-node network on the left and a four-node network on the right

The directed networks (i) and (ii) of {prf:ref}`ex-08-05-directed-terminology`.
```
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The definitions are those of {prf:ref}`def-04-01-path` and
{prf:ref}`def-04-01-tree`: a directed path follows
the arcs in their direction without repeating a node; a directed cycle is a
directed path that returns to its first node; a tree is a connected
subnetwork without cycles (directions ignored); a spanning tree is a tree
that contains every node.

**(a) Directed path.** (i) $1\to3\to4\to5$; (ii) $1\to2\to3\to4$. Each step
uses an arc in its own direction.

**(b) Directed cycle.** (i) $1\to3\to5\to1$ uses the arcs $(1,3)$, $(3,5)$,
$(5,1)$. Network (i) has four directed cycles in all: $1\to2\to5\to1$,
$1\to3\to5\to1$, $1\to3\to4\to5\to1$ and $1\to3\to4\to2\to5\to1$; each passes
through the arc $(5,1)$, the only arc entering node 1. (ii) There is **no**
directed cycle. Every arc goes from a smaller to a larger node number
($1<2<3<4$ is a topological order), so no directed path can come back.

**(c) Tree.** (i) the arcs $(1,3)$, $(3,4)$, $(3,5)$ on the nodes
$\{1,3,4,5\}$; (ii) the arcs $(1,3)$, $(1,2)$ on $\{1,2,3\}$. Both are
connected and have one arc fewer than nodes, hence no cycle.

**(d) Spanning tree.** (i) $(1,2)$, $(1,3)$, $(3,4)$, $(3,5)$: 4 arcs on all
5 nodes; (ii) $(1,3)$, $(1,2)$, $(2,4)$: 3 arcs on all 4 nodes. Each is
connected without cycles. (By the matrix-tree theorem of
{ref}`sec-04-01-trees`, network (i) has 45
spanning trees and network (ii) has 8, so many other answers are correct.)
:::

:::{prf:example} Drawing a network from its definition
:label: ex-08-05-draw-network

Draw the directed network with node set $\mathcal N=\{1,2,3,4,5\}$ and arc
set

$$
\mathcal A=\{(1,2),(1,5),(2,3),(2,4),(3,4),(3,5),(4,3),(4,5),(5,2)\}.
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Place the five nodes and draw one arrow per ordered pair, from the first node
to the second ({numref}`fig-08-05-draw-network`). The pairs $(3,4)$ and $(4,3)$
are two different arcs, drawn as two curved arrows between 3 and 4. A check:
the out-degrees are $2,2,2,2,1$ and the in-degrees $0,2,2,2,3$ for nodes
$1,\dots,5$; both add up to $\abs{\mathcal A}=9$. Node 1 has no entering arc,
so it can only be the start of a path.

```{figure} figures/08-05-draw-network.svg
:name: fig-08-05-draw-network
:width: 360px
:class: dark-light
:alt: Five-node directed network with nine arcs including the opposite pair between nodes 3 and 4

The network of {prf:ref}`ex-08-05-draw-network`.
```
:::

:::{prf:example} Placing numbers in squares with a network
:label: ex-08-05-squares-puzzle

Eight equal squares are arranged in three rows: two in the first row, four
in the second and two in the third, each row centered under the one above.
Name them A, B (row 1), C, D, E, F (row 2), G, H (row 3) from left to right,
so that A and B sit above D and E, and G and H below them. Write the numbers
$1,\dots,8$ into the squares, each once, so that no two squares that touch
horizontally, vertically or diagonally hold consecutive numbers. Use a
network to find the solutions systematically.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Network.** Let every square be a node and join two nodes by an edge when
the squares touch (including diagonally). This gives 17 edges:

A–B, A–C, A–D, A–E, B–D, B–E, B–F, C–D, C–G, D–E, D–G, D–H, E–F, E–G, E–H,
F–H, G–H.

The degrees are $\deg(D)=\deg(E)=6$, $\deg(A)=\deg(B)=\deg(G)=\deg(H)=4$,
$\deg(C)=\deg(F)=3$. The task is to label the nodes with $1,\dots,8$ so that
adjacent nodes never get consecutive labels.

**The two central squares.** D is adjacent to every square except F, and E to
every square except C. The number in D has one or two consecutive numbers
($k-1$ and $k+1$), and each must go into a square not adjacent to D; there is
only one such square, F. So the number in D has only one consecutive number:
it is 1 or 8. The same holds for E. Take $D=1$, $E=8$ (the case $D=8$, $E=1$
follows by replacing every number $k$ by $9-k$).

**Forced squares.** The number 2 must avoid all squares adjacent to D, so
$F=2$; the number 7 must avoid all squares adjacent to E, so $C=7$.

**The corners.** The numbers 3, 4, 5, 6 remain for A, B, G, H. A–B and G–H
are edges, and A, B are not adjacent to G, H, so $\{A,B\}$ and $\{G,H\}$ must
each be a nonconsecutive pair. The only split of $\{3,4,5,6\}$ into two
nonconsecutive pairs is $\{3,5\}$ and $\{4,6\}$. Further, B touches $F=2$,
so $B\ne3$; H touches $F=2$, so $H\ne3$; A and G touch $C=7$, so
$A\ne6$ and $G\ne6$. This leaves two possibilities:

- $\{A,B\}=\{3,5\}$: $A=3$, $B=5$, and then $G=4$, $H=6$;
- $\{A,B\}=\{4,6\}$: $A=4$, $B=6$, and then $G=3$, $H=5$.

**Answer.** With $D=1$, $E=8$ there are exactly two solutions, and the
mirror case gives two more, four in all (a complete enumeration of the
$8!$ placements confirms this). One of them is

```{list-table}
:header-rows: 0
:widths: 1 1 1 1

* -
  - 3
  - 5
  -
* - 7
  - 1
  - 8
  - 2
* -
  - 4
  - 6
  -
```

Check: the neighbors of 1 (squares A, C, E, G, H, B) hold 3, 7, 8, 4, 6, 5,
none equal to 2; the neighbors of 8 hold 3, 5, 1, 2, 4, 6, none equal to 7.
:::

(sec-08-05-mst)=
## Minimum spanning trees and clustering

:::{prf:example} Minimum spanning tree of an eight-node network
:label: ex-08-05-mst-eight

An undirected network has nodes $s, v_1,\dots,v_7$ and the edges (length in
km)

| Edge | $s$–$v_1$ | $s$–$v_2$ | $s$–$v_7$ | $v_1$–$v_2$ | $v_1$–$v_7$ | $v_2$–$v_3$ |
|---|---|---|---|---|---|---|
| Length | 3 | 6 | 1 | 2 | 2 | 1 |

| Edge | $v_3$–$v_4$ | $v_3$–$v_5$ | $v_3$–$v_6$ | $v_4$–$v_5$ | $v_5$–$v_6$ | $v_6$–$v_7$ |
|---|---|---|---|---|---|---|
| Length | 6 | 3 | 4 | 3 | 2 | 10 |

Find a minimum spanning tree with Prim's algorithm, starting from $s$, and
decide whether it is unique.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Prim's algorithm ({prf:ref}`alg-04-01-prim`) adds, at every step, the shortest
edge between the connected set $C$ and the other nodes.

| Step | Connected set $C$ before the step | Shortest edge leaving $C$ | Length | Total |
|---|---|---|---|---|
| 1 | $\{s\}$ | $s$–$v_7$ | 1 | 1 |
| 2 | $\{s,v_7\}$ | $v_7$–$v_1$ | 2 | 3 |
| 3 | $\{s,v_7,v_1\}$ | $v_1$–$v_2$ | 2 | 5 |
| 4 | $\{s,v_7,v_1,v_2\}$ | $v_2$–$v_3$ | 1 | 6 |
| 5 | $\ldots\cup\{v_3\}$ | $v_3$–$v_5$ | 3 | 9 |
| 6 | $\ldots\cup\{v_5\}$ | $v_5$–$v_6$ | 2 | 11 |
| 7 | $\ldots\cup\{v_6\}$ | $v_5$–$v_4$ | 3 | 14 |

For example, in step 2 the edges leaving $\{s,v_7\}$ are $s$–$v_1$ (3),
$s$–$v_2$ (6), $v_7$–$v_1$ (2) and $v_7$–$v_6$ (10); the shortest is
$v_7$–$v_1$. In step 7 the remaining node $v_4$ can be reached by
$v_5$–$v_4$ (3) or $v_3$–$v_4$ (6).

The minimum spanning tree has the 7 edges $s$–$v_7$, $v_7$–$v_1$,
$v_1$–$v_2$, $v_2$–$v_3$, $v_3$–$v_5$, $v_5$–$v_6$, $v_4$–$v_5$ and length
**14 km** ({numref}`fig-08-05-mst-eight`).

**Uniqueness.** Each of the five non-tree edges is strictly longer than every
tree edge on the cycle it closes: $s$–$v_1$ (3) against $1,2$; $s$–$v_2$ (6)
against $1,2,2$; $v_3$–$v_4$ (6) against $3,3$; $v_3$–$v_6$ (4) against
$3,2$; $v_6$–$v_7$ (10) against $2,2,1,3,2$. By {prf:ref}`cor-04-01-unique` the
tree is the only minimum spanning tree (an enumeration of all spanning trees finds
exactly one of length 14). Kruskal's algorithm gives the same tree.

```{figure} figures/08-05-mst-eight.svg
:name: fig-08-05-mst-eight
:width: 460px
:class: dark-light
:alt: Eight-node undirected network with the seven edges of the minimum spanning tree highlighted

The network of {prf:ref}`ex-08-05-mst-eight`; the minimum spanning tree (14 km)
is highlighted.
```
:::

:::{prf:example} Gas pipelines with two pressure groups
:label: ex-08-05-wellheads

Nine gas wellheads are to be linked by pipelines to a delivery point, which
is connected to wellhead 1. The possible links and their lengths in miles are

| Link | 1–2 | 1–3 | 1–4 | 1–5 | 1–8 | 1–9 | 2–3 | 3–4 | 3–5 | 4–5 |
|---|---|---|---|---|---|---|---|---|---|---|
| Length | 5 | 9 | 20 | 4 | 14 | 15 | 6 | 15 | 10 | 20 |

| Link | 4–6 | 4–7 | 5–6 | 5–7 | 5–8 | 5–9 | 7–8 | 8–9 | 3–7 |
|---|---|---|---|---|---|---|---|---|---|
| Length | 7 | 12 | 3 | 5 | 13 | 6 | 7 | 5 | 8 |

The wells form two groups by gas pressure: a high-pressure group
$\{2,3,4,7\}$ and a low-pressure group $\{5,6,8,9\}$. Wells of different
groups may not be linked directly, and each group must reach the delivery
point through wellhead 1. Find the shortest pipeline network.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Decomposition.** Remove the five links that join the two groups: 3–5, 4–5,
4–6, 5–7 and 7–8. What remains splits into two sub-networks that share only
node 1: the high-pressure part on $\{1,2,3,4,7\}$ and the low-pressure part on
$\{1,5,6,8,9\}$. A feasible pipeline network connects every well to node 1
without cross links, so it is the union of a spanning tree of each part, and
the shortest one consists of the two minimum spanning trees.

**High-pressure part** (links 1–2 5, 1–3 9, 1–4 20, 2–3 6, 3–4 15, 4–7 12,
3–7 8). Prim from node 1: 1–2 (5), 2–3 (6), 3–7 (8), 7–4 (12); total
**31 miles**.

**Low-pressure part** (links 1–5 4, 1–8 14, 1–9 15, 5–6 3, 5–8 13, 5–9 6,
8–9 5). Prim from node 1: 1–5 (4), 5–6 (3), 5–9 (6), 9–8 (5); total
**18 miles**.

**Answer.** Lay the links 1–2, 2–3, 3–7, 7–4 and 1–5, 5–6, 5–9, 9–8, in all
$31+18=49$ miles, plus the link from wellhead 1 to the delivery point
({numref}`fig-08-05-wellheads`). Without the pressure restriction the minimum
spanning tree would be 41 miles (it uses the cross links 4–6 and 5–7), so
the separation costs 8 miles of pipe.

```{figure} figures/08-05-wellheads.svg
:name: fig-08-05-wellheads
:width: 480px
:class: dark-light
:alt: Wellhead network with the two group trees highlighted and the forbidden cross links drawn dashed

The wellhead network of {prf:ref}`ex-08-05-wellheads` (lengths in miles,
delivery point P). Highlighted: the two minimum spanning trees, 49 miles in
total; dashed: links between the pressure groups, which may not be used.
```
:::

:::{prf:example} Machine cells from a minimum spanning tree
:label: ex-08-05-machine-cells

A plant makes 15 parts on 10 machines. Machine $i$ processes the set of
parts $\mathcal P_i$:

| Machine | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| Parts | 1, 6 | 2, 3, 7, 8, 9, 12, 13, 15 | 3, 5, 10, 14 | 2, 7, 8, 11, 12, 13 | 3, 5, 10, 11, 14 |

| Machine | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|
| Parts | 1, 4, 5, 9, 10 | 2, 5, 7, 8, 9, 10 | 3, 4, 15 | 4, 10 | 3, 8, 10, 14, 15 |

The machines are to be grouped into cells so that machines in the same cell
process similar parts. The dissimilarity of machines $i$ and $j$ is
$d_{ij}=1-n_{ij}/(n_{ij}+m_{ij})$, where $n_{ij}$ is the number of parts
processed on both machines and $m_{ij}$ the number processed on exactly one
of them.

(a) Formulate a network model (nodes, links, link attributes).
(b) Explain why the cells can be obtained from a minimum spanning tree.
(c) Find the two-cell and three-cell solutions.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** One node per machine, an undirected link between every pair of
machines (a complete network with 45 links), and the dissimilarity $d_{ij}$
as the length of link $i$–$j$. Since $n_{ij}+m_{ij}=\abs{\mathcal P_i\cup
\mathcal P_j}$, $d_{ij}=1-\abs{\mathcal P_i\cap\mathcal P_j}/\abs{\mathcal
P_i\cup\mathcal P_j}$. For machines 2 and 4, for example, the shared parts are
2, 7, 8, 12, 13 ($n=5$) and the union has 9 parts, so $d_{24}=1-5/9=0.4444$.
All lengths (4 decimals):

| | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| **1** | 1 | 1 | 1 | 1 | 0.8333 | 1 | 1 | 1 | 1 |
| **2** | | 0.9091 | 0.4444 | 0.9167 | 0.9167 | 0.6000 | 0.7778 | 1 | 0.7000 |
| **3** | | | 1 | 0.2000 | 0.7143 | 0.7500 | 0.8333 | 0.8000 | 0.5000 |
| **4** | | | | 0.9000 | 1 | 0.6667 | 1 | 1 | 0.9000 |
| **5** | | | | | 0.7500 | 0.7778 | 0.8571 | 0.8333 | 0.5714 |
| **6** | | | | | | 0.6250 | 0.8571 | 0.6000 | 0.8889 |
| **7** | | | | | | | 1 | 0.8571 | 0.7778 |
| **8** | | | | | | | | 0.7500 | 0.6667 |
| **9** | | | | | | | | | 0.8333 |

**(b)** For a partition into $k$ cells, call the smallest dissimilarity
between two machines in *different* cells its **spacing**; good cells have a
large spacing. Build a minimum spanning tree and delete its $k-1$ longest
links; the $k$ remaining components are the cells. This partition has the
largest possible spacing. Indeed, let $d^*$ be the length of the
$(k-1)$-th longest tree link. Two machines in different cells of the tree
partition are joined in the tree by a path that contains a deleted link, of
length at least $d^*$; by the optimality condition of a minimum spanning tree
({prf:ref}`thm-04-01-mst-optimality`),
the link between them is at least as long as every tree link on that path,
hence at least $d^*$. So the spacing is $d^*$ (the deleted link of length $d^*$
itself joins two cells). For any other partition
into $k$ cells, some two machines $p$, $q$ lie in one tree cell but in
different cells of the other partition. The tree path from $p$ to $q$ uses
only kept links, all of length at most $d^*$, and somewhere it steps from
one cell of the other partition to another; so that partition has spacing at
most $d^*$. Within each cell, the tree links connect the machines with the
least total dissimilarity. (This is single-linkage clustering.)

**(c)** Kruskal's algorithm ({prf:ref}`alg-04-01-kruskal`) accepts, in order,

| Link | 3–5 | 2–4 | 3–10 | 2–7 | 6–9 |
|---|---|---|---|---|---|
| $d_{ij}$ | 0.2000 | 0.4444 | 0.5000 | 0.6000 | 0.6000 |

| Link | 6–7 | 8–10 | 2–10 | 1–6 |
|---|---|---|---|---|
| $d_{ij}$ | 0.6250 | 0.6667 | 0.7000 | 0.8333 |

(links such as 5–10 (0.5714) and 4–7 (0.6667) are skipped because they would
close a cycle).
The minimum spanning tree has total length **5.1694**
({numref}`fig-08-05-machine-cells`).

- **Two cells:** delete the longest link 1–6 (0.8333): cells $\{1\}$ and
  $\{2,3,4,5,6,7,8,9,10\}$.
- **Three cells:** delete also the next longest, 2–10 (0.7000): cells
  $\{1\}$, $\{2,4,6,7,9\}$ and $\{3,5,8,10\}$.

Machine 1 shares only part 1 with machine 6 and nothing with the others, so
it forms a cell of its own in both solutions.

*(An earlier solution draws a spanning tree with the link 9–8 (0.75) in place
of 2–10 (0.70); its length is 5.2194, so it is not minimal. Its two- and
three-cell partitions happen to be the same.)*

```{figure} figures/08-05-machine-cells.svg
:name: fig-08-05-machine-cells
:width: 520px
:class: dark-light
:alt: Minimum spanning tree on ten machines with dissimilarities, two deleted links dashed, giving three cells

Minimum spanning tree of {prf:ref}`ex-08-05-machine-cells` with dissimilarities.
Deleting the dashed links 1–6 and 2–10 gives the three cells $\{1\}$,
$\{2,4,6,7,9\}$ and $\{3,5,8,10\}$.
```
:::

(sec-08-05-shortest-paths)=
## Shortest paths and shortest-path models

:::{prf:example} Dijkstra's algorithm on an undirected network
:label: ex-08-05-undirected-dijkstra

An undirected network has nodes $0,1,\dots,5$ and the edges (length in
minutes)

| Edge | 0–2 | 0–4 | 0–5 | 1–3 | 1–5 | 2–3 | 2–5 | 3–4 | 3–5 | 4–5 |
|---|---|---|---|---|---|---|---|---|---|---|
| Length | 13 | 16 | 8 | 6 | 10 | 14 | 11 | 5 | 17 | 7 |

Find the shortest paths from node 0 to all other nodes.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Dijkstra's algorithm ({prf:ref}`alg-04-02-dijkstra`) works on undirected networks
as well: every edge can be used in both directions with the same length. The
labels are $[V(j),P(j)]$ (shortest known length, predecessor); in every
iteration the temporary label with the smallest $V$ becomes permanent.

| Iteration | Made permanent | Temporary labels after scanning it |
|---|---|---|
| 0 | 0 $[0,-]$ | 2 $[13,0]$, 4 $[16,0]$, 5 $[8,0]$ |
| 1 | 5 $[8,0]$ | 1 $[18,5]$, 2 $[13,0]$, 3 $[25,5]$, 4 $[15,5]$ |
| 2 | 2 $[13,0]$ | 1 $[18,5]$, 3 $[25,5]$, 4 $[15,5]$ |
| 3 | 4 $[15,5]$ | 1 $[18,5]$, 3 $[20,4]$ |
| 4 | 1 $[18,5]$ | 3 $[20,4]$ |
| 5 | 3 $[20,4]$ | – |

In iteration 1, node 4 improves from 16 to $8+7=15$, while node 2 keeps 13
($8+11=19$ is longer). In iteration 2 the path via 2 to node 3 would cost
$13+14=27>25$. In iteration 3 node 3 improves to $15+5=20$, and in
iteration 4 the path via 1 ($18+6=24$) does not beat it.

**Answer.** The shortest lengths from node 0 are 8 min to node 5, 13 to node
2, 15 to node 4, 18 to node 1 and 20 to node 3. Following the predecessors
gives the shortest-path tree $0$–$5$, $0$–$2$, $5$–$4$, $5$–$1$, $4$–$3$
({numref}`fig-08-05-undirected-network`); for example, the path to node 3 is
$0\to5\to4\to3$ with $8+7+5=20$.

```{figure} figures/08-05-undirected-network.svg
:name: fig-08-05-undirected-network
:width: 420px
:class: dark-light
:alt: Six-node undirected network with final labels and the shortest-path tree from node 0 highlighted

Final labels $[V(j),P(j)]$ and the shortest-path tree of
{prf:ref}`ex-08-05-undirected-dijkstra`.
```
:::

:::{prf:example} Dijkstra's algorithm on a directed network
:label: ex-08-05-dijkstra-seven

A directed network has nodes $1,\dots,7$ and the arcs (length in km)

| Arc | (1,2) | (1,3) | (2,4) | (2,5) | (2,6) | (3,2) | (3,4) | (3,5) |
|---|---|---|---|---|---|---|---|---|
| Length | 5 | 1 | 7 | 1 | 6 | 2 | 6 | 7 |

| Arc | (4,3) | (4,6) | (4,7) | (5,4) | (5,6) | (5,7) | (6,2) | (6,7) |
|---|---|---|---|---|---|---|---|---|
| Length | 7 | 4 | 6 | 3 | 5 | 9 | 7 | 2 |

Use Dijkstra's algorithm to find the shortest routes from node 1 to every
other node.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Labels $[V(j),P(j)]$; ties between temporary labels are broken by the
smaller node number, and a label is changed only by a strictly shorter path.

| Iteration | Made permanent | Temporary labels after scanning it |
|---|---|---|
| 0 | 1 $[0,-]$ | 2 $[5,1]$, 3 $[1,1]$ |
| 1 | 3 $[1,1]$ | 2 $[3,3]$, 4 $[7,3]$, 5 $[8,3]$ |
| 2 | 2 $[3,3]$ | 4 $[7,3]$, 5 $[4,2]$, 6 $[9,2]$ |
| 3 | 5 $[4,2]$ | 4 $[7,3]$, 6 $[9,2]$, 7 $[13,5]$ |
| 4 | 4 $[7,3]$ | 6 $[9,2]$, 7 $[13,5]$ |
| 5 | 6 $[9,2]$ | 7 $[11,6]$ |
| 6 | 7 $[11,6]$ | – |

Notes on the steps: in iteration 1 node 2 improves from 5 to $1+2=3$; in
iteration 2 node 5 improves from 8 to $3+1=4$; in iteration 3 the paths via 5
to node 4 ($4+3=7$) and to node 6 ($4+5=9$) only tie with the existing
labels; in iteration 5 node 7 improves from 13 to $9+2=11$.

| Node | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| Distance (km) | 3 | 1 | 7 | 4 | 9 | 11 |
| Route | 1-3-2 | 1-3 | 1-3-4 | 1-3-2-5 | 1-3-2-6 | 1-3-2-6-7 |

Because of the ties, nodes 4, 6 and 7 have a second shortest route of the same
length: 1-3-2-5-4 (7 km), 1-3-2-5-6 (9 km) and 1-3-2-5-6-7 (11 km).
{numref}`fig-08-05-dijkstra-seven` shows the shortest-path tree.
*(An earlier solution lists the routes but not their lengths; the lengths
above complete it.)*

```{figure} figures/08-05-dijkstra-seven.svg
:name: fig-08-05-dijkstra-seven
:width: 460px
:class: dark-light
:alt: Seven-node directed network with arc lengths and the shortest-path tree from node 1 highlighted

The network of {prf:ref}`ex-08-05-dijkstra-seven` with the shortest-path tree
from node 1 highlighted and final distances next to the nodes.
```
:::

:::{prf:example} A purchasing plan as a shortest path
:label: ex-08-05-purchasing

A retailer needs 100, 140, 210 and 180 units of an item in the next four
months. The purchase price per unit is \$15, \$12, \$10 and \$14 in months
1 to 4, and every order costs \$200 regardless of its size. Stock may be
bought in a month to cover the demand of later months; each unit held from
one month to the next costs \$1.20. Demand must be met on time and no stock
is left at the end. Model the problem as a shortest-path problem and solve
it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Network.** Let node $i=1,\dots,4$ mean "start of month $i$, no stock on
hand" and node 5 the end of the horizon. An arc $(i,j)$ with $i<j$ means:
order at the start of month $i$ exactly the demand of months $i,\dots,j-1$.
Its cost is the order cost, the purchase cost, and the holding cost of the
units kept for later months:

$$
c_{ij}=200+p_i\sum_{k=i}^{j-1}d_k+1.2\sum_{k=i}^{j-1}(k-i)\,d_k,
$$

with $d=(100,140,210,180)$ and $p=(15,12,10,14)$. For example
$c_{13}=200+15\cdot240+1.2\cdot140=\$3\,968$ (the 140 units of month 2 are
held one month). All ten arcs:

| Arc | (1,2) | (1,3) | (1,4) | (1,5) | (2,3) | (2,4) | (2,5) | (3,4) | (3,5) | (4,5) |
|---|---|---|---|---|---|---|---|---|---|---|
| Units ordered | 100 | 240 | 450 | 630 | 140 | 350 | 530 | 210 | 390 | 180 |
| Cost (\$) | 1 700 | 3 968 | 7 622 | 10 970 | 1 880 | 4 652 | 7 244 | 2 300 | 4 316 | 2 720 |

A path from 1 to 5 is a purchasing plan (its arcs say when to order and for
which months), and its length is the plan's total cost.

**Solution.** The network is acyclic, so the shortest distances follow in
node order from $V(1)=0$ and Bellman's equations {eq}`eq-04-02-bellman`,
$V(j)=\min_{i<j}\{V(i)+c_{ij}\}$:

- $V(2)=1\,700$;
- $V(3)=\min\{3\,968,\ 1\,700+1\,880\}=3\,580$ (via 2);
- $V(4)=\min\{7\,622,\ 1\,700+4\,652,\ 3\,580+2\,300\}=5\,880$ (via 3);
- $V(5)=\min\{10\,970,\ 1\,700+7\,244,\ 3\,580+4\,316,\ 5\,880+2\,720\}=7\,896$
  (via 3).

The shortest path is $1\to2\to3\to5$ ({numref}`fig-08-05-purchasing`): buy 100
units in month 1, 140 units in month 2, and $210+180=390$ units in month 3,
holding 180 of them for one month. The minimum total cost is **\$7 896**.
It pays to buy ahead in month 3 because the month-3 price (\$10) plus one
month of holding (\$1.20) is far below the month-4 price (\$14), while buying
ahead in month 1 for month 2 would cost $15+1.20=\$16.20$ against \$12.
*(An earlier solution identifies this plan but does not state its cost.)*

```{figure} figures/08-05-purchasing.svg
:name: fig-08-05-purchasing
:width: 560px
:class: dark-light
:alt: Five nodes in a row with ten arcs labeled with costs, the path 1 to 2 to 3 to 5 highlighted

The purchasing network of {prf:ref}`ex-08-05-purchasing` (costs in dollars). The
shortest path $1\to2\to3\to5$ (\$7 896) is highlighted.
```
:::

:::{prf:example} A knapsack problem as a longest path
:label: ex-08-05-knapsack-path

A hiker's backpack holds 5 ft³. Three items are available, with volumes 2,
3 and 4 ft³ and values 30, 50 and 70 (on a scale of 0 to 100). Each item can
be taken at most once. Model the choice of items as a longest-path problem
on a network whose nodes are $[i,w]$: item $i$ is about to be considered and
$w$ ft³ of space remain. Find the best choice.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Network.** Start at $[1,5]$. From $[i,w]$ there are two arcs: "skip item
$i$" to $[i+1,w]$ with value 0, and, if the item fits, "take item $i$" to
$[i+1,w-\text{vol}_i]$ with the item's value. After the last item all nodes
lead to an end node E. The nodes that can be reached are

- stage 1: $[1,5]$;
- stage 2: $[2,5]$ (skip item 1) and $[2,3]$ (take it, value 30);
- stage 3: $[3,5]$, $[3,2]$ (from $[2,5]$, taking item 2 gives 50) and
  $[3,3]$, $[3,0]$ (from $[2,3]$, taking item 2 gives 50);
- end: from $[3,5]$ item 3 fits, and of the two possible arcs to E (skip it,
  value 0, or take it, value 70) only the better one is drawn; from $[3,2]$,
  $[3,3]$, $[3,0]$ it does not fit (value 0).

**Longest path.** The network ({numref}`fig-08-05-knapsack`) is acyclic, and
every path from $[1,5]$ to E is a feasible packing whose length is its total
value. The four paths have values

| Path | Items | Volume (ft³) | Value |
|---|---|---|---|
| $[1,5]\to[2,5]\to[3,5]\to E$ | 3 | 4 | 70 |
| $[1,5]\to[2,5]\to[3,2]\to E$ | 2 | 3 | 50 |
| $[1,5]\to[2,3]\to[3,3]\to E$ | 1 | 2 | 30 |
| $[1,5]\to[2,3]\to[3,0]\to E$ | 1, 2 | 5 | **80** |

so the hiker takes items 1 and 2 (5 ft³, value 80). A longest path in an
acyclic network is found by the same recursion as a shortest path, with
$\max$ in place of $\min$, or as a shortest path with negated values; the
recursion is the dynamic program of {ref}`sec-05-03-knapsack` (compare
{ref}`sec-04-02-longest`).

*(An earlier network for this problem allows two units of item 1, with an arc
of value 60; that is the integer rather than the 0–1 knapsack. Its optimum is
also 80, because two units of item 1 leave only 1 ft³.)*

```{figure} figures/08-05-knapsack.svg
:name: fig-08-05-knapsack
:width: 480px
:class: dark-light
:alt: Layered network of knapsack states from node 1,5 to an end node with item values on the arcs and the longest path highlighted

The knapsack network of {prf:ref}`ex-08-05-knapsack-path`; nodes $[i,w]$, arc
values, and the longest path (value 80) highlighted.
```
:::

:::{prf:example} Chains of contacts
:label: ex-08-05-hide-and-seek

Six children, Joe, Kay, Jim, Bob, Rae and Kim, play a version of hide and
seek in which each child knows the hiding places of only a few others:

- Joe knows where Bob, Jim and Rae hide;
- Kay knows where Bob, Jim and Rae hide;
- Jim and Bob each know only where Kay hides;
- Rae knows where Kim hides;
- Kim knows where Joe and Bob hide.

A child looking for another may ask the children whose places it knows, who
then ask the children they know, and so on. Using a shortest-path model, find
for every child a plan to find every other child with the fewest contacts.
Which search needs the most contacts?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Network.** One node per child and an arc $(a,b)$ of length 1 when $a$
knows where $b$ hides: 11 arcs ({numref}`fig-08-05-hide-and-seek`). A chain of
contacts from $a$ to $b$ is a directed path, and the number of contacts is
its number of arcs. With unit lengths, Dijkstra's algorithm reduces to a
breadth-first search from each child.

**Shortest chains.** Rows are the searching child, columns the child sought,
entries the fewest contacts:

| Searcher | Joe | Kay | Jim | Bob | Rae | Kim |
|---|---|---|---|---|---|---|
| **Joe** | – | 2 | 1 | 1 | 1 | 2 |
| **Kay** | 3 | – | 1 | 1 | 1 | 2 |
| **Jim** | 4 | 1 | – | 2 | 2 | 3 |
| **Bob** | 4 | 1 | 2 | – | 2 | 3 |
| **Rae** | 2 | 3 | 3 | 2 | – | 1 |
| **Kim** | 1 | 2 | 2 | 1 | 2 | – |

The plans (shortest-path trees) are:

- Joe: Bob, Jim, Rae directly; Kay via Bob (Joe–Bob–Kay); Kim via Rae.
- Kay: Bob, Jim, Rae directly; Kim via Rae; Joe via Rae and Kim
  (Kay–Rae–Kim–Joe).
- Jim: Kay; then Bob and Rae via Kay; Kim via Kay and Rae; Joe via Kay, Rae
  and Kim (Jim–Kay–Rae–Kim–Joe).
- Bob: the same as Jim, with Jim reached via Kay (Bob–Kay–Jim).
- Rae: Kim; Joe and Bob via Kim; Jim via Kim and Joe; Kay via Kim and Bob.
- Kim: Joe and Bob directly; Jim and Rae via Joe; Kay via Bob.

**Largest number of contacts.** The longest searches need **4 contacts**:
Jim looking for Joe and Bob looking for Joe, both along
$\cdot\to$Kay$\to$Rae$\to$Kim$\to$Joe. Joe is hard to find because only Kim
knows his place, and Kim is reached only through Rae.

*(An earlier drawing of this network joins Joe and Kim in both directions and
omits the arc from Joe to Jim; by the list above Joe does not know where Kim
hides, and he does know where Jim hides. The maximum of 4 contacts is the
same with the correct network.)*

```{figure} figures/08-05-hide-and-seek.svg
:name: fig-08-05-hide-and-seek
:width: 420px
:class: dark-light
:alt: Directed network of six children with an arc from each child to every child whose hiding place it knows

The "who knows where" network of {prf:ref}`ex-08-05-hide-and-seek`. The chain
Bob → Kay → Rae → Kim → Joe (4 contacts) is highlighted.
```
:::

(sec-08-05-max-flow)=
## Maximum flow and matching

:::{prf:example} Recovering a feasible flow from the final residual network
:label: ex-08-05-recover-flow

A network has source 1, sink 4 and arcs
$(1,2)$, $(1,3)$, $(2,3)$, $(2,4)$, $(3,4)$, each of capacity 5. (a) Run the
Ford–Fulkerson algorithm with the labeling procedure from the zero flow.
(b) Show that adding up the path flows does not give a feasible flow.
(c) Recover a feasible maximum flow from the final residual network, and give
a minimum cut.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** *Iteration 1.* From 1 both 2 and 3 have residual 5; the tie goes to
node 2, label $[5,1]$. From 2 the candidates 3 and 4 tie at 5; label 3 with
$[5,2]$. From 3 label 4 with $[5,3]$. Path $1\to2\to3\to4$, $\Delta=5$. Now
$r_{12}=r_{23}=r_{34}=0$ and $r_{21}=r_{32}=r_{43}=5$.

*Iteration 2.* From 1 only node 3 ($r_{13}=5$): label $[5,1]$. From 3 the only
candidate is 2, through the **reverse** residual arc $(3,2)$ with $r_{32}=5$:
label $[5,3]$. From 2 label 4 with $[5,2]$. Path $1\to3\to2\to4$, $\Delta=5$.
Now also $r_{13}=r_{32}=r_{24}=0$ and $r_{31}=r_{23}=r_{42}=5$.

*Iteration 3.* $r_{12}=r_{13}=0$: no arc leaves the source. Stop with
$v^*=5+5=10$.

The final residual network ({numref}`fig-08-05-recover-flow`(c)) has exactly
the arcs $(2,1)$, $(3,1)$, $(2,3)$, $(4,2)$ and $(4,3)$, each with residual
capacity 5.

**(b)** Adding the path flows arc by arc gives $x_{12}=x_{23}=x_{34}=5$ from
the first path and $x_{13}=x_{32}=x_{24}=5$ from the second. But the network
has no arc $(3,2)$: the second path did not send flow from 3 to 2, it
*canceled* the 5 units that the first path had put on $(2,3)$.

**(c)** Use $x_{ij}-x_{ji}=u_{ij}-r_{ij}$ {eq}`eq-04-03-net-flow` and keep the
positive direction {eq}`eq-04-03-recover`:

| Pair | $u_{ij}-r_{ij}$ | $u_{ji}-r_{ji}$ | Flow |
|---|---|---|---|
| (1,2) | $5-0=5$ | $0-5=-5$ | $x_{12}=5$ |
| (1,3) | $5-0=5$ | $0-5=-5$ | $x_{13}=5$ |
| (2,3) | $5-5=0$ | $0-0=0$ | $x_{23}=0$ |
| (2,4) | $5-0=5$ | $0-5=-5$ | $x_{24}=5$ |
| (3,4) | $5-0=5$ | $0-5=-5$ | $x_{34}=5$ |

The flow sends 5 units along $1\to2\to4$ and 5 along $1\to3\to4$; it is
feasible, and its value is 10. The last labeling reached only node 1, so
$S=\{1\}$ with capacity $u_{12}+u_{13}=10$ is a minimum cut; the cuts
$S=\{1,3\}$ and $S=\{1,2,3\}$ also have capacity 10.

*(An earlier drawing of the final residual network shows the arc between 3
and 4 as $3\to4$. Since $x_{34}=u_{34}=5$, the residual capacities are
$r_{34}=0$ and $r_{43}=5$, so the arc is $4\to3$; with the arc drawn the other
way, the recovery rule would wrongly give $x_{34}=0$. The maximum flow value
is 10.)*

```{figure} figures/08-05-recover-flow.svg
:name: fig-08-05-recover-flow
:width: 600px
:class: dark-light
:alt: Three panels, the original network, the residual network after the first augmentation and the final residual network

{prf:ref}`ex-08-05-recover-flow`: (a) the network with capacities; (b) the
residual network after sending 5 units along $1\to2\to3\to4$, with the second
augmenting path highlighted; (c) the final residual network.
```
:::

:::{prf:example} Maximum flow in an eight-node network
:label: ex-08-05-maxflow-eight

Find a maximum flow from $s$ to $t$ and a minimum cut in the network with
nodes $s,a,b,c,d,e,f,t$ and arcs (capacity in vehicles per minute)

| Arc | $(s,a)$ | $(s,b)$ | $(c,s)$ | $(a,d)$ | $(a,e)$ | $(b,d)$ |
|---|---|---|---|---|---|---|
| $u_{ij}$ | 7 | 2 | 3 | 3 | 3 | 1 |

| Arc | $(b,f)$ | $(c,f)$ | $(d,t)$ | $(e,c)$ | $(e,t)$ | $(f,t)$ |
|---|---|---|---|---|---|---|
| $u_{ij}$ | 3 | 4 | 5 | 3 | 3 | 4 |

Use the Ford–Fulkerson algorithm ({prf:ref}`alg-04-03-ford-fulkerson`) with the
labeling procedure ({prf:ref}`alg-04-03-labelling`; ties broken alphabetically,
$t$ last).
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Iteration 1.** From $s$: $a$ (7) beats $b$ (2), label $a$ $[7,s]$; from
$a$: $d$ and $e$ tie at 3, label $d$ $[3,a]$; from $d$: label $t$ $[5,d]$.
Path $s\to a\to d\to t$, $\Delta=3$.

**Iteration 2.** $a$ $[4,s]$; from $a$ only $e$ remains ($r_{ad}=0$), $[3,a]$;
from $e$: $c$ and $t$ tie at 3, label $c$ $[3,e]$; from $c$: $f$ $[4,c]$ ($s$
is labeled); from $f$: $t$ $[4,f]$. Path $s\to a\to e\to c\to f\to t$,
$\Delta=3$. Total 6.

**Iteration 3.** From $s$: $b$ (2) beats $a$ (1), label $b$ $[2,s]$; from
$b$: $f$ (3) beats $d$ (1), $[3,b]$; from $f$: the reverse arc $(f,c)$ has
residual 3 (the 3 units on $(c,f)$) and $(f,t)$ has 1, label $c$ $[3,f]$; from
$c$: the reverse arc $(c,e)$ has 3, label $e$ $[3,c]$; from $e$: the reverse
arc $(e,a)$ and $(e,t)$ tie at 3, label $a$ $[3,e]$; from $a$: nothing new
($r_{ad}=0$, $s$ and $e$ labeled), so **backtrack**: delete $a$, return to
$e$, label $t$ $[3,e]$. Path $s\to b\to f\to c\to e\to t$,
$\Delta=\min\{2,3,3,3,3\}=2$. Total 8. The path uses two reverse arcs: of
the 3 units that $a$ had routed $e\to c\to f\to t$, 2 now go from $e$
straight to $t$, and the capacity they free on $f\to t$ takes 2 new units
from $b$.

**Iteration 4.** From $s$ only $a$ ($r_{sa}=1$), $[1,s]$; from $a$ nothing
($r_{ad}=r_{ae}=0$): delete $a$; $s$ has no other candidate. Stop.

**Result.** $v^*=3+3+2=8$ vehicles per minute. The recovered flows are

| Arc | $(s,a)$ | $(s,b)$ | $(a,d)$ | $(a,e)$ | $(b,f)$ | $(d,t)$ |
|---|---|---|---|---|---|---|
| $x_{ij}$ | 6 | 2 | 3 | 3 | 2 | 3 |

| Arc | $(e,c)$ | $(e,t)$ | $(c,f)$ | $(f,t)$ | $(b,d)$ | $(c,s)$ |
|---|---|---|---|---|---|---|
| $x_{ij}$ | 1 | 2 | 1 | 3 | 0 | 0 |

Check: node $e$ receives 3 and sends $1+2$; node $f$ receives $2+1$ and sends
3; the sink receives $3+2+3=8$.

**Minimum cut.** The last search reached $S=\{s,a\}$. Its forward arcs are
$(s,b)$, $(a,d)$, $(a,e)$ with capacity $2+3+3=8=v^*$, and its backward arc
$(c,s)$ carries no flow ({numref}`fig-08-05-maxflow-eight`). It is the only
minimum cut. The bottleneck is node $a$: it receives up to 7 units but can
pass on only 6.

```{figure} figures/08-05-maxflow-eight.svg
:name: fig-08-05-maxflow-eight
:width: 520px
:class: dark-light
:alt: Eight-node flow network with flow over capacity labels and the source side of the minimum cut highlighted

Maximum flow of {prf:ref}`ex-08-05-maxflow-eight` with labels $x_{ij}/u_{ij}$.
The highlighted nodes $s$ and $a$ form the source side of the minimum cut.
```
:::

:::{prf:example} Assigning chores by maximum flow
:label: ex-08-05-chores

A parent wants five household chores (numbered 1 to 5) done by five children,
each child doing at most one chore and only a chore it accepts:

| Child | Rif | Mai | Ben | Kim | Ken |
|---|---|---|---|---|---|
| Accepted chores | 3, 4, 5 | 1 | 1, 2 | 1, 2, 5 | 2 |

Find the largest number of chores that can be done, and an assignment, with
the Ford–Fulkerson algorithm.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Network.** This is a bipartite matching problem ({ref}`sec-04-03-matching`):
a source $s$ with an arc of capacity 1 to every child, an arc of capacity 1
from each child to each chore it accepts, and an arc of capacity 1 from every
chore to the sink $t$ ({numref}`fig-08-05-chores`). An integral flow of value $v$
is an assignment of $v$ chores.

**Ford–Fulkerson.** All residual capacities are 0 or 1, so the labeling rule
breaks every tie by the order of the table (children), then chore numbers.

1. $s\to$Rif$\to3\to t$.
2. $s\to$Mai$\to1\to t$.
3. From Ben, chore 1 is tried first; it leads only back to Mai (reverse arc),
   a dead end, so Mai and then chore 1 are deleted. Chore 2 gives
   $s\to$Ben$\to2\to t$.
4. From Kim, chores 1 and 2 lead back to Mai and Ben, both dead ends; chore 5
   gives $s\to$Kim$\to5\to t$.
5. Only Ken is left at $s$; Ken → 2 → Ben → 1 → Mai is a dead end
   everywhere. No augmenting path: stop.

The maximum is **4 chores**: Rif does 3, Mai 1, Ben 2, Kim 5, and Ken has no
chore.

**Why not 5?** The last search reached $S=\{s,\text{Mai},\text{Ben},\text{Ken},1,2\}$.
The cut has forward arcs $s\to$Rif, $s\to$Kim, $1\to t$, $2\to t$, capacity
4. In words: Mai, Ben and Ken accept only chores 1 and 2, so at most two of
those three children can work, and at most $5-1=4$ chores get done.

**All maximum assignments.** Kim must take chore 5 (otherwise chores 3, 4, 5
can only go to Rif, and at most 3 chores are done); Rif takes 3 or 4; and
chores 1 and 2 go to Mai and Ben, Mai and Ken (Mai 1, Ken 2), or Ben and Ken
(Ben 1, Ken 2). That gives $2\times3=6$ maximum assignments.

*(An earlier answer lists the assignment Rif 3, Mai 1, Ben 2, Ken 5; it is
infeasible, because Ken accepts only chore 2.)*

```{figure} figures/08-05-chores.svg
:name: fig-08-05-chores
:width: 480px
:class: dark-light
:alt: Bipartite network from source through five children and five chores to the sink with the four assigned pairs highlighted

The matching network of {prf:ref}`ex-08-05-chores` (all capacities 1) with the
maximum assignment highlighted.
```
:::

(sec-08-05-min-cost-flow)=
## Minimum-cost flow, transportation and assignment

:::{prf:example} Transshipment through a warehouse with a dummy market
:label: ex-08-05-transshipment

Two factories F1 and F2 produce 60 and 50 units per week. Three stores S1,
S2, S3 need 30, 40 and 25 units per week. Goods may go through a warehouse W.
The shipping routes and their costs (\$ per unit) are

| Route | F1→W | F2→W | F1→S1 | W→S1 | W→S2 | W→S3 | F2→S3 |
|---|---|---|---|---|---|---|---|
| Cost (\$) | 2 | 3 | 7 | 3 | 4 | 6 | 5 |

The route F1→W can carry at most 40 units per week; the others are
unlimited.

(a) Formulate the problem as a balanced minimum-cost flow problem.
(b) Find the optimal shipping plan and prove its optimality with node
potentials.
(c) What is one more unit of capacity on F1→W worth?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** The factories supply $60+50=110$ units, the stores need
$30+40+25=95$. Add a dummy store D with $b_D=-15$ and arcs F1→D, F2→D of cost
0; the flow into D is the unused capacity. With $b=(60,50,0,-30,-40,-25,-15)$
for F1, F2, W, S1, S2, S3, D, the model is {eq}`eq-04-04-mcf` with nine arcs,
the bound $0\le x_{F1,W}\le40$ and all other $x\ge0$
({numref}`fig-08-05-transshipment`).

**(b)** The optimal plan (solving the linear program with HiGHS) is

| Arc | F1→W | F2→W | F1→S1 | W→S1 | W→S2 | W→S3 | F2→S3 | F1→D | F2→D |
|---|---|---|---|---|---|---|---|---|---|
| $x_{ij}$ | 40 | 25 | 5 | 25 | 40 | 0 | 25 | 15 | 0 |

with cost $2\cdot40+3\cdot25+7\cdot5+3\cdot25+4\cdot40+5\cdot25=\$550$ per
week. Factory F1 keeps 15 units of unused capacity.

*Optimality.* Take the potentials $\pi_{F1}=4$, $\pi_{F2}=3$, $\pi_W=0$,
$\pi_{S1}=-3$, $\pi_{S2}=-4$, $\pi_{S3}=-2$, $\pi_D=4$. The reduced costs
$\bar c_{ij}=c_{ij}-\pi_i+\pi_j$ are

| Arc | F1→W | F2→W | F1→S1 | W→S1 | W→S2 | W→S3 | F2→S3 | F1→D | F2→D |
|---|---|---|---|---|---|---|---|---|---|
| $\bar c_{ij}$ | $-2$ | 0 | 0 | 0 | 0 | 4 | 0 | 0 | 1 |

The six arcs with $\bar c_{ij}=0$ carry flow strictly inside their bounds;
F1→W has $\bar c<0$ and is at its upper bound 40; W→S3 and F2→D have
$\bar c>0$ and carry nothing. The conditions {eq}`eq-04-04-cs` hold, so the plan
is optimal, and since no arc at a bound has $\bar c=0$ it is the only optimal
plan. Check the balances: W receives $40+25$ and sends $25+40$; S1 receives
$5+25=30$.

**(c)** The reduced cost $\bar c_{F1,W}=-2$ says that each extra unit of
capacity on F1→W saves \$2 per week. The saving is visible in the network:
F1 sends one more unit to W (\$2), W forwards it to S1 (\$3), and F1 sends one
unit less directly to S1 (\$7), a change of $2+3-7=-\$2$. This rate holds
only for small increases (while F1 still ships directly to S1); without any
limit on F1→W the optimal cost is \$525 per week.

```{figure} figures/08-05-transshipment.svg
:name: fig-08-05-transshipment
:width: 460px
:class: dark-light
:alt: Transshipment network from two factories through a warehouse to three stores and a dummy store, with the arcs used in the optimal plan highlighted

The network of {prf:ref}`ex-08-05-transshipment`: unit costs on the arcs, the
bound $(0,40)$ on F1→W, net supplies in brackets, and the dummy store D.
Highlighted arcs carry flow in the optimal plan (\$550 per week).
```
:::

:::{prf:example} Assignment with a forbidden pair
:label: ex-08-05-assignment-forbidden

Four operators O1–O4 are to be assigned to four machines M1–M4, one operator
per machine. The setup times in minutes are below; operator O3 is not
qualified for machine M2.

| | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| O1 | 10 | 7 | 8 | 12 |
| O2 | 6 | 9 | 11 | 7 |
| O3 | 9 | – | 5 | 8 |
| O4 | 7 | 6 | 9 | 10 |

Find the assignment with the least total setup time and prove that it is
optimal.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Model.** A minimum-cost flow problem with $b=1$ for the operators, $b=-1$
for the machines, and one arc per allowed pair (the pair O3–M2 simply has no
arc). By total unimodularity the LP has an optimal 0–1 solution
({prf:ref}`thm-04-04-integrality`).

**Reduction.** Subtracting a constant from a row or a column of the cost
matrix changes every assignment's total by the same constant, because every
assignment uses each row and each column once. Subtract the row minima
$7,6,5,6$ and then the column minimum 1 of M4:

| | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| O1 | 3 | **0** | 1 | 4 |
| O2 | **0** | 3 | 5 | **0** |
| O3 | 4 | – | **0** | 2 |
| O4 | 1 | **0** | 3 | 3 |

Every assignment now costs at least 0, and the reductions add up to
$7+6+5+6+1=25$, so every assignment costs at least 25 minutes. But no
assignment uses only zeros: O1 and O4 both need M2. The zeros are covered by
three lines (column M2, row O2, column M3); the smallest uncovered entry is 1
(O4–M1). Subtract 1 from the three uncovered rows O1, O3, O4 and add 1 to the
two covered columns M2, M3 (net: $-1$ on uncovered entries, $+1$ where a
covered row meets a covered column):

| | M1 | M2 | M3 | M4 |
|---|---|---|---|---|
| O1 | 2 | **0** | 1 | 3 |
| O2 | **0** | 4 | 6 | **0** |
| O3 | 3 | – | **0** | 1 |
| O4 | **0** | **0** | 3 | 2 |

The lower bound rises to $25+3-2=26$. Now the zeros O1–M2, O2–M4, O3–M3,
O4–M1 form an assignment.

**Answer.** O1 → M2, O2 → M4, O3 → M3, O4 → M1, with total setup time
$7+7+5+7=26$ minutes, equal to the lower bound, hence optimal. In the
language of {ref}`sec-04-04-optimality` the reductions are node potentials:
$\pi_{O}=(8,6,6,7)$ for the operators, $\pi_{M}=(0,1,1,-1)$ for the machines,
and the last table lists the reduced costs $c_{ij}-\pi_i+\pi_j\ge0$, zero on
the chosen arcs; the dual objective is $\sum\pi_O-\sum\pi_M=27-1=26$.

**Uniqueness.** Every assignment costs 26 plus the sum of its entries in the
last table, so an optimal assignment uses only zeros. O1 has a single zero
(M2); then O4 must take M1, O3 must take M3 and O2 must take M4. The optimum is
unique; of the 18 assignments that respect the forbidden pair, the next best
costs 28 minutes.
:::

:::{prf:example} Cycle canceling on a small network
:label: ex-08-05-cycle-cancelling

Ten units must be sent from node 1 to node 4 ($b_1=10$, $b_4=-10$) in the
network with arcs (unit cost in \$, capacity in units)

| Arc | (1,2) | (1,3) | (2,3) | (2,4) | (3,4) |
|---|---|---|---|---|---|
| $c_{ij}$ (\$) | 2 | 5 | 1 | 6 | 2 |
| $u_{ij}$ | 8 | 10 | 5 | 6 | 10 |

Start from the feasible flow that sends 6 units along $1\to2\to4$ and 4 units
along $1\to3\to4$, and apply the cycle-canceling algorithm
({prf:ref}`alg-04-04-cycle-cancelling`), canceling the most negative cycle each
time.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Start.** $x_{12}=6$, $x_{24}=6$, $x_{13}=4$, $x_{34}=4$, $x_{23}=0$; cost
$2\cdot6+6\cdot6+5\cdot4+2\cdot4=\$76$.

**Residual network.** Forward arcs (cost $c_{ij}$) where $x_{ij}<u_{ij}$,
backward arcs (cost $-c_{ij}$) where $x_{ij}>0$:

| Residual arc | $1\to2$ | $2\to1$ | $1\to3$ | $3\to1$ | $2\to3$ | $4\to2$ | $3\to4$ | $4\to3$ |
|---|---|---|---|---|---|---|---|---|
| Cost (\$) | 2 | $-2$ | 5 | $-5$ | 1 | $-6$ | 2 | $-2$ |
| Residual capacity | 2 | 6 | 6 | 4 | 5 | 6 | 6 | 4 |

($2\to4$ is saturated, so only its backward arc appears.)

**Cycle 1.** The directed cycles with negative cost are
$2\to3\to4\to2$: $1+2-6=-3$; $1\to2\to3\to1$: $2+1-5=-2$; and
$1\to3\to4\to2\to1$: $5+2-6-2=-1$. Cancel the most negative,
$2\to3\to4\to2$, with $\Delta=\min\{5,6,6\}=5$: $x_{23}=5$, $x_{34}=9$,
$x_{24}=1$. Cost $76-3\cdot5=\$61$.

**Cycle 2.** Now $2\to3$ is saturated. The only negative cycle left is
$1\to3\to4\to2\to1$ with cost $-1$ and residual capacities $1\to3$: 6,
$3\to4$: 1, $4\to2$: 1, $2\to1$: 6, so $\Delta=1$: $x_{13}=5$, $x_{34}=10$,
$x_{24}=0$, $x_{12}=5$. Cost $61-1=\$60$.

**Stop.** The residual cycles are now $1\to2\to4\to3\to1$ ($2+6-2-5=1$),
$1\to3\to2\to1$ ($5-1-2=2$) and $2\to4\to3\to2$ ($6-2-1=3$), plus the
two-arc cycles of cost 0; none is negative, so the flow

$$
(x_{12},x_{13},x_{23},x_{24},x_{34})=(5,5,5,0,10)
$$

is optimal with cost **\$60** ({prf:ref}`thm-04-04-negative-cycle`). Potentials
confirm it: with $\pi=(5,3,0,-3)$ the reduced costs are $0$ on $(1,2)$ and
$(1,3)$ (flows strictly inside their bounds), $-2$ on $(2,3)$ and $-1$ on
$(3,4)$ (both at their upper bounds), and $0$ on $(2,4)$ (at its lower bound),
so {eq}`eq-04-04-cs` holds. Five units travel $1\to2\to3\to4$ (\$5 each) and
five travel $1\to3\to4$ (\$7 each): $25+35=\$60$.
:::

## References

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows:
  Theory, algorithms, and applications*. Prentice Hall. (Chapters 4–7, 9, 12,
  13.)
- Bellman, R. (1958). On a routing problem. *Quarterly of Applied
  Mathematics, 16*(1), 87–90. <https://doi.org/10.1090/qam/102435>
- Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
  *Numerische Mathematik, 1*(1), 269–271.
  <https://doi.org/10.1007/BF01386390>
- Ford, L. R., & Fulkerson, D. R. (1956). Maximal flow through a network.
  *Canadian Journal of Mathematics, 8*, 399–404.
  <https://doi.org/10.4153/CJM-1956-045-5>
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapters 9–10.)
- Kleinberg, J., & Tardos, É. (2006). *Algorithm design*. Pearson. (Section
  4.7: clustering with maximum spacing.)
- Kruskal, J. B. (1956). On the shortest spanning subtree of a graph and the
  traveling salesman problem. *Proceedings of the American Mathematical
  Society, 7*(1), 48–50. <https://doi.org/10.1090/S0002-9939-1956-0078686-7>
- Prim, R. C. (1957). Shortest connection networks and some generalizations.
  *Bell System Technical Journal, 36*(6), 1389–1401.
  <https://doi.org/10.1002/j.1538-7305.1957.tb01515.x>
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 6: network models, problem sets of Sections 6.1–6.4.)
