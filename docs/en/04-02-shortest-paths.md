# Shortest-Path Problems

## Overview

Finding the cheapest, fastest or most reliable route through a network is
the most frequently solved network optimization problem: route planners,
packet routing in communication networks and traffic-assignment models solve
millions of shortest-path problems a day, and many planning problems that do
not look like routing (equipment replacement, purchasing plans, project
schedules) turn into shortest paths once the right network is drawn. This
chapter states the problem and its variants, and derives everything from one
observation: a part of a shortest path is itself a shortest path. That
observation gives **Bellman's equations** for the shortest-path lengths.
**Dijkstra's algorithm** solves these equations for one origin when all arc
lengths are nonnegative, by making node labels permanent in order of
increasing distance; we prove that it is correct, count its work, recover the
paths from the predecessor labels, and show by a counterexample why a negative
arc breaks it. The **Floyd–Warshall algorithm** finds the shortest paths
between all pairs of nodes, also with negative arcs, and detects negative
cycles. Finally the shortest-path problem is written as a linear program,
whose dual variables are the node labels, and used as a modeling tool. The
chapter builds on {doc}`04-01-networks-mst` (networks, paths, cycles, trees and
the node–arc incidence matrix); {ref}`sec-04-02-lp` also uses the dual of a
linear program and complementary slackness from {doc}`03-05-lp-duality`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- state the single-pair, single-source and all-pairs shortest-path problems
  and the assumptions (nonnegative lengths, no negative cycles) under which
  each method applies;
- derive Bellman's equations from the optimality of subpaths;
- run Dijkstra's algorithm by hand with labels $[V(j),P(j)]$, recover the
  shortest paths and the shortest-path tree, and prove the algorithm correct
  for nonnegative lengths;
- run the Floyd–Warshall algorithm, recover paths from the predecessor
  matrix and detect negative cycles;
- formulate the shortest-path problem as a linear program, interpret its dual
  as node potentials and verify optimality with reduced lengths; and
- model decision problems such as equipment replacement as shortest-path
  problems.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $s$, $t$ | origin and destination node of a shortest-path problem |
| $c_{ij}$ | length (distance, time or cost) of arc $(i,j)$ |
| $c(Q)$ | length of a path $Q$, the sum of the lengths of its arcs |
| $V(j)$ | label of node $j$: length of the best path from $s$ to $j$ found so far; at termination the shortest-path length |
| $P(j)$ | predecessor of node $j$ on that path |
| $[V(j),P(j)]$ | the label of node $j$ in Dijkstra's algorithm |
| $S$, $T$ | sets of nodes with permanent and with temporary labels |
| $\mD_k=(d_{ij}^{(k)})$ | Floyd–Warshall distance matrix after iteration $k$ |
| $\mP_k=(p_{ij}^{(k)})$ | Floyd–Warshall predecessor matrix after iteration $k$ |
| $\bar c_{ij}=c_{ij}+V(i)-V(j)$ | reduced length of arc $(i,j)$ with respect to the labels $V$ |
| $\mathcal W$, $c(\mathcal W)$ | a directed cycle and its length |

(sec-04-02-problem)=
## The shortest-path problem

Let $G=(\mathcal N,\mathcal A)$ be a directed network in which every arc
$(i,j)$ has a length $c_{ij}$. The **length** of a directed path
$Q=(i_1,i_2,\dots,i_p)$ is

$$
c(Q)=\sum_{k=1}^{p-1}c_{i_ki_{k+1}}.
$$ (eq-04-02-path-length)

:::{prf:definition} Shortest-path problems
:label: def-04-02-shortest-path

Given an origin $s$ and a destination $t$, the **shortest-path problem** asks
for a directed path from $s$ to $t$ of least length. Its common variants are

- **single pair**: one origin $s$ and one destination $t$;
- **single source**: shortest paths from $s$ to every other node;
- **single destination**: shortest paths from every node to $t$ (reverse all
  arcs and solve a single-source problem from $t$);
- **all pairs**: shortest paths between every ordered pair of nodes.

The length of a shortest path from $s$ to $j$ is the **distance** from $s$ to
$j$; it is $+\infty$ when no directed path exists.
:::

The single-pair problem is usually solved as a single-source problem stopped
early, because the methods find the distances to the nodes near $s$ on the
way to $t$ anyway. Two remarks fix the scope of the chapter.

- **Undirected networks.** An edge that can be used in both directions is
  replaced by the two arcs $(i,j)$ and $(j,i)$ with the same length
  ({ref}`sec-04-01-graphs`). An edge of negative length then forms a cycle
  $i\to j\to i$ of negative length, so in undirected networks the methods
  below assume nonnegative lengths.
- **Negative cycles.** If a directed cycle of negative total length can be
  reached from $s$ and can reach $t$, going round it again and again makes
  routes from $s$ to $t$ arbitrarily short; the problem over routes is then
  unbounded, and the problem over paths (no repeated nodes) becomes very hard
  to solve. All methods of this chapter therefore assume that **no negative
  cycle exists**. Dijkstra's algorithm needs the stronger assumption
  $c_{ij}\ge0$ for all arcs; the Floyd–Warshall algorithm needs only the
  absence of negative cycles, and detects them when they exist.

:::{prf:example} Shortest route from node 1 to node 7
:label: ex-04-02-route

{numref}`fig-04-02-dijkstra-network` shows a directed network with seven nodes
and twelve arcs; the arc labels are lengths, for example travel times in
minutes. Find a shortest path from node 1 to node 7.

```{figure} figures/04-02-dijkstra-network.svg
:name: fig-04-02-dijkstra-network
:width: 540px
:class: dark-light
:alt: Seven nodes joined by twelve directed arcs with lengths between 1 and 11

The network of {prf:ref}`ex-04-02-route`. Arc labels are lengths $c_{ij}$ in
minutes; arcs $(3,1)$, $(3,2)$ and $(5,3)$ point back towards the origin.
```
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The most direct approach is **enumeration**: list every directed path from
node 1 to node 7 and keep the shortest. This network has seven such paths
({numref}`tab-04-02-paths`), and the shortest one is $1\to4\to3\to7$ with length
$4+3+5=12$ minutes.

```{list-table} All directed paths from node 1 to node 7, by length (minutes)
:name: tab-04-02-paths
:header-rows: 1

* - Path
  - Length
* - $1\to4\to3\to7$
  - $4+3+5=12$
* - $1\to2\to5\to7$
  - $2+5+6=13$
* - $1\to4\to6\to7$
  - $4+1+9=14$
* - $1\to2\to5\to3\to7$
  - $2+5+8+5=20$
* - $1\to4\to3\to6\to7$
  - $4+3+7+9=23$
* - $1\to4\to3\to2\to5\to7$
  - $4+3+11+5+6=29$
* - $1\to2\to5\to3\to6\to7$
  - $2+5+8+7+9=31$
```
:::

Enumeration, as in the solution of {prf:ref}`ex-04-02-route`, does not
scale. In a network with $q$ layers of $r$ nodes, each
node joined to every node of the next layer, there are $r^q$ paths from the
first to the last layer, already $10^{20}$ for 20 layers of 10 nodes. The
methods below exploit the structure of the problem instead.

(sec-04-02-bellman)=
## Optimality of subpaths and Bellman's equations

The enumeration in {numref}`tab-04-02-paths` wastes work: the paths
$1\to4\to3\to7$ and $1\to4\to3\to6\to7$ share the subpath $1\to4\to3$, whose
length is computed twice, and once it is known that $1\to4\to3$ is the best way
to reach node 3, no path that reaches node 3 in another way needs to be
considered further. The following lemma makes this precise.

:::{prf:lemma} Optimality of subpaths
:label: lem-04-02-subpath

Assume that $G$ has no negative cycle. If $Q=(s,\dots,i,\dots,j)$ is a
shortest path from $s$ to $j$, then its part from $s$ to $i$ is a shortest
path from $s$ to $i$. In particular, if $(i,j)$ is the last arc of $Q$, the
distance to $j$ equals the distance to $i$ plus $c_{ij}$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Write $Q$ as $Q_1$ (from $s$ to $i$) followed by $Q_2$ (from $i$ to $j$).
Suppose some path $R$ from $s$ to $i$ is shorter than $Q_1$. Then $R$ followed
by $Q_2$ is a route from $s$ to $j$ of length $c(R)+c(Q_2)<c(Q)$. If it
repeats a node, it consists of a path from $s$ to $j$ plus one or more cycles;
the cycles have nonnegative length, so removing them gives a path from $s$ to
$j$ no longer than the route, hence shorter than $Q$. This contradicts the
optimality of $Q$.
:::

This is the **principle of optimality** of dynamic programming in its
simplest form: the remaining part of an optimal decision sequence must be
optimal for the state it starts from ({prf:ref}`thm-05-01-principle`). Applied to
the last arc of every shortest path, the lemma gives a system of equations for
the distances.

:::{prf:theorem} Bellman's equations
:label: thm-04-02-bellman

Assume that $G$ has no negative cycle and that every node can be reached from
$s$. Then the distances $V(j)$ from $s$ satisfy

$$
\begin{aligned}
V(s)&=0,\\
V(j)&=\min_{i:\,(i,j)\in\mathcal A}\bigl\{V(i)+c_{ij}\bigr\}
\qquad\text{for all } j\ne s.
\end{aligned}
$$ (eq-04-02-bellman)

If, in addition, every cycle has positive length, the distances are the only
solution of {eq}`eq-04-02-bellman`.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

*The distances satisfy the equations.* $V(s)=0$ because the path consisting
of $s$ alone has length 0 and no route back to $s$ is negative. For $j\ne s$,
every path from $s$ to $j$ ends with some arc $(i,j)$ and is at least as long
as $V(i)+c_{ij}$, so $V(j)\ge\min_i\{V(i)+c_{ij}\}$. Conversely, for every arc
$(i,j)$ a shortest path to $i$ followed by $(i,j)$ is a route to $j$ of length
$V(i)+c_{ij}$, and removing cycles from it does not make it longer, so
$V(j)\le V(i)+c_{ij}$.

*Uniqueness.* Let $U$ be any solution. For every $j\ne s$ pick an arc
$(q_j,j)$ that attains the minimum, so $U(j)=U(q_j)+c_{q_jj}$. Following
these arcs backwards from $j$ never returns to a node already visited: along a
cycle $j_1\to\dots\to j_p\to j_1$ of chosen arcs the equations would add up to
a cycle length of 0, which is excluded. So the backward walk ends at $s$ and
traces a path from $s$ to $j$ of length $U(j)$; hence $U(j)\ge V(j)$. On the
other hand, along a shortest path $s=i_1,\dots,i_p=j$ the inequalities
$U(i_{k+1})\le U(i_k)+c_{i_ki_{k+1}}$ add up to $U(j)\le V(j)$. Hence $U=V$.
:::

For the network of {prf:ref}`ex-04-02-route` the equation for node 7 reads

$$
\begin{aligned}
V(7)&=\min\{V(3)+5,\ V(5)+6,\ V(6)+9\}\\
&=\min\{7+5,\ 7+6,\ 5+9\}=12,
\end{aligned}
$$

using the distances $V(3)=7$, $V(5)=7$ and $V(6)=5$ found in
{prf:ref}`ex-04-02-dijkstra`. The
difficulty is visible here: the equations are circular in general (the
equation for node 3 contains $V(5)$, and the equation for node 5 contains
$V(2)$, while node 3 has an arc back to node 2). They can be solved one node at
a time only in a suitable order. In an acyclic network a topological order
does the job ({ref}`sec-05-04-acyclic`); with nonnegative lengths,
Dijkstra's algorithm finds the order while it runs, namely the order of
increasing distance.

(sec-04-02-dijkstra)=
## Dijkstra's algorithm

The method published by Dijkstra (1959) is a **label-setting** algorithm.
Every node $j$ carries a label $[V(j),P(j)]$:
$V(j)$ is the length of the best path from $s$ to $j$ found so far and $P(j)$
the node before $j$ on that path. A label is **temporary** while it may still
improve and becomes **permanent** when it is known to be optimal. In each
iteration the most recent permanent node $i$ is **scanned**: every arc
$(i,j)$ to a node without a permanent label is used to try to improve the
label of $j$. Then the temporary label with the smallest $V$ is made
permanent.

:::{prf:algorithm} Dijkstra's algorithm
:label: alg-04-02-dijkstra

**Input:** a directed network $G=(\mathcal N,\mathcal A)$ with lengths
$c_{ij}\ge0$, an origin $s$ and, optionally, a destination $t$.
**Output:** labels $[V(j),P(j)]$: the distance from $s$ to $j$ and the
predecessor of $j$ on a shortest path, for every node $j$ reachable from $s$
(for the nodes up to $t$ if a destination is given).

1. **Initialize.** Give $s$ the permanent label $[V(s),P(s)]=[0,-]$. Set
   $S=\{s\}$, $T=\varnothing$ and $i=s$.
2. **Stopping test.** If $i=t$, or if every node has a permanent label,
   stop.
3. **Scan node $i$.** For every arc $(i,j)\in\mathcal A$ with $j\notin S$:
   - if $j$ has no label yet, give it the temporary label
     $[V(i)+c_{ij},\,i]$ and set $T\leftarrow T\cup\{j\}$;
   - if $j\in T$ and $V(i)+c_{ij}<V(j)$, set $V(j)\leftarrow V(i)+c_{ij}$ and
     $P(j)\leftarrow i$;
   - otherwise leave the label of $j$ unchanged.
4. **Make a label permanent.** If $T=\varnothing$, stop: the remaining nodes
   cannot be reached from $s$. Otherwise choose $r\in T$ with the smallest
   $V(r)$ (ties broken arbitrarily), set $S\leftarrow S\cup\{r\}$,
   $T\leftarrow T\setminus\{r\}$, $i\leftarrow r$, and go to Step 2.
:::

The path to any permanent node $j$ is read **backwards** from the predecessor
labels, $j\leftarrow P(j)\leftarrow P(P(j))\leftarrow\cdots\leftarrow s$.
The algorithm works for undirected networks too, with every edge scanned in
both directions.

:::{prf:example} Dijkstra's algorithm on the seven-node network
:label: ex-04-02-dijkstra

Apply {prf:ref}`alg-04-02-dijkstra` to the network of {prf:ref}`ex-04-02-route`,
from $s=1$ to $t=7$. When two temporary labels tie, make the node with the
smaller number permanent.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Iteration 0.** Node 1 receives the permanent label $[0,-]$; $S=\{1\}$,
$T=\varnothing$, $i=1$.

**Iteration 1.** Node 1 is not the destination. Scanning its arcs gives node 2
the label $[0+2,1]=[2,1]$ and node 4 the label $[0+4,1]=[4,1]$, so
$T=\{2,4\}$. The smallest temporary label is $V(2)=2$, and node 2 becomes
permanent.

**Iteration 2.** Scanning node 2: arc $(2,5)$ gives node 5 the label
$[2+5,2]=[7,2]$. Now $T=\{4,5\}$ with $V(4)=4<V(5)=7$, and node 4 becomes
permanent.

**Iteration 3.** Scanning node 4: node 3 receives $[4+3,4]=[7,4]$ and node 6
receives $[4+1,4]=[5,4]$. Among $T=\{3,5,6\}$ the smallest label is
$V(6)=5$, and node 6 becomes permanent.

**Iteration 4.** Scanning node 6: node 7 receives $[5+9,6]=[14,6]$. Now
$T=\{3,5,7\}$ with $V(3)=V(5)=7$; by the tie rule node 3 becomes permanent.

**Iteration 5.** Scanning node 3: arcs $(3,1)$, $(3,2)$ and $(3,6)$ lead to
permanent nodes and are skipped. Arc $(3,7)$ gives $V(3)+5=12<14$, so the
label of node 7 improves to $[12,3]$. Node 5 ($V(5)=7$) becomes permanent.

**Iteration 6.** Scanning node 5: arc $(5,3)$ leads to a permanent node. Arc
$(5,7)$ gives $7+6=13\ge12$, so the label of node 7 stays $[12,3]$. Node 7 is
the only temporary node and becomes permanent.

**Iteration 7.** The node just made permanent is the destination, so the
algorithm stops with $V(7)=12$.

{numref}`tab-04-02-dijkstra` summarizes the run. Backtracking from node 7 gives
$7\leftarrow P(7)=3\leftarrow P(3)=4\leftarrow P(4)=1$, so the shortest path is
$1\to4\to3\to7$, of length $4+3+5=12$ minutes, in agreement with
{numref}`tab-04-02-paths`. The two tied labels in iteration 4 do not affect the
result: making node 5 permanent first gives the same final labels.

```{list-table} Dijkstra's algorithm from node 1 (lengths in minutes); temporary labels are listed after the scan
:name: tab-04-02-dijkstra
:header-rows: 1

* - Iteration
  - Node scanned $i$
  - Label changes
  - Temporary labels ($T$)
  - Made permanent
* - 0
  - –
  - $1:[0,-]$
  - –
  - 1
* - 1
  - 1
  - $2:[2,1]$, $4:[4,1]$
  - $2:[2,1]$, $4:[4,1]$
  - 2
* - 2
  - 2
  - $5:[7,2]$
  - $4:[4,1]$, $5:[7,2]$
  - 4
* - 3
  - 4
  - $3:[7,4]$, $6:[5,4]$
  - $3:[7,4]$, $5:[7,2]$, $6:[5,4]$
  - 6
* - 4
  - 6
  - $7:[14,6]$
  - $3:[7,4]$, $5:[7,2]$, $7:[14,6]$
  - 3 (tie with 5)
* - 5
  - 3
  - $7:[14,6]\to[12,3]$
  - $5:[7,2]$, $7:[12,3]$
  - 5
* - 6
  - 5
  - none ($13\ge12$)
  - $7:[12,3]$
  - 7
* - 7
  - –
  - stop: $i=t=7$
  - –
  - –
```

*(The initial label of the origin is sometimes written $[0,0]$ or $(0,0)$; the
second entry is not a node, and here the dash in $[0,-]$ marks that the origin
has no predecessor.)*

```{figure} figures/04-02-dijkstra-tree.svg
:name: fig-04-02-dijkstra-tree
:width: 560px
:class: dark-light
:alt: The seven-node network with a label [V, P] at every node and six thick arcs forming a tree rooted at node 1

Final labels $[V(j),P(j)]$ of {prf:ref}`ex-04-02-dijkstra`. The thick arcs
$(P(j),j)$ form the shortest-path tree rooted at node 1; the filled nodes lie
on the shortest path $1\to4\to3\to7$ of length 12.
```
:::

As the solution of {prf:ref}`ex-04-02-dijkstra` shows, the run did more than
was asked: every node received its distance from node 1,
because node 7 happened to be the farthest. The arcs $(P(j),j)$ of all
permanent nodes are drawn thick in {numref}`fig-04-02-dijkstra-tree`; they form a
tree, which the next section explains.

(sec-04-02-dijkstra-correctness)=
## Correctness, the shortest-path tree and running time

The algorithm is greedy: it declares a label final as soon as it is the
smallest temporary one, and never revisits it. The following theorem shows
why this is safe when no arc is negative.

:::{prf:theorem} Correctness of Dijkstra's algorithm
:label: thm-04-02-dijkstra

Let $c_{ij}\ge0$ for all arcs. Whenever a node $r$ is made permanent, $V(r)$
is the distance from $s$ to $r$, and following the predecessor labels from
$r$ back to $s$ traces a shortest path.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Every label $[V(j),P(j)]$ describes an actual path: $V(j)=V(P(j))+c_{P(j)j}$
at the moment the label was set, and $P(j)$ was permanent then. So $V(j)$ is
never smaller than the distance $d(j)$ from $s$ to $j$, and we only have to
show $V(r)\le d(r)$.

Use induction on the order in which nodes become permanent; the claim holds
for $s$. Let $r$ be the next node to become permanent, and assume that every
node of $S$ has its correct distance. Take any path $Q$ from $s$ to $r$. It
starts in $S$ and ends outside $S$, so it has a first arc $(u,w)$ with $u\in S$
and $w\notin S$. When $u$ was scanned, $w$ received a label with
$V(w)\le V(u)+c_{uw}=d(u)+c_{uw}$, and labels only decrease afterwards. Because
the arcs of $Q$ after $w$ have nonnegative lengths,

$$
c(Q)\ \ge\ d(u)+c_{uw}\ \ge\ V(w)\ \ge\ V(r),
$$

where the last inequality holds because $w\in T$ and $r$ has the smallest
temporary label. This is true for every path $Q$, so $V(r)\le d(r)$.
:::

The nonnegativity assumption enters at exactly one point: the part of $Q$
after the first node outside $S$ cannot make $Q$ shorter. The argument also
shows that nodes become permanent in order of nondecreasing distance, so a
run with a destination $t$ can stop as soon as $t$ is permanent.

:::{prf:proposition} Shortest-path tree
:label: prop-04-02-tree

When Dijkstra's algorithm stops after all reachable nodes are permanent, the
arcs $(P(j),j)$, $j\ne s$, form a tree rooted at $s$ in which the unique path
from $s$ to every node $j$ is a shortest path. It is called a
**shortest-path tree**.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Every node $j\ne s$ has exactly one entering tree arc, from $P(j)$, and
$P(j)$ became permanent before $j$. Following predecessors therefore visits
nodes in strictly earlier order of permanence and must end at $s$, so the arcs
contain no cycle and connect every node to $s$: with $\abs{\mathcal N}-1$ arcs
on the reachable nodes they form a tree ({prf:ref}`thm-04-01-tree`). By
{prf:ref}`thm-04-02-dijkstra` the tree path to $j$ has length $V(j)=d(j)$.
:::

A single-source problem therefore has a compact answer, one predecessor per
node. In {numref}`fig-04-02-dijkstra-tree` the tree has the arcs $(1,2)$, $(1,4)$,
$(2,5)$, $(4,3)$, $(4,6)$ and $(3,7)$. Note the difference from a minimum
spanning tree ({doc}`04-01-networks-mst`): the shortest-path tree minimizes the
length of every root-to-node path, not the total length of the tree arcs.

**Running time.** Every node is scanned once and every arc is examined once
during the scans, which costs $O(\abs{\mathcal A})$ in total. The remaining
work is finding the smallest temporary label in Step 4, $\abs{\mathcal N}-1$
times. With the labels in an array each search costs $O(\abs{\mathcal N})$,
so the algorithm runs in $O(\abs{\mathcal N}^2)$ time, the best possible for
dense networks. With the temporary labels in a binary heap, each selection and
each label decrease costs $O(\log\abs{\mathcal N})$, giving
$O(\abs{\mathcal A}\log\abs{\mathcal N})$, which is much faster for sparse
networks such as road networks. Fibonacci heaps reduce this to
$O(\abs{\mathcal A}+\abs{\mathcal N}\log\abs{\mathcal N})$
(Ahuja et al., 1993, Section 4.7). The scheme is the same as Prim's
algorithm with labels ({ref}`sec-04-01-prim`); only the update differs.

:::{prf:example} A negative arc defeats Dijkstra's algorithm
:label: ex-04-02-negative-arc

A network has the arcs $(1,2)$ of length 2, $(1,3)$ of length 3 and $(3,2)$
of length $-2$, and no cycle at all. Run Dijkstra's algorithm from node 1 and
compare the label it gives node 2 with the true distance. Does adding a
constant to all lengths, to make them nonnegative, repair the method?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Dijkstra's algorithm from node 1 ({numref}`fig-04-02-negative-arc`) labels
node 2 with $[2,1]$ and node 3 with $[3,1]$, and makes node 2 permanent
because $2<3$. Scanning node 2 changes nothing. Node 3 becomes permanent with $[3,1]$, and when it is scanned the arc $(3,2)$ leads to
a permanent node and is ignored. The algorithm reports $V(2)=2$, but the path
$1\to3\to2$ has length $3-2=1$.

The proof of {prf:ref}`thm-04-02-dijkstra` fails at the arc after the first node
outside $S$: the path $1\to3\to2$ leaves $S=\{1\}$ through node 3 with length 3
and then becomes shorter. Adding a constant to all lengths to make them
nonnegative does not help either: adding 2 gives the lengths 4, 5 and 0, and
the shortest path becomes $1\to2$ (length 4 against 5), because a constant
added to every arc penalizes paths with many arcs more than paths with few.

```{figure} figures/04-02-negative-arc.svg
:name: fig-04-02-negative-arc
:width: 360px
:class: dark-light
:alt: Three nodes; arcs 1 to 2 of length 2, 1 to 3 of length 3 and 3 to 2 of length minus 2, with the labels Dijkstra assigns

The counterexample of {prf:ref}`ex-04-02-negative-arc`. Dijkstra's algorithm makes
the label $[2,1]$ of node 2 permanent too early; the thick path $1\to3\to2$ has
length 1.
```
:::

With negative arcs but no negative cycles, single-source problems are solved
by **label-correcting** methods, which allow a label to change after it has
been used. The best known is the Bellman–Ford method (Bellman, 1958; Ford, 1956), which
applies the update $V(j)\leftarrow\min\{V(j),V(i)+c_{ij}\}$ to all arcs in
rounds; at most $\abs{\mathcal N}-1$ rounds are needed, a total of
$O(\abs{\mathcal N}\abs{\mathcal A})$ work, and a further improvement in round
$\abs{\mathcal N}$ reveals a negative cycle. It is developed as a
successive-approximation procedure of dynamic programming in
{ref}`sec-05-04-successive-approximation`.

(sec-04-02-floyd)=
## All pairs: the Floyd–Warshall algorithm

The distances between all pairs of nodes can be found by running Dijkstra's
algorithm from every node, in $O(\abs{\mathcal N}^3)$ time with arrays, when
all lengths are nonnegative. The Floyd–Warshall algorithm (Floyd, 1962;
Warshall, 1962) needs the same time, is even simpler, and allows negative arcs.
It works on the **length matrix** $\mD_0=(d_{ij}^{(0)})$ of
{ref}`sec-04-01-representations`: $d_{ii}^{(0)}=0$, $d_{ij}^{(0)}=c_{ij}$ for
arcs and $d_{ij}^{(0)}=+\infty$ for missing arcs. Number the nodes $1,\dots,
\abs{\mathcal N}$ and let

$$
d_{ij}^{(k)}=\text{length of a shortest path from } i \text{ to } j
\text{ with all intermediate nodes in } \{1,\dots,k\}.
$$ (eq-04-02-floyd-def)

For $k=0$ no intermediate node is allowed, which gives the length matrix, and
for $k=\abs{\mathcal N}$ every path is allowed, which gives the distances.

:::{prf:theorem} Floyd–Warshall recurrence
:label: thm-04-02-floyd

If $G$ has no negative cycle, then for $k=1,\dots,\abs{\mathcal N}$ and all
$i$, $j$

$$
d_{ij}^{(k)}=\min\bigl\{d_{ij}^{(k-1)},\ d_{ik}^{(k-1)}+d_{kj}^{(k-1)}\bigr\}.
$$ (eq-04-02-floyd)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Call a path from $i$ to $j$ *admissible* if its intermediate nodes lie in
$\{1,\dots,k\}$, and consider a shortest admissible path. If it does not pass
through $k$, its intermediate nodes lie in $\{1,\dots,k-1\}$ and its length
is at least $d_{ij}^{(k-1)}$. If it passes through $k$, it does so once, and its
parts from $i$ to $k$ and from $k$ to $j$ have intermediate nodes in
$\{1,\dots,k-1\}$, so its length is at least $d_{ik}^{(k-1)}+d_{kj}^{(k-1)}$.
Hence $d_{ij}^{(k)}$ is at least the right-hand side of {eq}`eq-04-02-floyd`.
Conversely, both numbers on the right are lengths of admissible routes: the
first of an admissible path, the second of a route from $i$ via $k$ to $j$. If
the two parts of that route share a node, it contains a cycle, which has
nonnegative length and can be removed (the argument of
{prf:ref}`lem-04-02-subpath`), leaving an admissible path that is no longer. So
neither number is smaller than $d_{ij}^{(k)}$, and equality holds.
:::

The update {eq}`eq-04-02-floyd` is called the **triple operation** on $i$, $k$
and $j$. Row $k$ and column $k$ of the matrix do not change in iteration $k$
(because $d_{kk}^{(k-1)}=0$), so the matrix can be overwritten in place.
To recover the paths, a **predecessor matrix** $\mP_k=(p_{ij}^{(k)})$ stores
the node before $j$ on the current best path from $i$ to $j$, exactly like the
label $P(j)$ of Dijkstra's algorithm.

:::{prf:algorithm} Floyd–Warshall algorithm
:label: alg-04-02-floyd

**Input:** a directed network with nodes $1,\dots,\abs{\mathcal N}$ and
lengths $c_{ij}$ (any sign).
**Output:** the distance matrix $\mD$ and the predecessor matrix $\mP$, or a
negative cycle.

1. **Initialize.** Set $d_{ii}=0$; $d_{ij}=c_{ij}$ and $p_{ij}=i$ for every
   arc $(i,j)$; $d_{ij}=+\infty$ and $p_{ij}=-$ for every other pair.
2. **Iterate over the intermediate node.** For $k=1,\dots,\abs{\mathcal N}$
   (outer loop), and for every pair $i$, $j$ (inner loops): if
   $d_{ik}+d_{kj}<d_{ij}$, set $d_{ij}\leftarrow d_{ik}+d_{kj}$ and
   $p_{ij}\leftarrow p_{kj}$.
3. **Check for negative cycles.** If $d_{ii}<0$ for some node $i$ (this can
   be checked after every iteration $k$), stop: a negative cycle passes
   through $i$. Otherwise $\mD$ holds all distances.
4. **Recover a path.** The shortest path from $i$ to $j$ is read backwards:
   $j\leftarrow p_{ij}\leftarrow p_{i,p_{ij}}\leftarrow\cdots\leftarrow i$.
:::

The loop over $k$ must be the **outer** loop: the theorem builds
$\mD_k$ from the complete matrix $\mD_{k-1}$. The algorithm performs
$\abs{\mathcal N}$ triple operations for each of the $\abs{\mathcal N}^2$
pairs, $O(\abs{\mathcal N}^3)$ work and $O(\abs{\mathcal N}^2)$ storage,
independent of the number of arcs. Why a negative diagonal entry signals a
negative cycle: every entry $d_{ij}$ is always the length of some route from
$i$ to $j$, so $d_{ii}<0$ means a closed route of negative length through
$i$, and such a route contains a negative cycle. Conversely, by induction on
$k$, after iteration $k$ every entry $d_{uv}$ is at most the length of every
path from $u$ to $v$ with intermediate nodes in $\{1,\dots,k\}$. If a negative
cycle $\mathcal W$ exists, let $h$ be its node with the largest number and $g$
another of its nodes. The two parts of $\mathcal W$, from $g$ to $h$ and from $h$ to $g$, are
such paths for $k=h-1$, so iteration $h$ sets
$d_{gg}\le d_{gh}+d_{hg}\le c(\mathcal W)<0$ at the latest.

:::{prf:example} Floyd–Warshall on a four-node network
:label: ex-04-02-floyd

A courier operates the six legs of {numref}`fig-04-02-floyd-network` between four
depots. The arc labels are net costs in \$ hundred; the leg from depot 1 to
depot 4 carries a paid return load, so its net cost is $-3$ (a gain of
\$300). Find the cheapest route between every pair of depots.

```{figure} figures/04-02-floyd-network.svg
:name: fig-04-02-floyd-network
:width: 330px
:class: dark-light
:alt: Four nodes in a square with six directed arcs; the arc from 1 to 4 has length minus 3

The courier network of {prf:ref}`ex-04-02-floyd`. Arc labels are net costs in
\$ hundred; the negative arc $(1,4)$ is used by four of the twelve cheapest
routes.
```
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The network has one negative arc and the cycles $1\to4\to3\to1$ (cost
$-3+2+3=2$), $1\to2\to3\to1$ (11) and $3\to4\to3$ (9), none negative, so
{prf:ref}`alg-04-02-floyd` applies. In the matrices below, row $i$ and column
$j$ refer to the route from $i$ to $j$, and boxed entries changed in the
iteration.

$$
\begin{gathered}
\mD_0=\begin{pmatrix}0&5&\infty&-3\\ \infty&0&3&\infty\\ 3&\infty&0&7\\ \infty&\infty&2&0\end{pmatrix},\\
\mP_0=\begin{pmatrix}-&1&-&1\\ -&-&2&-\\ 3&-&-&3\\ -&-&4&-\end{pmatrix}.
\end{gathered}
$$

**Iteration $k=1$** (routes through depot 1). Only rows with $d_{i1}<\infty$
(row 3) and columns with $d_{1j}<\infty$ (columns 2 and 4) can improve:
$d_{32}=d_{31}+d_{12}=3+5=8<\infty$ and $d_{34}=3+(-3)=0<7$; both new paths end
with an arc out of depot 1, so $p_{32}=p_{34}=p_{12}=p_{14}=1$.

$$
\begin{gathered}
\mD_1=\begin{pmatrix}0&5&\infty&-3\\ \infty&0&3&\infty\\ 3&\boxed{8}&0&\boxed{0}\\ \infty&\infty&2&0\end{pmatrix},\\
\mP_1=\begin{pmatrix}-&1&-&1\\ -&-&2&-\\ 3&\boxed{1}&-&\boxed{1}\\ -&-&4&-\end{pmatrix}.
\end{gathered}
$$

**Iteration $k=2$.** Column 2 is finite in rows 1 and 3, row 2 only in column
3: $d_{13}=5+3=8<\infty$ with $p_{13}=p_{23}=2$; $d_{33}$ would become
$8+3=11>0$ and is not changed.

$$
\begin{gathered}
\mD_2=\begin{pmatrix}0&5&\boxed{8}&-3\\ \infty&0&3&\infty\\ 3&8&0&0\\ \infty&\infty&2&0\end{pmatrix},\\
\mP_2=\begin{pmatrix}-&1&\boxed{2}&1\\ -&-&2&-\\ 3&1&-&1\\ -&-&4&-\end{pmatrix}.
\end{gathered}
$$

**Iteration $k=3$.** Depot 3 is reached from every depot and reaches every
depot, and four entries improve: $d_{21}=3+3=6$, $d_{24}=3+0=3$,
$d_{41}=2+3=5$ and $d_{42}=2+8=10$, with $p_{21}=p_{41}=p_{31}=3$,
$p_{24}=p_{34}=1$ and $p_{42}=p_{32}=1$. Routes from depot 1 through depot 3
do not improve; for example $d_{13}+d_{32}=8+8=16>5=d_{12}$.

$$
\begin{gathered}
\mD_3=\begin{pmatrix}0&5&8&-3\\ \boxed{6}&0&3&\boxed{3}\\ 3&8&0&0\\ \boxed{5}&\boxed{10}&2&0\end{pmatrix},\\
\mP_3=\begin{pmatrix}-&1&2&1\\ \boxed{3}&-&2&\boxed{1}\\ 3&1&-&1\\ \boxed{3}&\boxed{1}&4&-\end{pmatrix}.
\end{gathered}
$$

**Iteration $k=4$.** Through depot 4 only one route improves:
$d_{13}=d_{14}+d_{43}=-3+2=-1<8$, with $p_{13}=p_{43}=4$.

$$
\begin{gathered}
\mD_4=\begin{pmatrix}0&5&\boxed{-1}&-3\\ 6&0&3&3\\ 3&8&0&0\\ 5&10&2&0\end{pmatrix},\\
\mP_4=\begin{pmatrix}-&1&\boxed{4}&1\\ 3&-&2&1\\ 3&1&-&1\\ 3&1&4&-\end{pmatrix}.
\end{gathered}
$$

The diagonal of $\mD_4$ is zero, which confirms that there is no negative
cycle, and $\mD_4$ holds all cheapest route costs. For example, the cheapest
route from depot 4 to depot 2 costs $d_{42}=10$, that is, \$1 000; backtracking
in row 4 of $\mP_4$ gives $2\leftarrow p_{42}=1\leftarrow p_{41}=3\leftarrow
p_{43}=4$, the route $4\to3\to1\to2$ with cost $2+3+5=10$. The route from
depot 2 to depot 4 is $2\to3\to1\to4$ with cost $3+3-3=3$ (\$300), which uses
the profitable leg; the route from depot 1 to depot 3 is $1\to4\to3$, with a
net gain of \$100. A check with an independent implementation gives the same
matrix and paths.

*Negative cycle.* If the return load on leg $(1,4)$ paid more, say
$c_{14}=-6$, the cycle $1\to4\to3\to1$ would have cost $-6+2+3=-1$. Rerunning
the algorithm, iteration 1 sets $d_{34}=3-6=-3$, and iteration 3 gives
$d_{44}=d_{43}+d_{34}=2-3=-1<0$: the negative diagonal entry reveals the
cycle $4\to3\to1\to4$, and no cheapest route exists between depots that can use
it.
:::

:::{prf:remark} Sequence matrix
:label: rem-04-02-sequence-matrix

Some textbooks (Taha, 2017, Section 6.3) record in place of
the predecessor the intermediate node: a **sequence matrix** starts with
$s_{ij}=j$ for every arc and sets $s_{ij}\leftarrow k$ whenever iteration $k$
improves $d_{ij}$. Then $s_{ij}=j$ means that the arc $(i,j)$ is used
directly, and otherwise $s_{ij}=k$ says that the route passes through $k$, so
the routes from $i$ to $k$ and from $k$ to $j$ are recovered in the same way.
In {prf:ref}`ex-04-02-floyd` the final sequence matrix has rows
$(-,2,4,4)$, $(3,-,3,3)$, $(1,1,-,1)$ and $(3,3,3,-)$; for instance
$s_{42}=3$, $s_{43}=3$ (direct) and $s_{32}=1$, $s_{31}=1$ and $s_{12}=2$
(direct) rebuild $4\to3\to1\to2$. Both matrices describe the same paths; the
predecessor matrix needs no recursion, which is why it is used here.
:::

(sec-04-02-lp)=
## Linear programming formulation and node potentials

The shortest-path problem is a special minimum-cost flow problem
({ref}`sec-04-04-special-cases`): send one unit of flow from $s$ to $t$ at the
least cost, where $x_{ij}$ is the flow on arc $(i,j)$:

$$
\begin{aligned}
\min_{\vx}\ & \sum_{(i,j)\in\mathcal A}c_{ij}x_{ij}\\
\st\ & \sum_{j:(i,j)\in\mathcal A}x_{ij}-\sum_{j:(j,i)\in\mathcal A}x_{ji}
=\begin{cases}1,& i=s,\\-1,& i=t,\\0,&\text{otherwise,}\end{cases}\\
& x_{ij}\ge0\qquad\text{for all }(i,j)\in\mathcal A.
\end{aligned}
$$ (eq-04-02-lp)

The left-hand sides are the rows of $\mA\vx$, where $\mA$ is the node–arc
incidence matrix ({numref}`tab-04-01-incidence`). Every directed path from $s$ to
$t$ gives a feasible solution with $x_{ij}=1$ on its arcs and cost equal to the
path length. Conversely, a basic optimal solution of {eq}`eq-04-02-lp` is
integral because the incidence matrix is totally unimodular
({prf:ref}`thm-04-04-tu`), and when no negative cycle exists it can be
taken to be such a path. (With a negative cycle the linear program is
unbounded: flow can circulate round the cycle at a negative cost.) To solve
the single-source problem in one linear program, give $s$ the supply
$\abs{\mathcal N}-1$ and every other node the demand 1; an optimal solution
sends the flow along a shortest-path tree.

Attach a free dual variable $\pi_i$ to the conservation equation of node $i$.
The column of $x_{ij}$ has $+1$ in row $i$ and $-1$ in row $j$, so the dual of
{eq}`eq-04-02-lp` ({doc}`03-05-lp-duality`) is
$\max\{\pi_s-\pi_t:\ \pi_i-\pi_j\le c_{ij}\ \text{for all }(i,j)\}$. Because
the rows of $\mA$ add up to zero, the $\pi_i$ are determined only up to a
common constant; writing $V(j)=\pi_s-\pi_j$ (so that $V(s)=0$) turns the dual
into

$$
\begin{aligned}
\max_{V}\ & V(t)-V(s)\\
\st\ & V(j)-V(i)\le c_{ij}\qquad\text{for all }(i,j)\in\mathcal A.
\end{aligned}
$$ (eq-04-02-dual)

The dual variables are **node potentials**: each arc limits how far the
potential may rise along it. A physical picture is a set of beads (nodes) tied
together by strings (arcs) of lengths $c_{ij}$: pulling $t$ away from $s$ as
far as possible, the strings that become tight form a shortest path. Define
the **reduced length** of an arc with respect to the potentials,

$$
\bar c_{ij}=c_{ij}+V(i)-V(j).
$$ (eq-04-02-reduced-length)

Dual feasibility means $\bar c_{ij}\ge0$ for every arc, and complementary
slackness ({prf:ref}`thm-03-05-complementary-slackness`) means $\bar c_{ij}=0$ on every arc that carries flow. This gives a
certificate of optimality that does not depend on how the labels were
computed.

:::{prf:theorem} Shortest-path optimality conditions
:label: thm-04-02-optimality

Let $V(s)=0$ and let $V(j)$ be the length of some directed path from $s$ to
$j$ for every node $j$. Then the $V(j)$ are the distances from $s$ if and only
if

$$
V(j)\le V(i)+c_{ij}\qquad\text{for all }(i,j)\in\mathcal A,
$$ (eq-04-02-optimality)

that is, if all reduced lengths {eq}`eq-04-02-reduced-length` are nonnegative.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If some arc has $V(j)>V(i)+c_{ij}$, the path to $i$ followed by $(i,j)$ is
shorter than $V(j)$ (after removing cycles, as before), so $V(j)$ is not the
distance. Conversely, for any path $s=i_1,\dots,i_p=j$, adding the
inequalities {eq}`eq-04-02-optimality` along its arcs gives
$V(j)-V(s)\le c(Q)$, so no path is shorter than $V(j)$.
:::

In the network of {prf:ref}`ex-04-02-route`, the labels of
{prf:ref}`ex-04-02-dijkstra` give the reduced lengths in
{numref}`tab-04-02-reduced-costs`. All of them are nonnegative, so the labels
are the distances, and the arcs of the path $1\to4\to3\to7$ have
$\bar c_{ij}=0$, as complementary slackness requires. Indeed, the linear
program {eq}`eq-04-02-lp` and its dual {eq}`eq-04-02-dual` both have the optimal
value 12. The zero reduced lengths are exactly the arcs of the shortest-path
tree. For the single-pair problem the optimal potentials are not unique:
holding the other labels fixed, any $V(6)$ between 3 and 5 satisfies
{eq}`eq-04-02-optimality`, because node 6 is not on the optimal path, so a solver
may return different potentials for such nodes.

```{list-table} Reduced lengths $\bar c_{ij}=c_{ij}+V(i)-V(j)$ with the labels $V=(0,2,7,4,7,5,12)$ of nodes 1 to 7
:name: tab-04-02-reduced-costs
:header-rows: 1

* - Arc
  - $c_{ij}$
  - $\bar c_{ij}$
  - Arc
  - $c_{ij}$
  - $\bar c_{ij}$
* - $(1,2)$
  - 2
  - 0
  - $(3,7)$
  - 5
  - 0
* - $(1,4)$
  - 4
  - 0
  - $(4,3)$
  - 3
  - 0
* - $(2,5)$
  - 5
  - 0
  - $(4,6)$
  - 1
  - 0
* - $(3,1)$
  - 10
  - 17
  - $(5,3)$
  - 8
  - 8
* - $(3,2)$
  - 11
  - 16
  - $(5,7)$
  - 6
  - 1
* - $(3,6)$
  - 7
  - 9
  - $(6,7)$
  - 9
  - 2
```

The reduced lengths also explain the label-correcting methods: Bellman–Ford
repairs arcs with $\bar c_{ij}<0$ until none is left, and replacing $c_{ij}$
by $\bar c_{ij}$ for some feasible potentials turns a network with negative
arcs (but no negative cycles) into one with nonnegative lengths and the same
shortest paths, since every path from $i$ to $j$ changes in length by the same
amount $V(i)-V(j)$. This is how all-pairs problems on sparse networks with
negative arcs are solved with one Bellman–Ford run followed by Dijkstra's
algorithm from every node.

(sec-04-02-modelling)=
## Modeling with shortest paths

Many sequential decision problems are shortest-path problems on a network
whose nodes are **states** and whose arcs are **decisions**. The network is
often acyclic because the states are ordered in time.

**Equipment replacement.** A machine is needed for the next $n$ years. Let
node $i$, $i=1,\dots,n+1$, represent the start of year $i$ with a new machine
just bought (node $n+1$ is the end of the horizon). An arc $(i,j)$ with
$i<j$ means "buy a machine at the start of year $i$, use it for $j-i$ years
and trade it in at the start of year $j$", and its length is

$$
c_{ij}=\text{purchase price}+\sum_{a=0}^{j-i-1}\text{operating cost at age }a
-\text{trade-in value at age }j-i .
$$ (eq-04-02-replacement)

Every replacement plan is a path from node 1 to node $n+1$, and its length is
the total cost of the plan, so the cheapest plan is a shortest path. The
network is acyclic, so its arcs may have any sign. {doc}`05-02-dp-equipment-replacement`
solves this model by dynamic programming with numerical data and draws the
same network ({ref}`sec-05-02-network`).

**Most reliable route.** If arc $(i,j)$ works with probability $q_{ij}\in(0,1]$,
independently of the others, a route works with probability
$\prod q_{ij}$. Maximizing this product is the same as minimizing
$\sum(-\ln q_{ij})$, a shortest-path problem with nonnegative lengths
$-\ln q_{ij}$, which Dijkstra's algorithm solves.

Further models of this kind are worked out in {doc}`08-05-examples-networks`:
a purchasing and inventory plan ({prf:ref}`ex-08-05-purchasing`), a knapsack
problem written as a path problem ({prf:ref}`ex-08-05-knapsack-path`), and
contact chains in a social network ({prf:ref}`ex-08-05-hide-and-seek`).

(sec-04-02-longest)=
## Longest paths in acyclic networks

A **longest path** problem asks for a path of maximum total length. It is a
shortest-path problem with lengths $-c_{ij}$, so in a network with cycles of
positive length it is a shortest-path problem with negative cycles, which is
hard in general. In an **acyclic** network there are no cycles at all, and the
problem is easy: number the nodes in a topological order (every arc leads from
a lower to a higher number) and solve Bellman's equations with $\max$ in place
of $\min$, $V(j)=\max_{i:(i,j)\in\mathcal A}\{V(i)+c_{ij}\}$, one node at a
time in that order, in $O(\abs{\mathcal A})$ time. Longest paths in acyclic
networks determine the duration of a project in the critical path method
(activities as arcs, events as nodes) and solve knapsack-type allocation
problems drawn as stage–state networks; the procedure is developed in
{ref}`sec-05-04-acyclic`.

## Summary and Common Mistakes

### Summary

- A shortest-path problem asks for a least-length directed path between
  given nodes ({prf:ref}`def-04-02-shortest-path`); it is well posed when no
  negative cycle exists. Undirected edges become pairs of opposite arcs.
- Subpaths of shortest paths are shortest paths ({prf:ref}`lem-04-02-subpath`),
  so the distances satisfy Bellman's equations {eq}`eq-04-02-bellman`
  ({prf:ref}`thm-04-02-bellman`).
- Dijkstra's algorithm ({prf:ref}`alg-04-02-dijkstra`) makes labels
  $[V(j),P(j)]$ permanent in order of increasing distance; it is correct for
  $c_{ij}\ge0$ ({prf:ref}`thm-04-02-dijkstra`), runs in
  $O(\abs{\mathcal N}^2)$ or $O(\abs{\mathcal A}\log\abs{\mathcal N})$ time,
  yields a shortest-path tree ({prf:ref}`prop-04-02-tree`), and fails with
  negative arcs ({prf:ref}`ex-04-02-negative-arc`). In
  {prf:ref}`ex-04-02-dijkstra` the shortest path from node 1 to node 7 is
  $1\to4\to3\to7$ with length 12 ({numref}`fig-04-02-dijkstra-tree`).
- The Floyd–Warshall algorithm ({prf:ref}`alg-04-02-floyd`) applies the triple
  operation {eq}`eq-04-02-floyd` for $k=1,\dots,\abs{\mathcal N}$ in the outer
  loop, in $O(\abs{\mathcal N}^3)$ time; it allows negative arcs and reveals a
  negative cycle by a negative diagonal entry ({prf:ref}`ex-04-02-floyd`).
- The shortest-path problem is the linear program {eq}`eq-04-02-lp`; its dual
  {eq}`eq-04-02-dual` has node potentials $V(j)$ as variables, and labels are
  optimal exactly when all reduced lengths {eq}`eq-04-02-reduced-length` are
  nonnegative ({prf:ref}`thm-04-02-optimality`, {numref}`tab-04-02-reduced-costs`).
- Equipment replacement and reliability problems are shortest-path models
  ({ref}`sec-04-02-modelling`); longest paths are easy in acyclic networks
  ({ref}`sec-04-02-longest`).

### Common mistakes

- **Using Dijkstra's algorithm with negative arcs, or "fixing" them by adding
  a constant.** The proof of {prf:ref}`thm-04-02-dijkstra` needs $c_{ij}\ge0$, and
  {prf:ref}`ex-04-02-negative-arc` shows the failure. A constant added to every arc
  changes which path is shortest, because paths have different numbers of
  arcs. Use Bellman–Ford or Floyd–Warshall instead.
- **Making a label permanent too early, or stopping when the destination is
  first labeled.** Only the smallest temporary label may become permanent. In
  {prf:ref}`ex-04-02-dijkstra` node 7 first receives $[14,6]$; stopping there
  would report the path $1\to4\to6\to7$ of length 14 instead of 12.
- **Writing the Floyd–Warshall loops in the wrong order.** The intermediate
  node $k$ must be the outermost loop ({prf:ref}`alg-04-02-floyd`). With $k$ as
  the innermost loop, the entry $d_{ij}$ is finalized after a single pass over
  $k$ using entries $d_{ik}$, $d_{kj}$ that are not yet optimal, and some
  distances come out too large.
- **Mixing up the indices of the triple operation.** The test is
  $d_{ik}+d_{kj}<d_{ij}$: the route goes from $i$ to $k$ and then from $k$ to
  $j$. Writing $d_{ik}+d_{jk}$ uses the route from $j$ to $k$, which in a
  directed network is a different quantity.
- **Updating the predecessor matrix with $k$.** In the predecessor matrix the
  update is $p_{ij}\leftarrow p_{kj}$ (the node before $j$ on the best route
  from $k$ to $j$). Storing $k$ gives the sequence matrix of
  {prf:ref}`rem-04-02-sequence-matrix`, which must be decoded recursively.
- **Reading a path forwards from the predecessor labels.** $P(j)$ points to
  the node *before* $j$; start at the destination and backtrack to the origin.
- **Confusing a shortest-path tree with a minimum spanning tree.** A
  minimum spanning tree minimizes the total length of its edges and ignores
  any origin; a shortest-path tree minimizes every path from the root
  ({prf:ref}`prop-04-02-tree`).

## References

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows: Theory,
  algorithms, and applications*. Prentice Hall. (Chapter 4: shortest paths with
  label-setting algorithms, acyclic networks, Dijkstra's algorithm and its heap
  implementations; Chapter 5: optimality conditions, label-correcting
  algorithms, negative cycles and the all-pairs problem.)
- Bellman, R. (1957). *Dynamic programming*. Princeton University Press.
  (Chapter 3: the principle of optimality.)
- Bellman, R. (1958). On a routing problem. *Quarterly of Applied
  Mathematics, 16*(1), 87–90. <https://doi.org/10.1090/qam/102435>
- Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
  *Numerische Mathematik, 1*(1), 269–271. <https://doi.org/10.1007/BF01386390>
- Floyd, R. W. (1962). Algorithm 97: Shortest path. *Communications of the
  ACM, 5*(6), 345. <https://doi.org/10.1145/367766.368168>
- Ford, L. R., Jr. (1956). *Network flow theory* (Report P-923). RAND
  Corporation.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Section 10.3: the shortest-path problem.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Section 6.3: shortest-route problem, examples of applications,
  Dijkstra's and Floyd's algorithms, linear programming formulation.)
- Warshall, S. (1962). A theorem on Boolean matrices. *Journal of the ACM,
  9*(1), 11–12.
