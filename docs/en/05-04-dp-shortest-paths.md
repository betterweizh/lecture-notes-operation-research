# Shortest Paths by Dynamic Programming

## Overview

The staged networks of {doc}`05-01-dp-fundamentals` had a built-in stage
variable: the column of a node. A general network has no columns, may
contain cycles, and may have arcs of negative length, yet the shortest-path
problem on it is still a dynamic program. This chapter shows how. The
shortest distances satisfy Bellman's equation, a recurrence that is
circular as soon as the network has a cycle; every shortest-path method
breaks the circle by choosing an *artificial stage variable*. In an acyclic
network the stage is the position of a node in a **topological order**,
which always exists (the proof is given), and one pass of the recurrence
suffices. In a network with nonnegative lengths the stage is the number of
nodes already labeled permanently, which gives Dijkstra's method. With
negative lengths the stage is the number of arcs allowed, which gives the
**successive-approximation** (Bellman–Ford) procedures; these need at most
$n-1$ iterations and detect negative cycles. For all pairs of nodes the
stage is the set of nodes allowed as intermediate stops, which gives the
recursions of Floyd and Dantzig. Each method is stated as a DP model with
its value function, recurrence and boundary conditions, proved correct, run
on a small example and costed by counting additions and comparisons.
{doc}`04-02-shortest-paths` presents Dijkstra's algorithm
({prf:ref}`alg-04-02-dijkstra`) and the Floyd–Warshall algorithm
({prf:ref}`alg-04-02-floyd`) as labeling procedures; this chapter gives the
dynamic-programming view and does not repeat the step-by-step labeling.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- write Bellman's equation for shortest distances and explain why it needs
  a stage variable when the network has cycles;
- prove that an acyclic network has a source node and a topological order,
  find the order, and solve the forward and backward one-pass recursions;
- describe Dijkstra's method as a DP whose stage is the number of permanent
  labels, and explain why it needs nonnegative lengths;
- run the successive-approximation procedures, prove that they stop within
  $n-1$ improving iterations, and use them to detect negative cycles;
- state and run the all-pairs recursions of Floyd and Dantzig and recover
  the paths from predecessor tables; and
- count the additions and comparisons of every method and choose a method
  for a given network.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $n$ | number of nodes; the nodes are numbered $1,\dots,n$, node 1 is the origin |
| $c_{ij}$ | length of arc $(i,j)$; $c_{ij}=\infty$ if there is no arc, $c_{jj}=0$ |
| $V(j)$, $P(j)$ | length of a shortest path from node 1 to node $j$; predecessor of $j$ on it |
| $\bar V(i)$ | length of a shortest path from node $i$ to node $n$ (backward recursion) |
| $\mathcal S_k$, $i_k$ | set of the $k$ nodes nearest to node 1; the node added at stage $k$ (Dijkstra) |
| $T^{(k)}(j)$ | tentative label of node $j$ at stage $k$ of Dijkstra's method |
| $V^{(k)}(j)$, $P^{(k)}(j)$ | shortest length from 1 to $j$ with at most $k$ arcs; its predecessor (procedure 1) |
| $\tilde V^{(k)}(j)$, $\hat V^{(k)}(j)$ | values of procedures 2 and 3 (paths with at most $k-1$ decreases, reversals) |
| $D^{(k)}(i,j)$, $P^{(k)}(i,j)$ | shortest length from $i$ to $j$ with intermediate nodes in $\{1,\dots,k\}$; predecessor of $j$ (Floyd) |
| $\tilde D^{(k)}(i,j)$ | the same, restricted to origins and destinations in $\{1,\dots,k\}$ (Dantzig) |

(sec-05-04-bellman-equation)=
## Shortest paths as a dynamic program

Let $G=(\mathcal N,\mathcal A)$ be a directed network with nodes
$1,\dots,n$ and arc lengths $c_{ij}$. We do not assume that
$c_{ij}=c_{ji}$, nor that the lengths satisfy the triangle inequality
$c_{ij}+c_{jk}\ge c_{ik}$, nor (until stated) that they are nonnegative.
A missing arc has length $\infty$, so every pair of nodes may be treated as
joined by an arc, and $c_{jj}=0$. As in DP generally, we solve a whole
family of problems at once: the shortest paths from node 1 to *every* node.

A **cycle** is a path that starts and ends at the same node; a **negative
cycle** is one whose arc lengths add up to less than zero. If a negative
cycle can be reached from node 1 and can reach node $j$, a walk from 1 to
$j$ can go round the cycle as often as it likes and has no shortest length.
The methods below therefore either assume that there is no negative cycle
or detect one.

:::{prf:proposition} Bellman's equation
:label: prop-05-04-bellman-equation

Suppose every node can be reached from node 1 and no cycle has negative
length. Then the shortest distances $V(j)$ satisfy (as in
{prf:ref}`thm-04-02-bellman`)

$$
\begin{aligned}
V(1)&=0,\\
V(j)&=\min_{i\ne j}\{V(i)+c_{ij}\}\quad (j=2,\dots,n).
\end{aligned}
$$ (eq-05-04-bellman)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Without negative cycles a shortest walk can be taken to be a path (removing
a cycle does not increase the length), so the shortest distances exist and
$V(1)=0$. Let $j\ne1$ and let $i$ be the node before $j$ on a shortest path
to $j$. The part of the path up to $i$ is a shortest path to $i$: a shorter
one could replace it (principle of optimality; if the shorter path passed
through $j$, cutting it at $j$ would give a path to $j$ that is shorter
still). Hence $V(j)=V(i)+c_{ij}\ge\min_{i'}\{V(i')+c_{i'j}\}$. Conversely,
for every $i'$ a shortest path to $i'$ followed by the arc $(i',j)$ is a walk
to $j$ of length $V(i')+c_{i'j}$, which is at least $V(j)$.
:::

Equation {eq}`eq-05-04-bellman` has the form of a DP recurrence with the node
as state and the predecessor $P(j)$ as decision, but it has no stage: $V(j)$
is expressed through values $V(i)$ that may in turn depend on $V(j)$ when
the network has a cycle through $j$. It is a system of equations, not a
recipe. Each method in this chapter turns it into a recipe by adding an
**artificial stage variable** chosen so that the values of one stage depend
only on values of earlier stages ({numref}`tab-05-04-methods` at the end of
the chapter collects the choices). In the language of
{prf:ref}`def-05-01-dp-model`, the node $j$ is the state and the stage is this
artificial variable, written $k$ as in the other DP chapters. Because $V(j)$
without an index already denotes the true distance, the stage is written as
a superscript, $V^{(k)}(j)$, like the iterates $\vx^{(k)}$ of an algorithm.

(sec-05-04-acyclic)=
## Acyclic networks

A network without cycles is called **acyclic**. Staged networks such as
the lattice of {doc}`05-01-dp-fundamentals`, precedence networks of projects
and time-expanded networks are acyclic. In such networks the nodes can be
numbered so that every arc points "forwards", and that number is a natural
stage variable.

### Source nodes and topological order

A **source node** is a node without incoming arcs.

:::{prf:lemma} Acyclic networks have a source node
:label: lem-05-04-source-node

Every acyclic network with at least one node has a source node.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Start at any node and walk *backwards* along incoming arcs: from the
current node move to the tail of one of its incoming arcs, as long as the
current node has one. If the walk never stops, then after at most $n$ steps
it visits some node twice, and the arcs between the two visits, traversed
forwards, form a cycle, which is impossible. So the walk stops, and it stops
only at a node without incoming arcs.
:::

:::{prf:theorem} Topological order
:label: thm-05-04-topological-order

A directed network is acyclic if and only if its nodes can be numbered
$1,\dots,n$ so that $i<j$ for every arc $(i,j)$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If such a numbering exists, the node numbers strictly increase along every
path, so no path returns to its first node: there is no cycle.

Conversely, let the network be acyclic. By {prf:ref}`lem-05-04-source-node` it
has a source node; give it the number 1 and delete it with its arcs.
Deleting nodes cannot create a cycle, so the remaining network is acyclic
and has a source node; give it the number 2, delete it, and continue until
all nodes are numbered. When node $j$ receives its number, all its
remaining incoming arcs have been deleted, so each of its incoming arcs
$(i,j)$ comes from a node $i$ numbered earlier: $i<j$.
:::

The proof is an algorithm. If at some point the remaining network has no
source node, the network contains a cycle, so the algorithm also tests
acyclicity. When several source nodes are available, any of them may be
chosen; the order is then not unique, and every order gives the same
shortest distances. With incoming-arc counters that are decreased as arcs
are deleted, the ordering needs only bookkeeping proportional to the number
of arcs, and no additions or comparisons of lengths.

### The one-pass recursion

Number the nodes in topological order. Every arc into node $j$ comes from a
node with a smaller number, so $V(j)$ depends only on $V(1),\dots,V(j-1)$:
the node number is a stage variable.

- Value function: $V(j)$, the length of a shortest path from node 1 to
  node $j$; policy $P(j)$, the node before $j$ on such a path.
- Recurrence:

  $$
  V(j)=\min_{i<j}\{V(i)+c_{ij}\}\qquad(j=2,3,\dots,n),
  $$ (eq-05-04-acyclic-forward)

  with $c_{ij}=\infty$ when $(i,j)\notin\mathcal A$.
- Boundary condition $V(1)=0$; answers $V(j)$ for $j=1,\dots,n$.

Here the "state" at stage $j$ is trivial (there is only one, node $j$
itself), so the value function depends on the stage alone.

:::{prf:theorem} Correctness of the one-pass recursion
:label: thm-05-04-acyclic

In an acyclic network numbered in topological order, {eq}`eq-05-04-acyclic-forward`
with $V(1)=0$ gives the shortest distances from node 1, whatever the signs
of the lengths.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

By induction on $j$. Node 1 is a source, so the only path from 1 to 1 is
the empty path and $V(1)=0$. Assume that $V(1),\dots,V(j-1)$ are correct.
Every path from 1 to $j$ enters $j$ by an arc $(i,j)$ with $i<j$, and by
the principle of optimality the best such path with last arc $(i,j)$ has
length $V(i)+c_{ij}$; choosing $i$ optimally gives
{eq}`eq-05-04-acyclic-forward`. The principle applies because, in an acyclic
network, a path to $i$ followed by $(i,j)$ is always a path (it cannot
return to $j$), so no negative cycle can spoil it.
:::

The **backward** recursion computes the distances *to* the last node
instead: with $\bar V(i)$ the length of a shortest path from node $i$ to
node $n$,

$$
\begin{aligned}
\bar V(n)&=0,\\
\bar V(i)&=\min_{j>i}\{c_{ij}+\bar V(j)\}\quad(i=n-1,\dots,1).
\end{aligned}
$$ (eq-05-04-acyclic-backward)

**Operation count.** In a network with an arc from every node to every
higher-numbered node, stage $j$ of {eq}`eq-05-04-acyclic-forward` needs $j-1$
additions and $j-2$ comparisons, so the whole recursion needs

$$
\begin{aligned}
\sum_{j=2}^{n}(j-1)&=\frac{n(n-1)}{2}\ \text{additions},\\
\sum_{j=2}^{n}(j-2)&=\frac{(n-1)(n-2)}{2}\ \text{comparisons},
\end{aligned}
$$ (eq-05-04-acyclic-count)

about $n^2$ operations, not counting the ordering. In a sparse network the
count is one addition per arc and, at each node, one comparison fewer than
its number of incoming arcs.

:::{prf:example} Topological order and one-pass recursion
:label: ex-05-04-acyclic

Find the shortest paths from node $D$ to all nodes of the acyclic network
whose arcs and lengths (km) are
$D\to A$ 4, $D\to C$ 2, $C\to A$ 1, $C\to E$ 6, $A\to E$ 2, $A\to F$ 5,
$E\to F$ 2, $E\to B$ 6 and $F\to B$ 1.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Topological order.** Only $D$ has no incoming arc: number it 1 and delete
it with $D\to A$ and $D\to C$. Now only $C$ has no incoming arc (A still
has $C\to A$): number 2. Deleting $C$ leaves $A$ as the only source (3),
then $E$ (4), then $F$ (5), then $B$ (6). The order $D,C,A,E,F,B$ happens
to be unique, because at each step exactly one source was available.

**Forward recursion.** In this order,

$$
\begin{aligned}
V(D)&=0,\\
V(C)&=V(D)+2=2,\\
V(A)&=\min\{V(D)+4,\ V(C)+1\}=\min\{4,\ 3\}=3\quad(P(A)=C),\\
V(E)&=\min\{V(C)+6,\ V(A)+2\}=\min\{8,\ 5\}=5\quad(P(E)=A),\\
V(F)&=\min\{V(A)+5,\ V(E)+2\}=\min\{8,\ 7\}=7\quad(P(F)=E),\\
V(B)&=\min\{V(E)+6,\ V(F)+1\}=\min\{11,\ 8\}=8\quad(P(B)=F).
\end{aligned}
$$

The shortest path from $D$ to $B$ is $D$–$C$–$A$–$E$–$F$–$B$, 8 km. The
predecessors form the shortest-path tree highlighted in
{numref}`fig-05-04-acyclic-network`; here it is a single path, because every
node lies on the shortest route to $B$. The recursion used 9 additions (one
per arc) and 4 comparisons.

**Backward recursion.** With $\bar V$ = distance to $B$:
$\bar V(B)=0$, $\bar V(F)=1$, $\bar V(E)=\min\{6+0,\ 2+1\}=3$,
$\bar V(A)=\min\{2+3,\ 5+1\}=5$, $\bar V(C)=\min\{1+5,\ 6+3\}=6$,
$\bar V(D)=\min\{4+5,\ 2+6\}=8$, the same 8 km.

**Check.** $V(j)+\bar V(j)=8$ at every node, as it must be for nodes on a
shortest $D$–$B$ path.

```{figure} figures/05-04-acyclic-network.svg
:name: fig-05-04-acyclic-network
:width: 520px
:class: dark-light
:alt: Six-node acyclic network with arc lengths, topological numbers next to the nodes and labels V and P, with the shortest-path tree from D highlighted

The acyclic network of {prf:ref}`ex-05-04-acyclic`. Gray digits give the
topological order $D,C,A,E,F,B$; the label $[V(j),P(j)]$ beside each node
is the result of the forward recursion; highlighted arcs form the
shortest-path tree from $D$.
```
:::

Replacing $\min$ by $\max$ in {eq}`eq-05-04-acyclic-forward` gives *longest*
paths, which are well defined in an acyclic network; this is the critical
path computation of project scheduling.

(sec-05-04-dijkstra)=
## Dijkstra's method as a dynamic program

When all lengths are nonnegative, the nodes can be ordered by their
distance from node 1, and Dijkstra's method ({doc}`04-02-shortest-paths`)
finds this order while it computes the distances. Its DP form makes the
stage variable explicit: **stage $k$ is the number of nodes whose distance
is known.**

Let $\mathcal S_k$ be a set of $k$ nodes nearest to node 1 (no node outside
$\mathcal S_k$ is closer to 1 than a node inside it), with
$\mathcal S_1=\{1\}$, and let $i_k$ be the node added at stage $k$. Define
$T^{(k)}(j)$ as the length of a shortest path from 1 to $j$ whose nodes
other than $j$ all lie in $\mathcal S_k$.

- **Decision at stage $k$** ($k=2,\dots,n$): choose
  $i_k\in\argmin_{j\notin\mathcal S_{k-1}}T^{(k-1)}(j)$ and set
  $\mathcal S_k=\mathcal S_{k-1}\cup\{i_k\}$.
- **Recurrence:**

  $$
  T^{(k)}(j)=\min\bigl\{T^{(k-1)}(j),\ T^{(k-1)}(i_k)+c_{i_k j}\bigr\}\qquad (j\notin\mathcal S_k),
  $$ (eq-05-04-dijkstra)

  while $T^{(k)}(j)=T^{(k-1)}(j)$ for $j\in\mathcal S_k$.
- **Boundary conditions:** $\mathcal S_1=\{1\}$, $i_1=1$, $T^{(1)}(j)=c_{1j}$.
- **Answers:** $V(i_k)=T^{(k-1)}(i_k)$ for $k=2,\dots,n$, and $V(1)=0$.

:::{prf:theorem} Correctness of Dijkstra's recursion
:label: thm-05-04-dijkstra

If $c_{ij}\ge0$ for all arcs, then at every stage $T^{(k-1)}(i_k)=V(i_k)$,
and $T^{(k)}(j)$ has the meaning stated above for every $j\notin\mathcal S_k$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Inductive step. Assume that the labels of stage $k-1$ have the stated
meaning and that the nodes of $\mathcal S_{k-1}$ carry their true
distances. Take any path from 1 to $i_k$ and let $j$ be its first node
outside $\mathcal S_{k-1}$. The part up to $j$ has length at least
$T^{(k-1)}(j)\ge T^{(k-1)}(i_k)$ by the choice of $i_k$, and the rest has
nonnegative length, so no path to $i_k$ is shorter than $T^{(k-1)}(i_k)$:
this is the true distance. Next, a path to $j\notin\mathcal S_k$ whose
other nodes lie in $\mathcal S_k$ either avoids $i_k$, and is counted in
$T^{(k-1)}(j)$, or passes through $i_k$. If it passes through $i_k$ but
$i_k$ is not the last node before $j$, the last node before $j$ is some
$m\in\mathcal S_{k-1}$, and the part up to $m$ can be replaced by a
shortest path to $m$ inside $\mathcal S_{k-1}$ (which does not use $i_k$,
because $m$ is at least as close to 1 as $i_k$), so it is again counted in
$T^{(k-1)}(j)$. The remaining paths end with the arc $(i_k,j)$, and the
best of them has length $T^{(k-1)}(i_k)+c_{i_kj}$. This is
{eq}`eq-05-04-dijkstra`.
:::

The proof uses $c\ge0$ in one place only, and that place really fails with
a negative arc. With arcs $1\to2$ of length 2, $1\to3$ of length 3 and
$3\to2$ of length $-2$, the recursion makes node 2 permanent with label 2,
although the path 1–3–2 has length 1.

**Operation count.** Stage $k$ ($k=2,\dots,n$) scans the $n-k+1$
unlabeled nodes, which takes $n-k$ comparisons, and then updates the
$n-k$ remaining labels with one addition and one comparison each. Summing,

$$
\begin{gathered}
\frac{(n-1)(n-2)}{2}\ \text{additions},\\
(n-1)(n-2)\ \text{comparisons},
\end{gathered}
$$ (eq-05-04-dijkstra-count)

about $\tfrac32n^2$ operations (the first stage only copies the arc
lengths). With a list-processing technique that keeps the unlabeled nodes
separate, Yen showed that this bound holds even when the bookkeeping is
included (Dreyfus and Law, 1977, Ch. 4, §3A). A worked run in this DP form is in
{doc}`08-06-examples-dp-fundamentals`.

(sec-05-04-successive-approximation)=
## Successive approximation and negative arcs

When some lengths are negative, the distance order of Dijkstra's method is
lost. Arcs with negative lengths occur naturally when arc numbers are costs
and some movements earn money, or in the residual networks of minimum-cost
flow algorithms ({doc}`04-04-minimum-cost-flow`). The remedy is a stage
variable that does not depend on the lengths at all: **the number of arcs a
path may use.** This is the method of successive approximation, proposed by
Ford (1956) and by Bellman (1958) and usually called the **Bellman–Ford**
algorithm. Dreyfus and Law (1977, Ch. 4, §3A) present it together with two
faster variants, here called procedures 1, 2 and 3.

### Procedure 1: at most $k$ arcs

- Value function: $V^{(k)}(j)$ is the length of a shortest walk from node 1
  to node $j$ that uses at most $k$ arcs.

  (A walk may repeat nodes; when no negative cycle can be reached, the
  shortest walks may be taken to be paths, so the word only matters for
  the detection of negative cycles.)

- Recurrence ($k=1,2,\dots$; $j=1,\dots,n$):

  $$
  V^{(k)}(j)=\min_{i\ne j}\bigl\{V^{(k-1)}(i)+c_{ij}\bigr\},
  $$ (eq-05-04-bellman-ford)

  except that for $j=1$ the term $i=1$ (the empty path, $V^{(k-1)}(1)+0$)
  is included as well.
- Boundary conditions: $V^{(0)}(1)=0$, $V^{(0)}(j)=\infty$ for
  $j=2,\dots,n$.
- Policy: $P^{(k)}(j)$, a minimizing $i$ in {eq}`eq-05-04-bellman-ford`, the
  predecessor of $j$ on the path.

The stage is the iteration counter $k$, the state is the node $j$, and the
decision is the last arc into $j$. Each iteration is an improved
approximation of $V$; hence the name.

:::{prf:theorem} Successive approximation
:label: thm-05-04-successive

1. $V^{(k)}(j)$, defined as above, is given by {eq}`eq-05-04-bellman-ford` and
   the boundary conditions, and $V^{(k)}(j)\le V^{(k-1)}(j)$.
2. If no negative cycle can be reached from node 1, then
   $V^{(n-1)}(j)=V(j)$ for all $j$.
3. If $V^{(k)}=V^{(k-1)}$ for some $k$ (**convergence**), then
   $V^{(k')}=V^{(k)}$ for all $k'\ge k$, no negative cycle can be reached
   from node 1, and $V^{(k)}(j)=V(j)$.
4. Consequently, convergence occurs by iteration $n$ if and only if no
   negative cycle can be reached from node 1; if $V^{(n)}\ne V^{(n-1)}$,
   such a cycle exists.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

1. A walk to $j\ne1$ with at most $k$ arcs ends with some arc $(i,j)$,
   $i\ne j$, and the part before it is a walk to $i$ with at most $k-1$
   arcs; by the principle of optimality its best length is
   $V^{(k-1)}(i)$, and choosing $i$ optimally gives
   {eq}`eq-05-04-bellman-ford`. For $j=1$ the empty walk must be allowed too.
   Walks with at most $k-1$ arcs also have at most $k$ arcs, so the values
   cannot increase.
2. Without reachable negative cycles a shortest walk to any reachable node
   can be taken to be a path, and a path visits at most $n$ nodes, so it
   has at most $n-1$ arcs.
3. Iteration $k+1$ applies the same formula to the same numbers as
   iteration $k$, so it produces the same numbers again, and so on for
   ever: $V^{(k')}=V^{(k)}$ for all $k'\ge k$. If a negative cycle could be
   reached from node 1, walks around it would make some $V^{(k')}(j)$ tend
   to $-\infty$ as $k'$ grows, a contradiction. Then part 2 gives
   $V^{(k)}=V^{(n-1)}=V$ (if $k<n-1$, use $V^{(n-1)}=V^{(k)}$).
4. If no negative cycle is reachable, part 2 gives $V^{(n-1)}=V$ and, since
   $V^{(n)}$ cannot be lower than the true distances, $V^{(n)}=V^{(n-1)}$.
   The converse is part 3.
:::

In practice one stops at the first $k\le n$ with $V^{(k)}=V^{(k-1)}$; if
the values still change in iteration $n$, a negative cycle exists, and
following the predecessors $P^{(n)}$ backwards from a node whose value
changed leads into it.

**Operation count.** Each iteration needs $n-1$ additions and $n-2$
comparisons for each of the nodes $j=2,\dots,n$, and $n$ additions and
$n-1$ comparisons for node 1:

$$
\begin{gathered}
n(n-1)+1\ \text{additions and}\\
n(n-2)+1\ \text{comparisons per iteration},
\end{gathered}
$$ (eq-05-04-bf-count)

about $2n^2$ operations per iteration and about $2n^3$ in the worst case of
$n$ iterations (plus the convergence test, $n$ comparisons per
iteration).

:::{prf:example} Successive approximation with a negative arc
:label: ex-05-04-bellman-ford

A network has arcs $1\to2$ of length 2, $1\to3$ of length 5, $2\to3$ of
length 6, $3\to2$ of length $-4$, $2\to4$ of length 3 and $3\to4$ of length
1 (lengths in \$100, negative = profit).
Find the shortest paths from node 1 with procedure 1.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Iteration 1.** From the boundary values $V^{(0)}=(0,\infty,\infty,\infty)$:

$$
\begin{aligned}
V^{(1)}(1)&=\min\{0+0,\ \infty,\ \infty,\ \infty\}=0,\\
V^{(1)}(2)&=\min\{0+2,\ \infty-4,\ \infty\}=2\quad(P=1),\\
V^{(1)}(3)&=\min\{0+5,\ \infty+6,\ \infty\}=5\quad(P=1),\\
V^{(1)}(4)&=\min\{\infty,\ \infty+3,\ \infty+1\}=\infty .
\end{aligned}
$$

**Iteration 2.**
$V^{(2)}(2)=\min\{0+2,\ V^{(1)}(3)-4\}=\min\{2,1\}=1$ (via 3),
$V^{(2)}(3)=\min\{0+5,\ V^{(1)}(2)+6\}=5$,
$V^{(2)}(4)=\min\{V^{(1)}(2)+3,\ V^{(1)}(3)+1\}=\min\{5,6\}=5$ (via 2).

**Iteration 3.** Only $V^{(3)}(4)=\min\{V^{(2)}(2)+3,\ V^{(2)}(3)+1\}=\min\{4,6\}=4$
changes. **Iteration 4** reproduces iteration 3, so the procedure has
converged ({numref}`tab-05-04-bellman-ford`).

```{table} Procedure 1 on the network with a negative arc: values $V^{(k)}(j)$ and predecessors $P^{(k)}(j)$.
:name: tab-05-04-bellman-ford

| $k$ | $V^{(k)}(1)$ | $V^{(k)}(2)$ | $V^{(k)}(3)$ | $V^{(k)}(4)$ | $P^{(k)}(2)$ | $P^{(k)}(3)$ | $P^{(k)}(4)$ |
|---|---|---|---|---|---|---|---|
| 0 | 0 | $\infty$ | $\infty$ | $\infty$ | – | – | – |
| 1 | 0 | 2 | 5 | $\infty$ | 1 | 1 | – |
| 2 | 0 | 1 | 5 | 5 | 3 | 1 | 2 |
| 3 | 0 | 1 | 5 | 4 | 3 | 1 | 2 |
| 4 | 0 | 1 | 5 | 4 | 3 | 1 | 2 |
```

**Answer.** $V=(0,1,5,4)$. Tracing predecessors from node 4: $P(4)=2$,
$P(2)=3$, $P(3)=1$, so the shortest path to node 4 is 1–3–2–4 with length
$5-4+3=4$ (\$400). It uses three arcs, which is why it appears only in
iteration 3. A shortest-path method that ignored the negative arc
(Dijkstra's) would have fixed node 2 at the label 2 and node 4 at 5.

**A negative cycle.** Change the length of $2\to3$ from 6 to 3. The cycle
2–3–2 now has length $3-4=-1$. Procedure 1 gives
$V^{(1)}=(0,2,5,\infty)$, $V^{(2)}=(0,1,5,5)$, $V^{(3)}=(0,1,4,4)$ and
$V^{(4)}=(0,0,4,4)$: the values still change in iteration $4=n$, so a
negative cycle exists, and following predecessors from node 2
($P^{(4)}(2)=3$, $P^{(4)}(3)=2$) finds it.

```{figure} figures/05-04-negative-arc.svg
:name: fig-05-04-negative-arc
:width: 400px
:class: dark-light
:alt: Four-node network with a negative arc from node 3 to node 2 and the shortest-path tree 1-3, 3-2, 2-4 highlighted

Network of {prf:ref}`ex-05-04-bellman-ford` and {prf:ref}`ex-05-04-floyd`. The
arc from 3 to 2 has length $-4$. Labels $[V(j),P(j)]$ give the shortest
distances from node 1 and the predecessors; the highlighted shortest-path
tree contains the path 1–3–2–4 of length 4.
```
:::

### Procedure 2: use the newest values

When procedure 1 computes $V^{(k)}(j)$, the new values $V^{(k)}(i)$ of the
nodes $i<j$ are already known and are no larger than the old ones. Using
them at once (as the Gauss–Seidel method does for linear equations) can only
help. Along a path $1,i_1,i_2,\dots,i_m,j$, say that a **decrease** occurs
at $i_p$ if $i_p<i_{p-1}$ (the node numbers go down). The value function of
procedure 2 is $\tilde V^{(k)}(j)$, the length of a shortest path from 1 to
$j$ with at most $k-1$ decreases, and its recurrence is

$$
\tilde V^{(k)}(j)=\min\Bigl\{\min_{i<j}\bigl\{\tilde V^{(k)}(i)+c_{ij}\bigr\},\ \min_{i>j}\bigl\{\tilde V^{(k-1)}(i)+c_{ij}\bigr\}\Bigr\}
\qquad(j=1,\dots,n),
$$ (eq-05-04-procedure-2)

where for $j=1$ the term $\tilde V^{(k-1)}(1)$ is included in the second
minimum, with boundary conditions $\tilde V^{(0)}(1)=0$ and
$\tilde V^{(0)}(j)=\infty$ otherwise. An arc into $j$ from a lower node
does not create a decrease, so it may use the current iteration's value; an
arc from a higher node does, so it uses the previous iteration's value. The
stopping rule is that of procedure 1. Every path with at most $k$ arcs has
at most $k-1$ decreases, so $\tilde V^{(k)}(j)\le V^{(k)}(j)$, and with the
same work per iteration procedure 2 never needs more iterations than
procedure 1. On the network of {prf:ref}`ex-05-04-bellman-ford` it converges at
$k=3$ instead of $k=4$; the proofs and the run are in
{doc}`08-06-examples-dp-fundamentals`.

### Procedure 3: alternate the direction of the sweep

Yen's procedure (Dreyfus and Law, 1977, Ch. 4, §3A) sweeps the nodes alternately
upwards ($j=1,\dots,n$) and downwards ($j=n,\dots,1$). A **reversal** occurs
at an intermediate node of a path where the node numbers change from going
up to going down or from going down to going up. With $\hat V^{(k)}(j)$ the
length of a shortest path from 1 to $j$ with at most $k-1$ reversals, the
recurrences are

$$
\begin{aligned}
\hat V^{(k)}(j)&=\min\Bigl\{\hat V^{(k-1)}(j),\ \min_{i<j}\bigl\{\hat V^{(k)}(i)+c_{ij}\bigr\}\Bigr\}
&&(k\text{ odd};\ j=1,2,\dots,n),\\
\hat V^{(k)}(j)&=\min\Bigl\{\hat V^{(k-1)}(j),\ \min_{i>j}\bigl\{\hat V^{(k)}(i)+c_{ij}\bigr\}\Bigr\}
&&(k\text{ even};\ j=n,n-1,\dots,1),
\end{aligned}
$$ (eq-05-04-procedure-3)

with the same boundary conditions. An upward sweep extends paths by a run of
arcs that go up in node number, a downward sweep by a run that goes down,
and every extra run adds one reversal. Each sweep looks only at the arcs in
one direction, so an iteration costs

$$
\begin{gathered}
\frac{n(n-1)}{2}\ \text{additions and}\\
\frac{n(n-1)}{2}\ \text{comparisons},
\end{gathered}
$$ (eq-05-04-procedure-3-count)

about half of {eq}`eq-05-04-bf-count`. A path without negative cycles has at
most $n-2$ reversals, so at most $n$ iterations are needed, and
non-convergence by iteration $n$ again signals a negative cycle. One can
show that $\hat V^{(2k-1)}(j)\le\tilde V^{(k)}(j)$, so procedure 3 needs at
most twice as many (half-price) iterations as procedure 2, and usually far
fewer ({doc}`08-06-examples-dp-fundamentals`).

**How much faster?** Dreyfus and Law (1977, Ch. 4, §3A) report the average
number of operations of the three procedures on 50 random networks with 20
nodes, in which 1% of the arcs were negative ({numref}`tab-05-04-procedures`).
When a negative cycle is present, every procedure runs its full $n=20$
iterations, and the counts are the worst-case values
$20\,[n(n-1)+1+n(n-2)+1]=14\,840$ and $20\cdot n(n-1)=7\,600$. Without
negative cycles procedure 3 needed about 40% of the work of procedure 1
($1\,959/4\,916=0.398$) and 56% of that of procedure 2.

```{table} Average number of operations (additions plus comparisons) of the successive-approximation procedures on 50 random 20-node networks, as reported by Dreyfus and Law (1977, Table 4.1).
:name: tab-05-04-procedures

| | Procedure 1 | Procedure 2 | Procedure 3 |
|---|---|---|---|
| No negative cycle (32 networks) | 4 916 | 3 501 | 1 959 |
| Negative cycle (18 networks) | 14 840 | 14 840 | 7 600 |
```

(sec-05-04-all-pairs)=
## All pairs: the recursions of Floyd and Dantzig

To find shortest paths between *every* pair of nodes one could run a
single-origin method $n$ times. With nonnegative lengths, $n$ runs of
Dijkstra's method cost about $\tfrac32n^3$ operations. With negative
lengths, $n$ runs of successive approximation could cost $2n^4$. Two DP
formulations do the job in about $2n^3$ operations for any lengths and
detect negative cycles. Their stage variable is the **set of nodes allowed
as intermediate stops**.

### Floyd's recursion

Let $D^{(k)}(i,j)$ be the length of a shortest path from node $i$ to node
$j$ whose intermediate nodes (all nodes except $i$ and $j$) belong to
$\{1,2,\dots,k\}$. Then

$$
D^{(k)}(i,j)=\min\bigl\{D^{(k-1)}(i,j),\ D^{(k-1)}(i,k)+D^{(k-1)}(k,j)\bigr\}
\qquad(k=1,\dots,n),
$$ (eq-05-04-floyd)

with boundary condition $D^{(0)}(i,j)=c_{ij}$ (so $D^{(0)}(i,i)=0$) and
answers $D^{(n)}(i,j)$ for all $i,j$. This is the Floyd–Warshall recursion (Floyd, 1962; Warshall, 1962)
{eq}`eq-04-02-floyd`, with $D^{(k)}(i,j)=d_{ij}^{(k)}$, read as a DP: the
decision at stage $k$ is whether node $k$ is used as an intermediate node. The case $i=j$ is included in
{eq}`eq-05-04-floyd`; it detects negative cycles: if $D^{(k)}(i,i)<0$ for some
$i$, a cycle through $i$ of negative length has been found, and the
procedure stops.

:::{prf:theorem} Correctness of Floyd's recursion
:label: thm-05-04-floyd

If the procedure has not stopped before stage $k$, then $D^{(k)}(i,j)$ given
by {eq}`eq-05-04-floyd` is the length of a shortest path from $i$ to $j$ with
intermediate nodes in $\{1,\dots,k\}$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Inductive step. A shortest path from $i$ to $j$ with intermediate nodes in
$\{1,\dots,k\}$ either does not use node $k$, and then its length is
$D^{(k-1)}(i,j)$, or uses node $k$ exactly once (a path that returned to
$k$ would contain a cycle through $k$, whose length is $\ge0$ because no
negative $D^{(k-1)}(k,k)$ was found, so the cycle can be removed). In the
second case it consists of a path from $i$ to $k$ and a path from $k$ to
$j$, both with intermediate nodes in $\{1,\dots,k-1\}$; by the principle of
optimality each is shortest, with lengths $D^{(k-1)}(i,k)$ and
$D^{(k-1)}(k,j)$.
:::

**Paths.** Let $P^{(k)}(i,j)$ be the predecessor of $j$ on the path found at
stage $k$, with $P^{(0)}(i,j)=i$ for every arc $(i,j)$. When the second
term of {eq}`eq-05-04-floyd` is strictly smaller, the path from $i$ to $j$ now
ends with the path from $k$ to $j$, so $P^{(k)}(i,j)=P^{(k-1)}(k,j)$;
otherwise $P^{(k)}(i,j)=P^{(k-1)}(i,j)$. The path from $i$ to $j$ is read
backwards: $j$, $P(i,j)$, $P(i,P(i,j))$, … until $i$ is reached.

**Operation count.** At stage $k$ the entries with $i=k$ or $j=k$ do not
change ($D^{(k-1)}(k,k)=0$), so $(n-1)^2$ entries need one addition and one
comparison each:

$$
\begin{gathered}
n(n-1)^2\ \text{additions and}\\
n(n-1)^2\ \text{comparisons},\\
2n(n-1)^2\approx2n^3\ \text{in total.}
\end{gathered}
$$ (eq-05-04-floyd-count)

### Dantzig's recursion

Dantzig's procedure (as presented by Dreyfus and Law, 1977, Ch. 4, §3B) builds
the distance matrix one node at a time. After stage $k$ it knows
$\tilde D^{(k)}(i,j)$, the shortest length from $i$ to $j$ with *all* nodes
of the path, ends included, in $\{1,\dots,k\}$, for $i,j\le k$. Adding node
$k$ at stage $k$ requires the distances from the old nodes to $k$ and from
$k$ to the old nodes, a test for negative cycles through $k$, and the
improvement of the old distances by paths through $k$:

$$
\begin{aligned}
\tilde D^{(k)}(i,k)&=\min_{l=1,\dots,k-1}\bigl\{\tilde D^{(k-1)}(i,l)+c_{lk}\bigr\}&&(i=1,\dots,k-1),\\
\tilde D^{(k)}(k,j)&=\min_{l=1,\dots,k-1}\bigl\{c_{kl}+\tilde D^{(k-1)}(l,j)\bigr\}&&(j=1,\dots,k-1),\\
\tilde D^{(k)}(k,k)&=\min\Bigl\{0,\ \min_{l=1,\dots,k-1}\bigl\{\tilde D^{(k)}(k,l)+\tilde D^{(k)}(l,k)\bigr\}\Bigr\},\\
\tilde D^{(k)}(i,j)&=\min\bigl\{\tilde D^{(k-1)}(i,j),\ \tilde D^{(k)}(i,k)+\tilde D^{(k)}(k,j)\bigr\}&&(i,j=1,\dots,k-1),
\end{aligned}
$$ (eq-05-04-dantzig)

for $k=2,\dots,n$, with boundary condition $\tilde D^{(1)}(1,1)=0$ and
answers $\tilde D^{(n)}(i,j)$. If $\tilde D^{(k)}(k,k)<0$, a negative cycle
through node $k$ exists. The first two lines are one-arc extensions of
known shortest paths, which is why they use $c$ directly; their correctness
follows as for Floyd's recursion.

Dantzig's recursion needs $\tfrac12n(n-1)(2n-1)$ additions and
$\tfrac12n(n-1)(2n-3)$ comparisons, again $2n(n-1)^2$ operations in total
(the derivation is in {doc}`08-06-examples-dp-fundamentals`). Its advantage
shows when the scans in the first two lines are pruned: if
$c_{lk}\ge\tilde D^{(k)}(l,k)$ as computed from the nodes before $l$, the
arc $(l,k)$ cannot improve any $\tilde D^{(k)}(i,k)$ and its whole column of
additions can be skipped (Tabourier's modification). On 50 random 20-node
networks Dreyfus and Law (1977, Ch. 4, §3B) report averages of 9 240 operations
for the modified procedure against 14 440 $=2n(n-1)^2$ for Floyd's and
Dantzig's recursions when there was no negative cycle.

:::{prf:example} All pairs on the network with a negative arc
:label: ex-05-04-floyd

Find the shortest distances between all pairs of nodes of the network of
{prf:ref}`ex-05-04-bellman-ford` with Floyd's recursion, and check them with
Dantzig's.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Boundary matrix.** Rows are origins $i$, columns destinations $j$:

$$
D^{(0)}=\begin{pmatrix}
0&2&5&\infty\\
\infty&0&6&3\\
\infty&-4&0&1\\
\infty&\infty&\infty&0
\end{pmatrix}.
$$

**Stage 1.** No arc enters node 1 ($D^{(0)}(i,1)=\infty$ for $i\ne1$), so
nothing changes: $D^{(1)}=D^{(0)}$.

**Stage 2** (paths through node 2): $D^{(2)}(1,4)=\min\{\infty,\ 2+3\}=5$
and $D^{(2)}(3,4)=\min\{1,\ -4+3\}=-1$; the diagonal entry
$D^{(1)}(3,2)+D^{(1)}(2,3)=-4+6=2\ge0$, so no negative cycle.

**Stage 3** (paths through node 3): $D^{(3)}(1,2)=\min\{2,\ 5-4\}=1$ and
$D^{(3)}(1,4)=\min\{5,\ 5+(-1)\}=4$; the diagonal entry
$D^{(2)}(2,3)+D^{(2)}(3,2)=2\ge0$.

**Stage 4.** No arc leaves node 4, so nothing changes.

$$
\begin{gathered}
D^{(4)}=\begin{pmatrix}
0&1&5&4\\
\infty&0&6&3\\
\infty&-4&0&-1\\
\infty&\infty&\infty&0
\end{pmatrix},\\
P^{(4)}=\begin{pmatrix}
\text{–}&3&1&2\\
\text{–}&\text{–}&2&2\\
\text{–}&3&\text{–}&2\\
\text{–}&\text{–}&\text{–}&\text{–}
\end{pmatrix}.
\end{gathered}
$$

Row 1 repeats the single-origin result of {prf:ref}`ex-05-04-bellman-ford`.
The path from 1 to 4 is recovered backwards: $P(1,4)=2$, $P(1,2)=3$,
$P(1,3)=1$, giving 1–3–2–4. The infinite entries are pairs with no path
(nothing reaches node 1, nothing leaves node 4).

**Dantzig's recursion.** Stage 2 gives the $2\times2$ matrix with rows
$(0,2)$ and $(\infty,0)$. Stage 3 adds node 3:
$\tilde D(1,3)=\min\{0+5,\ 2+6\}=5$, $\tilde D(2,3)=\min\{\infty+5,\ 0+6\}=6$,
$\tilde D(3,1)=\infty$, $\tilde D(3,2)=\min\{\infty+2,\ -4+0\}=-4$,
$\tilde D(3,3)=\min\{0,\ \infty,\ -4+6\}=0$, and then
$\tilde D(1,2)=\min\{2,\ 5-4\}=1$. Stage 4 adds node 4:
$\tilde D(1,4)=\min\{\infty,\ 1+3,\ 5+1\}=4$,
$\tilde D(2,4)=\min\{\infty,\ 0+3,\ 6+1\}=3$,
$\tilde D(3,4)=\min\{\infty,\ -4+3,\ 0+1\}=-1$, and node 4 has no outgoing
arcs. The final matrix equals $D^{(4)}$.
:::

(sec-05-04-comparison)=
## Choosing a method

All methods of this chapter are the same DP, Bellman's equation
{eq}`eq-05-04-bellman`, made computable by different stage variables.
{numref}`tab-05-04-methods` compares them.

```{table} Shortest-path methods as dynamic programs: stage variable, requirements and number of operations (additions plus comparisons) for a network with $n$ nodes and arcs between all pairs.
:name: tab-05-04-methods

| Method | Stage variable | Requirement | Operations |
|---|---|---|---|
| One-pass recursion {eq}`eq-05-04-acyclic-forward` | position in a topological order | acyclic network, any lengths | $\approx n^2$ (one per arc in sparse networks) |
| Dijkstra {eq}`eq-05-04-dijkstra` | number of permanent labels | $c_{ij}\ge0$ | $\approx\tfrac32n^2$ |
| Procedure 1 {eq}`eq-05-04-bellman-ford` | maximum number of arcs | detects negative cycles | $\approx2n^2$ per iteration, at most $n$ iterations |
| Procedures 2 and 3 | maximum number of decreases / reversals | detect negative cycles | as procedure 1 / half per iteration; fewer in practice |
| $n\times$ Dijkstra | (as Dijkstra) | $c_{ij}\ge0$; all pairs | $\approx\tfrac32n^3$ |
| Floyd {eq}`eq-05-04-floyd`, Dantzig {eq}`eq-05-04-dantzig` | set of allowed intermediate nodes | all pairs; detect negative cycles | $2n(n-1)^2\approx2n^3$ |
```

In short: order the nodes and use one pass if the network is acyclic; use
Dijkstra's method when all lengths are nonnegative; use successive
approximation (preferably procedure 3) when some lengths are negative or a
negative cycle must be detected; use Floyd's or Dantzig's recursion (or
Dijkstra's method $n$ times for nonnegative lengths) when distances between
all pairs are needed. Implementations are the subject of
{doc}`09-04-code-networks`.

## Summary and Common Mistakes

### Summary

- Shortest distances satisfy Bellman's equation
  ({prf:ref}`prop-05-04-bellman-equation`, {eq}`eq-05-04-bellman`), a recurrence
  without a stage; each method adds an artificial stage variable
  ({numref}`tab-05-04-methods`).
- Every acyclic network has a source node and a topological order
  ({prf:ref}`lem-05-04-source-node`, {prf:ref}`thm-05-04-topological-order`); in
  that order one pass of {eq}`eq-05-04-acyclic-forward` or
  {eq}`eq-05-04-acyclic-backward` gives all distances, for any signs of the
  lengths ({prf:ref}`thm-05-04-acyclic`, {prf:ref}`ex-05-04-acyclic`).
- Dijkstra's method is a DP whose stage is the number of permanent labels
  ({eq}`eq-05-04-dijkstra`, {prf:ref}`thm-05-04-dijkstra`); it needs
  nonnegative lengths and about $\tfrac32n^2$ operations.
- Successive approximation ({eq}`eq-05-04-bellman-ford`) handles negative
  lengths, converges within $n$ iterations exactly when no negative cycle
  is reachable, and detects such cycles ({prf:ref}`thm-05-04-successive`,
  {prf:ref}`ex-05-04-bellman-ford`); procedures 2 and 3
  ({eq}`eq-05-04-procedure-2`, {eq}`eq-05-04-procedure-3`) reach the same values
  with less work.
- Floyd's and Dantzig's recursions ({eq}`eq-05-04-floyd`, {eq}`eq-05-04-dantzig`)
  give all-pairs distances and predecessor tables in about $2n^3$
  operations and detect negative cycles through negative diagonal entries
  ({prf:ref}`thm-05-04-floyd`, {prf:ref}`ex-05-04-floyd`).

### Common mistakes

- **Running the one-pass recursion in the wrong order.** It is correct only
  if every node is processed after all its predecessors; number the nodes
  topologically first.
- **Using Dijkstra's method with negative lengths.** A label made permanent
  may be wrong; use successive approximation instead.
- **Stopping successive approximation too early or too late.** Stop at the
  first iteration with no change; if values still change in iteration $n$,
  report a negative cycle instead of distances.
- **Dropping the empty path at the origin.** In {eq}`eq-05-04-bellman-ford`
  the value $V^{(k-1)}(1)$ must be kept for $j=1$; otherwise $V^{(k)}(1)$
  becomes the length of a cycle through node 1.
- **Ignoring the diagonal in Floyd's recursion.** The entries $D^{(k)}(i,i)$
  are the negative-cycle test; they must be computed, not set to zero.
- **Reading Floyd paths forwards.** $P(i,j)$ is the predecessor of $j$;
  paths are traced from the destination back to the origin.

## References

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows:
  Theory, algorithms, and applications*. Prentice Hall. (Ch. 4: shortest
  paths in acyclic networks and Dijkstra's algorithm; Ch. 5: label-correcting
  algorithms, negative cycles and all-pairs shortest paths.)
- Bellman, R. (1958). On a routing problem. *Quarterly of Applied
  Mathematics, 16*(1), 87–90. <https://doi.org/10.1090/qam/102435>
- Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
  *Numerische Mathematik, 1*(1), 269–271.
  <https://doi.org/10.1007/BF01386390>
- Dreyfus, S. E., & Law, A. M. (1977). *The art and theory of dynamic
  programming*. Academic Press. (Ch. 4: acyclic networks, Dijkstra's
  procedure as a DP, procedures 1–3 for negative arcs, the procedures of
  Floyd and Dantzig and its modification, operation counts and Tables 4.1
  and 4.2.)
- Floyd, R. W. (1962). Algorithm 97: Shortest path. *Communications of the
  ACM, 5*(6), 345. <https://doi.org/10.1145/367766.368168>
- Ford, L. R., Jr. (1956). *Network flow theory* (Report P-923). RAND
  Corporation.
- Warshall, S. (1962). A theorem on Boolean matrices. *Journal of the ACM,
  9*(1), 11–12.
