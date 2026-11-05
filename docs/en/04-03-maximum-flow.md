# Maximum Flow

## Overview

Oil moves from refineries to terminals through pipelines, cars move from
suburbs to the city center through roads, data move through a computer
network. Each link can carry only so much per unit of time, and a planner
asks how much can be sent in total from one place to another. This is the
**maximum-flow problem**. The chapter states it as a linear program, shows
how to bring real networks (two-way links, several sources, node limits)
into that form, and then solves it without the simplex method. The key
device is the **residual network**, which records how much more flow each
arc can take and how much can be pushed back. An **augmenting path** in the
residual network increases the flow; the **Ford–Fulkerson algorithm**
repeats this until no augmenting path is left. Why the result is optimal is
explained by **cuts**: every cut gives an upper bound on the flow, and the
**max-flow min-cut theorem** says that the best bound is attained. The same
argument shows that integer capacities give an integer maximum flow, which
turns the method into a tool for combinatorial problems such as matching
workers to jobs. The chapter builds on {doc}`04-01-networks-mst` (networks,
directed paths, cycles and the node–arc incidence matrix) and
{doc}`04-02-shortest-paths` (labels and predecessors); the LP view uses
{doc}`03-01-lp-model-building` and the duality of {doc}`03-05-lp-duality`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- formulate a maximum-flow problem as a linear program and transform
  undirected links, several sources and sinks, and node capacities into the
  standard form;
- construct the residual network of a feasible flow and find augmenting
  paths in it with the labeling procedure;
- run the Ford–Fulkerson algorithm by hand and recover feasible arc flows
  from the final residual network;
- compute the capacity of a cut, prove weak duality and the max-flow
  min-cut theorem, and read a minimum cut from the last labeling;
- explain the integrality property, the running time of augmenting-path
  methods and the improvement of Edmonds and Karp; and
- model a bipartite matching problem as a maximum flow.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $u_{ij}$ | capacity of arc $(i,j)$, the largest flow it may carry (flow units per unit of time) |
| $x_{ij}$, $\vx$ | flow on arc $(i,j)$; flow vector |
| $s$, $t$ | source and sink |
| $v$ | value of a flow: net flow out of $s$ |
| $r_{ij}$ | residual capacity of $(i,j)$ with respect to a flow $\vx$ |
| $G(\vx)$ | residual network of the flow $\vx$ |
| $\mathcal P$, $\Delta$ | augmenting path and the amount sent along it |
| $[\gamma(j),P(j)]$ | label of node $j$ in the labeling procedure: residual capacity of the arc by which $j$ was reached, and the node it was reached from |
| $D$ | set of deleted nodes in the labeling procedure |
| $w_i$ | capacity of node $i$ (largest flow through it) |
| $[S,\bar S]$, $u[S,\bar S]$ | $s$–$t$ cut ($s\in S$, $t\in\bar S=\mathcal N\setminus S$) and its capacity |
| $\pi_i$, $y_{ij}$ | dual variables of the conservation equation of node $i$ and of the capacity constraint of arc $(i,j)$ |

(sec-04-03-problem)=
## The maximum-flow problem

A **capacitated network** is a directed network $G=(\mathcal N,\mathcal A)$
in which every arc $(i,j)$ has a capacity $u_{ij}\ge0$, the largest amount of
flow that may pass along it per unit of time. Two nodes are distinguished: a
**source** $s$, where flow enters, and a **sink** $t$, where it leaves. All
other nodes only pass flow on.

:::{prf:definition} Feasible flow and maximum-flow problem
:label: def-04-03-max-flow

A **feasible flow** is a vector $\vx=(x_{ij})_{(i,j)\in\mathcal A}$ that
satisfies

- the **capacity constraints** $0\le x_{ij}\le u_{ij}$ for every arc, and
- **flow conservation** at every node $i\ne s,t$: what flows in flows out.

Its **value** $v$ is the net flow out of the source, which equals the net
flow into the sink. The **maximum-flow problem** asks for a feasible flow of
largest value.
:::

Writing "flow out minus flow in" for each node gives the linear program

$$
\begin{aligned}
\max_{\vx,\,v}\ & v \\
\st\ & \sum_{j:(i,j)\in\mathcal A}x_{ij}-\sum_{j:(j,i)\in\mathcal A}x_{ji}=
\begin{cases} v, & i=s,\\ -v, & i=t,\\ 0, & i\ne s,t,\end{cases} \\
& 0\le x_{ij}\le u_{ij}, \quad (i,j)\in\mathcal A.
\end{aligned}
$$ (eq-04-03-maxflow-lp)

The left-hand side of the conservation constraint of node $i$ is row $i$ of
$\mA\vx$, where $\mA$ is the node–arc incidence matrix of
{ref}`sec-04-01-representations`; the model is therefore
$\max\{v:\mA\vx=v(\ve_s-\ve_t),\ \vzero\le\vx\le\vu\}$. Adding the
conservation equations of all nodes gives $0=0$ (each $x_{ij}$ appears once
with $+$ and once with $-$), so one equation is redundant; this is why the
net outflow of $s$ and the net inflow of $t$ are the same number $v$. The
model is a linear program and could be handed to the simplex method. Its
special structure, however, allows a much more direct method, and that
method also shows *why* the optimum is what it is.

(sec-04-03-modelling)=
### Bringing a network into standard form

Real networks rarely arrive in the form {eq}`eq-04-03-maxflow-lp`. Four
transformations are used constantly.

**Undirected links.** A pipe or a road that can be used in either direction
with capacity $u$ becomes two arcs $(i,j)$ and $(j,i)$, each with capacity
$u$. This allows $x_{ij}=x_{ji}=u$ at the same time, which the link could not
carry. It does no harm: if both are positive, subtracting
$\min\{x_{ij},x_{ji}\}$ from both changes neither the conservation equations
(only the difference $x_{ij}-x_{ji}$ enters them) nor the value, and leaves
flow in one direction only, at most $u$.

**Several sources or sinks.** Add a **super-source** $s$ with an arc
$(s,i)$ to every original source $i$ and a **super-sink** $t$ with an arc
$(j,t)$ from every original sink $j$. The capacity of $(s,i)$ is the supply
limit of source $i$ (infinite if there is none), and similarly for the sink
arcs. A flow from $s$ to $t$ in the enlarged network is a flow from the
original sources to the original sinks, and conversely.

**Node capacities.** If at most $w_i$ units may pass through node $i$ (a
pumping station), split $i$ into an entry node $i'$ and an exit node $i''$
joined by an arc $(i',i'')$ of capacity $w_i$; arcs that entered $i$ now
enter $i'$, arcs that left $i$ now leave $i''$.

**Infinite capacities.** An arc without a limit gets $u_{ij}=\infty$. If an
algorithm needs a finite number, any value larger than the total capacity of
the arcs leaving $s$ will do, since no flow can exceed that total.

:::{prf:example} Oil pipeline network
:label: ex-04-03-pipeline

Three refineries (nodes 1, 2, 3) send crude oil through pumping stations
(nodes 4, 5, 6) to two distribution terminals (nodes 7, 8). The capacities,
in thousand barrels per day, are:

- two-way pipelines: 1–4: 20, 3–5: 15, 4–5: 20, 4–6: 10, 5–6: 30;
- one-way pipelines: $2\to4$: 10, $2\to5$: 20, $4\to7$: 10, $5\to8$: 30,
  $6\to7$: 50, $6\to8$: 20.

Find the largest total amount that the refineries can deliver to the
terminals per day. The refineries have no production limit.

*(The capacity of $5\to8$ is assumed to be 30 here; with a smaller capacity
the answer changes, see the solution.)*
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Model.** Add a super-source $s=0$ with arcs $(0,1)$, $(0,2)$, $(0,3)$ and a
super-sink $t=9$ with arcs $(7,9)$, $(8,9)$, all of infinite capacity, and
replace every two-way pipeline by two opposite arcs of the same capacity
({numref}`fig-04-03-pipeline`). The network then has 10 nodes and
$11+2\cdot5=21$ arcs. The capacity constraints include
$0\le x_{14}\le20$, $0\le x_{41}\le20$, $0\le x_{24}\le10$ and
$x_{01}\ge0$ (no upper bound), and the conservation equations read, for
example,

$$
\begin{aligned}
x_{01}+x_{02}+x_{03} &= v &&\text{(source 0)}\\
x_{14}-(x_{01}+x_{41}) &= 0 &&\text{(node 1)}\\
(x_{53}+x_{54}+x_{56}+x_{58})-(x_{25}+x_{35}+x_{45}+x_{65}) &= 0 &&\text{(node 5)}\\
-(x_{79}+x_{89}) &= -v &&\text{(sink 9)}.
\end{aligned}
$$

```{figure} figures/04-03-pipeline.svg
:name: fig-04-03-pipeline
:width: 560px
:class: dark-light
:alt: Pipeline network from super-source 0 through refineries 1 to 3, pumping stations 4 to 6 and terminals 7 and 8 to super-sink 9, with capacities and a dashed cut line

The pipeline network of {prf:ref}`ex-04-03-pipeline` with capacities in thousand
barrels per day. Lines without arrow heads are two-way pipelines (two
opposite arcs of the same capacity). The dashed line separates
$S=\{0,1,2,3\}$ from the other nodes; the arcs it crosses from left to right
have total capacity $20+10+20+15=65$.
```

**An upper bound.** Whatever the flow, everything the terminals receive must
leave the refineries through the arcs $1\to4$, $2\to4$, $2\to5$ and $3\to5$,
whose capacities add up to $65$. So no flow can exceed 65. The terminals give
another bound: their incoming capacity is $10+50+30+20=110$. Bounds of this
kind come from cuts ({ref}`sec-04-03-cuts`), and the best of them turns out to
be exact.

**A flow of value 65.** The flow

| Arc | $(0,1)$ | $(0,2)$ | $(0,3)$ | $(1,4)$ | $(2,4)$ | $(2,5)$ | $(3,5)$ | $(4,5)$ |
|---|---|---|---|---|---|---|---|---|
| $x_{ij}$ | 20 | 30 | 15 | 20 | 10 | 20 | 15 | 10 |

| Arc | $(4,6)$ | $(4,7)$ | $(5,6)$ | $(5,8)$ | $(6,7)$ | $(7,9)$ | $(8,9)$ |
|---|---|---|---|---|---|---|---|
| $x_{ij}$ | 10 | 10 | 15 | 30 | 25 | 35 | 30 |

(all other arcs 0) is feasible: at node 4, $20+10$ in and $10+10+10$ out; at
node 5, $20+15+10$ in and $30+15$ out; at node 6, $10+15$ in and $25$ out; at
node 7, $10+25$ in and $35$ out. Its value is $20+30+15=65$, equal to the
upper bound, so it is a maximum flow: the refineries can deliver
65 thousand barrels per day, and the refinery pipelines, not the terminals,
are the bottleneck. Solving {eq}`eq-04-03-maxflow-lp` with an LP solver gives
the same value.

**Sensitivity to the assumed capacity.** With $5\to8$ of capacity
$u_{58}\le15$ the answer is $50+u_{58}$; for every $u_{58}\ge15$ the refinery
bound 65 is binding and the maximum flow is 65.
:::

The upper bound in the solution of {prf:ref}`ex-04-03-pipeline` has a second
side: why is this a *combinatorial* problem rather than just an LP? The answer
of the following sections is that optimal flows can be built from finitely many paths, that
optimality is certified by one of finitely many cuts, and that integer data
give integer flows. The algorithms exploit this and never need a simplex
tableau.

(sec-04-03-residual)=
## Residual networks

Suppose a feasible flow $\vx$ is known. How can it be improved? Along an arc
$(i,j)$ one may send $u_{ij}-x_{ij}$ more units, and one may also *cancel*
up to $x_{ij}$ units already sent, which has the same effect on conservation
as sending flow from $j$ to $i$. The residual network collects both
possibilities.

To keep the formulas short, assume that for every arc $(i,j)$ the reverse arc
$(j,i)$ also exists; where it does not, add it with $u_{ji}=0$. A zero
capacity forces $x_{ji}=0$, so the transformed network has the same
feasible flows, and hence the same maximum flow, as the original one.

:::{prf:definition} Residual network
:label: def-04-03-residual

Let $\vx$ be a feasible flow. The **residual capacity** of the ordered pair
$(i,j)$ is

$$
r_{ij}=u_{ij}-x_{ij}+x_{ji},
$$ (eq-04-03-residual)

the amount by which the flow from $i$ to $j$ can be increased: by filling
$(i,j)$ up to its capacity and by canceling the flow on $(j,i)$. The
**residual network** $G(\vx)$ has the same nodes and an arc $(i,j)$ with
capacity $r_{ij}$ for every pair with $r_{ij}>0$.
:::

Two cases are worth remembering.

- **Single arc** ($u_{ji}=0$, so $x_{ji}=0$): $r_{ij}=u_{ij}-x_{ij}$ is the
  unused capacity and $r_{ji}=x_{ij}$ is the flow that can be pushed back.
- **Opposite arcs** (both capacities positive): $r_{ij}=u_{ij}-x_{ij}+x_{ji}$
  and $r_{ji}=u_{ji}-x_{ji}+x_{ij}$, the formula {eq}`eq-04-03-residual` for both
  directions.

In both cases

$$
r_{ij}+r_{ji}=u_{ij}+u_{ji},
$$ (eq-04-03-residual-sum)

so pushing flow from $i$ to $j$ moves residual capacity from $(i,j)$ to
$(j,i)$ and never creates or destroys it. For the zero flow,
$G(\vzero)$ is the original network.

:::{prf:example} Residual network of a given flow
:label: ex-04-03-residual

A network has source 1, sink 4 and arcs $(1,2)$, $(1,3)$, $(3,2)$, $(2,4)$,
$(3,4)$ with capacities 4, 2, 3, 5 and 1. The flow $\vx^{(0)}$ sends 2 units along $1\to2\to4$:
$x_{12}=x_{24}=2$, $x_{13}=x_{32}=x_{34}=0$. Construct $G(\vx^{(0)})$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

{numref}`fig-04-03-residual-example`(a) shows the network and the flow. Every
arc here has no opposite arc, so the single-arc case applies:

- $(1,2)$: $r_{12}=4-2=2$ forward and $r_{21}=2$ backward;
- $(2,4)$: $r_{24}=5-2=3$ and $r_{42}=2$;
- $(1,3)$, $(3,2)$, $(3,4)$ carry no flow: $r_{13}=2$, $r_{32}=3$,
  $r_{34}=1$, and the backward residuals are 0.

$G(\vx^{(0)})$ therefore has the seven arcs drawn in
{numref}`fig-04-03-residual-example`(b). Arc $(4,2)$ says that up to 2 units
arriving at 4 via 2 could be diverted; arc $(2,1)$ that up to 2 units could
be returned to the source.

```{figure} figures/04-03-residual-example.svg
:name: fig-04-03-residual-example
:width: 560px
:class: dark-light
:alt: Left, the diamond network with flow over capacity labels; right, its residual network with forward and backward arcs

(a) The network of {prf:ref}`ex-04-03-residual` with labels $x_{ij}/u_{ij}$ for
the flow $\vx^{(0)}$. (b) The residual network $G(\vx^{(0)})$; the backward arcs
$(2,1)$ and $(4,2)$ come from the 2 units on $1\to2\to4$.
```
:::

(sec-04-03-augmenting)=
## Augmenting paths and the Ford–Fulkerson algorithm

:::{prf:definition} Augmenting path
:label: def-04-03-augmenting-path

An **augmenting path** (breakthrough path) for a feasible flow $\vx$ is a
directed path $\mathcal P$ from $s$ to $t$ in the residual network $G(\vx)$.
Its **residual capacity** is $\Delta=\min\{r_{ij}:(i,j)\in\mathcal P\}>0$.
:::

**Augmenting.** Sending $\Delta$ units along $\mathcal P$ means: for every
arc $(i,j)$ of $\mathcal P$, set $r_{ij}\leftarrow r_{ij}-\Delta$ and
$r_{ji}\leftarrow r_{ji}+\Delta$. In terms of flows, first cancel up to
$\Delta$ units on $(j,i)$ and put the rest on $(i,j)$. The result is again a
feasible flow: on each pair the change stays within the residual capacity,
so $0\le x_{ij}\le u_{ij}$ still holds; every intermediate node of the path
receives $\Delta$ more and sends $\Delta$ more; and the net outflow of $s$
grows by $\Delta$. So **the value increases by $\Delta>0$**.

### Finding a path: the labeling procedure

An augmenting path is a directed $s$–$t$ path in $G(\vx)$, so any search
will find one (breadth-first, depth-first, …). The procedure below is a
depth-first search that always follows the arc with the largest residual
capacity and backs up from dead ends. Each reached node $j$ receives a
**label** $[\gamma(j),P(j)]$: the residual capacity of the arc by which it
was reached and the node it was reached from.

:::{prf:algorithm} Labeling procedure for an augmenting path
:label: alg-04-03-labelling

**Input:** residual capacities $r_{ij}$, source $s$, sink $t$.
**Output:** an augmenting path, or the conclusion that none exists together
with the set $S$ of nodes reachable from $s$.

1. **Initialize.** Label $s$ with $[\infty,-]$. Set $i=s$ and the set of
   deleted nodes $D=\varnothing$.
2. **Scan.** Let $J_i$ be the set of nodes $j$ that are neither labeled nor
   deleted and have $r_{ij}>0$. If $J_i=\varnothing$, go to Step 4.
3. **Label.** Choose $k\in J_i$ with the largest $r_{ik}$ (ties: smallest
   index) and label it $[r_{ik},i]$. If $k=t$, stop: tracing the
   predecessors $P(\cdot)$ back from $t$ gives an augmenting path, whose
   residual capacity is the smallest $\gamma$ on it. Otherwise set $i=k$ and
   go to Step 2.
4. **Backtrack.** If $i=s$, stop: there is no augmenting path, and
   $S=\{s\}\cup D$ is the set of nodes reachable from $s$. Otherwise delete
   $i$ (remove its label, add it to $D$), set $i=P(i)$ and go to Step 2.
:::

A node is deleted only when every residual arc leaving it leads to labeled
or deleted nodes, so a deleted node cannot lead to $t$ during this search.
Each node is labeled at most once and deleted at most once. With adjacency
lists sorted by residual capacity, each arc is examined a bounded number of
times, so one search takes $O(\abs{\mathcal A}\log\abs{\mathcal N})$ time
including the sorting, and a plain depth-first search that takes the first
candidate instead of the largest one takes $O(\abs{\mathcal A})$. When the search fails, the
deleted nodes are exactly the nodes other than $s$ reachable from $s$, which
is what the optimality proof below needs.

### The algorithm

:::{prf:algorithm} Ford–Fulkerson algorithm
:label: alg-04-03-ford-fulkerson

**Input:** a network with capacities $u_{ij}$ (reverse arcs added with
capacity 0), source $s$, sink $t$.
**Output:** a maximum flow $\vx^*$, its value $v^*$ and a minimum cut.

1. **Initialize.** Set $v=0$ and $r_{ij}=u_{ij}$ for every arc (the zero
   flow).
2. **Find a path.** Apply {prf:ref}`alg-04-03-labelling` to the current
   residual capacities. If it finds no path, go to Step 4.
3. **Augment.** Let $\Delta$ be the smallest label $\gamma(j)$ on the path.
   For every arc $(i,j)$ of the path set $r_{ij}\leftarrow r_{ij}-\Delta$ and
   $r_{ji}\leftarrow r_{ji}+\Delta$; set $v\leftarrow v+\Delta$. Go to Step 2.
4. **Stop.** $v$ is the maximum flow value. Recover the arc flows by
   {eq}`eq-04-03-recover` and take the set $S$ returned by the last search: the
   cut $[S,\bar S]$ is a minimum cut.
:::

:::{prf:example} Maximum flow on the diamond network
:label: ex-04-03-diamond

Find a maximum flow from node 1 to node 4 in the network of
{prf:ref}`ex-04-03-residual` (capacities $u_{12}=4$, $u_{13}=2$, $u_{32}=3$,
$u_{24}=5$, $u_{34}=1$), starting from the zero flow, and find a minimum cut.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The residual networks of the iterations are drawn in
{numref}`fig-04-03-ff-iterations`.

**Iteration 1** (residual network $G(\vx^{(0)})$ = original network). From 1 the
candidates are 2 ($r_{12}=4$) and 3 ($r_{13}=2$); label 2 with $[4,1]$. From
2 the only candidate is 4 ($r_{24}=5$); label 4 with $[5,2]$. The path is
$1\to2\to4$ with $\Delta=\min\{4,5\}=4$. Update $r_{12}=0$, $r_{21}=4$,
$r_{24}=1$, $r_{42}=4$; $v=4$.
*(The residual capacities of $(1,2)$ and $(2,4)$ become 0 and 1, those of
$(2,1)$ and $(4,2)$ become 4 and 4; an earlier version of this example lists
$(2,4)$ twice.)*

**Iteration 2** (residual network $G(\vx^{(1)})$). From 1 only 3 is reachable
($r_{12}=0$): label 3 with $[2,1]$. From 3 the candidates are 2
($r_{32}=3$) and 4 ($r_{34}=1$); the larger residual wins, label 2 with
$[3,3]$. From 2 the only unlabeled candidate is 4 ($r_{24}=1$); label it
$[1,2]$. Path $1\to3\to2\to4$ with $\Delta=\min\{2,3,1\}=1$. Update
$r_{13}=1$, $r_{31}=1$, $r_{32}=2$, $r_{23}=1$, $r_{24}=0$, $r_{42}=5$;
$v=5$.

**Iteration 3** (residual network $G(\vx^{(2)})$). Label 3 with $[1,1]$; from 3
the candidates are 2 ($r_{32}=2$) and 4 ($r_{34}=1$), so label 2 with
$[2,3]$. From 2 nothing new is reachable: $r_{24}=0$, and $r_{21}=4$,
$r_{23}=1$ lead to labeled nodes. **Backtrack**: delete node 2 and return to
3. Now 4 is the only candidate: label it $[1,3]$. Path $1\to3\to4$ with
$\Delta=\min\{1,1\}=1$. Update $r_{13}=0$, $r_{31}=2$, $r_{34}=0$,
$r_{43}=1$; $v=6$.

**Iteration 4** (residual network $G(\vx^{(3)})$). Both $r_{12}=0$ and
$r_{13}=0$: the source has no residual arc, and the search stops at once with
$S=\{1\}$.

**Result.** The maximum flow value is $v^*=4+1+1=6$. The arc flows follow from
$x_{ij}=\max\{0,u_{ij}-r_{ij}\}$ ({ref}`sec-04-03-recover`): $x_{12}=4-0=4$,
$x_{13}=2-0=2$, $x_{32}=3-2=1$, $x_{24}=5-0=5$, $x_{34}=1-0=1$. Check:
node 2 receives $4+1$ and sends 5; node 3 receives 2 and sends $1+1$; the
sink receives $5+1=6$.

**Minimum cut.** The last search returned $S=\{1\}$. The cut
$[\{1\},\{2,3,4\}]$ consists of $(1,2)$ and $(1,3)$, capacity $4+2=6=v^*$, so
it is a minimum cut. It is not the only one: $[\{1,2,3\},\{4\}]$, made of
$(2,4)$ and $(3,4)$, also has capacity $5+1=6$. The other two cuts,
$S=\{1,2\}$ with $u_{13}+u_{24}=7$ and $S=\{1,3\}$ with
$u_{12}+u_{32}+u_{34}=8$, are larger ({numref}`fig-04-03-min-cut`).

```{figure} figures/04-03-ff-iterations.svg
:name: fig-04-03-ff-iterations
:width: 600px
:class: dark-light
:alt: Four residual networks of the diamond example with the augmenting path of each iteration highlighted

Residual networks of {prf:ref}`ex-04-03-diamond`. Each panel shows $G(\vx^{(k)})$
with the augmenting path found in it highlighted: $1\to2\to4$ ($\Delta=4$),
$1\to3\to2\to4$ ($\Delta=1$) and, after backtracking from the dead end 2,
$1\to3\to4$ ($\Delta=1$). In $G(\vx^{(3)})$ no arc leaves node 1.
```

```{figure} figures/04-03-min-cut.svg
:name: fig-04-03-min-cut
:width: 420px
:class: dark-light
:alt: Diamond network with the maximum flow as flow over capacity labels and two dashed cut lines each of capacity 6

Maximum flow of {prf:ref}`ex-04-03-diamond` with labels $x_{ij}/u_{ij}$, and
its two minimum cuts, $S=\{1\}$ (left dashed line) and $S=\{1,2,3\}$ (right
dashed line), each of capacity 6. Every arc crossing a minimum cut forwards is
saturated.
```
:::

(sec-04-03-cuts)=
## Cuts and weak duality

:::{prf:definition} $s$–$t$ cut
:label: def-04-03-cut

An **$s$–$t$ cut** is a partition of the nodes into $S$ and
$\bar S=\mathcal N\setminus S$ with $s\in S$ and $t\in\bar S$. Its
**forward arcs** go from $S$ to $\bar S$, its **backward arcs** from $\bar S$
to $S$. The **capacity** of the cut is the total capacity of the forward
arcs,

$$
u[S,\bar S]=\sum_{(i,j)\in\mathcal A:\ i\in S,\ j\in\bar S}u_{ij}.
$$ (eq-04-03-cut-capacity)
:::

Removing the forward arcs of a cut leaves no directed path from $s$ to $t$,
so every unit of flow must cross a forward arc. The next lemma makes this
precise.

:::{prf:lemma} Flow across a cut
:label: lem-04-03-flow-across-cut

For every feasible flow $\vx$ of value $v$ and every $s$–$t$ cut
$[S,\bar S]$,

$$
v=\sum_{i\in S,\ j\in\bar S}x_{ij}-\sum_{i\in\bar S,\ j\in S}x_{ij}.
$$ (eq-04-03-flow-across-cut)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Add the conservation equations {eq}`eq-04-03-maxflow-lp` of all nodes in $S$.
The right-hand sides add to $v$ (only $s$ contributes, because $t\notin S$).
On the left, an arc with both ends in $S$ appears once with $+$ (at its tail)
and once with $-$ (at its head) and cancels; an arc from $S$ to $\bar S$
appears once with $+$; an arc from $\bar S$ to $S$ once with $-$; arcs with
both ends in $\bar S$ do not appear.
:::

:::{prf:theorem} Weak duality for flows and cuts
:label: thm-04-03-weak-duality

For every feasible flow of value $v$ and every $s$–$t$ cut,
$v\le u[S,\bar S]$. Consequently, if $v=u[S,\bar S]$ for some flow and some
cut, the flow is maximum and the cut is minimum.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

In {eq}`eq-04-03-flow-across-cut` the first sum is at most $u[S,\bar S]$ because
$x_{ij}\le u_{ij}$, and the second sum is nonnegative because $x_{ij}\ge0$.
If equality holds, no flow can be larger than this cut capacity and no cut
can be smaller than this flow value.
:::

The upper bound 65 in {prf:ref}`ex-04-03-pipeline` is the capacity of the cut
$S=\{0,1,2,3\}$, and the bound 110 is that of the cut
$S=\mathcal N\setminus\{7,8,9\}$ (with $t=9$ alone outside, the arcs into 9
have infinite capacity and give no bound at all).

(sec-04-03-maxflow-mincut)=
## The max-flow min-cut theorem

:::{prf:theorem} Max-flow min-cut theorem
:label: thm-04-03-maxflow-mincut

For a feasible flow $\vx$ the following are equivalent.

1. $\vx$ is a maximum flow.
2. The residual network $G(\vx)$ has no augmenting path.
3. There is an $s$–$t$ cut with $u[S,\bar S]=v$.

Hence, when a maximum flow exists (for instance when all capacities are
finite), the maximum flow value equals the minimum capacity of an $s$–$t$
cut.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

(1 ⇒ 2) If there were an augmenting path, augmenting along it would increase
the value by $\Delta>0$, so $\vx$ would not be maximum.

(2 ⇒ 3) Let $S$ be the set of nodes reachable from $s$ by directed paths in
$G(\vx)$. Then $s\in S$, and $t\notin S$ because there is no augmenting path,
so $[S,\bar S]$ is an $s$–$t$ cut. Take $i\in S$ and $j\in\bar S$. Since $j$
is not reachable, $r_{ij}=u_{ij}-x_{ij}+x_{ji}=0$. Because $x_{ij}\le u_{ij}$
and $x_{ji}\ge0$, both terms must vanish: $x_{ij}=u_{ij}$ and $x_{ji}=0$.
Every forward arc of the cut is saturated and every backward arc is empty,
and {eq}`eq-04-03-flow-across-cut` gives $v=u[S,\bar S]$.

(3 ⇒ 1) This is {prf:ref}`thm-04-03-weak-duality`.

For the last statement, apply (1 ⇒ 3) to a maximum flow; by weak duality no
cut has smaller capacity.
:::

The proof of (2 ⇒ 3) is constructive: the set $S$ returned by the final,
unsuccessful labeling is the reachable set, so **Ford–Fulkerson delivers a
minimum cut for free**. In {prf:ref}`ex-04-03-diamond` it gave $S=\{1\}$. A
minimum cut is the planner's bottleneck: raising the capacity of an arc
that is not a forward arc of any minimum cut cannot increase the maximum
flow, and raising a forward arc of the unique minimum cut by one unit can
increase it by at most one unit.

(sec-04-03-integrality)=
## Integrality and termination

:::{prf:theorem} Integrality of maximum flows
:label: thm-04-03-integrality

If all capacities are integers, the Ford–Fulkerson algorithm stops after at
most $v^*$ augmentations, and it returns a maximum flow in which every
$x_{ij}$ is an integer.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Initially all residual capacities are integers. If they are integers before
an augmentation, then $\Delta$, a minimum of them, is a positive integer and
all residual capacities stay integers after it. So every augmentation raises
$v$ by at least 1. Since $v$ is bounded by the capacity of any cut, for
instance $\sum_j u_{sj}$, the number of augmentations is finite and at most
$v^*$. The recovered flows {eq}`eq-04-03-recover` are differences of integers.
By {prf:ref}`thm-04-03-maxflow-mincut` the final flow is maximum.
:::

The same argument works for rational capacities (multiply by a common
denominator). With irrational capacities, a badly chosen sequence of
augmenting paths may never stop and may even converge to a value below the
maximum; this is one reason for the path rule below. The integrality
property is not only a technicality. It says that the LP
{eq}`eq-04-03-maxflow-lp` with integer capacities always has an integer optimal
solution, so problems that need whole units (people, vehicles, matches) can
be solved as maximum flows without integer programming;
{prf:ref}`thm-04-04-tu` traces this to total unimodularity of the incidence
matrix.

(sec-04-03-edmonds-karp)=
### Running time and the Edmonds–Karp rule

A plain search costs $O(\abs{\mathcal A})$, so with integer capacities the
algorithm needs $O(\abs{\mathcal A}\,v^*)$ time. This bound depends on the
size of the numbers, and it can be attained. Take nodes $s,a,b,t$ with
$u_{sa}=u_{sb}=u_{at}=u_{bt}=100$ and $u_{ab}=1$. If the augmenting paths
alternate between $s\to a\to b\to t$ and $s\to b\to a\to t$ (the second uses
the reverse residual arc $(b,a)$ created by the first), each augmentation
sends a single unit, and 200 augmentations are needed, although the two paths
$s\to a\to t$ and $s\to b\to t$ reach $v^*=200$ in two steps. With capacities
$10^6$ the bad choice would take two million steps.

Edmonds and Karp (1972) removed this dependence: **always augment along a
shortest augmenting path, measured in number of arcs**, found by breadth-first
search. They proved that the length of the shortest augmenting path never
decreases and that it can stay the same for at most $\abs{\mathcal A}$
augmentations, which gives at most $O(\abs{\mathcal N}\abs{\mathcal A})$
augmentations and $O(\abs{\mathcal N}\abs{\mathcal A}^2)$ time, whatever the
capacities (even irrational ones). The largest-residual-first labeling of
{prf:ref}`alg-04-03-labelling` is a heuristic in the same spirit: on the bad
example it follows $s\to a\to t$ first (residual 100 rather than 1) and
finishes in two augmentations. Faster methods (capacity scaling,
preflow–push) are described in Ahuja et al. (1993, Chapter 7).

(sec-04-03-recover)=
## Recovering arc flows from the residual network

{prf:ref}`alg-04-03-ford-fulkerson` updates residual capacities only. It is
tempting to obtain the arc flows by adding up the amounts sent along each
path, but this "gross" bookkeeping can be wrong: a path may use a reverse
residual arc $(j,i)$, which *cancels* flow on $(i,j)$ rather than putting
flow on an arc $(j,i)$ that may not even exist.
{prf:ref}`ex-08-05-recover-flow` shows a case where the gross sums put 5 units on
such an arc. The correct flows follow from the definition
{eq}`eq-04-03-residual`, which can be read backwards:

$$
x_{ij}-x_{ji}=u_{ij}-r_{ij}.
$$ (eq-04-03-net-flow)

The residual capacities therefore determine the **net** flow between every
pair of nodes, and only the net flows matter for conservation and for the
value. Choose, for every pair $\{i,j\}$,

$$
\begin{aligned}
x_{ij}&=\max\{0,\ u_{ij}-r_{ij}\},\\
x_{ji}&=\max\{0,\ u_{ji}-r_{ji}\}.
\end{aligned}
$$ (eq-04-03-recover)

At most one of the two is positive. Indeed, by {eq}`eq-04-03-residual-sum`,

$$
(u_{ij}-r_{ij})+(u_{ji}-r_{ji})=0,
$$

so the two differences are negatives of each other, and both cannot be
positive. Three cases arise:

- **only $(i,j)$ exists** ($u_{ji}=0$): $u_{ji}-r_{ji}=-r_{ji}\le0$, so
  $x_{ji}=0$ and $x_{ij}=u_{ij}-r_{ij}$;
- **only $(j,i)$ exists** ($u_{ij}=0$): symmetric, $x_{ij}=0$ and
  $x_{ji}=u_{ji}-r_{ji}$;
- **both exist**: the positive difference (if any) gives the flow in its
  direction and the other direction gets 0.

The flows {eq}`eq-04-03-recover` are feasible: $0\le x_{ij}\le u_{ij}$ because
$r_{ij}\ge0$, the net flows equal those of a feasible flow, and therefore
conservation and the value are those of the final flow. In
{prf:ref}`ex-04-03-diamond` this rule gave the flows of {numref}`fig-04-03-min-cut`
directly from the residual capacities of $G(\vx^{(3)})$.

(sec-04-03-lp-dual)=
## The LP dual is the minimum-cut problem

The max-flow min-cut theorem is an instance of LP duality
({prf:ref}`thm-03-05-strong-duality`). Attach a multiplier $\pi_i$ (free) to the
conservation equation of node $i$ in {eq}`eq-04-03-maxflow-lp`, written as
$\sum_j x_{ji}-\sum_j x_{ij}+v[i=s]-v[i=t]=0$, and a multiplier $y_{ij}\ge0$
to each capacity constraint $x_{ij}\le u_{ij}$. The dual LP is

$$
\begin{aligned}
\min_{\vpi,\,\vy}\ & \sum_{(i,j)\in\mathcal A}u_{ij}y_{ij} \\
\st\ & y_{ij}\ge\pi_i-\pi_j, \quad (i,j)\in\mathcal A,\\
& \pi_s-\pi_t=1,\\
& y_{ij}\ge0, \quad (i,j)\in\mathcal A.
\end{aligned}
$$ (eq-04-03-dual-lp)

Every cut gives a dual solution: $\pi_i=1$ for $i\in S$, $\pi_i=0$ for
$i\in\bar S$, and $y_{ij}=1$ exactly on the forward arcs. Its objective value
is $u[S,\bar S]$. Conversely, one can show that the dual has an optimal
solution of this 0–1 form, so the dual LP *is* the minimum-cut problem, and
LP strong duality is the max-flow min-cut theorem. The flow–cut pair found by
Ford–Fulkerson is a primal–dual optimal pair, and complementary slackness
reads exactly as in the proof of {prf:ref}`thm-04-03-maxflow-mincut`: $y_{ij}>0$
(a forward arc) forces $x_{ij}=u_{ij}$, and an arc with
$\pi_i-\pi_j<0$ (a backward arc) has a positive reduced cost and carries no
flow.

(sec-04-03-matching)=
## Bipartite matching

Let $W$ be a set of workers, $J$ a set of jobs, and $E\subseteq W\times J$
the pairs (worker, job) that are acceptable. A **matching** is a set of pairs
of $E$ in which no worker and no job appears twice; the task is to find a
matching with as many pairs as possible. Build the network

- with arcs $(s,w)$ of capacity 1 for every worker $w$,
- arcs $(w,j)$ of capacity 1 for every acceptable pair $(w,j)\in E$, and
- arcs $(j,t)$ of capacity 1 for every job $j$.

An integer flow sends 0 or 1 unit through each arc. Because each worker
receives at most one unit from $s$ and each job passes at most one unit to
$t$, the arcs $(w,j)$ with $x_{wj}=1$ form a matching with $v$ pairs; and
every matching gives an integer flow of the same value. By
{prf:ref}`thm-04-03-integrality` Ford–Fulkerson returns an integer maximum flow,
hence a **maximum matching**, in at most $\min\{\abs{W},\abs{J}\}$
augmentations. An augmenting path here alternates between new pairs and
canceled ones: it may reassign a worker to another job to make room for a
new worker.

The minimum cut has a combinatorial meaning as well. Give the middle arcs
$(w,j)$ infinite capacity instead of 1; this changes no integer flow, because
a worker receives at most one unit anyway. A cut of finite capacity then has
no forward arc $(w,j)$, so every acceptable pair has its worker outside $S$
or its job inside $S$, and the capacity of the cut is the number of these
workers (arcs $(s,w)$ with $w\in\bar S$) plus the number of these jobs (arcs
$(j,t)$ with $j\in S$). Hence the size of a maximum matching equals the
smallest number of workers and jobs that together touch every acceptable
pair (König's theorem). {prf:ref}`ex-08-05-chores` works through
an assignment of chores to children in this way.

## Summary and Common Mistakes

### Summary

- The maximum-flow problem is the LP {eq}`eq-04-03-maxflow-lp`; undirected
  links, several sources or sinks, node capacities and unlimited arcs are
  brought into this form by the transformations of
  {ref}`sec-04-03-modelling` ({prf:ref}`ex-04-03-pipeline`, $v^*=65$).
- The residual network $G(\vx)$ ({prf:ref}`def-04-03-residual`) has capacities
  $r_{ij}=u_{ij}-x_{ij}+x_{ji}$ {eq}`eq-04-03-residual`; a directed $s$–$t$ path
  in it increases the flow by its smallest residual capacity.
- Ford–Fulkerson ({prf:ref}`alg-04-03-ford-fulkerson`) augments along paths
  found by the labeling procedure ({prf:ref}`alg-04-03-labelling`) until none is
  left; on the diamond network it sends 4, 1 and 1 units, $v^*=6$
  ({prf:ref}`ex-04-03-diamond`, {numref}`fig-04-03-ff-iterations`).
- Every flow is bounded by every cut capacity
  ({prf:ref}`thm-04-03-weak-duality`), and the maximum flow equals the minimum cut
  capacity ({prf:ref}`thm-04-03-maxflow-mincut`); the nodes reached by the last
  labeling form a minimum cut ({numref}`fig-04-03-min-cut`). The LP dual
  {eq}`eq-04-03-dual-lp` is the minimum-cut problem.
- Integer capacities give integer maximum flows
  ({prf:ref}`thm-04-03-integrality`). Augmenting along shortest paths
  (Edmonds–Karp) bounds the work by $O(\abs{\mathcal N}\abs{\mathcal A}^2)$.
- Arc flows are recovered from the final residual capacities by
  {eq}`eq-04-03-recover`, and bipartite matching is a maximum flow with unit
  capacities ({ref}`sec-04-03-matching`).

### Common mistakes

- **Using arc capacities instead of residual capacities** after the first
  iteration. The search runs on $G(\vx)$, including the reverse arcs; an
  optimal solution often needs them (iteration 2 of
  {prf:ref}`ex-04-03-diamond` did not, but {prf:ref}`ex-08-05-maxflow-eight` does).
- **Reading the flows from the gross path sums.** A path through a reverse
  residual arc cancels flow; recover the flows with {eq}`eq-04-03-recover`.
- **Counting backward arcs in a cut capacity.** Only arcs from $S$ to
  $\bar S$ count in {eq}`eq-04-03-cut-capacity`.
- **Forgetting both directions of an undirected link**, or giving the two
  arcs half the capacity each. Each direction gets the full capacity.
- **Stopping too early.** An empty candidate set at a node other than $s$ is
  a dead end, not the end: backtrack. Only an empty set at $s$ proves
  optimality.
- **Assuming the minimum cut is unique.** Different cuts can have the same
  minimum capacity, as in {prf:ref}`ex-04-03-diamond`; the final labeling gives
  the one closest to $s$.

## References

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows:
  Theory, algorithms, and applications*. Prentice Hall. (Chapter 6: flows and
  cuts, residual networks, generic augmenting-path and labeling algorithms,
  max-flow min-cut theorem, combinatorial implications; Chapter 7: shortest
  augmenting path and capacity scaling; Chapter 12: bipartite matching.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Section 7.5: the maximum flow problem.)
- Edmonds, J., & Karp, R. M. (1972). Theoretical improvements in algorithmic
  efficiency for network flow problems. *Journal of the ACM, 19*(2), 248–264.
  <https://doi.org/10.1145/321694.321699>
- Ford, L. R., & Fulkerson, D. R. (1956). Maximal flow through a network.
  *Canadian Journal of Mathematics, 8*, 399–404.
  <https://doi.org/10.4153/CJM-1956-045-5>
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Section 10.5: the maximum flow
  problem.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Section 6.4: maximal flow model, cuts and the maximal flow
  algorithm.)
