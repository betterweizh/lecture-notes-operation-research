# Minimum-Cost Flow

## Overview

The maximum-flow problem of {doc}`04-03-maximum-flow` asks *how much* can be
sent; it ignores what sending costs. Most logistics questions are the other
way round: given what each plant can supply and what each market needs, ship
the goods through a network of roads, depots and pipelines **at least
cost**, respecting arc capacities and, sometimes, minimum flows required by
contracts. This is the **minimum-cost flow problem**. It contains the
shortest-path, maximum-flow, transportation, assignment and transshipment
problems as special cases, and it is still a linear program with a very
special matrix: the node–arc incidence matrix. That matrix is **totally
unimodular**, so integer supplies and capacities give integer optimal flows.
The chapter states the model, shows how to handle unbalanced data and lower
bounds, proves the integrality property, derives the optimality conditions
in two equivalent forms (reduced costs with node potentials, and absence of
negative-cost cycles), and turns them into algorithms: successive shortest
paths, cycle canceling and an outline of the network simplex method. It
uses residual networks and augmenting paths from {doc}`04-03-maximum-flow`,
shortest paths with labels from {doc}`04-02-shortest-paths`, basic feasible
solutions from {doc}`03-02-lp-geometry` and the dual LP and complementary
slackness from {doc}`03-05-lp-duality`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- formulate a minimum-cost flow problem with supplies, demands, unit costs
  and lower and upper bounds as a linear program;
- detect infeasible or unbalanced data, add dummy nodes, and remove lower
  bounds by a change of variables;
- prove that the node–arc incidence matrix is totally unimodular and explain
  why integer data give integer optimal flows;
- write the shortest-path, maximum-flow, transportation, assignment and
  transshipment problems as minimum-cost flow problems;
- check the optimality of a flow with node potentials and reduced costs, or
  by searching the residual network for negative cycles; and
- solve small instances by successive shortest paths or cycle canceling and
  describe the network simplex method.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $b_i$ | net supply of node $i$: $b_i>0$ supply, $b_i<0$ demand, $b_i=0$ transshipment node (units) |
| $c_{ij}$ | cost per unit of flow on arc $(i,j)$ (\$ per unit) |
| $l_{ij}$, $u_{ij}$ | lower and upper bound on the flow $x_{ij}$ (units) |
| $z$ | total cost $\sum c_{ij}x_{ij}$ (\$) |
| $\pi_i$ | potential (dual variable) of node $i$ (\$ per unit) |
| $y_{ij}$ | dual variable of the upper bound $x_{ij}\le u_{ij}$ |
| $\bar c_{ij}=c_{ij}-\pi_i+\pi_j$ | reduced cost of arc $(i,j)$ |
| $G(\vx)$ | residual network of $\vx$, with cost $c_{ij}$ on forward and $-c_{ij}$ on backward residual arcs |
| $\mathcal W$, $c(\mathcal W)$ | a directed cycle in $G(\vx)$ and its cost |

(sec-04-04-model)=
## The minimum-cost flow problem

Let $G=(\mathcal N,\mathcal A)$ be a directed network. Each node $i$ has a
**net supply** $b_i$: it produces $b_i$ units if $b_i>0$, consumes $-b_i$
units if $b_i<0$, and only passes flow on if $b_i=0$. Each arc has a unit cost
$c_{ij}$ and bounds $l_{ij}\le x_{ij}\le u_{ij}$ (with $l_{ij}=0$ and
$u_{ij}=\infty$ when no bound is given).

:::{prf:definition} Minimum-cost flow problem
:label: def-04-04-mcf

The **minimum-cost flow problem** is the linear program

$$
\begin{aligned}
\min_{\vx}\ & z=\sum_{(i,j)\in\mathcal A}c_{ij}x_{ij} \\
\st\ & \sum_{j:(i,j)\in\mathcal A}x_{ij}-\sum_{j:(j,i)\in\mathcal A}x_{ji}=b_i, \quad i\in\mathcal N,\\
& l_{ij}\le x_{ij}\le u_{ij}, \quad (i,j)\in\mathcal A.
\end{aligned}
$$ (eq-04-04-mcf)

The equality constraints are **flow-balance** (conservation) equations: what
leaves node $i$ minus what enters it equals its net supply. In matrix form
they read $\mA\vx=\vb$ with the node–arc incidence matrix $\mA$ of
{ref}`sec-04-01-representations`, whose column for arc $(i,j)$ is $\ve_i-\ve_j$.
:::

:::{prf:example} A supply chain with bounded arcs
:label: ex-04-04-supply-chain

Three plants (nodes 1, 2, 3) supply 100, 200 and 50 units of a product;
three markets (nodes 4, 5, 6) need 150, 80 and 120 units. Node 3 is both a
plant and a depot through which other plants may ship. The arcs, their unit
costs in dollars and the bounds required by contracts and by the capacity of
the links are:

| Arc | $(1,2)$ | $(1,3)$ | $(1,4)$ | $(2,3)$ | $(2,5)$ |
|---|---|---|---|---|---|
| $c_{ij}$ (\$) | 3 | 4 | 1 | 5 | 6 |
| $[l_{ij},u_{ij}]$ | – | – | $[50,80]$ | – | – |

| Arc | $(3,4)$ | $(3,5)$ | $(4,6)$ | $(5,6)$ |
|---|---|---|---|---|
| $c_{ij}$ (\$) | 1 | 2 | 2 | 4 |
| $[l_{ij},u_{ij}]$ | $[70,120]$ | – | $[100,120]$ | – |

(a dash means $[0,\infty)$). Find the shipping plan of least cost.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With $\vb=(100,200,50,-150,-80,-120)$ the model {eq}`eq-04-04-mcf` is

$$
\begin{aligned}
\min\ & x_{14}+4x_{13}+3x_{12}+5x_{23}+6x_{25}+x_{34}+2x_{35}+2x_{46}+4x_{56}\\
\st\ & x_{12}+x_{13}+x_{14}=100\\
& x_{23}+x_{25}-x_{12}=200\\
& x_{34}+x_{35}-x_{13}-x_{23}=50\\
& x_{46}-x_{14}-x_{34}=-150\\
& x_{56}-x_{25}-x_{35}=-80\\
& -x_{46}-x_{56}=-120\\
& 50\le x_{14}\le80\\
& 70\le x_{34}\le120\\
& 100\le x_{46}\le120\\
& \text{all other } x_{ij}\ge0.
\end{aligned}
$$

*(The balance of node 2 subtracts the inflow $x_{12}$; a version of this
model with $x_{13}$ in that row does not describe the network.)*

These data turn out to be infeasible ({ref}`sec-04-04-feasibility`). With the
lower bound on $(4,6)$ dropped, the optimal plan of
{numref}`tab-04-04-supply-chain-solution` ships $x_{13}=20$, $x_{14}=80$,
$x_{23}=50$, $x_{25}=150$, $x_{34}=120$, $x_{46}=50$, $x_{56}=70$ and nothing
on $(1,2)$ and $(3,5)$, at a cost of \$1 810
({ref}`sec-04-04-supply-chain-solved`).
:::

Before the model of {prf:ref}`ex-04-04-supply-chain` is solved, the next
section checks whether a solution exists at all.

(sec-04-04-feasibility)=
## Feasibility, balance and dummy nodes

**Balance.** Adding all balance equations gives $0=\sum_i b_i$, because every
$x_{ij}$ appears once with $+$ and once with $-$. A problem with
$\sum_ib_i\ne0$ is therefore infeasible as stated, and it must be decided
what the imbalance means:

- **Excess supply** ($\sum_i b_i>0$): not everything produced can be sold.
  Add a **dummy demand node** $d$ with $b_d=-\sum_ib_i$ and an arc of cost 0
  (or the cost of storage or disposal) from every supply node to $d$. Flow
  into $d$ is the unused capacity of each plant.
- **Excess demand** ($\sum_ib_i<0$): some demand cannot be met. Add a
  **dummy supply node** with $b=-\sum_ib_i$ and arcs to every demand node with
  a cost equal to the penalty for each unit short (lost profit). Flow out of
  it is the shortage of each market.

In both cases the equations become consistent and $\mA$ has rank
$\abs{\mathcal N}-1$ on a connected network, so one balance equation is
redundant.

**Capacity.** Balanced data can still be infeasible because of the bounds.
For a set of nodes $S$, the flow leaving $S$ minus the flow entering it must
equal the net supply $\sum_{i\in S}b_i$; this is possible only if

$$
\sum_{i\in S}b_i\le\sum_{(i,j):\,i\in S,\,j\notin S}u_{ij}-\sum_{(j,i):\,j\notin S,\,i\in S}l_{ji}
\quad\text{for every } S\subseteq\mathcal N,
$$ (eq-04-04-gale)

that is, the net supply of $S$ must fit through the arcs leaving it after the
compulsory flows on the arcs entering it. Gale (1957) and Hoffman (1960)
showed that these conditions (together with $\vl\le\vu$ and balance) are
also sufficient. In
practice feasibility is checked with one maximum-flow computation: after
removing the lower bounds ({ref}`sec-04-04-lower-bounds`), connect a
super-source to every supply node $i$ by an arc of capacity $b_i$ and every
demand node $j$ to a super-sink by an arc of capacity $-b_j$; the problem is
feasible exactly when the maximum flow saturates all source arcs.

**The supply-chain example is infeasible.** Take $S=\mathcal N\setminus\{4\}$
in {eq}`eq-04-04-gale`: its net supply is $-b_4=150$, the arcs leaving it,
$(1,4)$ and $(3,4)$, can carry $80+120=200$, and the arc $(4,6)$ entering it
must carry at least 100, so the condition $150\le200-100$ fails. In words:
node 4 must deliver 150 units to its market and at least $l_{46}=100$ units
to node 6, so it needs an inflow of at least 250 units, but its two incoming
arcs can carry at most 200. An LP
solver confirms that {eq}`eq-04-04-mcf` has no feasible solution for these
data. Some bound must give way. The smallest change is in the lower bound of
$(4,6)$: any $l_{46}\le50$ makes the problem feasible, and the rest of this
chapter uses the repaired bounds $0\le x_{46}\le120$ (the optimum is the same
for every $l_{46}\le50$).

(sec-04-04-lower-bounds)=
## Removing lower bounds

Lower bounds are removed by measuring each flow above its minimum. Substitute

$$
\begin{aligned}
x_{ij}&=l_{ij}+x'_{ij},\\
0\le x'_{ij}&\le u_{ij}-l_{ij}.
\end{aligned}
$$ (eq-04-04-shift)

The compulsory $l_{ij}$ units leave node $i$ and arrive at node $j$ in any
case, so the net supplies change to

$$
b'_i=b_i-\sum_{j:(i,j)\in\mathcal A}l_{ij}+\sum_{j:(j,i)\in\mathcal A}l_{ji},
$$ (eq-04-04-shifted-supply)

and the cost becomes $z=\sum c_{ij}x'_{ij}+\sum c_{ij}l_{ij}$, a constant plus
a problem with zero lower bounds. Every algorithm below may therefore assume
$l_{ij}=0$.

For the repaired supply chain the bounded arcs are $(1,4)$ with $[50,80]$
and $(3,4)$ with $[70,120]$. Shifting them gives

- $b'_1=100-50=50$, $b'_3=50-70=-20$, $b'_4=-150+50+70=-30$, the other
  $b'_i$ unchanged;
- $0\le x'_{14}\le30$ and $0\le x'_{34}\le50$;
- the constant $1\cdot50+1\cdot70=\$120$.

Node 3 has turned from a supply node into a demand node: its own 50 units
and 20 more are already committed to the arc $(3,4)$.

(sec-04-04-integrality)=
## Total unimodularity and integral flows

:::{prf:definition} Totally unimodular matrix
:label: def-04-04-tu

A matrix is **totally unimodular** (TU) if the determinant of every square
submatrix is $0$, $+1$ or $-1$. In particular every entry is $0$, $\pm1$.
:::

:::{prf:theorem} The incidence matrix is totally unimodular
:label: thm-04-04-tu

The node–arc incidence matrix $\mA$ of a directed network is totally
unimodular. So is the matrix $\begin{pmatrix}\mA&\vzero\\ \mI&\mI\end{pmatrix}$
of the constraints $\mA\vx=\vb$, $\vx+\vw=\vu$ (upper bounds with slack
variables $\vw$).
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\mQ$ be a $k\times k$ submatrix of $\mA$; we use induction on $k$. For
$k=1$ the entry is $0$ or $\pm1$. For $k>1$ every column of $\mQ$ has at most
one $+1$ and at most one $-1$, because each column of $\mA$ does. Three cases
remain. If $\mQ$ has a zero column, $\det\mQ=0$. If some column has a single
nonzero entry, expand $\det\mQ$ along it: $\det\mQ=\pm\det\mQ'$ with a
$(k-1)\times(k-1)$ submatrix $\mQ'$, which is $0$ or $\pm1$ by induction.
Otherwise every column has exactly one $+1$ and one $-1$, the rows of $\mQ$
add up to the zero vector, and $\det\mQ=0$.

For the second matrix take a square submatrix. A column belonging to $\vw$
has a single 1 in the full matrix, so in the submatrix it is either zero
(determinant 0) or a unit column, along which we expand. Repeating this
removes all $\vw$-columns and leaves a square submatrix of
$\begin{pmatrix}\mA\\\mI\end{pmatrix}$; its rows from $\mI$ are zero or unit
rows, and expanding along them leaves a square submatrix of $\mA$.
:::

:::{prf:theorem} Integrality of minimum-cost flows
:label: thm-04-04-integrality

If all $b_i$, $l_{ij}$ and $u_{ij}$ are integers, every basic feasible
solution of {eq}`eq-04-04-mcf` is integral. Hence, if the problem has an optimal
solution, it has an integral one, and the simplex method finds it.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

After the shift {eq}`eq-04-04-shift`, which keeps the data integral, write the
upper bounds with slacks as in {prf:ref}`thm-04-04-tu` and drop one redundant
balance row. A basic feasible solution ({prf:ref}`def-03-02-basic-solution`) solves
$\mB\vx_B=\vb'$ with a nonsingular square submatrix $\mB$ of a TU matrix, so
$\det\mB=\pm1$. By Cramer's rule each component of $\vx_B$ is a ratio of an
integer determinant to $\pm1$, an integer. The nonbasic variables are 0.
An LP in standard form with an optimum has an optimal basic feasible
solution ({prf:ref}`thm-03-02-fundamental`).
:::

This is the reason network problems are "combinatorial although they are
LPs": the assignment problem below needs 0–1 decisions, the matching problem
of {ref}`sec-04-03-matching` needs whole pairs, and both are solved as plain
LPs. The costs $c_{ij}$ need not be integers.

(sec-04-04-special-cases)=
## Special cases

**Shortest path.** Send one unit from $s$ to $t$: $b_s=1$, $b_t=-1$, $b_i=0$
otherwise, $u_{ij}=\infty$, $c_{ij}$ the arc lengths. An integral optimal
flow is a 0–1 vector whose arcs form a shortest $s$–$t$ path (provided there
is no negative cycle). With $b_s=\abs{\mathcal N}-1$ and $b_i=-1$ for every
other node the solution is a shortest-path tree ({prf:ref}`prop-04-02-tree`).

**Maximum flow.** Set all costs to 0 and all $b_i=0$, and add an arc $(t,s)$
of infinite capacity and cost $-1$. A flow is now a circulation, and its cost
is minus the flow returned from $t$ to $s$, that is, minus the value of the
$s$–$t$ flow; the minimum cost is $-v^*$.

**Transportation.** The nodes are supply points (plants, $b_i>0$) and demand
points (markets, $b_j<0$), and every arc goes from a plant to a market. The
balance equations say that each plant ships its supply and each market
receives its demand:

$$
\begin{aligned}
\sum_j x_{ij}&=b_i &&\text{(plants)},\\
\sum_i x_{ij}&=-b_j &&\text{(markets)}.
\end{aligned}
$$ (eq-04-04-transportation)

**Assignment.** A transportation problem with $n$ plants (workers) and $n$
markets (jobs), every $b_i=1$ and every $b_j=-1$. By
{prf:ref}`thm-04-04-integrality` an optimal basic solution has $x_{ij}\in\{0,1\}$,
one 1 in every row and every column: a permutation. Specialized algorithms
(the Hungarian method of Kuhn, 1955) exist, but the LP is already exact.

**Transshipment.** Intermediate nodes with $b_i=0$ (depots, ports,
cross-docks) through which goods may pass. Any uncapacitated minimum-cost flow
problem is a transshipment problem; {prf:ref}`ex-08-05-transshipment` adds a
depot capacity and a dummy market.

:::{prf:example} A transportation problem
:label: ex-04-04-transportation

Two plants A and B produce 40 and 50 tonnes of a product per week; three
markets need 20, 25 and 45 tonnes. Shipping costs in dollars per tonne are

| | Market 1 | Market 2 | Market 3 | Supply |
|---|---|---|---|---|
| Plant A | 4 | 6 | 9 | 40 |
| Plant B | 2 | 7 | 5 | 50 |
| Demand | 20 | 25 | 45 | 90 |

Find the cheapest shipping plan.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The data are balanced ($40+50=20+25+45=90$); the network is the bipartite
network of {numref}`fig-04-04-transportation` with $b_A=40$, $b_B=50$,
$b_1=-20$, $b_2=-25$, $b_3=-45$ and six uncapacitated arcs. Successive
shortest paths ({ref}`sec-04-04-ssp`, {numref}`tab-04-04-ssp`) find the optimal plan
$x_{A1}=15$, $x_{A2}=25$, $x_{B1}=5$, $x_{B3}=45$, at \$445 per week.
:::

:::{prf:example} An assignment problem
:label: ex-04-04-assignment

Three workers W1, W2, W3 are to be given one job each among J1, J2, J3. The
time each worker needs for each job, in hours, is

| | J1 | J2 | J3 |
|---|---|---|---|
| W1 | 9 | 4 | 7 |
| W2 | 3 | 6 | 5 |
| W3 | 8 | 7 | 4 |

Minimize the total time.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

As a minimum-cost flow problem each worker is a node with $b=1$, each job a
node with $b=-1$, and the nine arcs (worker, job) cost the times above. The
LP relaxation, solved without any integrality constraint, returns
$x_{W1,J2}=x_{W2,J1}=x_{W3,J3}=1$ and all other $x=0$: W1 does J2, W2 does
J1 and W3 does J3, in $4+3+4=11$ hours. Listing the $3!=6$ permutations
confirms it (their totals are 11, 17, 17, 19, 21 and 21 hours). The LP did not
need to be told that the answer must be a permutation; total unimodularity
guarantees it. The potentials of the next section certify the result: with
$\pi_{W1}=1$, $\pi_{W2}=3$, $\pi_{W3}=4$, $\pi_{J1}=0$, $\pi_{J2}=-3$,
$\pi_{J3}=0$, every reduced cost $c_{ij}-\pi_i+\pi_j$ is nonnegative (for
example $9-1+0=8$ for W1–J1) and the three chosen arcs have reduced cost 0.
:::

(sec-04-04-optimality)=
## Optimality conditions

### Reduced costs and node potentials

Assume $\vl=\vzero$ (after {eq}`eq-04-04-shift`) and attach a dual variable
$\pi_i$ to the balance equation of node $i$ and $y_{ij}\ge0$ to the upper
bound $x_{ij}\le u_{ij}$. The dual of {eq}`eq-04-04-mcf` ({doc}`03-05-lp-duality`)
is

$$
\begin{aligned}
\max_{\vpi,\,\vy}\ & \sum_i b_i\pi_i-\sum_{(i,j)}u_{ij}y_{ij}\\
\st\ & \pi_i-\pi_j-y_{ij}\le c_{ij}, \quad (i,j)\in\mathcal A,\\
& y_{ij}\ge0,
\end{aligned}
$$ (eq-04-04-dual)

with $\vpi$ free. The $\pi_i$ are called **node potentials**, and

$$
\bar c_{ij}=c_{ij}-\pi_i+\pi_j
$$ (eq-04-04-reduced-cost)

is the **reduced cost** of arc $(i,j)$: the reduced cost
$c_j-\vpi^\T\mA_{\cdot j}$ of the simplex method for the column
$\ve_i-\ve_j$. Adding a constant to all potentials changes no reduced cost
(one balance equation is redundant), so one potential can be fixed, for
example at 0.

:::{prf:theorem} Reduced-cost optimality conditions
:label: thm-04-04-reduced-cost

A feasible flow $\vx$ of {eq}`eq-04-04-mcf` is optimal if and only if there are
potentials $\vpi$ such that for every arc

$$
\begin{aligned}
\bar c_{ij}>0 &\ \Rightarrow\ x_{ij}=l_{ij},\\
\bar c_{ij}<0 &\ \Rightarrow\ x_{ij}=u_{ij},\\
l_{ij}<x_{ij}<u_{ij} &\ \Rightarrow\ \bar c_{ij}=0.
\end{aligned}
$$ (eq-04-04-cs)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

With $\vl=\vzero$, choose $y_{ij}=\max\{0,-\bar c_{ij}\}$, the smallest $y_{ij}$ that
makes $(\vpi,\vy)$ dual feasible. Complementary slackness for the pair
{eq}`eq-04-04-mcf`–{eq}`eq-04-04-dual` requires $x_{ij}>0\Rightarrow
\pi_i-\pi_j-y_{ij}=c_{ij}$ and $y_{ij}>0\Rightarrow x_{ij}=u_{ij}$, which with
this $\vy$ are exactly the three implications {eq}`eq-04-04-cs`. By LP duality
({prf:ref}`thm-03-05-complementary-slackness`) a feasible primal solution is
optimal if and only if a dual solution complementary to it exists. With
$\vl\ne\vzero$ apply the argument to the shifted
problem; reduced costs do not change under the shift.
:::

The conditions have an economic reading. Think of $\pi_i$ as the value of one
unit of the product at node $i$ (high at plants that have it, low at markets
that need it). Shipping from $i$ to $j$ moves a unit worth $\pi_i$ to a place
where it is worth $\pi_j$ at a cost $c_{ij}$. If $\bar c_{ij}>0$ the trip
costs more than it gains and the arc is used as little as allowed; if
$\bar c_{ij}<0$ it gains and the arc is filled; arcs in between are exactly
break-even.

### Negative cycles

The residual network of {doc}`04-03-maximum-flow` extends to costs. For a flow
$\vx$ with $\vl\le\vx\le\vu$, $G(\vx)$ contains

- a forward arc $(i,j)$ with residual capacity $u_{ij}-x_{ij}$ and cost
  $c_{ij}$ when $x_{ij}<u_{ij}$, and
- a backward arc $(j,i)$ with residual capacity $x_{ij}-l_{ij}$ and cost
  $-c_{ij}$ when $x_{ij}>l_{ij}$ (pushing flow back saves $c_{ij}$ per unit).

A directed cycle $\mathcal W$ in $G(\vx)$ can carry
$\Delta=\min$ (residual capacity on $\mathcal W$) $>0$ units without changing
any net supply, and the total cost changes by $\Delta\,c(\mathcal W)$, where
$c(\mathcal W)$ is the sum of the costs around the cycle.

:::{prf:theorem} Negative-cycle optimality conditions
:label: thm-04-04-negative-cycle

A feasible flow $\vx$ is optimal if and only if its residual network $G(\vx)$
contains no directed cycle of negative cost.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If $G(\vx)$ has a negative cycle, sending flow around it gives a feasible
flow of lower cost, so $\vx$ is not optimal. Conversely, suppose there is no
negative cycle. Add a node $r$ with zero-cost arcs to every node; then the
shortest-path lengths $d_i$ from $r$ in $G(\vx)$ are well defined and satisfy
$d_j\le d_i+(\text{cost of }(i,j))$ for every residual arc
({prf:ref}`thm-04-02-optimality`). Put $\pi_i=-d_i$. For a forward residual arc
($x_{ij}<u_{ij}$) this gives $c_{ij}-\pi_i+\pi_j=c_{ij}+d_i-d_j\ge0$, so
$\bar c_{ij}\ge0$; for a backward residual arc ($x_{ij}>l_{ij}$) it gives
$-c_{ij}+d_j-d_i\ge0$, so $\bar c_{ij}\le0$. Hence $\bar c_{ij}>0$ forces
$x_{ij}=l_{ij}$, $\bar c_{ij}<0$ forces $x_{ij}=u_{ij}$, and an arc strictly
between its bounds has $\bar c_{ij}=0$: these are the conditions
{eq}`eq-04-04-cs`, and $\vx$ is optimal by {prf:ref}`thm-04-04-reduced-cost`.
:::

(sec-04-04-supply-chain-solved)=
### The supply chain solved

The repaired supply-chain problem has the optimal flow of
{numref}`tab-04-04-supply-chain-solution`, at a cost of

$$
\begin{aligned}
z^*&=1\cdot80+4\cdot20+3\cdot0+5\cdot50+6\cdot150+1\cdot120+2\cdot0+2\cdot50+4\cdot70\\
&=\$1\,810.
\end{aligned}
$$

The potentials in the table (normalized to $\pi_6=0$) are $\pi=(9,10,5,2,4,0)$
for nodes 1 to 6. Every arc satisfies {eq}`eq-04-04-cs`: the five arcs strictly
inside their bounds have $\bar c_{ij}=0$; $(1,2)$ and $(3,5)$ are at their
lower bound 0 with $\bar c_{ij}>0$; $(1,4)$ and $(3,4)$ are at their upper
bounds with $\bar c_{ij}<0$. So the flow is optimal, and since no arc at a
bound has reduced cost 0, it is the only optimal flow. The flows are
integers, as {prf:ref}`thm-04-04-integrality` predicts. The balance checks are
immediate; for example, node 3 receives $20+50$ and sends $120+0$, a net supply of
$50$.

```{list-table} Optimal flow of the repaired supply chain, with potentials $\pi=(9,10,5,2,4,0)$
:name: tab-04-04-supply-chain-solution
:header-rows: 1
:stub-columns: 1

* - Arc
  - $c_{ij}$ (\$)
  - $[l_{ij},u_{ij}]$
  - $x^*_{ij}$
  - $\bar c_{ij}=c_{ij}-\pi_i+\pi_j$
  - status
* - $(1,2)$
  - 3
  - $[0,\infty)$
  - 0
  - 4
  - at lower bound
* - $(1,3)$
  - 4
  - $[0,\infty)$
  - 20
  - 0
  - basic
* - $(1,4)$
  - 1
  - $[50,80]$
  - 80
  - $-6$
  - at upper bound
* - $(2,3)$
  - 5
  - $[0,\infty)$
  - 50
  - 0
  - basic
* - $(2,5)$
  - 6
  - $[0,\infty)$
  - 150
  - 0
  - basic
* - $(3,4)$
  - 1
  - $[70,120]$
  - 120
  - $-2$
  - at upper bound
* - $(3,5)$
  - 2
  - $[0,\infty)$
  - 0
  - 1
  - at lower bound
* - $(4,6)$
  - 2
  - $[0,120]$
  - 50
  - 0
  - basic
* - $(5,6)$
  - 4
  - $[0,\infty)$
  - 70
  - 0
  - basic
```

```{figure} figures/04-04-supply-chain.svg
:name: fig-04-04-supply-chain
:width: 600px
:class: dark-light
:alt: Supply-chain network with supplies and demands, arc costs and bounds, and on the right the optimal flows with the five basic arcs highlighted

(a) The supply chain of {prf:ref}`ex-04-04-supply-chain`: net supplies in
brackets, arc labels "cost (lower, upper)" for the bounded arcs. The lower
bound 100 on $(4,6)$ makes the data infeasible and is dropped. (b) The
optimal flows $x^*_{ij}$; the highlighted arcs are strictly between their
bounds and form a spanning tree, $(1,4)$ and $(3,4)$ are at their upper
bounds, $(1,2)$ and $(3,5)$ carry nothing.
```

The reduced costs are also prices. The two binding upper bounds have
$\bar c_{14}=-6$ and $\bar c_{34}=-2$: one more unit of capacity on $(1,4)$
would lower the cost by \$6, on $(3,4)$ by \$2. The potential difference
$\pi_1-\pi_6=9$ is the marginal cost of one more unit made at plant 1 and
sold at market 6: it would travel $1\to3$ (\$4), displace a unit on
$2\to3$ (saving \$5), and that unit would go $2\to5\to6$ (\$6 + \$4), in all
$4-5+6+4=\$9$. Re-solving the LP with these changes gives exactly \$1 804,
\$1 808 and \$1 819.

(sec-04-04-ssp)=
## Successive shortest paths

The idea of the first algorithm is to keep the cost optimal at all times and
work towards feasibility. Add a super-source $s$ with an arc $(s,i)$ of
capacity $b_i$ to each supply node and a super-sink $t$ with an arc $(j,t)$ of
capacity $-b_j$ from each demand node, both of cost 0, and assume $\vl=\vzero$.
Starting from the zero flow, which has no negative cycle when all
$c_{ij}\ge0$, repeatedly send flow along a *cheapest* augmenting path.

:::{prf:algorithm} Successive shortest paths
:label: alg-04-04-ssp

**Input:** a balanced problem {eq}`eq-04-04-mcf` with $\vl=\vzero$, $\vc\ge\vzero$, extended
by $s$ and $t$ as above.
**Output:** a minimum-cost flow, or the conclusion that none is feasible.

1. **Initialize.** $\vx=\vzero$, $z=0$.
2. **Shortest path.** In the residual network $G(\vx)$ (costs $c_{ij}$
   forward, $-c_{ij}$ backward) find a shortest directed path $\mathcal P$
   from $s$ to $t$. If there is none, go to Step 4.
3. **Augment.** Let $\Delta$ be the smallest residual capacity on
   $\mathcal P$. Send $\Delta$ units along $\mathcal P$ (increase forward
   arcs, cancel flow on backward arcs), add $\Delta\cdot c(\mathcal P)$ to
   $z$, and go to Step 2.
4. **Stop.** If all arcs $(s,i)$ are saturated, $\vx$ restricted to the
   original arcs is optimal with cost $z$; otherwise the problem is
   infeasible.
:::

**Why it works.** Suppose $G(\vx)$ has no negative cycle and $\mathcal P$ is a
shortest $s$–$t$ path. After augmenting, the new residual arcs are reverses of
arcs of $\mathcal P$. With $d$ the shortest-path lengths before the
augmentation, every arc of $\mathcal P$ satisfies $d_j=d_i+c_{ij}$, so its
reverse also has cost $-c_{ij}=d_i-d_j$, and every residual arc after the
augmentation still satisfies $c_{ij}+d_i-d_j\ge0$. Summing this around any
cycle shows that no cycle has negative cost: optimality
({prf:ref}`thm-04-04-negative-cycle`) is preserved. With integer data each
augmentation sends at least one unit, so at most $\sum_{b_i>0}b_i$
augmentations are needed. Because reverse arcs have negative costs, Step 2
needs a method that allows them (Bellman–Ford,
{ref}`sec-05-04-successive-approximation`);
alternatively the distances $d$ define potentials $\pi=-d$ with which all
reduced costs are nonnegative, and Dijkstra's algorithm can be run on the
reduced costs.

**The transportation example.** {numref}`tab-04-04-ssp` lists the iterations
for {prf:ref}`ex-04-04-transportation`.

```{list-table} Successive shortest paths for the transportation example
:name: tab-04-04-ssp
:header-rows: 1

* - Iteration
  - Shortest path
  - Path cost (\$/t)
  - $\Delta$ (t)
  - Cumulative cost (\$)
* - 1
  - $s\to B\to1\to t$
  - 2
  - 20
  - 40
* - 2
  - $s\to B\to3\to t$
  - 5
  - 30
  - 190
* - 3
  - $s\to A\to2\to t$
  - 6
  - 25
  - 340
* - 4
  - $s\to A\to1\to B\to3\to t$
  - $4-2+5=7$
  - 15
  - 445
```

1. The cheapest route of all is B → 1 at \$2; market 1 needs 20 t, so
   $\Delta=20$.
2. Market 1 is served. The cheapest remaining route is B → 3 at \$5; B has
   30 t left, so $\Delta=30$ and plant B is exhausted.
3. Only plant A has product. Its routes cost \$6 to market 2, \$9 to market 3,
   and \$7 to market 3 by the residual path $A\to1\to B\to3$. The cheapest is
   $A\to2$; $\Delta=25$ fills market 2.
4. Market 3 still needs 15 t. The direct arc $A\to3$ costs \$9, but the path
   $A\to1\to B\to3$ costs $4-2+5=\$7$: plant A takes over 15 t of market 1
   from plant B (the backward arc $1\to B$ cancels 15 t of $x_{B1}$), and B
   sends those 15 t to market 3. $\Delta=\min\{15,20,15\}=15$.

All supply arcs are saturated after 90 t. The optimal plan is $x_{A1}=15$,
$x_{A2}=25$, $x_{B1}=5$, $x_{B3}=45$ and $x_{A3}=x_{B2}=0$, with cost $4\cdot15+6\cdot25+2\cdot5+5\cdot45=\$445$ per week, the last
cumulative cost of the table. The potentials $\pi_A=7$, $\pi_B=5$, $\pi_1=3$,
$\pi_2=1$, $\pi_3=0$ give reduced cost 0 on the four used arcs and
$\bar c_{A3}=9-7+0=2$, $\bar c_{B2}=7-5+1=3$ on the unused ones, which proves
optimality ({prf:ref}`thm-04-04-reduced-cost`); both are positive, so the plan is
unique. Iteration 4 shows why the residual network is needed: a greedy rule
that never revises a shipment would have used $A\to3$ and paid \$30 more.

```{figure} figures/04-04-transportation.svg
:name: fig-04-04-transportation
:width: 420px
:class: dark-light
:alt: Bipartite network with plants A and B on the left, markets 1 to 3 on the right, unit costs on the arcs, and the four arcs of the optimal plan highlighted

The transportation network of {prf:ref}`ex-04-04-transportation` with unit costs
(\$ per tonne) and net supplies in brackets. Highlighted arcs carry flow in
the optimal plan: $A\to1$ 15 t, $A\to2$ 25 t, $B\to1$ 5 t, $B\to3$ 45 t.
```

(sec-04-04-cycle-cancelling)=
## Cycle canceling

The second algorithm works the other way round: it starts from any feasible
flow and improves its cost until {prf:ref}`thm-04-04-negative-cycle` certifies
optimality.

:::{prf:algorithm} Cycle canceling
:label: alg-04-04-cycle-cancelling

**Input:** a minimum-cost flow problem {eq}`eq-04-04-mcf`.
**Output:** a minimum-cost flow.

1. **Feasible flow.** Find a feasible flow $\vx$ (for instance with one
   maximum-flow computation, {ref}`sec-04-04-feasibility`). If there is none,
   stop.
2. **Detect.** Search $G(\vx)$ for a directed cycle $\mathcal W$ with
   $c(\mathcal W)<0$ (Bellman–Ford or Floyd–Warshall detect one). If there is
   none, stop: $\vx$ is optimal.
3. **Cancel.** Let $\Delta$ be the smallest residual capacity on
   $\mathcal W$; send $\Delta$ units around $\mathcal W$, which lowers the cost
   by $\Delta\abs{c(\mathcal W)}$. Go to Step 2.
:::

With integer data each cancellation lowers the cost by at least 1, so the
method terminates. Choosing the cycle with the most negative *mean* cost per
arc makes it polynomial (Goldberg & Tarjan, 1989). The method also explains the
"pushing back" in iteration 4 above: the plan before that iteration plus 15 t
on $A\to3$ would leave the residual cycle $A\to1\to B\to3\to A$ of cost
$4-2+5-9=-2$, and canceling it gives the optimum. A full run is worked in
{prf:ref}`ex-08-05-cycle-cancelling`.

(sec-04-04-network-simplex)=
## The network simplex method (outline)

The simplex method applied to {eq}`eq-04-04-mcf` can be carried out on the
network itself, without a tableau, because bases have a graph meaning.

:::{prf:theorem} Bases are spanning trees
:label: thm-04-04-spanning-tree-basis

For a connected network, a set of $\abs{\mathcal N}-1$ columns of $\mA$ (with
one redundant row removed) is a basis if and only if the corresponding arcs
form a spanning tree (directions ignored).
:::

The reason: the columns of the arcs of a cycle are linearly dependent (add
them with signs $\pm1$ according to the direction of travel and they cancel),
so a basis contains no cycle; and $\abs{\mathcal N}-1$ arcs without a cycle
form a spanning tree ({prf:ref}`thm-04-01-tree`). A basic solution therefore
has every non-tree arc at its lower or upper bound, and the tree flows are
determined by the balance equations, computed from the leaves inwards. One
iteration of the **network simplex method** is:

1. **Potentials.** Fix $\pi_r=0$ at a root and solve
   $\bar c_{ij}=c_{ij}-\pi_i+\pi_j=0$ along the tree arcs, from the root
   outwards.
2. **Pricing.** Compute $\bar c_{ij}$ for the non-tree arcs. If every arc at
   its lower bound has $\bar c_{ij}\ge0$ and every arc at its upper bound has
   $\bar c_{ij}\le0$, the conditions {eq}`eq-04-04-cs` hold: stop.
3. **Entering arc.** Otherwise choose a violating arc. Together with the tree
   it forms exactly one cycle.
4. **Ratio test and pivot.** Push flow around that cycle in the direction
   that lowers the cost until some arc of the cycle reaches a bound; it
   leaves the tree, the entering arc joins it, and the potentials are
   updated.

All arithmetic is addition and subtraction along tree paths, which is why the
method is fast in practice and why it preserves integrality. Degenerate
pivots (zero-length steps) can cycle; the rule of "strongly feasible trees"
of Cunningham (1976) prevents this. For the supply chain, the optimal tree is the
highlighted set of {numref}`fig-04-04-supply-chain`(b), $\{(1,3),(2,3),(2,5),
(5,6),(4,6)\}$; with $\pi_6=0$ step 1 gives $\pi_5=4$, $\pi_2=10$,
$\pi_3=5$, $\pi_1=9$, $\pi_4=2$, and step 2 stops at once, as
{numref}`tab-04-04-supply-chain-solution` shows.

## Summary and Common Mistakes

### Summary

- The minimum-cost flow problem {eq}`eq-04-04-mcf` minimizes
  $\sum c_{ij}x_{ij}$ subject to flow balance $\mA\vx=\vb$ and bounds
  $\vl\le\vx\le\vu$ ({prf:ref}`def-04-04-mcf`).
- It needs $\sum_ib_i=0$ (dummy nodes otherwise) and enough capacity across
  every node set {eq}`eq-04-04-gale`; the supply chain of
  {prf:ref}`ex-04-04-supply-chain` fails this at node 4 and is repaired by
  dropping one lower bound. Lower bounds are removed by
  {eq}`eq-04-04-shift`–{eq}`eq-04-04-shifted-supply`.
- The incidence matrix is totally unimodular ({prf:ref}`thm-04-04-tu`), so
  integer data give integer basic optimal flows
  ({prf:ref}`thm-04-04-integrality`).
- Shortest path, maximum flow, transportation, assignment and transshipment
  are special cases ({ref}`sec-04-04-special-cases`).
- A flow is optimal if and only if there are potentials satisfying
  {eq}`eq-04-04-cs` ({prf:ref}`thm-04-04-reduced-cost`), equivalently if its
  residual network has no negative cycle ({prf:ref}`thm-04-04-negative-cycle`).
  The supply chain costs \$1 810 at the optimum
  ({numref}`tab-04-04-supply-chain-solution`).
- Successive shortest paths ({prf:ref}`alg-04-04-ssp`) solve the transportation
  example in four augmentations at \$445 ({numref}`tab-04-04-ssp`); cycle
  canceling ({prf:ref}`alg-04-04-cycle-cancelling`) and the network simplex
  method (bases = spanning trees, {prf:ref}`thm-04-04-spanning-tree-basis`) are
  the other standard methods.

### Common mistakes

- **Unbalanced data.** If $\sum b_i\ne0$ the LP is infeasible whatever the
  network; add a dummy node and decide what its arcs cost.
- **Forgetting the lower bounds in the supplies.** After
  $x_{ij}=l_{ij}+x'_{ij}$, both end nodes change their net supply
  {eq}`eq-04-04-shifted-supply`, and the constant $\sum c_{ij}l_{ij}$ must be
  added back to the cost.
- **Sign of the reduced cost.** With $\bar c_{ij}=c_{ij}-\pi_i+\pi_j$, an arc
  at its *upper* bound is optimal when $\bar c_{ij}\le0$, not $\ge0$.
- **Ignoring backward arcs** in the residual network. They have negative cost,
  so Dijkstra's algorithm cannot be applied to them without potentials, and
  without them the successive-shortest-path method can stop at a non-optimal
  plan.
- **Assuming integrality is automatic.** It needs integer $\vb$, $\vl$ and $\vu$ and
  a *basic* optimal solution; an interior-point solver may return a
  fractional point on an optimal face.

## References

- Ahuja, R. K., Magnanti, T. L., & Orlin, J. B. (1993). *Network flows:
  Theory, algorithms, and applications*. Prentice Hall. (Chapter 9:
  optimality conditions, duality, cycle-canceling and successive
  shortest path algorithms; Chapter 11: network simplex algorithm; Chapter
  12: assignment problems.)
- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear programming
  and network flows* (4th ed.). Wiley. (Chapter 9: minimal-cost network flow
  problems, total unimodularity and the network simplex method; Chapter 10:
  the transportation and assignment problems.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Chapter 7: network flow problems, the
  network simplex algorithm and integrality.)
- Cunningham, W. H. (1976). A network simplex method. *Mathematical
  Programming, 11*(1), 105–116.
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapters 14–17 and 19–20: the transportation, assignment
  and transshipment problems, networks, maximal flows and the primal–dual
  method for transportation problems.)
- Gale, D. (1957). A theorem on flows in networks. *Pacific Journal of
  Mathematics, 7*(2), 1073–1082.
- Goldberg, A. V., & Tarjan, R. E. (1989). Finding minimum-cost circulations
  by canceling negative cycles. *Journal of the ACM, 36*(4), 873–886.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 9: transportation and
  assignment problems; Sections 10.6–10.7: minimum cost flow problem and the
  network simplex method.)
- Hoffman, A. J. (1960). Some recent applications of the theory of linear
  inequalities to extremal combinatorial analysis. In R. Bellman & M. Hall
  (Eds.), *Combinatorial analysis* (pp. 113–127). American Mathematical
  Society.
- Kuhn, H. W. (1955). The Hungarian method for the assignment problem.
  *Naval Research Logistics Quarterly, 2*(1–2), 83–97.
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 5: transportation model and its variants; Chapter 6:
  network models.)
