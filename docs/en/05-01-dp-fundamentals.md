# Dynamic Programming Fundamentals

## Overview

Many planning problems ask for a *sequence* of decisions: which road to take
at each junction, when to replace a machine, how much of a budget to give to
each project. The decisions interact, because each one changes the situation
in which the next one is made, so the problem cannot be solved one decision
at a time. **Dynamic programming** (DP) solves such problems by splitting
them into a family of smaller subproblems, one for every possible situation,
and solving the subproblems in an order in which each one needs only answers
that are already known. The justification is Bellman's **principle of
optimality**: whatever happened before, the remaining decisions of an
optimal plan must be optimal for the situation they start from. This chapter
develops the method on a single network with sixteen nodes. It states and
proves the principle of optimality for additive costs and shows when it
fails. It then names the parts of every DP model (stage, state, decision,
transition, stage cost, value function, policy function, recurrence,
boundary condition and answer cell) and solves the network by backward and
by forward recursion. A turn penalty shows how the state must sometimes be
enlarged (state augmentation). The chapter ends with operation counts, the
curse of dimensionality and a short comparison of deterministic and
stochastic DP. The applications follow in {doc}`05-02-dp-equipment-replacement`,
{doc}`05-03-dp-resource-allocation` and {doc}`05-04-dp-shortest-paths`. Only
the vocabulary of networks (nodes, arcs, paths) from {doc}`04-01-networks-mst`
is assumed.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- explain the principle of optimality, prove it for additive costs and
  recognize objectives and state definitions for which it fails;
- formulate a DP model by naming its stage, state, decision, transition,
  stage cost, value function, policy function, recurrence, boundary
  condition and answer cell;
- solve a staged network by backward and by forward recursion, record the
  value and policy tables, and trace the optimal decisions;
- enlarge the state when the cost of a decision depends on the past, and
  write the new recurrence and boundary conditions; and
- count the additions and comparisons of a DP and of complete enumeration
  and explain the curse of dimensionality.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $k=0,1,\dots,N$ | stage index (number of decisions already made); $N$ = number of stages |
| $y$, $\mathcal Y_k$ | state, set of states at stage $k$ |
| $d$, $\mathcal D_k(y)$ | decision, set of decisions allowed in state $y$ at stage $k$ |
| $t_k(y,d)$ | transition function: state at stage $k+1$ after decision $d$ in state $y$ at stage $k$ |
| $c_k(y,d)$ | stage cost of decision $d$ in state $y$ at stage $k$ |
| $V_k(y)$ | optimal value function: least cost from state $y$ at stage $k$ to the end (backward recursion) |
| $P_k(y)$ | optimal policy function: a decision that attains $V_k(y)$ |
| $F_k(y)$ | least cost from the initial state to state $y$ at stage $k$ (forward recursion) |
| $(k,y)$ | vertex of the lattice network: column $k$, level $y$ |
| $c^{\uparrow}(k,y)$, $c^{\downarrow}(k,y)$ | cost of the arc from $(k,y)$ up to $(k+1,y+1)$ and down to $(k+1,y-1)$; $\infty$ if there is no such arc |
| $\mathrm U$, $\mathrm D$ | the decisions "up" and "down" |
| $\omega\in\{\mathrm U,\mathrm D\}$ | direction in which a vertex is left (augmented state) |
| $V_k(y,\omega)$ | value function of the turn-penalty problem |

(sec-05-01-multistage)=
## Multistage decision problems

A **multistage decision problem** is an optimization problem in which
decisions are taken one after another, each decision changes the situation
of the decision maker, and the total cost depends on the whole sequence.
Three ingredients make such problems hard to solve by inspection: the number
of decision sequences grows exponentially with their length, a decision that
looks cheap now may lead into an expensive situation later, and the
decisions cannot be separated because each one determines what is possible
next. The simplest multistage problems are path problems on a network whose
arcs all point "forwards", and they are used throughout this chapter.

:::{prf:example} Fastest route through a staged network
:label: ex-05-01-diamond

A traveler must go from node $A$ to node $B$ through the network of
{numref}`fig-05-01-diamond-network`. Every arc points to the right, and the
number on an arc is its travel time in minutes. Each node lies in one of
seven columns, and every arc joins a node in column $k$ to a node in column
$k+1$, one level up or one level down. Find a path from $A$ to $B$ of least
total travel time.

```{figure} figures/05-01-diamond-network.svg
:name: fig-05-01-diamond-network
:width: 640px
:class: dark-light
:alt: A diamond-shaped lattice network of sixteen nodes from A to B with travel times on its arcs, stage lines k = 0 to 6, and the fastest path A-C-F-J-M-P-B highlighted

The staged network of {prf:ref}`ex-05-01-diamond`. Each node $(k,y)$ sits in
column (stage) $k$ at level $y$; every arc goes one column to the right and
one level up or down. Black numbers are travel times in minutes; the bold
orange number above a node is its value $V_k(y)$, the least time from that
node to $B$, computed in {ref}`sec-05-01-backward`. The fastest path
$A$–$C$–$F$–$J$–$M$–$P$–$B$ (13 min) is highlighted.
```
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Complete enumeration.** Every path from $A$ to $B$ uses six arcs, three up
and three down, in some order, and every such order stays inside the
diamond. There are therefore $\binom{6}{3}=20$ paths. Evaluating one path
takes 5 additions (six numbers are added), so all paths take
$20\times5=100$ additions, and choosing the smallest of 20 totals takes 19
comparisons. This is manageable here, but a lattice with 30 columns already
has $\binom{30}{15}\approx1.55\times10^{8}$ paths
({ref}`sec-05-01-effort`).

The fastest path is $A$–$C$–$F$–$J$–$M$–$P$–$B$, of 13 minutes. Dynamic
programming finds it without listing the paths, by backward recursion
({prf:ref}`ex-05-01-backward`) or by forward recursion
({prf:ref}`ex-05-01-forward`).
:::

**A better idea.** Complete enumeration, as in the solution of
{prf:ref}`ex-05-01-diamond`, evaluates every path. Look instead at the first
decision, at $A$: go up to $C$ (1 min) or down to $D$ (0 min). Suppose we already knew the least time from $C$ to
$B$ and the least time from $D$ to $B$. Then the best route from $A$ would
follow at once: compare $1+(\text{best time from }C)$ with
$0+(\text{best time from }D)$. Inferior paths from $C$ or from $D$ to $B$ are
irrelevant, because whichever way we leave $A$, we should continue in the
best possible way from where we arrive. The best time from $C$ depends in
the same way on the best times from $E$ and $F$, the best time from $D$ on
those from $F$ and $G$, and so on to the nodes next to $B$, whose best times
are obvious. The work is done once per *node* instead of once per *path*,
and nodes such as $F$, which are shared by many paths, are evaluated only
once. This is the whole idea of dynamic programming; the next two sections
make it precise.

(sec-05-01-principle)=
## The principle of optimality

Bellman (1957, Ch. III, §3) stated the principle that justifies the
argument above:

> An optimal policy has the property that whatever the initial state and
> initial decision are, the remaining decisions must constitute an optimal
> policy with regard to the state resulting from the first decision.

For a path problem it reads: *the best path from $A$ to $B$ has the property
that, whatever its first arc, the rest of the path, from the node it reaches
to $B$, is a best path from that node to $B$.* An equivalent form: every
tail of an optimal trajectory is optimal for the problem that starts at the
beginning of that tail. Writing $V(X)$ for the least time from node $X$ to
$B$, the principle gives at once

$$
\begin{aligned}
V(A)&=\min\{1+V(C),\ 0+V(D)\},\\
V(C)&=\min\{5+V(E),\ 4+V(F)\},\\
V(D)&=\min\{7+V(F),\ 3+V(G)\},
\end{aligned}
$$

a set of relations in which every value depends only on values one column
further right. Once $V(P)$ and $V(Q)$ are known, the relations can be worked
backwards to $V(A)$.

### Why it holds for additive costs

The principle is not a law of nature. It holds because of the way the cost
of a decision sequence is built, and it is worth seeing exactly where.

:::{prf:theorem} Principle of optimality for additive costs
:label: thm-05-01-principle

Consider a multistage problem with stages $k=0,\dots,N-1$, in which the
state $y_{k+1}=t_k(y_k,d_k)$ after a decision depends only on the current
state and decision, the set $\mathcal D_k(y_k)$ of allowed decisions depends
only on the current state, and the total cost is the sum
$\sum_{k=0}^{N-1}c_k(y_k,d_k)$ of stage costs that depend only on the
current state and decision. Let $(d_0^*,\dots,d_{N-1}^*)$ be an optimal
decision sequence from the initial state $y_0$, with states
$y_0^*=y_0,y_1^*,\dots,y_N^*$. Then for every $j$, the tail
$(d_j^*,\dots,d_{N-1}^*)$ is an optimal decision sequence for the problem
that starts in state $y_j^*$ at stage $j$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Split the total cost of the optimal sequence into a head and a tail,

$$
\sum_{k=0}^{N-1}c_k(y_k^*,d_k^*)
=\underbrace{\sum_{k=0}^{j-1}c_k(y_k^*,d_k^*)}_{\text{head}}
+\underbrace{\sum_{k=j}^{N-1}c_k(y_k^*,d_k^*)}_{\text{tail}}.
$$

Suppose the tail were not optimal from $(j,y_j^*)$: some feasible sequence
$(d_j',\dots,d_{N-1}')$ from state $y_j^*$ at stage $j$ has a smaller tail
cost. Because feasibility and costs from stage $j$ on depend only on the
state $y_j^*$, and not on how it was reached, the head followed by the new
tail is a feasible sequence from $y_0$. Its cost is the same head cost plus
a smaller tail cost, which is less than the optimal cost, a contradiction.
:::

The proof uses three properties, and each is necessary:

1. **The state summarizes the past.** Given the state at stage $j$, the
   allowed future decisions, their costs and the resulting states do not
   depend on earlier states or decisions. (This is the *Markov* or
   *no-memory* property of the state.)
2. **Separability.** The total cost is the head cost combined with a tail
   cost that depends only on the state at stage $j$ and the later
   decisions.
3. **Monotonicity.** The combination is nondecreasing in the tail cost, so
   a better tail never makes the total worse. Addition has this property;
   so do $\max\{a,\cdot\}$ and multiplication by a positive number, which is
   why minimax and product objectives can also be handled by DP
   ({doc}`08-06-examples-dp-fundamentals`).

### When it fails

The principle fails, or more precisely the *chosen definition of the state*
fails, as soon as one of the three properties is lost.

* **The state forgets something that matters.** Add a time penalty of 3 min
  for every change between an up arc and a down arc of
  {numref}`fig-05-01-diamond-network`. Take the node as state. From $F$, the
  best continuation now depends on how $F$ was reached: a traveler who
  arrived moving down can leave down without penalty, while one who arrived
  moving up pays 3 min to leave down. With the values computed in
  {ref}`sec-05-01-augmentation`, the least remaining time from $F$ is 12 min
  for a traveler arriving downwards but 14 min for one arriving upwards.
  "The best path from $F$" is no longer a well-defined subproblem, and a
  recurrence on nodes alone gives wrong answers. The cure is to put the
  missing information into the state.
* **The objective is not separable.** Suppose the traveler wants to arrive
  as close as possible to 10 min after leaving, so the cost of a path of
  length $T$ is $\abs{T-10}$. The best tail from a node depends on the time
  already spent, so again the node alone is not a state; the elapsed time
  must be added to it.
* **The combination is not monotone.** If the cost of a path is the product
  of arc numbers that may be negative, a tail with the *largest* product may
  be what the whole path needs (a negative head turns the largest tail into
  the smallest total). The remedy, carrying both the smallest and the
  largest tail product, is worked out in {doc}`08-06-examples-dp-fundamentals`.

In each case the principle is rescued by a richer state. This is the
general lesson: the principle of optimality is a property of a *model*, and
building a DP model means choosing a state that makes it true.

(sec-05-01-model)=
## Elements of a DP model

The diamond example has six features that every DP model shares:

1. The problem divides into **stages**, with a decision at each stage. Here
   the stages are the columns of the network and the decision is where to
   go next.
2. Each stage has a number of **states**, the possible situations at that
   stage. Here the state is the node reached.
3. A decision transforms the current state into a state of the next stage:
   the decision "up" at $(k,y)$ leads to $(k+1,y+1)$.
4. Given the current state, the best remaining decisions do not depend on
   the previous states or decisions: it does not matter *how* a node was
   reached, only *that* it was reached.
5. A recursive relation gives the optimal decision for every state of stage
   $k$ once stage $k+1$ has been solved.
6. The last stage can be solved by itself.

:::{prf:definition} Elements of a DP model
:label: def-05-01-dp-model

A (deterministic, finite-horizon) **DP model** consists of

- a **stage** index $k=0,1,\dots,N$, which counts the decisions already
  made;
- for every stage a set $\mathcal Y_k$ of **states** $y$, describing
  everything about the past that matters for the future;
- for every state a set $\mathcal D_k(y)$ of allowed **decisions** $d$;
- a **transition function** $t_k(y,d)\in\mathcal Y_{k+1}$, the state that
  decision $d$ leads to;
- a **stage cost** $c_k(y,d)$, the cost incurred by decision $d$;
- the **optimal value function** $V_k(y)$, the least cost of completing
  the process from state $y$ at stage $k$ (its *arguments* $(k,y)$ name a
  subproblem);
- the **optimal policy function** $P_k(y)$, a decision that achieves
  $V_k(y)$;
- the **recurrence relation** that expresses $V_k$ through $V_{k+1}$;
- the **boundary conditions**, the values of $V$ that follow from the
  problem statement without computation; and
- the **answer cell**, the argument $(k,y)$ whose value solves the original
  problem.
:::

With these elements the principle of optimality becomes the recurrence

$$
V_k(y)=\min_{d\in\mathcal D_k(y)}\bigl\{c_k(y,d)+V_{k+1}\bigl(t_k(y,d)\bigr)\bigr\},
\qquad y\in\mathcal Y_k,\ k=N-1,\dots,0,
$$ (eq-05-01-recurrence)

with boundary condition and answer

$$
\begin{gathered}
V_N(y)=g(y)\qquad (y\in\mathcal Y_N),\\
\text{answer: } V_0(y_0),
\end{gathered}
$$ (eq-05-01-boundary)

where $g(y)$ is the terminal cost of ending in state $y$ (zero in a path
problem that ends at a fixed node). The policy function records a
minimizer: $P_k(y)\in\argmin_{d}\{c_k(y,d)+V_{k+1}(t_k(y,d))\}$.
{numref}`fig-05-01-stage-state` shows one stage of the recurrence.

```{figure} figures/05-01-stage-state.svg
:name: fig-05-01-stage-state
:width: 520px
:class: dark-light
:alt: Schematic of one DP stage: a state y at stage k with two decisions leading to two states at stage k plus 1, labeled with stage costs and values

One stage of a DP model. Each decision $d$ in state $y$ costs $c_k(y,d)$ and
leads to the state $t_k(y,d)$ of the next stage, whose value $V_{k+1}$ is
already known. The recurrence {eq}`eq-05-01-recurrence` compares
$c_k(y,d)+V_{k+1}(t_k(y,d))$ over the decisions; the better decision
(highlighted) is the policy $P_k(y)$.
```

The recurrence is a theorem, not a definition: {prf:ref}`thm-05-01-principle`
shows that a best completion from $(k,y)$ consists of some first decision
$d$ followed by a best completion from $t_k(y,d)$, and trying all first
decisions gives {eq}`eq-05-01-recurrence`. A model is correct exactly when the
state has the no-memory property 4 above.

:::{prf:algorithm} Building and solving a DP model
:label: alg-05-01-dp-steps

1. **Arguments.** Choose the stage variable and the state variables, so
   that each pair $(k,y)$ describes a subproblem and the state summarizes
   everything about the past that matters for the future.
2. **Value function.** Define $V_k(y)$ in words, precisely enough that the
   principle of optimality applies to it.
3. **Recurrence.** Write {eq}`eq-05-01-recurrence` for this problem: the
   decisions, their costs and the transitions.
4. **Boundary conditions.** Write the values that need no computation, and
   name the answer cell.
5. **Solve.** Starting from the boundary conditions, evaluate the
   recurrence stage by stage, recording $V_k(y)$ and $P_k(y)$ for every
   state.
6. **Trace.** Starting from the answer cell, follow the policy function
   forward to read off the optimal decisions.
:::

Choosing the stage and state variables is the creative step; the rest is
routine. Two rules of thumb guide it. If *too few* arguments are used, the
state forgets something that matters and no correct recurrence can be
written (the turn penalty of {ref}`sec-05-01-augmentation` is the standard
example). If *too many* are used, a correct recurrence exists, but the
table becomes needlessly large ({ref}`sec-05-01-effort`).

**The lattice in coordinates.** Give vertex $X$ of
{numref}`fig-05-01-diamond-network` the coordinates $(k,y)$: $k$ is its column
and $y$ its level, so $A=(0,0)$, $C=(1,1)$, $D=(1,-1)$, $H=(3,3)$,
$K=(3,-3)$ and $B=(6,0)$. Write $c^{\uparrow}(k,y)$ for the time of the arc
from $(k,y)$ to $(k+1,y+1)$ and $c^{\downarrow}(k,y)$ for the time of the
arc from $(k,y)$ to $(k+1,y-1)$, with the value $\infty$ when the arc does
not exist (for instance $c^{\uparrow}(3,3)=\infty$ at the top vertex $H$).
The DP model is:

- stage $k=0,\dots,6$ (the column), state $y$ (the level), decision
  $d\in\{\mathrm U,\mathrm D\}$, transition $t_k(y,\mathrm U)=y+1$,
  $t_k(y,\mathrm D)=y-1$, stage cost $c^{\uparrow}(k,y)$ or
  $c^{\downarrow}(k,y)$;
- value function $V_k(y)$ = least time from vertex $(k,y)$ to $B=(6,0)$;
- recurrence

  $$
  V_k(y)=\min\bigl\{c^{\uparrow}(k,y)+V_{k+1}(y+1),\ c^{\downarrow}(k,y)+V_{k+1}(y-1)\bigr\};
  $$ (eq-05-01-lattice-backward)

- boundary condition $V_6(0)=0$ (equivalently $V_5(1)=2$ and
  $V_5(-1)=1$, the only arcs into $B$); answer cell $V_0(0)$.

The value $\infty$ makes a missing arc lose every comparison, so one
formula covers interior and boundary vertices alike.

(sec-05-01-backward)=
## Backward recursion

The recurrence {eq}`eq-05-01-lattice-backward` is solved from the last stage to
the first, which is why it is called **backward recursion**: it answers the
question "what is the best route from $(k,y)$ to the end?" for every
vertex, and the whole problem is solved when $(k,y)$ is the origin.

:::{prf:example} Backward solution of the staged network
:label: ex-05-01-backward

Solve {prf:ref}`ex-05-01-diamond` with {eq}`eq-05-01-lattice-backward`.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Stage 5.** $V_5(1)=c^{\downarrow}(5,1)+V_6(0)=2$ at $P$ and
$V_5(-1)=c^{\uparrow}(5,-1)+V_6(0)=1$ at $Q$.

**Stage 4.** At $M=(4,0)$ both arcs exist:
$V_4(0)=\min\{2+V_5(1),\ 8+V_5(-1)\}=\min\{4,\ 9\}=4$, so the best decision
is up, to $P$. At $L$ and $N$ only one arc exists.

The remaining stages follow the same pattern. {numref}`tab-05-01-backward`
lists every option and the result.

```{table} Backward recursion on the staged network: for each vertex the options $c+V_{k+1}$ (up; down), the value $V_k(y)$ and the policy $P_k(y)$ with the node it leads to.
:name: tab-05-01-backward

| Stage $k$ | Node $(k,y)$ | Up: $c^{\uparrow}+V_{k+1}(y+1)$ | Down: $c^{\downarrow}+V_{k+1}(y-1)$ | $V_k(y)$ | $P_k(y)$ |
|---|---|---|---|---|---|
| 6 | $B\ (6,0)$ | – | – | 0 | – |
| 5 | $P\ (5,1)$ | – | $2+0=2$ | 2 | D ($B$) |
| 5 | $Q\ (5,-1)$ | $1+0=1$ | – | 1 | U ($B$) |
| 4 | $L\ (4,2)$ | – | $5+2=7$ | 7 | D ($P$) |
| 4 | $M\ (4,0)$ | $2+2=4$ | $8+1=9$ | 4 | U ($P$) |
| 4 | $N\ (4,-2)$ | $4+1=5$ | – | 5 | U ($Q$) |
| 3 | $H\ (3,3)$ | – | $3+7=10$ | 10 | D ($L$) |
| 3 | $I\ (3,1)$ | $3+7=10$ | $4+4=8$ | 8 | D ($M$) |
| 3 | $J\ (3,-1)$ | $2+4=6$ | $2+5=7$ | 6 | U ($M$) |
| 3 | $K\ (3,-3)$ | $2+5=7$ | – | 7 | U ($N$) |
| 2 | $E\ (2,2)$ | $2+10=12$ | $1+8=9$ | 9 | D ($I$) |
| 2 | $F\ (2,0)$ | $1+8=9$ | $2+6=8$ | 8 | D ($J$) |
| 2 | $G\ (2,-2)$ | $5+6=11$ | $4+7=11$ | 11 | U or D ($J$ or $K$) |
| 1 | $C\ (1,1)$ | $5+9=14$ | $4+8=12$ | 12 | D ($F$) |
| 1 | $D\ (1,-1)$ | $7+8=15$ | $3+11=14$ | 14 | D ($G$) |
| 0 | $A\ (0,0)$ | $1+12=13$ | $0+14=14$ | 13 | U ($C$) |
```

**Answer and trace.** $V_0(0)=13$ min. Following the policy from the
answer cell: $P_0(0)=\mathrm U$ to $C$, $P_1(1)=\mathrm D$ to $F$,
$P_2(0)=\mathrm D$ to $J$, $P_3(-1)=\mathrm U$ to $M$,
$P_4(0)=\mathrm U$ to $P$, $P_5(1)=\mathrm D$ to $B$. The fastest path is
$A$–$C$–$F$–$J$–$M$–$P$–$B$ with times $1+4+2+2+2+2=13$ min, as
highlighted in {numref}`fig-05-01-diamond-network`.

**Check.** Enumerating all 20 paths gives the same minimum, 13 min, on the
same path only; the next best paths take 14 min
($A$–$C$–$F$–$I$–$M$–$P$–$B$ and $A$–$C$–$F$–$J$–$N$–$Q$–$B$).

**Work.** Each arc is used in exactly one addition: 24 additions for 24
arcs (22 if the two values of stage 5 are taken as boundary conditions).
A comparison is needed at each vertex with two outgoing arcs, and there are
9 such vertices ($A$, $C$, $D$, $E$, $F$, $G$, $I$, $J$, $M$): 9
comparisons. Complete enumeration needed 100 additions and 19 comparisons.
:::

The backward table {numref}`tab-05-01-backward` holds more than the answer.
Its policy column is an optimal
**feedback rule**: for *every* vertex it says which way to go, so a
traveler pushed off the optimal path (by a closed road, say) can continue
optimally from any vertex without a new computation. The tie at $G$
shows that the policy function may be set-valued; any choice gives an
optimal continuation.

(sec-05-01-forward)=
## Forward recursion

The recursion can also be organized from the origin outwards. Instead of
"what is the best route from $(k,y)$ to the end?", ask "what is the best
route from the start to $(k,y)$?". Let $F_k(y)$ be the least time from
$A=(0,0)$ to vertex $(k,y)$. A best path to $(k,y)$ arrives either from
$(k-1,y-1)$ by an up arc or from $(k-1,y+1)$ by a down arc, and by the
principle of optimality its part up to that predecessor is a best path to
the predecessor. Hence

$$
F_k(y)=\min\bigl\{F_{k-1}(y-1)+c^{\uparrow}(k-1,y-1),\ F_{k-1}(y+1)+c^{\downarrow}(k-1,y+1)\bigr\},
$$ (eq-05-01-lattice-forward)

for $k=1,\dots,6$, with boundary condition $F_0(0)=0$ and answer cell
$F_6(0)$. The policy now records the *predecessor* of each vertex on a best
path to it.

:::{prf:example} Forward solution of the staged network
:label: ex-05-01-forward

Solve {prf:ref}`ex-05-01-diamond` with {eq}`eq-05-01-lattice-forward`.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

{numref}`tab-05-01-forward` gives the result stage by stage.

```{table} Forward recursion on the staged network: options $F_{k-1}+c$ from the lower-left and upper-left predecessors, the value $F_k(y)$ and the best predecessor.
:name: tab-05-01-forward

| Stage $k$ | Node | From below-left: $F_{k-1}(y-1)+c^{\uparrow}$ | From above-left: $F_{k-1}(y+1)+c^{\downarrow}$ | $F_k(y)$ | Predecessor |
|---|---|---|---|---|---|
| 0 | $A$ | – | – | 0 | – |
| 1 | $C$ | $0+1=1$ | – | 1 | $A$ |
| 1 | $D$ | – | $0+0=0$ | 0 | $A$ |
| 2 | $E$ | $1+5=6$ | – | 6 | $C$ |
| 2 | $F$ | $0+7=7$ | $1+4=5$ | 5 | $C$ |
| 2 | $G$ | – | $0+3=3$ | 3 | $D$ |
| 3 | $H$ | $6+2=8$ | – | 8 | $E$ |
| 3 | $I$ | $5+1=6$ | $6+1=7$ | 6 | $F$ |
| 3 | $J$ | $3+5=8$ | $5+2=7$ | 7 | $F$ |
| 3 | $K$ | – | $3+4=7$ | 7 | $G$ |
| 4 | $L$ | $6+3=9$ | $8+3=11$ | 9 | $I$ |
| 4 | $M$ | $7+2=9$ | $6+4=10$ | 9 | $J$ |
| 4 | $N$ | $7+2=9$ | $7+2=9$ | 9 | $J$ or $K$ |
| 5 | $P$ | $9+2=11$ | $9+5=14$ | 11 | $M$ |
| 5 | $Q$ | $9+4=13$ | $9+8=17$ | 13 | $N$ |
| 6 | $B$ | $13+1=14$ | $11+2=13$ | 13 | $P$ |
```

The answer is $F_6(0)=13$ min. Tracing predecessors backwards from $B$:
$B\leftarrow P\leftarrow M\leftarrow J\leftarrow F\leftarrow C\leftarrow A$,
the same path as before.

**Check.** For every vertex, $F_k(y)+V_k(y)$ is the length of the best
$A$–$B$ path *through* that vertex. On the optimal path the sum is 13 at
every vertex ($A$: $0+13$, $C$: $1+12$, $F$: $5+8$, $J$: $7+6$, $M$: $9+4$,
$P$: $11+2$, $B$: $13+0$); at every other vertex it is larger (for example
$H$: $8+10=18$), which confirms both tables.
:::

Both recursions do the same amount of work and give the same optimal value,
because the forward recursion is the backward recursion applied to the
network with all arcs reversed and the roles of $A$ and $B$ exchanged. They
answer different families of questions, however. The backward table gives
the best route *to the fixed terminal* from every vertex; the forward table
gives the best route *from the fixed origin* to every vertex. Backward
recursion is the natural choice when the end of the process is fixed and
the start may vary (a policy that must work from any state); forward
recursion when the start is fixed and several possible ends are of interest,
or when the data of later stages become known only as time goes on.

(sec-05-01-augmentation)=
## State augmentation

The no-memory property of the state is easily lost when the cost of a
decision depends on earlier decisions. The remedy is to **augment** the
state with whatever part of the past the future depends on.

:::{prf:example} Staged network with a turn penalty
:label: ex-05-01-turn-penalty

In the network of {prf:ref}`ex-05-01-diamond`, every change of direction costs
an extra 3 min: whenever the traveler arrives at a vertex on an up arc and
leaves it on a down arc, or arrives on a down arc and leaves on an up arc.
Find the fastest route from $A$ to $B$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Choice of state.** As {ref}`sec-05-01-principle` showed, the vertex alone is
not enough: the cost of leaving $F$ downwards is $c^{\downarrow}(2,0)$ or
$c^{\downarrow}(2,0)+3$ depending on how $F$ was entered. For a backward
recursion it is convenient to record the direction in which a vertex is
*left*: when the next vertex is processed, the direction of arrival is
then known. The state becomes $(y,\omega)$, with
$\omega\in\{\mathrm U,\mathrm D\}$.

**Value function.** $V_k(y,\omega)$ is the least total time (travel times
plus turn penalties) from vertex $(k,y)$ to $B$, given that the traveler
leaves $(k,y)$ in direction $\omega$. Penalties for turns at $(k,y)$ itself
are *not* included; they belong to the decision made on arrival.

**Recurrence.** Leaving upwards takes the traveler to $(k+1,y+1)$, where
the traveler either continues upwards (no penalty) or turns down (penalty 3), and
symmetrically for leaving downwards:

$$
\begin{aligned}
V_k(y,\mathrm U)&=c^{\uparrow}(k,y)+\min\bigl\{V_{k+1}(y+1,\mathrm U),\ V_{k+1}(y+1,\mathrm D)+3\bigr\},\\
V_k(y,\mathrm D)&=c^{\downarrow}(k,y)+\min\bigl\{V_{k+1}(y-1,\mathrm U)+3,\ V_{k+1}(y-1,\mathrm D)\bigr\},
\end{aligned}
$$ (eq-05-01-turn)

for $k=4,3,\dots,0$.

**Boundary conditions.** At stage 5 the traveler leaves on the last arc,
into $B$, where no further decision (and no turn) occurs:
$V_5(1,\mathrm D)=c^{\downarrow}(5,1)=2$, $V_5(-1,\mathrm U)=1$, and
$V_5(1,\mathrm U)=V_5(-1,\mathrm D)=\infty$ (no such arcs).

**Answer cell.** At $A$ there is no arrival direction and no turn, so the
answer is $\min\{V_0(0,\mathrm U),V_0(0,\mathrm D)\}$.

**Computation.** {numref}`tab-05-01-turn` lists the values; the letter in
brackets is the direction chosen at the *next* vertex. Two entries as an
illustration:

$$
\begin{aligned}
V_4(0,\mathrm U)&=c^{\uparrow}(4,0)+\min\{V_5(1,\mathrm U),\ V_5(1,\mathrm D)+3\}\\
&=2+\min\{\infty,\ 5\}=7\ (\mathrm D),\\
V_3(-1,\mathrm D)&=c^{\downarrow}(3,-1)+\min\{V_4(-2,\mathrm U)+3,\ V_4(-2,\mathrm D)\}\\
&=2+\min\{8,\ \infty\}=10\ (\mathrm U).
\end{aligned}
$$

```{table} Turn-penalty recursion: $V_k(y,\mathrm U)$ and $V_k(y,\mathrm D)$ for every vertex, with the direction chosen at the next vertex in brackets; "–" means that the vertex has no arc in that direction ($V=\infty$).
:name: tab-05-01-turn

| Stage $k$ | Node $(k,y)$ | $V_k(y,\mathrm U)$ | $V_k(y,\mathrm D)$ |
|---|---|---|---|
| 5 | $P\ (5,1)$ | – | 2 (end) |
| 5 | $Q\ (5,-1)$ | 1 (end) | – |
| 4 | $L\ (4,2)$ | – | 7 (D) |
| 4 | $M\ (4,0)$ | 7 (D) | 12 (U) |
| 4 | $N\ (4,-2)$ | 5 (U) | – |
| 3 | $H\ (3,3)$ | – | 10 (D) |
| 3 | $I\ (3,1)$ | 13 (D) | 14 (U) |
| 3 | $J\ (3,-1)$ | 9 (U) | 10 (U) |
| 3 | $K\ (3,-3)$ | 7 (U) | – |
| 2 | $E\ (2,2)$ | 15 (D) | 15 (D) |
| 2 | $F\ (2,0)$ | 14 (U) | 12 (D) |
| 2 | $G\ (2,-2)$ | 14 (U) | 14 (U) |
| 1 | $C\ (1,1)$ | 20 (U) | 16 (D) |
| 1 | $D\ (1,-1)$ | 21 (U) | 17 (D) |
| 0 | $A\ (0,0)$ | 20 (D) | 17 (D) |
```

**Answer and trace.** $\min\{V_0(0,\mathrm U),V_0(0,\mathrm D)\}=\min\{20,17\}=17$
min: leave $A$ downwards. The brackets then give D at $D$, D at $G$, U at
$K$ (a turn), U at $N$ and U at $Q$: the path is
$A$–$D$–$G$–$K$–$N$–$Q$–$B$ ({numref}`fig-05-01-turn-penalty-path`), with
travel time $0+3+4+2+4+1=14$ min and one turn, at $K$, so the total is
$14+3=17$ min.

**Check.** Enumerating the 20 paths with their turns gives the same optimum,
unique, followed by two paths of 20 min
($A$–$C$–$F$–$J$–$N$–$Q$–$B$ and $A$–$D$–$G$–$J$–$M$–$P$–$B$). The path
that was optimal without the penalty turns three times (at $C$, $J$ and
$P$) and now costs $13+9=22$ min.

```{figure} figures/05-01-turn-penalty-path.svg
:name: fig-05-01-turn-penalty-path
:width: 640px
:class: dark-light
:alt: The staged network with pairs of values beside each node and the path A-D-G-K-N-Q-B highlighted, with one turn at K

Turn-penalty problem of {prf:ref}`ex-05-01-turn-penalty`. Beside each node is
the pair $[V_k(y,\mathrm U),\ V_k(y,\mathrm D)]$ ("–" when there is no arc
in that direction). The optimal path $A$–$D$–$G$–$K$–$N$–$Q$–$B$ goes
straight down to $K$, turns once (penalty 3) and goes straight up to $B$:
$14+3=17$ min.
```
:::

Three remarks on the formulation.

* **Leaving or arriving.** One could equally record the direction in which a
  vertex was *entered*. The recurrence then charges the penalty at the
  vertex being processed rather than at the next one. Both models are
  correct; what matters is that the definition of $V$ and the recurrence
  use the same convention.
* **A boundary at $B$.** Some texts write the boundary conditions as
  $V_6(0,\mathrm U)=V_6(0,\mathrm D)=0$ and apply {eq}`eq-05-01-turn` at stage
  5 as well. This gives $B$ a fictitious "leaving direction", but the
  minimum in {eq}`eq-05-01-turn` then simply picks the branch without a
  penalty (for example $V_5(1,\mathrm D)=2+\min\{0+3,\ 0\}=2$), so the
  numbers are the same. Stating the boundary at stage 5, as in the solution of
  {prf:ref}`ex-05-01-turn-penalty`, avoids
  the artificial direction.
* **The price of augmentation.** The augmented model has two states per
  vertex instead of one, so it needs about twice the work: 15 vertices with
  outgoing arcs give 30 states (of which 6 are infinite), each evaluated
  with one addition for the arc plus up to two for the continuation.
  Augmentation always multiplies the size of the state space; it should add
  exactly the information the future needs and nothing more.

(sec-05-01-effort)=
## Computational effort and the curse of dimensionality

**Lattice networks.** Consider the diamond lattice with $N$ columns of arcs
($N$ even), from $A=(0,0)$ to $B=(N,0)$, so that each side of the diamond
has $m=N/2$ arcs. A path consists of $m$ up arcs and $m$ down arcs in some
order, so there are $\binom{N}{m}$ paths. Complete enumeration needs $N-1$
additions per path and one comparison per path but one:

$$
\begin{aligned}
\text{enumeration: }&(N-1)\binom{N}{N/2}\ \text{additions},\\
&\binom{N}{N/2}-1\ \text{comparisons}.
\end{aligned}
$$ (eq-05-01-count-enumeration)

Backward recursion uses each arc in one addition and makes one comparison at
every vertex with two outgoing arcs. The lattice has $(m+1)^2$ vertices,
$2m(m+1)$ arcs, and $m^2$ vertices with two outgoing arcs, so

$$
\begin{aligned}
\text{DP: }&2m(m+1)=\tfrac12N(N+2)\ \text{additions},\\
&m^2=\tfrac14N^2\ \text{comparisons}.
\end{aligned}
$$ (eq-05-01-count-dp)

For $N=6$ these are 24 additions and 9 comparisons, as counted in
{prf:ref}`ex-05-01-backward`. {numref}`tab-05-01-counts` shows how differently
the two counts grow: enumeration grows like $2^N/\sqrt N$, DP like $N^2$.

```{table} Additions and comparisons for the diamond lattice with $N$ stages.
:name: tab-05-01-counts

| $N$ | Paths | Enumeration: additions | Enumeration: comparisons | DP: additions | DP: comparisons |
|---|---|---|---|---|---|
| 6 | 20 | 100 | 19 | 24 | 9 |
| 10 | 252 | 2 268 | 251 | 60 | 25 |
| 20 | 184 756 | 3 510 364 | 184 755 | 220 | 100 |
| 30 | 155 117 520 | 4 498 408 080 | 155 117 519 | 480 | 225 |
```

**General models.** In the general recurrence {eq}`eq-05-01-recurrence`, each
state of each stage needs one addition per decision and one comparison per
decision but one. The total work is therefore about

$$
\sum_{k=0}^{N-1}\ \sum_{y\in\mathcal Y_k}\abs{\mathcal D_k(y)}
\ \approx\ N\cdot(\text{states per stage})\cdot(\text{decisions per state}),
$$ (eq-05-01-work)

and the memory needed is one value and one policy entry per state. The work
grows only *linearly* in the number of stages, which is the great strength
of DP. Its weakness is the number of states per stage. If the state consists
of $r$ variables, each taking $q$ values, there are $q^{r}$ states per stage,
and the work grows exponentially in $r$. Bellman called this the **curse of
dimensionality**. The turn penalty doubled the state space; recording in
addition the number of turns made so far, as needed when the $j$th turn
costs $j$ ({doc}`08-06-examples-dp-fundamentals`), multiplies it by up to $N$.
Allocating two resources instead of one ({doc}`05-03-dp-resource-allocation`)
or tracking the age at the last overhaul of a machine
({doc}`05-02-dp-equipment-replacement`) adds a whole dimension. A good
formulation keeps the state as small as the principle of optimality allows:
this is why "too many arguments" in {prf:ref}`alg-05-01-dp-steps` is a real,
if less serious, error.

(sec-05-01-stochastic)=
## Deterministic and stochastic dynamic programming

All models in this part of the book are **deterministic**: a decision in a
known state leads to a known next state at a known cost. In a **stochastic**
DP the transition $t_k(y,d)$ is replaced by a probability distribution over
next states (and the cost may be random), the value function is the least
*expected* cost-to-go, and the recurrence becomes
$V_k(y)=\min_d\mathbb E\bigl[c_k(y,d,\xi)+V_{k+1}(t_k(y,d,\xi))\bigr]$ for a
random disturbance $\xi$. The principle of optimality still holds, for the
same three reasons (the state summarizes the past, the expected cost is
separable, the expectation is monotone), but the answer is no longer a fixed
sequence of decisions: it is a policy that prescribes a decision for each
state that might occur. This is why even in the deterministic case it pays
to record the whole policy table rather than only the optimal path.

## Summary and Common Mistakes

### Summary

- DP solves a multistage decision problem by solving one subproblem per
  state and combining the answers through a recurrence
  ({eq}`eq-05-01-recurrence`, {eq}`eq-05-01-boundary`).
- The principle of optimality ({prf:ref}`thm-05-01-principle`) holds when the
  state summarizes the past, the cost is separable and the combination is
  monotone; when it seems to fail, the state is too small.
- A DP model consists of stage, state, decision, transition, stage cost,
  value function, policy function, recurrence, boundary conditions and
  answer cell ({prf:ref}`def-05-01-dp-model`, {prf:ref}`alg-05-01-dp-steps`).
- Backward recursion ({eq}`eq-05-01-lattice-backward`, {numref}`tab-05-01-backward`)
  gives the best decision from every state to the end; forward recursion
  ({eq}`eq-05-01-lattice-forward`, {numref}`tab-05-01-forward`) gives the best
  path from the origin to every state. On the staged network both give
  13 min on $A$–$C$–$F$–$J$–$M$–$P$–$B$
  ({numref}`fig-05-01-diamond-network`).
- State augmentation adds the part of the past that the future depends on;
  with a turn penalty of 3 min the state is (vertex, direction) and the
  optimum is 17 min on $A$–$D$–$G$–$K$–$N$–$Q$–$B$ ({eq}`eq-05-01-turn`,
  {numref}`tab-05-01-turn`, {numref}`fig-05-01-turn-penalty-path`).
- DP work grows linearly with the number of stages but exponentially with
  the number of state variables, the curse of dimensionality
  ({eq}`eq-05-01-count-dp`, {eq}`eq-05-01-work`, {numref}`tab-05-01-counts`).

### Common mistakes

- **A state that forgets the past.** If the cost or the feasibility of a
  decision depends on earlier decisions (turns, coupons used, resources
  spent, time elapsed), that information must be part of the state.
  Otherwise the recurrence is wrong even though it can be evaluated.
- **Mixing conventions.** Defining $V$ with the direction of *arrival* and
  writing the recurrence with the direction of *leaving* (or the reverse)
  gives a model that may happen to work on one example and fail on another.
  Write the definition of $V$ in words first and derive the recurrence from
  it.
- **Missing boundary cases.** Arcs that do not exist need the value
  $\infty$; states at the last stage need explicit values; a recurrence that
  jumps two stages ahead needs care at the second-to-last stage.
- **Forgetting the answer cell.** State which value solves the original
  problem, for example a minimum over the augmented states at the origin.
- **Recording only the value.** Without the policy function the optimal
  decisions cannot be traced; keep $P_k(y)$ (all ties) for every state.
- **Confusing the two recursions.** $V_k(y)$ (backward) is the cost from
  $(k,y)$ to the end; $F_k(y)$ (forward) is the cost from the start to
  $(k,y)$. They are different functions with the same answer at the answer
  cells.

## References

- Bellman, R. (1957). *Dynamic programming*. Princeton University Press.
  (Ch. III, §3: the principle of optimality; Ch. I: multistage decision
  processes.)
- Dreyfus, S. E., & Law, A. M. (1977). *The art and theory of dynamic
  programming*. Academic Press. (Ch. 1: elementary path problems, the
  staged network, backward and forward recursion, turn penalties and
  operation counts.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Ch. 12: §12.1 recursive nature of DP computations, §12.2
  forward and backward recursion, §12.4 the problem of dimensionality.)
