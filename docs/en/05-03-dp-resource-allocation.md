# Resource Allocation and Capital Budgeting

## Overview

A limited resource (money, staff, equipment hours, water) has to be divided
among several activities, each of which turns its share into a return. When
the returns are arbitrary tables, with no linearity or concavity, linear and
nonlinear programming do not apply directly, but dynamic programming does:
the activities are taken one at a time, as if the allocation were a sequence
of decisions, and the amount of resource still unallocated is the state.
This chapter develops the **resource allocation** model with one resource,
solves it with stage tables and reads it as a longest path, contrasts
integer and continuous allocations, and shows that the knapsack problem is a
special case. With two resources the state becomes a pair, which shows how
quickly the work of DP grows with the number of state variables, and a
Lagrange multiplier can remove one of them. The chapter ends with **capital
budgeting**, where a single budget buys two priced resources, solved by two
DP formulations whose effort is compared. It builds on
{doc}`05-01-dp-fundamentals` (stages, states, value and policy functions,
recurrences and boundary conditions); the Lagrange multipliers of
{ref}`sec-05-03-lagrange` are introduced in
{doc}`02-03-convex-optimization-duality`, and only the idea is used.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- formulate the allocation of one resource among $N$ activities as a DP and
  solve it with stage tables, including all alternative optima;
- interpret the allocation DP as a longest path in a layered network and
  count the work it needs;
- distinguish integer from continuous allocations and solve simple
  continuous cases in closed form;
- write the knapsack problem as an allocation DP;
- formulate allocation problems with two resources, count their states, and
  use a Lagrange multiplier to reduce the state dimension; and
- formulate capital budgeting with a money state, reduce the stage decision
  to a single variable, and compare the effort of the two formulations.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $N$, $k$ | number of activities; activity (stage) index $k=1,\dots,N$ |
| $X$ | total units of the resource |
| $y$ | units not yet allocated when activity $k$ is considered (the state) |
| $x_k$, $x$ | units allocated to activity $k$ (the decision); $x$ when the stage is clear |
| $r_k(x)$ | return from allocating $x$ units to activity $k$ |
| $V_k(y)$, $P_k(y)$ | maximum return from activities $k,\dots,N$ with $y$ units left; optimal allocation(s) to activity $k$ there |
| $a_k$, $h_k$ | size and worth of item $k$ in the knapsack problem |
| $X_1$, $X_2$ | totals of resources 1 and 2; $(y_1,y_2)$ units left, $(x_1,x_2)$ units allocated at a stage |
| $\lambda$, $q(\lambda)$ | Lagrange multiplier (price) of resource 2; Lagrange dual function (an upper bound) |
| $B$, $b$ | budget and dollars left (capital budgeting) |
| $c_1$, $c_2$ | unit prices of resources 1 and 2 |
| $\hat r_k(\beta)$ | best return of activity $k$ from at most $\beta$ dollars |

(sec-05-03-problem)=
## The allocation problem

There are $X$ units of a resource and $N$ activities. Allocating $x$ units to
activity $k$ returns $r_k(x)$, a function given by a table. The units are
indivisible, so $x_k\in\{0,1,\dots,X\}$. The problem is

$$
\begin{aligned}
\max_{x_1,\dots,x_N}\ & \sum_{k=1}^{N}r_k(x_k)\\
\st\ & \sum_{k=1}^{N}x_k=X,\\
& x_k\in\{0,1,\dots,X\},\quad k=1,\dots,N.
\end{aligned}
$$ (eq-05-03-model)

The returns of different activities are independent of one another: the
return of activity $k$ depends on $x_k$ only. This **separability** is what
the DP needs; the shape of the $r_k$ does not matter. If unused units are
allowed, the constraint becomes $\sum_k x_k\le X$; when every $r_k$ is
nondecreasing, the two versions have the same optimal value, since leftover
units can always be given to some activity without loss.

Nothing in {eq}`eq-05-03-model` happens in time. The DP view is a device:
imagine that the activities are served one after another. First an amount
$x_1$ is given to activity 1; then $x_2$ of the remaining $X-x_1$ units goes
to activity 2; and so on. The activities become the **stages**, and when
activity $k$ is served, the only thing that matters about the earlier
decisions is how much resource they left.

(sec-05-03-one-resource)=
## The one-resource DP

The elements of {prf:ref}`def-05-01-dp-model` are as follows.

**Stages.** Activity $k$, $k=1,\dots,N$.

**States.** The number of units $y\in\{0,1,\dots,X\}$ not yet allocated
when activity $k$ is served. At stage 1 the state is $X$.

**Decisions and transitions.** Allocate $x\in\{0,1,\dots,y\}$ units to
activity $k$; the return is $r_k(x)$ and the next state is $y-x$.

:::{prf:definition} Value function of the allocation problem
:label: def-05-03-value

$V_k(y)$ is the maximum total return obtainable from activities
$k,k+1,\dots,N$ when $y$ units remain to be allocated to them. $P_k(y)$ is
the set of allocations $x$ that attain it.
:::

By the principle of optimality, whatever is given to activity $k$, the
remaining $y-x$ units must be allocated optimally among activities
$k+1,\dots,N$:

$$
V_k(y)=\max_{x=0,1,\dots,y}\bigl\{r_k(x)+V_{k+1}(y-x)\bigr\},\qquad y=0,1,\dots,X.
$$ (eq-05-03-recurrence)

With the equality constraint of {eq}`eq-05-03-model`, the last activity must
take whatever is left:

$$
V_N(y)=r_N(y),\qquad y=0,1,\dots,X.
$$ (eq-05-03-boundary)

(With the inequality version one writes instead $V_{N+1}(y)=0$ for all $y$
and applies {eq}`eq-05-03-recurrence` also at $k=N$.) The answer cell is $V_1(X)$; at
stage 1 only the state $y=X$ is needed. The allocation is recovered forward:
$x_1\in P_1(X)$, then $x_2\in P_2(X-x_1)$, and so on.

:::{prf:example} Eight units among four activities
:label: ex-05-03-one-resource

Eight units of a resource are to be allocated among four activities with
the returns of {numref}`tab-05-03-returns`. Find an allocation of maximum
total return.

```{list-table} Returns $r_k(x)$ of the four activities
:name: tab-05-03-returns
:header-rows: 1
:stub-columns: 1

* - $x$
  - 0
  - 1
  - 2
  - 3
  - 4
  - 5
  - 6
  - 7
  - 8
* - $r_1(x)$
  - 0
  - 3
  - 7
  - 10
  - 12
  - 13
  - 14
  - 14
  - 14
* - $r_2(x)$
  - 0
  - 1
  - 2
  - 4
  - 8
  - 13
  - 17
  - 19
  - 20
* - $r_3(x)$
  - 0
  - 2
  - 4
  - 6
  - 8
  - 10
  - 12
  - 14
  - 16
* - $r_4(x)$
  - 0
  - 1
  - 3
  - 6
  - 9
  - 12
  - 14
  - 16
  - 17
```
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

```{figure} figures/05-03-return-functions.svg
:name: fig-05-03-returns
:width: 520px
:class: dark-light
:alt: Return functions of the four activities as connected points, with the optimal allocation marked

Return functions $r_k(x)$ of {prf:ref}`ex-05-03-one-resource`. Activity 1 has
mostly decreasing marginal returns, activity 3 constant ones, and
activities 2 and 4 increasing and then decreasing ones. The
rings mark the optimal allocation $x=(2,6,0,0)$, worth $7+17=24$.
```

{numref}`fig-05-03-returns` plots the returns. They are not concave (activity
2 returns little for its first units and a lot for units 4 to 6), so greedy
allocation by marginal return would fail; the DP handles any shape. By {eq}`eq-05-03-boundary`, $V_4=r_4$. The
recurrence at stage 3 compares $y+1$ options for each state, for example

$$
\begin{aligned}
V_3(3)&=\max\{0+V_4(3),\ 2+V_4(2),\ 4+V_4(1),\ 6+V_4(0)\}\\
&=\max\{6,\ 5,\ 5,\ 6\}=6\qquad(x=0\text{ or }3),\\
V_3(8)&=\max\{0+17,\ 2+16,\ 4+14,\ 6+12,\ 8+9,\ 10+6,\ 12+3,\ 14+1,\ 16+0\}\\
&=18\qquad(x=1,2\text{ or }3).
\end{aligned}
$$

{numref}`tab-05-03-one-resource` lists all values and maximizers. At stage 1
only $y=8$ is needed:

$$
\begin{aligned}
V_1(8)&=\max\{0+21,\ 3+19,\ 7+17,\ 10+13,\ 12+9,\ 13+6,\ 14+4,\ 14+2,\ 14+0\}\\
&=\max\{21,22,24,23,21,19,18,16,14\}=24\qquad(x=2).
\end{aligned}
$$

```{list-table} Stage table of the one-resource example
:name: tab-05-03-one-resource
:header-rows: 1
:stub-columns: 1

* - $y$
  - 0
  - 1
  - 2
  - 3
  - 4
  - 5
  - 6
  - 7
  - 8
* - $V_4(y)=r_4(y)$
  - 0
  - 1
  - 3
  - 6
  - 9
  - 12
  - 14
  - 16
  - 17
* - $V_3(y)$
  - 0
  - 2
  - 4
  - 6
  - 9
  - 12
  - 14
  - 16
  - 18
* - $P_3(y)$
  - 0
  - 1
  - 2
  - 0, 3
  - 0
  - 0
  - 0, 1
  - 0, 1, 2
  - 1, 2, 3
* - $V_2(y)$
  - 0
  - 2
  - 4
  - 6
  - 9
  - 13
  - 17
  - 19
  - 21
* - $P_2(y)$
  - 0
  - 0
  - 0
  - 0
  - 0
  - 5
  - 6
  - 6, 7
  - 6, 7
* - $V_1(y)$
  -
  -
  -
  -
  -
  -
  -
  -
  - 24
* - $P_1(y)$
  -
  -
  -
  -
  -
  -
  -
  -
  - 2
```

The maximum return is $V_1(8)=24$. Forward: $x_1=2$ leaves $y=6$ at stage 2,
where $P_2(6)=6$, so $x_2=6$; then $y=0$ and $x_3=x_4=0$. The optimal
allocation is $x=(2,6,0,0)$ with return $r_1(2)+r_2(6)=7+17=24$. Although
the stage tables contain many ties, the path from $y=8$ meets none of them,
so this allocation is the unique optimum; enumerating the 165 allocations
of 8 units confirms it.
:::

**Effort.** Each stage has $X+1$ states, and state $y$ compares $y+1$
options, so a full stage costs

$$
\sum_{y=0}^{X}(y+1)=\frac{(X+1)(X+2)}{2}
$$ (eq-05-03-stage-count)

additions and slightly fewer comparisons. Stages $2,\dots,N-1$ are full, the
boundary stage costs nothing and stage 1 needs only $X+1$ options, so the
total is about $(N-2)X^2/2$: 99 additions in the example. Enumeration must
consider all $\binom{X+N-1}{N-1}$ allocations, which grows like $X^{N-1}$;
the DP is polynomial in both $N$ and $X$.

**Longest-path view.** Draw a node for every pair (stage $k$, units left $y$)
and an arc from $(k,y)$ to $(k+1,y-x)$ of length $r_k(x)$ for every
$x\le y$, plus a terminal arc of length $r_N(y)$ out of each node of stage
$N$. Then every allocation is a path from $(1,X)$ to the terminal node, its
return is the path length, and the problem is a **longest path** in an
acyclic layered network; {eq}`eq-05-03-recurrence` is the acyclic
longest-path recursion. The network for this example with $X=5$ is drawn and
solved in {prf:ref}`ex-08-07-l9-longest-path`.

(sec-05-03-continuous)=
## Integer and continuous allocations

The DP above assumes integer units. When the resource is divisible, the
decision $x$ ranges over the interval $[0,y]$ and the recurrence becomes

$$
V_k(y)=\max_{0\le x\le y}\bigl\{r_k(x)+V_{k+1}(y-x)\bigr\},\qquad 0\le y\le X,
$$ (eq-05-03-continuous)

an equation between functions rather than tables. Two ways to use it:

- **Analytically**, when the $r_k$ are simple. For example, with
  $r_k(x)=\alpha_k\sqrt x$ and $\alpha_k>0$, the last stage gives
  $V_N(y)=\alpha_N\sqrt y$. Suppose $V_{k+1}(y)=\gamma\sqrt y$ with $\gamma>0$.
  The function $\alpha_k\sqrt x+\gamma\sqrt{y-x}$ is concave in $x$, and its
  derivative vanishes where $\alpha_k^2(y-x)=\gamma^2x$, that is at
  $x=y\alpha_k^2/(\alpha_k^2+\gamma^2)$, with value
  $\sqrt{\alpha_k^2+\gamma^2}\,\sqrt y$. By induction
  $V_k(y)=\bigl(\sum_{l=k}^{N}\alpha_l^2\bigr)^{1/2}\sqrt y$, and the optimal
  allocation is proportional to the squared coefficients:
  $x_k=X\alpha_k^2/\sum_l\alpha_l^2$. With $\alpha=(1,2,3)$ and $X=14$ this
  gives $x=(1,4,9)$ and a total return of $\sqrt{14\cdot14}=14$; a numerical
  solver returns the same allocation.
- **Numerically**, by restricting $y$ and $x$ to a grid of step $\Delta$ and
  solving the integer problem in units of $\Delta$. The error shrinks with
  $\Delta$, at the price of more states: halving $\Delta$ doubles the states
  and multiplies the work per stage by about four.

When every $r_k$ is concave and differentiable, the continuous problem is a
convex program and the KKT conditions ({prf:ref}`thm-06-03-convex-kkt`)
apply: at an optimum the marginal returns $r_k'(x_k)$ of all activities with
$x_k>0$ are equal (to a multiplier $\lambda$), and no activity with $x_k=0$
has a larger marginal return at zero. The square-root example satisfies
this: $\alpha_k/(2\sqrt{x_k})=\sqrt{\sum_l\alpha_l^2}/(2\sqrt X)$ for every
$k$. The DP does not need concavity, which is its advantage when returns
have fixed set-up effects, thresholds or tables like
{numref}`tab-05-03-returns`. An integer solution cannot, in general, be found by
rounding a continuous one. Other objectives fit the same scheme: a continuous
min–max allocation is solved stage by stage in {prf:ref}`ex-08-07-a8-q5`.

(sec-05-03-knapsack)=
## The knapsack problem

A hiker can carry a load of at most $X$ kg. Item $k$ weighs $a_k$ kg and is
worth $h_k$; each item is taken at most once (the **0–1 knapsack problem**)
or in any integer number of copies (the **integer knapsack problem**). Both
are allocation problems in which the decision at stage $k$ is the number of
copies $x_k$ of item $k$, the return is $h_kx_k$, and the resource used is
$a_kx_k$:

$$
\begin{aligned}
V_k(y)&=\max_{x=0,1,\dots,\lfloor y/a_k\rfloor}\bigl\{h_kx+V_{k+1}(y-a_kx)\bigr\},\\
V_{N+1}(y)&=0,
\end{aligned}
$$ (eq-05-03-knapsack)

with $x\in\{0,1\}$ only (and $a_k\le y$ for $x=1$) in the 0–1 case. In the
language of {eq}`eq-05-03-recurrence`, $u$ units of capacity given to item $k$
return $r_k(u)=h_k\lfloor u/a_k\rfloor$ in the integer case, and $h_k$ if
$u\ge a_k$ (0 otherwise) in the 0–1 case; the inequality version of the
boundary condition is used because capacity may be left unused. A 0–1 knapsack with five housing applications is solved in
{prf:ref}`ex-08-07-a8-q4`, and a course-selection problem with a lower bound on
every allocation in {prf:ref}`ex-08-07-a8-q3`.

(sec-05-03-two-resources)=
## Two resources

Now each activity uses two resources: there are $X_1$ units of resource 1 and
$X_2$ units of resource 2, and allocating $x_1$ units of resource 1 and $x_2$
units of resource 2 to activity $k$ returns $r_k(x_1,x_2)$. The integer
program is

$$
\begin{aligned}
\max\ & \sum_{k=1}^{N}r_k(x_{k1},x_{k2})\\
\st\ & \sum_{k=1}^{N}x_{k1}=X_1,\\
& \sum_{k=1}^{N}x_{k2}=X_2,\\
& x_{k1}\in\{0,\dots,X_1\},\quad k=1,\dots,N,\\
& x_{k2}\in\{0,\dots,X_2\},\quad k=1,\dots,N.
\end{aligned}
$$ (eq-05-03-two-model)

To decide what to give activity $k$ one must know how much of **each**
resource is left, so the state is the pair $(y_1,y_2)$ and

$$
V_k(y_1,y_2)=\max_{\substack{x_1=0,\dots,y_1\\ x_2=0,\dots,y_2}}
\bigl\{r_k(x_1,x_2)+V_{k+1}(y_1-x_1,\,y_2-x_2)\bigr\},
$$ (eq-05-03-two-resources)

$$
V_N(y_1,y_2)=r_N(y_1,y_2),
$$ (eq-05-03-two-boundary)

for $0\le y_1\le X_1$, $0\le y_2\le X_2$. The answer cell is $V_1(X_1,X_2)$.

:::{prf:example} Two resources among three activities
:label: ex-05-03-two-resources

Three units of each of two resources ($X_1=X_2=3$) are to be allocated among
three activities with the returns of {numref}`tab-05-03-two-returns` (rows: units
$x_1$ of resource 1; columns: units $x_2$ of resource 2). Find an allocation
of maximum total return.

```{list-table} Returns $r_k(x_1,x_2)$ of the three activities
:name: tab-05-03-two-returns
:header-rows: 2
:stub-columns: 1

* -
  - $r_1$
  -
  -
  -
  - $r_2$
  -
  -
  -
  - $r_3$
  -
  -
  -
* - $x_1\backslash x_2$
  - 0
  - 1
  - 2
  - 3
  - 0
  - 1
  - 2
  - 3
  - 0
  - 1
  - 2
  - 3
* - 0
  - 0
  - 1
  - 3
  - 6
  - 0
  - 2
  - 4
  - 6
  - 0
  - 3
  - 5
  - 8
* - 1
  - 4
  - 5
  - 6
  - 7
  - 1
  - 4
  - 6
  - 7
  - 2
  - 5
  - 7
  - 9
* - 2
  - 5
  - 6
  - 7
  - 8
  - 4
  - 6
  - 8
  - 9
  - 4
  - 7
  - 9
  - 11
* - 3
  - 6
  - 7
  - 8
  - 9
  - 6
  - 8
  - 10
  - 11
  - 6
  - 9
  - 11
  - 13
```
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The boundary is $V_3=r_3$. At stage 2 each of the 16 states compares all
pairs $(x_1,x_2)\le(y_1,y_2)$; for example, writing the options as
$(x_1,x_2)\!:\ r_2(x_1,x_2)+V_3(\text{rest})$,

$$
\begin{aligned}
V_2(1,1)&=\max\{(0,0)\!:0+5,\ (0,1)\!:2+2,\ (1,0)\!:1+3,\ (1,1)\!:4+0\}=5,\\
V_2(3,1)&=\max\{(0,0)\!:0+9,\ (0,1)\!:2+6,\ (1,0)\!:1+7,\ (1,1)\!:4+4,\\
&\qquad\quad (2,0)\!:4+5,\ (2,1)\!:6+2,\ (3,0)\!:6+3,\ (3,1)\!:8+0\}=9.
\end{aligned}
$$

{numref}`tab-05-03-two-stage` gives the whole of $V_2$ and one optimal decision
per state.

```{list-table} Stage 2 of the two-resource example: $V_2(y_1,y_2)$ and an optimal decision $P_2(y_1,y_2)$
:name: tab-05-03-two-stage
:header-rows: 1
:stub-columns: 1

* - $y_1\backslash y_2$
  - 0
  - 1
  - 2
  - 3
* - 0
  - 0 $\;(0,0)$
  - 3 $\;(0,0)$
  - 5 $\;(0,0)$
  - 8 $\;(0,0)$
* - 1
  - 2 $\;(0,0)$
  - 5 $\;(0,0)$
  - 7 $\;(0,0)$
  - 9 $\;(0,0)$
* - 2
  - 4 $\;(0,0)$
  - 7 $\;(0,0)$
  - 9 $\;(0,0)$
  - 12 $\;(2,0)$
* - 3
  - 6 $\;(0,0)$
  - 9 $\;(0,0)$
  - 11 $\;(0,0)$
  - 14 $\;(3,0)$
```

Many states have ties: the decision is unique only in states $(0,0)$,
$(0,1)$, $(0,3)$, $(1,0)$, $(1,1)$, $(2,3)$ and $(3,3)$; in $(3,2)$, for
instance, seven decisions give 11. At stage 1 only $(3,3)$ is needed:

$$
\begin{aligned}
V_1(3,3)&=\max_{(x_1,x_2)}\{r_1(x_1,x_2)+V_2(3-x_1,3-x_2)\}\\
&=r_1(1,0)+V_2(2,3)=4+12=16,
\end{aligned}
$$

where the 16 options are, row by row in $x_1=0,1,2,3$ and $x_2=0,1,2,3$:
$14,12,12,12$; $16,14,13,11$; $14,13,12,10$; $14,12,11,9$. The maximizer
$(1,0)$ is unique. Forward: activity 1 gets $(1,0)$, leaving $(2,3)$;
$P_2(2,3)=(2,0)$ is unique, leaving $(0,3)$ for activity 3. The optimal
allocation is

$$
\begin{aligned}
(x_{11},x_{12})&=(1,0),\\
(x_{21},x_{22})&=(2,0),\\
(x_{31},x_{32})&=(0,3),
\end{aligned}
$$

with return $4+4+8=16$. Enumeration over all allocations confirms that it is
the unique optimum, also when unused units are allowed.
:::

**State growth.** With one resource each stage had $X+1$ states and
$(X+1)(X+2)/2$ options in total. With two resources it has $(X_1+1)(X_2+1)$
states, and by the same count the options number

$$
\sum_{y_1=0}^{X_1}\sum_{y_2=0}^{X_2}(y_1+1)(y_2+1)
=\frac{(X_1+1)(X_1+2)}{2}\cdot\frac{(X_2+1)(X_2+2)}{2},
$$ (eq-05-03-two-count)

which is 100 for the example. With $m$ resources of $X$ units each, a stage
has $(X+1)^m$ states and about $(X^2/2)^m$ options: the work grows
exponentially in the number of state variables. This is Bellman's **curse of
dimensionality**, and it limits DP to problems with few state variables
(three resources are discussed in {prf:ref}`ex-08-07-t2-q10`).

(sec-05-03-lagrange)=
## Removing a resource with a Lagrange multiplier

One way to fight the curse is to put a **price** $\lambda\ge0$ on resource 2
and drop its constraint. For fixed $\lambda$, consider

$$
q(\lambda)=\max\Bigl\{\sum_{k}\bigl[r_k(x_{k1},x_{k2})-\lambda x_{k2}\bigr]
:\ \sum_k x_{k1}\le X_1,\ x_{k2}\in\{0,\dots,X_2\}\Bigr\}+\lambda X_2.
$$ (eq-05-03-lagrange)

Every allocation feasible for the original problem is feasible here and its
objective is not smaller (the added term $\lambda(X_2-\sum_kx_{k2})$ is
nonnegative), so the **dual function** $q(\lambda)$ is an **upper bound** on
the optimal value for every $\lambda\ge0$, the analog for a maximization
of weak duality ({prf:ref}`thm-02-03-weak-duality`). The relaxed problem has
one resource only: each activity
first chooses its best $x_2$ for every $x_1$,

$$
\tilde r_k(x_1)=\max_{x_2=0,\dots,X_2}\bigl\{r_k(x_1,x_2)-\lambda x_2\bigr\},
$$

and the one-resource DP {eq}`eq-05-03-recurrence` is run on the $\tilde r_k$. If
the relaxed optimum happens to use exactly $X_2$ units of resource 2 (and all
$X_1$ units of resource 1), it is feasible for the original problem and its
value equals the bound, so it is
**optimal**. The multiplier is adjusted (raised if resource 2 is
over-used, lowered if under-used) until this happens or the bound is as small
as possible.

In {prf:ref}`ex-05-03-two-resources`, $\lambda=9/4$ gives
$\tilde r_1=(0,4,5,6)$, $\tilde r_2=(0,\tfrac74,4,6)$ and
$\tilde r_3=(\tfrac54,\tfrac{11}4,\tfrac{19}4,\tfrac{27}4)$ for
$x_1=0,1,2,3$. The one-resource DP gives $x_1=(1,2,0)$ with value
$4+4+\tfrac54=\tfrac{37}4$ and, behind $\tilde r_3(0)$, $x_2=(0,0,3)$: exactly
three units of resource 2. The bound is $\tfrac{37}4+3\cdot\tfrac94=16$, the
optimum found in the solution of the example. A poorly chosen price gives only a bound: at
$\lambda=3/2$ the relaxed solution uses five units of resource 2 and
$q=17$. Because the problem is discrete, there may be no price at which the
relaxed solution uses exactly $X_2$ (a **duality gap**); the multiplier then
gives a bound and a good heuristic, not a guarantee. The connection with LP
duality is developed in {doc}`02-03-convex-optimization-duality` and
{doc}`03-05-lp-duality`.

(sec-05-03-capital-budgeting)=
## Capital budgeting

In **capital budgeting** the activities do not receive the resources
directly; they receive money, with which resources are bought. Let the budget
be $B$ dollars, and let one unit of resource 1 cost $c_1$ and one unit of
resource 2 cost $c_2$. If activity $k$ buys $x_{k1}$ and $x_{k2}$ units, it
returns $r_k(x_{k1},x_{k2})$:

$$
\begin{aligned}
\max\ & \sum_{k=1}^{N}r_k(x_{k1},x_{k2})\\
\st\ & c_1\sum_{k=1}^{N}x_{k1}+c_2\sum_{k=1}^{N}x_{k2}\le B,\\
& x_{k1},x_{k2}\ \text{nonnegative integers},\quad k=1,\dots,N.
\end{aligned}
$$ (eq-05-03-budget-model)

The two resources are now linked through one budget constraint, so a single
state variable, the money left, suffices.

**Formulation 1: money as the state, two decisions per stage.** With $b$
dollars left when activity $k$ is served, any pair with
$c_1x_1+c_2x_2\le b$ can be bought:

$$
\begin{aligned}
V_k(b)&=\max_{\substack{x_1,x_2\ \text{integer}\ \ge0\\ c_1x_1+c_2x_2\le b}}
\bigl\{r_k(x_1,x_2)+V_{k+1}(b-c_1x_1-c_2x_2)\bigr\},\\
V_{N+1}(b)&=0,
\end{aligned}
$$ (eq-05-03-budget-1)

for $b=0,1,\dots,B$ (prices and budget in whole dollars). The answer cell is
$V_1(B)$. The joint condition $c_1x_1+c_2x_2\le b$ belongs in the maximum;
the separate bounds $x_1\le\lfloor b/c_1\rfloor$ and $x_2\le\lfloor b/c_2\rfloor$
alone would allow overspending.

**Formulation 2: collapse each stage to one decision.** The pair $(x_1,x_2)$
matters to later stages only through the money it costs. So first solve, for
every activity and every amount $\beta=0,1,\dots,B$, the small problem

$$
\hat r_k(\beta)=\max\bigl\{r_k(x_1,x_2):\ c_1x_1+c_2x_2\le\beta\bigr\},
$$ (eq-05-03-hat-r)

the best return of activity $k$ from at most $\beta$ dollars. Then allocate
money as a single resource:

$$
\begin{aligned}
V_k(b)&=\max_{\beta=0,1,\dots,b}\bigl\{\hat r_k(\beta)+V_{k+1}(b-\beta)\bigr\},\\
V_{N+1}(b)&=0.
\end{aligned}
$$ (eq-05-03-budget-2)

The two formulations give the same value function: in {eq}`eq-05-03-budget-1`,
group the pairs by the amount $\beta$ they cost at most; for a fixed $\beta$
the best is $\hat r_k(\beta)$, and $V_{k+1}$ is nonincreasing in the money
spent, so nothing is lost. Formulation 2 is the one-resource DP of
{ref}`sec-05-03-one-resource` applied to money, after a preprocessing step.

:::{prf:example} A six-dollar budget for four activities
:label: ex-05-03-capital-budgeting

A budget of $B=\$6$ is available. Resource 1 costs $c_1=\$3$ per unit and
resource 2 costs $c_2=\$1$ per unit. The returns of the four activities are
given in {numref}`tab-05-03-budget-returns` for every affordable pair
($3x_1+x_2\le6$). Find all optimal plans.

```{list-table} Returns $r_k(x_1,x_2)$ of the capital-budgeting example (pairs with $3x_1+x_2\le6$)
:name: tab-05-03-budget-returns
:header-rows: 1
:stub-columns: 2

* - $k$
  - $x_1$
  - $x_2=0$
  - 1
  - 2
  - 3
  - 4
  - 5
  - 6
* - 1
  - 0
  - 0
  - 5
  - 9
  - 11
  - 12
  - 13
  - 14
* - 1
  - 1
  - 5
  - 9
  - 13
  - 18
  -
  -
  -
* - 1
  - 2
  - 10
  -
  -
  -
  -
  -
  -
* - 2
  - 0
  - 0
  - 3
  - 5
  - 6
  - 7
  - 8
  - 9
* - 2
  - 1
  - 2
  - 5
  - 9
  - 12
  -
  -
  -
* - 2
  - 2
  - 7
  -
  -
  -
  -
  -
  -
* - 3
  - 0
  - 0
  - 1
  - 3
  - 6
  - 10
  - 12
  - 13
* - 3
  - 1
  - 4
  - 6
  - 8
  - 11
  -
  -
  -
* - 3
  - 2
  - 11
  -
  -
  -
  -
  -
  -
* - 4
  - 0
  - 0
  - 1
  - 3
  - 6
  - 10
  - 13
  - 14
* - 4
  - 1
  - 7
  - 10
  - 14
  - 17
  -
  -
  -
* - 4
  - 2
  - 13
  -
  -
  -
  -
  -
  -
```

*(The tables are filled exactly for the pairs with $3x_1+x_2\le6$, that is,
for a price of \$1 per unit of resource 2. The same tables with a price of
\$2 define a different problem, solved in {prf:ref}`ex-08-07-l9-budget-dearer`;
its optimum is 13.)*
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Preprocessing.** For activity 1 and $\beta=6$, the affordable pairs are
$(0,0),\dots,(0,6)$, $(1,0),\dots,(1,3)$ and $(2,0)$; the best is $r_1(1,3)=18$.
For $\beta=5$ the best is $13$, reached by $(0,5)$ and by $(1,2)$. Doing this
for all activities and amounts gives {numref}`tab-05-03-hat-r`. Because these
return tables increase in $x_2$, the best pair for $\beta$ always spends the
whole amount: $\hat r_k(\beta)=\max_{x_1}r_k(x_1,\beta-3x_1)$.

```{list-table} Best return $\hat r_k(\beta)$ from at most $\beta$ dollars, with the pairs $(x_1,x_2)$ that reach it
:name: tab-05-03-hat-r
:header-rows: 1

* - $\beta$
  - $\hat r_1$
  - $\hat r_2$
  - $\hat r_3$
  - $\hat r_4$
* - 0
  - 0 $\;(0,0)$
  - 0 $\;(0,0)$
  - 0 $\;(0,0)$
  - 0 $\;(0,0)$
* - 1
  - 5 $\;(0,1)$
  - 3 $\;(0,1)$
  - 1 $\;(0,1)$
  - 1 $\;(0,1)$
* - 2
  - 9 $\;(0,2)$
  - 5 $\;(0,2)$
  - 3 $\;(0,2)$
  - 3 $\;(0,2)$
* - 3
  - 11 $\;(0,3)$
  - 6 $\;(0,3)$
  - 6 $\;(0,3)$
  - 7 $\;(1,0)$
* - 4
  - 12 $\;(0,4)$
  - 7 $\;(0,4)$
  - 10 $\;(0,4)$
  - 10 $\;(0,4),(1,1)$
* - 5
  - 13 $\;(0,5),(1,2)$
  - 9 $\;(1,2)$
  - 12 $\;(0,5)$
  - 14 $\;(1,2)$
* - 6
  - 18 $\;(1,3)$
  - 12 $\;(1,3)$
  - 13 $\;(0,6)$
  - 17 $\;(1,3)$
```

**Money allocation.** Now {eq}`eq-05-03-budget-2` runs backwards from
$V_5\equiv0$, so $V_4=\hat r_4$. At stage 3, for example,

$$
\begin{aligned}
V_3(4)&=\max\{\hat r_3(0)+V_4(4),\ \hat r_3(1)+V_4(3),\ \hat r_3(2)+V_4(2),\ \hat r_3(3)+V_4(1),\ \hat r_3(4)+V_4(0)\}\\
&=\max\{10,\ 8,\ 6,\ 7,\ 10\}=10\qquad(\beta=0\text{ or }4).
\end{aligned}
$$

*(The option $\beta=3$ here is $6+V_4(1)=6+1=7$; a value $6+3$ sometimes
printed for it is a slip that does not change the result.)* All values are
in {numref}`tab-05-03-budget-stages`, where the decision is the amount $\beta$
given to activity $k$.

```{list-table} Stage table of the money allocation (values in dollars of return; decision $\beta$ in dollars)
:name: tab-05-03-budget-stages
:header-rows: 1

* - $b$
  - $V_4(b)$
  - $V_3(b)$
  - $\beta_3$
  - $V_2(b)$
  - $\beta_2$
  - $V_1(b)$
  - $\beta_1$
* - 0
  - 0
  - 0
  - 0
  - 0
  - 0
  -
  -
* - 1
  - 1
  - 1
  - 0, 1
  - 3
  - 1
  -
  -
* - 2
  - 3
  - 3
  - 0, 2
  - 5
  - 2
  -
  -
* - 3
  - 7
  - 7
  - 0
  - 7
  - 0
  -
  -
* - 4
  - 10
  - 10
  - 0, 4
  - 10
  - 0, 1
  -
  -
* - 5
  - 14
  - 14
  - 0
  - 14
  - 0
  -
  -
* - 6
  - 17
  - 17
  - 0
  - 17
  - 0, 1
  - 19
  - 1, 2
```

The maximum return is $V_1(6)=19$. Following the ties forward gives every
optimal plan:

- $\beta_1=1$: activity 1 buys $(0,1)$; then $b=5$, $\beta_2=0$, $\beta_3=0$,
  and activity 4 gets \$5 and buys $(1,2)$. Return $5+14=19$.
- $\beta_1=2$: activity 1 buys $(0,2)$; then $b=4$, where $\beta_2\in\{0,1\}$.
  - $\beta_2=1$: activity 2 buys $(0,1)$, $b=3$, $\beta_3=0$, activity 4
    buys $(1,0)$. Return $9+3+7=19$.
  - $\beta_2=0$: $b=4$ at stage 3, where $\beta_3\in\{0,4\}$. With $\beta_3=4$
    activity 3 buys $(0,4)$: return $9+10=19$. With $\beta_3=0$ activity 4
    gets \$4 and buys $(0,4)$ or $(1,1)$: return $9+10=19$ in both cases.

There are thus **five** optimal plans (listed as the pairs of activities
1–4):

| Plan | Activity 1 | Activity 2 | Activity 3 | Activity 4 | Return |
|---|---|---|---|---|---|
| 1 | $(0,1)$ | $(0,0)$ | $(0,0)$ | $(1,2)$ | 19 |
| 2 | $(0,2)$ | $(0,1)$ | $(0,0)$ | $(1,0)$ | 19 |
| 3 | $(0,2)$ | $(0,0)$ | $(0,4)$ | $(0,0)$ | 19 |
| 4 | $(0,2)$ | $(0,0)$ | $(0,0)$ | $(0,4)$ | 19 |
| 5 | $(0,2)$ | $(0,0)$ | $(0,0)$ | $(1,1)$ | 19 |

Every plan spends the whole \$6, and enumeration of all affordable
combinations finds exactly these five. Formulation 1 gives the same value
function $V_k(b)$ at every stage, as it must.
:::

**Effort of the two formulations.** Per stage, formulation 1 examines, for
every $b$, all pairs with $c_1x_1+c_2x_2\le b$; there are about
$b^2/(2c_1c_2)$ of them, so a stage costs about $B^3/(6c_1c_2)$ evaluations.
Formulation 2 examines each affordable pair once to build $\hat r_k$ (about
$B^2/(2c_1c_2)$ pairs, plus a running maximum over $\beta$) and then
$(B+1)(B+2)/2$ money splits. Its cost grows like $B^2$ instead of $B^3$.
{numref}`tab-05-03-budget-effort` shows the counts per stage for $c_1=3$,
$c_2=1$.

```{list-table} Evaluations per stage of the two capital-budgeting formulations ($c_1=\$3$, $c_2=\$1$)
:name: tab-05-03-budget-effort
:header-rows: 1

* - Budget $B$
  - Formulation 1
  - Formulation 2 (pairs + splits)
* - 6
  - 39
  - 12 + 28 = 40
* - 12
  - 185
  - 35 + 91 = 126
* - 30
  - 1 991
  - 176 + 496 = 672
* - 60
  - 13 881
  - 651 + 1 891 = 2 542
* - 120
  - 103 361
  - 2 501 + 7 381 = 9 882
```

For the small example the two are about equal; for realistic budgets the
two-step formulation is an order of magnitude cheaper. The same idea, solving
an inner problem once and reusing it, appears whenever a stage decision has
several components that affect the future only through one aggregate.

## Summary and Common Mistakes

### Summary

- Allocation of one resource among $N$ activities with separable returns is
  a DP with stages = activities and state = units left
  ({eq}`eq-05-03-recurrence`, {eq}`eq-05-03-boundary`); each stage costs
  $(X+1)(X+2)/2$ evaluations ({eq}`eq-05-03-stage-count`), and the problem is a
  longest path in a layered network.
- In {prf:ref}`ex-05-03-one-resource` eight units are best split as
  $(2,6,0,0)$, return 24 ({numref}`tab-05-03-one-resource`); the returns need
  not be concave ({numref}`fig-05-03-returns`).
- Continuous allocations satisfy {eq}`eq-05-03-continuous`; concave returns lead
  to equal marginal returns, and the knapsack problem is the allocation DP
  {eq}`eq-05-03-knapsack`.
- Two resources make the state a pair ({eq}`eq-05-03-two-resources`) and
  multiply the work ({eq}`eq-05-03-two-count`); the example's unique optimum is
  $(1,0),(2,0),(0,3)$ with return 16. A Lagrange multiplier on one resource
  ({eq}`eq-05-03-lagrange`) gives an upper bound and, at $\lambda=9/4$, a
  certificate of optimality.
- Capital budgeting has the money left as its state; collapsing each stage
  with $\hat r_k$ ({eq}`eq-05-03-hat-r`, {eq}`eq-05-03-budget-2`) turns it into a
  one-resource allocation and lowers the work from order $B^3$ to $B^2$ per
  stage. The six-dollar example has five optimal plans with return 19.

### Common mistakes

- **Wrong boundary condition for the version of the constraint.** With
  $\sum_kx_k=X$ use $V_N(y)=r_N(y)$; with $\sum_kx_k\le X$ use
  $V_{N+1}(y)=0$. Mixing them silently changes the problem when some $r_k$
  decreases.
- **Computing $V_1$ for every state.** Only $V_1(X)$ is needed; the earlier
  stages need all states.
- **Missing alternative optima.** Follow every tie on the forward pass, as in
  the capital-budgeting example; ties off the optimal path do not matter.
- **Separate bounds instead of a joint budget.** In capital budgeting the
  pair $(x_1,x_2)$ must satisfy $c_1x_1+c_2x_2\le b$ jointly.
- **Assuming the whole amount is spent.** Writing
  $\hat r_k(\beta)=\max_{x_1}r_k(x_1,\beta-c_1x_1)$ is valid only when the
  returns increase in $x_2$; definition {eq}`eq-05-03-hat-r` is always valid.
- **Rounding a continuous allocation** to get an integer one; solve the
  integer DP instead.

## References

- Bellman, R. (1957). *Dynamic programming*. Princeton University Press.
  (Preface and Chapter 1: the multistage allocation process and the curse of
  dimensionality.)
- Dreyfus, S. E., & Law, A. M. (1977). *The art and theory of dynamic
  programming*. Academic Press. (Chapter 3: resource allocation, including
  problems with several resources and the Lagrange multiplier procedure.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 11: dynamic programming,
  Section 11.3: distributing effort and resources among activities.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 12: deterministic dynamic programming, Section 12.3:
  the knapsack model and other applications.)
