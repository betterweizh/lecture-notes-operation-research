# Equipment Replacement

## Overview

Machines, vehicles and pumps grow more expensive to run as they age, and
their resale value falls. An owner who must keep one machine in service over
a planning horizon therefore faces a sequence of linked decisions: at the
start of every year, keep the current machine for another year or trade it
for a new one. Each decision changes the age of the machine that the next
decision inherits, which is exactly the structure that dynamic programming
(DP) exploits. This chapter builds the classical **equipment replacement**
model. It first takes the year as the stage and the age of the machine as
the state, which gives a two-argument value function $V_k(i)$; the model is
solved in full on a five-year example and read as a shortest path in an
acyclic network. It then derives a second, one-argument formulation that
uses the fact that every replacement brings the machine back to age 1, and
compares the computational effort of the two. The chapter closes with the
choice of boundary conditions at the end of the horizon, an extension in
which an overhaul is a third option and the state needs a second
coordinate, and a list of variants that are solved in
{doc}`08-07-examples-dp-applications`. It builds on the DP vocabulary of
{doc}`05-01-dp-fundamentals` (stages, states, value and policy functions,
recurrences and boundary conditions); the network view uses shortest paths in
acyclic networks from {doc}`04-02-shortest-paths` and
{doc}`05-04-dp-shortest-paths`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- formulate the equipment replacement problem as a DP with year as stage,
  age as state and keep/buy as decision, and state its recurrence and
  boundary conditions;
- solve the model with a stage table, recover all optimal policies, and
  check the result by enumeration or as a shortest path in an acyclic
  network;
- derive the one-argument formulation, explain why it is valid, and compare
  the number of states and operations of the two formulations;
- choose boundary conditions for different end-of-horizon rules (salvage,
  no salvage, a required final state) and discuss long-horizon policies; and
- extend the state when the cost depends on more of the history, as in the
  overhaul model, and count the resulting states.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $N$ | planning horizon in years; year $k=1,\dots,N$, and "year $N+1$" is the end of the horizon |
| $i$, $i_0$ | age of the machine (in years) at the start of a year; initial age at the start of year 1 |
| $c(i)$ | cost of operating for one year a machine that is of age $i$ at the start of that year |
| $p$ | price of a new machine (age 0) |
| $t(i)$ | trade-in value of a machine of age $i$ traded at the start of a year |
| $s(i)$ | salvage value at the end of year $N$ of a machine that has just turned age $i$ |
| $V_k(i)$, $P_k(i)$ | minimum cost of years $k,\dots,N$ starting year $k$ with a machine of age $i$; optimal decision there, $\mathrm{B}$ (buy) or $\mathrm{K}$ (keep) |
| $U_k$ | minimum cost of years $k,\dots,N$ starting year $k$ with a one-year-old machine (one-argument formulation), $U_k=V_k(1)$ |
| $j$ | age at the last overhaul ($j=0$: never overhauled), in {ref}`sec-05-02-overhaul` |
| $e$, $o$ | exchange cost and overhaul cost in the overhaul model |
| $A(m)$ | average cost per year of the stationary policy that trades every $m$ years |
| $\alpha$ | discount factor per year, $0<\alpha<1$ |

All money amounts in the examples are in thousands of dollars (\$1 000).

(sec-05-02-problem)=
## The replacement problem

A firm needs one machine in service in each of the next $N$ years. At the
start of year 1 it owns a machine of age $i_0$. At the start of every year
it either **keeps** the machine it has, or **buys** a new one and trades the
old one in. The data are:

- the operating cost $c(i)$ of a year in which the machine starts at age
  $i$; a new machine has age 0, so a machine bought at the start of a year
  costs $c(0)$ to run in that year;
- the price $p$ of a new machine and the trade-in value $t(i)$ of a machine
  of age $i\ge1$;
- the salvage value $s(i)$ received at the end of year $N$ for the machine
  then in hand, which has just turned age $i$.

The question is when, if ever, to replace the current machine, when to
replace its replacement, and so on, so that the total cost over the $N$
years, net of the final salvage value, is as small as possible. Revenue is
the same whatever the decisions (the machine is always in service), so it
does not enter the model.

The assumptions behind the model are worth stating. Costs and values depend
on the age only, not on the calendar year (this is relaxed in
{prf:ref}`ex-08-07-t2-q2`); decisions are made only at the start of a year; a
purchase and a trade-in happen together; and money is not discounted (see
{ref}`sec-05-02-horizon` for discounting).

(sec-05-02-two-argument)=
## Formulation with year and age

The DP model follows the steps of {prf:ref}`alg-05-01-dp-steps` and uses the
elements of {prf:ref}`def-05-01-dp-model`.

**Stages.** One decision is made at the start of each year, so year $k$ is
stage $k$, $k=1,\dots,N$. The end of the horizon is labeled stage $N+1$; no
decision is made there, but the age of the machine at that point fixes the
salvage value.

**States.** Everything that matters for the future is the age of the machine
in hand: the future operating costs, trade-in values and the salvage value
all depend on the age, and not on how the machine came to have it. The state
at stage $k$ is therefore the age $i$ at the start of year $k$. Which ages
can occur? If the original machine has never been replaced, it is of age
$i_0+k-1$. Otherwise the last purchase took place at the start of some year
$k'$ with $1\le k'\le k-1$, and the machine is now $k-k'$ years old, which
is one of $1,2,\dots,k-1$. The set of possible states at stage $k$ is

$$
\mathcal I_k=\{1,2,\dots,k-1\}\cup\{i_0+k-1\},
$$ (eq-05-02-states)

which has exactly $k$ elements (the two parts never overlap because
$i_0\ge1$). Stage 1 has the single state $i_0$, and the end of the horizon
has $N+1$ possible ages.

**Decisions and transitions.** In state $i$ the decision is $\mathrm{B}$ (buy a new
machine, trade in the old one) or $\mathrm{K}$ (keep). Buying costs $p-t(i)$ now plus
$c(0)$ for running the new machine, and the next state is age 1. Keeping
costs $c(i)$ and the next state is age $i+1$.

:::{prf:definition} Value and policy functions of the replacement model
:label: def-05-02-value

For $k=1,\dots,N+1$ and $i\in\mathcal I_k$, $V_k(i)$ is the minimum total
cost of years $k$ through $N$, net of the final salvage value, when year $k$
starts with a machine of age $i$. $P_k(i)\in\{\mathrm{B},\mathrm{K}\}$ is a decision that
attains the minimum.
:::

The principle of optimality gives the recurrence: whatever is decided in
year $k$, the remaining years must be run optimally from the state that
decision produces, so

$$
V_k(i)=\min\left\{
\begin{aligned}
&\mathrm{B}:\ p-t(i)+c(0)+V_{k+1}(1)\\
&\mathrm{K}:\ c(i)+V_{k+1}(i+1)
\end{aligned}
\right\},\qquad k=N,N-1,\dots,1,\ i\in\mathcal I_k.
$$ (eq-05-02-recurrence)

At the end of the horizon the machine is sold for its salvage value, a
negative cost:

$$
V_{N+1}(i)=-s(i),\qquad i\in\mathcal I_{N+1}.
$$ (eq-05-02-boundary)

The answer cell is $V_1(i_0)$, and an optimal replacement plan is read forward
from $P_1(i_0)$, following the transitions.

:::{prf:example} Five-year replacement plan
:label: ex-05-02-replacement

A firm plans for $N=5$ years. Its machine is $i_0=2$ years old at the start
of year 1, and a new machine costs $p=50$. The other data are in
{numref}`tab-05-02-data`.

```{list-table} Data of the replacement example (thousands of dollars)
:name: tab-05-02-data
:header-rows: 1
:stub-columns: 1

* - Age $i$
  - 0
  - 1
  - 2
  - 3
  - 4
  - 5
  - 6
  - 7
* - Operating cost $c(i)$
  - 10
  - 13
  - 20
  - 40
  - 70
  - 100
  - 100
  -
* - Trade-in value $t(i)$
  -
  - 32
  - 21
  - 11
  - 5
  - 0
  - 0
  -
* - Salvage value $s(i)$
  -
  - 25
  - 17
  - 8
  - 0
  - 0
  - 0
  - 0
```

Find a replacement plan of minimum total cost.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The cost of a buy decision depends only on the age traded in,
$p-t(i)+c(0)=60-t(i)$, which is 28, 39, 49, 55, 60 and 60 for ages 1 to 6.
By {eq}`eq-05-02-states` the possible ages are $\{2\}$ in year 1, $\{1,3\}$ in
year 2, $\{1,2,4\}$ in year 3, $\{1,2,3,5\}$ in year 4, $\{1,2,3,4,6\}$ in
year 5 and $\{1,2,3,4,5,7\}$ at the end. The boundary condition
{eq}`eq-05-02-boundary` gives

$$
\begin{aligned}
V_6(1)&=-25,\\
V_6(2)&=-17,\\
V_6(3)&=-8,\\
V_6(4)&=V_6(5)=V_6(7)=0.
\end{aligned}
$$

The recurrence {eq}`eq-05-02-recurrence` is then applied year by year,
backwards. For example, in year 5 with a machine of age 3,

$$
\begin{aligned}
V_5(3)&=\min\{\mathrm{B}:\ 49+V_6(1),\ \mathrm{K}:\ 40+V_6(4)\}\\
&=\min\{24,\ 40\}=24\quad(\mathrm{B}),
\end{aligned}
$$

and in year 2 with a machine of age 1 the two options tie:

$$
\begin{aligned}
V_2(1)&=\min\{\mathrm{B}:\ 28+V_3(1),\ \mathrm{K}:\ 13+V_3(2)\}\\
&=\min\{28+48,\ 13+63\}=76\quad(\mathrm{B}\text{ or }\mathrm{K}).
\end{aligned}
$$

{numref}`tab-05-02-two-arg` collects all fifteen evaluations.

```{list-table} Stage table of the two-argument formulation (thousands of dollars)
:name: tab-05-02-two-arg
:header-rows: 1

* - Year $k$
  - Age $i$
  - $\mathrm{B}$: $60-t(i)+V_{k+1}(1)$
  - $\mathrm{K}$: $c(i)+V_{k+1}(i+1)$
  - $V_k(i)$
  - $P_k(i)$
* - 5
  - 1
  - $28-25=3$
  - $13-17=-4$
  - $-4$
  - $\mathrm{K}$
* - 5
  - 2
  - $39-25=14$
  - $20-8=12$
  - 12
  - $\mathrm{K}$
* - 5
  - 3
  - $49-25=24$
  - $40+0=40$
  - 24
  - $\mathrm{B}$
* - 5
  - 4
  - $55-25=30$
  - $70+0=70$
  - 30
  - $\mathrm{B}$
* - 5
  - 6
  - $60-25=35$
  - $100+0=100$
  - 35
  - $\mathrm{B}$
* - 4
  - 1
  - $28-4=24$
  - $13+12=25$
  - 24
  - $\mathrm{B}$
* - 4
  - 2
  - $39-4=35$
  - $20+24=44$
  - 35
  - $\mathrm{B}$
* - 4
  - 3
  - $49-4=45$
  - $40+30=70$
  - 45
  - $\mathrm{B}$
* - 4
  - 5
  - $60-4=56$
  - $100+35=135$
  - 56
  - $\mathrm{B}$
* - 3
  - 1
  - $28+24=52$
  - $13+35=48$
  - 48
  - $\mathrm{K}$
* - 3
  - 2
  - $39+24=63$
  - $20+45=65$
  - 63
  - $\mathrm{B}$
* - 3
  - 4
  - $55+24=79$
  - $70+56=126$
  - 79
  - $\mathrm{B}$
* - 2
  - 1
  - $28+48=76$
  - $13+63=76$
  - 76
  - $\mathrm{B}$, $\mathrm{K}$
* - 2
  - 3
  - $49+48=97$
  - $40+79=119$
  - 97
  - $\mathrm{B}$
* - 1
  - 2
  - $39+76=115$
  - $20+97=117$
  - 115
  - $\mathrm{B}$
```

*(The keep option in year 5 at age 6 is $c(6)+V_6(7)=100+0=100$; a value of
70 sometimes printed for it is a slip and does not change the decision.)*

The minimum cost is $V_1(2)=115$, that is, \$115 000. The optimal plan is
read forward through the policy column. In year 1 the two-year-old machine
is replaced ($P_1(2)=\mathrm{B}$), so year 2 starts at age 1, where both decisions are
optimal:

- **$\mathrm{B}$ in year 2**: age 1 in year 3, $P_3(1)=\mathrm{K}$; age 2 in year 4,
  $P_4(2)=\mathrm{B}$; age 1 in year 5, $P_5(1)=\mathrm{K}$. The plan is $\mathrm{BBKBK}$.
- **$\mathrm{K}$ in year 2**: age 2 in year 3, $P_3(2)=\mathrm{B}$; age 1 in year 4,
  $P_4(1)=\mathrm{B}$; age 1 in year 5, $P_5(1)=\mathrm{K}$. The plan is $\mathrm{BKBBK}$.

Check of the first plan: $39+28+13+39+13-s(2)=132-17=115$. Both plans end
with a two-year-old machine. Enumerating all $2^5=32$ keep/buy sequences
confirms that these two are the only optimal plans; the next best,
$\mathrm{BKBKK}$, costs 116, and never replacing ($\mathrm{KKKKK}$) costs
$20+40+70+100+100=330$.
:::

Two features of the solution of {prf:ref}`ex-05-02-replacement` are typical. Near the end of the horizon
keeping becomes more attractive (in year 5 a young machine is kept) because
a replacement bought then would be sold after one year at a large loss
($p-s(1)=25$). And ties in the stage table produce alternative optimal
plans, which the policy function records as sets.

(sec-05-02-network)=
## Replacement as a shortest path

Every state $(k,i)$ of the model can be drawn as a node, and every decision
as an arc to the state it leads to: a keep arc $(k,i)\to(k+1,i+1)$ of length
$c(i)$ and a buy arc $(k,i)\to(k+1,1)$ of length $p-t(i)+c(0)$. A final arc
from each end node $(N+1,i)$ to a single terminal node carries the length
$-s(i)$. A replacement plan is then a path from $(1,i_0)$ to the terminal
node, its cost is the path length, and $V_k(i)$ is the length of a shortest
path from node $(k,i)$ to the terminal node.

```{figure} figures/05-02-replacement-network.svg
:name: fig-05-02-replacement-network
:width: 620px
:class: dark-light
:alt: Age-versus-year network of the replacement example with keep and buy arcs and the two optimal plans highlighted

Decision network of {prf:ref}`ex-05-02-replacement`. Nodes are (year, age)
pairs; rising arcs keep the machine, arcs down to age 1 buy a new one. Arc
lengths depend on the age only and are listed once per age on the left
(column K: $c(i)$, column B: $60-t(i)$); the numbers right of year 6 are the
terminal lengths $-s(i)$. The two optimal plans $\mathrm{BBKBK}$ and $\mathrm{BKBBK}$, of
length 115, are highlighted with their arc lengths; they split at year 2
and meet again at year 5.
```

{numref}`fig-05-02-replacement-network` shows the network of the example: 21
nodes and 30 decision arcs. The network is **acyclic** because every arc
goes from year $k$ to year $k+1$, so ordering the nodes by year is a
topological order, and the backward recursion {eq}`eq-05-02-recurrence` is the
shortest-path recursion {eq}`eq-05-04-acyclic-backward` for acyclic networks
run from the terminal node. A shortest-path code (see
{doc}`04-02-shortest-paths`) returns length 115 and the two highlighted paths.
The terminal arcs have negative lengths; this does no harm in an acyclic
network, but it rules out Dijkstra's algorithm unless a constant (here
25) is first added to all terminal arcs, which changes every path length by
the same amount.

The network view also explains the effort. Brute force evaluates
$2^N$ plans, each with $N$ additions: 160 additions and 31 comparisons for
$N=5$, and the count doubles with every extra year. The DP visits each node
once and each arc once.

(sec-05-02-one-argument)=
## A formulation with one argument

{numref}`fig-05-02-replacement-network` reveals a special structure: apart
from the nodes on the line of the original machine, every node lies on a
path that returned to age 1 at some earlier year. Each time a new machine
is bought, the process starts afresh from age 1, and what happens until
the next purchase is completely determined by how long the machine is kept.
This suggests a value function of the year alone. (The network of
{eq}`eq-04-02-replacement`, whose arcs are the periods between two purchases,
is built on the same idea.)

**Stages and states.** Let stage $k$ ($k=2,\dots,N+1$) be the start of a year
in which the firm holds a one-year-old machine; there is only one state per
stage.

**Decisions.** At such a stage the firm decides how long to keep the
machine: either until the end of the horizon, or for $j-k$ years
($j=k,k+1,\dots,N$), trading it at the start of year $j$, when it is of age
$j-k+1$, and running the new machine in year $j$. The next decision is then
taken at the start of year $j+1$ with a one-year-old machine.

:::{prf:definition} One-argument value function
:label: def-05-02-one-argument

For $k=2,\dots,N+1$, $U_k$ is the minimum total cost of years $k$ through
$N$, net of the final salvage value, when year $k$ starts with a
one-year-old machine. Hence $U_k=V_k(1)$.
:::

Keeping the machine throughout costs the operating costs of ages
$1,\dots,N+1-k$ less the salvage value at age $N+2-k$. Keeping it for
$j-k$ years and buying at the start of year $j$ costs the operating costs of
ages $1,\dots,j-k$, then $p-t(j-k+1)+c(0)$, then $U_{j+1}$. Taking the best
of these options,

$$
U_k=\min\left\{
\begin{aligned}
&\sum_{a=1}^{N+1-k}c(a)-s(N+2-k)\\
&\min_{j=k,\dots,N}\Bigl[\sum_{a=1}^{j-k}c(a)+p-t(j-k+1)+c(0)+U_{j+1}\Bigr]
\end{aligned}
\right\},\qquad k=N,\dots,2,
$$ (eq-05-02-one-arg)

with the boundary condition

$$
U_{N+1}=-s(1).
$$ (eq-05-02-one-arg-boundary)

The empty sum (for $j=k$) is zero: the machine is traded at once, at age 1.
The start of year 1 is handled separately because the original machine has
age $i_0$, not 1. Keeping it for $j-1$ years and trading at age $i_0+j-1$,

$$
V_1(i_0)=\min\left\{
\begin{aligned}
&\sum_{a=0}^{N-1}c(i_0+a)-s(i_0+N)\\
&\min_{j=1,\dots,N}\Bigl[\sum_{a=0}^{j-2}c(i_0+a)+p-t(i_0+j-1)+c(0)+U_{j+1}\Bigr]
\end{aligned}
\right\}.
$$ (eq-05-02-one-arg-start)

:::{prf:example} The five-year plan with one state per year
:label: ex-05-02-one-argument

Solve {prf:ref}`ex-05-02-replacement` with {eq}`eq-05-02-one-arg` to
{eq}`eq-05-02-one-arg-start`.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The boundary value is $U_6=-s(1)=-25$. The options of every stage are
listed in {numref}`tab-05-02-one-arg`; "B$j$" means keep the one-year-old
machine until the start of year $j$ and buy there.

```{list-table} Options of the one-argument formulation (thousands of dollars)
:name: tab-05-02-one-arg
:header-rows: 1

* - Stage
  - Keep to the end
  - B$k$
  - B$(k+1)$
  - B$(k+2)$
  - B$(k+3)$
  - B$(k+4)$
  - Value
  - Decision
* - $k=5$
  - $13-17=-4$
  - $28+U_6=3$
  -
  -
  -
  -
  - $U_5=-4$
  - keep
* - $k=4$
  - $33-8=25$
  - $28+U_5=24$
  - $13+39+U_6=27$
  -
  -
  -
  - $U_4=24$
  - B4
* - $k=3$
  - $73-0=73$
  - $28+U_4=52$
  - $13+39+U_5=48$
  - $33+49+U_6=57$
  -
  -
  - $U_3=48$
  - B4
* - $k=2$
  - $143-0=143$
  - $28+U_3=76$
  - $13+39+U_4=76$
  - $33+49+U_5=78$
  - $73+55+U_6=103$
  -
  - $U_2=76$
  - B2, B3
* - start (age 2)
  - $330-0=330$
  - $39+U_2=115$
  - $20+49+U_3=117$
  - $60+55+U_4=139$
  - $130+60+U_5=186$
  - $230+60+U_6=265$
  - $V_1(2)=115$
  - B1
```

In the start row the columns B$k$, …, B$(k+4)$ are B1, …, B5. For example,
B3 keeps the original machine through years 1 and 2 (ages 2 and 3, cost
$20+40=60$), trades it at age 4 ($p-t(4)+c(0)=55$) and continues with
$U_4=24$.

The answer is again 115 with a purchase in year 1. The plan is recovered by
following the decisions: after B1 the firm is at stage 2, where B2 and B3
tie. B2 leads to stage 3, where B4 is optimal, then to stage 5, where the
machine is kept to the end: $\mathrm{BBKBK}$. B3 leads to stage 4, where B4 is
optimal, then to stage 5: $\mathrm{BKBBK}$. The values agree with the two-argument
solution, $U_k=V_k(1)$ for $k=2,\dots,6$ (compare {numref}`tab-05-02-two-arg`),
as they must.
:::

(sec-05-02-effort)=
## Computational effort of the two formulations

The effort of a DP is the number of states times the work per state.

- **Two arguments.** Stage $k$ has $k$ states by {eq}`eq-05-02-states`, so
  stages $1,\dots,N$ have $1+2+\dots+N=N(N+1)/2$ states. Each state needs 2
  additions and 1 comparison, for a total of $N(N+1)$ additions and
  $N(N+1)/2$ comparisons: 30 additions and 15 comparisons in the example.
- **One argument.** There are only $N$ values to compute ($U_N,\dots,U_2$
  and the start), but stage $k$ has $N-k+2$ options and the start has
  $N+1$. The comparisons number
  $\sum_{k=2}^{N}(N-k+1)+N=N(N+1)/2$, the same as before. If the running
  sums $\sum_{a=1}^{m}c(a)$ are computed once (in $N$ additions), every option
  costs about two additions, so the additions are again of order $N^2$ (35
  in the example).

Both formulations therefore need $O(N^2)$ operations, far fewer than the
$N2^N$ of enumeration. The one-argument formulation stores $N$ numbers
instead of $N(N+1)/2$, and its options are long-range moves rather than
one-year steps; the two-argument formulation is easier to extend, because
any new feature of the problem (a forced replacement age, year-dependent
data, an overhaul) only changes the arcs out of each state. When the age is
capped at a maximum $M$ (a machine of age $M$ must be replaced), the
two-argument model has at most $M$ states per stage and the one-argument
model at most $M$ options per stage, and both need $O(NM)$ operations.

(sec-05-02-horizon)=
## Boundary conditions at the end of the horizon

The boundary condition {eq}`eq-05-02-boundary` encodes what happens to the
machine when the planning horizon ends. Other end-of-horizon rules need only
a different boundary condition, not a different recurrence:

- **No salvage value** (the machine is scrapped): $V_{N+1}(i)=0$ for all $i$.
- **A required final state**, for example a machine no older than $m$ at the
  end: $V_{N+1}(i)=-s(i)$ for $i\le m$ and $V_{N+1}(i)=+\infty$ for $i>m$.
  An infinite boundary value excludes every plan that ends in that state.
- **A terminal valuation** other than resale, such as the value of the
  machine to the next planning period: $V_{N+1}(i)=-v(i)$ with the chosen
  value $v(i)$.

The choice matters because the last few decisions are driven by it, as the
year-5 decisions of {prf:ref}`ex-05-02-replacement` show. This **end-of-horizon
effect** is a reason to plan over a horizon longer than the decisions that
will actually be implemented, and to re-plan every year.

For a very long horizon with stationary data, the optimal policy settles
into a cycle: buy a new machine, keep it for $m$ years, trade it at age $m$,
repeat. One cycle costs $p-t(m)+\sum_{a=0}^{m-1}c(a)$ and lasts $m$ years,
so the best cycle length minimizes the average cost per year

$$
A(m)=\frac{1}{m}\Bigl[p-t(m)+\sum_{a=0}^{m-1}c(a)\Bigr].
$$ (eq-05-02-average-cost)

With the data of {numref}`tab-05-02-data`, $A(1)=28$, $A(2)=52/2=26$,
$A(3)=82/3\approx27.33$, $A(4)=128/4=32$, so trading every two years is the
best stationary policy. The five-year plans $\mathrm{BBKBK}$ and $\mathrm{BKBBK}$ are close to
it: after the first purchase they keep machines for one or two years.
When money is discounted by a factor $\alpha<1$ per year, every future value
in {eq}`eq-05-02-recurrence` is multiplied by $\alpha$:
$V_k(i)=\min\{p-t(i)+c(0)+\alpha V_{k+1}(1),\ c(i)+\alpha V_{k+1}(i+1)\}$,
and the cycle criterion becomes a ratio of discounted sums.

(sec-05-02-overhaul)=
## Overhaul as a third option

Suppose that, besides keeping and replacing, the machine can be
**overhauled** at the start of a year. An overhauled machine performs better
than one that was not overhauled, but not as well as a new one. Its
performance depends on its actual age and on the time since its last
overhaul, but not on earlier overhauls. The data may now also depend on the
year $k$:

- $c(k,i,j)$, the cost of operating during year $k$ a machine of age $i$
  last overhauled at age $j$ ($j=0$ means never overhauled);
- $e(k,i,j)$, the cost of exchanging such a machine for a new one at the
  start of year $k$ (price less trade-in);
- $o(k,i)$, the cost of overhauling a machine of age $i$ at the start of
  year $k$;
- $s(i,j)$, the salvage value at the end of year $N$ of a machine that has
  just turned age $i$ and was last overhauled at age $j$.

To decide well, one must know the year, the age of the machine and the age
at which it was last overhauled; because earlier overhauls do not matter,
nothing else is needed. The state is the pair $(i,j)$ with $0\le j<i$ (or
$(i,j)=(0,0)$ for a brand-new machine). The three decisions move the state
as follows:

- **buy** a new machine ($\mathrm{B}$): $(i,j)\to(1,0)$, cost $e(k,i,j)+c(k,0,0)$;
- **keep** ($\mathrm{K}$): $(i,j)\to(i+1,j)$, cost $c(k,i,j)$;
- **overhaul** ($\mathrm{O}$): $(i,j)\to(i+1,i)$, cost $o(k,i)+c(k,i,i)$, since the
  machine runs in year $k$ as one of age $i$ just overhauled at age $i$.

With $V_k(i,j)$ the minimum cost of years $k,\dots,N$ from state $(i,j)$,

$$
V_k(i,j)=\min\left\{
\begin{aligned}
&\mathrm{B}:\ e(k,i,j)+c(k,0,0)+V_{k+1}(1,0)\\
&\mathrm{K}:\ c(k,i,j)+V_{k+1}(i+1,j)\\
&\mathrm{O}:\ o(k,i)+c(k,i,i)+V_{k+1}(i+1,i)
\end{aligned}
\right\},
$$ (eq-05-02-overhaul)

$$
V_{N+1}(i,j)=-s(i,j),
$$ (eq-05-02-overhaul-boundary)

and the answer cell is $V_1(i_0,j_0)$ for a machine of age $i_0$ last overhauled at
age $j_0$.

**How many states?** With an initial state $(i_0,j_0)$, the age at the start
of year $k$ is $i_0+k-1$ if the machine was never replaced, and otherwise one
of $1,\dots,k-1$, exactly as in {eq}`eq-05-02-states`. For each age, the age at
the last overhaul ranges over $j_0$ (if no overhaul has happened since year
1, for the original machine) or 0 (for a replacement), and over every
earlier age at which an overhaul could have been made. A machine of age
$m\le k-1$ bought in year $k-m$ can have $j\in\{0,1,\dots,m-1\}$, which is $m$
values, and the original machine, now of age $i_0+k-1$, has up to $k$ values
($j_0$, or any of the ages $i_0,\dots,i_0+k-2$ it had at the start of years
$1,\dots,k-1$). The number of states in year $k$ is therefore up to
$k+(1+2+\dots+(k-1))=k(k+1)/2$ instead of $k$, and the total work grows like
$N^3$ instead of $N^2$: each extra state coordinate multiplies the work by the
number of values it can take.

:::{prf:example} Replacement with overhaul
:label: ex-05-02-overhaul

Take the data of {numref}`tab-05-02-data` over a horizon of $N=3$ years,
starting with a two-year-old machine never overhauled, $(i_0,j_0)=(2,0)$.
An overhaul costs $o=9$ at any age. After an overhaul at age $j$, the
machine runs like a never-overhauled machine of age $i-j+1$, so
$c(i,j)=c(i-j+1)$ for $j\ge1$ and $c(i,0)=c(i)$. Trade-in and salvage values
depend on the actual age only: $e(i,j)=p-t(i)$ and $s(i,j)=s(i)$. No datum
depends on the year. Find the optimal plan and compare it with the best
plan without overhauls.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The reachable states are $(2,0)$ in year 1; $(1,0)$, $(3,0)$, $(3,2)$ in year
2; $(1,0)$, $(2,0)$, $(2,1)$, $(4,0)$, $(4,2)$, $(4,3)$ in year 3; and ten
states at the end, whose values $-s(i)$ are $-25$ for age 1, $-17$ for age 2,
$-8$ for age 3 and 0 for age 5. The replacement option costs
$p-t(i)+c(0)+V_{k+1}(1,0)$ and the overhaul option
$9+c(1)+V_{k+1}(i+1,i)=22+V_{k+1}(i+1,i)$.

```{list-table} Stage table of the overhaul example (thousands of dollars)
:name: tab-05-02-overhaul
:header-rows: 1

* - Year
  - State $(i,j)$
  - $\mathrm{B}$
  - $\mathrm{K}$
  - $\mathrm{O}$
  - $V_k(i,j)$
  - Decision
* - 3
  - $(1,0)$
  - $28-25=3$
  - $13-17=-4$
  - $22-17=5$
  - $-4$
  - $\mathrm{K}$
* - 3
  - $(2,0)$
  - $39-25=14$
  - $20-8=12$
  - $22-8=14$
  - 12
  - $\mathrm{K}$
* - 3
  - $(2,1)$
  - $39-25=14$
  - $20-8=12$
  - $22-8=14$
  - 12
  - $\mathrm{K}$
* - 3
  - $(4,0)$
  - $55-25=30$
  - $70+0=70$
  - $22+0=22$
  - 22
  - $\mathrm{O}$
* - 3
  - $(4,2)$
  - $55-25=30$
  - $40+0=40$
  - $22+0=22$
  - 22
  - $\mathrm{O}$
* - 3
  - $(4,3)$
  - $55-25=30$
  - $20+0=20$
  - $22+0=22$
  - 20
  - $\mathrm{K}$
* - 2
  - $(1,0)$
  - $28-4=24$
  - $13+12=25$
  - $22+12=34$
  - 24
  - $\mathrm{B}$
* - 2
  - $(3,0)$
  - $49-4=45$
  - $40+22=62$
  - $22+20=42$
  - 42
  - $\mathrm{O}$
* - 2
  - $(3,2)$
  - $49-4=45$
  - $20+22=42$
  - $22+20=42$
  - 42
  - $\mathrm{K}$, $\mathrm{O}$
* - 1
  - $(2,0)$
  - $39+24=63$
  - $20+42=62$
  - $22+42=64$
  - 62
  - $\mathrm{K}$
```

The minimum cost is $V_1(2,0)=62$. The plan keeps the machine in year 1
(state $(3,0)$ in year 2), overhauls it at age 3 at the start of year 2
(state $(4,3)$ in year 3) and keeps it in year 3; it ends with a five-year-old
machine of no salvage value. Check: $20+(9+13)+20-0=62$. Without the
overhaul option the same DP (drop column $\mathrm{O}$) gives 63, attained by
replacing in years 1 and 2 and keeping in year 3
($39+28+13-17=63$), so the option saves \$1 000. Enumerating all $3^3=27$
plans confirms that $\mathrm{KOK}$ is the unique optimum. The number of states grows
from $1,2,3$ per year without overhauls to $1,3,6$ with them, in line with
the count above.
:::

(sec-05-02-variants)=
## Variants

The replacement model adapts to many practical rules by changing the state,
the decision set or the boundary conditions. The following variants are
formulated and solved in {doc}`08-07-examples-dp-applications`:

- trading for a used machine of any age, with a forced replacement at a
  maximum age ({prf:ref}`ex-08-07-t2-q1`), and year-dependent data such as
  inflation ({prf:ref}`ex-08-07-t2-q2`);
- a leasing option, which needs a second value function for the state
  "no machine owned" ({prf:ref}`ex-08-07-t2-q3`);
- a rule that each machine must be traded no older than its predecessor,
  which adds the previous trade-in age to the state ({prf:ref}`ex-08-07-t2-q4`);
- overhauls with revenues instead of costs, overhauls that take two years,
  and a nonnegative-capital constraint, where a forward recursion needs one
  state variable fewer than the backward one ({prf:ref}`ex-08-07-t2-q5`,
  {prf:ref}`ex-08-07-t2-q6`, {prf:ref}`ex-08-07-t2-q7`); and
- a small business with rising running costs and an inflating purchase price
  ({prf:ref}`ex-08-07-a8-q7`).

## Summary and Common Mistakes

### Summary

- The replacement problem has stages = years, state = age of the machine and
  decisions keep or buy; its recurrence is {eq}`eq-05-02-recurrence` with the
  salvage boundary condition {eq}`eq-05-02-boundary` and $k$ possible states in
  year $k$ ({eq}`eq-05-02-states`).
- The five-year example costs \$115 000 with the two optimal plans $\mathrm{BBKBK}$
  and $\mathrm{BKBBK}$ ({numref}`tab-05-02-two-arg`); it is a shortest path in an
  acyclic age–year network ({numref}`fig-05-02-replacement-network`).
- Because every purchase restarts the machine at age 1, a one-argument value
  function $U_k=V_k(1)$ suffices after the first purchase
  ({eq}`eq-05-02-one-arg` to {eq}`eq-05-02-one-arg-start`); it stores $N$ values
  instead of $N(N+1)/2$, and both formulations need $O(N^2)$ operations.
- End-of-horizon rules are boundary conditions; for long horizons the best
  stationary cycle minimizes the average cost {eq}`eq-05-02-average-cost`.
- An overhaul option requires the age at the last overhaul as a second state
  coordinate ({eq}`eq-05-02-overhaul`), which raises the number of states per
  year from $k$ to about $k^2/2$ ({prf:ref}`ex-05-02-overhaul`).

### Common mistakes

- **Mixing up the timing of costs.** In this model a machine bought at the
  start of a year is run in that same year at cost $c(0)$, and the next
  state is age 1, not age 0.
- **Forgetting the salvage value**, or giving it the wrong sign: at the end
  it is a revenue, so $V_{N+1}(i)=-s(i)$ in a cost-minimization model.
- **Listing impossible states or missing possible ones.** In year $k$ the
  original machine has age $i_0+k-1$, any replacement an age between 1 and
  $k-1$; nothing else can occur.
- **Reporting one optimal plan when the stage table has ties.** A tie in the
  policy column on the optimal path means alternative optimal plans; follow
  every branch.
- **Too few state variables.** In the overhaul model, the age alone does not
  determine the future costs; a recurrence in $V_k(i)$ only would not satisfy
  the principle of optimality.

## References

- Bellman, R. (1957). *Dynamic programming*. Princeton University Press.
  (Chapter I: multistage decision processes; Chapter III, §3: the principle
  of optimality.)
- Dreyfus, S. E., & Law, A. M. (1977). *The art and theory of dynamic
  programming*. Academic Press. (Chapter 2: equipment replacement, including
  the one-argument formulation and the overhaul model.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 11: dynamic programming,
  Sections 11.1–11.3.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 12: deterministic dynamic programming, Section 12.3:
  the equipment replacement model.)
