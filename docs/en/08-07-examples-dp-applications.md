# DP Applications

This chapter works through the examples and problems of
{doc}`05-02-dp-equipment-replacement` and of
{doc}`05-03-dp-resource-allocation`, two applications of dynamic programming
(DP). The notation is that of those chapters and of {doc}`05-01-dp-fundamentals`:
$k$ is the stage and $N$ the number of stages, $V_k$ is the optimal value
function and $P_k$ the optimal policy of stage $k$, and decisions are named
by upright letters ($\mathrm{B}$ buy, $\mathrm{K}$ keep, $\mathrm{O}$
overhaul). In replacement problems the stage is the year and the state is the
age $i$ of the machine; in allocation problems the stage is the activity, the
state is the amount $y$ of resource left and the decision is the amount $x$
allocated. Every solution states the stages, states, decisions,
recurrence, boundary condition and answer cell (the value that solves the
problem); numeric
problems add their stage tables and the optimal policy. All numbers were
recomputed, and every optimum was checked by enumerating all feasible
decision sequences.

(sec-08-07-replacement)=
## Replacement variants

:::{prf:example} Trading for a machine of any age
:label: ex-08-07-t2-q1

In the replacement model of {ref}`sec-05-02-two-argument`, suppose that a
machine reaching age $M$ must be replaced, and that at the start of a year the
current machine of age $i$ can be traded for a machine of any age
$j\in\{0,1,\dots,M-1\}$, new or used. The cost of the exchange is $u(i,j)$,
with $u(i,0)=p-t(i)$ (a new machine) and $u(i,i)=0$ (keeping the machine).
Give a backward DP and count its operations.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Stages.** Years $k=1,\dots,N$; the end of the horizon is stage $N+1$.

**States.** The age $i$ of the machine at the start of year $k$. Since any
age can be bought, every $i\in\{1,\dots,M\}$ is possible from year 2 on;
year 1 has the single state $i_0$.

**Decisions.** The age $j\in\{0,\dots,M-1\}$ of the machine operated in year
$k$. Choosing $j=i$ keeps the machine (allowed only for $i\le M-1$, which the
range of $j$ ensures); the machine is age $j+1$ a year later.

**Recurrence.** With $V_k(i)$ the minimum cost of years $k,\dots,N$ from age $i$,

$$
V_k(i)=\min_{j=0,\dots,M-1}\bigl\{u(i,j)+c(j)+V_{k+1}(j+1)\bigr\}.
$$

**Boundary condition.** $V_{N+1}(i)=-s(i)$ for $i=1,\dots,M$.

**Answer cell.** $V_1(i_0)$; the optimal ages are read forward from the minimizers.

**Operation count.** Years $2,\dots,N$ have $M$ states each and year 1 has
one. Each state compares $M$ options, each needing two additions. The total
is $2M\bigl(M(N-1)+1\bigr)$ additions and $(M-1)\bigl(M(N-1)+1\bigr)$
comparisons, about $2M^2N$ and $M^2N$: linear in the horizon, quadratic in the
maximum age.
:::

:::{prf:example} Year-dependent costs
:label: ex-08-07-t2-q2

Suppose that in {prf:ref}`ex-08-07-t2-q1` all data depend on the year, for
instance because of inflation: $u_k(i,j)$ is the cost of trading an age-$i$
machine for an age-$j$ machine at the start of year $k$, and $c_k(j)$ the
cost of operating an age-$j$ machine during year $k$. Give the backward DP.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The stages, states and decisions are those of {prf:ref}`ex-08-07-t2-q1`; the
year is already the stage, so it only has to be carried into the data:

$$
\begin{aligned}
V_k(i)&=\min_{j=0,\dots,M-1}\bigl\{u_k(i,j)+c_k(j)+V_{k+1}(j+1)\bigr\},\\
V_{N+1}(i)&=-s(i).
\end{aligned}
$$

The answer cell is $V_1(i_0)$. Many more data are needed (one table per year), but
the number of states, options and operations is the same as before. Data that
depend on the stage never enlarge the state: the stage index is always known.
:::

:::{prf:example} Leasing as an option
:label: ex-08-07-t2-q3

Add a third option to the basic replacement model of
{ref}`sec-05-02-two-argument`: sell the owned machine for its salvage value $s(i)$
and lease a new machine for the year at cost $l$ (excluding operating cost).
A firm that leased in the previous year owns no machine; it can lease again or
buy a new machine at price $p$. The firm starts year 1 with an owned, new
machine. Give a backward DP.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**States.** Two kinds of situation occur at the start of a year: an owned
machine of age $i$, or no owned machine (after a lease). Each gets its own
value function:

- $V_k(i)$: minimum cost of years $k,\dots,N$ starting year $k$ with an owned
  machine of age $i$;
- $W_k$: the same when year $k$ starts without an owned machine.

**Recurrences.** An owned machine can be traded for a new one ($\mathrm{B}$), kept
($\mathrm{K}$) or sold so that a machine is leased ($\mathrm{L}$); a leased machine is always
new, so it runs at cost $c(0)$:

$$
V_k(i)=\min\left\{
\begin{aligned}
&\mathrm{B}:\ p-t(i)+c(0)+V_{k+1}(1)\\
&\mathrm{K}:\ c(i)+V_{k+1}(i+1)\\
&\mathrm{L}:\ -s(i)+l+c(0)+W_{k+1}
\end{aligned}
\right\},
$$

$$
W_k=\min\left\{
\begin{aligned}
&\mathrm{B}:\ p+c(0)+V_{k+1}(1)\\
&\mathrm{L}:\ l+c(0)+W_{k+1}
\end{aligned}
\right\}.
$$

For the new machine of year 1 ($i=0$) only $\mathrm{K}$ and $\mathrm{L}$ make sense.

**Boundary conditions.** $V_{N+1}(i)=-s(i)$ and $W_{N+1}=0$ (a leased machine
is simply returned).

**Answer cell.** $V_1(0)$. The two value functions call each other, but always at
the next stage, so they are computed together, year by year backwards.
:::

:::{prf:example} Trade-in ages that never increase
:label: ex-08-07-t2-q4

Give a backward DP for the basic replacement model of
{ref}`sec-05-02-two-argument` with the extra rule that, after the original
machine has been traded, every machine must be traded at or before the age at which its predecessor was traded: if
$a_m$ is the age at trade-in of the $m$-th machine traded, then
$a_{m+1}\le a_m$. The original machine has age $i_0$ at the start of year 1.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**States.** To decide, one must know the year, the age $i$ of the current
machine, and the latest age $z$ at which it may be traded, which is the age at
which the previous machine was traded. Before the first trade there is no
limit, written $z=\infty$. The state is $(i,z)$ with $i\le z$.

**Decisions and recurrence.** Buying sets the new limit to the age just
traded; keeping is allowed only while $i<z$:

$$
V_k(i,z)=\begin{cases}
\min\bigl\{\,p-t(i)+c(0)+V_{k+1}(1,i),\ \ c(i)+V_{k+1}(i+1,z)\bigr\}, & i<z,\\[2pt]
p-t(i)+c(0)+V_{k+1}(1,i), & i=z.
\end{cases}
$$

**Boundary condition.** $V_{N+1}(i,z)=-s(i)$ (the final sale is not a trade).

**Answer cell.** $V_1(i_0,\infty)$.

**Range of the states.** In year $k$, $i\in\{1,\dots,k-1\}\cup\{i_0+k-1\}$ as
in {eq}`eq-05-02-states`. A finite $z$ is the age of a machine traded at the
start of an earlier year $k'\le k-1$, so $z\le i_0+k'-1\le i_0+k-2$; and
$z\ge i$. A safe range is therefore $z\in\{i,\dots,i_0+k-2\}\cup\{\infty\}$.
Computing $V$ for states that cannot occur does no harm.

*(A range $z=i,\dots,k-1$ that is sometimes suggested is too small: when the
original machine is traded at an age above $k-1$, later states have such a
$z$. On random instances the DP with that range misses reachable states,
while the range above always covers them.)*
:::

:::{prf:example} Replacing a lawn mower
:label: ex-08-07-a8-q7

A 13-year-old runs a lawn-mowing business with 10 customers, each mowed 3
times a season at \$50 per mowing. A new mower has just been bought for \$200.
Its operating and maintenance cost is \$120 in its first year of service and
grows by 20% per year after that. A one-year-old mower can be sold for \$150,
and the resale value falls by 10% for every further year of age. The price of
a new mower rises by 10% per year. The business will run for four seasons
(through the season in which the owner is 16). The owner believes that buying
a new mower every two years is the most economical policy. Formulate a DP and
check the claim.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The revenue, $10\cdot3\cdot\$50=\$1\,500$ per season, does not depend on the
decisions, and the \$200 already paid is sunk, so the DP minimizes the net
cost of the mowers. As in {ref}`sec-05-02-two-argument`, $c(i)$ is the cost of a
season that starts with a mower of age $i$; with compound growth
$c(i)=120\cdot1.2^{\,i}$ (so $c(0)=120$ for a new mower), the resale value at
age $i\ge1$ is $s(i)=150\cdot0.9^{\,i-1}$, and a new mower costs
$p_k=200\cdot1.1^{\,k-1}$ in season $k$.

**Stages.** Seasons $k=1,2,3,4$; stage 5 is the end.

**States.** The age $i$ of the mower at the start of season $k$; $i=0$ in
season 1.

**Decisions.** Keep ($\mathrm{K}$): pay $c(i)$, next age $i+1$. Buy
($\mathrm{B}$, for $i\ge1$): sell the mower for $s(i)$, buy a new one at
$p_k$ and run it at $c(0)$; next age 1.

**Recurrence.**
$V_k(i)=\min\{\mathrm{K}:\ c(i)+V_{k+1}(i+1),\ \mathrm{B}:\ p_k-s(i)+c(0)+V_{k+1}(1)\}$.

**Boundary condition.** $V_5(i)=-s(i)$: the mower is sold after the last
season.

**Answer cell.** $V_1(0)$.

| Season $k$ | Age $i$ | $\mathrm{K}$ | $\mathrm{B}$ | $V_k(i)$ | Decision |
|---|---|---|---|---|---|
| 4 | 1 | $144-135=9.00$ | $266.20-150+120-150=86.20$ | 9.00 | $\mathrm{K}$ |
| 4 | 2 | $172.80-121.50=51.30$ | $266.20-135+120-150=101.20$ | 51.30 | $\mathrm{K}$ |
| 4 | 3 | $207.36-109.35=98.01$ | $266.20-121.50+120-150=114.70$ | 98.01 | $\mathrm{K}$ |
| 3 | 1 | $144+51.30=195.30$ | $242-150+120+9=221.00$ | 195.30 | $\mathrm{K}$ |
| 3 | 2 | $172.80+98.01=270.81$ | $242-135+120+9=236.00$ | 236.00 | $\mathrm{B}$ |
| 2 | 1 | $144+236=380.00$ | $220-150+120+195.30=385.30$ | 380.00 | $\mathrm{K}$ |
| 1 | 0 | $120+380=500.00$ | — | 500.00 | $\mathrm{K}$ |

The minimum net cost is \$500.00: keep the mower in seasons 1 and 2, replace
it at the start of season 3 (age 2, price \$242), keep the new one in season 4
and sell it at age 2. Check: $120+144+(242-135+120)+144-135=500$. This is
exactly "a new mower every two years", so the claim is right for a four-season
horizon. The next best plans cost \$505.30 (replace after one year) and
\$534.81 (never replace). The profit over the four seasons is
$4\cdot1\,500-500=\$5\,500$.

*(An earlier version of this solution reports \$458 for the same plan; its
table mixes two readings of the percentage changes and adds $144+239$ as 338
instead of 383. With compound rates the minimum is \$500.00; with linear
rates, +\$24 cost and −\$20 resale value per year and +\$20 price per year,
the same plan is optimal at \$508.)* If the business ran for only three
seasons, never replacing would be best (\$315.30), so the conclusion depends
on the horizon.
:::

(sec-08-07-overhaul)=
## Overhaul and capital constraints

The problems of this section use the overhaul model of
{ref}`sec-05-02-overhaul`: state $(i,j)$ = age and age at the last overhaul,
exchange cost $e(k,i,j)$, overhaul cost $o(k,i)$ and salvage value $s(i,j)$.

:::{prf:example} Overhaul with revenues
:label: ex-08-07-t2-q5

In the overhaul model, suppose that a machine produces a positive net revenue
$n(k,i,j)$ during year $k$ (when of age $i$ and last overhauled at age $j$)
instead of a cost. Net revenue plus salvage value, minus exchange and
overhaul costs, is to be maximized. Give a DP formulation.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The stages, states and transitions are those of {eq}`eq-05-02-overhaul`; costs
become negative terms and the minimum a maximum. With $V_k(i,j)$ the maximum
net return of years $k,\dots,N$ from state $(i,j)$:

$$
V_k(i,j)=\max\left\{
\begin{aligned}
&\mathrm{B}:\ -e(k,i,j)+n(k,0,0)+V_{k+1}(1,0)\\
&\mathrm{K}:\ n(k,i,j)+V_{k+1}(i+1,j)\\
&\mathrm{O}:\ -o(k,i)+n(k,i,i)+V_{k+1}(i+1,i)
\end{aligned}
\right\},
$$

with boundary condition $V_{N+1}(i,j)=s(i,j)$ (salvage is now a revenue).
The answer cell is $V_1(i_0,j_0)$ for the initial machine.
:::

:::{prf:example} An overhaul that takes two years
:label: ex-08-07-t2-q6

Modify the formulation of {prf:ref}`ex-08-07-t2-q5` for an overhaul that
takes two years, during which the machine ages but earns no revenue. An overhaul cannot be started at the
start of year $N$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Only the overhaul transition changes: an overhaul started at the start of
year $k$ on a machine of age $i$ costs $o(k,i)$, earns nothing in years $k$
and $k+1$, and leads to state $(i+2,i)$ at the start of year $k+2$:

$$
V_k(i,j)=\max\left\{
\begin{aligned}
&\mathrm{B}:\ -e(k,i,j)+n(k,0,0)+V_{k+1}(1,0)\\
&\mathrm{K}:\ n(k,i,j)+V_{k+1}(i+1,j)\\
&\mathrm{O}:\ -o(k,i)+V_{k+2}(i+2,i)
\end{aligned}
\right\}.
$$

The boundary conditions are $V_{N+1}(i,j)=s(i,j)$ and
$V_{N+2}(i,j)=-\infty$. The second one makes the overhaul option worthless in
year $N$, which enforces the rule without a special case; an overhaul started
in year $N-1$ ends exactly at the horizon. The recurrence now reaches two
stages ahead, so the backward pass must have $V_{k+1}$ and $V_{k+2}$ at hand.
The answer cell is $V_1(i_0,j_0)$.
:::

:::{prf:example} Overhaul with a nonnegative capital position
:label: ex-08-07-t2-q7

Consider the revenue model of {prf:ref}`ex-08-07-t2-q5`, starting year 1 with a
new machine and zero capital. The capital position $C_k$ at the start of year
$k$ evolves as $C_1=0$ and $C_{k+1}=C_k+{}$(net revenue of year $k$)${}-{}$(exchange
or overhaul cost paid at the start of year $k$). Exchange and overhaul costs
are paid before any revenue of the year comes in, and the capital position
must never be negative. Maximize the final capital $C_{N+1}$ plus the salvage
value. Give a backward and a forward DP.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Backward procedure.** Whether an exchange or overhaul is affordable depends
on the capital in hand, so the capital becomes a state variable. Let
$V_k(i,j,C)$ be the maximum final capital plus salvage value when year $k$
starts with a machine in state $(i,j)$ and capital $C$:

$$
V_k(i,j,C)=\max\left\{
\begin{aligned}
&\mathrm{B}\ (\text{if }C\ge e(k,i,j)):\ V_{k+1}\bigl(1,0,\ C-e(k,i,j)+n(k,0,0)\bigr)\\
&\mathrm{K}:\ V_{k+1}\bigl(i+1,j,\ C+n(k,i,j)\bigr)\\
&\mathrm{O}\ (\text{if }C\ge o(k,i)):\ V_{k+1}\bigl(i+1,i,\ C-o(k,i)+n(k,i,i)\bigr)
\end{aligned}
\right\},
$$

with $V_{N+1}(i,j,C)=C+s(i,j)$. The answer cell is $V_1(0,0,0)$. The capital can
take many values, so this DP has a large state space.

*(Equivalently, the year's net revenue can be written outside $V$ if $V$
counts only what is earned from year $k$ on, with boundary $s(i,j)$; the two
conventions must not be combined.)*

**Forward procedure.** As in the forward recursion of {ref}`sec-05-01-forward`,
write $F_k$ for a value accumulated from the start. Let $F_k(i,j)$ be the
**largest** capital that can be
held at the start of year $k$ in machine state $(i,j)$, over all feasible
histories ($-\infty$ if the state cannot be reached). The key observation is
that more capital is never worse: every constraint has the form $C\ge$ cost, so
any decision available with less capital is also available with more, and the
final objective increases with $C$. Among all histories reaching $(i,j)$ in
year $k$, only the richest one matters, so the capital need not be a state
variable. Starting from $F_1(0,0)=0$, each state passes its capital forward
along every affordable decision:

- keep: $(i,j)\to(i+1,j)$, candidate $F_k(i,j)+n(k,i,j)$;
- overhaul ($i\ge1$, if $F_k(i,j)\ge o(k,i)$): $(i,j)\to(i+1,i)$, candidate
  $F_k(i,j)-o(k,i)+n(k,i,i)$;
- buy (if $F_k(i,j)\ge e(k,i,j)$): $(i,j)\to(1,0)$, candidate
  $F_k(i,j)-e(k,i,j)+n(k,0,0)$;

and $F_{k+1}$ of each state is the largest candidate arriving there. In
particular

$$
\begin{aligned}
F_{k+1}(1,0)=\max\Bigl\{&F_k(0,0)+n(k,0,0),\\
&\max_{(i,j):\,F_k(i,j)\ge e(k,i,j)}\bigl[F_k(i,j)-e(k,i,j)+n(k,0,0)\bigr]\Bigr\},
\end{aligned}
$$

where the first term (keeping the new machine) occurs only in year 1. The
answer is $\max_{(i,j)}\bigl[F_{N+1}(i,j)+s(i,j)\bigr]$.

*(A version of this recursion that obtains $F_{k+1}(1,0)$ only from
replacements misses the first year, in which the new machine is simply kept;
with the first term included, the forward and backward procedures agree on
every test instance.)*

**Why it matters.** The forward procedure has one state variable fewer than
the backward one, because the value of the optimal value function itself
plays the role of the capital. The reduction rests on the monotonicity above.
If there were also an upper limit on capital (say, purchases forbidden when
capital exceeds some amount), a history with less than the maximum capital
might be the only one that allows a later purchase, and the reduction would
no longer be valid.
:::

(sec-08-07-one-resource)=
## Allocation: one resource and knapsack

:::{prf:example} Allocation as a longest path
:label: ex-08-07-l9-longest-path

Take the four return functions of {numref}`tab-05-03-returns` and allocate
$X=5$ units, every unit being used. Model the problem as a layered network:
node $(k,y)$ means that activity $k$ is next and $y$ units are left, and the
arc to $(k+1,y')$ allocates $y-y'$ units to activity $k$. Give the arc
lengths, solve the longest-path problem by DP and state all optimal
allocations.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

```{figure} figures/08-07-allocation-network.svg
:name: fig-08-07-allocation-network
:width: 560px
:class: dark-light
:alt: Layered network of the allocation of five units among four activities with the two longest paths highlighted

Longest-path network of the allocation of 5 units. Arcs out of stage 1 carry
$r_1(x)$; the numbers right of the last column are the terminal lengths
$r_4(y)$. The two longest paths, of length 14, are highlighted with their
arc lengths.
```

**Arc lengths.** In the network of {numref}`fig-08-07-allocation-network`,
arc $(k,y)\to(k+1,y')$ has length $r_k(y-y')$ for $y'\le y$; the last activity takes the rest, so node $(4,y)$ has terminal
length $r_4(y)$. There are $6+21+21=48$ arcs between the columns.

**Recurrence.** $V_k(y)$ is the length of a longest path from node $(k,y)$
to the end: $V_k(y)=\max_{x\le y}\{r_k(x)+V_{k+1}(y-x)\}$, with $V_4=r_4$;
the answer cell is $V_1(5)$. This is {eq}`eq-05-03-recurrence` with $X=5$.

| $y$ | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| $V_4(y)=r_4(y)$ | 0 | 1 | 3 | 6 | 9 | 12 |
| $V_3(y)$ | 0 | 2 | 4 | 6 | 9 | 12 |
| $P_3(y)$ | 0 | 1 | 2 | 0, 3 | 0 | 0 |
| $V_2(y)$ | 0 | 2 | 4 | 6 | 9 | 13 |
| $P_2(y)$ | 0 | 0 | 0 | 0 | 0 | 5 |

At stage 1,
$V_1(5)=\max\{0+13,\ 3+9,\ 7+6,\ 10+4,\ 12+2,\ 13+0\}=14$, attained by
$x_1=3$ and by $x_1=4$.

- $x_1=3$: $y=2$ at stage 2, $P_2(2)=0$; $y=2$ at stage 3, $P_3(2)=2$; the
  allocation is $(3,0,2,0)$ with return $10+4=14$.
- $x_1=4$: $y=1$, $P_2(1)=0$, $P_3(1)=1$; the allocation is $(4,0,1,0)$ with
  return $12+2=14$.

The longest paths have length 14 and there are exactly two of them
(confirmed by enumeration). Note that with eight units the best allocation
gives nothing to activity 3 and six units to activity 2
({prf:ref}`ex-05-03-one-resource`): the optimal allocation of a smaller amount is
not a part of the optimal allocation of a larger one when returns are not
concave.
:::

:::{prf:example} Choosing elective courses
:label: ex-08-07-a8-q3

A student must take ten elective courses from four departments, at least one
from each department. The student rates the knowledge gained from $n$ courses
of each department on a 100-point scale:

| Department | $n=1$ | 2 | 3 | 4 | 5 | 6 | $\ge7$ |
|---|---|---|---|---|---|---|---|
| I | 25 | 50 | 60 | 80 | 100 | 100 | 100 |
| II | 20 | 70 | 90 | 100 | 100 | 100 | 100 |
| III | 40 | 60 | 80 | 100 | 100 | 100 | 100 |
| IV | 10 | 20 | 30 | 40 | 50 | 60 | 70 |

Formulate a DP that maximizes total knowledge and solve it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Let $b_k(n)$ be the knowledge from $n$ courses of department $k$ (I to IV
$=1$ to 4).

**Stages.** Departments $k=1,2,3,4$.

**States.** The number $y$ of courses still to be chosen when department $k$
is considered. Departments $k,\dots,4$ need at least one course each and
departments $1,\dots,k-1$ took at least one each, so
$5-k\le y\le 11-k$; at stage 1, $y=10$.

**Decisions.** $x\in\{1,\dots,y-(4-k)\}$ courses from department $k$, leaving
at least one for each later department.

**Recurrence and boundary.**
$V_k(y)=\max_{1\le x\le y-(4-k)}\{b_k(x)+V_{k+1}(y-x)\}$, $V_4(y)=b_4(y)$.

**Answer cell.** $V_1(10)$.

| $y$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| $V_4(y)$ | 10 | 20 | 30 | 40 | 50 | 60 | 70 | | | |
| $V_3(y)$ ($x$) | | 50 (1) | 70 (2) | 90 (3) | 110 (4) | 120 (4) | 130 (4) | 140 (4) | | |
| $V_2(y)$ ($x$) | | | 70 (1) | 120 (2) | 140 (2, 3) | 160 (2, 3) | 180 (2, 3) | 200 (3) | 210 (3, 4) | |
| $V_1(y)$ ($x$) | | | | | | | | | | 250 (2) |

At stage 1 the options $x_1=1,\dots,7$ give $235, 250, 240, 240, 240, 220,
170$, so $V_1(10)=250$ with $x_1=2$. Forward: $y=8$ at stage 2, where
$x_2=3$ ($90+110=200$); $y=5$ at stage 3, where $x_3=4$ ($100+10=110$); and
$x_4=1$. The student should take **2, 3, 4 and 1** courses from departments
I to IV, for $50+90+100+10=250$ points; this plan is the unique optimum.
:::

:::{prf:example} Selecting housing applications
:label: ex-08-07-a8-q4

A volunteer organization that builds homes for families in need has five
applications. It can count on at most 23 volunteers over the next six months.
Each application has a need score (higher means more need) and requires a
number of volunteers:

| Application | House size (ft²) | Score | Volunteers |
|---|---|---|---|
| 1 | 1 200 | 78 | 7 |
| 2 | 1 000 | 64 | 4 |
| 3 | 1 100 | 68 | 6 |
| 4 | 1 000 | 62 | 5 |
| 5 | 1 200 | 85 | 8 |

Which applications should be approved to maximize the total score? Formulate
and solve a DP.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

This is a 0–1 knapsack problem ({eq}`eq-05-03-knapsack`) with worth $h_k$ =
score and size $a_k$ = volunteers; the house sizes play no role.

**Stages.** Applications $k=1,\dots,5$.

**States.** The number $y$ of volunteers still available, $0\le y\le23$.

**Decisions.** Reject ($x=0$) or approve ($x=1$, only if $a_k\le y$).

**Recurrence and boundary.**
$V_k(y)=\max\{V_{k+1}(y),\ h_k+V_{k+1}(y-a_k)\}$, the second term only if
$a_k\le y$; $V_6(y)=0$.

**Answer cell.** $V_1(23)$.

Only the states that can be reached from $y=23$ need to be listed. Each cell
gives reject / approve / $V_k(y)$:

| Stage | Reachable states $y$: reject / approve / $V_k(y)$ |
|---|---|
| 5 | 23, 19, 18, 17, 16, 14, 13, 12, 11, 10, 8: 0 / 85 / 85; 7, 6, 5, 1: 0 / – / 0 |
| 4 | 23, 19, 17, 16, 13: 85 / 147 / 147; 12, 10: 85 / 62 / 85; 6: 0 / 62 / 62 |
| 3 | 23: 147 / 215 / 215; 19: 147 / 215 / 215; 16: 147 / 153 / 153; 12: 85 / 130 / 130 |
| 2 | 23: 215 / 279 / 279; 16: 153 / 194 / 194 |
| 1 | 23: 279 / 272 / 279 |

$V_1(23)=279$ with application 1 rejected. Forward: approve application 2
($y=23\to19$), approve 3 ($19\to13$), approve 4 ($13\to8$), approve 5
($8\to0$). The organization should **approve applications 2, 3, 4 and 5**,
using all 23 volunteers, for a total score of $64+68+62+85=279$. Approving
application 1 instead costs 7 points ($V$ = 272).
:::

:::{prf:example} A continuous min–max allocation
:label: ex-08-07-a8-q5

Formulate a DP for

$$
\begin{aligned}
\min\ & \max\{\varphi_1(x_1),\dots,\varphi_N(x_N)\}\\
\st\ & x_1+\dots+x_N=X,\\
& x_k\ge0,\quad k=1,\dots,N,
\end{aligned}
$$

where the $x_k$ are real numbers. Solve it for $N=3$, $X=10$,
$\varphi_1(x)=x+5$, $\varphi_2(x)=5x+3$ and $\varphi_3(x)=x-2$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Stages.** Variables $k=1,\dots,N$.

**States.** The amount $y\ge0$ still to be distributed among
$x_k,\dots,x_N$.

**Decisions.** $x_k\in[0,y]$.

**Recurrence.** Let $V_k(y)$ be the smallest possible value of
$\max\{\varphi_k(x_k),\dots,\varphi_N(x_N)\}$ with $x_k+\dots+x_N=y$. The
maximum splits like a sum does, so the principle of optimality holds:

$$
\begin{aligned}
V_k(y)&=\min_{0\le x\le y}\max\bigl\{\varphi_k(x),\ V_{k+1}(y-x)\bigr\},\\
V_N(y)&=\varphi_N(y).
\end{aligned}
$$

**Answer cell.** $V_1(X)$.

**Stage 3.** $V_3(y)=y-2$ with $x_3=y$.

**Stage 2.** $V_2(y)=\min_{0\le x\le y}\max\{5x+3,\ y-x-2\}$. The first term
increases and the second decreases in $x$, so the minimum is where they meet,
$5x+3=y-x-2$, that is, $x=(y-5)/6$, if that point is feasible ($y\ge5$);
otherwise at $x=0$:

$$
V_2(y)=\begin{cases}3, & 0\le y\le5\ \ (x_2=0),\\[2pt] \dfrac{5y-7}{6}, & y>5\ \ \bigl(x_2=\tfrac{y-5}{6}\bigr).\end{cases}
$$

*(An earlier version writes the second branch in terms of the decision,
$(5x_2-7)/6$; it is a function of the state $y$.)*

**Stage 1.** $V_1(y)=\min_{0\le x\le y}\max\{x+5,\ V_2(y-x)\}$. At $x=0$ the
value is $\max\{5,V_2(y)\}=5$ as long as $V_2(y)\le5$, that is, $y\le37/5$.
Otherwise the increasing $x+5$ meets the decreasing $(5(y-x)-7)/6$ at
$x=(5y-37)/11$:

$$
V_1(y)=\begin{cases}5, & y\le\tfrac{37}{5},\\[2pt] \dfrac{5y+18}{11}, & y>\tfrac{37}{5}\ \ \bigl(x_1=\tfrac{5y-37}{11}\bigr).\end{cases}
$$

**Solution.** $V_1(10)=68/11\approx6.182$ with $x_1=13/11$. Then
$y=97/11$ at stage 2 and $x_2=(97/11-5)/6=7/11$; finally
$x_3=90/11$. Check: $\varphi_1=13/11+5$, $\varphi_2=35/11+3$ and
$\varphi_3=90/11-2$ all equal $68/11$, and $13/11+7/11+90/11=10$. Equalizing
the terms is typical of min–max allocations; the same optimum follows from
the linear program (LP) $\min\{u:\ u\ge\varphi_k(x_k),\ \sum_k x_k=10,\ \vx\ge\vzero\}$.
:::

(sec-08-07-several-resources)=
## Several resources and side conditions

:::{prf:example} Two resources
:label: ex-08-07-t2-q8

$X_1$ units of resource 1 and $X_2$ units of resource 2 are to be allocated
among $N$ activities; allocating $x_1$ and $x_2$ units to activity $k$
returns $r_k(x_1,x_2)$. Write the DP formulation.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Stages are the activities; the state is the pair $(y_1,y_2)$ of units left;
the decision is the pair $(x_1,x_2)$ with $x_1\le y_1$, $x_2\le y_2$. The
recurrence and boundary condition are {eq}`eq-05-03-two-resources` and
{eq}`eq-05-03-two-boundary`,

$$
\begin{aligned}
V_k(y_1,y_2)&=\max_{\substack{x_1=0,\dots,y_1\\ x_2=0,\dots,y_2}}
\bigl\{r_k(x_1,x_2)+V_{k+1}(y_1-x_1,y_2-x_2)\bigr\},\\
V_N(y_1,y_2)&=r_N(y_1,y_2),
\end{aligned}
$$

and the answer cell is $V_1(X_1,X_2)$. A numeric instance with
$X_1=X_2=3$ is solved in {prf:ref}`ex-05-03-two-resources`.
:::

:::{prf:example} Three resources and the dimension of the problem
:label: ex-08-07-t2-q10

Suppose there are three resources, with initial amounts $X_1$, $X_2$, $X_3$
and returns $r_k(x_1,x_2,x_3)$, $k=1,\dots,N$. Give the DP formulation. What
is the dimension of the problem, and how much work does it need?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The state is the triple $(y_1,y_2,y_3)$ of units left:

$$
V_k(y_1,y_2,y_3)=\max_{\substack{x_l=0,\dots,y_l\\ l=1,2,3}}
\bigl\{r_k(x_1,x_2,x_3)+V_{k+1}(y_1-x_1,\,y_2-x_2,\,y_3-x_3)\bigr\},
$$

with $V_N=r_N$ and answer $V_1(X_1,X_2,X_3)$.

**Dimension.** The problem is **three-dimensional**: the value function has
three state variables besides the stage index. Each stage has
$(X_1+1)(X_2+1)(X_3+1)$ states, and, by the count of {eq}`eq-05-03-two-count`
extended to three factors, it needs

$$
\prod_{l=1}^{3}\frac{(X_l+1)(X_l+2)}{2}
$$

evaluations of $r_k+V_{k+1}$. With ten units of each resource that is 1 331
states and $66^3=287\,496$ evaluations per stage, against 121 states and
4 356 evaluations with two resources and 66 with one: each resource adds a
factor of about $X^2/2$ to the work. Storage grows the same way, since all of
$V_{k+1}$ must be kept.
:::

:::{prf:example} A minimum total for some activities
:label: ex-08-07-t2-q11

In the one-resource allocation problem {eq}`eq-05-03-model`, require in addition
that a specified group of $M$ of the $N$ activities ($M<N$) receives at least
$Y$ units in total ($Y<X$). Give a DP.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Number the activities so that the specified ones are $1,\dots,M$. The units
given to them are then exactly the units used in stages $1,\dots,M$, and the
condition says that **at most $X-Y$ units** may be left after stage $M$.

- For $k=M+1,\dots,N$: the usual recurrence {eq}`eq-05-03-recurrence` with
  boundary $V_N(y)=r_N(y)$, computed only for $y=0,1,\dots,X-Y$.
- For $k=M$: if $y\le X-Y$ the usual recurrence applies; if $y>X-Y$ the
  allocation must bring the state down to $X-Y$ or less:

  $$
  V_M(y)=\max_{x=y-(X-Y),\dots,y}\bigl\{r_M(x)+V_{M+1}(y-x)\bigr\}.
  $$

- For $k<M$: the usual recurrence for $y=0,\dots,X$.

The answer cell is $V_1(X)$. The side condition costs nothing extra: it only
restricts the states after stage $M$ and the decisions at stage $M$.
:::

:::{prf:example} At most $M$ activities used
:label: ex-08-07-t2-q12

In the allocation problem {eq}`eq-05-03-model`, at most $M$ of the $N$
allocations ($M<N$) may be nonzero. Give a DP.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**States.** Besides the units left $y$, one must know how many activities
have already been used. Let the state at stage $k$ be $(y,l)$, where $l\le M$
is the number of activities among $1,\dots,k-1$ with a nonzero allocation.

**Recurrence.**

$$
\begin{aligned}
V_k(y,l)=\max\Bigl\{&r_k(0)+V_{k+1}(y,l),\\
&\max_{x=1,\dots,y}\bigl[r_k(x)+V_{k+1}(y-x,\,l+1)\bigr]\Bigr\},
\end{aligned}
$$

where $V_{k+1}(\cdot,l+1)=-\infty$ when $l+1>M$ (so with $l=M$ only a zero
allocation is possible).

**Boundary condition** (the last activity takes the rest):

$$
V_N(y,l)=\begin{cases}r_N(0), & y=0,\ l\le M,\\ r_N(y), & y>0,\ l<M,\\ -\infty, & y>0,\ l=M.\end{cases}
$$

**Answer cell.** $V_1(X,0)$. The extra state variable multiplies the number of
states by at most $M+1$.
:::

:::{prf:example} A charge on the number of activities used
:label: ex-08-07-t2-q13

In the allocation problem {eq}`eq-05-03-model`, a charge $h(l)$ is deducted from
the total return when exactly $l$ of the $N$ activities receive a nonzero
allocation ($h(0)=0$). Give an efficient DP.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

As in {prf:ref}`ex-08-07-t2-q12`, the state is $(y,l)$ with $l$ the number of
activities used so far; now $l$ is unrestricted ($l\le k-1$). The charge
depends only on the final count, so it is applied once, in the boundary
condition:

$$
\begin{aligned}
V_k(y,l)=\max\Bigl\{&r_k(0)+V_{k+1}(y,l),\\
&\max_{x=1,\dots,y}\bigl[r_k(x)+V_{k+1}(y-x,\,l+1)\bigr]\Bigr\},
\end{aligned}
$$

$$
V_N(y,l)=\begin{cases}r_N(y)-h(l+1), & y>0,\\ r_N(0)-h(l), & y=0.\end{cases}
$$

The answer cell is $V_1(X,0)$. Charging $h$ at the end keeps the recurrence as
simple as in {prf:ref}`ex-08-07-t2-q12`; equivalently one could deduct
$h(l+1)-h(l)$ whenever an activity is used. Stage $k$ has only $k$ values of
$l$, so the work is about $N/2$ times that of the basic problem.
:::

:::{prf:example} Profit with a goodwill target
:label: ex-08-07-t2-q14

A company has $N$ divisions. Allocating $x$ units of a resource to division
$k$ produces a dollar return $r_k(x)$ and a goodwill score $\gamma_k(x)$ (an
integer, possibly negative: a very profitable division may harm the company's
image). $X$ units are available and need not all be used. Maximize the dollar
return subject to a total goodwill of at least $\Gamma$. Give a DP
formulation, including the range of the goodwill state.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**States.** The units left $y$ and the goodwill $G$ accumulated by divisions
$1,\dots,k-1$.

**Recurrence.** With $V_k(y,G)$ the maximum dollar return from divisions
$k,\dots,N$,

$$
V_k(y,G)=\max_{x=0,\dots,y}\bigl\{r_k(x)+V_{k+1}\bigl(y-x,\ G+\gamma_k(x)\bigr)\bigr\}.
$$

**Boundary condition.** $V_{N+1}(y,G)=0$ if $G\ge\Gamma$ and $-\infty$ if
$G<\Gamma$: an ending with too little goodwill is infeasible. Unused units are
allowed, so no condition on $y$ is needed.

**Answer cell.** $V_1(X,0)$.

**Range of $G$.** Let $\bar\gamma_l=\max_x\gamma_l(x)$ and
$\underline\gamma_l=\min_x\gamma_l(x)$. At stage $k$, $G$ lies between
$\sum_{l=1}^{k-1}\underline\gamma_l$ and $\sum_{l=1}^{k-1}\bar\gamma_l$, a
conservative range. *(The sums run over the divisions $l=1,\dots,k-1$; a
version that sums the term of division $k$ itself is a slip.)* If all
goodwill scores are nonnegative, every $G\ge\Gamma$ can be replaced by
$\Gamma$, which caps the range; with negative scores this is not allowed,
because goodwill above $\Gamma$ may be needed to absorb later losses.
:::

(sec-08-07-capital-budgeting)=
## Capital budgeting

:::{prf:example} Capital budgeting with a dearer second resource
:label: ex-08-07-l9-budget-dearer

Solve the capital-budgeting problem of {prf:ref}`ex-05-03-capital-budgeting`
(budget \$6, four activities, the returns of {numref}`tab-05-03-budget-returns`)
when resource 2 costs \$2 per unit instead of \$1. Resource 1 still costs
\$3.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The model is {eq}`eq-05-03-budget-model` with $c_1=3$, $c_2=2$, $B=6$. Only the
pairs with $3x_1+2x_2\le6$ are affordable: $(0,0)$, $(0,1)$, $(0,2)$, $(0,3)$,
$(1,0)$, $(1,1)$, $(2,0)$.

**Formulation.** State $b$ = dollars left; collapse each stage with
$\hat r_k(\beta)=\max\{r_k(x_1,x_2):\ 3x_1+2x_2\le\beta\}$ ({eq}`eq-05-03-hat-r`)
and allocate money with {eq}`eq-05-03-budget-2`, $V_5\equiv0$. The answer is
$V_1(6)$.

| $\beta$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| $\hat r_1$ | 0 | 0 | 5 | 5 | 9 | 9 | 11 |
| $\hat r_2$ | 0 | 0 | 3 | 3 | 5 | 5 | 7 |
| $\hat r_3$ | 0 | 0 | 1 | 4 | 4 | 6 | 11 |
| $\hat r_4$ | 0 | 0 | 1 | 7 | 7 | 10 | 13 |

Now an odd dollar cannot be spent on resource 2 alone, and the best use of
\$6 in activities 3 and 4 is two units of resource 1 ($r_3(2,0)=11$,
$r_4(2,0)=13$).

| $b$ | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| $V_4(b)=\hat r_4(b)$ | 0 | 0 | 1 | 7 | 7 | 10 | 13 |
| $V_3(b)$ | 0 | 0 | 1 | 7 | 7 | 10 | 13 |
| $V_2(b)$ | 0 | 0 | 3 | 7 | 7 | 10 | 13 |
| $V_1(b)$ | | | | | | | 13 |

At $b=6$ the best split at every stage gives \$0 to activities 1, 2 and 3:
$V_1(6)=\max\{0+13,\ 0+10,\ 5+7,\ 5+7,\ 9+3,\ 9+0,\ 11+0\}=13$ for
$\beta_1=0,\dots,6$, attained only by $\beta_1=0$, and likewise at stages 2
and 3. The unique optimal plan is to
give the whole budget to **activity 4, which buys two units of resource 1**,
for a return of $r_4(2,0)=13$ (enumeration confirms it). Raising the price of
resource 2 from \$1 to \$2 lowers the best return from 19 to 13.
:::

(sec-08-07-sequential)=
## Other sequential decisions: rental and livestock

:::{prf:example} Renting cars for tours
:label: ex-08-07-a8-q6

A tour operator needs 7, 4, 7 and 8 rental cars in the next four weeks. A
local dealer rents cars at \$220 per car per week, plus a flat fee of \$500
for every rental transaction (each time additional cars are picked up).
Cars that are no longer needed can be returned free of charge, or kept for a
later week at the weekly rent. Formulate a DP that minimizes the rental cost
and solve it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Let $d=(7,4,7,8)$ be the weekly needs and $d_{\max}=8$. Holding more than 8
cars never helps.

**Stages.** Weeks $k=1,2,3,4$.

**States.** The number $y$ of cars held in the previous week ($y=0$ before
week 1, and $d_{k-1}\le y\le8$ afterwards).

**Decisions.** The number $x\in\{d_k,\dots,8\}$ of cars held in week $k$. If
$x>y$, new cars are picked up and the fee is due.

**Recurrence.**
$V_k(y)=\min_{d_k\le x\le8}\bigl\{220x+500\,[x>y]+V_{k+1}(x)\bigr\}$, where
$[x>y]$ is 1 if $x>y$ and 0 otherwise.

**Boundary condition.** $V_5(y)=0$.

**Answer cell.** $V_1(0)$.

| Week $k$ | State $y$ | $V_k(y)$ | Decision $x$ |
|---|---|---|---|
| 4 | 7 | $1\,760+500=2\,260$ | 8 |
| 4 | 8 | $1\,760$ | 8 |
| 3 | 4, 5, 6 | $\min\{1\,540+500+2\,260,\ 1\,760+500+1\,760\}=4\,020$ | 8 |
| 3 | 7 | $\min\{1\,540+2\,260,\ 1\,760+500+1\,760\}=3\,800$ | 7 |
| 3 | 8 | $\min\{1\,540+2\,260,\ 1\,760+1\,760\}=3\,520$ | 8 |
| 2 | 7 | $\min\{880+4\,020,\ 1\,100+4\,020,\ 1\,320+4\,020,\ 1\,540+3\,800,\ 2\,260+3\,520\}=4\,900$ | 4 |
| 2 | 8 | $\min\{880+4\,020,\ \dots,\ 1\,540+3\,800,\ 1\,760+3\,520\}=4\,900$ | 4 |
| 1 | 0 | $\min\{1\,540+500+4\,900,\ 1\,760+500+4\,900\}=6\,940$ | 7 |

The minimum cost is **\$6 940**: rent 7 cars in week 1; return 3 and keep 4
in week 2; pick up 4 more in week 3, holding 8 although only 7 are needed; keep
the 8 cars in week 4. Check:
$(1\,540+500)+880+(1\,760+500)+1\,760=6\,940$. Holding the eighth car one
week early (\$220) avoids a second fee (\$500) in week 4. The next best plans
cost \$7 160.

*(An earlier version of this solution totals \$6 900; the correct minimum is
\$6 940. Its recurrence is also written with \$200 per car and a \$550 fee,
while the data, and its own tables, use \$220 and \$500.)*
:::

:::{prf:example} Selling sheep that multiply
:label: ex-08-07-a8-q8

A farmer owns $Y$ sheep. At the end of each year the farmer decides how
many sheep to sell and how many to keep; a sheep sold in year $k$ brings a
profit $p_k$, and every sheep kept becomes two by the end of the next year.
All sheep will be sold at the end of year $N$.

(a) Develop a DP model.

(b) Solve it for $N=3$, $Y=2$, $p_1=\$100$, $p_2=\$130$, $p_3=\$120$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** **Stages.** Years $k=1,\dots,N$ (decisions at the end of each year).
**States.** The number $y$ of sheep available at the end of year $k$; $y=Y$
in year 1, and since the flock at most doubles each year,
$0\le y\le2^{\,k-1}Y$. **Decisions.** The number $x\in\{0,\dots,y\}$ sold; the
$y-x$ sheep kept become $2(y-x)$. **Recurrence.**

$$
V_k(y)=\max_{x=0,\dots,y}\bigl\{p_kx+V_{k+1}\bigl(2(y-x)\bigr)\bigr\},\qquad k=1,\dots,N-1,
$$

**Boundary condition.** $V_N(y)=p_Ny$ (everything is sold).
**Answer cell.** $V_1(Y)$.

*(A state bound $y\le2^{\,k}Y$ that is sometimes given is twice too large;
it is harmless but wasteful.)*

**(b)** $V_3(y)=120y$ for $y=0,\dots,8$. At stage 2, for $y=0,\dots,4$,
$V_2(y)=\max_x\{130x+240(y-x)\}=240y$ with $x=0$ (keeping is worth
$2\cdot120=240$ per sheep, selling 130). At stage 1, with $y=2$:

$$
\begin{aligned}
V_1(2)&=\max\{0+V_2(4),\ 100+V_2(2),\ 200+V_2(0)\}\\
&=\max\{960,\ 580,\ 200\}=960\quad(x=0).
\end{aligned}
$$

Keep both sheep in year 1 (4 sheep in year 2), keep all 4 in year 2 (8 sheep
in year 3), and sell all 8 at the end of year 3 for $8\cdot\$120=\$960$.

**A closed form.** By induction $V_k(y)=m_ky$ with $m_N=p_N$ and
$m_k=\max\{p_k,\ 2m_{k+1}\}$: in each year it is optimal to sell all sheep if
$p_k\ge2m_{k+1}$ and to keep all otherwise. Here $m_3=120$, $m_2=240$,
$m_1=480$ and $V_1(2)=960$.
:::

(sec-08-07-lp-by-dp)=
## Linear, integer and nonlinear programs by DP

:::{prf:example} Three small programs solved by DP
:label: ex-08-07-a8-q9

Solve by DP, taking one variable per stage and the remaining right-hand sides
as the state:

(a) $\max\ 4x_1+14x_2$ s.t. $2x_1+7x_2\le21$, $7x_1+2x_2\le21$,
$x_1,x_2\ge0$;

(b) $\max\ 8x_1+7x_2$ s.t. $2x_1+x_2\le8$, $5x_1+2x_2\le15$, $x_1,x_2\ge0$
and integer;

(c) $\max\ 7x_1^2+6x_1+5x_2^2$ s.t. $x_1+2x_2\le10$, $x_1-3x_2\le9$,
$x_1,x_2\ge0$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Each problem has two constraints, so the state is the pair
$(\beta_1,\beta_2)$ of right-hand sides still available; a stage fixes one
variable and reduces both. The last stage is solved in closed form, and the
first stage becomes a one-variable problem.

**(a)** Stage 2 (variable $x_2$):
$V_2(\beta_1,\beta_2)=\max\{14x_2:\ 7x_2\le\beta_1,\ 2x_2\le\beta_2,\ x_2\ge0\}
=14\min\{\beta_1/7,\ \beta_2/2\}$. Stage 1 (variable $x_1\in[0,3]$):

$$
V_1(21,21)=\max_{0\le x_1\le3}\Bigl\{4x_1+14\min\Bigl\{\frac{21-2x_1}{7},\ \frac{21-7x_1}{2}\Bigr\}\Bigr\}.
$$

The first term of the minimum is the smaller one exactly when
$x_1\le7/3$. On $[0,7/3]$ the objective is $4x_1+2(21-2x_1)=42$; on
$[7/3,3]$ it is $4x_1+7(21-7x_1)=147-45x_1$, which decreases. Hence
$V_1(21,21)=42$ for **every** $x_1\in[0,7/3]$, with
$x_2=3-\tfrac27x_1$: the optimal solutions form the segment from $(0,3)$ to
$(\tfrac73,\tfrac73)$. The reason is that the objective
$4x_1+14x_2=2(2x_1+7x_2)$ is parallel to the first constraint. An LP solver
returns one end of the segment, $(0,3)$, with value 42.

**(b)** Stage 2: $V_2(\beta_1,\beta_2)=7\min\{\lfloor\beta_1\rfloor,\lfloor\beta_2/2\rfloor\}$
with $x_2$ equal to that minimum. Stage 1: $x_1\le\min\{\lfloor8/2\rfloor,\lfloor15/5\rfloor\}=3$ and

| $x_1$ | $\beta=(8-2x_1,\ 15-5x_1)$ | $x_2$ | $8x_1+7x_2$ |
|---|---|---|---|
| 0 | (8, 15) | $\min\{8,7\}=7$ | 49 |
| 1 | (6, 10) | $\min\{6,5\}=5$ | 43 |
| 2 | (4, 5) | $\min\{4,2\}=2$ | 30 |
| 3 | (2, 0) | $\min\{2,0\}=0$ | 24 |

The optimum is $(x_1,x_2)=(0,7)$ with value **49**; check:
$2\cdot0+7=7\le8$ and $5\cdot0+14=14\le15$.

**(c)** Here it is convenient to fix $x_2$ first (stage 1) and $x_1$ last
(stage 2). After $x_2$ is chosen, the right-hand sides left for $x_1$ are
$\beta_1=10-2x_2$ and $\beta_2=9+3x_2$, so $x_1\le m=\min\{\beta_1,\beta_2\}$.
Stage 2: $7x_1^2+6x_1$ increases for $x_1\ge0$, so
$V_2(\beta_1,\beta_2)=7m^2+6m$ with $x_1=m$. Stage 1 ($0\le x_2\le5$):

$$
\begin{aligned}
V_1&=\max_{0\le x_2\le5}\bigl\{5x_2^2+7m^2+6m\bigr\},\\
m&=\min\{10-2x_2,\ 9+3x_2\}.
\end{aligned}
$$

The two expressions for $m$ are equal at $x_2=1/5$. For $x_2\le1/5$,
$m=9+3x_2$ and the objective is $68x_2^2+396x_2+621$, increasing, with value
$702.92$ at $x_2=1/5$. For $x_2\ge1/5$, $m=10-2x_2$ and the objective is
$33x_2^2-292x_2+760$, a convex function whose maximum on $[1/5,5]$ is at an
end point: $702.92$ at $x_2=1/5$ and $125$ at $x_2=5$. Hence the optimum is
$x_2=1/5$, $x_1=m=48/5=9.6$, with value $17\,573/25=702.92$. It is the vertex
of the feasible region where both constraints are tight; the other vertices
give $0$ at $(0,0)$, $621$ at $(9,0)$ and $125$ at $(0,5)$, as expected for the
maximum of a convex function, which lies at a vertex.

*(An earlier version writes the stage-1 expression with $\beta-x_2$ in place of
$10-2x_2$ and $9+3x_2$; its expanded polynomials and the final answer are
correct.)*
:::

## References

- Dreyfus, S. E., & Law, A. M. (1977). *The art and theory of dynamic
  programming*. Academic Press. (Chapter 2: equipment replacement problems;
  Chapter 3: resource allocation problems.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 12: deterministic dynamic programming, problems of
  Sections 12.2–12.4.)
