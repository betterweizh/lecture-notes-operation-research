# DP Fundamentals and Shortest Paths

This chapter works through the examples and problems of
{doc}`05-01-dp-fundamentals` and of {doc}`05-04-dp-shortest-paths`. The first two
sections practice dynamic programming (DP) formulations: two staged decision problems, then eight
variants of the lattice path problem in which the cost structure forces a
richer state. The remaining sections follow {doc}`05-04-dp-shortest-paths`:
acyclic networks and topological order, successive approximation with
negative arcs, all-pairs recursions with their operation counts, and
label-setting methods with special objectives. The notation is that of the
two theory chapters: $k$ for the stage and $N$ for the number of stages,
$V$ for optimal values, $P$ for policies and predecessors, $c_{ij}$ for arc
lengths; in the shortest-path recursions the stage is written as a
superscript, $V^{(k)}(j)$. Every number was recomputed by
dynamic programming and checked by complete enumeration of the paths or with
a network library; every formulation of the lattice problems was also
checked against enumeration on random lattices.

(sec-08-06-formulations)=
## Principle of optimality and stage-wise formulations

:::{prf:example} Five-day hiking schedule
:label: ex-08-06-hike

Two hikers spend five days in a mountain area with three peaks, W, J and A.
The trail from the base of a mountain to its peak is 6 miles long for W,
4 miles for J and 5 miles for A, and each climb goes up and back down the
same trail. The bases are joined by trails of 3 miles (W–J), 2 miles (J–A)
and 5 miles (W–A). The hikers start on day 1 at the base of W and must be
back at the base of W at the end of day 5. Every day they climb exactly one
mountain, and in the evening they walk to the base of the mountain they will
climb the next day and camp there; they never climb the same mountain on
two consecutive days.

(a) Formulate a DP model that maximizes the total distance hiked and solve
it.

(b) Solve the problem again when W may not be climbed on day 5.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a) Model.**

- *Stage:* the day $k=1,\dots,5$ (more precisely, the start of day $k$).
- *State:* the base $i\in\{\mathrm W,\mathrm J,\mathrm A\}$ where the hikers
  are at the start of day $k$; it is the mountain climbed on day $k$.
- *Decision:* the base $j\ne i$ to walk to in the evening (the mountain of
  day $k+1$).
- *Data:* a climb of mountain $i$ is a round trip of $t_i$ miles,
  $t_{\mathrm W}=12$, $t_{\mathrm J}=8$, $t_{\mathrm A}=10$; $d_{ij}$ is the
  base-to-base distance, with $d_{\mathrm{WW}}=0$.
- *Value function:* $V_k(i)$ = the largest number of miles that can be hiked
  from the start of day $k$ at base $i$ to the end of day 5 at the base of
  W.
- *Recurrence:*

  $$
  V_k(i)=t_i+\max_{j\ne i}\bigl\{d_{ij}+V_{k+1}(j)\bigr\}\qquad(k=1,\dots,4).
  $$

- *Boundary condition:* on day 5 the hikers climb mountain $i$ and walk back
  to the base of W: $V_5(i)=t_i+d_{i\mathrm W}$, so
  $V_5(\mathrm W)=12$, $V_5(\mathrm J)=11$, $V_5(\mathrm A)=15$.
- *Answer cell:* $V_1(\mathrm W)$.

**Solution.** The table gives $V_k(i)$ with the best next base in brackets.
For example
$V_4(\mathrm W)=12+\max\{3+V_5(\mathrm J),\ 5+V_5(\mathrm A)\}=12+\max\{14,20\}=32$.

| Day $k$ | $V_k(\mathrm W)$ | $V_k(\mathrm J)$ | $V_k(\mathrm A)$ |
|---|---|---|---|
| 5 | 12 | 11 | 15 |
| 4 | 32 (A) | 25 (A) | 27 (W) |
| 3 | 44 (A) | 43 (W) | 47 (W) |
| 2 | 64 (A) | 57 (A) | 59 (W) |
| 1 | 76 (A) | – | – |

The optimum is $V_1(\mathrm W)=76$ miles with the schedule W–A–W–A–W:
climbs $12+10+12+10+12=56$ miles and four walks between W and A,
$4\times5=20$ miles. Enumerating all 16 admissible schedules confirms that
it is the only one with 76 miles.

**(b)** Now $V_5(\mathrm W)=-\infty$ (forbidden), and the recursion gives

| Day $k$ | $V_k(\mathrm W)$ | $V_k(\mathrm J)$ | $V_k(\mathrm A)$ |
|---|---|---|---|
| 5 | $-\infty$ | 11 | 15 |
| 4 | 32 (A) | 25 (A) | 23 (J) |
| 3 | 40 (J or A) | 43 (W) | 47 (W) |
| 2 | 64 (A) | 57 (A) | 55 (W or J) |
| 1 | 72 (J or A) | – | – |

The optimum drops to 72 miles, and the ties give four optimal schedules:
W–A–W–A–J, W–A–W–J–A, W–J–A–W–A and W–A–J–W–A.

*(The rule must forbid repeating a mountain on consecutive days, as the
recurrence $j\ne i$ assumes. A boundary condition stated at "day 6" with
$V_6(\mathrm W)=0$ and $V_6(\mathrm J)=V_6(\mathrm A)=-\infty$ does not
work together with $j\ne i$, because it would forbid returning to W after a
climb of W on day 5; the day-5 boundary used here avoids the problem.)*
:::

:::{prf:example} Shortest route through seven cities
:label: ex-08-06-staged-route

A road network joins seven cities; its arcs and their lengths (km) are $1\to2$ 5, $1\to3$ 14, $1\to4$ 6,
$2\to3$ 6, $2\to5$ 12, $2\to6$ 12, $3\to5$ 2, $3\to6$ 3, $4\to3$ 7,
$4\to6$ 9, $5\to7$ 4 and $6\to7$ 4. Formulate a DP model for the shortest
route from city 1 to city 7 and solve it.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

```{figure} figures/08-06-staged-route.svg
:name: fig-08-06-staged-route
:width: 720px
:class: dark-light
:alt: Left, a seven-node road network with distances and the route 1-2-3-5-7 highlighted; right, the same network drawn in five stages with a dashed zero-length arc from city 7 to city 7

(a) The road network of {prf:ref}`ex-08-06-staged-route` with the shortest
route 1–2–3–5–7 (17 km). (b) The same network drawn in stages: stage $k$
holds the cities that can be reached with $k-1$ arcs and that still lead to
city 7; the dashed arc of length 0 keeps a route that has already arrived
at city 7.
```

**Model.** The network of {numref}`fig-08-06-staged-route`(a) is not drawn in columns, so a stage variable must be made.
Take as stage the number of arcs traveled plus one ({numref}`fig-08-06-staged-route`(b)):
city 1 is the only state of stage 1, cities 2, 3, 4 form stage 2, cities
3, 5, 6 stage 3, cities 5, 6, 7 stage 4 and city 7 stage 5. A route with
three arcs reaches city 7 at stage 4; the dummy arc $7\to7$ of length 0
carries it to stage 5, so that every route ends in the same state.

- *State:* city $i$ at stage $k$; *decision:* the next city $j$.
- *Value function:* $V_k(i)$ = the length of a shortest route from city
  $i$ at stage $k$ to city 7 at stage 5.
- *Recurrence:* $V_k(i)=\min_{j}\{c_{ij}+V_{k+1}(j)\}$, the minimum over the
  arcs $(i,j)$ between stage $k$ and stage $k+1$ (including $c_{77}=0$).
- *Boundary condition:* $V_5(7)=0$; *answer cell:* $V_1(1)$.

**Solution.**

| Stage | Values (best next city) |
|---|---|
| 4 | $V_4(5)=4$ (7), $V_4(6)=4$ (7), $V_4(7)=0$ (7) |
| 3 | $V_3(3)=\min\{2+4,\ 3+4\}=6$ (5), $V_3(5)=4$ (7), $V_3(6)=4$ (7) |
| 2 | $V_2(2)=\min\{6+6,\ 12+4,\ 12+4\}=12$ (3), $V_2(3)=6$ (5), $V_2(4)=\min\{7+6,\ 9+4\}=13$ (3 or 6) |
| 1 | $V_1(1)=\min\{5+12,\ 14+6,\ 6+13\}=17$ (2) |

The shortest route is 1–2–3–5–7 with $5+6+2+4=17$ km. Of the nine routes
from 1 to 7, the next best is 1–2–3–6–7 with 18 km.

**Another stage variable.** The network is acyclic, so the one-pass
recursion of {doc}`05-04-dp-shortest-paths` also applies, but the city numbers
are not a topological order (the arc $4\to3$ points backwards). The order
1, 2, 4, 3, 5, 6, 7 is topological, and in it
$V(2)=5$, $V(4)=6$, $V(3)=\min\{14,\ 5+6,\ 6+7\}=11$, $V(5)=\min\{17,\ 11+2\}=13$,
$V(6)=\min\{17,\ 11+3,\ 15\}=14$ and $V(7)=\min\{13+4,\ 14+4\}=17$, the
same answer without the dummy arc.
:::

(sec-08-06-lattice)=
## Lattice variants: turns, coupons, minimax and products

The eight problems of this section use the lattice network of
{doc}`05-01-dp-fundamentals` with $N$ stages: vertex $(k,y)$ is in column $k$
at level $y$, $A=(0,0)$, $B=(N,0)$, and from each vertex an up arc of cost
$c^{\uparrow}(k,y)$ leads to $(k+1,y+1)$ and a down arc of cost
$c^{\downarrow}(k,y)$ leads to $(k+1,y-1)$; a cost $\infty$ means that the
arc does not exist. Every path from $A$ to $B$ has exactly $N$ arcs. A
**turn** at a vertex is a change from an up arc to a down arc or from a down
arc to an up arc. In each problem, formulate a DP (stage, state, decision,
value function, recurrence, boundary conditions, answer cell) for the
cheapest path from $A$ to $B$ under the stated cost structure. Each
solution ends with a numerical check on the 16-vertex network of
{numref}`fig-05-01-diamond-network` ($N=6$), whose 20 paths can be enumerated.

As in {doc}`05-01-dp-fundamentals`, the direction variable
$\omega\in\{\mathrm U,\mathrm D\}$ is the direction in which a vertex is
**left**, and every value function is defined with that convention.

:::{prf:example} Growing turn penalties
:label: ex-08-06-lattice-a

The cost of a path is the sum of its arc costs plus a turn penalty: the
first turn costs 1, the second 2, the $j$th turn costs $j$. (A path with
three turns pays $1+2+3=6$.) Give a backward and a forward formulation.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The penalty of the next turn depends on how many turns were made before,
so the state must count them.

**Backward formulation.**

- *Stage* $k=0,\dots,N-1$; *state* $(y,\omega,\ell)$: the vertex level, the
  direction in which the vertex is left, and the number $\ell$ of turns
  made at the vertices before $(k,y)$.
- *Value function:* $V_k(y,\omega,\ell)$ = least cost from $(k,y)$ to $B$
  when $(k,y)$ is left in direction $\omega$ and $\ell$ turns have been
  made before $(k,y)$.
- *Decision:* the direction in which the next vertex is left (a turn there
  is the $(\ell+1)$th turn and costs $\ell+1$).
- *Recurrence* ($k=N-2,\dots,0$; $\ell=0,\dots,k-1$, or $\ell=0$ at $k=0$):

  $$
  \begin{aligned}
  V_k(y,\mathrm U,\ell)&=c^{\uparrow}(k,y)+\min\bigl\{V_{k+1}(y+1,\mathrm U,\ell),\ (\ell+1)+V_{k+1}(y+1,\mathrm D,\ell+1)\bigr\},\\
  V_k(y,\mathrm D,\ell)&=c^{\downarrow}(k,y)+\min\bigl\{(\ell+1)+V_{k+1}(y-1,\mathrm U,\ell+1),\ V_{k+1}(y-1,\mathrm D,\ell)\bigr\}.
  \end{aligned}
  $$

- *Boundary conditions:* the last arc enters $B$, where nothing more is
  paid: $V_{N-1}(y,\omega,\ell)=c^{\omega}(N-1,y)$ for every $\ell$, where
  $c^{\mathrm U}=c^{\uparrow}$ and $c^{\mathrm D}=c^{\downarrow}$.
- *Answer cell:* $\min\{V_0(0,\mathrm U,0),\ V_0(0,\mathrm D,0)\}$.

**Forward formulation.** Let $F_k(y,\omega,\ell)$ be the least cost from
$A$ to $(k,y)$ when $(k,y)$ is *entered* in direction $\omega$ and $\ell$
turns were made at the vertices before $(k,y)$. A turn at the previous
vertex $(k-1,y\mp1)$ is then the $\ell$th turn and costs $\ell$:

$$
\begin{aligned}
F_k(y,\mathrm U,\ell)&=c^{\uparrow}(k-1,y-1)+\min\bigl\{F_{k-1}(y-1,\mathrm U,\ell),\ \ell+F_{k-1}(y-1,\mathrm D,\ell-1)\bigr\},\\
F_k(y,\mathrm D,\ell)&=c^{\downarrow}(k-1,y+1)+\min\bigl\{F_{k-1}(y+1,\mathrm D,\ell),\ \ell+F_{k-1}(y+1,\mathrm U,\ell-1)\bigr\},
\end{aligned}
$$

for $k=2,\dots,N$ (a term with $\ell-1<0$ is $\infty$), with boundary
conditions $F_1(1,\mathrm U,0)=c^{\uparrow}(0,0)$,
$F_1(-1,\mathrm D,0)=c^{\downarrow}(0,0)$ (all other $F_1=\infty$) and answer
$\min_{\omega,\ \ell}F_N(0,\omega,\ell)$, the minimum over both directions and
$\ell=0,\dots,N-1$.

Both formulations need four arguments (stage, level, direction, number of
turns); the number of states is about $N$ times that of the plain turn
penalty of {prf:ref}`ex-05-01-turn-penalty`.

**Check.** On the network of {numref}`fig-05-01-diamond-network` both
recursions give 15, on $A$–$D$–$G$–$K$–$N$–$Q$–$B$: travel time 14 and one
turn, which costs 1. The path that is fastest without penalties,
$A$–$C$–$F$–$J$–$M$–$P$–$B$, turns three times and costs
$13+1+2+3=19$. Enumeration of the 20 paths agrees.

*(A formulation that defines the direction as the direction of **arrival**
at $(k,y)$ but uses the recurrence above, with $c^{\uparrow}$ in both rows
of $V_k(y,\mathrm U,\ell)$, mixes the two conventions: the recurrence is
correct only for the leaving direction. With the arrival convention the
penalty of a turn at $(k,y)$ itself must be charged when $(k,y)$ is left in
the other direction.)*
:::

:::{prf:example} At most a given number of turns
:label: ex-08-06-lattice-b

Turns cost nothing, but only paths with at most $L$ turns are allowed.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The state must count the turns made so far, and a value function with too
many turns must be infeasible.

- *Stage* $k$; *state* $(y,\omega,\ell)$ with $\ell=0,\dots,L$ the number
  of turns made before $(k,y)$; *decision:* the leaving direction at the
  next vertex.
- *Value function:* $V_k(y,\omega,\ell)$ = least cost from $(k,y)$ to $B$,
  leaving $(k,y)$ in direction $\omega$, when $\ell$ turns have been made
  before $(k,y)$ and at most $L-\ell$ more are allowed.
- *Recurrence* ($k=N-2,\dots,0$):

  $$
  \begin{aligned}
  V_k(y,\mathrm U,\ell)&=c^{\uparrow}(k,y)+\min\bigl\{V_{k+1}(y+1,\mathrm U,\ell),\ V_{k+1}(y+1,\mathrm D,\ell+1)\bigr\},\\
  V_k(y,\mathrm D,\ell)&=c^{\downarrow}(k,y)+\min\bigl\{V_{k+1}(y-1,\mathrm U,\ell+1),\ V_{k+1}(y-1,\mathrm D,\ell)\bigr\},
  \end{aligned}
  $$

  with the convention $V_{k+1}(\cdot,\cdot,L+1)=\infty$, which removes the
  turning branch when $\ell=L$.
- *Boundary conditions:* $V_{N-1}(y,\omega,\ell)=c^{\omega}(N-1,y)$ for
  $\ell=0,\dots,L$.
- *Answer cell:* $\min\{V_0(0,\mathrm U,0),\ V_0(0,\mathrm D,0)\}$.

**Check.** With $L=0$ the answer is $\infty$: no path from $A$ to $B$ is
straight. With $L=1$ only the two "bent" paths along the edges of the
diamond are allowed, and the answer is 14 on $A$–$D$–$G$–$K$–$N$–$Q$–$B$
(the other one, along the top, takes 18). With $L=2$ the answer is still 14,
now attained by three paths: $A$–$D$–$G$–$K$–$N$–$Q$–$B$,
$A$–$C$–$F$–$J$–$N$–$Q$–$B$ and $A$–$D$–$G$–$J$–$M$–$P$–$B$. With $L=3$
the unrestricted optimum 13 becomes feasible.
:::

:::{prf:example} The cheapest $N-1$ arcs
:label: ex-08-06-lattice-c

The cost of a path is the sum of its $N-1$ smallest arc costs.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Every path has exactly $N$ arcs, so the sum of its $N-1$ smallest arc
costs is its total cost minus its largest arc cost. Equivalently, the
traveler holds one coupon that makes one arc of his choice free, and he
chooses the arc optimally; the state records whether the coupon is still
held.

- *Stage* $k$; *state* $(y,q)$ with $q\in\{0,1\}$ the number of coupons
  still held; *decision:* the direction, and (if $q=1$) whether to use the
  coupon on this arc.
- *Value function:* $V_k(y,q)$ = least cost from $(k,y)$ to $B$ holding $q$
  coupons.
- *Recurrence:*

  $$
  \begin{aligned}
  V_k(y,0)&=\min\bigl\{c^{\uparrow}(k,y)+V_{k+1}(y+1,0),\ c^{\downarrow}(k,y)+V_{k+1}(y-1,0)\bigr\},\\
  V_k(y,1)&=\min\bigl\{c^{\uparrow}(k,y)+V_{k+1}(y+1,1),\ V_{k+1}(y+1,0),\\
  &\qquad\quad c^{\downarrow}(k,y)+V_{k+1}(y-1,1),\ V_{k+1}(y-1,0)\bigr\},
  \end{aligned}
  $$

  where the second and fourth terms use the coupon on the arc (an arc that
  does not exist cannot be used, with or without a coupon).
- *Boundary conditions:* $V_N(0,0)=0$ and $V_N(0,1)=\infty$: the coupon
  must have been used, so that exactly $N-1$ arcs are paid.
- *Answer cell:* $V_0(0,1)$.

Since the minimum over the arc that gets the coupon is the same as dropping
the largest arc, the DP finds $\min_{\text{paths}}(\text{total}-\text{largest arc})$.

**Check.** The answer is 8, on $A$–$D$–$F$–$J$–$M$–$P$–$B$ with arcs
$0,7,2,2,2,2$: total 15, the arc $D$–$F$ of cost 7 is free. The path that is
fastest in total (13) has largest arc 4 and costs 9 here.

*(When all arc costs are nonnegative, the boundary condition $V_N(0,1)=0$,
which lets the coupon go unused, gives the same answer, because using the
coupon never hurts. With negative costs it does not: if every arc of a path
is negative, dropping the largest one increases the sum, and only
$V_N(0,1)=\infty$ enforces "exactly $N-1$ arcs".)*
:::

:::{prf:example} Minimize the largest arc
:label: ex-08-06-lattice-d

The cost of a path is its largest arc cost.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

- *Stage* $k$; *state* $y$; *decision* U or D.
- *Value function:* $V_k(y)$ = the smallest possible largest arc cost on a
  path from $(k,y)$ to $B$.
- *Recurrence:*

  $$
  V_k(y)=\min\Bigl\{\max\bigl\{c^{\uparrow}(k,y),\ V_{k+1}(y+1)\bigr\},\ \max\bigl\{c^{\downarrow}(k,y),\ V_{k+1}(y-1)\bigr\}\Bigr\}.
  $$

  The cost of going up is the larger of the arc met now and the best
  possible largest arc from the next vertex on.
- *Boundary condition:* $V_N(0)=-\infty$, the largest element of an empty
  set of arcs. When all costs are nonnegative, $V_N(0)=0$ gives the same
  values.
- *Answer cell:* $V_0(0)$.

The principle of optimality holds because the cost is combined by
$\max\{c,\cdot\}$, which is nondecreasing in the tail value
({prf:ref}`thm-05-01-principle`, monotonicity). Minimax objectives are
common: a storage facility must be as large as the largest stock it will
hold, a pipe as wide as the largest flow, and a controller may minimize the
largest deviation from a target trajectory.

**Check.** The answer is 4, attained by four paths:
$A$–$C$–$F$–$I$–$M$–$P$–$B$, $A$–$C$–$F$–$J$–$M$–$P$–$B$,
$A$–$C$–$F$–$J$–$N$–$Q$–$B$ and $A$–$D$–$G$–$K$–$N$–$Q$–$B$. No path avoids
an arc of cost 4 or more: from $A$ one must use $C$–$E$ (5), $C$–$F$ (4),
$D$–$F$ (7) or $D$–$G$ (3) followed by $G$–$J$ (5) or $G$–$K$ (4).
:::

:::{prf:example} Coupons that reduce arc costs
:label: ex-08-06-lattice-e

The traveler starts at $A$ with $Q$ coupons. Spending $\nu$ coupons on an
up arc from $(k,y)$ reduces its cost by $\gamma^{\uparrow}(k,y,\nu)$, and
spending $\nu$ coupons on a down arc reduces its cost by
$\gamma^{\downarrow}(k,y,\nu)$; both reductions increase with $\nu$, and
$\gamma(\cdot,\cdot,0)=0$. The cost of a path is the sum of its reduced arc
costs.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

- *Stage* $k$; *state* $(y,q)$, $q=0,\dots,Q$ the coupons still held;
  *decision:* the direction $\omega$ and the number $\nu\in\{0,\dots,q\}$ of
  coupons spent on the arc.
- *Value function:* $V_k(y,q)$ = least cost from $(k,y)$ to $B$ with $q$
  coupons in hand.
- *Recurrence:*

  $$
  \begin{aligned}
  V_k(y,q)=\min\Bigl\{&\min_{\nu=0,\dots,q}\bigl\{c^{\uparrow}(k,y)-\gamma^{\uparrow}(k,y,\nu)+V_{k+1}(y+1,q-\nu)\bigr\},\\
  &\min_{\nu=0,\dots,q}\bigl\{c^{\downarrow}(k,y)-\gamma^{\downarrow}(k,y,\nu)+V_{k+1}(y-1,q-\nu)\bigr\}\Bigr\}.
  \end{aligned}
  $$

- *Boundary conditions:* $V_N(0,q)=0$ for $q=0,\dots,Q$ (coupons left over
  are worthless).
- *Answer cell:* $V_0(0,Q)$.

Each state now needs up to $2(q+1)$ evaluations, so the work grows like
$Q^2$ per vertex: a second dimension of the state and a decision set that
grows with it.

**Check.** Take $Q=2$ and the same reductions on every arc: one coupon
saves 2, two coupons on the same arc save 3, and an arc can never cost less
than 0 (the reduction is at most the arc cost). The DP gives 9: on the
path $A$–$C$–$F$–$J$–$M$–$P$–$B$ (13), one coupon on each of two arcs of cost
at least 2 saves 4 (ten such allocations tie). Putting both coupons on one
arc saves only 3, and no other path reaches 9. The assumption
$\gamma(\cdot,\cdot,0)=0$ matters: it makes "spend nothing" the plain arc
cost.
:::

:::{prf:example} No turns at two consecutive vertices
:label: ex-08-06-lattice-f

The cost is the usual sum of arc costs, but a path that turns at a vertex
may not turn at the next vertex. Give two formulations: one with a flag in
the state, and one in which the stage advances by two after a turn.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Formulation with a flag.**

- *State* $(y,\omega,\beta)$: level, leaving direction, and $\beta=1$ if
  the path turned *at* $(k,y)$ (so it may not turn at the next vertex),
  $\beta=0$ otherwise.
- *Value function:* $V_k(y,\omega,\beta)$ = least cost from $(k,y)$ to $B$,
  leaving in direction $\omega$, with turn flag $\beta$.
- *Recurrence:*

  $$
  \begin{aligned}
  V_k(y,\mathrm U,\beta)&=c^{\uparrow}(k,y)+\min\bigl\{V_{k+1}(y+1,\mathrm U,0),\ V_{k+1}(y+1,\mathrm D,1)\ \text{if }\beta=0\bigr\},\\
  V_k(y,\mathrm D,\beta)&=c^{\downarrow}(k,y)+\min\bigl\{V_{k+1}(y-1,\mathrm D,0),\ V_{k+1}(y-1,\mathrm U,1)\ \text{if }\beta=0\bigr\}.
  \end{aligned}
  $$

- *Boundary conditions:* $V_{N-1}(y,\omega,\beta)=c^{\omega}(N-1,y)$;
  *answer cell:* $\min_{\omega}V_0(0,\omega,0)$ (leaving $A$ is not a turn).

**Formulation with a two-stage jump.** If the path turns at $(k+1,\cdot)$,
it must go straight at $(k+2,\cdot)$, so both arcs can be taken in one
decision. Let $W_k(y,\omega)$ be the least cost from $(k,y)$ to $B$ when
$(k,y)$ is left in direction $\omega$ and a turn is allowed at the next
vertex. Then

$$
\begin{aligned}
W_k(y,\mathrm U)&=c^{\uparrow}(k,y)+\min\bigl\{W_{k+1}(y+1,\mathrm U),\ c^{\downarrow}(k+1,y+1)+W_{k+2}(y,\mathrm D)\bigr\},\\
W_k(y,\mathrm D)&=c^{\downarrow}(k,y)+\min\bigl\{W_{k+1}(y-1,\mathrm D),\ c^{\uparrow}(k+1,y-1)+W_{k+2}(y,\mathrm U)\bigr\},
\end{aligned}
$$

with boundary conditions $W_N(0,\mathrm U)=W_N(0,\mathrm D)=0$ (the
direction at $B$ is a dummy) and $c^{\uparrow}(N,\cdot)=c^{\downarrow}(N,\cdot)=\infty$,
and answer $\min_\omega W_0(0,\omega)$. The second branch moves the stage
from $k$ to $k+2$; the backward computation still works because $W_{k+2}$
is known when $W_k$ is computed. In fact $W_k(y,\omega)=V_k(y,\omega,0)$:
the jump simply substitutes the forced straight move of state $\beta=1$.

**Check.** Both formulations give 13 on the path
$A$–$C$–$F$–$J$–$M$–$P$–$B$, which turns at $C$, $J$ and $P$ (stages 1, 3
and 5, never consecutive), so the restriction costs nothing here.

*(A version of the jump formulation that records the direction of
**arrival** and, when it turns at $(k,y)$ itself, jumps from stage $k$ to
stage $k+2$, must treat $k=N-1$ separately: there a turn leads straight
into $B$ and no second arc exists. Without this special case every path
whose last turn is at the vertex before $B$ is lost; on the check network
that version returns 14 instead of 13, because the optimal path turns at
$P$. In random tests it was wrong in 88 of 200 lattices. The flag
formulation also needs the value function for $\beta=1$, written out
above.)*
:::

:::{prf:example} Product of positive costs
:label: ex-08-06-lattice-g

All arc costs are positive, and the cost of a path is the product of its
arc costs. (Such costs arise as failure or survival factors in fault trees.)
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

- *Stage* $k$; *state* $y$; *value function* $V_k(y)$ = least product of
  arc costs from $(k,y)$ to $B$.
- *Recurrence:*

  $$
  V_k(y)=\min\bigl\{c^{\uparrow}(k,y)\cdot V_{k+1}(y+1),\ c^{\downarrow}(k,y)\cdot V_{k+1}(y-1)\bigr\}.
  $$

- *Boundary condition:* $V_N(0)=1$, the empty product (not 0).
- *Answer cell:* $V_0(0)$.

Multiplication by a positive number preserves order, so the combination is
monotone and the principle of optimality holds. Products can be handled as
easily as sums; alternatively, since $\ln$ is increasing and
$\ln\prod c=\sum\ln c$, minimizing the product is the ordinary shortest-path
problem with arc costs $\ln c^{\uparrow}$, $\ln c^{\downarrow}$.

**Check.** The network of {numref}`fig-05-01-diamond-network` has an arc of
cost 0, so add 1 to every arc cost. The least product is then 600, on
$A$–$D$–$G$–$K$–$N$–$Q$–$B$ with factors $1\cdot4\cdot5\cdot3\cdot5\cdot2$.
:::

:::{prf:example} Product of signed costs
:label: ex-08-06-lattice-h

Some arc costs are positive and some negative; the cost of a path is the
product of its arc costs.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Multiplication by a negative number reverses order: if the product of the
arcs up to a vertex is negative, the *largest* product of the remaining
arcs gives the smallest total. So the minimum alone is not enough, and the
principle of optimality fails for a value function that records only the
least product. Record both extremes.

- *Value functions:* $V^{\min}_k(y)$ and $V^{\max}_k(y)$ = the least and
  the greatest product of arc costs over the paths from $(k,y)$ to $B$.
  (Equivalently, a state variable records whether the remaining product is
  to be minimized or maximized.)
- *Recurrence:* for each existing arc from $(k,y)$ with cost $c$ to the next
  vertex $(k+1,y')$, its candidates are

  $$
  \begin{aligned}
  \text{for the minimum: }&\begin{cases}c\,V^{\min}_{k+1}(y'),&c\ge0,\\ c\,V^{\max}_{k+1}(y'),&c<0,\end{cases}\\[4pt]
  \text{for the maximum: }&\begin{cases}c\,V^{\max}_{k+1}(y'),&c\ge0,\\ c\,V^{\min}_{k+1}(y'),&c<0;\end{cases}
  \end{aligned}
  $$

  $V^{\min}_k(y)$ is the least of the candidates of the (one or two) arcs
  and $V^{\max}_k(y)$ the greatest. Arcs that do not exist are left out
  (here $\infty$ cannot be used, because $c\cdot\infty$ has no meaning when
  $c$ is negative or zero).
- *Boundary conditions:* $V^{\min}_N(0)=V^{\max}_N(0)=1$.
- *Answer cell:* $V^{\min}_0(0)$.

**Check.** Subtract 3 from every arc cost of
{numref}`fig-05-01-diamond-network`, which gives costs between $-3$ and 5.
The DP gives $V^{\min}_0(0)=-240$, on $A$–$D$–$F$–$I$–$M$–$Q$–$B$ with
factors $(-3)\cdot4\cdot(-2)\cdot1\cdot5\cdot(-2)$, and
$V^{\max}_0(0)=120$. Enumeration of the 20 paths gives the same extremes.
A recursion that kept only $V^{\min}$ would discard, at the vertex $I$,
the largest continuation that the negative arcs before it need.
:::

(sec-08-06-acyclic)=
## Acyclic networks and topological order

:::{prf:example} Every acyclic network has a source node
:label: ex-08-06-source-node

Show that every acyclic network (with at least one node) has a node with
only outgoing arcs, that is, a node that no arc enters.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Among all paths of the network choose one with the largest number of arcs;
there are finitely many paths (a path visits each node at most once), and a
single node counts as a path with no arcs. Let $s$ be its first node.
Suppose an arc $(u,s)$ enters $s$. If $u$ is not on the path, then $u$
followed by the path is a longer path, contradicting the choice. If $u$ is
on the path, then the path from $s$ to $u$ followed by the arc $(u,s)$ is a
cycle, contradicting acyclicity. So no arc enters $s$.

(A second proof, walking backwards along incoming arcs until the walk must
stop, is given in {prf:ref}`lem-05-04-source-node`.) The statement is the base
of the topological ordering: repeatedly removing a source node numbers all
nodes so that every arc points from a lower to a higher number
({prf:ref}`thm-05-04-topological-order`).
:::

:::{prf:example} Numbering the nodes of an acyclic network
:label: ex-08-06-topological-order

Number the nodes of the network of {numref}`fig-08-06-acyclic-network` so that
$i<j$ for every arc $(i,j)$.
:::

```{figure} figures/08-06-acyclic-network.svg
:name: fig-08-06-acyclic-network
:width: 380px
:class: dark-light
:alt: Six-node acyclic network with nodes named a to f and arc lengths

An acyclic network with unnumbered nodes $a,\dots,f$ and arc lengths, used
in {prf:ref}`ex-08-06-topological-order`, {prf:ref}`ex-08-06-acyclic-forward` and
{prf:ref}`ex-08-06-acyclic-backward`.
```

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The arcs are $e\to a$, $e\to d$, $e\to f$, $a\to d$, $a\to b$, $d\to b$,
$d\to f$, $d\to c$, $f\to c$ and $c\to b$. Remove source nodes one at a
time:

| Step | Source nodes of the remaining network | Numbered |
|---|---|---|
| 1 | $e$ | $e=1$ |
| 2 | $a$ ($d$ is still entered by $a\to d$) | $a=2$ |
| 3 | $d$ | $d=3$ |
| 4 | $f$ ($c$ is still entered by $f\to c$) | $f=4$ |
| 5 | $c$ | $c=5$ |
| 6 | $b$ | $b=6$ |

The numbering $e,a,d,f,c,b\mapsto1,\dots,6$ turns the arcs into
$1\to2$, $1\to3$, $1\to4$, $2\to3$, $2\to6$, $3\to6$, $3\to4$, $3\to5$,
$4\to5$, $5\to6$: every arc goes from a lower to a higher number. At every
step exactly one source node was available, so this is the only valid
numbering of this network.
:::

:::{prf:example} One-pass forward recursion
:label: ex-08-06-acyclic-forward

Use the forward recursion {eq}`eq-05-04-acyclic-forward` to find the shortest
distances from node $e$ to all nodes of {numref}`fig-08-06-acyclic-network`.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

In the topological order of {prf:ref}`ex-08-06-topological-order`:

$$
\begin{aligned}
V(e)&=0,\\
V(a)&=V(e)+4=4\quad(P=e),\\
V(d)&=\min\{V(e)+3,\ V(a)+2\}=\min\{3,\ 6\}=3\quad(P=e),\\
V(f)&=\min\{V(e)+5,\ V(d)+1\}=\min\{5,\ 4\}=4\quad(P=d),\\
V(c)&=\min\{V(d)+5,\ V(f)+2\}=\min\{8,\ 6\}=6\quad(P=f),\\
V(b)&=\min\{V(a)+6,\ V(d)+7,\ V(c)+3\}=\min\{10,\ 10,\ 9\}=9\quad(P=c).
\end{aligned}
$$

The shortest path to $b$ is $e$–$d$–$f$–$c$–$b$ with length
$3+1+2+3=9$. The recursion used one addition per arc (10) and 5
comparisons; the general bound {eq}`eq-05-04-acyclic-count` for six nodes, 15
additions and 10 comparisons, assumes arcs between all pairs.
:::

:::{prf:example} Backward recursion for an acyclic network
:label: ex-08-06-acyclic-backward

Give a backward DP formulation for the shortest distances from every node
to node $n$ of an acyclic network, prove it, and apply it to
{numref}`fig-08-06-acyclic-network` with $n=b$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Formulation.** Number the nodes topologically. Let $\bar V(i)$ be the
length of a shortest path from node $i$ to node $n$, and $\bar P(i)$ the
node after $i$ on it. Then

$$
\bar V(i)=\min_{j>i}\bigl\{c_{ij}+\bar V(j)\bigr\}\qquad(i=n-1,n-2,\dots,1),
$$

with $c_{ij}=\infty$ for missing arcs, boundary condition $\bar V(n)=0$,
and answers $\bar V(i)$ for all $i$.

**Proof.** By downward induction on $i$. Node $n$ has no outgoing arc to a
higher node, so the only path from $n$ to $n$ is the empty one. Assume
$\bar V(i+1),\dots,\bar V(n)$ correct. Every path from $i$ to $n$ leaves $i$
by an arc $(i,j)$ with $j>i$, and its remainder is a path from $j$ to $n$,
at best of length $\bar V(j)$; choosing $j$ optimally gives the recurrence.
The work is the same as for the forward recursion.

**Application.** In the order $e,a,d,f,c,b$:

$$
\begin{aligned}
\bar V(b)&=0,\\
\bar V(c)&=3,\\
\bar V(f)&=2+\bar V(c)=5,\\
\bar V(d)&=\min\{7+0,\ 1+5,\ 5+3\}=6\quad(\bar P=f),\\
\bar V(a)&=\min\{6+0,\ 2+6\}=6\quad(\bar P=b),\\
\bar V(e)&=\min\{4+6,\ 3+6,\ 5+5\}=9\quad(\bar P=d).
\end{aligned}
$$

$\bar V(e)=9=V(b)$, as it must. The sums $V(j)+\bar V(j)$ equal 9 at
$e,d,f,c,b$ (the nodes of the shortest path) and 10 at $a$.
:::

(sec-08-06-successive)=
## Successive approximation and negative cycles

The problems of this section use procedure 2 ({eq}`eq-05-04-procedure-2`) and
procedure 3 ({eq}`eq-05-04-procedure-3`) of {doc}`05-04-dp-shortest-paths`. Along
a walk $1,i_1,\dots,i_m,j$, a *decrease* occurs at a node whose number is
smaller than that of the node before it (the end node $j$ included), and a
*reversal* occurs at an intermediate node where the node numbers change from
increasing to decreasing or back. As in procedure 1, values are lengths of
shortest *walks* of the stated kind.

:::{prf:example} Correctness of procedure 2
:label: ex-08-06-proc2-proof

Show that $\tilde V^{(k)}(j)$, the length of a shortest walk from node 1 to
node $j$ with at most $k-1$ decreases, is given by the recurrence
{eq}`eq-05-04-procedure-2` and the boundary conditions
$\tilde V^{(0)}(1)=0$, $\tilde V^{(0)}(j)=\infty$ ($j\ne1$).
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Read $\tilde V^{(0)}$ as "at most $-1$ decreases": only the empty walk at
node 1 qualifies, which is the boundary condition. Fix $k\ge1$, assume that
$\tilde V^{(k-1)}$ is correct, and proceed by induction on $j=1,\dots,n$,
assuming $\tilde V^{(k)}(i)$ correct for $i<j$.

A walk to $j$ with at most $k-1$ decreases is either the empty walk (only
when $j=1$), or ends with an arc $(i,j)$, $i\ne j$.

- If $i<j$, the last arc is not a decrease, so the walk up to $i$ has at
  most $k-1$ decreases; its best length is $\tilde V^{(k)}(i)$, known
  because $i<j$. Conversely, appending an arc $(i,j)$ with $i<j$ to any
  walk with at most $k-1$ decreases gives such a walk to $j$.
- If $i>j$, the last arc is a decrease, so the walk up to $i$ has at most
  $k-2$ decreases; its best length is $\tilde V^{(k-1)}(i)$, and appending
  $(i,j)$ adds exactly one decrease.

Choosing the best $i$ in each group gives {eq}`eq-05-04-procedure-2`. For
$j=1$ there is no $i<1$, and the empty walk is represented by the term
$\tilde V^{(k-1)}(1)$, which is 0 as long as no negative cycle passes
through node 1. Because the nodes are processed in increasing order, every
value on the right-hand side is known when it is needed.
:::

:::{prf:example} Procedure 2 on the network with a negative arc
:label: ex-08-06-proc2-run

Find the shortest paths from node 1 in the network of
{numref}`fig-05-04-negative-arc` with procedure 2, and compare the number of
iterations with procedure 1 ({prf:ref}`ex-05-04-bellman-ford`).
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The arcs are $1\to2$ (2), $1\to3$ (5), $2\to3$ (6), $3\to2$ ($-4$),
$2\to4$ (3), $3\to4$ (1). Within an iteration, nodes are processed in the
order 1, 2, 3, 4; arcs from lower nodes use the new values, arcs from
higher nodes the old ones.

**Iteration 1.** $\tilde V^{(1)}(1)=0$;
$\tilde V^{(1)}(2)=\min\{0+2,\ \tilde V^{(0)}(3)-4\}=2$;
$\tilde V^{(1)}(3)=\min\{0+5,\ \tilde V^{(1)}(2)+6\}=5$;
$\tilde V^{(1)}(4)=\min\{\tilde V^{(1)}(2)+3,\ \tilde V^{(1)}(3)+1\}=\min\{5,6\}=5$.

**Iteration 2.**
$\tilde V^{(2)}(2)=\min\{0+2,\ \tilde V^{(1)}(3)-4\}=\min\{2,1\}=1$;
$\tilde V^{(2)}(3)=\min\{0+5,\ 1+6\}=5$;
$\tilde V^{(2)}(4)=\min\{1+3,\ 5+1\}=4$.

**Iteration 3** repeats iteration 2, so the procedure has converged:

| $k$ | $\tilde V^{(k)}(1)$ | $\tilde V^{(k)}(2)$ | $\tilde V^{(k)}(3)$ | $\tilde V^{(k)}(4)$ |
|---|---|---|---|---|
| 0 | 0 | $\infty$ | $\infty$ | $\infty$ |
| 1 | 0 | 2 | 5 | 5 |
| 2 | 0 | 1 | 5 | 4 |
| 3 | 0 | 1 | 5 | 4 |

The distances $V=(0,1,5,4)$ and the path 1–3–2–4 are those of procedure 1.
The path has one decrease (at node 2), so it is found in iteration 2, and
convergence is detected in iteration 3, one iteration earlier than with
procedure 1 (which needs iteration 3 for the three arcs of the path and
iteration 4 to detect convergence). Already in iteration 1, procedure 2
reaches node 4 with a finite value, because it uses the new value of node 2
at once.
:::

:::{prf:example} Procedure 2 is never worse than procedure 1
:label: ex-08-06-proc2-vs-proc1

Show that $\tilde V^{(k)}(j)\le V^{(k)}(j)$ for all $j$ and $k\ge1$, and
conclude that procedure 2 never needs more operations than procedure 1.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

A decrease occurs at the head of an arc that goes from a higher to a lower
node number. The first arc of a walk leaves node 1, the lowest number, so
it can never be a decrease; hence a walk with at most $k$ arcs has at most
$k-1$ decreases. The walks admitted in $V^{(k)}(j)$ are therefore among
those admitted in $\tilde V^{(k)}(j)$, and the minimum over the larger set
is no larger: $\tilde V^{(k)}(j)\le V^{(k)}(j)$. (For the empty walk at
node 1 both values are 0.)

Now suppose that procedure 1 converges at iteration $k$, so that
$V^{(k-1)}=V$ (no negative cycle; {prf:ref}`thm-05-04-successive`). Every
$\tilde V^{(k-1)}(j)$ is the length of some walk, so it is at least $V(j)$;
together with $\tilde V^{(k-1)}\le V^{(k-1)}=V$ this gives
$\tilde V^{(k-1)}=V$. Iteration $k$ of procedure 2 then cannot change
anything, so procedure 2 has converged by iteration $k$ as well. Both
procedures do $n(n-1)+1$ additions and $n(n-2)+1$ comparisons per
iteration ({eq}`eq-05-04-bf-count`), so procedure 2 never needs more
operations. If a negative cycle exists, both run their $n$ iterations.
:::

:::{prf:example} Correctness of procedure 3
:label: ex-08-06-proc3-proof

Show that $\hat V^{(k)}(j)$, the length of a shortest walk from node 1 to
node $j$ with at most $k-1$ reversals, is given by the recurrences
{eq}`eq-05-04-procedure-3` and the boundary conditions
$\hat V^{(0)}(1)=0$, $\hat V^{(0)}(j)=\infty$ ($j\ne1$).
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Split a walk into **runs**: maximal pieces along which the node numbers
increase, or decrease. Consecutive runs alternate, a reversal separates two
runs, and the first run increases (it starts at node 1, the lowest number).
So a walk with $r$ runs has $r-1$ reversals, and "at most $k-1$ reversals"
means "at most $k$ runs". When the walk has exactly $k$ runs, its last run
increases if $k$ is odd and decreases if $k$ is even. Read $\hat V^{(0)}$ as
"at most 0 runs": only the empty walk at node 1, which is the boundary
condition.

Let $k$ be odd, assume $\hat V^{(k-1)}$ correct, and proceed by induction on
$j=1,2,\dots,n$. Take a walk to $j$ with at most $k$ runs.

- If it has at most $k-1$ runs, its length is at least
  $\hat V^{(k-1)}(j)$, the first term of {eq}`eq-05-04-procedure-3`.
- If it has exactly $k$ runs, its last run increases, so its last arc is
  $(i,j)$ with $i<j$, and the walk up to $i$ has at most $k$ runs. Its
  length is at least $\hat V^{(k)}(i)+c_{ij}$, the second term.

Conversely, every candidate on the right-hand side is the length of a walk
with at most $k$ runs: this is clear for the first term; for the second,
append an increasing arc $(i,j)$ to a walk to $i$ with at most $k$ runs. If
that walk has exactly $k$ runs, its last run increases ($k$ odd) and the new
arc extends it; if it has fewer, the new arc adds at most one run. Hence
{eq}`eq-05-04-procedure-3` holds for odd $k$; the values $\hat V^{(k)}(i)$,
$i<j$, are known because the sweep goes upwards. For even $k$ the argument
is the same with "decreasing", $i>j$ and a downward sweep.
:::

:::{prf:example} Convergence of procedure 3 and negative cycles
:label: ex-08-06-proc3-convergence

Show that if $\hat V^{(k)}=\hat V^{(k-1)}$ for some $k$ with $2\le k\le n$,
then $\hat V^{(k)}(j)=V(j)$ for all $j$. Show also that if convergence does
not occur by iteration $n$, a negative cycle exists, so procedure 3 never
needs more than $n$ iterations.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**After a sweep, its own arcs are satisfied.** After an upward sweep $k$,
every arc $(i,j)$ with $i<j$ satisfies $\hat V^{(k)}(j)\le\hat V^{(k)}(i)+c_{ij}$,
because the minimum in {eq}`eq-05-04-procedure-3` includes this term. After a
downward sweep the same holds for every arc with $i>j$.

**Convergence.** Suppose $\hat V^{(k)}=\hat V^{(k-1)}$ with $k\ge2$. The
values satisfy the inequalities of the arcs of sweep $k-1$ (they are the
values $\hat V^{(k-1)}$) and those of sweep $k$, one family going up and the
other down, so

$$
\begin{aligned}
\hat V^{(k)}(j)&\le\hat V^{(k)}(i)+c_{ij}\quad\text{for every arc }(i,j),\\
\hat V^{(k)}(1)&\le0 .
\end{aligned}
$$

Every later sweep therefore finds nothing to improve, and the values stay
fixed. Take any walk $1=i_0,i_1,\dots,i_m=j$. Adding the inequalities
along it gives $\hat V^{(k)}(j)\le\hat V^{(k)}(1)+\sum c\le$ (length of the
walk). Hence no walk from 1 to $j$ is shorter than $\hat V^{(k)}(j)$, so no
negative cycle can be reached from node 1 (otherwise walks would become
arbitrarily short) and $\hat V^{(k)}(j)\le V(j)$. Since every finite
$\hat V^{(k)}(j)$ is the length of some walk, $\hat V^{(k)}(j)\ge V(j)$ as
well, and the two are equal. Nodes that cannot be reached keep the value
$\infty$, which is also correct.

**At most $n$ iterations.** Without a reachable negative cycle a shortest
walk can be taken to be a path. A path has at most $n$ nodes, hence at most
$n-2$ intermediate nodes and at most $n-2$ reversals, so it is counted in
$\hat V^{(n-1)}$: $\hat V^{(n-1)}=V$. Iteration $n$ then changes nothing and
convergence is detected at the latest at $k=n$. Conversely, if the values
still change in iteration $n$, a negative cycle must exist.
:::

:::{prf:example} Procedure 3 on the network with a negative arc
:label: ex-08-06-proc3-run

Find the shortest paths from node 1 in the network of
{numref}`fig-05-04-negative-arc` with procedure 3 and compare the work with
procedures 1 and 2.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Odd iterations sweep upwards ($j=1,2,3,4$, using arcs from lower nodes),
even iterations downwards ($j=4,3,2,1$, using arcs from higher nodes).

**Iteration 1 (up).** $\hat V^{(1)}(2)=\min\{\infty,\ 0+2\}=2$;
$\hat V^{(1)}(3)=\min\{\infty,\ 0+5,\ 2+6\}=5$;
$\hat V^{(1)}(4)=\min\{\infty,\ 2+3,\ 5+1\}=5$.

**Iteration 2 (down).** $\hat V^{(2)}(4)=5$ (no arc from a higher node);
$\hat V^{(2)}(3)=\min\{5,\ \hat V^{(2)}(4)+c_{43}\}=5$ ($c_{43}=\infty$);
$\hat V^{(2)}(2)=\min\{2,\ \hat V^{(2)}(3)-4\}=1$; $\hat V^{(2)}(1)=0$.

**Iteration 3 (up).** $\hat V^{(3)}(2)=1$, $\hat V^{(3)}(3)=\min\{5,\ 0+5,\ 1+6\}=5$,
$\hat V^{(3)}(4)=\min\{5,\ 1+3,\ 5+1\}=4$.

**Iteration 4 (down)** changes nothing: converged.

| $k$ | sweep | $\hat V^{(k)}(1)$ | $\hat V^{(k)}(2)$ | $\hat V^{(k)}(3)$ | $\hat V^{(k)}(4)$ |
|---|---|---|---|---|---|
| 0 | – | 0 | $\infty$ | $\infty$ | $\infty$ |
| 1 | up | 0 | 2 | 5 | 5 |
| 2 | down | 0 | 1 | 5 | 5 |
| 3 | up | 0 | 1 | 5 | 4 |
| 4 | down | 0 | 1 | 5 | 4 |

The path 1–3–2–4 goes up, down, up: two reversals, so it appears at $k=3$.

**Work** ($n=4$). Procedure 3: 4 iterations of $n(n-1)/2=6$ additions and 6
comparisons, 48 operations. Procedure 2: 3 iterations of $n(n-1)+1=13$
additions and $n(n-2)+1=9$ comparisons, 66 operations. Procedure 1: 4 such
iterations, 88 operations. Procedure 3 needs one more iteration than
procedure 2 but each costs about half as much.
:::

:::{prf:example} Procedure 3 is never worse than procedure 2
:label: ex-08-06-proc3-vs-proc2

Show that $\hat V^{(2k-1)}(j)\le\tilde V^{(k)}(j)$ for all $j$ and
$k\ge1$, and that if procedure 2 converges in $k$ iterations, procedure 3
converges in at most $2k-2$ iterations.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Inequality.** Take a walk with at most $k-1$ decreases. Every decreasing
run contains at least one arc that goes down, and every such arc is a
decrease, so the walk has at most $k-1$ decreasing runs. Runs alternate and
start with an increasing run, so the walk has at most $2(k-1)+1=2k-1$
runs, that is, at most $2k-2$ reversals. Every walk admitted in
$\tilde V^{(k)}(j)$ is therefore admitted in $\hat V^{(2k-1)}(j)$, and
$\hat V^{(2k-1)}(j)\le\tilde V^{(k)}(j)$.

**Iterations.** If procedure 2 converges in iteration $k$ ($k\ge2$), then
$\tilde V^{(k-1)}=V$ and no negative cycle exists. By the inequality,
$\hat V^{(2k-3)}\le\tilde V^{(k-1)}=V$, and since every value is the length
of a walk, $\hat V^{(2k-3)}=V$. The next iteration, $2k-2$, cannot change
these values, so procedure 3 has converged by iteration $2k-2$.

**Consequence.** Procedure 3 needs at most $2k-2$ iterations of $n(n-1)$
operations ({eq}`eq-05-04-procedure-3-count`), procedure 2 needs $k$
iterations of $n(n-1)+1+n(n-2)+1=2n^2-3n+2$ operations
({eq}`eq-05-04-bf-count`). Since $(2k-2)\,n(n-1)<k\,(2n^2-3n+2)$ is equivalent
to $k(n-2)<2n(n-1)$, which holds for every $k\le n$, procedure 3 always needs
fewer operations, and in practice far fewer ({numref}`tab-05-04-procedures`);
the worst cases of the two procedures occur for different networks. On the
four-node example, procedure 2 converged at $k=3$ and procedure 3 at
$k=4=2\cdot3-2$, the bound.
:::

(sec-08-06-all-pairs)=
## All-pairs methods and operation counts

:::{prf:example} Floyd's recursion on a mixed network
:label: ex-08-06-floyd-run

Find the shortest paths between all pairs of nodes of the network in
{numref}`fig-08-06-all-pairs-network` with Floyd's recursion {eq}`eq-05-04-floyd`.
Edges drawn without arrow heads can be used in both directions with the
same length; the pairs of arcs between nodes 2 and 4 and between nodes 3
and 5 are directed: $2\to4$ 3, $4\to2$ $-2$, $3\to5$ $-5$, $5\to3$ 6.
:::

```{figure} figures/08-06-all-pairs-network.svg
:name: fig-08-06-all-pairs-network
:width: 560px
:class: dark-light
:alt: Six-node network with undirected edges and two pairs of opposite arcs, one arc of each pair with a negative length

Network for the all-pairs problems ({prf:ref}`ex-08-06-floyd-run`,
{prf:ref}`ex-08-06-dantzig-run`). Edges without arrow heads have the same
length in both directions; the curved arcs between 2 and 4 and between 3
and 5 are one-way, and two of them are negative.
```

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Boundary matrix** (row = origin, column = destination):

$$
D^{(0)}=\begin{pmatrix}
0&1&2&\infty&\infty&\infty\\
1&0&4&3&9&\infty\\
2&4&0&8&-5&\infty\\
\infty&-2&8&0&7&5\\
\infty&9&6&7&0&3\\
\infty&\infty&\infty&5&3&0
\end{pmatrix}.
$$

**Stages.** At stage $k$ every entry $(i,j)$ with $i,j\ne k$ is compared
with $D^{(k-1)}(i,k)+D^{(k-1)}(k,j)$. The entries that improve:

| Stage $k$ | Improved entries $(i,j)$: old $\to$ new |
|---|---|
| 1 | $(2,3)$: $4\to3$; $(3,2)$: $4\to3$ |
| 2 | $(1,4)$: $\infty\to4$; $(1,5)$: $\infty\to10$; $(3,4)$: $8\to6$; $(4,1)$: $\infty\to-1$; $(4,3)$: $8\to1$; $(5,1)$: $\infty\to10$ |
| 3 | $(1,5)$: $10\to-3$; $(2,5)$: $9\to-2$; $(4,5)$: $7\to-4$; $(5,1)$: $10\to8$ |
| 4 | $(1,6)$: $\infty\to9$; $(2,6)$: $\infty\to8$; $(3,6)$: $\infty\to11$; $(5,1)$: $8\to6$; $(5,2)$: $9\to5$; $(6,1)$: $\infty\to4$; $(6,2)$: $\infty\to3$; $(6,3)$: $\infty\to6$; $(6,5)$: $3\to1$ |
| 5 | $(1,6)$: $9\to0$; $(2,6)$: $8\to1$; $(3,1)$: $2\to1$; $(3,2)$: $3\to0$; $(3,4)$: $6\to2$; $(3,6)$: $11\to-2$; $(4,6)$: $5\to-1$ |
| 6 | none |

The diagonal stays 0 at every stage (for instance at stage 2,
$D^{(1)}(4,2)+D^{(1)}(2,4)=-2+3=1\ge0$, and at stage 3,
$D^{(2)}(5,3)+D^{(2)}(3,5)=6-5=1\ge0$), so there is no negative cycle.

**Result.** Distances $D^{(6)}$ and predecessors $P^{(6)}(i,j)$ of $j$:

$$
D^{(6)}=\begin{pmatrix}
0&1&2&4&-3&0\\
1&0&3&3&-2&1\\
1&0&0&2&-5&-2\\
-1&-2&1&0&-4&-1\\
6&5&6&7&0&3\\
4&3&6&5&1&0
\end{pmatrix},
$$

$$
P^{(6)}=\begin{pmatrix}
\text{–}&1&1&2&3&5\\
2&\text{–}&1&2&3&5\\
2&4&\text{–}&5&3&5\\
2&4&1&\text{–}&3&5\\
2&4&5&5&\text{–}&5\\
2&4&1&6&3&\text{–}
\end{pmatrix}.
$$

**Paths.** Read backwards with the predecessor row of the origin. From 1 to
6: $P(1,6)=5$, $P(1,5)=3$, $P(1,3)=1$, so 1–3–5–6 with $2-5+3=0$. From 3
to 1: $P(3,1)=2$, $P(3,2)=4$, $P(3,4)=5$, $P(3,5)=3$: 3–5–4–2–1 with
$-5+7-2+1=1$. From 6 to 3: 6–4–2–1–3 with $5-2+1+2=6$. All 36 distances
agree with an independent Bellman–Ford computation from every origin.
:::

:::{prf:example} Operation count of Dantzig's recursion
:label: ex-08-06-dantzig-count

Show that Dantzig's recursion {eq}`eq-05-04-dantzig` needs
$\tfrac12n(n-1)(2n-1)$ additions and $\tfrac12n(n-1)(2n-3)$ comparisons,
$2n(n-1)^2$ operations in all. (When $i=j$, the last line of
{eq}`eq-05-04-dantzig` needs no computation, because
$\tilde D^{(k)}(i,i)=\tilde D^{(k-1)}(i,i)=0$ when no negative cycle has
been found.)
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Count stage $k$ ($k=2,\dots,n$) line by line; let $m=k-1$ be the number of
old nodes.

| Line of {eq}`eq-05-04-dantzig` | Entries | Additions | Comparisons |
|---|---|---|---|
| $\tilde D^{(k)}(i,k)$ | $m$ | $m$ each | $m-1$ each |
| $\tilde D^{(k)}(k,j)$ | $m$ | $m$ each | $m-1$ each |
| $\tilde D^{(k)}(k,k)$ | 1 | $m$ | $m$ (with the 0) |
| $\tilde D^{(k)}(i,j)$, $i\ne j$ | $m(m-1)$ | 1 each | 1 each |

Additions at stage $k$: $2m^2+m+m(m-1)=3m^2$. Comparisons:
$2m(m-1)+m+m(m-1)=3m^2-2m$. Summing over $m=1,\dots,n-1$ with
$\sum m=\tfrac12n(n-1)$ and $\sum m^2=\tfrac16(n-1)n(2n-1)$:

$$
\begin{aligned}
\text{additions}&=3\sum_{m=1}^{n-1}m^2=\tfrac12n(n-1)(2n-1),\\
\text{comparisons}&=3\sum_{m=1}^{n-1}m^2-2\sum_{m=1}^{n-1}m
=\tfrac12n(n-1)(2n-1)-n(n-1)=\tfrac12n(n-1)(2n-3).
\end{aligned}
$$

The total is $\tfrac12n(n-1)\bigl[(2n-1)+(2n-3)\bigr]=2n(n-1)^2$, exactly
the count {eq}`eq-05-04-floyd-count` of Floyd's recursion. For $n=4$, 6, 10
and 20 this gives 42 + 30, 165 + 135, 855 + 765 and 7 410 + 7 030
operations, which agree with an instrumented run of the recursion.
:::

:::{prf:example} Dantzig's recursion on the mixed network
:label: ex-08-06-dantzig-run

Find all shortest paths of the network in {numref}`fig-08-06-all-pairs-network`
with Dantzig's recursion {eq}`eq-05-04-dantzig`.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Stage $k$ adds node $k$; the matrices hold $\tilde D^{(k)}(i,j)$ for
$i,j\le k$ (rows = origins).

**Stage 2.** $\tilde D(1,2)=c_{12}=1$, $\tilde D(2,1)=c_{21}=1$:
$\begin{pmatrix}0&1\\1&0\end{pmatrix}$.

**Stage 3.** $\tilde D(1,3)=\min\{0+2,\ 1+4\}=2$,
$\tilde D(2,3)=\min\{1+2,\ 0+4\}=3$, $\tilde D(3,1)=\min\{2+0,\ 4+1\}=2$,
$\tilde D(3,2)=\min\{2+1,\ 4+0\}=3$, $\tilde D(3,3)=\min\{0,\ 2+2,\ 3+3\}=0$;
no old entry improves:
$\begin{pmatrix}0&1&2\\1&0&3\\2&3&0\end{pmatrix}$.

**Stage 4.** $\tilde D(1,4)=\min\{\infty,\ 1+3,\ 2+8\}=4$,
$\tilde D(2,4)=\min\{\infty,\ 0+3,\ 3+8\}=3$,
$\tilde D(3,4)=\min\{\infty,\ 3+3,\ 0+8\}=6$,
$\tilde D(4,1)=\min\{\infty,\ -2+1,\ 8+2\}=-1$,
$\tilde D(4,2)=\min\{\infty,\ -2+0,\ 8+3\}=-2$,
$\tilde D(4,3)=\min\{\infty,\ -2+3,\ 8+0\}=1$;
$\tilde D(4,4)=\min\{0,\ -1+4,\ -2+3,\ 1+6\}=0$; no old entry improves
(for example $\tilde D(2,1)$: $3-1=2>1$).

**Stage 5.** $\tilde D(i,5)$ for $i=1,\dots,4$: $-3$, $-2$, $-5$, $-4$
(all via the arc $3\to5$; for example, $\tilde D(1,5)=\tilde D(1,3)-5=-3$);
$\tilde D(5,j)$ for $j=1,\dots,4$: $6$, $5$, $6$, $7$ (for example,
$\tilde D(5,1)=\min\{9+1,\ 6+2,\ 7-1\}=6$ via node 4);
$\tilde D(5,5)=\min\{0,\ 6-3,\ 5-2,\ 6-5,\ 7-4\}=0$. Old entries improved
through node 5: $\tilde D(3,1)=-5+6=1$, $\tilde D(3,2)=-5+5=0$,
$\tilde D(3,4)=-5+7=2$.

**Stage 6.** $\tilde D(i,6)$ for $i=1,\dots,5$: $0$, $1$, $-2$, $-1$, $3$
(via the edge 5–6 for $i\le4$, for example $\tilde D(1,6)=-3+3=0$);
$\tilde D(6,j)$ for $j=1,\dots,5$: $4$, $3$, $6$, $5$, $1$ (via node 4 or
node 5, for example $\tilde D(6,1)=c_{64}+\tilde D(4,1)=5-1=4$);
$\tilde D(6,6)=0$; no old entry improves.

The final matrix is the matrix $D^{(6)}$ of {prf:ref}`ex-08-06-floyd-run`. The
two recursions differ only in the order in which the same comparisons are
made: Floyd's lets every pair use node $k$ at stage $k$, Dantzig's builds
the full distance matrix of the first $k$ nodes.
:::

(sec-08-06-label-setting)=
## Label-setting methods and special variants

:::{prf:example} Dijkstra's method as a dynamic program
:label: ex-08-06-dijkstra-dp

Find the shortest paths from node 1 in the network whose arcs are $1\to2$ 5, $1\to3$ 1,
$2\to3$ 4, $2\to4$ 6, $3\to4$ 7, $3\to5$ 2, $5\to2$ 1, $5\to4$ 3, $4\to6$ 2
and $5\to6$ 8, with the DP form {eq}`eq-05-04-dijkstra` of Dijkstra's method.
Record the set $\mathcal S_k$, the node $i_k$ and the labels $T^{(k)}$ at
every stage.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

```{figure} figures/08-06-dijkstra-network.svg
:name: fig-08-06-dijkstra-network
:width: 560px
:class: dark-light
:alt: Six-node directed network with a cycle, nonnegative arc lengths, labels V and P, and the shortest-path tree from node 1 highlighted

Network of {prf:ref}`ex-08-06-dijkstra-dp`. It contains the cycle 2–3–5–2, so
there is no topological order. Labels $[V(j),P(j)]$ and the highlighted
shortest-path tree are the result of the example.
```

Stage 1: $\mathcal S_1=\{1\}$, $T^{(1)}(j)=c_{1j}$. At each later stage the
unlabeled node with the smallest $T$ is added, and the labels of the
remaining nodes are updated through it ($\infty$ = no path yet; bold =
permanent):

| Stage $k$ | $i_k$ | $T(1)$ | $T(2)$ | $T(3)$ | $T(4)$ | $T(5)$ | $T(6)$ |
|---|---|---|---|---|---|---|---|
| 1 | 1 | **0** | 5 | 1 | $\infty$ | $\infty$ | $\infty$ |
| 2 | 3 | **0** | 5 | **1** | $\min\{\infty,1+7\}=8$ | $\min\{\infty,1+2\}=3$ | $\infty$ |
| 3 | 5 | **0** | $\min\{5,3+1\}=4$ | **1** | $\min\{8,3+3\}=6$ | **3** | $\min\{\infty,3+8\}=11$ |
| 4 | 2 | **0** | **4** | **1** | $\min\{6,4+6\}=6$ | **3** | 11 |
| 5 | 4 | **0** | **4** | **1** | **6** | **3** | $\min\{11,6+2\}=8$ |
| 6 | 6 | **0** | **4** | **1** | **6** | **3** | **8** |

The answers are $V(i_k)=T^{(k-1)}(i_k)$: $V=(0,4,1,6,3,8)$, found in the
order 1, 3, 5, 2, 4, 6, which is the order of increasing distance. The
predecessors (the node through which each final label was last improved)
are $P(3)=1$, $P(5)=3$, $P(2)=5$, $P(4)=5$, $P(6)=4$: the tree highlighted
in {numref}`fig-08-06-dijkstra-network`. For example, the shortest path to 6 is
1–3–5–4–6 of length $1+2+3+2=8$.

**Work.** Stage $k$ scans $n-k+1$ unlabeled nodes ($n-k$ comparisons) and
updates $n-k$ labels: with $n=6$, $(n-1)(n-2)/2=10$ additions and
$(n-1)(n-2)=20$ comparisons in the count {eq}`eq-05-04-dijkstra-count`, which
treats missing arcs as arcs of length $\infty$.
:::

:::{prf:example} Minimax paths
:label: ex-08-06-minimax-path

In a network with $n$ nodes and arc numbers $c_{ij}$ of any sign, the
value of a path is the largest $c_{ij}$ on it. Give an efficient procedure
for paths of least value from node 1 to all other nodes.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Let $V(j)$ be the least value of a path from 1 to $j$. A path is extended by
an arc $(i,j)$ with the combination $\max\{V(i),c_{ij}\}$, which never
decreases when the path is extended, *whatever the signs of the numbers*.
This is exactly the property that Dijkstra's method needs (in its
correctness proof, the nonnegativity of the lengths was used only to show
that extending a path cannot make it better). So Dijkstra's method applies
with addition replaced by maximum:

1. $\mathcal S=\{1\}$, $V(1)=-\infty$ (the empty path has no arc), and
   $T(j)=c_{1j}$ for $j\ne1$.
2. Make the unlabeled node $k$ with the smallest $T(k)$ permanent:
   $V(k)=T(k)$.
3. For every unlabeled $j$: $T(j)\leftarrow\min\{T(j),\ \max\{V(k),c_{kj}\}\}$,
   and record $k$ as predecessor when the label improves. Repeat from step 2
   until all nodes are labeled.

Cycles do no harm: going round a cycle can never lower the maximum, so there
is no analog of a negative cycle. The work is that of Dijkstra's method,
about $\tfrac32n^2$ operations (maxima instead of additions). This is the
directed version of the minimax problem on the lattice
({prf:ref}`ex-08-06-lattice-d`); in an undirected network the minimax paths
are also the paths in a minimum spanning tree ({doc}`04-01-networks-mst`).

**Illustration.** In {numref}`fig-05-04-negative-arc`: $T=(\cdot,2,5,\infty)$;
node 2 becomes permanent with $V(2)=2$; then
$T(3)=\min\{5,\max\{2,6\}\}=5$ and $T(4)=\max\{2,3\}=3$; node 4 is next
with $V(4)=3$ (path 1–2–4), and node 3 last with $V(3)=5$ (the arc
$1\to3$).
:::

:::{prf:example} Most reliable route
:label: ex-08-06-reliable-path

A message is sent from node 1 to node $n$. If it is routed from node $i$ to
node $j$, it arrives at $j$ with known probability $p_{ij}$ and is lost
otherwise, independently of what happens on other arcs. Give an efficient
procedure for the route with the largest probability of arrival.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

By independence, the arrival probability of a route is the product of its
$p_{ij}$, and we want to maximize it. Two equivalent procedures:

* **Transform.** With $c_{ij}=-\ln p_{ij}\ge0$ (and $c_{ij}=\infty$ when
  $p_{ij}=0$), maximizing $\prod p_{ij}$ is minimizing
  $\sum c_{ij}=-\ln\prod p_{ij}$, an ordinary shortest-path problem with
  nonnegative lengths. Dijkstra's method solves it in about
  $\tfrac32n^2$ operations, and the arrival probability is
  $e^{-V(n)}$.
* **Direct.** Run Dijkstra's method with "largest" in place of "smallest"
  and multiplication in place of addition: $R(1)=1$; repeatedly make
  permanent the unlabeled node with the *largest* tentative probability and
  update $R(j)\leftarrow\max\{R(j),\ R(k)\,p_{kj}\}$. It is valid because
  a product of probabilities can only decrease as a route is extended, the
  mirror image of nonnegative lengths.

As in the product problem on the lattice ({prf:ref}`ex-08-06-lattice-g`), the
boundary value is 1 (the empty product), and only positive factors keep the
principle of optimality intact. One may stop as soon as node $n$ becomes
permanent.
:::

:::{prf:example} Time-dependent travel times
:label: ex-08-06-time-dependent

In a network with $n$ nodes, $c_{ij}(t)$ is the time needed to go from node
$i$ to node $j$ when node $i$ is left at time $t$. Give an efficient
procedure for the earliest possible arrival time at every node when node 1
is left at time 0.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Let $V(j)$ be the earliest arrival time at node $j$, and
$a_{ij}(t)=t+c_{ij}(t)$ the arrival time at $j$ when $i$ is left at $t$.

**Assumption (no overtaking).** $a_{ij}(t)$ is nondecreasing in $t$:
leaving later never means arriving earlier. This holds for most traffic
models ("first in, first out"). Under it, arriving at a node earlier is
never worse, so the earliest arrival time is a correct state, and the
principle of optimality holds: an earliest route to $j$ through $i$ reaches
$i$ as early as possible.

**Procedure** (Dijkstra's method with arrival times as labels):

1. $V(1)=0$; $T(j)=a_{1j}(0)$ for $j\ne1$.
2. Make permanent the unlabeled node $k$ with the smallest $T(k)$:
   $V(k)=T(k)$.
3. For every unlabeled $j$: $T(j)\leftarrow\min\{T(j),\ a_{kj}(V(k))\}$.
   Repeat from step 2.

Correctness follows as for Dijkstra's method, because $a_{ij}(t)\ge t$
(travel times are nonnegative) plays the role of nonnegative lengths, and
monotonicity of $a_{ij}$ ensures that leaving $k$ at its earliest time is
best. The work is about $\tfrac32n^2$ operations plus $\tfrac12n^2$
evaluations of the functions $c_{ij}$.

**Without the assumption.** If waiting at nodes is allowed, replace $a_{ij}$
by $\bar a_{ij}(t)=\min_{s\ge t}a_{ij}(s)$ (wait for the best departure),
which is nondecreasing, and apply the same procedure. If waiting is not
allowed and overtaking is possible, arriving earlier can be worse, so the
state must be the pair (node, time), and the problem becomes a shortest
path in the time-expanded network, which is acyclic in time and can be
solved by the one-pass recursion of {doc}`05-04-dp-shortest-paths`.
:::

:::{prf:example} Dijkstra's method with a turn penalty
:label: ex-08-06-turn-dijkstra

Nine nodes form a diamond grid: node 1 in the first column, nodes 2
(upper) and 3 (lower) in the second, nodes 4, 5, 6 (top to bottom) in the
third, nodes 7 (upper) and 8 (lower) in the fourth and node 9 in the fifth.
Every arc goes one column to the right, up or down: $1\to2$ 3, $1\to3$ 1, $2\to4$ 4, $2\to5$ 2,
$3\to5$ 5, $3\to6$ 4, $4\to7$ 4, $5\to7$ 2, $5\to8$ 3, $6\to8$ 5,
$7\to9$ 3 and $8\to9$ 4. Each turn (a change between an up arc and a down
arc) costs an extra 2. Modify Dijkstra's method to find the cheapest path
from node 1 to node 9. (Hint: give each node one label for every arc that
enters it.)
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

```{figure} figures/08-06-turn-grid.svg
:name: fig-08-06-turn-grid
:width: 400px
:class: dark-light
:alt: Nine-node diamond grid with arc lengths and the path 1-3-5-7-9 highlighted

Network of {prf:ref}`ex-08-06-turn-dijkstra`. With a penalty of 2 per turn
the cheapest path is 1–3–5–7–9 (length 11, two turns, total 15); without
the penalty it would be 1–2–5–7–9 (length 10).
```

**Labels.** The cost of leaving a node depends on how it was entered, so the
node alone is not a state. Give node $j$ one label $T(j,a)$ for each arc $a$
that enters it: the least cost of a path from 1 that enters $j$ by arc $a$.
The labels are the nodes of an expanded network in which "enter $j$ by $a$,
then leave by $a'$" costs $c_{a'}$ plus 2 if $a$ and $a'$ point in different
directions. Dijkstra's method on the expanded network is valid (all costs
are nonnegative). Making one label of a node permanent does not settle the
others: a path that enters the node by another arc may be more expensive
now but cheaper later, because it needs no turn.

**Run.** Labels are made permanent in this order (each with the path that
gives it):

| Order | Label | Value | Path |
|---|---|---|---|
| 1 | node 1 | 0 | 1 |
| 2 | 3 by $1\to3$ | 1 | 1–3 |
| 3 | 2 by $1\to2$ | 3 | 1–2 |
| 4 | 6 by $3\to6$ | $1+4=5$ | 1–3–6 |
| 5 | 4 by $2\to4$ | $3+4=7$ | 1–2–4 |
| 6 | 5 by $2\to5$ | $3+2+2=7$ | 1–2–5 (turn at 2) |
| 7 | 5 by $3\to5$ | $1+5+2=8$ | 1–3–5 (turn at 3) |
| 8 | 7 by $5\to7$ | $8+2=10$ | 1–3–5–7 |
| 9 | 8 by $5\to8$ | $7+3=10$ | 1–2–5–8 |
| 10 | 8 by $6\to8$ | $5+5+2=12$ | 1–3–6–8 (turn at 6) |
| 11 | 7 by $4\to7$ | $7+4+2=13$ | 1–2–4–7 (turn at 4) |
| 12 | 9 by $7\to9$ | $10+3+2=15$ | 1–3–5–7–9 (turn at 7) |

When the first label of node 9 becomes permanent the procedure stops: the
cheapest path is 1–3–5–7–9 with arc lengths $1+5+2+3=11$ and two turns
(at 3 and 7), total 15 ({numref}`fig-08-06-turn-grid`). (The other label of node 9, by $8\to9$, would be
16.) Checking all six paths: 1–2–4–7–9 costs $14+2=16$, 1–2–5–7–9
$10+6=16$, 1–2–5–8–9 $12+4=16$, 1–3–5–7–9 $11+4=15$, 1–3–5–8–9 $13+6=19$,
1–3–6–8–9 $14+2=16$. Note that 1–3–5–7–9 enters node 5 by its *more*
expensive label (8 rather than 7), which is why both labels must be kept.

If the lines other than $1\to2$ and $1\to3$ are taken as two-way roads
(a path may then also move to the left), the answer does not change: any
move to the left adds length and turns. The work is that of Dijkstra's
method on the expanded network, which has one label per arc plus the
origin: 13 labels instead of 9 nodes.
:::

:::{prf:example} Shortest path between two sets of nodes
:label: ex-08-06-set-to-set

In a network with $n$ nodes and nonnegative lengths, $\mathcal N_A$ and
$\mathcal N_B$ are disjoint sets of nodes with $n_A$ and $n_B$ elements. Give an
efficient Dijkstra-like procedure for the shortest path that starts at some
node of $\mathcal N_A$ and ends at some node of $\mathcal N_B$ (both chosen
optimally). How many additions and comparisons does it need?
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Idea.** Add a super-source $0$ joined to every node of $\mathcal N_A$ by an
arc of length 0. A shortest path from 0 to a node of $\mathcal N_B$, with its
first arc removed, is a shortest $\mathcal N_A$–$\mathcal N_B$ path. In
Dijkstra's method, the super-source and then all of $\mathcal N_A$ become
permanent with label 0 at once, so they need not be created explicitly.

**Procedure.**

1. Make every node of $\mathcal N_A$ permanent with $V=0$. For every other
   node $j$ set $T(j)=\min_{i\in\mathcal N_A}c_{ij}$ and remember the
   minimizing $i$.
2. Make permanent the unlabeled node $k$ with the smallest $T(k)$. If
   $k\in\mathcal N_B$, stop: $T(k)$ is the answer, and the path is traced back
   through the predecessors to its starting node in $\mathcal N_A$.
3. Update $T(j)\leftarrow\min\{T(j),\ V(k)+c_{kj}\}$ for every unlabeled
   $j$ and return to step 2.

The first node of $\mathcal N_B$ to become permanent is the closest one to the
set $\mathcal N_A$, by the correctness of Dijkstra's method.

**Operation count.** Step 1 needs no additions and $n_A-1$ comparisons for
each of the $n-n_A$ other nodes. Steps 2–3 are Dijkstra's method on a network
whose origin is the merged set $\mathcal N_A$, that is, with $n-n_A+1$ nodes; by
{eq}`eq-05-04-dijkstra-count` they need at most

$$
\begin{aligned}
&\frac{(n-n_A)(n-n_A-1)}{2}\ \text{additions},\\
&(n-n_A)(n-n_A-1)\ \text{comparisons},
\end{aligned}
$$

plus the $(n_A-1)(n-n_A)$ comparisons of step 1. Stopping at the first node of
$\mathcal N_B$ saves at least the last $n_B-1$ iterations, so these are upper
bounds; for $n_A=n_B=1$ they reduce to the count of the single-pair problem.
:::

## References

- Bellman, R. (1958). On a routing problem. *Quarterly of Applied
  Mathematics, 16*(1), 87–90. <https://doi.org/10.1090/qam/102435>
- Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
  *Numerische Mathematik, 1*(1), 269–271.
  <https://doi.org/10.1007/BF01386390>
- Dreyfus, S. E., & Law, A. M. (1977). *The art and theory of dynamic
  programming*. Academic Press. (Ch. 1: path problems on lattices with turn
  penalties, coupons, minimax and product costs; Ch. 4: the general
  shortest-path problem, its procedures and Problems 4.1–4.19, restated
  here.)
- Floyd, R. W. (1962). Algorithm 97: Shortest path. *Communications of the
  ACM, 5*(6), 345. <https://doi.org/10.1145/367766.368168>
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Ch. 12: §12.1–12.3, the hiking-schedule and staged-route
  problems.)
