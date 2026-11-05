# Sensitivity Analysis

## Overview

The data of a linear program (LP) are rarely exact: costs and prices are
forecasts, resource limits are budgets that may be renegotiated, and new
products or new restrictions appear after the model has been solved.
**Sensitivity analysis** (post-optimal analysis) answers, from the optimal
basis alone and without solving again, how far the data may move before
that basis stops being optimal, how the optimal value changes inside that
range, and how to re-optimize cheaply outside it. Everything follows from
one observation: a basis is optimal exactly when it is primal feasible and
dual feasible, and each kind of data change touches only one of the two
conditions. Changes of costs and new variables affect dual feasibility and
are repaired by primal simplex pivots; changes of right-hand sides and new
constraints affect primal feasibility and are repaired by dual simplex pivots
({doc}`03-06-revised-dual-simplex`). The chapter derives the ranging rules for
nonbasic and basic costs and for right-hand sides, explains shadow prices
and their range of validity, the 100 percent rule for simultaneous changes,
new activities and new constraints, and the parametric right-hand side, all
on one running example with three constraints and six variables, organized
around the bordered inverse $\mT$ {eq}`eq-03-06-bordered-inverse`. Dual
values and their economic meaning come from {doc}`03-05-lp-duality`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- state which data changes affect primal feasibility and which affect dual
  feasibility of an optimal basis;
- compute the range of a nonbasic cost, a basic cost and a right-hand side
  over which the optimal basis stays optimal, and the optimal value inside
  the range;
- interpret shadow prices, their validity range and the 100 percent rule;
- price a new activity and add a new constraint using the $\mT$-matrix;
- re-optimize with primal or dual simplex pivots when a range is left; and
- trace the optimal value as a piecewise-linear function of one right-hand
  side.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\mB$, $\vx_B^*$, $\vpi$ | optimal basis matrix, optimal basic solution and its multipliers, $\vpi^\T=\vc_B^\T\mB^{-1}$ |
| $\delta$ | change of a single coefficient: $c_r\to c_r+\delta$ or $b_r\to b_r+\delta$ |
| $\bar a_{ij}$ | entries of the optimal tableau $\mB^{-1}\mA$ |
| $\mB^{-1}\ve_r$ | column $r$ of $\mB^{-1}$ |
| $z^*(\cdot)$ | optimal value as a function of the changing data |
| $t$ | scalar parameter of a parametric right-hand side $\vb+t\vd$ |

(sec-03-07-optimality-conditions)=
## What keeps a basis optimal

Throughout, the LP is $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ with
$\mA$ of full row rank, solved to an optimal basis $\mB$. By
{prf:ref}`thm-03-03-optimality` and {prf:ref}`prop-03-05-basis-optimal`, $\mB$ is
optimal when two conditions hold; for a nondegenerate basis they are also
necessary, and sensitivity analysis keeps track of exactly these two
conditions.

:::{prf:theorem} Optimality conditions of a basis
:label: thm-03-07-basis-optimal

A basis $\mB$ is optimal for $\min\{\vc^\T\vx:\mA\vx=\vb,\vx\ge\vzero\}$ if

$$
\mB^{-1}\vb\ge\vzero\qquad\text{(primal feasibility)}
$$ (eq-03-07-primal-feasible)

and

$$
\bar c_j=c_j-\vc_B^\T\mB^{-1}\mA_{\cdot j}\ge0\ \text{ for all } j
\qquad\text{(dual feasibility)}.
$$ (eq-03-07-dual-feasible)

For a maximization, {eq}`eq-03-07-dual-feasible` reads $\bar c_j\le0$.
:::

Condition {eq}`eq-03-07-primal-feasible` involves $\mB$ and $\vb$ only;
condition {eq}`eq-03-07-dual-feasible` involves $\mB$, $\vc$ and the columns
$\mA_{\cdot j}$ only. Hence:

* a change of $\vc$ leaves {eq}`eq-03-07-primal-feasible` intact and can only
  break {eq}`eq-03-07-dual-feasible`; the same holds for a new variable or a
  changed nonbasic column;
* a change of $\vb$ leaves {eq}`eq-03-07-dual-feasible` intact and can only
  break {eq}`eq-03-07-primal-feasible`; the same holds for a new constraint,
  whose slack enters the basis.

Within the range where both conditions hold, the optimal **basis** stays
the same, although the optimal **values** may move: $\vx_B^*=\mB^{-1}\vb$
moves with $\vb$, and $z^*=\vc_B^\T\vx_B^*$ moves with $\vc_B$ and $\vb$.
Outside the range, the old basis is still a good starting point: a basis
that has lost dual feasibility is primal feasible, so the primal simplex
method continues from it; a basis that has lost primal feasibility is dual
feasible, so the dual simplex method continues from it. Usually a few
pivots suffice. This chapter studies one change at a time; simultaneous
changes are treated in {ref}`sec-03-07-hundred-percent` and
{ref}`sec-03-07-parametric`.

(sec-03-07-geometry)=
## A geometric view

Two pictures make the algebra concrete. They use the two-product LP of
{prf:ref}`ex-08-04-a6-q1-hidec`,

$$
\begin{aligned}
\max_{\vx}\ & 3x_1+4x_2 \\
\st\ & 2x_1+3x_2\le1200\quad\text{(resistors)},\\
& 2x_1+x_2\le1000\quad\text{(capacitors)},\\
& 4x_2\le800\quad\text{(chips)},\\
& x_1,x_2\ge0,
\end{aligned}
$$

whose optimum is $\vx^*=(450,100)$ with revenue \$1750, at the corner where
the resistor and capacitor constraints bind.

**Changing a right-hand side moves a constraint parallel to itself.**
In {numref}`fig-03-07-rhs-ranging` the resistor limit $b_1$ varies. The optimal
corner stays the intersection of the same two lines, hence the same basis,
while the resistor line moves between $b_1=1000$, where the corner reaches
the $x_1$-axis, and $b_1=1400$, where it reaches the chip line. Inside this
range the optimum slides along the capacitor line and the revenue grows
linearly, by the dual price \$1.25 per resistor: from \$1500 at $b_1=1000$
to \$2000 at $b_1=1400$. Beyond either end another constraint becomes
binding, the basis changes, and so does the slope.

```{figure} figures/03-07-rhs-ranging.svg
:name: fig-03-07-rhs-ranging
:width: 520px
:class: dark-light
:alt: Feasible region of a two-variable LP with the resistor constraint drawn at three right-hand sides; the optimal corner slides along the capacitor line

Right-hand-side ranging. The resistor line $2x_1+3x_2=b_1$ is drawn for
$b_1=1000$, $1200$ (the current value, shaded region) and $1400$. Over this
range the optimal corner (dots) moves along the capacitor line
$2x_1+x_2=1000$ from $(500,0)$ to $(400,200)$; the same two constraints
bind, so the basis is unchanged.
```

**Changing a cost rotates the objective contour.** In
{numref}`fig-03-07-cost-ranging` the price $c_1$ of the first product varies
with $c_2=4$ fixed. The contour $c_1x_1+4x_2=z$ pivots about the optimal
corner. The corner stays optimal while the contour's slope $-c_1/4$ lies
between the slopes of the two binding constraints, $-2/3$ (resistors) and
$-2$ (capacitors), that is for $8/3\le c_1\le8$. At either end the contour
lies along a binding constraint, a whole edge is optimal, and beyond it the
optimum jumps to the neighboring corner, $(300,200)$ or $(500,0)$.

```{figure} figures/03-07-cost-ranging.svg
:name: fig-03-07-cost-ranging
:width: 520px
:class: dark-light
:alt: Feasible region with three objective contours through the optimal corner, for three values of the first price

Cost ranging. Objective contours $c_1x_1+4x_2=z$ through the optimal corner
$(450,100)$ for $c_1=8/3$ (the contour coincides with the resistor line),
$c_1=3$ (current price) and $c_1=8$ (it coincides with the capacitor line);
the constraint lines are gray. For every $c_1$ strictly between $8/3$ and 8
the contour touches the region only at this corner; at the two ends it
touches a whole edge.
```

The two pictures are the geometric content of the two optimality conditions:
the right-hand side decides *where* the corner is (primal feasibility), the
costs decide *which* corner is best (dual feasibility). With more than two
variables pictures fail, and the algebra below takes over.

(sec-03-07-running)=
## The running example and the $\mT$-matrix

All rules below are applied to one LP.

:::{prf:example} Running example
:label: ex-03-07-running-lp

$$
\begin{aligned}
\min_{\vx}\ & 3x_1+2x_2-3x_3-6x_4+10x_5-5x_6 \\
\st\ & x_1+2x_2+x_4-6x_6=11,\\
& x_2+x_3+3x_4-2x_5-x_6=6,\\
& x_1+2x_2+x_3+3x_4-x_5-5x_6=13,\\
& x_1,\dots,x_6\ge0 .
\end{aligned}
$$ (eq-03-07-running-lp)

Show that $(x_1,x_2,x_3)$ is an optimal basis: give $\mB^{-1}$, the basic
solution, the multipliers $\vpi$, the optimal value and the optimal tableau.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The optimal basis is $(x_1,x_2,x_3)$ with

$$
\begin{gathered}
\mB=\begin{pmatrix}1&2&0\\0&1&1\\1&2&1\end{pmatrix},\\
\mB^{-1}=\begin{pmatrix}-1&-2&2\\1&1&-1\\-1&0&1\end{pmatrix},
\end{gathered}
$$

so that $\vx_B^*=\mB^{-1}\vb=(3,4,2)^\T$, the multipliers are
$\vpi^\T=(3,2,-3)\mB^{-1}=(2,-4,1)$ and $z^*=\vpi^\T\vb=22-24+13=11$. The
optimal tableau is {numref}`tab-03-07-optimal-tableau`; every reduced cost is
nonnegative and every basic value is positive, so the basis is optimal and
nondegenerate. A solver (HiGHS) confirms $\vx^*=(3,4,2,0,0,0)$, $z^*=11$
and the duals $(2,-4,1)$.

```{table} Optimal tableau of the running example.
:name: tab-03-07-optimal-tableau

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 0 | 0 | −1 | 2 | −2 | 3 |
| $x_2$ | 0 | 1 | 0 | 1 | −1 | −2 | 4 |
| $x_3$ | 0 | 0 | 1 | 2 | −1 | 1 | 2 |
| $-z$ | 0 | 0 | 0 | 1 | 3 | 8 | −11 |
```
:::

For sensitivity analysis it is convenient to keep the optimal basis
$(x_1,x_2,x_3)$ of {prf:ref}`ex-03-07-running-lp`, with
$\mB^{-1}$ and $\vpi^\T=(2,-4,1)$ from its solution, in bordered form, as the matrix $\mT$ of {eq}`eq-03-06-bordered-inverse`:

$$
\begin{aligned}
\mT=\begin{pmatrix}\mB&\vzero\\\vc_B^\T&1\end{pmatrix}^{-1}
&=\begin{pmatrix}\mB^{-1}&\vzero\\-\vpi^\T&1\end{pmatrix}\\
&=\begin{pmatrix}-1&-2&2&0\\1&1&-1&0\\-1&0&1&0\\-2&4&-1&1\end{pmatrix}.
\end{aligned}
$$ (eq-03-07-T-running)

By {eq}`eq-03-06-T-column`, $\mT$ turns any column of data into a column of the
optimal tableau with its reduced cost:
$\mT(\mA_{\cdot j};c_j)=(\bar\mA_{\cdot j};\bar c_j)$ and
$\mT(\vb;0)=(\vx_B^*;-z^*)=(3,4,2,-11)^\T$. A new or modified column is
therefore priced and updated with one matrix–vector product, without
touching the rest of the tableau. If the LP was solved by the bordered
revised method of {ref}`sec-03-06-two-phase-revised`, $\mT$ is available at the
end at no cost; in a tableau that started from slack variables it is the
block under the slack columns.

(sec-03-07-nonbasic-cost)=
## Ranging a nonbasic cost

Let $x_r$ be nonbasic and let only $c_r$ change. The multipliers
$\vpi^\T=\vc_B^\T\mB^{-1}$ do not involve $c_r$, so the only reduced cost
that changes is $\bar c_r=c_r-\vpi^\T\mA_{\cdot r}$ itself, and it changes one
for one with $c_r$. The basis stays optimal if and only if

$$
c_r\ \ge\ \vpi^\T\mA_{\cdot r}\qquad\text{(minimization)},
$$ (eq-03-07-nonbasic-range)

with $\le$ for a maximization. The quantity $\vpi^\T\mA_{\cdot r}$ is the
**break-even cost** of activity $r$: the value, at the current shadow prices,
of the resources one unit of $x_r$ would use. The optimal solution and value
do not change at all inside the range, because $x_r^*=0$. Below the range
$\bar c_r<0$, and $x_r$ enters the basis; the primal simplex method
continues from the old basis.

:::{prf:example} Cost of a nonbasic variable
:label: ex-03-07-nonbasic-cost

In the running example, find the range of $c_5$ (now 10) over which the
basis $(x_1,x_2,x_3)$ stays optimal, and re-optimize for $c_5=6$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The break-even cost is
$\vpi^\T\mA_{\cdot5}=(2,-4,1)(0,-2,-1)^\T=8-1=7$, so
$\bar c_5=c_5-7\ge0$ and the range is $c_5\ge7$. As long as producing $x_5$
costs at least 7 per unit, it stays at zero.

For $c_5=6$, the column of $x_5$ in the optimal tableau, with its new reduced
cost, is

$$
\mT\begin{pmatrix}\mA_{\cdot5}\\c_5\end{pmatrix}
=\begin{pmatrix}-1&-2&2&0\\1&1&-1&0\\-1&0&1&0\\-2&4&-1&1\end{pmatrix}
\begin{pmatrix}0\\-2\\-1\\6\end{pmatrix}
=\begin{pmatrix}2\\-1\\-1\\-1\end{pmatrix}.
$$

So $\bar c_5=-1<0$ and $x_5$ enters. Only row $x_1$ has a positive entry,
the ratio is $3/2$, and $x_1$ leaves. After the pivot:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_5$ | 1/2 | 0 | 0 | −1/2 | 1 | −1 | 3/2 |
| $x_2$ | 1/2 | 1 | 0 | 1/2 | 0 | −3 | 11/2 |
| $x_3$ | 1/2 | 0 | 1 | 3/2 | 0 | 0 | 7/2 |
| $-z$ | 1/2 | 0 | 0 | 1/2 | 0 | 7 | −19/2 |

All reduced costs are nonnegative: the new optimum is
$\vx^*=(0,\tfrac{11}{2},\tfrac72,0,\tfrac32,0)$ with
$z^*=2\cdot\tfrac{11}{2}-3\cdot\tfrac72+6\cdot\tfrac32=\tfrac{19}{2}$, lower
than 11, as it must be for a cheaper activity. One pivot replaced a new
solve; HiGHS gives the same optimum. At exactly $c_5=7$ the reduced cost is
0, and both bases are optimal (alternative optima).
:::

(sec-03-07-basic-cost)=
## Ranging a basic cost

Now let $x_r$ be basic, in row $i$ of the tableau, and change $c_r$ to
$c_r+\delta$. Then $\vc_B$ changes by $\delta\ve_i$, the multipliers change
by $\delta\,\ve_i^\T\mB^{-1}$, and every nonbasic reduced cost changes:

$$
\begin{aligned}
\bar c_j(\delta)&=c_j-(\vpi+\delta\,\mB^{-\T}\ve_i)^\T\mA_{\cdot j}\\
&=\bar c_j-\delta\,\bar a_{ij}\qquad(j\text{ nonbasic}),
\end{aligned}
$$ (eq-03-07-basic-cost)

where $\bar a_{ij}$ is the entry of row $i$ (the row of $x_r$) of the
optimal tableau. Requiring $\bar c_j(\delta)\ge0$ for all nonbasic $j$
gives the range

$$
\max_{j:\ \bar a_{ij}<0}\frac{\bar c_j}{\bar a_{ij}}
\ \le\ \delta\ \le\
\min_{j:\ \bar a_{ij}>0}\frac{\bar c_j}{\bar a_{ij}},
$$ (eq-03-07-basic-range)

with an empty maximum read as $-\infty$ and an empty minimum as $+\infty$.
For a maximization the inequalities $\bar c_j(\delta)\le0$ give the same
formula with the roles of the two bounds exchanged. Inside the range the
solution $\vx^*$ does not move and $z^*$ changes by $\delta x_r^*$.
Outside, the nonbasic variable whose reduced cost turned negative enters,
and the primal simplex method continues.

:::{prf:example} Cost of a basic variable
:label: ex-03-07-basic-cost

In the running example, find the range of $c_1$ (now 3) over which the basis
stays optimal, and re-optimize for $c_1=5$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With $c_1$ as a symbol,

$$
\vpi^\T=(c_1,2,-3)\mB^{-1}=(5-c_1,\ 2-2c_1,\ 2c_1-5),
$$

and the reduced costs of $x_4,x_5,x_6$ are

$$
\begin{aligned}
\bar c_4&=c_1-2,\\
\bar c_5&=9-2c_1,\\
\bar c_6&=2+2c_1 .
\end{aligned}
$$

All three are nonnegative exactly when $2\le c_1\le\tfrac92$. The same range
follows from {eq}`eq-03-07-basic-range` with row 1 of
{numref}`tab-03-07-optimal-tableau`, $(\bar a_{14},\bar a_{15},\bar a_{16})=(-1,2,-2)$
and $(\bar c_4,\bar c_5,\bar c_6)=(1,3,8)$:
$\max\{1/(-1),\,8/(-2)\}=-1\le\delta\le3/2$. Inside the range
$z^*(c_1)=3c_1+2\cdot4-3\cdot2=2+3c_1$.

For $c_1=5$ we get $\vpi^\T=(0,-8,5)$ and
$(\bar c_4,\bar c_5,\bar c_6)=(3,-1,12)$, so $x_5$ enters. The last row of
the $\mT$-matrix is now $(-\vpi^\T,1)=(0,8,-5,1)$, while the cost of $x_5$ is
still 10, so

$$
\mT\begin{pmatrix}\mA_{\cdot5}\\c_5\end{pmatrix}
=\begin{pmatrix}-1&-2&2&0\\1&1&-1&0\\-1&0&1&0\\0&8&-5&1\end{pmatrix}
\begin{pmatrix}0\\-2\\-1\\10\end{pmatrix}
=\begin{pmatrix}2\\-1\\-1\\-1\end{pmatrix}.
$$

The ratio test again removes $x_1$, and the pivot gives the same basis
$(x_5,x_2,x_3)$ as in {prf:ref}`ex-03-07-nonbasic-cost`, with reduced costs
$(\tfrac12,0,0,\tfrac52,0,11)$ and
$\vx^*=(0,\tfrac{11}{2},\tfrac72,0,\tfrac32,0)$,
$z^*=\tfrac{31}{2}$. Making $x_1$ expensive pushes it out of the solution,
and $x_5$ takes over its role. HiGHS confirms $z^*=15.5$; at $c_1=\tfrac92$
both bases give $z=15.5$, and at $c_1=2$ the basis $(x_1,x_2,x_4)$ ties.

*(Earlier versions of this example multiply $\mT$ by the column
$(0,-2,-1,5)^\T$; the cost of $x_5$ is unchanged at 10, and only with 10
does the product equal $(2,-1,-1,-1)^\T$.)*
:::

Applying {eq}`eq-03-07-nonbasic-range` and {eq}`eq-03-07-basic-range` to every
variable produces the cost-ranging part of a solver's sensitivity report,
{numref}`tab-03-07-ranges`.

```{table} Ranges of the running example: each coefficient alone may vary in the interval shown without changing the optimal basis $(x_1,x_2,x_3)$.
:name: tab-03-07-ranges

| Coefficient | Current value | Range | Effect inside the range |
|---|---|---|---|
| $c_1$ (basic) | 3 | $[2,\ 9/2]$ | $z^*$ changes by $3\delta$ |
| $c_2$ (basic) | 2 | $[-1,\ 3]$ | $z^*$ changes by $4\delta$ |
| $c_3$ (basic) | −3 | $[-6,\ -5/2]$ | $z^*$ changes by $2\delta$ |
| $c_4$ (nonbasic) | −6 | $[-7,\ \infty)$ | none |
| $c_5$ (nonbasic) | 10 | $[7,\ \infty)$ | none |
| $c_6$ (nonbasic) | −5 | $[-13,\ \infty)$ | none |
| $b_1$ | 11 | $[7,\ 13]$ | $z^*$ changes by $\pi_1\delta=2\delta$ |
| $b_2$ | 6 | $[2,\ 15/2]$ | $z^*$ changes by $\pi_2\delta=-4\delta$ |
| $b_3$ | 13 | $[23/2,\ 17]$ | $z^*$ changes by $\pi_3\delta=\delta$ |
```

(sec-03-07-rhs)=
## Ranging a right-hand side and shadow prices

Change $b_r$ to $b_r+\delta$. The reduced costs do not depend on $\vb$, so
the basis stays dual feasible. Its basic values become

$$
\vx_B(\delta)=\mB^{-1}(\vb+\delta\ve_r)=\vx_B^*+\delta\,\mB^{-1}\ve_r ,
$$ (eq-03-07-rhs-values)

where $\mB^{-1}\ve_r$ is column $r$ of $\mB^{-1}$, read off the $\mT$-matrix.
The basis stays optimal while $\vx_B(\delta)\ge\vzero$:

$$
\max_{i:\ (\mB^{-1})_{ir}>0}\frac{-x_{B(i)}^*}{(\mB^{-1})_{ir}}
\ \le\ \delta\ \le\
\min_{i:\ (\mB^{-1})_{ir}<0}\frac{-x_{B(i)}^*}{(\mB^{-1})_{ir}} .
$$ (eq-03-07-rhs-range)

This rule is the same for minimization and maximization. Inside the range,

$$
\begin{aligned}
z^*(\delta)&=\vc_B^\T\vx_B(\delta)=z^*+\delta\,\vc_B^\T\mB^{-1}\ve_r\\
&=z^*+\pi_r\,\delta .
\end{aligned}
$$ (eq-03-07-shadow-price)

:::{prf:definition} Shadow price
:label: def-03-07-shadow-price

The **shadow price** of constraint $r$ is the rate of change of the optimal
value with respect to its right-hand side, $\partial z^*/\partial b_r$. When
the optimal basis is nondegenerate it equals the multiplier $\pi_r$, and
{eq}`eq-03-07-shadow-price` holds on the whole range {eq}`eq-03-07-rhs-range`, the
**validity range** of the shadow price.
:::

The shadow price is what one more unit of resource $r$ is worth at the
margin ({prf:ref}`thm-03-05-shadow-price`): the most one should pay for it, or the least
one should accept to give it up. Its validity range says how many units may
be bought or sold at that price; beyond it the basis changes and so does the
price, usually for the worse (diminishing returns, see
{ref}`sec-03-07-parametric`). At a degenerate optimum one end of the range can
be $\delta=0$, and the rates of change to the left and to the right may
then differ; the multiplier $\pi_r$ is only one of them.

:::{prf:example} Right-hand side of the first constraint
:label: ex-03-07-rhs

In the running example, (a) find the range of $b_1$ (now 11) over which the
basis stays optimal and the optimal value inside it; (b) re-optimize for
$b_1=15$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** Column 1 of $\mB^{-1}$ is $(-1,1,-1)^\T$, so by {eq}`eq-03-07-rhs-values`

$$
\vx_B(b_1)=\mB^{-1}\begin{pmatrix}b_1\\6\\13\end{pmatrix}
=\begin{pmatrix}14-b_1\\b_1-7\\13-b_1\end{pmatrix}\ge\vzero ,
$$

that is, $7\le b_1\le13$.

Inside the range $z^*(b_1)=3(14-b_1)+2(b_1-7)-3(13-b_1)=2b_1-11$, with slope
$\pi_1=2$: each extra unit of $b_1$ raises the minimum cost by 2. At the ends
of the range $x_2$ (at $b_1=7$) or $x_3$ (at $b_1=13$) drops to zero.

**(b)** For $b_1=15$, $\mT(\vb;0)=(-1,8,-2,-19)^\T$: the basis is dual
feasible but $x_1=-1$ and $x_3=-2$ are negative. The dual simplex method
({prf:ref}`alg-03-06-dual-simplex`) removes the most negative, $x_3$. Its row
of $\mB^{-1}\mA$ is row 3 of $\mB^{-1}$ times $\mA$,

$$
(-1,0,1)\,\mA=(0,0,1,2,-1,1),
$$

whose only negative entry is $\bar a_{35}=-1$; the dual ratio is
$\bar c_5/1=3$, so $x_5$ enters. After the pivot:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | RHS |
|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 0 | 2 | 3 | 0 | 0 | −5 |
| $x_2$ | 0 | 1 | −1 | −1 | 0 | −3 | 10 |
| $x_5$ | 0 | 0 | −1 | −2 | 1 | −1 | 2 |
| $-z$ | 0 | 0 | 3 | 7 | 0 | 11 | −25 |

Now $x_1=-5$ leaves, but its row $(1,0,2,3,0,0)$ has no negative entry: it
states $x_1+2x_3+3x_4=-5$, impossible for $\vx\ge\vzero$. **The LP with
$b_1=15$ has no feasible solution**, and HiGHS agrees. (The parametric
analysis of {ref}`sec-03-07-parametric` shows that the LP is feasible exactly
for $b_1\le40/3$.)

The new $\mT$-matrix is obtained as in {prf:ref}`thm-03-06-update`: the
bordered column of $x_5$ is $(2,-1,-1,3)^\T$ (its tableau column and
reduced cost), the pivot is in row 3, and the $4\times4$ eta matrix
{eq}`eq-03-06-eta` multiplies $\mT$ on the left:

$$
\begin{aligned}
\mT_{\text{new}}&=
\begin{pmatrix}1&0&2&0\\0&1&-1&0\\0&0&-1&0\\0&0&3&1\end{pmatrix}
\begin{pmatrix}-1&-2&2&0\\1&1&-1&0\\-1&0&1&0\\-2&4&-1&1\end{pmatrix}\\
&=\begin{pmatrix}-3&-2&4&0\\2&1&-2&0\\1&0&-1&0\\-5&4&2&1\end{pmatrix},
\end{aligned}
$$

which is indeed the inverse of the bordered basis of $(x_1,x_2,x_5)$.

*(Earlier versions of this example head the ratio row with
$\bar a_{2j}$; the pivot row is the row of $x_3$, row 3. They also stop
after choosing $x_5$; the continuation above shows that the problem is
infeasible.)*
:::

Doing the same as in part (a) of the solution of {prf:ref}`ex-03-07-rhs` for
all three rows gives the right-hand-side part of {numref}`tab-03-07-ranges`.
Note the sign of $\pi_2=-4$: raising $b_2$ lowers the minimum cost.

(sec-03-07-hundred-percent)=
## Simultaneous changes and the 100 percent rule

The ranges of {numref}`tab-03-07-ranges` hold for one coefficient at a time.
When several right-hand sides change at once, the exact test is simply
$\mB^{-1}(\vb+\Delta\vb)\ge\vzero$. A quick sufficient test, the 100 percent
rule of Bradley et al. (1977), needs only the table.

:::{prf:theorem} The 100 percent rule
:label: thm-03-07-hundred-percent

For each right-hand side $b_r$ that changes by $\delta_r$, let $\Delta_r>0$
be the allowed increase (if $\delta_r>0$) or the allowed decrease (if
$\delta_r<0$) from its range. If

$$
\sum_r\frac{\abs{\delta_r}}{\Delta_r}\le1,
$$ (eq-03-07-hundred-percent)

the optimal basis stays optimal, and
$z^*$ changes by $\sum_r\pi_r\delta_r$. Likewise, if several costs $c_j$
change by $\delta_j$ and the fractions of their allowed changes (from the
cost ranges) add up to at most one, the basis stays optimal, and $z^*$
changes by $\sum_j\delta_jx_j^*$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\kappa_r=\abs{\delta_r}/\Delta_r$ and $\kappa_0=1-\sum_r\kappa_r\ge0$.
The change vector is $\Delta\vb=\sum_r\kappa_r\,(\pm\Delta_r\ve_r)$, a
convex combination, with weights $\kappa_0,\kappa_1,\dots$, of the zero
change and of the single changes $\pm\Delta_r\ve_r$ to the ends of the
ranges. Each of these keeps $\mB^{-1}(\vb+\cdot)\ge\vzero$, and that set of
changes is convex (it is defined by linear inequalities), so it contains
$\Delta\vb$. (A change whose allowed amount is unlimited, $\Delta_r=\infty$,
counts as $\kappa_r=0$: the set contains the whole half-line
$\pm t\ve_r$, $t\ge0$, from every one of its points, so this change can be
added last.) The basis is unchanged, so
$z^*=\vc_B^\T\mB^{-1}(\vb+\Delta\vb)$ changes by
$\vpi^\T\Delta\vb=\sum_r\pi_r\delta_r$. For the costs, the same argument
applies to the set of cost vectors with $\bar c_j\ge0$ for all $j$, which is
defined by inequalities linear in $\vc$; $\vx^*$ does not move, so $z^*$
changes by $\sum_j\delta_jx_j^*$.
:::

The rule is sufficient, not necessary: a combination above 100 percent may
still keep the basis optimal. In the running example, raising $b_1$ from 11
to 12 uses $1/2$ of its allowed increase and raising $b_3$ from 13 to 14 uses
$1/4$ of its allowed increase: 75 percent in total, so the basis stays
optimal, with $\vx_B=\mB^{-1}(12,6,14)^\T=(4,4,2)^\T$ and
$z^*=11+2\cdot1+1\cdot1=14$ (HiGHS agrees). A uniform scaling of all costs
by a positive factor never changes the optimal basis, whatever the rule
says, because all reduced costs scale by the same factor
({prf:ref}`ex-08-04-a6-q2-electra`).

(sec-03-07-new-activity)=
## A new activity

A new variable $x_{n+1}$ with column $\mA_{\cdot n+1}$ and cost $c_{n+1}$ is
proposed: a new product, a new route, a new technology. The old optimal basis
is still a basis of the enlarged LP, with the new variable nonbasic at zero,
so it is primal feasible. It stays optimal if and only if the new reduced cost
is nonnegative,

$$
\bar c_{n+1}=c_{n+1}-\vpi^\T\mA_{\cdot n+1}\ge0 ,
$$ (eq-03-07-new-activity)

that is, if the activity costs at least the value of the resources it would
consume at the current shadow prices. This one inner product is how new
products are screened and priced in practice. If $\bar c_{n+1}<0$, the new
variable enters; its tableau column and reduced cost come from the
$\mT$-matrix, $\mT(\mA_{\cdot n+1};c_{n+1})$, and the primal simplex method
continues.

:::{prf:example} Two new activities
:label: ex-03-07-new-activity

A seventh variable is added to the running example. (a) Should it be used if
$\mA_{\cdot7}=(1,2,-3)^\T$ and $c_7=-7$? (b) Re-optimize for
$\mA_{\cdot7}=(3,-1,1)^\T$ and $c_7=4$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** The break-even cost is $\vpi^\T\mA_{\cdot7}=2-8-3=-9$, so
$\bar c_7=-7-(-9)=2>0$. The activity is not worth introducing; the optimal
solution and $z^*=11$ are unchanged.

**(b)** Now $\vpi^\T\mA_{\cdot7}=6+4+1=11$ and $\bar c_7=4-11=-7<0$: the new
activity is profitable and enters. Its bordered column is

$$
\mT\begin{pmatrix}3\\-1\\1\\4\end{pmatrix}
=\begin{pmatrix}-3+2+2\\3-1-1\\-3+0+1\\-6-4-1+4\end{pmatrix}
=\begin{pmatrix}1\\1\\-2\\-7\end{pmatrix}.
$$

The ratios are $3/1$ (row $x_1$) and $4/1$ (row $x_2$), so $x_1$ leaves. After
the pivot, with the new column appended:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | RHS |
|---|---|---|---|---|---|---|---|---|
| $x_7$ | 1 | 0 | 0 | −1 | 2 | −2 | 1 | 3 |
| $x_2$ | −1 | 1 | 0 | 2 | −3 | 0 | 0 | 1 |
| $x_3$ | 2 | 0 | 1 | 0 | 3 | −3 | 0 | 8 |
| $-z$ | 7 | 0 | 0 | −6 | 17 | −6 | 0 | 10 |

The objective dropped from 11 to $-10$, but two reduced costs are negative.
The column of $x_6$ is $(-2,0,-3)^\T\le\vzero$: increasing $x_6$ by $\theta$
and adjusting $x_7=3+2\theta$, $x_3=8+3\theta$ keeps every constraint
satisfied and lowers the cost by $6\theta$ without limit. **With this new
activity the LP is unbounded.** (One can check the ray directly: the
direction $\vd$ with $d_3=3$, $d_6=1$, $d_7=2$ and all other components zero
satisfies $\mA\vd=\vzero$ and has cost
$-9-5+8=-6<0$.) HiGHS also reports unboundedness. In a real model this
signals that the new activity, combined with $x_6$, creates a money machine,
which usually means a missing constraint or a wrong sign in its data.

*(Earlier versions of this example stop after the first pivot choice; the
continuation shows that the modified LP is unbounded.)*
:::

(sec-03-07-new-constraint)=
## A new constraint

A new constraint $\va^\T\vx\le\beta$ is handled as in
{ref}`sec-03-06-adding-constraint`. If the current optimum satisfies it, nothing
changes. Otherwise the constraint's slack enters the basis with the negative
value $\beta-\va^\T\vx^*$, the reduced costs stay as they were, and dual
simplex pivots restore feasibility. A new constraint can never improve the
optimal value.

:::{prf:example} Adding a constraint
:label: ex-03-07-new-constraint

Add $x_1+x_2\le6$ to the running example and re-optimize.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

At $\vx^*=(3,4,2,0,0,0)$ the left-hand side is $7>6$, so the optimum is cut
off. With the slack $x_7$, the new row is $x_1+x_2+x_7=6$. Subtracting the
rows of $x_1$ and $x_2$ of {numref}`tab-03-07-optimal-tableau` puts it in
canonical form:

$$
\begin{aligned}
&(1,1,0,0,0,0,1\mid6)\\
&\quad-(1,0,0,-1,2,-2,0\mid3)\\
&\quad-(0,1,0,1,-1,-2,0\mid4)\\
&=(0,0,0,0,-1,4,1\mid-1).
\end{aligned}
$$

The tableau is dual feasible with $x_7=-1$. Row $x_7$ leaves; its only
negative entry is in column $x_5$, with dual ratio $3/1=3$, so $x_5$ enters:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | RHS |
|---|---|---|---|---|---|---|---|---|
| $x_1$ | 1 | 0 | 0 | −1 | 0 | 6 | 2 | 1 |
| $x_2$ | 0 | 1 | 0 | 1 | 0 | −6 | −1 | 5 |
| $x_3$ | 0 | 0 | 1 | 2 | 0 | −3 | −1 | 3 |
| $x_5$ | 0 | 0 | 0 | 0 | 1 | −4 | −1 | 1 |
| $-z$ | 0 | 0 | 0 | 1 | 0 | 20 | 3 | −14 |

The tableau is primal and dual feasible: $\vx^*=(1,5,3,0,1,0)$ with
$z^*=3+10-9+10=14$. The new constraint is binding ($x_1+x_2=6$) and costs 3
units of objective, its shadow price, the reduced cost of its slack. HiGHS
gives the same optimum.
:::

(sec-03-07-coefficient-changes)=
### Changes of constraint coefficients

A change in a column $\mA_{\cdot j}$ of a **nonbasic** variable affects only
$\bar c_j$; it is handled exactly like a new activity with the modified
column, {eq}`eq-03-07-new-activity`. A change in the column of a **basic**
variable changes $\mB$ itself, and with it $\mB^{-1}$, $\vx_B$ and $\vpi$, so
both optimality conditions may fail at once. A practical way out is to add
the modified column as a new variable $x_j'$, priced with $\mT$, and to drive
the old $x_j$ out of the basis by giving it a prohibitive cost (or by
treating it as an artificial variable); the primal simplex method then
switches to $x_j'$. If both feasibilities are lost, a Phase I computation or
the artificial-constraint technique of {ref}`sec-03-06-artificial-constraint`
restores one of them. Changes of a coefficient in a new constraint are
covered by the procedure above.

(sec-03-07-parametric)=
## Parametric right-hand side

Ranging studies one change up to the first basis change. **Parametric
analysis** follows the optimal solution over the whole line
$\vb(t)=\vb+t\vd$, $t\in\R$, where $\vd$ is a fixed direction of change. On
each interval of $t$ one basis is optimal, $\vx_B(t)=\mB^{-1}\vb+t\mB^{-1}\vd$
is affine in $t$, and $z^*(t)=\vpi^\T\vb+t\,\vpi^\T\vd$ is affine with slope
$\vpi^\T\vd$. At an end point some basic variable reaches zero; one dual
simplex pivot on its row gives the basis of the next interval, or proves
that the LP becomes infeasible there.

:::{prf:theorem} The optimal value is convex in the right-hand side
:label: thm-03-07-value-convex

Let $\mA$ have full row rank. On the set of right-hand sides $\vb$ for which
the LP $\min\{\vc^\T\vx:\mA\vx=\vb,\vx\ge\vzero\}$ has a finite optimum, the
optimal value $z^*(\vb)$ is a convex, piecewise-linear function of $\vb$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Where the optimum is finite, strong duality gives
$z^*(\vb)=\max\{\vb^\T\vpi:\mA^\T\vpi\le\vc\}$. The dual feasible set does not
depend on $\vb$; because $\mA^\T$ has full column rank it contains no line,
so it has finitely many vertices $\vpi^1,\dots,\vpi^K$, and a finite maximum
of a linear function over it is attained at a vertex. Hence
$z^*(\vb)=\max_k\vb^\T\vpi^k$. A maximum of finitely many linear functions
is convex and piecewise linear.
:::

So along any line the slopes $\vpi^\T\vd$ increase from piece to piece: in a
minimization each further unit of a resource helps less (or hurts more) than
the one before. For a maximization $z^*(\vb)$ is concave. In the same way
(by the primal instead of the dual) the optimal value of a minimization is a
concave, piecewise-linear function of the cost vector $\vc$.

:::{prf:example} Optimal value as a function of $b_1$
:label: ex-03-07-parametric

In the running example, find $z^*(b_1)$ for all values of $b_1$, the other
data fixed.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Start from the optimal basis at $b_1=11$ and move in both directions.

* **$7\le b_1\le13$, basis $(x_1,x_2,x_3)$.** From {prf:ref}`ex-03-07-rhs`,
  $\vx_B=(14-b_1,\,b_1-7,\,13-b_1)$ and $z^*=2b_1-11$.
* **Above 13.** At $b_1=13$, $x_3$ reaches zero. For $b_1>13$ it would be
  negative, so the dual simplex method pivots on row $x_3$, whose only
  negative entry is in column $x_5$ (as in {prf:ref}`ex-03-07-rhs`). The new
  basis $(x_1,x_2,x_5)$ has $\vx_B=(40-3b_1,\,2b_1-20,\,b_1-13)$,
  $\vpi=(5,-4,-2)$ and $z^*=5b_1-50$, valid for $13\le b_1\le40/3$. At
  $b_1=40/3$, $x_1$ reaches zero, and its row $(1,0,2,3,0,0)$ has no
  negative entry: for $b_1>40/3$ the LP is infeasible.
* **Below 7.** At $b_1=7$, $x_2$ reaches zero. Its row
  $(0,1,0,1,-1,-2)$ has negative entries for $x_5$ (ratio $3/1=3$) and $x_6$
  (ratio $8/2=4$), so $x_5$ enters. The basis $(x_1,x_3,x_5)$ has
  $\vx_B=(b_1,\,20-2b_1,\,7-b_1)$, $\vpi=(-1,-7,4)$ and $z^*=10-b_1$, valid
  for $0\le b_1\le7$.
* **Below 0.** At $b_1=0$, $x_1$ reaches zero; in this basis its row is
  $(1,2,0,1,0,-6)$, whose only negative entry brings in $x_6$ (ratio
  $2/6$). The basis $(x_3,x_5,x_6)$ has
  $\vx_B=(20-\tfrac32b_1,\,7-\tfrac23b_1,\,-\tfrac16b_1)$,
  $\vpi=(-\tfrac43,-7,4)$ and $z^*=10-\tfrac43b_1$, feasible for every
  $b_1\le0$.

```{table} Optimal value of the running example as a function of $b_1$.
:name: tab-03-07-parametric

| Interval of $b_1$ | Optimal basis | Shadow price $\pi_1$ | $z^*(b_1)$ |
|---|---|---|---|
| $b_1\le0$ | $(x_3,x_5,x_6)$ | $-4/3$ | $10-\tfrac43b_1$ |
| $0\le b_1\le7$ | $(x_1,x_3,x_5)$ | $-1$ | $10-b_1$ |
| $7\le b_1\le13$ | $(x_1,x_2,x_3)$ | $2$ | $2b_1-11$ |
| $13\le b_1\le40/3$ | $(x_1,x_2,x_5)$ | $5$ | $5b_1-50$ |
| $b_1>40/3$ | none | — | infeasible |
```

The pieces join continuously ($z^*=10$ at $b_1=0$, 3 at 7, 15 at 13,
$50/3$ at $40/3$) and their slopes $-\tfrac43<-1<2<5$ increase, as
{prf:ref}`thm-03-07-value-convex` predicts; {numref}`fig-03-07-parametric-value`
shows the graph. Each slope is the shadow price $\pi_1$ of its basis, so the
shadow price of the first constraint is valid only on its own piece. HiGHS
reproduces $z^*$ at sample points of every piece (for example $z^*=7$ at
$b_1=3$, $9$ at $b_1=10$) and reports infeasibility for $b_1=14$.

```{figure} figures/03-07-parametric-value.svg
:name: fig-03-07-parametric-value
:width: 520px
:class: dark-light
:alt: Piecewise-linear convex graph of the optimal value against the first right-hand side, with four pieces and an infeasible region on the right

Optimal value $z^*(b_1)$ of the running example. Four bases are optimal on
four intervals; the slopes (shadow prices $\pi_1$) increase from $-4/3$ to
$5$, so the function is convex. For $b_1>40/3$ the LP is infeasible (shaded).
The dot marks the current value $b_1=11$, $z^*=11$.
```
:::

The slopes of $z^*(b_1)$ in {numref}`tab-03-07-parametric` are the shadow
prices $\pi_1$ of the successive optimal bases, so a shadow price is valid
only on its own piece of the optimal-value function.

A parametric analysis with all right-hand sides moving together is worked
in {prf:ref}`ex-08-04-a6-q6-parametric`.

(sec-03-07-reoptimisation)=
## Re-optimization when the basis changes

{numref}`tab-03-07-reoptimisation` collects the tests and the repair steps. In
every case the old optimal basis is the starting point, and $\mT$ provides
the new columns.

```{table} Post-optimal analysis of a minimization: what to compute after each kind of change, and how to re-optimize when the test fails.
:name: tab-03-07-reoptimisation

| Change | What to recompute | Basis stays optimal if | Otherwise |
|---|---|---|---|
| cost $c_r$, $x_r$ nonbasic | $\bar c_r$ | $c_r\ge\vpi^\T\mA_{\cdot r}$ | $x_r$ enters; primal simplex |
| cost $c_r$, $x_r$ basic | $\vpi$, all $\bar c_j$ | $\delta$ in {eq}`eq-03-07-basic-range` | primal simplex |
| right-hand side $b_r$ | $\vx_B=\mB^{-1}\vb$ | $\delta$ in {eq}`eq-03-07-rhs-range` | dual simplex |
| new activity | $\bar c_{n+1}$, then $\mT(\mA_{\cdot n+1};c_{n+1})$ | $\bar c_{n+1}\ge0$ | new variable enters; primal simplex |
| new constraint | $\va^\T\vx^*$ | $\va^\T\vx^*\le\beta$ | add row; dual simplex |
| column of a nonbasic variable | $\bar c_j$ | $\bar c_j\ge0$ | as a new activity |
| column of a basic variable | $\mB^{-1}$, $\vx_B$, $\vpi$ | both conditions | new variable plus primal or dual steps |
```

Two caveats apply. First, a re-optimized solution is optimal for the new
data, but the ranges computed for the old basis no longer apply; a new
analysis starts from the new basis. Second, all statements assume a unique,
nondegenerate optimal basis. With alternative optima (a zero reduced cost)
or degeneracy (a zero basic value), ranges and shadow prices depend on which
optimal basis a solver reports, and should be read with care.

## Summary and Common Mistakes

### Summary

- A basis is optimal when it is primal feasible {eq}`eq-03-07-primal-feasible`
  and dual feasible {eq}`eq-03-07-dual-feasible`
  ({prf:ref}`thm-03-07-basis-optimal`); costs and new columns affect only the
  second, right-hand sides and new rows only the first.
- A nonbasic cost may rise freely and fall to its break-even value
  $\vpi^\T\mA_{\cdot r}$ ({eq}`eq-03-07-nonbasic-range`); a basic cost moves all
  reduced costs, with the range {eq}`eq-03-07-basic-range` read from its tableau
  row ({prf:ref}`ex-03-07-nonbasic-cost`, {prf:ref}`ex-03-07-basic-cost`).
- A right-hand side may vary while $\vx_B^*+\delta\mB^{-1}\ve_r\ge\vzero$
  ({eq}`eq-03-07-rhs-range`); inside this validity range the optimal value
  changes at the shadow price $\pi_r$ ({prf:ref}`def-03-07-shadow-price`,
  {prf:ref}`ex-03-07-rhs`).
- Simultaneous changes are safe when their fractions of the allowed changes
  add up to at most one ({prf:ref}`thm-03-07-hundred-percent`).
- The $\mT$-matrix {eq}`eq-03-07-T-running` prices and updates any new or
  modified column with one product; new activities and constraints are
  handled by one test and, if needed, primal or dual simplex pivots
  ({prf:ref}`ex-03-07-new-activity`, {prf:ref}`ex-03-07-new-constraint`).
- The optimal value is convex and piecewise linear in $\vb$
  ({prf:ref}`thm-03-07-value-convex`); following it with dual simplex pivots
  gives the parametric analysis of {numref}`fig-03-07-parametric-value`.
- Geometrically, a right-hand side moves a constraint parallel to itself and
  a cost rotates the objective contour ({numref}`fig-03-07-rhs-ranging`,
  {numref}`fig-03-07-cost-ranging`).

### Common mistakes

- Treating a cost range as a range for the optimal *value*: inside the range
  of a basic cost the solution is fixed but $z^*$ moves by $\delta x_r^*$.
- Using the shadow price outside its validity range, for example valuing 500
  extra units of a resource whose range allows 200
  ({prf:ref}`ex-08-04-a6-q1-hidec`).
- Forgetting that changing a basic cost changes $\vpi$ and *every* nonbasic
  reduced cost, not only one.
- Updating the last row of $\mT$ but multiplying it by the old cost of an
  unchanged variable (or the new cost of the wrong one).
- Applying single-coefficient ranges to simultaneous changes without the
  100 percent rule or a direct check.
- Reading signs with the wrong convention: for a maximization optimality is
  $\bar c_j\le0$, and the inequalities of the cost ranges flip.
- Re-solving from scratch when one primal or dual simplex pivot from the old
  basis would do.

## References

- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear
  programming and network flows* (4th ed.). Wiley. (Chapter 6: sensitivity
  analysis, adding constraints and activities, parametric analysis.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Chapter 5: local sensitivity analysis,
  global dependence on the right-hand side and on the cost vector,
  parametric programming.)
- Bradley, S. P., Hax, A. C., & Magnanti, T. L. (1977). *Applied
  mathematical programming*. Addison-Wesley. (Chapter 3: sensitivity
  analysis and the 100 percent rule.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 7: sensitivity analysis and
  the 100 percent rule; Chapter 8: parametric linear programming.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Section 3.6: graphical and simplex-based sensitivity analysis;
  Section 4.5: post-optimal analysis; Section 7.5: parametric linear
  programming.)
