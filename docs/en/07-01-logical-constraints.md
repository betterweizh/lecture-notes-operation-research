# Logical Constraints and the Big-M Method

## Overview

Linear programs (LPs) describe decisions that can be scaled continuously: produce
3.4 tonnes, ship 12.5 trucks' worth, build 0.3 of a bridge. Many planning
questions are not like that. A bridge is built or not; a plant that opens
pays its fixed cost in full; a pump runs either not at all or above its
minimum rate; of two technical rules, *at least one* must hold. This chapter
shows how such yes/no decisions and logical conditions are written with
**binary variables** $\delta_i\in\{0,1\}$ and linear inequalities, turning
the model into a **mixed-integer linear program** (MILP). The central device
is the **big-$M$ method**: a binary variable switches a constraint on or off
by adding a constant $M$ to its right-hand side. We derive how large $M$ must
be (and why it should be no larger), and then build the standard logical
patterns from it: either–or constraints, if–then constraints, "$k$ of $n$"
conditions, fixed charges, semi-continuous variables and logical relations
between binary variables. The geometry of the either–or constraint explains
why some correct formulations are better than others. Two infrastructure
models, a project portfolio with precedence rules and a facility-opening
problem, close the chapter. It builds on LP modeling
({doc}`03-01-lp-model-building`), polyhedral geometry
({doc}`03-02-lp-geometry`) and convex hulls
({prf:ref}`def-02-01-convex-hull`).

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- explain how a binary variable and a constant $M$ switch an inequality on
  and off, and compute a valid and tight $M$ from variable bounds;
- explain why an unnecessarily large $M$ weakens the LP relaxation and
  causes numerical trouble;
- formulate either–or, if–then, "at least / at most / exactly $k$ of $m$",
  fixed-charge and semi-continuous conditions with binaries;
- translate logical statements about yes/no decisions into linear
  inequalities and check them by enumeration;
- compare formulations of the same condition by the strength of their LP
  relaxations, using the convex hull as the yardstick; and
- model project-selection and facility-opening problems with logical side
  conditions.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\delta_i\in\{0,1\}$ | binary (indicator) variable |
| $\Omega_0$ | set of points allowed by all *other* constraints of the model (usually a box $\vl\le\vx\le\vu$ or a bounded polyhedron) |
| $g_i(\vx)\le 0$ | constraint that is switched on or off by a binary variable |
| $M$, $M_i$ | big-$M$ constants: upper bounds of $g_i(\vx)$ on $\Omega_0$ |
| $\underline g$, $\underline g_i$ | lower bounds of $g(\vx)$, $g_i(\vx)$ on $\Omega_0$ |
| $\varepsilon$ | small positive tolerance used to model strict inequalities |
| $l_j$, $u_j$ | lower and upper bound of a variable $x_j$ |
| $F$, $F_i$ | fixed cost (set-up or opening cost), in currency units |
| $\conv S$ | convex hull of a set $S$ |
| $\omega$ | binary defined as a logical function (AND, OR) of other binaries |

(sec-07-01-milp)=
## Binary variables and mixed-integer programs

A **mixed-integer linear program** has a linear objective, linear
constraints, and some variables restricted to integer values:

$$
\begin{aligned}
\min_{\vx,\vdelta}\ & \vc^\T\vx+\vq^\T\vdelta \\
\st\ & \mA\vx+\mD\vdelta\le\vb, \\
& \vx\ge\vzero, \\
& \vdelta\in\{0,1\}^p .
\end{aligned}
$$ (eq-07-01-milp)

Here $\vdelta=(\delta_1,\dots,\delta_p)$ collects the binary
variables and $\vq$ their objective coefficients. Its feasible set
$\Omega$ is a finite union of polyhedra, one for each value of
$\vdelta$, and is in general not convex. Dropping the integrality
requirement, that is, replacing $\delta_i\in\{0,1\}$ by
$0\le\delta_i\le1$, gives the **LP relaxation**. Because its feasible set
contains $\Omega$, the LP relaxation of a minimization problem has an optimal
value no larger than that of the MILP: it is a **lower bound** (an upper
bound for maximization). Branch-and-bound solvers compute such bounds at
every node of their search tree, so the closer the relaxation is to the
MILP, the faster the search. Two formulations with the same integer
solutions can therefore behave very differently, and much of this chapter
is about writing the stronger one.

A binary variable has two natural meanings: a **decision** ("build the
bridge", $\delta=1$) or an **indicator** of whether a constraint is
enforced. The big-$M$ method links the two.

(sec-07-01-big-m)=
## The big-M method

Let $\Omega_0$ be the set of points allowed by all constraints of the model
except one, $g(\vx)\le0$, and suppose that $g$ is bounded above on
$\Omega_0$ by a known number $M$:

$$
g(\vx)\le M\quad\text{for all }\vx\in\Omega_0 .
$$ (eq-07-01-valid-m)

Then the constraint $g(\vx)\le M$ is **redundant** on $\Omega_0$: it cuts
nothing off. Introduce a binary variable $\delta$ and replace $g(\vx)\le0$
by

$$
\begin{gathered}
g(\vx)\le M(1-\delta),\\
\delta\in\{0,1\},\quad \vx\in\Omega_0 .
\end{gathered}
$$ (eq-07-01-bigm-relax)

For $\delta=1$ the right-hand side is $0$ and the constraint is enforced;
for $\delta=0$ it reads $g(\vx)\le M$, which every point of $\Omega_0$
satisfies, so the constraint is **relaxed** (it may or may not hold). In
words, $\delta=1\Rightarrow g(\vx)\le0$. The implication runs one way only:
$\delta=0$ does not force $g(\vx)>0$.

A constraint in the other direction, $g(\vx)\ge0$, is first rewritten as
$-g(\vx)\le0$:

$$
\begin{aligned}
-g(\vx)&\le M(1-\delta),\\
M&\ge\max_{\vx\in\Omega_0}\bigl(-g(\vx)\bigr) .
\end{aligned}
$$ (eq-07-01-bigm-ge)

An **equality** $g(\vx)=0$ cannot be switched off by a single big-$M$ term,
because no constant on the right-hand side turns an equality into a
redundant constraint. It is split into the two inequalities $g(\vx)\le0$ and
$-g(\vx)\le0$, each relaxed with the *same* binary and its own constant:
$g(\vx)\le M^+(1-\delta)$ and $-g(\vx)\le M^-(1-\delta)$, with $M^+$ and
$M^-$ upper bounds of $g$ and $-g$ on $\Omega_0$.

(sec-07-01-choosing-m)=
### Choosing M

Condition {eq}`eq-07-01-valid-m` says what a **valid** $M$ is: any number at
least $\max_{\vx\in\Omega_0}g(\vx)$. Too small a value is an error, not a
mere weakness: if $M<g(\bar\vx)$ for some $\bar\vx\in\Omega_0$, then the
point $\bar\vx$ is cut off even when the constraint is supposed to be
relaxed, and the model may lose its true optimum.

For a linear constraint $g(\vx)=\va^\T\vx-b$ over a box $\vl\le\vx\le\vu$
the maximum is attained term by term: $a_jx_j$ is largest at $u_j$ when
$a_j>0$ and at $l_j$ when $a_j<0$. Hence the smallest valid constant is

$$
M=\sum_{j:\,a_j>0}a_ju_j+\sum_{j:\,a_j<0}a_jl_j-b .
$$ (eq-07-01-m-from-bounds)

When $\Omega_0$ is a general polyhedron, $\max_{\Omega_0}g$ is itself an LP;
the box bound {eq}`eq-07-01-m-from-bounds` is valid but may be larger. If a
variable has no finite bound, no valid $M$ exists. Bounds are then derived
from the problem's physics (no plant ships more than total demand, no flow
exceeds a pipe's capacity), or the condition is handed to the solver as an
**indicator constraint**: most current MILP solvers accept a statement
"$\delta=1\Rightarrow\va^\T\vx\le b$" directly and enforce it by branching,
without any $M$.

Why not simply take $M=10^6$ and be safe? There are two reasons.

1. **Weak relaxation.** In the LP relaxation $\delta$ may be fractional. The
   constraint {eq}`eq-07-01-bigm-relax` is then satisfied by any $\vx$ with
   $\delta\le1-g(\vx)/M$: the larger $M$, the closer $\delta$ can stay to
   $1$ while the constraint is violated, or, read the other way, the
   smaller the fraction of $\delta$ that has to be "paid" to switch the
   constraint off. With a huge $M$ the relaxation behaves almost as if the
   logical condition were absent, its bound is poor, and branch and bound
   must explore many more nodes. {prf:ref}`ex-07-01-either-or` and the
   fixed-charge model of {ref}`sec-07-01-fixed-charge` show this numerically.
2. **Numerical error.** Solvers accept a value as integer when it lies
   within a tolerance of an integer, typically $10^{-6}$. With
   $M=10^{6}$ the value $\delta=1-10^{-6}$ counts as $1$, yet the constraint
   $g(\vx)\le M(1-\delta)=1$ then allows a violation of a whole unit. Large
   coefficients also make the LP bases ill-conditioned.

The rule is therefore: **take $M$ as the smallest valid bound**, computed
from {eq}`eq-07-01-m-from-bounds` or from a tighter argument, and use a
separate constant $M_i$ for each constraint.

(sec-07-01-reverse)=
### Forcing the binary: the reverse implication

Formulation {eq}`eq-07-01-bigm-relax` models $\delta=1\Rightarrow g(\vx)\le0$.
Sometimes the converse is needed: the binary must *record* whether the
constraint holds, $g(\vx)\le0\Rightarrow\delta=1$, or equivalently
$\delta=0\Rightarrow g(\vx)>0$. Strict inequalities cannot be imposed in an
LP, so a tolerance $\varepsilon>0$ is used: $\delta=0\Rightarrow
g(\vx)\ge\varepsilon$. With a lower bound $\underline g\le g(\vx)$ on $\Omega_0$
($\underline g<0$), this is

$$
g(\vx)\ge\varepsilon-(\varepsilon-\underline g)\,\delta .
$$ (eq-07-01-reverse)

For $\delta=0$ it gives $g(\vx)\ge\varepsilon$; for $\delta=1$ it gives
$g(\vx)\ge\underline g$, which always holds. Together, {eq}`eq-07-01-bigm-relax` and
{eq}`eq-07-01-reverse` make $\delta=1$ exactly when $g(\vx)\le0$, except that
points with $0<g(\vx)<\varepsilon$ are excluded. When $g$ takes only integer values on
$\Omega_0$, $\varepsilon=1$ is exact; otherwise $\varepsilon$ must exceed the solver's
feasibility tolerance.

(sec-07-01-either-or)=
## Either–or constraints

Two constraints $g_1(\vx)\le0$ and $g_2(\vx)\le0$ are **alternatives** when
at least one of them must hold. Typical cases: a machine processes job $A$
before job $B$ or job $B$ before job $A$; a building is placed either west of
a road or north of a river. With one binary per constraint,

$$
\begin{aligned}
g_1(\vx)&\le M_1(1-\delta_1),\\
g_2(\vx)&\le M_2(1-\delta_2),\\
\delta_1+\delta_2&\ge1,\qquad \delta_1,\delta_2\in\{0,1\} .
\end{aligned}
$$

At least one of the constraints is enforced. Since enforcing both is never
required, $\delta_1+\delta_2=1$ may be used instead, and then
$\delta_2=1-\delta_1$ can be eliminated:

$$
\begin{aligned}
g_1(\vx)&\le M_1(1-\delta),\\
g_2(\vx)&\le M_2\,\delta,\qquad \delta\in\{0,1\} .
\end{aligned}
$$ (eq-07-01-either-or)

With $\delta=1$ the first constraint holds, with $\delta=0$ the second.
A constraint that is a system of several inequalities (a polyhedron)
receives the same binary in each of its rows, each row with its own
constant.

Geometrically, the feasible set of an either–or condition is the **union**
of two sets, which is usually not convex. Any LP relaxation is convex, so
the best it can do is the **convex hull** of the union; how close it comes
depends on the constants.

:::{prf:example} Either a triangle or a square
:label: ex-07-01-either-or

On the box $\Omega_0=[0,4]^2$ the point $\vx=(x_1,x_2)$ must lie in the
triangle $T=\{x_1+x_2\le2\}$ or in the square $S=\{x_1\ge3,\ x_2\ge3\}$.
Maximize $-x_1+2x_2$. Formulate the condition with big-$M$ constraints,
find the integer optimum, and compare the bounds given by the LP relaxation
with the smallest valid constants and with $M=10$ for every row.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Formulation.* Let $\delta=1$ select $T$ and $\delta=0$ select $S$. The
triangle's row $x_1+x_2-2\le0$ is relaxed when $\delta=0$; the square's rows
$3-x_1\le0$ and $3-x_2\le0$ are relaxed when $\delta=1$:

$$
\begin{aligned}
x_1+x_2-2&\le M_T(1-\delta),\\
3-x_1&\le M_1\delta,\\
3-x_2&\le M_2\delta .
\end{aligned}
$$

By {eq}`eq-07-01-m-from-bounds` the smallest valid constants are
$M_T=4+4-2=6$ and $M_1=M_2=3-0=3$.

*Integer solution.* With $\delta=1$ the best point of $T$ is $(0,2)$ with
value $4$; with $\delta=0$ the best point of $S$ is $(3,4)$ with value $5$.
The optimum is $\vx^*=(3,4)$, $\delta^*=0$, value $5$.

*LP relaxations.* {numref}`fig-07-01-either-or` shows the projections onto
$(x_1,x_2)$ of the LP relaxations. The convex hull of $T\cup S$ is the
hexagon with vertices $(0,0)$, $(2,0)$, $(4,3)$, $(4,4)$, $(3,4)$, $(0,2)$
and area $10$. The relaxation with the tight constants is the larger hexagon
$(0,0)$, $(2,0)$, $(4,2)$, $(4,4)$, $(2,4)$, $(0,2)$ of area $12$, and its
LP optimum is $6$ at $\vx=(2,4)$ with $\delta=\tfrac13$: the point satisfies
$x_1+x_2=6\le2+6\cdot\tfrac23$ and $3-x_1=1\le3\cdot\tfrac13$, so one third
of a "triangle decision" lets it escape both sets. With $M_T=M_1=M_2=10$
the relaxation is the whole box (area $16$): the LP optimum is $8$ at
$(0,4)$, for example with $\delta=0.3$. Both MILPs return the same integer
optimum, $5$; only the bounds differ: $6$ (gap $20\%$) against $8$ (gap
$60\%$).

```{figure} figures/07-01-either-or.svg
:name: fig-07-01-either-or
:width: 620px
:class: dark-light
:alt: Union of a triangle and a square in a box, with its convex hull and the projections of two big-M relaxations

Either–or constraint of {prf:ref}`ex-07-01-either-or`. Left: the feasible set
$T\cup S$ (dark) and its convex hull (dashed outline). Right: the projection
of the LP relaxation with the tight constants $M_T=6$, $M_1=M_2=3$ (shaded
hexagon, area 12) and with $M=10$ (the whole box, area 16), with $T\cup S$
and its hull repeated for comparison. The optimum
$\vx^*=(3,4)$ of $\max\,-x_1+2x_2$ and the relaxed optimum $(2,4)$ are
marked.
```
:::

As the solution of {prf:ref}`ex-07-01-either-or` shows, even the tight
constants leave the relaxation strictly larger than the convex hull of
$T\cup S$. The convex hull itself can be reached by a larger formulation. Write
$\vx=\vy+\vz$, require $\vy\in\delta\,T$ and $\vz\in(1-\delta)S$, that is,
scale the right-hand sides (and the bounds) of each set by its binary:
$y_1+y_2\le2\delta$, $\vzero\le\vy\le4\delta\vone$, $z_1\ge3(1-\delta)$,
$z_2\ge3(1-\delta)$, $\vz\le4(1-\delta)\vone$. For $\delta\in\{0,1\}$ one of
the copies is zero and the other lies in its set; for fractional $\delta$
the point $\vx$ is a convex combination of a point of $T$ and a point of $S$.
This **disaggregated** (or convex-hull) formulation of disjunctive programming
(Balas, 1985) has the hull as its relaxation: its LP optimum in
{prf:ref}`ex-07-01-either-or` is $5$, equal to the integer optimum. It needs a copy of $\vx$ per
alternative, so it is used when the bound matters more than the model size
(Williams, 2013, Chapters 9–10).

(sec-07-01-if-then)=
## If–then constraints

An **if–then** condition reads: *if* $g_1(\vx)>0$, *then* $g_2(\vx)\le0$.
When $g_1(\vx)\le0$ nothing is required of $g_2$. Logically,
$A\Rightarrow B$ is the same as $(\text{not }A)\text{ or }B$, so the condition
is the either–or constraint

$$
g_1(\vx)\le0\quad\text{or}\quad g_2(\vx)\le0,
$$

and {eq}`eq-07-01-either-or` gives

$$
\begin{aligned}
g_1(\vx)&\le M_1(1-\delta),\\
g_2(\vx)&\le M_2\,\delta,\qquad \delta\in\{0,1\} .
\end{aligned}
$$ (eq-07-01-if-then)

Check of the logic: if $g_1(\vx)>0$, the first inequality excludes
$\delta=1$, so $\delta=0$ and the second inequality enforces
$g_2(\vx)\le0$. If $g_1(\vx)\le0$, the choice $\delta=1$ satisfies the first
inequality and relaxes the second, so $g_2$ is free. In set terms, the
condition describes

$$
\begin{aligned}
&\{\vx: g_1(\vx)>0,\ g_2(\vx)\le0\}\cup\{\vx: g_1(\vx)\le0\}\\
&\quad=\{\vx: g_2(\vx)\le0\}\cup\{\vx: g_1(\vx)\le0\},
\end{aligned}
$$

the union of two half-planes when $g_1,g_2$ are affine. For instance, with
$g_1=x_1+x_2-1$ and $g_2=1-x_1+x_2$ ("if $x_1+x_2>1$ then $x_1-x_2\ge1$"),
the feasible set is the plane minus the open wedge where $x_1+x_2>1$ and
$x_1-x_2<1$.

:::{prf:example} An if–then condition between binary decisions
:label: ex-07-01-if-then-binaries

Four yes/no decisions $\delta_1,\dots,\delta_4$ are linked by: if
$\delta_1=1$, then $\delta_2=\delta_3=\delta_4=0$ (choosing option 1
excludes options 2–4). Write the condition as linear constraints, first
from the general if–then rule {eq}`eq-07-01-if-then` and then as one row per
pair of excluded options, and compare the LP relaxations of the two
systems.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*From the general rule.* The premise $\delta_1=1$ is $g_1=\delta_1>0$ and
the conclusion is $g_2=\delta_2+\delta_3+\delta_4\le0$. On
$\Omega_0=\{0,1\}^4$ the bounds are $M_1=1$ and $M_2=3$, and
{eq}`eq-07-01-if-then` with an auxiliary binary $y$ reads
$\delta_1\le1-y$, $\delta_2+\delta_3+\delta_4\le3y$. Eliminating $y$ (the
two inequalities admit some $y\in\{0,1\}$ exactly when
$\delta_2+\delta_3+\delta_4\le3(1-\delta_1)$) gives the **aggregated**
constraint

$$
\delta_2+\delta_3+\delta_4\le3(1-\delta_1) .
$$

*(The premise must be written $\delta_1>0$, or $\delta_1\ge1$; a binary is
never greater than $1$, so a premise "$\delta_1>1$" would never fire.)*

*Disaggregated form.* The same condition is "at most one of $\delta_1$ and
$\delta_j$" for each $j$:

$$
\delta_1+\delta_j\le1,\qquad j=2,3,4 .
$$

Enumerating the $16$ binary vectors confirms that both systems accept
exactly the vectors with $\delta_1=0$ or $\delta_2=\delta_3=\delta_4=0$.
Their relaxations differ: the fractional point
$(\delta_1,\dots,\delta_4)=(\tfrac23,1,0,0)$ satisfies the aggregated
constraint ($1\le3\cdot\tfrac13$) but violates $\delta_1+\delta_2\le1$. Every
point of the disaggregated relaxation satisfies the aggregated one (add the
three inequalities).
:::

As the solution of {prf:ref}`ex-07-01-if-then-binaries` shows, the
disaggregated form is never weaker than the aggregated one and is preferred,
although it uses three rows instead of one.

(sec-07-01-k-of-n)=
## At least, at most and exactly k of m constraints

Let $g_i(\vx)\le0$, $i=1,\dots,m$, be constraints with valid constants
$M_i$. Requiring that **at least $k$** of them hold is a direct
generalization of the either–or case:

$$
\begin{aligned}
g_i(\vx)&\le M_i(1-\delta_i),\quad i=1,\dots,m,\\
\sum_{i=1}^m\delta_i&\ge k,\\
\delta_i&\in\{0,1\},\quad i=1,\dots,m .
\end{aligned}
$$ (eq-07-01-k-of-m)

(Either–or is $k=1$, $m=2$.) A single constant $M=\max_iM_i$ is valid for
every row but, by {ref}`sec-07-01-choosing-m`, individual constants are better.

"**At most $k$** hold" and "**exactly $k$** hold" are different in kind.
Replacing $\ge k$ by $\le k$ or $=k$ in {eq}`eq-07-01-k-of-m` limits only the
number of constraints that are *enforced*; the other constraints are relaxed,
not violated, and may hold anyway. So $\sum_i\delta_i\le k$ is no
restriction at all ($\vdelta=\vzero$ is allowed), and
$\sum_i\delta_i=k$ still means "at least $k$". To count the constraints that
are satisfied, each binary must also be forced to $1$ when its constraint
holds, through the reverse implication {eq}`eq-07-01-reverse`:

$$
\begin{aligned}
g_i(\vx)&\le M_i(1-\delta_i),\quad i=1,\dots,m,\\
g_i(\vx)&\ge\varepsilon-(\varepsilon-\underline g_i)\,\delta_i,\quad i=1,\dots,m,\\
\sum_{i=1}^m\delta_i&\le k\quad(\text{or }=k),\\
\delta_i&\in\{0,1\},\quad i=1,\dots,m .
\end{aligned}
$$ (eq-07-01-at-most-k)

Now $\delta_i=1$ exactly when $g_i(\vx)\le0$ (up to the tolerance
$\varepsilon$), and the cardinality constraint counts satisfied
constraints. Conditions of the "at most" kind are rare in practice (they ask
the model to violate rules on purpose); they arise, for example, when at
most $k$ service standards may be breached.

(sec-07-01-fixed-charge)=
## Fixed charges and semi-continuous variables

A **fixed-charge** cost has a set-up part $F>0$ that is paid as soon as an
activity is used, and a variable part $c\ge0$ per unit:

$$
C(x)=\begin{cases}0, & x=0,\\ F+cx, & 0<x\le u .\end{cases}
$$

$C$ is discontinuous at $0$ and concave on $[0,u]$, so it cannot appear in an
LP. With a binary $\delta$ ("the activity is used") the cost becomes
$F\delta+cx$, and the link between $x$ and $\delta$ is the big-$M$ constraint
$x\le M\delta$ (if $\delta=0$ then $x\le0$). The smallest valid constant is
the bound $u$ itself, which merges with the capacity constraint:

$$
\begin{aligned}
\min_{x,\delta}\ & F\delta+cx\\
\st\ & 0\le x\le u\,\delta,\\
& \delta\in\{0,1\} .
\end{aligned}
$$ (eq-07-01-fixed-charge)

No constraint forces $\delta=0$ when $x=0$: minimization does that, because
$\delta=1$ costs $F$ for nothing.

In the LP relaxation, $\delta$ only needs to reach $x/u$, so the relaxed cost
of producing $x$ is $(F/u+c)\,x$: the fixed charge is spread over the
capacity as if it were a variable cost. {numref}`fig-07-01-fixed-charge` shows
this for $F=\$500$, $c=\$20$ per tonne and $u=100$ tonnes. Producing $x=60$
tonnes truly costs $500+20\cdot60=\$1700$; the relaxation charges
$(5+20)\cdot60=\$1500$. With the needlessly large constant $M=1000$ the
relaxation charges only $(0.5+20)\cdot60=\$1230$: the lower bound loses
$94\%$ of the fixed charge instead of $40\%$.

```{figure} figures/07-01-fixed-charge.svg
:name: fig-07-01-fixed-charge
:width: 520px
:class: dark-light
:alt: Fixed-charge cost with a jump at zero and the two straight lines of its LP relaxations

Fixed-charge cost $C(x)$ with $F=\$500$, $c=\$20$/t and capacity
$u=100$ t (solid; filled dot: $C(0)=0$; open circle: the value $F$
that $C$ approaches as $x\downarrow0$ but does not take). The LP relaxation with $M=u$
replaces it by the line $25x$, the convex envelope of $C$ on $[0,u]$
(dashed); with $M=1000$ the line is $20.5x$ (dotted). At $x=60$ t the three
values are \$1700, \$1500 and \$1230.
```

A variable is **semi-continuous** if it is either $0$ or between a positive
lower bound and an upper bound: $x\in\{0\}\cup[l,u]$ with $0<l\le u$. A
pump either is off or delivers between $l=20$ and $u=100$ liters per second;
a supplier accepts either no order or at least its minimum order quantity.
With $\delta=1$ meaning "on",

$$
l\,\delta\le x\le u\,\delta,\qquad \delta\in\{0,1\} .
$$ (eq-07-01-semicontinuous)

For $\delta=0$ both bounds are $0$; for $\delta=1$ they are $l$ and $u$. A
fixed charge and a minimum level often appear together, and
{eq}`eq-07-01-semicontinuous` then carries both. Some solvers accept
semi-continuous variables directly and branch on the disjunction
$x=0$ or $x\ge l$ without an explicit binary.

(sec-07-01-logic)=
## Logical relations between binary variables

Once every yes/no statement has a binary variable, logical relations between
the statements become linear inequalities between the binaries. Let
$\delta_A,\delta_B,\delta_C$ indicate statements $A$, $B$, $C$ (decisions,
or constraints made to hold through {eq}`eq-07-01-bigm-relax`). Three rules
cover everything:

* negation: "not $A$" is $1-\delta_A$;
* a disjunction (OR) of literals holds when their sum is at least $1$;
* a conjunction (AND) of conditions is the list of their constraints.

Any logical condition can be put in **conjunctive normal form**, an AND of
OR-clauses, and each clause $\ell_1\lor\dots\lor\ell_r$ becomes
$\sum_{\text{positive }\ell}\delta+\sum_{\text{negated }\ell}(1-\delta)\ge1$.
For example, $A\Rightarrow B$ is the clause $\lnot A\lor B$, that is
$(1-\delta_A)+\delta_B\ge1$, or $\delta_B\ge\delta_A$. {numref}`tab-07-01-logic`
lists the common cases; every row was checked against the truth table by
enumerating the $8$ binary triples.

```{list-table} Logical conditions on binary variables and their linear forms
:name: tab-07-01-logic
:header-rows: 1
:widths: 40 60

* - Condition
  - Linear constraint(s)
* - if $A$ then $B$
  - $\delta_B\ge\delta_A$
* - if not $A$ then $B$ ($A$ or $B$)
  - $\delta_A+\delta_B\ge1$
* - at most one of $A$, $B$
  - $\delta_A+\delta_B\le1$
* - exactly one of $A$, $B$
  - $\delta_A+\delta_B=1$
* - $A$ if and only if $B$
  - $\delta_A=\delta_B$
* - if ($A$ or $C$) then $B$
  - $\delta_B\ge\delta_A$ and $\delta_B\ge\delta_C$ (weaker single row: $2\delta_B\ge\delta_A+\delta_C$)
* - if ($A$ and $C$) then $B$
  - $\delta_B\ge\delta_A+\delta_C-1$
* - if $A$ then ($B$ or $C$)
  - $\delta_B+\delta_C\ge\delta_A$
* - if $A$ then ($B$ and $C$)
  - $\delta_B\ge\delta_A$ and $\delta_C\ge\delta_A$ (weaker single row: $\delta_B+\delta_C\ge2\delta_A$)
```

*(A single row $\delta_B\ge\delta_A+\delta_C$ is sometimes given for "if
($A$ or $C$) then $B$"; it is wrong, because it forbids $\delta_A=\delta_C=1$
altogether: the right-hand side would be $2$.)* Where the table offers a
single row as an alternative, it has the same integer solutions but a weaker
relaxation, as in {prf:ref}`ex-07-01-if-then-binaries`.

A new binary can also be *defined* as a logical function of others. For
$\omega=\delta_A\land\delta_B$ (AND),

$$
\begin{aligned}
\omega&\le\delta_A,\\
\omega&\le\delta_B,\\
\omega&\ge\delta_A+\delta_B-1,
\end{aligned}
$$ (eq-07-01-and)

for $\omega=\delta_A\lor\delta_B$ (OR), $\omega\ge\delta_A$, $\omega\ge\delta_B$,
$\omega\le\delta_A+\delta_B$; and for the exclusive or, $\omega\ge\delta_A-\delta_B$,
$\omega\ge\delta_B-\delta_A$, $\omega\le\delta_A+\delta_B$, $\omega\le2-\delta_A-\delta_B$.
Enumeration of the four pairs $(\delta_A,\delta_B)$ confirms each definition.
The AND variable is the product $\delta_A\delta_B$, the starting point of
{doc}`07-02-linearization-techniques`.

(sec-07-01-infrastructure)=
## Infrastructure applications

The two examples below combine the patterns of this chapter. Their optimal
solutions were found with a MILP solver and confirmed by enumerating all
candidate decisions.

:::{prf:example} Selecting transport projects
:label: ex-07-01-project-selection

A city has a capital budget of \$120 million for six candidate projects
({numref}`tab-07-01-projects`), with benefits measured as present values of
net benefits.

```{list-table} Candidate projects (amounts in \$ million)
:name: tab-07-01-projects
:header-rows: 1

* - $j$
  - Project
  - Cost
  - Benefit
* - 1
  - river bridge, design A
  - 40
  - 70
* - 2
  - river bridge, design B
  - 30
  - 52
* - 3
  - ring road
  - 35
  - 50
* - 4
  - interchange on the ring road
  - 15
  - 28
* - 5
  - tram line
  - 45
  - 75
* - 6
  - bus rapid transit (BRT) line
  - 20
  - 33
```

The rules are: (i) at most one bridge design is built; (ii) the interchange
can be built only together with the ring road; (iii) at least two of the
three corridor projects (ring road, tram, BRT) are built; (iv) if both the
tram and the BRT line are built, a bridge must be built to carry one of
them. Which projects maximize the total benefit? Compare the answer with
the selection the budget alone would allow and with the LP relaxation.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Formulation.* With $\delta_j=1$ if project $j$ is selected,

$$
\begin{aligned}
\max_{\vdelta}\ & 70\delta_1+52\delta_2+50\delta_3+28\delta_4\\
& \quad+75\delta_5+33\delta_6\\
\st\ & 40\delta_1+30\delta_2+35\delta_3+15\delta_4\\
& \quad+45\delta_5+20\delta_6\le120, \\
& \delta_1+\delta_2\le1, \\
& \delta_4\le\delta_3, \\
& \delta_3+\delta_5+\delta_6\ge2, \\
& \delta_1+\delta_2\ge\delta_5+\delta_6-1, \\
& \vdelta\in\{0,1\}^6 .
\end{aligned}
$$

Rule (ii) is a **precedence** constraint ("if 4 then 3", first row of
{numref}`tab-07-01-logic`); rule (iii) is "at least $2$ of $3$"; rule (iv) is
"if (5 and 6) then (1 or 2)", obtained from the AND row with
$\delta_1+\delta_2$ in place of $\delta_B$ (rule (i) keeps this sum binary).

*Solution.* Of the $64$ binary vectors, $12$ satisfy all constraints. The
best selects bridge A, the ring road and the tram,
$\vdelta^*=(1,0,1,0,1,0)$, at a cost of exactly \$120 million and
a benefit of \$195 million. Without rules (i)–(iv) the budget alone would
allow bridge A, the interchange, the tram and the BRT line (benefit
\$206 million), which builds an interchange on a road that does not exist.
The LP relaxation gives \$201.4 million, for instance with
$\vdelta=(1,0,0.3,0.3,1,1)$: fractional projects fill the budget
and satisfy rule (iii) with $0.3+1+1\ge2$.
:::

:::{prf:example} Opening recycling facilities
:label: ex-07-01-facility

Four candidate sites can host a recycling facility serving five districts.
Site $i$ has an annual fixed cost $F_i$ (\$ thousand per year) and a
capacity $u_i$ (thousand tonnes per year); district $j$ produces $d_j$
thousand tonnes per year; moving waste from district $j$ to site $i$ costs
$c_{ij}$ dollars per tonne ({numref}`tab-07-01-facility`). Which sites should
open, and how strong is the LP relaxation of the model?

```{list-table} Data of the facility-opening problem
:name: tab-07-01-facility
:header-rows: 1

* - Site $i$
  - $c_{i1}$
  - $c_{i2}$
  - $c_{i3}$
  - $c_{i4}$
  - $c_{i5}$
  - $F_i$
  - $u_i$
* - 1
  - 4
  - 6
  - 9
  - 10
  - 8
  - 300
  - 90
* - 2
  - 7
  - 3
  - 5
  - 9
  - 10
  - 250
  - 70
* - 3
  - 10
  - 8
  - 4
  - 3
  - 6
  - 280
  - 80
* - 4
  - 9
  - 10
  - 7
  - 5
  - 3
  - 220
  - 60
* - $d_j$
  - 30
  - 25
  - 35
  - 20
  - 30
  -
  -
```
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Formulation.* Let $x_{ij}$ (thousand tonnes per year) be the waste sent from
district $j$ to site $i$ and $\delta_i=1$ if site $i$ opens. A unit cost of
\$1 per tonne on a thousand tonnes is \$1 thousand, so the objective is in
\$ thousand per year:

$$
\begin{aligned}
\min_{\vx,\vdelta}\ & \sum_{i=1}^4F_i\delta_i+\sum_{i=1}^4\sum_{j=1}^5c_{ij}x_{ij}\\
\st\ & \sum_{i=1}^4x_{ij}=d_j, & j&=1,\dots,5,\\
& \sum_{j=1}^5x_{ij}\le u_i\,\delta_i, & i&=1,\dots,4,\\
& x_{ij}\ge0, & i&=1,\dots,4,\ j=1,\dots,5,\\
& \delta_i\in\{0,1\}, & i&=1,\dots,4 .
\end{aligned}
$$

The capacity row is the fixed-charge link {eq}`eq-07-01-fixed-charge`: a closed
site handles nothing.

*Solution.* Total waste is $140$ thousand tonnes per year, more than any
single site can take, so at least two sites open. Solving the transport LP
for each feasible set of open sites gives the totals below (fixed plus
transport cost, \$ thousand per year); sites $\{2,4\}$ have only $130$ of
capacity.

| Open sites | {1,2} | {1,3} | {1,4} | {2,3} | {3,4} | {1,2,3} | {1,2,4} | {1,3,4} | {2,3,4} | all |
|---|---|---|---|---|---|---|---|---|---|---|
| Total cost | 1350 | 1240 | 1275 | **1200** | 1260 | 1410 | 1330 | 1360 | 1325 | 1535 |

The optimum opens sites 2 and 3: fixed cost $250+280=530$, transport cost
$670$, total \$1 200 thousand per year. Site 2 takes districts 1 and 2
($30$ and $25$) and $5$ of district 3; site 3 takes the remaining $30$ of
district 3 and districts 4 and 5, which fills it to its capacity of $80$.
The MILP solver returns the same solution.

*Strength of three formulations.* The LP relaxation bound depends on how the
link between $x_{ij}$ and $\delta_i$ is written:

1. with the valid but loose constant $M=140$ (total waste) in
   $\sum_jx_{ij}\le140\,\delta_i$, plus the capacity rows
   $\sum_jx_{ij}\le u_i$, the LP bound is $751.07$;
2. with $M=u_i$ as in the model above, it is $976.79$;
3. adding the **disaggregated** links $x_{ij}\le\min(d_j,u_i)\,\delta_i$
   (if a site is closed, each single flow is zero) raises it to $1168.46$.

The gaps to the optimum $1200$ are $37\%$, $19\%$ and $2.6\%$. All three
formulations have the same integer solutions; the third is the one to give a
solver.
:::

## Summary and Common Mistakes

### Summary

- A binary $\delta$ switches $g(\vx)\le0$ on and off through
  $g(\vx)\le M(1-\delta)$ ({eq}`eq-07-01-bigm-relax`), valid when
  $M\ge\max_{\Omega_0}g$ ({eq}`eq-07-01-valid-m`); for linear $g$ over a box the
  tight value is {eq}`eq-07-01-m-from-bounds`. Equalities are split into two
  inequalities; the reverse implication needs {eq}`eq-07-01-reverse`.
- Larger $M$ gives weaker LP relaxations and numerical errors; indicator
  constraints avoid $M$ when no bound is known.
- Either–or and if–then conditions use one binary
  ({eq}`eq-07-01-either-or`, {eq}`eq-07-01-if-then`); "at least $k$ of $m$" is
  {eq}`eq-07-01-k-of-m`, while "at most / exactly $k$" needs
  {eq}`eq-07-01-at-most-k`.
- Fixed charges and semi-continuous variables are modeled by
  $x\le u\delta$ and $l\delta\le x\le u\delta$ ({eq}`eq-07-01-fixed-charge`,
  {eq}`eq-07-01-semicontinuous`); the relaxation spreads the fixed charge as
  $F/u$ per unit ({numref}`fig-07-01-fixed-charge`).
- Logical relations between binaries follow from negation, clauses and
  conjunctions ({numref}`tab-07-01-logic`, {eq}`eq-07-01-and`).
- The feasible set of a disjunction is a union of polyhedra; its convex hull
  is the ideal relaxation ({numref}`fig-07-01-either-or`), reached by the
  disaggregated formulation; tight constants and disaggregated rows move
  big-$M$ models toward it ({prf:ref}`ex-07-01-facility`).

### Common mistakes

- **$M$ too small.** A constant below $\max_{\Omega_0}g$ silently cuts off
  feasible points. Derive $M$ from bounds with {eq}`eq-07-01-m-from-bounds`;
  never guess it.
- **$M$ too large.** It is valid but weakens the relaxation and invites
  tolerance errors; {prf:ref}`ex-07-01-either-or` and
  {prf:ref}`ex-07-01-facility` show bounds degrading as $M$ grows.
- **Reading the implication backwards.** $g(\vx)\le M(1-\delta)$ means
  $\delta=1\Rightarrow g(\vx)\le0$ only. If the model needs
  $g(\vx)\le0\Rightarrow\delta=1$, add {eq}`eq-07-01-reverse`.
- **"At most $k$" with $\sum\delta_i\le k$ alone.** Without the reverse
  implication it limits enforced constraints, not satisfied ones.
- **Strict inequalities.** $g(\vx)>0$ cannot be imposed; use
  $g(\vx)\ge\varepsilon$ with $\varepsilon$ above the solver tolerance (or
  $1$ for integer data).
- **Aggregating logical rows.** One row for several implications is correct
  but weaker; prefer the disaggregated rows.

## References

- Balas, E. (1985). Disjunctive programming and a hierarchy of relaxations
  for discrete optimization problems. *SIAM Journal on Algebraic and Discrete
  Methods, 6*(3), 466–486.
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Sections 10.1–10.2: modeling with
  integer variables, fixed charges, disjunctive constraints, facility
  location, strength of formulations.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Sections 12.2–12.3: binary decision
  variables, either–or, $K$ out of $N$, fixed charges.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Section 9.1: integer-programming applications, fixed-charge and
  either–or constraints.)
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapters 9–10: building integer programming models,
  logical conditions and zero–one variables, disjunctions, alternative
  formulations and the strength of their relaxations.)
