# Geometry of Linear Programs

## Overview

This chapter looks at linear programs (LPs) as geometric objects. With two
variables the feasible region can be drawn, and the **graphical method**
({prf:ref}`alg-03-02-graphical`) finds the optimum by sliding a line of
constant objective value across it; applied to the product mix of
{doc}`03-01-lp-model-building` ({prf:ref}`ex-03-02-graphical`) it shows that the
optimum sits at a corner, that some constraints are binding and others
have slack, and that rounding a fractional optimum is not the same as
solving the integer problem. The four possible outcomes of an LP (unique
optimum, alternative optima, unbounded, infeasible) are seen in
{numref}`fig-03-02-graphical-cases`. To go beyond two variables the chapter
turns to algebra. Every LP can be written in **standard form**
({prf:ref}`prop-03-02-standard-form`), its feasible set is a polyhedron, and
the corners of that polyhedron are exactly the **basic feasible solutions**
of the equation system ({prf:ref}`thm-03-02-extreme-bfs`). The **fundamental
theorem of linear programming** ({prf:ref}`thm-03-02-fundamental`) states that
an LP with an optimal solution has an optimal basic feasible solution, so
the search can be restricted to finitely many points. Finally, writing the
LP in terms of a basis gives the **reduced form**
{eq}`eq-03-02-reduced-form` and a first optimality test
({prf:ref}`thm-03-02-optimality-test`), the starting point of the simplex
method in {doc}`03-03-simplex-method`. The chapter uses linear
independence, rank and inverses from {doc}`01-01-linear-algebra` and
convex sets, polyhedra and extreme points from {doc}`02-01-convex-sets`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- solve a two-variable LP graphically and recognize the cases of a unique
  optimum, alternative optima, unboundedness and infeasibility;
- compute slacks and surpluses and identify binding constraints;
- convert any LP into standard form;
- explain why the feasible set of an LP is a convex polyhedron;
- compute the basic solutions of $\mA\vx=\vb$, classify them as feasible
  or infeasible, and recognize degenerate ones;
- prove that extreme points and basic feasible solutions coincide, and
  state and prove the fundamental theorem of LP; and
- write an LP in reduced form with respect to a basis and apply the
  reduced-cost optimality test.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $P$ | feasible set (polyhedron) of an LP, $P=\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ in standard form |
| $s_i$ | slack ($\le$ row) or surplus ($\ge$ row) of constraint $i$ |
| $\mB$, $\mN$ | basis matrix ($m$ linearly independent columns of $\mA$) and the matrix of the remaining columns |
| $\vx_B$, $\vx_N$; $\vc_B$, $\vc_N$ | basic and nonbasic variables; their objective coefficients |
| $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ | reduced cost of $x_j$, with $\vpi^\T=\vc_B^\T\mB^{-1}$ |
| $J(\vx)$ | support of $\vx$, the index set $\{j:x_j>0\}$ |

(sec-03-02-graphical)=
## The graphical method

An LP in two variables $x_1,x_2$ can be solved by drawing. Each linear
inequality $a_1x_1+a_2x_2\le b$ holds on one side of the line
$a_1x_1+a_2x_2=b$, so the feasible region is the polygon (possibly
unbounded) where all these half-planes overlap. The points with a given
objective value $z=c_1x_1+c_2x_2$ form a line, a **contour** (iso-profit
or iso-cost line); changing $z$ moves the line parallel to itself, in the
direction of the vector $(c_1,c_2)$ as $z$ increases.

:::{prf:algorithm} Graphical method for a two-variable LP
:label: alg-03-02-graphical

**Input:** an LP in the variables $x_1,x_2$.

1. **Draw the constraints.** For each constraint, draw the line where it
   holds with equality and determine the side on which it holds (test a
   point such as the origin).
2. **Shade the feasible region**, the intersection of all these
   half-planes and the nonnegativity quadrant. If it is empty, stop: the
   LP is infeasible.
3. **Draw a contour** $c_1x_1+c_2x_2=z_0$ for some convenient $z_0$.
4. **Slide the contour** parallel to itself in the improving direction
   ($+(c_1,c_2)$ for max, $-(c_1,c_2)$ for min) as far as it still touches
   the region. If it can be moved indefinitely, stop: the LP is unbounded.
5. **Read the optimum.** The last points touched are optimal. Find their
   coordinates by solving the two equations of the lines that meet there,
   and evaluate $z$.
:::

:::{prf:example} Graphical solution of the product mix
:label: ex-03-02-graphical

Solve the product-mix LP {eq}`eq-03-01-product-mix` of
{prf:ref}`ex-03-01-product-mix`,

$$
\begin{aligned}
\max\ & z=500x_1+450x_2 \\
\st\ & 6x_1+5x_2\le60 \\
& 10x_1+20x_2\le150 \\
& x_1\le8 \\
& x_1,x_2\ge0,
\end{aligned}
$$

with $x_1,x_2$ in hundreds of cases per week. Find the slack of every
constraint at the optimum, and compare the LP optimum with the best plan in
whole hundreds of cases.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Region.** The machine line $6x_1+5x_2=60$ meets the axes at $(10,0)$ and
$(0,12)$; the storage line $10x_1+20x_2=150$ at $(15,0)$ and $(0,7.5)$; the
demand line is $x_1=8$. All three constraints hold at the origin, so the
region lies on the origin's side of each line. It is the pentagon with
vertices

$$
\begin{gathered}
(0,0),\quad (8,0),\quad (8,\tfrac{12}{5}),\\
(\tfrac{45}{7},\tfrac{30}{7}),\quad (0,\tfrac{15}{2}),
\end{gathered}
$$

shown in {numref}`fig-03-02-product-mix-region`. The vertex
$(\tfrac{45}{7},\tfrac{30}{7})$ is the intersection of the machine and
storage lines: multiplying $6x_1+5x_2=60$ by 4 and subtracting
$10x_1+20x_2=150$ gives $14x_1=90$, so $x_1=\tfrac{45}{7}$ and
$x_2=(60-6\cdot\tfrac{45}{7})/5=\tfrac{30}{7}$. The lines $x_1=8$ and
$6x_1+5x_2=60$ meet at $(8,\tfrac{12}{5})$.

**Contours.** The contours $500x_1+450x_2=z$ have slope $-\tfrac{10}{9}$.
The contour $z=2250$ crosses the region; $z=4000$ passes through the vertex
$(8,0)$. Moving further, the last point of the region that a contour
touches is $(\tfrac{45}{7},\tfrac{30}{7})$. The contour slope
$-\tfrac{10}{9}$ lies between the slopes of the two lines that meet there,
$-\tfrac65$ (machine) and $-\tfrac12$ (storage), which is why this corner,
and no other, is touched last. Hence

$$
\begin{aligned}
x_1^*&=\tfrac{45}{7}=6\tfrac37,\\
x_2^*&=\tfrac{30}{7}=4\tfrac27,\\
z^*&=500\cdot\tfrac{45}{7}+450\cdot\tfrac{30}{7}\\
&=\tfrac{36000}{7}=5142\tfrac67 .
\end{aligned}
$$

As a check, the objective values at the other vertices are 0, 4000, 5080
and 3375, all smaller.

**Slack.** At the optimum the machine is used $6\cdot\tfrac{45}{7}+5\cdot
\tfrac{30}{7}=60$ h and the store holds $10\cdot\tfrac{45}{7}+20\cdot
\tfrac{30}{7}=150$ m³: both constraints are binding, with zero slack. The
demand constraint has slack $8-\tfrac{45}{7}=\tfrac{11}{7}=1\tfrac47$
hundred cases, about 157 cases of juice glasses per week that could still
be sold.

*(An earlier version of this example states the slack as $1\tfrac47$ cases
per week; because $x_1$ is measured in hundreds of cases, the slack is
$1\tfrac47$ hundred cases, or $157\tfrac17$ cases, per week.)*

**Whole hundreds of cases.** If the plan had to be in whole hundreds of
cases, rounding the optimum down gives $(6,4)$, which is feasible with
$z=4800$. It is not the best integer plan: the integer points of the region
with $z\ge4800$ are $(8,2)$ with $z=4900$, $(7,3)$ with $z=4850$ and $(6,4)$
with $z=4800$, so the integer optimum is $(8,2)$, far from the LP optimum
in the $x_1$ direction. Rounding the LP solution is not the same as solving
the integer program.

```{figure} figures/03-02-product-mix-region.svg
:name: fig-03-02-product-mix-region
:width: 540px
:class: dark-light
:alt: Pentagon-shaped feasible region bounded by three constraint lines and the axes, with three parallel dashed objective contours; the highest one touches the region only at the vertex (45/7, 30/7)

Feasible region of the product-mix LP. The dashed contours
$500x_1+450x_2=z$ for $z=2250$, $4000$ and $\tfrac{36000}{7}$ move up and
to the right as $z$ grows; the last one touches the region only at the
optimal vertex $\vx^*=(\tfrac{45}{7},\tfrac{30}{7})$.
```
:::

The solution of {prf:ref}`ex-03-02-graphical` uses the vocabulary of
slack and binding constraints, which is general.

:::{prf:definition} Slack, surplus and binding constraints
:label: def-03-02-slack

At a feasible point $\vx$, the **slack** of a constraint
$\mA_{i\cdot}\vx\le b_i$ is $s_i=b_i-\mA_{i\cdot}\vx\ge0$, the unused amount
of the resource; the **surplus** of a constraint $\mA_{i\cdot}\vx\ge b_i$
is $s_i=\mA_{i\cdot}\vx-b_i\ge0$, the amount by which the requirement is
exceeded. A constraint with zero slack or surplus at $\vx$ is **binding**
(active) at $\vx$; otherwise it is **non-binding**.
:::

A binding constraint is a scarce resource: relaxing it would let the
objective improve (in {doc}`03-05-lp-duality` this improvement per unit is
the shadow price, {ref}`sec-03-05-shadow-prices`). A non-binding
constraint could be tightened up to its slack without changing the
optimum.

(sec-03-02-graphical-cases)=
### The four possible outcomes

The product mix has a unique optimal solution. Three other outcomes can
occur, and {numref}`fig-03-02-graphical-cases` shows one small LP of each kind.

- **Unique optimum.** For $\max\ x_1+x_2$ subject to $x_1+2x_2\le6$,
  $2x_1+x_2\le6$, $\vx\ge\vzero$, the last contour touches the region only
  at $(2,2)$, with $z^*=4$.
- **Alternative optima.** With the objective $\max\ x_1+2x_2$ over the same
  region, the contours are parallel to the edge $x_1+2x_2=6$, and the last
  contour, $z^*=6$, contains the whole edge from $(0,3)$ to $(2,2)$. Every
  point $\theta(0,3)+(1-\theta)(2,2)$, $0\le\theta\le1$, is optimal. The
  optimal set is still a face of the region, and it still contains corners.
- **Unbounded LP.** For $\max\ x_1+x_2$ subject to $-x_1+x_2\le2$,
  $x_1-2x_2\le2$, $\vx\ge\vzero$, the region is unbounded and so is the
  objective: from $(0,0)$ the points $t\,(2,1)$, $t\ge0$, are all feasible
  ($-2t+t\le0\le2$ and $2t-2t=0\le2$) and their objective value $3t$ grows
  without limit. An unbounded *region* alone does not make an LP unbounded:
  $\min\ x_1+x_2$ over the same region has the optimum $(0,0)$.
- **Infeasible LP.** The constraints $x_1+x_2\le2$ and $x_1+2x_2\ge6$ with
  $\vx\ge\vzero$ have no common point, since $x_1+2x_2\le2(x_1+x_2)\le4$
  for nonnegative $\vx$.

```{figure} figures/03-02-graphical-cases.svg
:name: fig-03-02-graphical-cases
:width: 640px
:class: dark-light
:alt: Four small panels: a region with a single optimal corner; the same region with an optimal edge; an unbounded region with an arrow along which the objective grows; two half-planes that do not overlap

The four outcomes of an LP. (a) Unique optimum at $(2,2)$. (b) Alternative
optima: the contour $x_1+2x_2=6$ contains a whole edge. (c) Unbounded:
the objective $x_1+x_2$ grows along the ray from the origin in the
direction $(2,1)$. (d) Infeasible: the half-planes $x_1+x_2\le2$ and
$x_1+2x_2\ge6$ do not meet in the nonnegative quadrant.
```

In every case with an optimum, a corner of the region is optimal. This is
the observation on which the simplex method rests, and the rest of the
chapter proves it for any number of variables
({prf:ref}`thm-03-02-fundamental`). The proof needs two things: an algebraic
description of the corners, which drawing cannot supply in $\R^n$, and a
uniform way of writing every LP.

(sec-03-02-standard-form)=
## Standard form

:::{prf:definition} Standard form
:label: def-03-02-standard-form

An LP is in **standard form** if it reads

$$
\begin{aligned}
\min_{\vx}\ & z=\vc^\T\vx \\
\st\ & \mA\vx=\vb, \\
& \vx\ge\vzero,
\end{aligned}
$$ (eq-03-02-standard-form)

with $\mA\in\R^{m\times n}$, $\vb\in\R^m$ and $\vc\in\R^n$: all
constraints are equations and all variables are nonnegative. When a
starting basis is to be built ({doc}`03-03-simplex-method`), one also asks
for $\vb\ge\vzero$.
:::

The theory below is stated for minimization; a maximization with the same
constraints is handled by {eq}`eq-03-01-max-min` or by reversing the
inequalities in the optimality test, and the running examples of this
chapter are maximizations.

Every LP can be brought into standard form, which is why the theory and
the simplex method need to deal with this one form only.

:::{prf:proposition} Conversion to standard form
:label: prop-03-02-standard-form

Every LP can be rewritten as an LP in standard form by the following
operations, which change neither the set of optimal plans (after the
substitutions are undone) nor, apart from the sign in step 1, the optimal
value.

1. **Objective sense.** Replace $\max\ \vc^\T\vx$ by $\min\ (-\vc)^\T\vx$.
2. **$\le$ constraints.** Replace $\mA_{i\cdot}\vx\le b_i$ by
   $\mA_{i\cdot}\vx+s_i=b_i$ with a new **slack variable** $s_i\ge0$.
3. **$\ge$ constraints.** Replace $\mA_{i\cdot}\vx\ge b_i$ by
   $\mA_{i\cdot}\vx-s_i=b_i$ with a new **surplus variable** $s_i\ge0$.
4. **Nonpositive variables.** If $x_j\le0$, substitute $x_j=-x_j'$ with
   $x_j'\ge0$.
5. **Free variables.** If $x_j$ has no sign restriction, substitute
   $x_j=x_j^+-x_j^-$ with $x_j^+,x_j^-\ge0$.
6. **Bounds.** If $x_j\ge l_j$, substitute $x_j=l_j+x_j'$ with $x_j'\ge0$;
   an upper bound $x_j\le u_j$ is an ordinary $\le$ constraint (step 2).
7. **Negative right-hand sides.** Multiply an equation with $b_i<0$ by $-1$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Step 1 is {eq}`eq-03-01-max-min`. In step 2, the map
$\vx\mapsto(\vx,\ s_i=b_i-\mA_{i\cdot}\vx)$ is a one-to-one correspondence
between the points satisfying the inequality and the points $(\vx,s_i)$
satisfying the equation with $s_i\ge0$; the objective does not involve
$s_i$, so it takes the same value at corresponding points. Step 3 is the
same with $s_i=\mA_{i\cdot}\vx-b_i$. Steps 4 and 6 are changes of variable
that map $\{x_j\le0\}$ onto $\{x_j'\ge0\}$ and $\{x_j\ge l_j\}$ onto
$\{x_j'\ge0\}$ one to one, and the objective changes only by the constant
$c_jl_j$ in step 6. In step 5 every real $x_j$ can be written as
$x_j^+-x_j^-$ with $x_j^\pm\ge0$ (for example $x_j^+=\max\{x_j,0\}$,
$x_j^-=\max\{-x_j,0\}$), and every such pair gives back a real $x_j$; the
map is onto though not one to one, and the objective depends only on the
difference, so optimal values agree and every optimal $(x_j^+,x_j^-)$
gives an optimal $x_j$. Step 7 does not change the solution set of the
equation.
:::

Each conversion is a modeling device, not an approximation. Note two
consequences. A slack or surplus variable has the meaning of
{prf:ref}`def-03-02-slack` at every point, so it reports directly which
constraints bind. And a free variable doubles into two columns that are
negatives of each other; such columns are linearly dependent, so (by the
results below) at most one of $x_j^+$, $x_j^-$ is positive at a basic
feasible solution, and the split recovers $x_j^+=\max\{x_j,0\}$,
$x_j^-=\max\{-x_j,0\}$.

For the product mix, with slack variables named $x_3,x_4,x_5$ (unused
machine hours, unused m³ of storage, unused demand in hundreds of cases),
the standard form is

$$
\begin{aligned}
\max\ & z=500x_1+450x_2 \\
\st\ & 6x_1+5x_2+x_3=60 \\
& 10x_1+20x_2+x_4=150 \\
& x_1+x_5=8 \\
& x_1,\dots,x_5\ge0,
\end{aligned}
$$ (eq-03-02-product-mix-standard)

so that

$$
\begin{gathered}
\mA=\begin{pmatrix}6&5&1&0&0\\10&20&0&1&0\\1&0&0&0&1\end{pmatrix},\\
\vb=\begin{pmatrix}60\\150\\8\end{pmatrix},\\
\vc=(500,\,450,\,0,\,0,\,0)^\T .
\end{gathered}
$$

(The maximization is kept here, and in the simplex method of
{doc}`03-03-simplex-method` the rules for max are given alongside those for
min.) The plan with $x_1=\tfrac{45}{7}$, $x_2=\tfrac{30}{7}$ corresponds to
$\vx=(\tfrac{45}{7},\tfrac{30}{7},0,0,\tfrac{11}{7})$: the zero slacks
$x_3=x_4=0$ mark the binding constraints.

(sec-03-02-polyhedra)=
## Half-spaces and polyhedra

For a nonzero $\va\in\R^n$ and $b\in\R$, the set $\{\vx:\va^\T\vx=b\}$ is a
**hyperplane** and the sets $\{\vx:\va^\T\vx\ge b\}$ and
$\{\vx:\va^\T\vx\le b\}$ are the two closed **half-spaces** it bounds.
In the plane, a hyperplane is a line and a half-space the closed half-plane
on one side of it ({numref}`fig-03-02-half-space`). Every linear inequality
constraint of an LP describes a half-space, and every nonnegativity
condition $x_j\ge0$ is the half-space with $\va=\ve_j$, $b=0$.

```{figure} figures/03-02-half-space.svg
:name: fig-03-02-half-space
:width: 420px
:class: dark-light
:alt: The line x1 + x2 = 1 through (0,1) and (1,0) with the region above and to its right shaded and labeled x1 + x2 at least 1

The closed half-space $x_1+x_2\ge1$ (shaded) and its boundary, the
hyperplane $x_1+x_2=1$ through $(0,1)$ and $(1,0)$. The normal vector
$(1,1)$ points into the half-space.
```

What about equality constraints? An equation $\va^\T\vx=b$ holds exactly
when both $\va^\T\vx\le b$ and $\va^\T\vx\ge b$ hold, so it is the
intersection of two half-spaces. Hence the feasible set of every LP,
whatever the mix of $\le$, $\ge$ and $=$ rows, is the intersection of
finitely many closed half-spaces, that is, a **polyhedron**
({prf:ref}`def-02-01-polyhedron`); a bounded polyhedron is a **polytope**.

:::{prf:proposition} The feasible set of an LP is a convex polyhedron
:label: prop-03-02-convex-polyhedron

The feasible set of an LP is a closed convex polyhedron. In particular,
$P=\{\vx\in\R^n:\mA\vx=\vb,\ \vx\ge\vzero\}$ is convex: if $\vy,\vw\in P$
and $0\le\theta\le1$, then $\theta\vy+(1-\theta)\vw\in P$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

It is an intersection of finitely many closed half-spaces, as just shown.
Directly: $\mA(\theta\vy+(1-\theta)\vw)=\theta\vb+(1-\theta)\vb=\vb$, and
$\theta\vy+(1-\theta)\vw\ge\vzero$ because both terms are nonnegative.
Closedness holds because each half-space is the preimage of a closed
interval under a continuous linear function.
:::

Convexity is what makes LP tractable. A point of a convex set that is not
optimal can always be improved by moving toward a better point along a
segment that stays feasible, and a linear objective changes linearly along
the segment; consequently a *local* optimum of an LP is a global one
({prf:ref}`thm-02-03-local-global`). The corners of the region are
captured by the notion of an extreme point, defined in
{prf:ref}`def-02-01-extreme-point` and restated here.

:::{prf:definition} Extreme point
:label: def-03-02-extreme-point

A point $\vx$ of a convex set $C$ is an **extreme point** (corner point,
vertex) of $C$ if it cannot be written as $\vx=\theta\vy+(1-\theta)\vw$
with $\vy,\vw\in C$, $\vy\ne\vw$ and $0<\theta<1$; that is, $\vx$ does not
lie strictly inside any segment of $C$.
:::

The extreme points of the product-mix region are the five vertices listed
in the solution of {prf:ref}`ex-03-02-graphical`. A point on an edge, such as
$(8,1)=\tfrac12(8,0)+\tfrac12(8,2)$, or an interior point, is not extreme.
To find extreme points in $\R^n$ without a picture, they must be described
algebraically.

(sec-03-02-basic-solutions)=
## Basic solutions

From now on the LP is in standard form {eq}`eq-03-02-standard-form` and its
feasible set is $P=\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$. The following
assumption is made throughout the simplex chapters.

:::{prf:assumption} Full row rank
:label: asm-03-02-full-rank

$\mA\in\R^{m\times n}$ has rank $m$, so its $m$ rows are linearly
independent and $m\le n$.
:::

The assumption loses nothing. If a row of $\mA$ is a linear combination of
other rows, then either the system $\mA\vx=\vb$ is inconsistent (and the LP
is infeasible) or the same combination holds for $\vb$ and the row is a
**redundant** equation that can be deleted. (Phase I of the simplex method
detects such rows, {prf:ref}`thm-03-04-phase-one`.)

The equations $\mA\vx=\vb$ say that $\vb$ is a combination of the columns
of $\mA$ with coefficients $x_j$:

$$
x_1\mA_{\cdot1}+x_2\mA_{\cdot2}+\dots+x_n\mA_{\cdot n}=\vb .
$$ (eq-03-02-column-view)

A feasible point is a way of writing $\vb$ as a *nonnegative* combination
of columns. The simplest such representations use as few columns as
possible.

:::{prf:definition} Basis, basic solution, basic feasible solution
:label: def-03-02-basic-solution

Under {prf:ref}`asm-03-02-full-rank`:

1. A set of $m$ linearly independent columns of $\mA$ is a **basis**. The
   $m\times m$ matrix $\mB$ formed by these columns (the **basis matrix**)
   is nonsingular. The other columns form the matrix $\mN$.
2. The variables of the columns in $\mB$ are the **basic variables**
   $\vx_B$; the others are the **nonbasic variables** $\vx_N$.
3. The **basic solution** of the basis is the unique solution of
   $\mA\vx=\vb$ with $\vx_N=\vzero$, namely
   $$
   \begin{aligned}
   \vx_B&=\mB^{-1}\vb,\\
   \vx_N&=\vzero .
   \end{aligned}
   $$
4. A basic solution with $\vx_B\ge\vzero$ is a **basic feasible solution**
   (BFS).
:::

A basic solution is computed in three steps: choose $m$ columns, check that
they are independent ($\det\mB\ne0$), solve $\mB\vx_B=\vb$. There are at
most $\binom{n}{m}$ ways to choose the columns, hence at most
$\binom{n}{m}$ basic solutions.

:::{prf:example} Enumerating basic solutions of a small system
:label: ex-03-02-basic-solutions

Find all basic solutions of
$$
\begin{gathered}
2x_1+3x_2+x_3=5,\\
x_1+2x_2-x_3=3,\\
x_1,x_2,x_3\ge0,
\end{gathered}
$$

and decide which are feasible.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Here $m=2$, $n=3$ and
$\mA=\begin{pmatrix}2&3&1\\1&2&-1\end{pmatrix}$, $\vb=\begin{pmatrix}5\\3\end{pmatrix}$.
The rows are linearly independent (neither is a multiple of the other), so
{prf:ref}`asm-03-02-full-rank` holds. There are $\binom32=3$ candidate bases.

- $\{x_1,x_2\}$: $\mB=\begin{pmatrix}2&3\\1&2\end{pmatrix}$,
  $\det\mB=1$, $\mB^{-1}=\begin{pmatrix}2&-3\\-1&2\end{pmatrix}$,
  $\vx_B=\mB^{-1}\vb=(1,1)^\T$. The basic solution $(1,1,0)$ is feasible.
- $\{x_1,x_3\}$: $\mB=\begin{pmatrix}2&1\\1&-1\end{pmatrix}$,
  $\det\mB=-3$, $\mB^{-1}=\begin{pmatrix}\tfrac13&\tfrac13\\
  \tfrac13&-\tfrac23\end{pmatrix}$, $\vx_B=(\tfrac83,-\tfrac13)^\T$. The
  basic solution $(\tfrac83,0,-\tfrac13)$ is infeasible ($x_3<0$).
- $\{x_2,x_3\}$: $\mB=\begin{pmatrix}3&1\\2&-1\end{pmatrix}$,
  $\det\mB=-5$, $\mB^{-1}=\begin{pmatrix}\tfrac15&\tfrac15\\
  \tfrac25&-\tfrac35\end{pmatrix}$, $\vx_B=(\tfrac85,\tfrac15)^\T$. The
  basic solution $(0,\tfrac85,\tfrac15)$ is feasible.

Each solution satisfies both equations, for example
$2\cdot0+3\cdot\tfrac85+\tfrac15=5$ and $0+2\cdot\tfrac85-\tfrac15=3$. The
system has two basic feasible solutions, $(1,1,0)$ and $(0,\tfrac85,\tfrac15)$.
The feasible set is the segment joining them; {prf:ref}`ex-08-01-feasible-segment`
draws it.
:::

For the product mix in standard form {eq}`eq-03-02-product-mix-standard`
($m=3$, $n=5$) there are $\binom53=10$ candidate bases.
{numref}`tab-03-02-product-mix-basic` lists them, and
{numref}`fig-03-02-product-mix-basic-solutions` locates each basic solution in
the $(x_1,x_2)$ plane: every basic solution is an intersection of two of
the five lines $6x_1+5x_2=60$, $10x_1+20x_2=150$, $x_1=8$, $x_1=0$, $x_2=0$,
because the two nonbasic variables are set to zero and each variable being
zero means "on that line".

```{list-table} The ten candidate bases of the product mix in standard form. The nonbasic variables are zero; a dash marks a singular basis.
:name: tab-03-02-product-mix-basic
:header-rows: 1

* - Basic variables
  - $(x_1,x_2,x_3,x_4,x_5)$
  - Feasible?
  - $z$
* - $x_1,x_2,x_3$
  - $(8,\tfrac72,-\tfrac{11}{2},0,0)$
  - no
  - 5575
* - $x_1,x_2,x_4$
  - $(8,\tfrac{12}{5},0,22,0)$
  - yes
  - 5080
* - $x_1,x_2,x_5$
  - $(\tfrac{45}{7},\tfrac{30}{7},0,0,\tfrac{11}{7})$
  - yes
  - $\tfrac{36000}{7}$
* - $x_1,x_3,x_4$
  - $(8,0,12,70,0)$
  - yes
  - 4000
* - $x_1,x_3,x_5$
  - $(15,0,-30,0,-7)$
  - no
  - 7500
* - $x_1,x_4,x_5$
  - $(10,0,0,50,-2)$
  - no
  - 5000
* - $x_2,x_3,x_4$
  - —
  - singular
  -
* - $x_2,x_3,x_5$
  - $(0,\tfrac{15}{2},\tfrac{45}{2},0,8)$
  - yes
  - 3375
* - $x_2,x_4,x_5$
  - $(0,12,0,-90,8)$
  - no
  - 5400
* - $x_3,x_4,x_5$
  - $(0,0,60,150,8)$
  - yes
  - 0
```

```{figure} figures/03-02-product-mix-basic-solutions.svg
:name: fig-03-02-product-mix-basic-solutions
:width: 540px
:class: dark-light
:alt: The product-mix region with five filled dots at its corners, the basic feasible solutions, and four hollow dots outside it at (8, 3.5), (10, 0), (15, 0) and (0, 12), the infeasible basic solutions

The basic solutions of the product mix, drawn in the $(x_1,x_2)$ plane.
Filled dots are the five basic feasible solutions, which are exactly the
corners of the feasible region; hollow dots are the four infeasible basic
solutions, intersections of constraint lines that violate some other
constraint, such as $(8,\tfrac72)$, which needs 65.5 machine hours.
```

The table confirms what the figure suggests. The five basic feasible
solutions are the five corners of the region, and the largest $z$ among
them, $\tfrac{36000}{7}$, is the optimum found graphically. The infeasible
basic solutions can have larger $z$ (7500 at $(15,0)$), which is why
feasibility must be checked. The basis $\{x_2,x_3,x_4\}$ is singular: its
columns have zeros in the third row, and indeed the lines $x_1=0$ and
$x_1=8$ are parallel and never meet.

A basic solution can have basic variables equal to zero. This seemingly
minor case matters for the simplex method.

:::{prf:definition} Degenerate basic solution
:label: def-03-02-degenerate

A basic solution is **degenerate** if at least one basic variable is zero,
that is, if it has fewer than $m$ positive components. Otherwise it is
**nondegenerate**.
:::

A degenerate BFS belongs to more than one basis: the zero basic variable
can be exchanged for a suitable nonbasic variable without changing the
point. Geometrically, in $\R^2$, more than two of the lines that bound the
region pass through the corner. In {doc}`03-03-simplex-method`, degeneracy
causes pivots that change the basis but not the point, and in rare cases
cycling ({ref}`sec-03-03-degeneracy`). {prf:ref}`ex-08-01-bfs-test` contains a degenerate BFS.

(sec-03-02-independence-test)=
### Testing linear independence by pivoting

Deciding whether a set of columns is a basis, or whether the columns used
by a given solution are independent, requires a test of linear
independence. Besides computing a determinant (only for square matrices)
or a rank, a simple pivoting test also returns a dependence relation when
there is one.

:::{prf:algorithm} Linear-independence test by pivoting
:label: alg-03-02-independence

**Input:** vectors $\vv_1,\dots,\vv_k\in\R^m$.

1. **Set up** the tableau $[\,\mI_k\mid\mV\,]$, where row $r$ of $\mV$ is
   $\vv_r^\T$ (the left part is the $k\times k$ identity).
2. **Examine** rows $r=1,\dots,k$ in turn.
   - If the right part of row $r$ is zero, **stop**: the vectors are
     linearly dependent, and the left part $(\gamma_1,\dots,\gamma_k)$ of
     row $r$ gives the dependence $\sum_i\gamma_i\vv_i=\vzero$.
   - Otherwise choose a nonzero entry of the right part of row $r$ as the
     pivot, divide row $r$ by it, and subtract multiples of row $r$ from
     all other rows to make the rest of the pivot column zero.
3. **Conclude.** If all $k$ rows have been pivoted, the vectors are
   linearly independent.
:::

Why does a zero right part mean dependence? Every row of the tableau is
a combination of the initial rows, which are $[\ve_r^\T\mid\vv_r^\T]$.
A row equal to $\sum_i\gamma_i[\ve_i^\T\mid\vv_i^\T]$ has left part
$(\gamma_1,\dots,\gamma_k)$ and right part $\sum_i\gamma_i\vv_i^\T$, and
row operations preserve this structure. So a row whose right part is zero
states $\sum_i\gamma_i\vv_i=\vzero$, and its left part is not zero
because the left block stays nonsingular under row operations. Hence the
vectors are dependent, with exactly these coefficients. Conversely, if $k$
pivots succeed, the right part has been reduced to $k$ independent rows,
and since the left block is nonsingular the original vectors span the same
$k$-dimensional space. {prf:ref}`ex-08-01-independence-pivoting` applies the
test to four vectors in $\R^5$ and finds
$-\vv_1-\vv_2-\vv_3+\vv_4=\vzero$.

(sec-03-02-extreme-bfs)=
## Extreme points are basic feasible solutions

The table and figure of the product mix showed that its corners and its
basic feasible solutions coincide. This is true for every LP.

:::{prf:theorem} Extreme points and basic feasible solutions
:label: thm-03-02-extreme-bfs

Let $P=\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ with $\mA$ of full row rank. A
point $\vx\in P$ is an extreme point of $P$ if and only if it is a basic
feasible solution. Moreover, $\vx\in P$ is a BFS if and only if the columns
$\{\mA_{\cdot j}:j\in J(\vx)\}$ of its positive components are linearly
independent.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

*BFS $\Rightarrow$ extreme point.* Let $\vx$ be the BFS of a basis $\mB$,
and suppose $\vx=\theta\vy+(1-\theta)\vw$ with $\vy,\vw\in P$ and
$0<\theta<1$. For every nonbasic index $j$, $0=x_j=\theta y_j+(1-\theta)w_j$
with $y_j,w_j\ge0$ forces $y_j=w_j=0$. Then $\mB\vy_B=\mA\vy=\vb$, so
$\vy_B=\mB^{-1}\vb=\vx_B$; likewise $\vw_B=\vx_B$. Hence $\vy=\vw=\vx$, and
$\vx$ is extreme.

*Positive columns independent $\Rightarrow$ BFS.* Let $J=J(\vx)$ and suppose
the columns $\mA_{\cdot j}$, $j\in J$, are independent. Then $\abs{J}\le m$,
and since $\mA$ has rank $m$ these columns can be completed by $m-\abs{J}$
further columns of $\mA$ to a basis $\mB$. The point $\vx$ satisfies
$\mA\vx=\vb$ and is zero outside $J$, hence outside the basis; by the
uniqueness in item 3 of {prf:ref}`def-03-02-basic-solution` it is the basic solution of
$\mB$, and it is feasible.

*Extreme point $\Rightarrow$ positive columns independent.* Suppose the
columns $\mA_{\cdot j}$, $j\in J$, are dependent: there is $\vd\ne\vzero$
with $d_j=0$ for $j\notin J$ and $\mA\vd=\sum_{j\in J}d_j\mA_{\cdot j}=\vzero$.
Because $x_j>0$ for $j\in J$, the points $\vx+\varepsilon\vd$ and
$\vx-\varepsilon\vd$ are nonnegative for small $\varepsilon>0$, and both
satisfy $\mA(\vx\pm\varepsilon\vd)=\vb$. So both lie in $P$, they are
distinct, and $\vx=\tfrac12(\vx+\varepsilon\vd)+\tfrac12(\vx-\varepsilon\vd)$
is not extreme.

The three implications prove both statements: a BFS has independent
positive columns (they are among the columns of $\mB$), and a point with
independent positive columns is a BFS.
:::

The theorem turns a geometric property into an algebraic one that can be
checked: a feasible point is a corner exactly when the columns it uses
are independent. Two consequences follow at once.

:::{prf:corollary} Finitely many extreme points
:label: cor-03-02-finite

The polyhedron $P$ has at most $\binom{n}{m}$ extreme points. A
nondegenerate BFS corresponds to exactly one basis; a degenerate BFS may
correspond to several.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Each extreme point is the basic solution of at least one of the at most
$\binom nm$ bases, and a basis has only one basic solution. A
nondegenerate BFS has $m$ positive components, whose columns must form its
basis; for a degenerate BFS the columns of the positive components can be
completed to a basis in more than one way.
:::

(sec-03-02-fundamental-theorem)=
## The fundamental theorem of linear programming

Graphically, an LP with an optimum has an optimal corner. The algebraic
statement is the following theorem, whose proof reduces any feasible or
optimal point to a basic one by moving along directions in which the
constraints stay satisfied.

:::{prf:theorem} Fundamental theorem of linear programming
:label: thm-03-02-fundamental

Consider the LP {eq}`eq-03-02-standard-form` with $\mA$ of full row rank.

1. If the LP has a feasible solution, it has a basic feasible solution.
2. If the LP has an optimal solution, it has an optimal basic feasible
   solution.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

We prove part 2; part 1 is the same argument without the objective (take
$\vc=\vzero$, for which every feasible point is optimal).

Among all optimal solutions choose one, $\vx$, with the smallest number of
positive components, and let $J=J(\vx)$. If the columns $\mA_{\cdot j}$,
$j\in J$, are linearly independent, $\vx$ is a BFS by
{prf:ref}`thm-03-02-extreme-bfs`, and we are done. Suppose they are dependent:
there is $\vd\ne\vzero$ with $d_j=0$ for $j\notin J$ and $\mA\vd=\vzero$.
Replacing $\vd$ by $-\vd$ if necessary, we may assume that $\vd$ has at
least one negative component.

For every $\varepsilon$ with $\abs{\varepsilon}$ small enough, $\vx+\varepsilon\vd$
is feasible, as in the proof of {prf:ref}`thm-03-02-extreme-bfs`. Its
objective value is $\vc^\T\vx+\varepsilon\,\vc^\T\vd$. If $\vc^\T\vd\ne0$,
a small $\varepsilon$ of the sign opposite to $\vc^\T\vd$ would give a
feasible point with a smaller objective value than the optimal $\vx$,
which is impossible. Hence $\vc^\T\vd=0$, and $\vx+\varepsilon\vd$ is
optimal whenever it is feasible.

Now increase $\varepsilon$ from 0 as far as nonnegativity allows:

$$
\varepsilon^*=\min\left\{\frac{x_j}{-d_j}:\ d_j<0\right\}>0 .
$$

The point $\vx'=\vx+\varepsilon^*\vd$ satisfies $\mA\vx'=\vb$, is
nonnegative by the choice of $\varepsilon^*$, and is optimal since
$\vc^\T\vx'=\vc^\T\vx$. Its components outside $J$ are still zero, and the
component $j$ that attains the minimum has become zero. So $\vx'$ is an
optimal solution with fewer positive components than $\vx$, contradicting
the choice of $\vx$. Hence the columns of $J$ are independent and $\vx$ is
an optimal BFS.
:::

Together with {prf:ref}`thm-03-02-extreme-bfs`, the theorem says: *if an LP
has an optimal solution, some extreme point of its feasible set is
optimal.* It does not say that every optimal solution is an extreme point:
with alternative optima ({numref}`fig-03-02-graphical-cases`(b)) the points
inside an optimal edge are optimal too, but the edge ends at optimal
corners.

The theorem also explains when an optimum exists. A nonempty polyhedron
in standard form has a BFS by part 1, and {doc}`03-03-simplex-method` shows
that, starting from a BFS, the simplex method with an anticycling rule
either finds an optimal BFS or a ray along which the objective improves
without limit ({prf:ref}`thm-03-03-bland`). So every LP
ends in exactly one of three ways: it is **infeasible**, it is
**unbounded**, or it has an **optimal BFS**. In particular, an LP that is
feasible and whose objective is bounded on the feasible set always attains
its optimum, unlike a general optimization problem.

**Why not simply enumerate?** The theorem gives a finite algorithm:
compute all basic solutions, discard the infeasible ones, and keep the best.
For the product mix that means ten small systems
({numref}`tab-03-02-product-mix-basic`), and {doc}`08-01-examples-lp-modeling`
solves several LPs this way ({prf:ref}`ex-08-01-basic-solutions-2d`,
{prf:ref}`ex-08-01-enumerate-bases`).
But the number of candidate bases grows explosively: an LP with $m=25$
equations and $n=50$ variables has $\binom{50}{25}\approx1.26\times10^{14}$
of them. The simplex method instead moves from one BFS to a neighboring,
better one, and it needs a test that tells when a BFS is optimal without
comparing it with all the others. The reduced form supplies that test.

(sec-03-02-reduced-form)=
## Bases and the reduced form

Fix a basis. After reordering the columns, write $\mA=[\,\mB\ \mN\,]$,
$\vx=(\vx_B,\vx_N)$ and $\vc=(\vc_B,\vc_N)$. The constraints and the
objective become

$$
\begin{aligned}
\mB\vx_B+\mN\vx_N&=\vb,\\
z&=\vc_B^\T\vx_B+\vc_N^\T\vx_N .
\end{aligned}
$$

Since $\mB$ is invertible, the constraints can be solved for the basic
variables, $\vx_B=\mB^{-1}\vb-\mB^{-1}\mN\vx_N$, and substituting into $z$
gives

$$
z=\vc_B^\T\mB^{-1}\vb+\bigl(\vc_N^\T-\vc_B^\T\mB^{-1}\mN\bigr)\vx_N .
$$

Every solution of $\mA\vx=\vb$ is therefore determined by its nonbasic
part $\vx_N$, and the LP is equivalent to the **reduced form** in the
nonbasic variables alone:

$$
\begin{aligned}
\min_{\vx_N}\ & z=\vc_B^\T\mB^{-1}\vb+\bar\vc_N^\T\vx_N \\
\st\ & \mB^{-1}\vb-\mB^{-1}\mN\vx_N\ge\vzero, \\
& \vx_N\ge\vzero,
\end{aligned}
$$ (eq-03-02-reduced-form)

where the first constraint says $\vx_B\ge\vzero$ and

$$
\begin{gathered}
\bar\vc_N^\T=\vc_N^\T-\vc_B^\T\mB^{-1}\mN,\\
\text{that is,}\quad \bar c_j=c_j-\vpi^\T\mA_{\cdot j},\\
\text{with}\quad\vpi^\T=\vc_B^\T\mB^{-1}.
\end{gathered}
$$ (eq-03-02-reduced-costs)

The number $\bar c_j$ is the **reduced cost** (relative cost) of $x_j$: the
rate at which $z$ changes when the nonbasic variable $x_j$ is increased
from zero and the basic variables adjust to keep $\mA\vx=\vb$. For a basic
variable the same formula gives $\bar c_j=c_j-\vc_B^\T\mB^{-1}\mB\ve_j=0$.
The vector $\vpi$ will reappear in {doc}`03-05-lp-duality` as the vector of
simplex multipliers or shadow prices ({prf:ref}`thm-03-05-shadow-price`). At the basic solution
($\vx_N=\vzero$) the objective value is $z=\vc_B^\T\mB^{-1}\vb=\vc_B^\T\vx_B$.

:::{prf:theorem} Optimality test for a basic feasible solution
:label: thm-03-02-optimality-test

Let $\vx$ be the BFS of a basis $\mB$, with reduced costs
{eq}`eq-03-02-reduced-costs`.

1. For a **minimization**, if $\bar c_j\ge0$ for all nonbasic $j$, then
   $\vx$ is optimal.
2. For a **maximization**, if $\bar c_j\le0$ for all nonbasic $j$, then
   $\vx$ is optimal.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\vy$ be any feasible solution. It satisfies $\mA\vy=\vb$, so by the
derivation of the reduced form its objective value is
$\vc^\T\vy=\vc_B^\T\mB^{-1}\vb+\bar\vc_N^\T\vy_N$, while
$\vc^\T\vx=\vc_B^\T\mB^{-1}\vb$ because $\vx_N=\vzero$. Hence

$$
\vc^\T\vy-\vc^\T\vx=\sum_{j\ \text{nonbasic}}\bar c_j\,y_j .
$$

Every $y_j\ge0$. In case 1 every $\bar c_j\ge0$, so the sum is
nonnegative and $\vc^\T\vy\ge\vc^\T\vx$; in case 2 every $\bar c_j\le0$, so
$\vc^\T\vy\le\vc^\T\vx$. Since $\vy$ was arbitrary, $\vx$ is optimal.
:::

The condition is **sufficient**. It is also necessary for a nondegenerate
BFS: if, in a minimization, some $\bar c_j<0$ and all basic variables are
positive, then increasing $x_j$ slightly keeps $\vx_B>\vzero$ and lowers
$z$, so $\vx$ is not optimal. For a degenerate BFS it is not necessary:
in $\min\{-x_1:\ x_1+x_2=0,\ x_1,x_2\ge0\}$ the only feasible point
$(0,0)$ is optimal, but with the basis $\{x_2\}$ the reduced cost of $x_1$ is
$\bar c_1=-1<0$. Some other basis of the same degenerate point passes the
test (here $\{x_1\}$, with $\bar c_2=1$), and the simplex method finds it.

:::{prf:example} Reduced form and optimality of a basis
:label: ex-03-02-reduced-form

Consider

$$
\begin{aligned}
\max\ & z=3x_1+x_2+3x_3 \\
\st\ & 2x_1+x_2+x_3+x_4=5 \\
& x_1+2x_2+3x_3+x_5=5 \\
& 2x_1+2x_2+x_3+x_6=6 \\
& x_1,\dots,x_6\ge0,
\end{aligned}
$$

where $x_4,x_5,x_6$ are slack variables. Take the basis of columns
$\mA_{\cdot1},\mA_{\cdot3},\mA_{\cdot6}$. Write the LP in reduced form,
find the BFS and decide whether it is optimal.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Partition.** With $\vx_B=(x_1,x_3,x_6)$ and $\vx_N=(x_2,x_4,x_5)$,

$$
\begin{gathered}
\mB=\begin{pmatrix}2&1&0\\1&3&0\\2&1&1\end{pmatrix},\\
\mN=\begin{pmatrix}1&1&0\\2&0&1\\2&0&0\end{pmatrix},\\
\vc_B=\begin{pmatrix}3\\3\\0\end{pmatrix},\\
\vc_N=\begin{pmatrix}1\\0\\0\end{pmatrix}.
\end{gathered}
$$

$\det\mB=5\ne0$, so the columns form a basis, and

$$
\begin{gathered}
\mB^{-1}=\begin{pmatrix}\tfrac35&-\tfrac15&0\\-\tfrac15&\tfrac25&0\\
-1&0&1\end{pmatrix},\\
\mB^{-1}\vb=\begin{pmatrix}2\\1\\1\end{pmatrix},\\
\mB^{-1}\mN=\begin{pmatrix}\tfrac15&\tfrac35&-\tfrac15\\
\tfrac35&-\tfrac15&\tfrac25\\1&-1&0\end{pmatrix}.
\end{gathered}
$$

**BFS.** $\vx_B=\mB^{-1}\vb=(2,1,1)\ge\vzero$, so the basic solution
$\vx=(2,0,1,0,0,1)$ is feasible, with $z=3\cdot2+3\cdot1+0\cdot1=9$. It
satisfies the three equations: $4+1+0=5$, $2+3+0=5$, $4+1+1=6$.

**Reduced costs.** $\vpi^\T=\vc_B^\T\mB^{-1}=(\tfrac65,\tfrac35,0)$, so

$$
\begin{gathered}
\bar c_2=1-\vpi^\T\mA_{\cdot2}=1-\bigl(\tfrac65+\tfrac65+0\bigr)=-\tfrac75,\\
\bar c_4=0-\tfrac65=-\tfrac65,\\
\bar c_5=0-\tfrac35=-\tfrac35 .
\end{gathered}
$$

**Reduced form.**

$$
\begin{aligned}
\max\ & z=9-\tfrac75x_2-\tfrac65x_4-\tfrac35x_5 \\
\st\ & x_1=2-\tfrac15x_2-\tfrac35x_4+\tfrac15x_5\ge0 \\
& x_3=1-\tfrac35x_2+\tfrac15x_4-\tfrac25x_5\ge0 \\
& x_6=1-x_2+x_4\ge0 \\
& x_2,x_4,x_5\ge0 .
\end{aligned}
$$

**Optimality.** All three reduced costs are negative, so by
{prf:ref}`thm-03-02-optimality-test` (maximization) the BFS is optimal: any
feasible point has $x_2,x_4,x_5\ge0$, and each unit of any of them lowers
$z$ below 9. The optimum is $x_1^*=2$, $x_3^*=1$, the first two
constraints binding and the third with slack $x_6^*=1$, and $z^*=9$, which
agrees with a solver. Of the $\binom63=20$ candidate bases of this
system, 19 are nonsingular and 8 of their basic solutions are feasible; the
reduced-cost test identified the optimal one without looking at the others.

*(An earlier version of this example regroups the first and third
constraints with $x_2+x_3$ and $2x_2+x_5$ as the nonbasic parts; with the
basis $\{x_1,x_3,x_6\}$ the nonbasic parts are $x_2+x_4$ and $2x_2$, and
the index range of the nonnegativity condition is $1,\dots,6$.)*
:::

When the test fails, the reduced form shows what to do. In a maximization,
a nonbasic variable with $\bar c_j>0$ improves $z$ at the rate $\bar c_j$
when it is increased from zero; it can be increased until one of the basic
variables, which change according to the constraint rows of
{eq}`eq-03-02-reduced-form`, falls to zero. That basic variable leaves the
basis, $x_j$ enters, and a new BFS with a better (or, when degenerate,
equal) objective value is reached. Repeating this step is the **simplex
method** of {doc}`03-03-simplex-method`; geometrically it walks from corner to
adjacent corner of $P$ along edges that improve the objective.

## Summary and Common Mistakes

### Summary

- A two-variable LP is solved graphically by shading the feasible region
  and sliding an objective contour ({prf:ref}`alg-03-02-graphical`); the
  product mix has its optimum at the corner $(\tfrac{45}{7},\tfrac{30}{7})$
  with $z^*=\tfrac{36000}{7}$ ({prf:ref}`ex-03-02-graphical`,
  {numref}`fig-03-02-product-mix-region`). Slack, surplus and binding
  constraints describe how resources are used ({prf:ref}`def-03-02-slack`).
- An LP has a unique optimum, alternative optima, is unbounded or is
  infeasible ({numref}`fig-03-02-graphical-cases`).
- Every LP can be converted to standard form $\min\{\vc^\T\vx:\mA\vx=\vb,\
  \vx\ge\vzero\}$ ({prf:ref}`prop-03-02-standard-form`); its feasible set is a
  convex polyhedron ({prf:ref}`prop-03-02-convex-polyhedron`).
- A basis is a set of $m$ independent columns; its basic solution is
  $\vx_B=\mB^{-1}\vb$, $\vx_N=\vzero$, and it is a BFS if $\vx_B\ge\vzero$
  ({prf:ref}`def-03-02-basic-solution`); a BFS with a zero basic variable is
  degenerate ({prf:ref}`def-03-02-degenerate`).
- Extreme points are exactly the BFS, and there are at most $\binom nm$ of
  them ({prf:ref}`thm-03-02-extreme-bfs`, {prf:ref}`cor-03-02-finite`).
- If an LP has an optimal solution, it has an optimal BFS
  ({prf:ref}`thm-03-02-fundamental`).
- In reduced form {eq}`eq-03-02-reduced-form` the objective is
  $\vc_B^\T\mB^{-1}\vb+\bar\vc_N^\T\vx_N$ with reduced costs
  $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$; a BFS with $\bar\vc_N\ge\vzero$ (min)
  or $\bar\vc_N\le\vzero$ (max) is optimal ({prf:ref}`thm-03-02-optimality-test`,
  {prf:ref}`ex-03-02-reduced-form`).

### Common mistakes

- Reading the optimum off a graph without computing the intersection of
  the two binding lines exactly; graphs give the corner, not its
  coordinates.
- Concluding that an LP is unbounded because its feasible region is
  unbounded; the objective must also improve along a ray.
- Forgetting to check $\vx_B\ge\vzero$: a basic solution can lie outside
  the feasible region and have a better objective value than the optimum.
- Choosing $m$ columns without checking that they are independent; a
  singular set of columns is not a basis.
- Treating the reduced-cost condition as necessary at a degenerate BFS, or
  using the min rule ($\bar c_j\ge0$) for a max problem.
- Reporting slacks in the wrong units when the variables are scaled
  (hundreds of cases, thousands of barrels).

## References

- Bazaraa, M. S., Jarvis, J. J., & Sherali, H. D. (2010). *Linear
  programming and network flows* (4th ed.). Wiley. (Ch. 2, §§2.4–2.6:
  convex sets, polyhedral sets, extreme points; ch. 3, §§3.1–3.4: basic
  feasible solutions, key to the simplex method.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Ch. 1, §1.4: graphical
  representation; ch. 2, §§2.1–2.6: polyhedra, extreme points, vertices and
  basic feasible solutions, standard form, degeneracy, existence and
  optimality of extreme points.)
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Ch. 7: the geometry of linear programs.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Ch. 3, §3.1: graphical solution;
  ch. 5, §5.1: foundations of the simplex method.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Ch. 2, §2.2: graphical LP solution; ch. 3, §§3.1–3.2: LP in
  equation form, from graphical to algebraic solution; ch. 7: from extreme
  points to basic solutions.)
