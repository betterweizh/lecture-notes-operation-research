# The Simplex Method

## Overview

The fundamental theorem of linear programming ({doc}`03-02-lp-geometry`) says
that a linear program (LP) with an optimal solution has an optimal **basic feasible
solution** (BFS), that is, an optimal extreme point. The simplex method,
devised by G. B. Dantzig in 1947, turns this into an algorithm: start at a
BFS, and as long as some neighboring BFS is better, move to it. This chapter
develops the method from the algebra of bases. It first writes an LP in
**canonical form** with respect to a basis, which shows at a glance the
current BFS, its objective value and the **reduced costs** of the nonbasic
variables. The reduced costs give an optimality test; when the test fails,
one nonbasic variable is increased, the **minimum-ratio test** decides how
far, and one basic variable leaves the basis. The same step is then carried
out on a **simplex tableau** by one Gauss–Jordan **pivot**. The whole method
is run on the product-mix LP and drawn on its feasible region. The chapter
ends with what can go wrong or be special at the end: unboundedness,
alternative optima, degeneracy and **cycling**, and **Bland's rule**, which
prevents cycling. Starting bases for LPs without an obvious one are the
subject of {doc}`03-04-two-phase-big-m`. Inverse matrices and
Gauss–Jordan elimination from {doc}`01-01-linear-algebra` are used
throughout.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- write an LP in canonical form with respect to a basis, in equations, in
  matrix form and as a tableau, and read off the BFS and its objective value;
- compute reduced costs and simplex multipliers and apply the optimality
  test for minimization and maximization problems;
- derive the minimum-ratio test from the requirement that the next solution
  stays feasible, and explain why the result is again a BFS;
- carry out pivots on a tableau with exact fractions and follow the path of
  the method on a two-variable feasible region;
- prove that the method terminates when no BFS is degenerate, and recognize
  unboundedness, alternative optima and degeneracy in a tableau; and
- explain how cycling can occur and apply Bland's rule to prevent it.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\mathcal B$ | index set of the basic variables, $\abs{\mathcal B}=m$ |
| $B(i)$ | index of the basic variable of row $i$ of a tableau, so $\vx_B=(x_{B(1)},\dots,x_{B(m)})^\T$ |
| $\bar{\va}_j=\mB^{-1}\mA_{\cdot j}$ | updated column of $x_j$ with entries $\bar a_{ij}$ |
| $s$, $r$ | index of the entering variable, row of the leaving variable (pivot row) |
| $\bar a_{rs}$ | pivot element |
| $\vd$ | edge direction of a simplex step, $\vd_B=-\bar{\va}_s$, $d_s=1$ |
| $\theta$ | step length of a pivot, the minimum ratio |

The other symbols ($\mA$, $\vb$, $\vc$, $\mB$, $\mN$, $\vx_B$, $\vx_N$,
$\bar{\vb}$, $\bar c_j$, $\vpi$, $z$) are those of the book's notation page.

(sec-03-03-canonical-form)=
## Bases, basic feasible solutions and canonical form

### Standard form and bases

The simplex method works on an LP in **standard form**

$$
\begin{aligned}
\min_{\vx}\ & z=\vc^\T\vx\\
\st\ & \mA\vx=\vb,\\
& \vx\ge\vzero,
\end{aligned}
$$ (eq-03-03-standard-form)

with $\mA\in\R^{m\times n}$ of full row rank $m$ and $\vb\ge\vzero$.
{doc}`03-02-lp-geometry` shows how any LP is brought to this form: inequality
constraints receive slack or surplus variables, free variables are split into
two nonnegative parts, lower bounds are shifted to zero, a constraint with a
negative right-hand side is multiplied by $-1$, and $\max\vc^\T\vx$ is the
same as $-\min(-\vc)^\T\vx$. Worked transformations are in
{prf:ref}`ex-08-02-slack-surplus` and {prf:ref}`ex-08-02-bounds-free`. Maximization
problems are kept as they are in this book; the rules for them are stated
next to the rules for minimization.

A **basis** is a set $\mathcal B$ of $m$ column indices whose columns form a
nonsingular $m\times m$ matrix $\mB=[\mA_{\cdot B(1)},\dots,\mA_{\cdot B(m)}]$.
The remaining columns form the nonbasic matrix $\mN$, and the variables are
split accordingly into $\vx_B$ and $\vx_N$, the costs into $\vc_B$ and
$\vc_N$. Setting $\vx_N=\vzero$ and solving $\mB\vx_B=\vb$ gives the
**basic solution** $\vx_B=\bar{\vb}=\mB^{-1}\vb$, $\vx_N=\vzero$. It is a
**basic feasible solution** when $\bar{\vb}\ge\vzero$, and it is
**degenerate** when some basic variable is zero. BFSs are exactly the extreme
points of the feasible polyhedron ({prf:ref}`thm-03-02-extreme-bfs`), and two bases that differ in one index
are called **adjacent**; their BFSs are equal or joined by an edge of the
polyhedron ({doc}`03-02-lp-geometry`).

### Canonical form

Partition $\mA\vx=\vb$ as $\mB\vx_B+\mN\vx_N=\vb$ and multiply by $\mB^{-1}$:

$$
\vx_B+\mB^{-1}\mN\vx_N=\mB^{-1}\vb=\bar{\vb}.
$$ (eq-03-03-canonical-constraints)

Every basic variable now appears in exactly one row, with coefficient $+1$.
Substituting $\vx_B=\bar{\vb}-\mB^{-1}\mN\vx_N$ into the objective gives

$$
\begin{aligned}
z&=\vc_B^\T\vx_B+\vc_N^\T\vx_N\\
&=\vc_B^\T\mB^{-1}\vb+\left(\vc_N^\T-\vc_B^\T\mB^{-1}\mN\right)\vx_N\\
&=\vc_B^\T\bar{\vb}+\sum_{j\notin\mathcal B}\bar c_j x_j ,
\end{aligned}
$$ (eq-03-03-objective-reduced)

where the **simplex multipliers** and the **reduced costs** (relative cost
coefficients) are

$$
\begin{aligned}
\vpi^\T&=\vc_B^\T\mB^{-1},\\
\bar c_j&=c_j-\vpi^\T\mA_{\cdot j}, \qquad j=1,\dots,n.
\end{aligned}
$$ (eq-03-03-reduced-cost)

For the basic variable of row $i$, $j=B(i)$, the column is
$\mA_{\cdot j}=\mB\ve_i$, so $\bar c_j=c_j-\vc_B^\T\mB^{-1}\mB\ve_i=c_j-c_{B(i)}=0$. The
multipliers $\vpi$ are the dual variables of {doc}`03-05-lp-duality`
({ref}`sec-03-05-reading-duals`); here they are only a convenient way to compute
all $\bar c_j$ at once. Part 1 of the next theorem is
{prf:ref}`thm-03-02-optimality-test`, restated for the canonical form.

:::{prf:definition} Canonical form
:label: def-03-03-canonical-form

An LP in standard form is in **canonical form** with respect to a basis
$\mathcal B$ when it is written as

$$
\begin{aligned}
\vx_B+\mB^{-1}\mN\vx_N&=\bar{\vb},\\
-z+\bar{\vc}_N^\T\vx_N&=-\vc_B^\T\bar{\vb},
\end{aligned}
$$

with $\bar{\vb}=\mB^{-1}\vb\ge\vzero$: every constraint contains one basic
variable with coefficient $+1$ that appears in no other constraint and not
in the objective row, and the right-hand sides are nonnegative.
:::

The objective has been written as one more equation, $-z+\vc^\T\vx=0$, in
which $-z$ plays the part of a basic variable that never leaves. In canonical
form the current BFS can be read off without any computation: nonbasic
variables are zero, the basic variable of row $i$ equals $\bar b_i$, and the
right-hand side of the objective row equals $-z$.

:::{prf:example} The product-mix LP in canonical form
:label: ex-03-03-product-mix

The product-mix LP of {prf:ref}`ex-03-01-product-mix` plans the weekly output of
a molding shop: $x_1$ and $x_2$ are the hundreds of cases of juice glasses
and of cocktail glasses produced per week, which contribute \$500 and \$450
per hundred cases. One hundred cases take 6 h and 5 h of machine time
(60 h available), and 10 m³ and 20 m³ of storage space (150 m³ available),
and at most 8 hundred cases of juice glasses can be sold:

$$
\begin{aligned}
\max_{\vx}\ & z=500x_1+450x_2\\
\st\ & 6x_1+5x_2\le 60,\\
& 10x_1+20x_2\le 150,\\
& x_1\le 8,\\
& x_1,x_2\ge 0.
\end{aligned}
$$

Write this LP in standard form and in canonical form with respect to the
slack basis $\mathcal B=\{3,4,5\}$ and with respect to the basis
$\mathcal B=\{1,3,4\}$. For each basis give the BFS, its objective value
and the reduced costs of the nonbasic variables.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

With slack variables $x_3$ (idle machine hours), $x_4$ (unused storage, m³)
and $x_5$ (unmet juice-glass demand, hundreds of cases) the standard form
{eq}`eq-03-02-product-mix-standard` has

$$
\begin{gathered}
\mA=\begin{pmatrix}6&5&1&0&0\\10&20&0&1&0\\1&0&0&0&1\end{pmatrix},\\
\vb=\begin{pmatrix}60\\150\\8\end{pmatrix},\\
\vc=(500,450,0,0,0)^\T .
\end{gathered}
$$

The slack columns form an identity matrix, so the system is already in
canonical form with respect to $\mathcal B=\{3,4,5\}$:

$$
\begin{aligned}
6x_1+5x_2+x_3&=60,\\
10x_1+20x_2+x_4&=150,\\
x_1+x_5&=8,\\
-z+500x_1+450x_2&=0 .
\end{aligned}
$$

Its BFS is $\vx=(0,0,60,150,8)^\T$ with $z=0$: produce nothing, and all
resources are idle. Here $\vc_B=\vzero$, so $\vpi=\vzero$ and the reduced
costs are simply $\bar c_1=500$ and $\bar c_2=450$.

A second basis, used in {ref}`sec-03-03-matrix-step`, is
$\mathcal B=\{1,3,4\}$ with $\vx_B=(x_1,x_3,x_4)^\T$:

$$
\begin{gathered}
\mB=\begin{pmatrix}6&1&0\\10&0&1\\1&0&0\end{pmatrix},\\
\mB^{-1}=\begin{pmatrix}0&0&1\\1&0&-6\\0&1&-10\end{pmatrix},\\
\bar{\vb}=\mB^{-1}\vb=\begin{pmatrix}8\\12\\70\end{pmatrix}\ge\vzero ,
\end{gathered}
$$

so this basis is feasible. With $\vc_B=(500,0,0)^\T$ the multipliers are
$\vpi^\T=\vc_B^\T\mB^{-1}=(0,0,500)$ and, for the nonbasic columns
$\mN=[\mA_{\cdot2},\mA_{\cdot5}]$,

$$
\begin{gathered}
\mB^{-1}\mN=\begin{pmatrix}0&1\\5&-6\\20&-10\end{pmatrix},\\
\bar c_2=450-\vpi^\T\mA_{\cdot2}=450,\\
\bar c_5=0-\vpi^\T\mA_{\cdot5}=-500,
\end{gathered}
$$

and $z=\vc_B^\T\bar{\vb}=500\cdot8=4000$. The canonical form with respect to
this basis reads

$$
\begin{aligned}
x_1+x_5&=8,\\
5x_2+x_3-6x_5&=12,\\
20x_2+x_4-10x_5&=70,\\
-z+450x_2-500x_5&=-4000 :
\end{aligned}
$$

the BFS $\vx=(8,0,12,70,0)^\T$ produces 8 hundred cases of juice glasses,
leaves 12 machine hours and 70 m³ of storage unused, and earns \$4000 per
week. {prf:ref}`ex-08-02-product-mix-bfs` obtains the same equations by
Gauss–Jordan elimination.
:::

(sec-03-03-matrix-step)=
## The simplex step in matrix form

### Optimality test

Equation {eq}`eq-03-03-objective-reduced` holds for **every** solution of
$\mA\vx=\vb$, not only for the BFS. It says how the objective changes when the
nonbasic variables are raised from zero: $\bar c_j$ is the rate of change of
$z$ per unit of $x_j$ when the basic variables adjust to keep $\mA\vx=\vb$.

:::{prf:theorem} Optimality test
:label: thm-03-03-optimality

Let $\mathcal B$ be a feasible basis of the minimization problem
{eq}`eq-03-03-standard-form` with BFS $\bar\vx$.

1. If $\bar c_j\ge 0$ for all $j$, then $\bar\vx$ is optimal.
2. If $\bar\vx$ is nondegenerate and $\bar c_s<0$ for some $s$, then
   $\bar\vx$ is not optimal.

For a maximization problem, read $\bar c_j\le 0$ in part 1 and
$\bar c_s>0$ in part 2.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

1. Let $\vx$ be any feasible solution. By
   {eq}`eq-03-03-objective-reduced`, $\vc^\T\vx=\vc_B^\T\bar{\vb}+\sum_{j\notin\mathcal B}\bar c_jx_j\ge\vc_B^\T\bar{\vb}=\vc^\T\bar\vx$,
   because $\bar c_j\ge0$ and $x_j\ge 0$.
2. Increase $x_s$ to a small $\theta>0$, keep the other nonbasic variables
   at zero, and set $\vx_B=\bar{\vb}-\theta\bar{\va}_s$. This satisfies
   $\mA\vx=\vb$, and $\vx_B\ge\vzero$ for $\theta$ small enough because
   $\bar{\vb}>\vzero$. Its cost is $\vc^\T\bar\vx+\theta\bar c_s<\vc^\T\bar\vx$.

For a maximization problem apply the result to $\min(-\vc)^\T\vx$, whose
reduced costs are $-\bar c_j$.
:::

Part 1 is a *sufficient* condition. It is not necessary at a degenerate BFS:
then the step of part 2 may have length zero, and an optimal BFS can be
represented by a basis whose reduced costs fail the test.
{prf:ref}`ex-08-02-a3-q16` shows such a basis. The optimal point is the same;
only its basis has to be changed, by pivots of length zero, until the test
holds. (Some texts state the test as "if and only if"; that is true only
under nondegeneracy.)

### Moving along an edge

Suppose the test fails: $\bar c_s<0$ for a nonbasic $x_s$ (for a
maximization problem, $\bar c_s>0$). The simplex method changes the BFS in
the simplest possible way. It raises $x_s$ from zero to some $\theta\ge0$
and keeps every other nonbasic variable at zero. The basic variables then
follow from {eq}`eq-03-03-canonical-constraints`:

$$
\begin{gathered}
\vx_B(\theta)=\bar{\vb}-\theta\,\bar{\va}_s,\\
x_s(\theta)=\theta,\\
z(\theta)=z+\theta\,\bar c_s .
\end{gathered}
$$ (eq-03-03-edge)

In vector form $\vx(\theta)=\bar\vx+\theta\vd$ with $d_s=1$,
$\vd_B=-\bar{\va}_s$ and the other components zero. Geometrically, $\vd$
points along an edge of the feasible polyhedron that leaves $\bar\vx$; the
objective improves at rate $\abs{\bar c_s}$ per unit of $\theta$.

How far can we go? The objective keeps improving, so the only limit is
feasibility, $\vx_B(\theta)\ge\vzero$, that is,
$\bar b_i-\theta\bar a_{is}\ge 0$ for every row $i$. Rows with
$\bar a_{is}\le 0$ impose no limit: their basic variable stays the same or
grows. Rows with $\bar a_{is}>0$ require $\theta\le\bar b_i/\bar a_{is}$.
The largest feasible step is therefore given by the **minimum-ratio test**

$$
\theta=\min\left\{\frac{\bar b_i}{\bar a_{is}}:\ \bar a_{is}>0\right\}
=\frac{\bar b_r}{\bar a_{rs}} .
$$ (eq-03-03-ratio-test)

At $\theta$ the basic variable $x_{B(r)}$ of a minimizing row $r$ has fallen
to zero. It **leaves** the basis and $x_s$ **enters**: the new basis is
$\mathcal B'=\mathcal B\setminus\{B(r)\}\cup\{s\}$. Two questions remain.

*Is the new point again a BFS?* The point $\vx(\theta)$ is feasible by
construction and has $x_{B(r)}=0$ and all nonbasic variables other than
$x_s$ at zero, so it is the basic solution of $\mathcal B'$ provided that
$\mathcal B'$ is a basis. It is: $\mB^{-1}\mB'$ equals the identity matrix
with column $r$ replaced by $\bar{\va}_s$, whose determinant is
$\bar a_{rs}\neq 0$ ({prf:ref}`prop-01-01-basis-exchange`). So $\mB'$ is nonsingular, and $\mathcal B'$ is a
feasible basis adjacent to $\mathcal B$.

*What if no entry of $\bar{\va}_s$ is positive?* Then $\vx(\theta)$ is
feasible for every $\theta\ge0$ while $z(\theta)=z+\theta\bar c_s$ improves
without bound: the LP is **unbounded**, and $\{\bar\vx+\theta\vd:\theta\ge0\}$
is a ray of feasible solutions along which this happens
({ref}`sec-03-03-special-cases`).

The answer to the question "which variable leaves?" is thus forced. Any
other choice either does not give a basis or gives an infeasible basic
solution, as the solution of {prf:ref}`ex-03-03-matrix-step` shows.

:::{prf:example} One simplex step by matrix algebra
:label: ex-03-03-matrix-step

Continue {prf:ref}`ex-03-03-product-mix` at the basis $\mathcal B=\{1,3,4\}$,
BFS $(x_1,x_2)=(8,0)$, $z=\$4000$. The reduced costs are $\bar c_2=450>0$
and $\bar c_5=-500$. Carry out one simplex step by matrix algebra: choose
the entering variable, examine each basic variable that could leave the
basis and compare with the ratio test {eq}`eq-03-03-ratio-test`, and compute
the new basis inverse, multipliers, reduced costs and objective value.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The problem is a maximization, so $x_2$ is the entering variable: each additional hundred cases of cocktail glasses raises the
contribution by \$450 while the basic variables adjust. Its updated column
and the current values are

$$
\begin{gathered}
\bar{\va}_2=\mB^{-1}\mA_{\cdot2}=\begin{pmatrix}0\\5\\20\end{pmatrix},\\
\bar{\vb}=\begin{pmatrix}8\\12\\70\end{pmatrix} \quad\text{(rows of }x_1,\ x_3,\ x_4).
\end{gathered}
$$

Three variables could in principle leave the basis.

- **$x_1$ leaves.** The columns of $x_2,x_3,x_4$ are $(5,20,0)^\T$,
  $(1,0,0)^\T$ and $(0,1,0)^\T$; their third components are all zero, so
  $\mB$ is singular and $\{2,3,4\}$ is not a basis. This is the case
  $\bar a_{12}=0$: raising $x_2$ does not change $x_1$ at all.
- **$x_4$ leaves.** The basis $\{1,3,2\}$ gives
  $(x_1,x_3,x_2)=(8,-\tfrac{11}{2},\tfrac72)$. It is a basic solution but not
  feasible: the point $(x_1,x_2)=(8,\tfrac72)$ lies on the storage line and
  on the demand line, and needs $6\cdot8+5\cdot\tfrac72=65.5$ h of machine
  time, more than the 60 h available. The step $\theta=\tfrac72$ was too
  long.
- **$x_3$ leaves.** The basis $\{1,2,4\}$ gives
  $(x_1,x_2,x_4)=(8,\tfrac{12}{5},22)\ge\vzero$, a BFS.

The ratio test {eq}`eq-03-03-ratio-test` finds the third case directly: the
ratios are $12/5$ (row of $x_3$) and $70/20=7/2$ (row of $x_4$), and the row
of $x_1$ has no ratio because $\bar a_{12}=0$. The minimum
$\theta=\tfrac{12}{5}$ belongs to $x_3$. Along the edge
{eq}`eq-03-03-edge`, $x_3=12-5\theta$ and $x_4=70-20\theta$; machine time runs
out first.

For the new basis, in the order $(x_1,x_2,x_4)$,

$$
\begin{gathered}
\mB'^{-1}=\begin{pmatrix}0&0&1\\\tfrac15&0&-\tfrac65\\-4&1&14\end{pmatrix},\\
\vpi^\T=(500,450,0)\,\mB'^{-1}=(90,0,-40),
\end{gathered}
$$

so $\bar c_3=0-90=-90$ and $\bar c_5=0-(-40)=40$, and
$z=500\cdot8+450\cdot\tfrac{12}{5}=5080$. As {eq}`eq-03-03-edge` predicts,
$z$ grew by $\theta\bar c_2=\tfrac{12}{5}\cdot450=1080$. The new BFS
$(8,\tfrac{12}{5})$ is not yet optimal: $\bar c_5=40>0$ says that giving up
juice-glass production in favor of more cocktail glasses still pays.
:::

Recomputing $\mB^{-1}$ at every step is wasteful. The tableau below updates
all the quantities of the canonical form with one elimination step instead;
the revised simplex method of {doc}`03-06-revised-dual-simplex` returns to the
matrix form and updates $\mB^{-1}$ in product form
({ref}`sec-03-06-product-form`, {eq}`eq-01-01-product-form`).

(sec-03-03-tableau)=
## The simplex tableau

### Layout and sign conventions

A **simplex tableau** stores the canonical form
({prf:ref}`def-03-03-canonical-form`) as a table of detached coefficients. In
this book it has one row per constraint, labeled by its basic variable,
and the objective row last, labeled $-z$:

| Basic | $\vx_B$ | $\vx_N$ | RHS |
|:--|:-:|:-:|:-:|
| $\vx_B$ | $\mI$ | $\mB^{-1}\mN$ | $\bar{\vb}=\mB^{-1}\vb$ |
| $-z$ | $\vzero^\T$ | $\bar{\vc}_N^\T=\vc_N^\T-\vpi^\T\mN$ | $-\vc_B^\T\bar{\vb}$ |

The entries of the constraint rows are the $\bar a_{ij}$, the right-hand
sides (column RHS) are the $\bar b_i$, the objective row holds the reduced costs
$\bar c_j$ (zero in basic columns), and its right-hand side is $-z$. The
label $-z$ reminds us that this row is the equation $-z+\sum_j\bar c_jx_j=-z_0$,
where $z_0$ is the current objective value.

The rules depend only on the sign of $\bar c_j$:

```{list-table} Rules of the tableau simplex method
:name: tab-03-03-rules
:header-rows: 1
:stub-columns: 1

* - 
  - Minimization
  - Maximization
* - Optimal when
  - all $\bar c_j\ge0$
  - all $\bar c_j\le0$
* - Entering variable (Dantzig's rule)
  - most negative $\bar c_j$
  - most positive $\bar c_j$
* - Unbounded when
  - $\bar c_s<0$ and $\bar a_{is}\le0$ for all $i$
  - $\bar c_s>0$ and $\bar a_{is}\le0$ for all $i$
* - Leaving variable
  - minimum ratio {eq}`eq-03-03-ratio-test`
  - minimum ratio {eq}`eq-03-03-ratio-test`
* - Change of $z$ in a pivot
  - $\theta\bar c_s\le0$
  - $\theta\bar c_s\ge0$
```

Some textbooks, Taha (2017) among them, use a **$z$-row** instead: the
objective is written $z-\vc^\T\vx=0$, so the row holds $-\bar c_j$ and the
right-hand side $+z$, and every sign test above is reversed (for
maximization, optimal when all $z$-row entries are $\ge0$, entering
variable most negative). Both conventions give the same pivots; this book
uses the $-z$ row throughout.

### The pivot operation

Moving from basis $\mathcal B$ to $\mathcal B'$ means making the column of
$x_s$ a unit vector with its $1$ in row $r$, while keeping the other basic
columns unit vectors. This is one step of Gauss–Jordan elimination
({prf:ref}`alg-01-01-pivot`) on the **pivot element** $\bar a_{rs}$:

$$
\begin{aligned}
\text{row } r:&\quad \bar a'_{rj}=\frac{\bar a_{rj}}{\bar a_{rs}},\\
&\quad \bar b'_r=\frac{\bar b_r}{\bar a_{rs}}=\theta,\\
\text{row } i\ne r:&\quad \bar a'_{ij}=\bar a_{ij}-\bar a_{is}\,\bar a'_{rj},\\
&\quad \bar b'_i=\bar b_i-\bar a_{is}\,\theta,\\
\text{objective row}:&\quad \bar c'_j=\bar c_j-\bar c_s\,\bar a'_{rj},\\
&\quad -z'=-z-\bar c_s\,\theta .
\end{aligned}
$$ (eq-03-03-pivot)

The objective row is updated by the same row operation as the other rows.
The result is again the canonical form: the row operations replace the
equations by equivalent ones, the columns of the basic variables are unit
vectors, and the objective row is zero in them. Since the canonical form
with respect to a basis is unique (it is $\mB'^{-1}\mA\vx=\mB'^{-1}\vb$ and
{eq}`eq-03-03-objective-reduced`), the updated objective row holds exactly the
new reduced costs, and its right-hand side is $-z'$ with

$$
z'=z+\theta\,\bar c_s .
$$ (eq-03-03-z-change)

Since $\theta\ge0$, a pivot never makes the objective worse, and it
improves it strictly when $\theta>0$. The value $\theta=0$ occurs exactly
when the leaving variable is already zero, at a degenerate BFS
({ref}`sec-03-03-degeneracy`).

The pivot also produces $\mB'^{-1}$ for free: the columns that formed the
identity matrix in the initial tableau (here the slack columns) hold
$\mB^{-1}$ in every later tableau, because they record the row operations
done so far. The same columns of the objective row hold $-\vpi^\T$ when the
initial basic variables have zero cost.

(sec-03-03-algorithm)=
## The simplex algorithm

:::{prf:algorithm} Simplex method (tableau form)
:label: alg-03-03-simplex

**Input:** an LP in standard form and a tableau in canonical form with
respect to a feasible basis $\mathcal B$.
**Output:** an optimal BFS, or a ray along which the objective is unbounded.

1. **Optimality test.** If $\bar c_j\ge0$ for all $j$ (maximization:
   $\bar c_j\le0$), stop: the current BFS is optimal
   ({prf:ref}`thm-03-03-optimality`).
2. **Entering variable.** Choose a nonbasic $x_s$ with $\bar c_s<0$
   (maximization: $\bar c_s>0$). Dantzig's rule takes the most negative
   (most positive) $\bar c_s$; Bland's rule ({prf:ref}`alg-03-03-bland`) takes
   the smallest index.
3. **Unboundedness test.** If $\bar a_{is}\le0$ for every row $i$, stop: the
   LP is unbounded along the ray $\bar\vx+\theta\vd$, $\theta\ge0$, of
   {eq}`eq-03-03-edge`.
4. **Ratio test.** Compute $\theta=\min\{\bar b_i/\bar a_{is}:\bar a_{is}>0\}$
   and choose a row $r$ that attains it (ties: the row whose basic variable
   has the smallest index). The variable $x_{B(r)}$ leaves.
5. **Pivot.** Pivot on $\bar a_{rs}$ by {eq}`eq-03-03-pivot`, label row $r$ with
   $x_s$, and go to Step 1.
:::

Step 2 leaves a choice. Dantzig's rule is the classical one and usually
needs few iterations, but it does not always minimize the number of pivots
({prf:ref}`ex-08-02-a3-q8` compares it with the least-positive rule). Any rule
that picks an improving variable is correct; the choice matters for speed
and, at degenerate vertices, for termination.

:::{prf:example} The product-mix LP solved by the simplex method
:label: ex-03-03-simplex-run

Solve {prf:ref}`ex-03-03-product-mix` from the slack basis with Dantzig's rule.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The pivot element of each tableau is printed in bold.

```{table} Initial tableau of the product-mix LP (basis $\{x_3,x_4,x_5\}$, vertex $(0,0)$).
:name: tab-03-03-pm-0

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $6$ | $5$ | $1$ | $0$ | $0$ | $60$ |
| $x_4$ | $10$ | $20$ | $0$ | $1$ | $0$ | $150$ |
| $x_5$ | $\boldsymbol{1}$ | $0$ | $0$ | $0$ | $1$ | $8$ |
| $-z$ | $500$ | $450$ | $0$ | $0$ | $0$ | $0$ |
```

**Iteration 1.** The largest reduced cost is $\bar c_1=500$, so $x_1$
enters. The ratios are $60/6=10$, $150/10=15$ and $8/1=8$; the minimum
$\theta=8$ is in the row of $x_5$, so $x_5$ leaves: the demand limit is
reached before machine time or storage run out. Pivoting on $1$ subtracts
6 times, 10 times and 500 times the $x_5$ row from the other rows, and
$z$ rises by $8\cdot500=4000$.

```{table} Tableau 1 (basis $\{x_3,x_4,x_1\}$, vertex $(8,0)$).
:name: tab-03-03-pm-1

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_3$ | $0$ | $\boldsymbol{5}$ | $1$ | $0$ | $-6$ | $12$ |
| $x_4$ | $0$ | $20$ | $0$ | $1$ | $-10$ | $70$ |
| $x_1$ | $1$ | $0$ | $0$ | $0$ | $1$ | $8$ |
| $-z$ | $0$ | $450$ | $0$ | $0$ | $-500$ | $-4000$ |
```

This is the canonical form found by matrix algebra in
{prf:ref}`ex-03-03-product-mix`, with the rows in another order.

**Iteration 2.** Now $\bar c_2=450$ is the only positive entry, so $x_2$
enters. The ratios are $12/5$ and $70/20=7/2$ (the row of $x_1$ has
$\bar a_{12}=0$ and no ratio); $x_3$ leaves with $\theta=\tfrac{12}{5}$, as
in {prf:ref}`ex-03-03-matrix-step`, and $z$ rises by
$\tfrac{12}{5}\cdot450=1080$ to 5080.

```{table} Tableau 2 (basis $\{x_2,x_4,x_1\}$, vertex $(8,\tfrac{12}{5})$).
:name: tab-03-03-pm-2

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $0$ | $1$ | $\tfrac{1}{5}$ | $0$ | $-\tfrac{6}{5}$ | $\tfrac{12}{5}$ |
| $x_4$ | $0$ | $0$ | $-4$ | $1$ | $\boldsymbol{14}$ | $22$ |
| $x_1$ | $1$ | $0$ | $0$ | $0$ | $1$ | $8$ |
| $-z$ | $0$ | $0$ | $-90$ | $0$ | $40$ | $-5080$ |
```

**Iteration 3.** The slack $x_5$ of the demand limit has $\bar c_5=40>0$
and enters: the shop should produce *fewer* juice glasses. The column of
$x_5$ is $(-\tfrac65,14,1)^\T$; the first entry is negative (as $x_5$ grows,
$x_2$ grows too), and the ratios are $22/14=\tfrac{11}{7}$ and $8/1=8$. So
$x_4$ leaves with $\theta=\tfrac{11}{7}$ and $z$ rises by
$\tfrac{11}{7}\cdot40=\tfrac{440}{7}$.

```{table} Tableau 3, optimal (basis $\{x_2,x_5,x_1\}$, vertex $(\tfrac{45}{7},\tfrac{30}{7})$).
:name: tab-03-03-pm-3

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_2$ | $0$ | $1$ | $-\tfrac{1}{7}$ | $\tfrac{3}{35}$ | $0$ | $\tfrac{30}{7}$ |
| $x_5$ | $0$ | $0$ | $-\tfrac{2}{7}$ | $\tfrac{1}{14}$ | $1$ | $\tfrac{11}{7}$ |
| $x_1$ | $1$ | $0$ | $\tfrac{2}{7}$ | $-\tfrac{1}{14}$ | $0$ | $\tfrac{45}{7}$ |
| $-z$ | $0$ | $0$ | $-\tfrac{550}{7}$ | $-\tfrac{20}{7}$ | $0$ | $-\tfrac{36000}{7}$ |
```

**Optimality.** All reduced costs are $\le0$, so Tableau 3 is optimal:
$x_1^*=\tfrac{45}{7}\approx6.429$ and $x_2^*=\tfrac{30}{7}\approx4.286$
hundred cases per week, with contribution
$z^*=\tfrac{36000}{7}\approx\$5142.86$ per week. Machine time and storage
are used up ($x_3=x_4=0$), and the demand limit has
$x_5=\tfrac{11}{7}\approx1.571$ hundred cases to spare. The nonbasic
reduced costs are strictly negative, so this optimum is unique
({ref}`sec-03-03-special-cases`).

*(An earlier version of this tableau shows the entry in row $x_2$, column
$x_3$ as $\tfrac17$; the correct entry is $-\tfrac17$, which is also the
only value consistent with the objective row:
$0-\left(450\cdot(-\tfrac17)+500\cdot\tfrac27\right)=-\tfrac{550}{7}$.)*

**Checks.** The solution satisfies $6\cdot\tfrac{45}{7}+5\cdot\tfrac{30}{7}=60$,
$10\cdot\tfrac{45}{7}+20\cdot\tfrac{30}{7}=150$ and $\tfrac{45}{7}\le8$, and
$500\cdot\tfrac{45}{7}+450\cdot\tfrac{30}{7}=\tfrac{36000}{7}$; the steps
$4000+1080+\tfrac{440}{7}=\tfrac{36000}{7}$ add up. In the slack columns of
the objective row stand $-\pi_1=-\tfrac{550}{7}$ and $-\pi_2=-\tfrac{20}{7}$:
one more machine hour would be worth $\tfrac{550}{7}\approx\$78.57$ and
one more cubic meter of storage $\tfrac{20}{7}\approx\$2.86$ per week. These
shadow prices are the subject of {doc}`03-05-lp-duality`
({prf:ref}`ex-03-05-product-mix-dual`, {prf:ref}`thm-03-05-shadow-price`) and
{doc}`03-07-sensitivity-analysis`.

```{figure} figures/03-03-simplex-path.svg
:name: fig-03-03-simplex-path
:width: 520px
:class: dark-light
:alt: Feasible region of the product-mix LP with the simplex path from the origin along the x1-axis to (8,0), up to (8,12/5) and on to the optimum (45/7,30/7)

Path of the simplex method on the product-mix LP. Each tableau is a
vertex: $\vx^{(0)}=(0,0)$, $\vx^{(1)}=(8,0)$, $\vx^{(2)}=(8,\tfrac{12}{5})$
and the optimum $\vx^{(3)}=(\tfrac{45}{7},\tfrac{30}{7})$, where the dashed
contour $500x_1+450x_2=\tfrac{36000}{7}$ touches the region. Every pivot
moves along an edge to an adjacent vertex. The hollow point
$(8,\tfrac72)$ is the infeasible basic solution obtained if $x_4$ instead
of $x_3$ leaves in the second step ({prf:ref}`ex-03-03-matrix-step`).
```
:::

On the path of {numref}`fig-03-03-simplex-path`, Dantzig's rule picked $x_1$
first because its reduced cost was larger, and the method needed three
pivots. Had $x_2$ entered first, the path $(0,0)\to(0,\tfrac{15}{2})\to(\tfrac{45}{7},\tfrac{30}{7})$
would have reached the optimum in two: the largest reduced cost measures
the improvement per unit of the entering variable, not per pivot.

(sec-03-03-special-cases)=
## Termination, unboundedness and alternative optima

### Finite termination

:::{prf:theorem} Finite termination without degeneracy
:label: thm-03-03-finite

If every BFS of the LP is nondegenerate, then
{prf:ref}`alg-03-03-simplex`, with any rule for choosing an improving entering
variable, stops after finitely many iterations, either with an optimal
basis or with a ray along which the objective is unbounded.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Take a minimization problem. In every iteration that does not stop,
$\bar c_s<0$ and, by nondegeneracy, $\theta=\bar b_r/\bar a_{rs}>0$, so by
{eq}`eq-03-03-z-change` the objective strictly decreases. The objective value
of a BFS is determined by its basis, so no basis can occur twice. There are
at most $\binom{n}{m}$ bases, hence at most that many iterations, and the
algorithm can only stop in Step 1 (optimal) or Step 3 (unbounded). When it
stops in Step 1, the BFS is optimal by {prf:ref}`thm-03-03-optimality`.
:::

With an initial feasible basis ({doc}`03-04-two-phase-big-m` shows how to find
one), this proves constructively the fundamental theorem of LP
({prf:ref}`thm-03-02-fundamental`) for nondegenerate problems: an LP that is feasible and bounded has an optimal
BFS. The bound $\binom{n}{m}$ is enormous, and examples are known (Klee &
Minty, 1972) on which Dantzig's rule visits exponentially many vertices; in practice the number of
pivots is usually a small multiple of $m$.

### Unboundedness

Step 3 of the algorithm stops when the entering column has no positive
entry. By {eq}`eq-03-03-edge` the points $\bar\vx+\theta\vd$ are feasible for
all $\theta\ge0$, and the objective changes by $\theta\bar c_s$, without
limit. The feasible region contains the whole ray, and the LP has no optimal
solution. Two remarks:

- An unbounded *feasible region* does not imply an unbounded *objective*.
  Unboundedness of the LP is a property of the region and the cost vector
  together, and the tableau detects it only for the entering column it
  examines; another improving column may still have a positive entry.
- In a modeling context an unbounded LP almost always means a missing
  constraint (an unlimited resource, a forgotten demand limit).

{prf:ref}`ex-08-02-unbounded-ray` and {prf:ref}`ex-08-02-a3-q15` work through
unbounded LPs.

### Alternative optima

In an optimal tableau of a minimization problem, suppose some nonbasic
$x_j$ has $\bar c_j=0$. Entering $x_j$ does not change $z$ by
{eq}`eq-03-03-edge`, so every point reached is also optimal:

- if the ratio test gives $\theta>0$, the pivot leads to a different
  optimal BFS, and the whole edge between the two vertices is optimal;
- if $\theta=0$ (degeneracy), the pivot changes the basis but not the point;
- if $\bar{\va}_j\le\vzero$, the ray $\bar\vx+\theta\vd$, $\theta\ge0$,
  consists of optimal solutions, none of them a vertex except
  $\bar\vx$ itself ({prf:ref}`ex-08-02-a3-q14`).

Conversely, if every nonbasic $\bar c_j$ is strictly positive (strictly
negative for maximization), the optimum is **unique**: by
{eq}`eq-03-03-objective-reduced`, any feasible $\vx$ with $\vx_N\neq\vzero$
costs more than $\bar\vx$, and $\vx_N=\vzero$ forces $\vx_B=\bar{\vb}$. At a
degenerate optimum the converse fails: there the optimum can be unique even
though some $\bar c_j=0$.

(sec-03-03-degeneracy)=
## Degeneracy and cycling

### Degenerate pivots

A BFS is degenerate when a basic variable is zero. Geometrically, more
constraints are active at the vertex than are needed to define it: in the
plane, three or more constraint lines (counting the axes) pass through it,
and several bases describe the same point. Degeneracy arises in the
algorithm through a **tie** in the ratio test: when two rows attain the
minimum ratio, only one of their basic variables leaves, and the other stays
basic with value zero.

At a degenerate BFS the next ratio test may give $\theta=0$: the
leaving variable is a basic variable already at zero, and the pivot changes
the basis without moving the point or changing $z$. Such a **degenerate
pivot** is not wasted. It may be needed to reach a basis that passes the
optimality test, or to find an edge along which the method can move again.
{prf:ref}`ex-08-02-a3-q13` shows such *temporary degeneracy*, and
{prf:ref}`ex-08-02-a3-q12` counts the bases that can describe one degenerate
vertex.

### Cycling

The proof of {prf:ref}`thm-03-03-finite` breaks down at degenerate vertices:
if $\theta=0$ the objective does not decrease, and nothing prevents the
method from returning to a basis it has already visited. Once that happens
with a deterministic rule, the same sequence repeats forever. This is
**cycling**. It is rare in practice but possible, as the following classic
example shows; it is the version of Beale's (1955) example given by Bertsimas
and Tsitsiklis (1997, Example 3.6).

:::{prf:example} Cycling under Dantzig's rule
:label: ex-03-03-cycling

Consider

$$
\begin{aligned}
\min_{\vx}\ & z=-\tfrac34x_4+20x_5-\tfrac12x_6+6x_7\\
\st\ & x_1+\tfrac14x_4-8x_5-x_6+9x_7=0,\\
& x_2+\tfrac12x_4-12x_5-\tfrac12x_6+3x_7=0,\\
& x_3+x_6=1,\\
& x_1,\dots,x_7\ge0,
\end{aligned}
$$

which is in canonical form with respect to $\{x_1,x_2,x_3\}$:

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_1$ | $1$ | $0$ | $0$ | $\boldsymbol{\tfrac{1}{4}}$ | $-8$ | $-1$ | $9$ | $0$ |
| $x_2$ | $0$ | $1$ | $0$ | $\tfrac{1}{2}$ | $-12$ | $-\tfrac{1}{2}$ | $3$ | $0$ |
| $x_3$ | $0$ | $0$ | $1$ | $0$ | $0$ | $1$ | $0$ | $1$ |
| $-z$ | $0$ | $0$ | $0$ | $-\tfrac{3}{4}$ | $20$ | $-\tfrac{1}{2}$ | $6$ | $0$ |

The BFS $(0,0,1,0,0,0,0)$ is degenerate. Use Dantzig's rule (most negative
$\bar c_j$ enters) and break ties in the ratio test by the smallest index of
the leaving variable. Show that the simplex method then cycles.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Exact pivoting gives the following sequence:

```{list-table} Six degenerate pivots that return to the starting basis
:name: tab-03-03-cycling
:header-rows: 1

* - Iteration
  - Entering ($\bar c_s$)
  - Rows in the tie
  - Leaving
  - Pivot element
  - Basis after the pivot
* - 1
  - $x_4$ ($-\tfrac34$)
  - $x_1,x_2$
  - $x_1$
  - $\tfrac14$
  - $\{x_4,x_2,x_3\}$
* - 2
  - $x_5$ ($-4$)
  - $x_2$
  - $x_2$
  - $4$
  - $\{x_4,x_5,x_3\}$
* - 3
  - $x_6$ ($-2$)
  - $x_4,x_5$
  - $x_4$
  - $8$
  - $\{x_6,x_5,x_3\}$
* - 4
  - $x_7$ ($-3$)
  - $x_5$
  - $x_5$
  - $\tfrac{3}{16}$
  - $\{x_6,x_7,x_3\}$
* - 5
  - $x_1$ ($-1$)
  - $x_6,x_7$
  - $x_6$
  - $2$
  - $\{x_1,x_7,x_3\}$
* - 6
  - $x_2$ ($-2$)
  - $x_7$
  - $x_7$
  - $\tfrac13$
  - $\{x_1,x_2,x_3\}$
```

Every pivot has $\theta=0$, the point and the value $z=0$ never change,
and after six pivots the tableau is again the initial one. The method would
repeat these six pivots forever, although the LP has the optimal solution
$\vx^*=(\tfrac34,0,0,1,0,1,0)$ with $z^*=-\tfrac54$.
:::

### Bland's rule

Cycling can be prevented by a careful choice of the entering and leaving
variables. The simplest such choice is due to Bland (1977).

:::{prf:algorithm} Bland's smallest-subscript rule
:label: alg-03-03-bland

In Steps 2 and 4 of {prf:ref}`alg-03-03-simplex`:

1. **Entering variable.** Among the nonbasic variables with
   $\bar c_j<0$ (maximization: $\bar c_j>0$), choose the one with the
   **smallest index** $j$.
2. **Leaving variable.** Among the rows that attain the minimum ratio,
   choose the one whose basic variable has the **smallest index**.
:::

:::{prf:theorem} Bland (1977)
:label: thm-03-03-bland

With Bland's rule the simplex method never cycles; it stops after finitely
many iterations with an optimal basis or with an unbounded ray, whether or
not the LP is degenerate.
:::

The proof (Bland, 1977; Bertsimas & Tsitsiklis, 1997, §3.4) assumes a
cycle and looks at the variable with the largest index among those that
enter and leave during the cycle. Comparing the tableau in which it leaves
with the tableau in which it enters again leads to a sign contradiction in
the reduced costs. Applied to {prf:ref}`ex-03-03-cycling`, Bland's rule makes
the same first three pivots. In the fourth tableau, however, both $x_1$
($\bar c_1=-2$) and $x_7$ ($\bar c_7=-3$) are improving, and Bland's rule
takes $x_1$ where Dantzig's rule took $x_7$:

| Iteration | Entering | Leaving | $\theta$ | $z$ after the pivot |
|:-:|:-:|:-:|:-:|:-:|
| 4 | $x_1$ | $x_5$ | $0$ | $0$ |
| 5 | $x_2$ | $x_3$ | $\tfrac12$ | $-\tfrac12$ |
| 6 | $x_4$ | $x_2$ | $1$ | $-\tfrac54$ |

The sixth tableau,

| Basic | $x_1$ | $x_2$ | $x_3$ | $x_4$ | $x_5$ | $x_6$ | $x_7$ | RHS |
|:--|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| $x_6$ | $0$ | $0$ | $1$ | $0$ | $0$ | $1$ | $0$ | $1$ |
| $x_1$ | $1$ | $-\tfrac{1}{2}$ | $\tfrac{3}{4}$ | $0$ | $-2$ | $0$ | $\tfrac{15}{2}$ | $\tfrac{3}{4}$ |
| $x_4$ | $0$ | $2$ | $1$ | $1$ | $-24$ | $0$ | $6$ | $1$ |
| $-z$ | $0$ | $\tfrac{3}{2}$ | $\tfrac{5}{4}$ | $0$ | $2$ | $0$ | $\tfrac{21}{2}$ | $\tfrac{5}{4}$ |

has all reduced costs $\ge0$ and gives the optimum
$\vx^*=(\tfrac34,0,0,1,0,1,0)$, $z^*=-\tfrac54$ (confirmed by an LP solver).

Bland's rule is simple but often slow, because it ignores the size of the
reduced costs. Two other remedies are classical. The **lexicographic rule**
breaks ties in the ratio test by comparing whole rows of
$[\bar{\vb}\ \ \mB^{-1}]$, divided by the pivot-column entry,
lexicographically; it is equivalent to Charnes's (1952) idea of
**perturbing** the right-hand sides to $b_i+\varepsilon^i$ for a tiny
$\varepsilon>0$, which makes every BFS nondegenerate (Dantzig, 1963, ch. 10;
Bertsimas & Tsitsiklis, 1997, §3.4). Commercial solvers use small random
perturbations of $\vb$ or of the bounds, which have the same effect in
floating-point arithmetic.

## Summary and Common Mistakes

### Summary

- Relative to a feasible basis, the LP takes the canonical form
  ({prf:ref}`def-03-03-canonical-form`): $\vx_B=\bar{\vb}-\mB^{-1}\mN\vx_N$ and
  $z=\vc_B^\T\bar{\vb}+\sum_{j\notin\mathcal B}\bar c_jx_j$, with reduced
  costs $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ and $\vpi^\T=\vc_B^\T\mB^{-1}$
  ({eq}`eq-03-03-objective-reduced`, {eq}`eq-03-03-reduced-cost`).
- A BFS is optimal for a minimization if all $\bar c_j\ge0$ (for a
  maximization, all $\bar c_j\le0$); the test is also necessary at a
  nondegenerate BFS ({prf:ref}`thm-03-03-optimality`).
- An improving variable $x_s$ is raised along the edge {eq}`eq-03-03-edge`;
  the minimum-ratio test {eq}`eq-03-03-ratio-test` gives the longest feasible
  step and the leaving variable, the new basis is adjacent, and
  $z'=z+\theta\bar c_s$ ({eq}`eq-03-03-z-change`).
- In the tableau (objective row last, labeled $-z$, holding $\bar c_j$
  and $-z$), one Gauss–Jordan pivot {eq}`eq-03-03-pivot` performs the whole
  step ({prf:ref}`alg-03-03-simplex`, {numref}`tab-03-03-rules`).
- The product-mix LP is solved in three pivots along the path of
  {numref}`fig-03-03-simplex-path`, with optimum $(\tfrac{45}{7},\tfrac{30}{7})$
  and $z^*=\tfrac{36000}{7}\approx\$5142.86$ per week.
- Without degeneracy the method terminates ({prf:ref}`thm-03-03-finite`). An
  improving column without positive entries proves unboundedness; zero
  reduced costs of nonbasic variables at the optimum signal alternative
  optima.
- Degenerate pivots ($\theta=0$) can make the method cycle
  ({prf:ref}`ex-03-03-cycling`); Bland's smallest-subscript rule
  ({prf:ref}`alg-03-03-bland`, {prf:ref}`thm-03-03-bland`) and lexicographic or
  perturbation rules prevent it.

### Common mistakes

- **Wrong sign test.** In the $-z$ row of this book, a *maximization* is
  optimal when all $\bar c_j\le0$ and the entering variable has the most
  *positive* $\bar c_j$; with Taha's $z$-row the signs are reversed. Mixing
  the two conventions in one problem leads to pivots that make $z$ worse.
- **Ratios with nonpositive entries.** Rows with $\bar a_{is}\le0$ are
  left out of the ratio test. A negative ratio, or a ratio with a zero
  denominator, does not bound the step.
- **Choosing the leaving variable by the largest ratio** or by the largest
  pivot. Any row other than a minimizing one gives an infeasible basic
  solution ({prf:ref}`ex-03-03-matrix-step`).
- **Forgetting the objective row** when pivoting, or updating it with a
  different multiplier. The objective row is just one more row of the
  elimination step {eq}`eq-03-03-pivot`.
- **Reading the value from the objective row with the wrong sign.** The
  right-hand side of the $-z$ row is $-z$.
- **Declaring a degenerate tableau non-optimal** because some
  $\bar c_j$ fails the test, or declaring the optimum unique because the
  test holds strictly. Both conclusions need nondegeneracy.
- **Rounding.** Decimal rounding in a tableau accumulates quickly; keep
  fractions exact when working by hand.

## References

- Beale, E. M. L. (1955). Cycling in the dual simplex algorithm. *Naval
  Research Logistics Quarterly, 2*(4), 269–275. (The cycling example.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Sections 3.1–3.3: optimality
  conditions, development of the simplex method, the full tableau and
  Example 3.6 on cycling; Section 3.4: lexicographic and Bland's rules;
  Section 3.7: computational efficiency.)
- Bland, R. G. (1977). New finite pivoting rules for the simplex method.
  *Mathematics of Operations Research, 2*(2), 103–107.
  <https://doi.org/10.1287/moor.2.2.103>
- Charnes, A. (1952). Optimality and degeneracy in linear programming.
  *Econometrica, 20*(2), 160–170. <https://doi.org/10.2307/1907845>
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapter 5: the simplex method; Chapter 10: resolving
  degeneracy.)
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 4: the simplex method, its
  algebra and tableau form, tie breaking.)
- Klee, V., & Minty, G. J. (1972). How good is the simplex algorithm? In
  O. Shisha (Ed.), *Inequalities III* (pp. 159–175). Academic Press.
  (Exponentially many pivots under Dantzig's rule.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Sections 3.1–3.3: equation form, from the graphical to the
  algebraic solution, the simplex method with the $z$-row convention;
  Section 3.5: degeneracy, alternative optima and unboundedness.)
