# Convex Sets

## Overview

The feasible set $\Omega$ of an optimization problem is a subset of
$\R^n$, and its shape decides much of what can be said about the problem.
This chapter collects the facts about sets that the rest of the book uses.
The first part is topological: open, closed, bounded and compact sets,
sequences and their limit points, and Weierstrass' theorem
({prf:ref}`thm-02-01-weierstrass`), which guarantees that a continuous
function attains its minimum on a nonempty compact set. The second part is
geometric: affine and convex sets, convex hulls and cones, hyperplanes and
half-spaces, polyhedra, balls and ellipsoids, the cone of positive
semidefinite matrices, the operations that preserve convexity, and extreme
points. Seven example sets ({prf:ref}`ex-02-01-seven-sets`) are examined
against all these properties. The chapter ends with the separating and
supporting hyperplane theorems ({prf:ref}`thm-02-01-separation-point`,
{prf:ref}`thm-02-01-supporting-hyperplane`), the geometric facts behind
duality and optimality conditions.
It uses the vectors, norms, inner products and definite matrices of
{doc}`01-01-linear-algebra` and the feasible sets and optimal values of
{doc}`01-03-introduction`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- decide whether a set in $\R^n$ is open, closed, bounded or compact, and
  find the limit points of a sequence;
- apply Weierstrass' theorem to show that an optimization problem has an
  optimal solution;
- recognize affine sets, convex sets, cones, hyperplanes, half-spaces,
  polyhedra, balls, ellipsoids and norm cones, and prove that they are
  convex;
- form convex hulls and conic hulls, and identify extreme points;
- use intersections and affine maps to prove that a set is convex; and
- state the separating and supporting hyperplane theorems and explain the
  proof of the first.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\vx,\vy,\vz$; $\vx_1,\dots,\vx_k$ | points of $\R^n$ (a bold letter with a subscript is a point, not a component) |
| $\vx^{(k)}$ | $k$-th term of a sequence of points |
| $B(\vx,r)$, $\bar B(\vx,r)$ | open ball $\{\vy:\norm{\vy-\vx}<r\}$ and closed ball $\{\vy:\norm{\vy-\vx}\le r\}$ |
| $\intr\Omega$, $\operatorname{cl}\Omega$, $\operatorname{bd}\Omega$ | interior, closure and boundary of $\Omega$ |
| $\theta$, $\theta_i$ | weights of an affine, convex or conic combination |
| $\aff S$, $\conv S$, $\cone S$ | affine hull, convex hull and conic hull of a set $S$ |
| $C$, $K$ | a convex set; a cone |
| $H$, $P$ | a hyperplane; a polyhedron |
| $\R^n_+$ | nonnegative orthant $\{\vx\in\R^n:\vx\ge\vzero\}$ |
| $\mathcal E$ | an ellipsoid |
| $\Sym^n$, $\Sym^n_+$, $\Sym^n_{++}$ | symmetric, positive semidefinite and positive definite $n\times n$ matrices |
| $\operatorname{proj}_C(\vy)$ | the point of $C$ closest to $\vy$ |

Unless stated otherwise, $\norm{\cdot}$ is the Euclidean norm
$\norm{\vx}=\sqrt{\vx^\T\vx}$.

(sec-02-01-topology)=
## Open, closed, bounded and compact sets

Optimization problems ask for the best point of a set, and whether a best
point exists depends on whether the set contains its "edge" and whether it
extends to infinity. Balls make these ideas precise. The **open ball** of
radius $r>0$ around $\vx$ is $B(\vx,r)=\{\vy\in\R^n:\norm{\vy-\vx}<r\}$, and
the **closed ball** is $\bar B(\vx,r)=\{\vy:\norm{\vy-\vx}\le r\}$.

:::{prf:definition} Interior, open and closed sets
:label: def-02-01-open-closed

Let $\Omega\subseteq\R^n$.

1. A point $\vx\in\Omega$ is an **interior point** of $\Omega$ if some ball
   $B(\vx,\varepsilon)$, $\varepsilon>0$, lies entirely in $\Omega$. The
   set of interior points is the **interior** $\intr\Omega$.
2. $\Omega$ is **open** if every point of $\Omega$ is an interior point,
   $\Omega=\intr\Omega$.
3. $\Omega$ is **closed** if its complement $\R^n\setminus\Omega$ is open.
4. A point $\vx$ (in $\Omega$ or not) is a **boundary point** of $\Omega$
   if every ball around $\vx$ contains points of $\Omega$ and points not in
   $\Omega$. The set of boundary points is the **boundary**
   $\operatorname{bd}\Omega$, and $\operatorname{cl}\Omega=\Omega\cup
   \operatorname{bd}\Omega$ is the **closure** of $\Omega$.
:::

An open set contains none of its boundary points, a closed set contains all
of them. Most sets are neither: the interval $[0,1)$ contains the boundary
point $0$ but not the boundary point $1$. "Closed" is therefore not the
opposite of "open". The empty set and $\R^n$ are both open and closed. In
optimization, sets defined by finitely many inequalities $g_i(\vx)\le0$
with continuous $g_i$ are closed, while strict inequalities $g_i(\vx)<0$
define open sets.

### Sequences and limit points

A **sequence** in $\R^n$ is an infinite list of points
$\vx^{(0)},\vx^{(1)},\vx^{(2)},\dots$, written $\{\vx^{(k)}\}$. Iterative
algorithms ({doc}`06-04-line-search`, {doc}`06-05-unconstrained-algorithms`)
produce such sequences, and their analysis asks whether and where the
sequence settles.

:::{prf:definition} Convergence, subsequence, limit point
:label: def-02-01-limit-point

1. $\{\vx^{(k)}\}$ **converges** to $\vx^*$, written $\vx^{(k)}\to\vx^*$,
   if $\norm{\vx^{(k)}-\vx^*}\to0$ as $k\to\infty$.
2. If $k_1<k_2<k_3<\cdots$ are increasing indices, the sequence
   $\vx^{(k_1)},\vx^{(k_2)},\dots$ is a **subsequence** of $\{\vx^{(k)}\}$;
   for example the terms with even index.
3. $\vx^*$ is a **limit point** of $\{\vx^{(k)}\}$ if some subsequence
   converges to $\vx^*$.
:::

A convergent sequence has exactly one limit point, its limit. A sequence can
also have several limit points, or none.

:::{prf:example} A sequence with two limit points
:label: ex-02-01-limit-points

Find the limit points of the sequence in $\R^2$ with terms

$$
\vx^{(k)}=\left(\frac{1}{k+1},\ (-1)^k\,2+\frac{1}{k^2-0.5}\right)^\T,
\qquad k=0,1,2,\dots
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The first component tends to $0$. The second depends on the parity of $k$.
For even $k=2s$ it equals $2+1/(4s^2-0.5)$, which tends to $2$; for odd
$k=2s+1$ it equals $-2+1/((2s+1)^2-0.5)$, which tends to $-2$. Hence the
subsequence of even-indexed terms converges to $(0,2)$ and the subsequence
of odd-indexed terms converges to $(0,-2)$: for example
$\vx^{(10)}=(0.0909,\ 2.0101)$ and $\vx^{(11)}=(0.0833,\ -1.9917)$.

There are no other limit points. Any convergent subsequence contains
infinitely many even-indexed terms or infinitely many odd-indexed terms,
and these terms approach $(0,2)$ or $(0,-2)$, so the limit of the
subsequence must be one of these two points. The limit points are
therefore $(0,2)$ and $(0,-2)$. Because there are two of them, the sequence
itself does not converge.
:::

Closed sets can be described by sequences, and this description is often
easier to use than the definition.

:::{prf:theorem} Closed sets contain their limits
:label: thm-02-01-closed-sequences

A set $\Omega\subseteq\R^n$ is closed if and only if the limit of every
convergent sequence of points of $\Omega$ belongs to $\Omega$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\Omega$ be closed, $\vx^{(k)}\in\Omega$ and $\vx^{(k)}\to\vx^*$. If
$\vx^*\notin\Omega$, then $\vx^*$ lies in the open complement, so some ball
$B(\vx^*,\varepsilon)$ contains no point of $\Omega$. But
$\norm{\vx^{(k)}-\vx^*}<\varepsilon$ for all large $k$, a contradiction.

Conversely, suppose $\Omega$ contains the limits of its convergent
sequences, and let $\vy\notin\Omega$. If every ball $B(\vy,1/k)$ contained a
point $\vx^{(k)}\in\Omega$, then $\vx^{(k)}\to\vy$ and $\vy$ would belong to
$\Omega$. So some ball around $\vy$ misses $\Omega$, that is, the complement
is open.
:::

### Bounded and compact sets

:::{prf:definition} Bounded and compact sets
:label: def-02-01-bounded-compact

A set $\Omega\subseteq\R^n$ is **bounded** if there is a number $M$ with
$\norm{\vx}\le M$ for all $\vx\in\Omega$ (equivalently, $\abs{x_j}\le M'$ for
all components and some $M'$). It is **compact** if it is closed and
bounded.
:::

The key property of compact sets is that sequences in them cannot escape.
By the **Bolzano–Weierstrass theorem**, every bounded sequence in $\R^n$ has
at least one limit point; if the sequence lies in a compact set $\Omega$,
the limit point lies in $\Omega$ by {prf:ref}`thm-02-01-closed-sequences`.
The sequence of {prf:ref}`ex-02-01-limit-points` is bounded, and indeed it has
limit points; the sequence $\vx^{(k)}=(k,0)$ is unbounded and has none.

(sec-02-01-weierstrass)=
### Existence of a minimum

:::{prf:theorem} Weierstrass
:label: thm-02-01-weierstrass

Let $\Omega\subseteq\R^n$ be nonempty and compact and let
$f:\Omega\to\R$ be continuous. Then $f$ attains a minimum and a maximum on
$\Omega$: there are points $\vx^*,\vx^\sharp\in\Omega$ with
$f(\vx^*)\le f(\vx)\le f(\vx^\sharp)$ for all $\vx\in\Omega$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $p^*=\inf\{f(\vx):\vx\in\Omega\}$, possibly $-\infty$. Choose points
$\vx^{(k)}\in\Omega$ with $f(\vx^{(k)})\to p^*$. Because $\Omega$ is
compact, a subsequence converges to some $\vx^*\in\Omega$, and continuity
of $f$ gives $f(\vx^*)=\lim f(\vx^{(k_s)})=p^*$ along that subsequence. So
$p^*$ is finite and attained. Apply the same argument to $-f$ for the
maximum.
:::

Each hypothesis is needed, as three small examples show. On the unbounded
set $[1,\infty)$ the function $1/x$ has infimum $0$ but never reaches it
({doc}`01-03-introduction`). On the set $(0,1]$, which is not closed, $f(x)=x$
has infimum $0$, not attained. On the compact set $[0,1]$, the
discontinuous function with $f(0)=1$ and $f(x)=x$ for $x>0$ again has
infimum $0$, not attained.

Many feasible sets in practice are closed but unbounded, for example
$\{\vx:\vx\ge\vzero\}$. The theorem still applies if the objective grows
without bound far away. Call $f$ **coercive** on $\Omega$ if
$f(\vx)\to\infty$ whenever $\norm{\vx}\to\infty$ with $\vx\in\Omega$.

:::{prf:corollary} Coercive functions attain a minimum
:label: cor-02-01-coercive

If $\Omega$ is nonempty and closed and $f$ is continuous and coercive on
$\Omega$, then $f$ attains a minimum on $\Omega$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Pick $\vx^{(0)}\in\Omega$. The set $\{\vx\in\Omega:f(\vx)\le
f(\vx^{(0)})\}$ is nonempty, closed (a continuous inequality within a
closed set) and bounded (by coercivity), hence compact. A minimizer of $f$
on this set, which exists by {prf:ref}`thm-02-01-weierstrass`, minimizes $f$
on all of $\Omega$, because the points left out have larger values.
:::

For instance, the squared distance $\norm{\vx-\vy}^2$ to a fixed point
$\vy$ is coercive, so every nonempty closed set contains a point closest to
$\vy$. This fact is used in {ref}`sec-02-01-separation`.

(sec-02-01-affine-sets)=
## Lines, segments and affine sets

Let $\vx\ne\vy$ be two points of $\R^n$. The points

$$
\vz=\theta\vx+(1-\theta)\vy=\vy+\theta(\vx-\vy),\qquad\theta\in\R,
$$ (eq-02-01-line)

form the **line** through $\vx$ and $\vy$: starting at $\vy$ ($\theta=0$)
one moves along the direction $\vx-\vy$ and reaches $\vx$ at $\theta=1$.
Restricting $\theta$ to $[0,1]$ gives the **line segment** between $\vx$ and
$\vy$, written $[\vx,\vy]$.

:::{prf:definition} Affine set
:label: def-02-01-affine-set

A set $A\subseteq\R^n$ is **affine** if it contains the line through any two
distinct points of $A$. A point $\theta_1\vx_1+\dots+\theta_k\vx_k$ with
$\theta_1+\dots+\theta_k=1$ is an **affine combination** of
$\vx_1,\dots,\vx_k$.
:::

An affine set contains every affine combination of its points (by induction
on $k$, as in {prf:ref}`thm-02-01-convex-combinations` below). The most
important example is the solution set of a system of linear equations.

:::{prf:proposition} Solution sets of linear equations are affine
:label: prop-02-01-linear-equations-affine

For $\mA\in\R^{m\times n}$ and $\vb\in\R^m$, the set
$\{\vx:\mA\vx=\vb\}$ is affine.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If $\mA\vx=\vb$, $\mA\vy=\vb$ and $\theta\in\R$, then
$\mA(\theta\vx+(1-\theta)\vy)=\theta\vb+(1-\theta)\vb=\vb$.
:::

The converse also holds: every affine set can be written as
$\{\vx:\mA\vx=\vb\}$ for some $\mA$ and $\vb$. If $A$ is affine and
$\vx_0\in A$, then $V=A-\vx_0=\{\vx-\vx_0:\vx\in A\}$ is a linear subspace,
and $A=\vx_0+V$ is a translated subspace; its **dimension** is that of
$V$. In the language of {doc}`01-01-linear-algebra`, $V$ is the null space of
$\mA$ and has dimension $n-\rank\mA$. The **affine hull** $\aff S$ of a set
$S$ is the set of all affine combinations of points of $S$, the smallest
affine set containing $S$. The affine hull of two distinct points is the
line through them; that of three points not on a line is a plane.

(sec-02-01-convex-sets)=
## Convex sets

:::{prf:definition} Convex set
:label: def-02-01-convex-set

A set $C\subseteq\R^n$ is **convex** if it contains the line segment between
any two of its points:

$$
\vx,\vy\in C,\ \theta\in[0,1]\ \Longrightarrow\ \theta\vx+(1-\theta)\vy\in C .
$$ (eq-02-01-convex-set)
:::

Every affine set is convex, since it contains whole lines and not only
segments. The empty set, a single point, a ray $\{\vx_0+t\vd:t\ge0\}$ and
$\R^n$ are convex. {numref}`fig-02-01-convex-nonconvex` shows a convex set and
three sets that are not convex: for each of the latter, two points of the
set can be found whose segment leaves the set.

```{figure} figures/02-01-convex-nonconvex.svg
:name: fig-02-01-convex-nonconvex
:width: 640px
:class: dark-light
:alt: A shaded hexagon with a segment inside it, and three nonconvex sets (a kidney shape, a ring, two hexagons touching at one corner) each with a segment that leaves the set

(a) A convex hexagon: the segment between any two of its points stays in
it. (b) A kidney-shaped set, (c) a ring and (d) two hexagons that touch at
one vertex are not convex: the dashed parts of the segments lie outside the
set.
```

### Convex combinations and the convex hull

A point $\theta_1\vx_1+\dots+\theta_k\vx_k$ with $\theta_i\ge0$ and
$\theta_1+\dots+\theta_k=1$ is a **convex combination** of
$\vx_1,\dots,\vx_k$: a weighted average of the points. For $k=2$ these are
the points of the segment $[\vx_1,\vx_2]$; for three points not on a line
they fill the triangle with these vertices.

:::{prf:theorem} Convex sets contain all convex combinations
:label: thm-02-01-convex-combinations

A set $C$ is convex if and only if every convex combination of finitely
many points of $C$ belongs to $C$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

If $C$ contains all convex combinations, it contains those of two points,
so it is convex. Conversely let $C$ be convex; we show by induction on $k$
that convex combinations of $k$ points of $C$ lie in $C$. For $k=1$ there is
nothing to show, and $k=2$ is the definition. Let
$\vz=\sum_{i=1}^k\theta_i\vx_i$ with $k\ge3$. If $\theta_k=1$, then
$\vz=\vx_k\in C$. Otherwise

$$
\vz=(1-\theta_k)\sum_{i=1}^{k-1}\frac{\theta_i}{1-\theta_k}\vx_i+\theta_k\vx_k .
$$

The weights $\theta_i/(1-\theta_k)$, $i<k$, are nonnegative and sum to
$1$, so the first sum is a convex combination of $k-1$ points and lies in
$C$ by the induction hypothesis. Then $\vz$ is a convex combination of two
points of $C$ and lies in $C$.
:::

:::{prf:definition} Convex hull
:label: def-02-01-convex-hull

The **convex hull** $\conv S$ of a set $S\subseteq\R^n$ is the set of all
convex combinations of finitely many points of $S$.
:::

$\conv S$ is convex, since a convex combination of two convex combinations
is again a convex combination of points of $S$, and by
{prf:ref}`thm-02-01-convex-combinations` every convex set that contains $S$
contains $\conv S$. Hence $\conv S$ is the *smallest* convex set containing
$S$. For a finite set $S$ it is a polygon (in $\R^2$) or a polytope (in
general) whose corners are points of $S$; see {numref}`fig-02-01-hulls`(a).
**Carathéodory's theorem** adds that in $\R^n$ every point of $\conv S$ is a
convex combination of at most $n+1$ points of $S$; in the plane, every point
of the hull lies in a triangle with corners in $S$.

```{figure} figures/02-01-hulls.svg
:name: fig-02-01-hulls
:width: 600px
:class: dark-light
:alt: Left, nine points and the shaded polygon that is their convex hull, with its corner points marked. Right, two vectors from the origin and the shaded wedge between them that is their conic hull

(a) The convex hull of nine points is the shaded polygon; its corners
(filled) are the extreme points, the other points (open) lie inside or on
an edge. (b) The conic hull of $\vx_1$ and $\vx_2$ is the wedge of all
combinations $\theta_1\vx_1+\theta_2\vx_2$ with $\theta_1,\theta_2\ge0$.
```

(sec-02-01-cones)=
## Cones

:::{prf:definition} Cone and convex cone
:label: def-02-01-cone

A set $K\subseteq\R^n$ is a **cone** if $\vx\in K$ and $\theta\ge0$ imply
$\theta\vx\in K$. It is a **convex cone** if, in addition, it is convex;
equivalently, $\theta_1\vx_1+\theta_2\vx_2\in K$ for all
$\vx_1,\vx_2\in K$ and $\theta_1,\theta_2\ge0$.
:::

A point $\theta_1\vx_1+\dots+\theta_k\vx_k$ with all $\theta_i\ge0$ (but no
condition on their sum) is a **conic combination** of the $\vx_i$. The
**conic hull** $\cone S$ is the set of all conic combinations of points of
$S$, the smallest convex cone containing $S$ ({numref}`fig-02-01-hulls`(b)).
Geometrically, a cone is a union of rays from the origin; a convex cone is
a cone without "gaps".

Examples:

- the nonnegative orthant $\R^n_+=\{\vx:\vx\ge\vzero\}$;
- a half-space through the origin, $\{\vx:\va^\T\vx\le0\}$;
- the cone generated by the columns of a matrix,
  $\{\mA\vw:\vw\ge\vzero\}$, which is the conic hull of the columns
  $\mA_{\cdot1},\dots,\mA_{\cdot n}$. Cones of this kind appear in Farkas'
  lemma ({prf:ref}`thm-03-05-farkas`) and in the geometric form
  {eq}`eq-06-03-kkt-cone` of the Karush–Kuhn–Tucker (KKT) conditions, where the negative gradient of the
  objective must lie in the cone generated by the gradients of the active
  constraints;
- the norm cones and the positive semidefinite cone of the next sections.

The union of the two half-axes $\{(t,0):t\ge0\}\cup\{(0,t):t\ge0\}$ is a cone
but not a convex cone: $(1,0)+(0,1)$ is missing.

(sec-02-01-hyperplanes)=
## Hyperplanes and half-spaces

:::{prf:definition} Hyperplane and half-space
:label: def-02-01-hyperplane

For $\va\in\R^n$, $\va\ne\vzero$, and $b\in\R$, the set
$H=\{\vx:\va^\T\vx=b\}$ is a **hyperplane** with **normal vector** $\va$. The
sets $\{\vx:\va^\T\vx\le b\}$ and $\{\vx:\va^\T\vx\ge b\}$ are the two
closed **half-spaces** bounded by $H$; with strict inequalities they are
open half-spaces.
:::

In $\R^2$ a hyperplane is a line, in $\R^3$ a plane. If $\vx_0$ is any
point of $H$, so that $\va^\T\vx_0=b$, then

$$
\begin{aligned}
H&=\{\vx:\va^\T(\vx-\vx_0)=0\},\\
\{\vx:\va^\T\vx\le b\}&=\{\vx:\va^\T(\vx-\vx_0)\le0\}.
\end{aligned}
$$ (eq-02-01-hyperplane-offset)

The first set consists of the points $\vx$ for which $\vx-\vx_0$ is
orthogonal to $\va$: $H$ is the affine set $\vx_0+\{\va\}^\perp$ of
dimension $n-1$. The half-space $\va^\T\vx\le b$ consists of the points for
which $\vx-\vx_0$ makes an obtuse or right angle with $\va$, that is, the
points on the side of $H$ *away* from $\va$ ({numref}`fig-02-01-hyperplane`).
Moving from $H$ in the direction of $\va$ increases $\va^\T\vx$; this is
the geometric content of the statement that the gradient $\va$ of the
linear function $\va^\T\vx$ points in the direction of increase
({prf:ref}`thm-01-02-steepest-ascent`).

```{figure} figures/02-01-hyperplane.svg
:name: fig-02-01-hyperplane
:width: 460px
:class: dark-light
:alt: A line in the plane with a normal vector a drawn at a point x0 on it; the side away from a is labeled a-transpose x at most b, the other side at least b

The hyperplane $\va^\T\vx=b$ through $\vx_0$ with normal vector $\va$. For
a point $\vx$ in the shaded half-space $\va^\T\vx\le b$, the vector
$\vx-\vx_0$ makes an obtuse angle with $\va$.
```

A hyperplane is affine and hence convex. A half-space is convex: if
$\va^\T\vx\le b$, $\va^\T\vy\le b$ and $\theta\in[0,1]$, then
$\va^\T(\theta\vx+(1-\theta)\vy)=\theta\va^\T\vx+(1-\theta)\va^\T\vy\le\theta
b+(1-\theta)b=b$. Half-spaces are the building blocks of the feasible sets of linear programs (LPs):
each linear inequality constraint defines one.

(sec-02-01-polyhedra)=
## Polyhedra

:::{prf:definition} Polyhedron and polytope
:label: def-02-01-polyhedron

A **polyhedron** is the intersection of finitely many closed half-spaces,

$$
P=\{\vx\in\R^n:\va_i^\T\vx\le b_i,\ i=1,\dots,m\}=\{\vx:\mA\vx\le\vb\},
$$ (eq-02-01-polyhedron)

where the rows of $\mA$ are $\va_i^\T$ and $\le$ between vectors holds
component by component. A bounded polyhedron is a **polytope**.
:::

Equality constraints are allowed, because each equation is a pair of
inequalities: $\va^\T\vx=b$ holds exactly when $\va^\T\vx\le b$ and
$-\va^\T\vx\le-b$. So $\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$, the feasible set
of an LP in standard form, is a polyhedron, and so is every set described by
finitely many linear equations and inequalities. Other examples are the
nonnegative orthant $\R^n_+$, a box $\{\vx:l_j\le x_j\le u_j\}$, a half-space,
an affine set, and the **unit simplex** $\{\vx:\vx\ge\vzero,\
\vone^\T\vx=1\}$ (the set of probability vectors).

Every polyhedron is convex, as an intersection of convex sets
({ref}`sec-02-01-operations`), and closed, as an intersection of closed sets. A
polytope is therefore compact, and by {prf:ref}`thm-02-01-weierstrass` a
continuous function, in particular a linear one, attains its minimum on a
nonempty polytope.

Polyhedra also have an "inner" description. The **Minkowski–Weyl theorem**
states that a set is a polytope if and only if it is the convex hull of
finitely many points, and that every polyhedron is the sum of the convex
hull of finitely many points and the conic hull of finitely many
directions. The fundamental theorem of LP ({prf:ref}`thm-03-02-fundamental`)
shows that an LP in standard form with an optimal solution has an optimal
corner point.

(sec-02-01-balls)=
## Balls, ellipsoids and norm cones

The closed Euclidean ball with center $\vx_c$ and radius $r>0$ is

$$
\bar B(\vx_c,r)=\{\vx:\norm{\vx-\vx_c}\le r\}=\{\vx_c+r\vu:\norm{\vu}\le1\}.
$$ (eq-02-01-ball)

It is convex: if $\norm{\vx-\vx_c}\le r$, $\norm{\vy-\vx_c}\le r$ and
$\theta\in[0,1]$, the triangle inequality gives

$$
\begin{aligned}
\norm{\theta\vx+(1-\theta)\vy-\vx_c}
&=\norm{\theta(\vx-\vx_c)+(1-\theta)(\vy-\vx_c)}\\
&\le\theta\norm{\vx-\vx_c}+(1-\theta)\norm{\vy-\vx_c}\\
&\le r .
\end{aligned}
$$

The proof used only the triangle inequality and homogeneity, so the ball
$\{\vx:\norm{\vx-\vx_c}\le r\}$ of *any* norm is convex. In $\R^2$ the unit
ball of $\norm{\cdot}_1$ is a diamond and that of $\norm{\cdot}_\infty$ a
square ({numref}`fig-02-01-balls`(a)). Both are polytopes: the $\ell_1$ ball is
the convex hull of the $2n$ points $\pm\ve_j$, the $\ell_\infty$ ball is the
box $\{-1\le x_j\le1\}$ with $2^n$ corners. The Euclidean ball is not a
polyhedron: its boundary is curved, and finitely many half-spaces cannot
reproduce it.

:::{prf:definition} Ellipsoid
:label: def-02-01-ellipsoid

For $\vx_c\in\R^n$ and a symmetric positive definite matrix
$\mP\in\Sym^n_{++}$, the set

$$
\mathcal E=\{\vx:(\vx-\vx_c)^\T\mP^{-1}(\vx-\vx_c)\le1\}
$$ (eq-02-01-ellipsoid)

is an **ellipsoid** with center $\vx_c$.
:::

With $\mP=r^2\mI$ the ellipsoid is the ball $\bar B(\vx_c,r)$. In general,
let $\mP=\mV\diag(\lambda_1(\mP),\dots,\lambda_n(\mP))\mV^\T$ be the
spectral decomposition of $\mP$ ({doc}`01-01-linear-algebra`). The semi-axes of
$\mathcal E$ point along the eigenvectors (columns of $\mV$) and have
lengths $\sqrt{\lambda_i(\mP)}$ ({numref}`fig-02-01-balls`(b)). Equivalently,
$\mathcal E=\{\vx_c+\mP^{1/2}\vu:\norm{\vu}\le1\}$ is the image of the unit
ball under an affine map, which proves its convexity by
{prf:ref}`thm-02-01-affine-maps`. Ellipsoids describe, for example, the
uncertainty region of estimated demands, and they are the sublevel sets of
convex quadratic functions ({prf:ref}`thm-02-02-sublevel`).

```{figure} figures/02-01-balls.svg
:name: fig-02-01-balls
:width: 600px
:class: dark-light
:alt: Left, the unit balls of the 1-norm (diamond), 2-norm (disc) and infinity-norm (square) drawn on top of each other. Right, an ellipse with its center and two perpendicular semi-axes

(a) Unit balls of $\norm{\cdot}_1$ (diamond), $\norm{\cdot}_2$ (disc) and
$\norm{\cdot}_\infty$ (square) in $\R^2$. (b) The ellipsoid
{eq}`eq-02-01-ellipsoid` with $\vx_c=(2,1)$ and $\mP$ with eigenvalues $4$ and
$1$; the semi-axes are $2\vv_1$ and $\vv_2$, where $\vv_1$ and $\vv_2$ are
unit eigenvectors of $\mP$ for these eigenvalues.
```

Adding the radius as a variable gives a cone. The **norm cone** of a norm
$\norm{\cdot}$ is

$$
K=\{(\vx,t)\in\R^{n+1}:\norm{\vx}\le t\}.
$$ (eq-02-01-norm-cone)

It is a convex cone: scaling $(\vx,t)$ by $\theta\ge0$ preserves
$\norm{\vx}\le t$, and convexity follows from the triangle inequality as for
the ball. For the Euclidean norm, $K$ is the **second-order cone** (or
"ice-cream cone"); in $\R^3$ it is the solid circular cone
$\sqrt{x_1^2+x_2^2}\le t$. Optimization problems with constraints of the form
$\norm{\mA\vx+\vb}\le\vc^\T\vx+d$ (second-order cone programs) can be solved
almost as efficiently as LPs (Boyd & Vandenberghe, 2004, §4.4.2).

(sec-02-01-psd-cone)=
## The positive semidefinite cone

Sets of matrices can be convex too. Let $\Sym^n$ be the set of symmetric
$n\times n$ matrices, $\Sym^n_+$ the positive semidefinite ones
($\mA\succeq0$, that is, $\vz^\T\mA\vz\ge0$ for all $\vz$) and $\Sym^n_{++}$
the positive definite ones ({doc}`01-01-linear-algebra`).

$\Sym^n_+$ is a convex cone: if $\mA,\mB\succeq0$ and
$\theta_1,\theta_2\ge0$, then for every $\vz$

$$
\vz^\T(\theta_1\mA+\theta_2\mB)\vz=\theta_1\vz^\T\mA\vz+\theta_2\vz^\T\mB\vz\ge0 .
$$

For $n=2$, the matrix $\begin{pmatrix}u&v\\v&w\end{pmatrix}$ is positive
semidefinite if and only if $u\ge0$, $w\ge0$ and $uw-v^2\ge0$. In the
coordinates $(u,v,w)\in\R^3$ this set is a solid cone with a curved
boundary, so $\Sym^2_+$ is a convex cone but not a polyhedron. Constraints
of the form "a matrix depending affinely on $\vx$ is positive semidefinite"
define semidefinite programs (Boyd & Vandenberghe, 2004, §4.6.2); in this
book the cone appears mainly through the second-order test for convex
functions ({prf:ref}`thm-02-02-second-order`).

(sec-02-01-operations)=
## Operations that preserve convexity

Checking definition {eq}`eq-02-01-convex-set` directly is often awkward. It is
usually easier to build a set from simple convex sets by operations that
preserve convexity.

:::{prf:theorem} Intersections
:label: thm-02-01-intersection

The intersection of any collection (finite or infinite) of convex sets is
convex.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\vx,\vy$ lie in every set of the collection and $\theta\in[0,1]$. Each
set is convex, so each contains $\theta\vx+(1-\theta)\vy$, and so does the
intersection.
:::

Polyhedra are intersections of half-spaces and are therefore convex. The
theorem also covers infinite intersections:
$\Sym^n_+=\bigcap_{\vz\in\R^n}\{\mA\in\Sym^n:\vz^\T\mA\vz\ge0\}$ is an
intersection of infinitely many half-spaces of $\Sym^n$, one for each
$\vz$. The feasible set $\{\vx:g_i(\vx)\le0,\ i=1,\dots,m\}$ of a problem
is convex whenever each set $\{\vx:g_i(\vx)\le0\}$ is, which is the case
for convex functions $g_i$ ({prf:ref}`thm-02-02-sublevel`).

A **union** of convex sets is in general *not* convex; set $\Omega_3$ of
{prf:ref}`ex-02-01-seven-sets` is the union of two convex polygons and is not
convex.

:::{prf:theorem} Affine maps
:label: thm-02-01-affine-maps

Let $\phi(\vx)=\mA\vx+\vb$ with $\mA\in\R^{m\times n}$, $\vb\in\R^m$. If
$C\subseteq\R^n$ is convex, then its image $\phi(C)=\{\mA\vx+\vb:\vx\in C\}$
is convex. If $D\subseteq\R^m$ is convex, then its preimage
$\phi^{-1}(D)=\{\vx:\mA\vx+\vb\in D\}$ is convex.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

An affine map preserves convex combinations:
$\phi(\theta\vx+(1-\theta)\vy)=\theta\phi(\vx)+(1-\theta)\phi(\vy)$ because
$\theta\vb+(1-\theta)\vb=\vb$. If $\phi(\vx),\phi(\vy)\in\phi(C)$ with
$\vx,\vy\in C$, the combination $\theta\phi(\vx)+(1-\theta)\phi(\vy)$ is
the image of the point $\theta\vx+(1-\theta)\vy\in C$. If
$\phi(\vx),\phi(\vy)\in D$, then
$\phi(\theta\vx+(1-\theta)\vy)\in D$ because $D$ is convex.
:::

Consequences:

- translations $C+\vx_0$, scalings $\alpha C$ and rotations of convex sets
  are convex;
- the projection of a convex set onto some of its coordinates,
  $\{\vx_1:(\vx_1,\vx_2)\in C\}$, is convex (image under a linear map);
- the sum $C_1+C_2=\{\vx+\vy:\vx\in C_1,\vy\in C_2\}$ of two convex sets is
  convex (image of the product $C_1\times C_2$, which is convex, under
  $(\vx,\vy)\mapsto\vx+\vy$);
- the set $\{\vx:\norm{\mA\vx-\vb}\le r\}$ is convex (preimage of a ball);
- the ellipsoid {eq}`eq-02-01-ellipsoid` is convex (image of the unit ball).

Two further operations, the perspective map and linear-fractional maps,
also preserve convexity (Boyd & Vandenberghe, 2004, §2.3.3); they are not
needed in this book.

(sec-02-01-seven-sets)=
## Seven example sets

The following example applies the definitions of this chapter to a family
of sets, five of them in the plane ({numref}`fig-02-01-example-sets`).

```{figure} figures/02-01-example-sets.svg
:name: fig-02-01-example-sets
:width: 720px
:class: dark-light
:alt: Five small panels showing the quadrilateral Omega 1, the disc Omega 2, the nonconvex quadrilateral Omega 3 with a dashed segment outside it, Omega 4 with two dashed open edges, and the half-plane Omega 5

The sets $\Omega_1$ to $\Omega_5$ of {prf:ref}`ex-02-01-seven-sets`. Dashed
edges of $\Omega_4$ do not belong to the set. In $\Omega_3$ the segment
from $(0,0)$ to $(3,0)$ leaves the set, so $\Omega_3$ is not convex.
$\Omega_5$ continues without bound to the upper right.
```

:::{prf:example} Properties of seven sets
:label: ex-02-01-seven-sets

Decide for each set whether it is open, closed, bounded, compact, convex and
a polyhedron.

$$
\begin{aligned}
\Omega_1&=\{\vx\in\R^2:x_1+x_2\le2,\ x_1-x_2\le1,\ x_1\ge0,\ x_2\ge0\},\\
\Omega_2&=\{\vx\in\R^2:(x_1-1)^2+(x_2-1)^2\le1\},\\
\Omega_3&=\{\vx\in\R^2:x_1+x_2\le3,\ g(\vx)\le0,\ x_1\ge0,\ x_2\ge0\},\\
\Omega_4&=\{\vx\in\R^2:x_1+x_2<2,\ x_1-x_2<1,\ x_1\ge0,\ x_2\ge0\},\\
\Omega_5&=\{\vx\in\R^2:x_1+x_2\ge2\},\\
\Omega_6&=\R^n_+=\{\vx\in\R^n:x_j\ge0,\ j=1,\dots,n\},\\
\Omega_7&=\R^n,
\end{aligned}
$$

where $g(\vx)=x_1-x_2$ for $0\le x_1\le1$ and $g(\vx)=3-x_1-2x_2$ for
$x_1\ge1$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**$\Omega_1$** is the quadrilateral with corners $(0,0)$, $(1,0)$,
$(\tfrac32,\tfrac12)$ and $(0,2)$; the corner $(\tfrac32,\tfrac12)$ solves
$x_1+x_2=2$, $x_1-x_2=1$. It is the intersection of four closed half-spaces,
hence a polyhedron, closed and convex. It is bounded ($\norm{\vx}\le2$ on
it), so it is a polytope and compact. It is not open: the boundary point
$(0,0)$ belongs to it.

**$\Omega_2$** is the closed disc with center $(1,1)$ and radius $1$,
touching both axes. It is closed, bounded, compact and convex
({ref}`sec-02-01-balls`), not open, and not a polyhedron because its boundary
is a circle.

**$\Omega_3$** is bounded by the line $x_1+x_2=3$ above and by the broken
line from $(0,0)$ to $(1,1)$ to $(3,0)$ below (where $g(\vx)=0$): it is the
quadrilateral with corners $(0,0)$, $(1,1)$, $(3,0)$, $(0,3)$, whose
interior angle at $(1,1)$ is about $252^\circ$. It is closed (a union of
two closed polygons, the parts with $x_1\le1$ and with $x_1\ge1$) and
bounded, hence compact. It is **not convex**: $(0,0)$ and $(3,0)$ belong to
$\Omega_3$, but their midpoint $(1.5,0)$ does not, since at $x_1=1.5$ the
set requires $x_2\ge(3-1.5)/2=0.75$. Consequently it is not a polyhedron.
It is not open.

**$\Omega_4$** is $\Omega_1$ without the two slanted edges. It is convex,
as the intersection of two open and two closed half-spaces, and bounded. It
is neither open nor closed. It is not open because the point $(0.5,0)$
belongs to it but every ball around it contains points with $x_2<0$. It is
not closed because the points $(\tfrac32-\tfrac1k,\tfrac12)$,
$k=2,3,\dots$, belong to $\Omega_4$ while their limit $(\tfrac32,\tfrac12)$
does not ({prf:ref}`thm-02-01-closed-sequences`). Hence it is not compact and
not a polyhedron (polyhedra are closed). Its closure is $\Omega_1$.

**$\Omega_5$** is a closed half-plane: closed, convex, a polyhedron, not
open. It is unbounded (it contains $(t,t)$ for all $t\ge1$), so it is not
compact.

**$\Omega_6$**, the nonnegative orthant, is closed, convex, a polyhedron and
a convex cone; it is unbounded and not open.

**$\Omega_7=\R^n$** is both open and closed, convex, and unbounded. It is a
polyhedron, since it has the form {eq}`eq-02-01-polyhedron` with $m=0$
inequalities (or with the single row $\vzero^\T\vx\le0$), and it is also an
affine set and a cone.

{numref}`tab-02-01-seven-sets` collects the answers.

```{list-table} Properties of the seven example sets $\Omega_1$ to $\Omega_7$
:name: tab-02-01-seven-sets
:header-rows: 1

* - Set
  - Open
  - Closed
  - Bounded
  - Compact
  - Convex
  - Polyhedron
* - $\Omega_1$
  - no
  - yes
  - yes
  - yes
  - yes
  - yes (polytope)
* - $\Omega_2$
  - no
  - yes
  - yes
  - yes
  - yes
  - no
* - $\Omega_3$
  - no
  - yes
  - yes
  - yes
  - no
  - no
* - $\Omega_4$
  - no
  - no
  - yes
  - no
  - yes
  - no
* - $\Omega_5$
  - no
  - yes
  - no
  - no
  - yes
  - yes
* - $\Omega_6$
  - no
  - yes
  - no
  - no
  - yes
  - yes
* - $\Omega_7$
  - yes
  - yes
  - no
  - no
  - yes
  - yes
```
:::

By {prf:ref}`thm-02-01-weierstrass` and the answers in
{numref}`tab-02-01-seven-sets`, every continuous objective attains its
minimum on $\Omega_1$, $\Omega_2$ and $\Omega_3$. On $\Omega_4$ it may not:
$f(\vx)=-x_1$ approaches $-\tfrac32$ but never reaches it. On $\Omega_5$ to
$\Omega_7$ existence depends on the objective, as in
{prf:ref}`cor-02-01-coercive`.

(sec-02-01-extreme-points)=
## Extreme points

:::{prf:definition} Extreme point
:label: def-02-01-extreme-point

A point $\vx$ of a convex set $C$ is an **extreme point** of $C$ if it is
not a convex combination of two other points of $C$: whenever
$\vx=\theta\vy+(1-\theta)\vz$ with $\vy,\vz\in C$ and $\theta\in(0,1)$, then
$\vy=\vz=\vx$.
:::

An extreme point does not lie strictly inside any segment of $C$. The
extreme points of a polygon are its corners; those of the convex hull of a
finite set are the corners of the hull ({numref}`fig-02-01-hulls`(a)). Every
boundary point of a disc is an extreme point, so $\Omega_2$ has infinitely
many. A half-space, a line and $\R^n$ have none. The orthant $\R^n_+$ has
exactly one, the origin. In {prf:ref}`ex-02-01-seven-sets`, $\Omega_1$ has the
four extreme points $(0,0)$, $(1,0)$, $(\tfrac32,\tfrac12)$, $(0,2)$, while
$\Omega_4$ keeps only $(0,0)$: the other three corners lie on the removed
edges.

Extreme points matter because they are enough to describe a compact convex
set and to find the minimum of a linear function on it. The theorem of
Minkowski (in infinite dimensions, Krein and Milman) states that a compact
convex set in $\R^n$ is the convex hull of its extreme points. Consequently,
if $\vx^*$ minimizes $\vc^\T\vx$ over such a set and
$\vx^*=\sum_i\theta_i\vv_i$ is a convex combination of extreme points
$\vv_i$, then

$$
\vc^\T\vx^*=\sum_i\theta_i\,\vc^\T\vv_i\ \ge\ \min_i\vc^\T\vv_i ,
$$

so some extreme point is also optimal. For LP feasible sets, which are
polyhedra and possibly unbounded, {doc}`03-02-lp-geometry` proves the
corresponding result ({prf:ref}`thm-03-02-fundamental`) and identifies the
extreme points with the basic feasible solutions on which the simplex
method works ({prf:ref}`thm-03-02-extreme-bfs`).

(sec-02-01-separation)=
## Separating and supporting hyperplanes

A convex set and a point outside it can always be separated by a
hyperplane. This simple-looking fact underlies duality theory
({prf:ref}`thm-02-03-slater`), Farkas' lemma ({prf:ref}`thm-03-05-farkas`) and
the KKT conditions ({prf:ref}`thm-06-03-kkt`).
The proof uses the closest point of a set.

:::{prf:theorem} Projection onto a closed convex set
:label: thm-02-01-projection

Let $C\subseteq\R^n$ be nonempty, closed and convex, and let $\vy\in\R^n$.
There is exactly one point $\vx^*\in C$ closest to $\vy$, written
$\vx^*=\operatorname{proj}_C(\vy)$, and a point $\vx^*\in C$ is this point
if and only if

$$
(\vy-\vx^*)^\T(\vx-\vx^*)\le0\quad\text{for all }\vx\in C .
$$ (eq-02-01-projection)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

*Existence.* The function $\norm{\vx-\vy}^2$ is continuous and coercive, so
it attains a minimum on the closed set $C$ by {prf:ref}`cor-02-01-coercive`.

*Characterization.* Let $\vx^*$ be a closest point, $\vx\in C$ and
$t\in(0,1]$. The point $\vx^*+t(\vx-\vx^*)$ lies in $C$ by convexity, so

$$
\begin{aligned}
\norm{\vy-\vx^*}^2
&\le\norm{\vy-\vx^*-t(\vx-\vx^*)}^2\\
&=\norm{\vy-\vx^*}^2-2t(\vy-\vx^*)^\T(\vx-\vx^*)+t^2\norm{\vx-\vx^*}^2 .
\end{aligned}
$$

Canceling, dividing by $2t$ and letting $t\to0$ gives
{eq}`eq-02-01-projection`. Conversely, if $\vx^*\in C$ satisfies
{eq}`eq-02-01-projection`, then for every $\vx\in C$

$$
\begin{aligned}
\norm{\vy-\vx}^2&=\norm{\vy-\vx^*}^2+\norm{\vx^*-\vx}^2-2(\vy-\vx^*)^\T(\vx-\vx^*)\\
&\ge\norm{\vy-\vx^*}^2 ,
\end{aligned}
$$

so $\vx^*$ is a closest point.

*Uniqueness.* If $\vx^*$ and $\vx^{**}$ both satisfy
{eq}`eq-02-01-projection`, apply it to $\vx^*$ with $\vx=\vx^{**}$ and to
$\vx^{**}$ with $\vx=\vx^*$ and add: $\norm{\vx^{**}-\vx^*}^2\le0$.
:::

Condition {eq}`eq-02-01-projection` says that the angle between $\vy-\vx^*$ and
every direction $\vx-\vx^*$ into the set is at least $90^\circ$. It already
contains the separating hyperplane.

:::{prf:theorem} Separation of a point from a closed convex set
:label: thm-02-01-separation-point

Let $C\subseteq\R^n$ be nonempty, closed and convex, and let $\vy\notin C$.
Then there are $\va\ne\vzero$ and $b$ such that

$$
\va^\T\vx<b<\va^\T\vy\quad\text{for all }\vx\in C .
$$ (eq-02-01-strict-separation)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\vx^*=\operatorname{proj}_C(\vy)$ and $\va=\vy-\vx^*$, which is
nonzero because $\vy\notin C$. Take the hyperplane through the midpoint of
$\vx^*$ and $\vy$, $b=\va^\T(\vx^*+\vy)/2$. Then
$\va^\T\vy=b+\tfrac12\norm{\va}^2>b$, and $\va^\T\vx^*=b-\tfrac12\norm{\va}^2$.
For $\vx\in C$, {eq}`eq-02-01-projection` gives $\va^\T(\vx-\vx^*)\le0$, so
$\va^\T\vx\le\va^\T\vx^*=b-\tfrac12\norm{\va}^2<b$.
:::

{numref}`fig-02-01-separation`(a) shows the construction. A consequence is that
a closed convex set is the intersection of all closed half-spaces that
contain it: every point outside the set is cut off by one of them. This is
the "outer" description of convex sets, of which polyhedra are the case
with finitely many half-spaces.

```{figure} figures/02-01-separation.svg
:name: fig-02-01-separation
:width: 640px
:class: dark-light
:alt: Left, a disc C, a point y outside it, the closest point x-star on the circle, and a line between them perpendicular to y minus x-star. Right, a convex polygon with a line touching it at a boundary point x0 and the normal vector a pointing away from the set

(a) Separating hyperplane for the disc $C=\bar B((1,1),1)$ and the point
$\vy=(3,2)$: the normal is $\va=\vy-\operatorname{proj}_C(\vy)$ and the
hyperplane passes through the midpoint of the dashed segment. (b) A
supporting hyperplane at the boundary point $\vx_0$: the whole set lies in
the half-space $\va^\T\vx\le\va^\T\vx_0$.
```

Two disjoint convex sets can also be separated.

:::{prf:theorem} Separating hyperplane theorem
:label: thm-02-01-separating-hyperplane

If $C$ and $D$ are nonempty, convex and disjoint, there are $\va\ne\vzero$
and $b$ with $\va^\T\vx\le b$ for all $\vx\in C$ and $\va^\T\vx\ge b$ for all
$\vx\in D$. If, in addition, $C$ is closed and $D$ is compact, the
inequalities can be made strict.
:::

The strict case follows from {prf:ref}`thm-02-01-separation-point`: the set
$C-D=\{\vx-\vy:\vx\in C,\vy\in D\}$ is convex
({prf:ref}`thm-02-01-affine-maps`), closed (because $D$ is compact) and does
not contain $\vzero$, so $\vzero$ can be separated from it. The general
case needs a limiting argument (Boyd & Vandenberghe, 2004, §2.5.1). Strict
separation can fail for closed sets that are not compact: the closed convex
sets $\{\vx:x_2\le0\}$ and $\{\vx:x_1>0,\ x_1x_2\ge1\}$ are disjoint, but
they come arbitrarily close, and the only separating hyperplane is
$x_2=0$, which touches neither.

:::{prf:theorem} Supporting hyperplane theorem
:label: thm-02-01-supporting-hyperplane

Let $C$ be convex and $\vx_0\in\operatorname{bd}C$. Then there is
$\va\ne\vzero$ with $\va^\T\vx\le\va^\T\vx_0$ for all $\vx\in C$. The
hyperplane $\{\vx:\va^\T\vx=\va^\T\vx_0\}$ is a **supporting hyperplane** of
$C$ at $\vx_0$.
:::

:::{dropdown} Proof sketch
:class-container: proof
:animate: fade-in

Because $\vx_0$ is a boundary point of the convex set $C$, there are points
$\vy^{(k)}\notin\operatorname{cl}C$ that converge to $\vx_0$. Separate each
from $\operatorname{cl}C$ by {prf:ref}`thm-02-01-separation-point` with a
normal $\va^{(k)}$ of length 1. The bounded sequence $\{\va^{(k)}\}$ has a
limit point $\va$ (Bolzano–Weierstrass), $\norm{\va}=1$, and passing to the
limit along the corresponding subsequence in
$\va^{(k)\T}\vx<\va^{(k)\T}\vy^{(k)}$ gives $\va^\T\vx\le\va^\T\vx_0$ for all
$\vx\in C$.
:::

Geometrically, a supporting hyperplane touches $C$ at $\vx_0$ and has $C$
entirely on one side ({numref}`fig-02-01-separation`(b)). Where the boundary is
smooth, it is the tangent plane; at a corner there are many supporting
hyperplanes. In optimization, the level set $\{\vx:\vc^\T\vx=\vc^\T\vx^*\}$
through a minimizer $\vx^*$ of $\vc^\T\vx$ over $C$ is a supporting
hyperplane of $C$ at $\vx^*$ with $\va=-\vc$, and the optimality condition
of {prf:ref}`thm-02-03-optimality` generalizes this picture to
nonlinear objectives.

## Summary and Common Mistakes

### Summary

- Open, closed, bounded and compact sets are defined through balls
  ({prf:ref}`def-02-01-open-closed`, {prf:ref}`def-02-01-bounded-compact`); a set
  is closed if and only if it contains the limits of its convergent
  sequences ({prf:ref}`thm-02-01-closed-sequences`). A sequence's limit points
  are the limits of its convergent subsequences
  ({prf:ref}`ex-02-01-limit-points`).
- A continuous function attains its minimum on a nonempty compact set
  ({prf:ref}`thm-02-01-weierstrass`), and on a closed set if it is coercive
  ({prf:ref}`cor-02-01-coercive`).
- Affine sets contain lines, convex sets contain segments
  ({prf:ref}`def-02-01-convex-set`), convex cones contain conic combinations
  ({prf:ref}`def-02-01-cone`); $\conv S$ is the smallest convex set
  containing $S$.
- Hyperplanes, half-spaces, polyhedra $\{\vx:\mA\vx\le\vb\}$, norm balls,
  ellipsoids, norm cones and the cone $\Sym^n_+$ are convex.
- Intersections and affine images and preimages preserve convexity
  ({prf:ref}`thm-02-01-intersection`, {prf:ref}`thm-02-01-affine-maps`); unions do
  not.
- The seven sets of {prf:ref}`ex-02-01-seven-sets` illustrate every property
  ({numref}`tab-02-01-seven-sets`).
- Extreme points ({prf:ref}`def-02-01-extreme-point`) generate compact convex
  sets, and a linear function attains its minimum over such a set at one of
  them.
- A point outside a closed convex set can be strictly separated from it
  ({prf:ref}`thm-02-01-separation-point`), using the projection
  ({prf:ref}`thm-02-01-projection`); every boundary point of a convex set has
  a supporting hyperplane ({prf:ref}`thm-02-01-supporting-hyperplane`).

### Common mistakes

- **Taking "closed" as "not open".** Most sets, such as $\Omega_4$ or
  $[0,1)$, are neither; $\R^n$ is both.
- **Forgetting that unions do not preserve convexity.** Intersections do;
  unions, such as $\Omega_3$, generally do not.
- **Calling every convex set with corners a polyhedron.** A polyhedron is an
  intersection of *finitely many closed* half-spaces; a disc is not one,
  and neither is $\Omega_4$, which is not closed.
- **Applying Weierstrass' theorem to unbounded or open sets.** A continuous
  function on a closed but unbounded set needs an additional argument such
  as coercivity ({prf:ref}`cor-02-01-coercive`).
- **Confusing boundary points with extreme points.** Every extreme point
  of a compact convex set is a boundary point, but points in the middle of
  an edge are boundary points that are not extreme.
- **Confusing limit point and limit.** A sequence can have several limit
  points ({prf:ref}`ex-02-01-limit-points`); it converges only if it has
  exactly one and is bounded.

## References

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 2: convex
  hulls, closure and interior, Weierstrass' theorem, separation and support
  of sets, convex cones, polyhedral sets and extreme points; Appendix A,
  mathematical review.)
- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Sections 2.1–2.2, polyhedra and extreme
  points.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Chapter 2,
  Sections 2.1–2.5: affine and convex sets, important examples, operations
  that preserve convexity, separating and supporting hyperplanes; Appendix
  A.2, topology.)
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Appendix B, convex sets.)
