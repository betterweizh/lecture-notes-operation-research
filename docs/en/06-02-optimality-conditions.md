# Optimality Conditions

## Overview

How can one recognize a minimum of a nonlinear program without comparing
it with every other feasible point? This chapter answers the question with
derivatives. Two kinds of directions organize the theory: **feasible
directions**, along which one can move a little without leaving $\Omega$,
and **descent directions**, along which $f$ decreases. At a local minimum
no direction can be both ({prf:ref}`thm-06-02-first-order`). From this
principle follow the first-order necessary conditions (a zero gradient at an
interior point, a variational inequality on a convex set), the second-order
necessary conditions ({prf:ref}`thm-06-02-second-order-necessary`) and the
second-order sufficient condition for unconstrained problems
({prf:ref}`thm-06-02-second-order-sufficient`). For convex programs the
first-order condition is also sufficient for a global minimum
({prf:ref}`thm-06-02-convex-sufficient`). The conditions here are stated in
terms of directions; {doc}`06-03-kkt-conditions` translates them into
conditions on the constraint functions, the Karush–Kuhn–Tucker (KKT)
conditions.
The chapter builds on the nonlinear programming (NLP) model and local
minima of
{doc}`06-01-nlp-models`, on gradients, Hessians and Taylor's theorem from
{doc}`01-02-matrix-calculus`, and on convex sets and functions from
{doc}`02-01-convex-sets` and {doc}`02-02-convex-functions`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- decide whether a direction is feasible at a point of a polyhedral or
  convex set, and whether it is a descent direction of $f$;
- state and prove the first-order necessary condition
  $\nabla f(\vx^*)^\T\vd\ge0$ for all feasible directions, and derive from
  it $\nabla f(\vx^*)=\vzero$ at interior points;
- state and prove the second-order necessary conditions and the
  second-order sufficient condition for unconstrained minima;
- classify stationary points as local minima, local maxima or saddle
  points with the Hessian; and
- show that for a convex program a point satisfying the first-order
  condition is a global minimum.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\vd$ | a direction in $\R^n$; $d_j$ its components |
| $\alpha$, $\bar\alpha$ | step length along $\vd$; a step length up to which $\vx+\alpha\vd$ stays feasible |
| $D(\vx)$ | set of feasible directions of $\Omega$ at $\vx$ |
| $\ve_j$ | $j$-th unit vector of $\R^n$ |
| $\vartheta\in(0,1)$ | position of the intermediate point in Taylor's theorem |
| $\sigma$ | smallest eigenvalue $\lambda_{\min}(\nabla^2f(\vx^*))$ in the proof of {prf:ref}`thm-06-02-second-order-sufficient` |

Throughout, $f$ is continuously differentiable on an open set containing
$\Omega$, and twice continuously differentiable whenever its Hessian
appears. The problem is $\min\{f(\vx):\vx\in\Omega\}$ with
$\Omega\subseteq\R^n$; the description of $\Omega$ by constraint functions
is not needed until {doc}`06-03-kkt-conditions`.

(sec-06-02-feasible-directions)=
## Feasible directions

From a feasible point one can move only in certain directions without
leaving $\Omega$. The following definition makes this precise; it requires
the whole segment from $\vx$ to $\vx+\bar\alpha\vd$ to be feasible, not only
its end point.

:::{prf:definition} Feasible direction
:label: def-06-02-feasible-direction

Let $\vx\in\Omega$. A vector $\vd\in\R^n$ is a **feasible direction** of
$\Omega$ at $\vx$ if there is an $\bar\alpha>0$ such that

$$
\vx+\alpha\vd\in\Omega\qquad\text{for all }0\le\alpha\le\bar\alpha .
$$ (eq-06-02-feasible-direction)

The set of feasible directions at $\vx$ is denoted $D(\vx)$.
:::

If $\vd\in D(\vx)$ then $t\vd\in D(\vx)$ for every $t>0$ (use
$\bar\alpha/t$), so $D(\vx)$ is a cone, and only the direction of $\vd$
matters, not its length. Three situations cover most uses.

:::{prf:proposition} Feasible directions in three common cases
:label: prop-06-02-feasible-directions

1. If $\vx$ is an interior point of $\Omega$, every $\vd\in\R^n$ is a
   feasible direction at $\vx$.
2. If $\Omega$ is convex and $\vx,\vy\in\Omega$, then $\vd=\vy-\vx$ is a
   feasible direction at $\vx$ (with $\bar\alpha=1$).
3. If $\Omega=\{\vx:\va_i^\T\vx\le b_i,\ i=1,\dots,m\}$ is a polyhedron,
   then $\vd$ is a feasible direction at $\vx\in\Omega$ if and only if
   $\va_i^\T\vd\le0$ for every $i$ with $\va_i^\T\vx=b_i$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

We treat the three cases in turn.

1. Some ball $B(\vx,\varepsilon)$ lies in $\Omega$. For $\vd\ne\vzero$ take
   $\bar\alpha=\varepsilon/(2\norm{\vd})$; then
   $\norm{\alpha\vd}\le\varepsilon/2<\varepsilon$ for $0\le\alpha\le\bar\alpha$,
   so $\vx+\alpha\vd\in B(\vx,\varepsilon)\subseteq\Omega$.
2. $\vx+\alpha(\vy-\vx)=(1-\alpha)\vx+\alpha\vy$ is a convex combination of
   two points of $\Omega$ for $0\le\alpha\le1$.
3. For a constraint with $\va_i^\T\vx=b_i$ (an *active* constraint),
   $\va_i^\T(\vx+\alpha\vd)=b_i+\alpha\va_i^\T\vd$, which is $\le b_i$ for
   some $\alpha>0$ exactly when $\va_i^\T\vd\le0$. A constraint with
   $\va_i^\T\vx<b_i$ has slack $s_i=b_i-\va_i^\T\vx>0$ and stays satisfied
   for $\alpha\le s_i/\va_i^\T\vd$ if $\va_i^\T\vd>0$, and for all $\alpha\ge0$
   otherwise. Taking $\bar\alpha$ as the smallest of these finitely many
   bounds proves the claim.
:::

Part 3 says that at a point of a polyhedron only the active constraints
restrict the directions, and they do so through a linear test. The same
test for nonlinear constraints is the subject of {doc}`06-03-kkt-conditions`,
where it turns out to be more delicate.

:::{prf:remark} Feasible directions on curved sets
:label: rem-06-02-curved

On a curved boundary there may be far fewer feasible directions than one
expects. At the point $(1,0)$ of the unit disc
$\{\vx:x_1^2+x_2^2\le1\}$, the direction $\vd=(d_1,d_2)$ is feasible
exactly when $d_1<0$ or $\vd=\vzero$: moving along the tangent $(0,1)$
leaves the disc at once, since $1+\alpha^2>1$. On the unit *circle*
$\{\vx:x_1^2+x_2^2=1\}$, a set defined by an equation, the only feasible
direction at any point is $\vd=\vzero$, because no segment lies on a
circle. Then {prf:ref}`thm-06-02-first-order` holds trivially and says
nothing. For such sets one must move along feasible *curves* rather than
straight lines; {doc}`06-03-kkt-conditions` does this with the tangent cone
and obtains conditions that are informative for curved constraints as
well. The results of the present chapter are sharp for interior points,
for polyhedra and, through {prf:ref}`prop-06-02-feasible-directions`, for
convex sets.
:::

(sec-06-02-descent-directions)=
## Descent directions

:::{prf:definition} Descent direction
:label: def-06-02-descent-direction

A vector $\vd$ is a **descent direction** of $f$ at $\vx$ if
$\nabla f(\vx)^\T\vd<0$. A **feasible descent direction** at $\vx\in\Omega$
is a direction that is both feasible and a descent direction.
:::

The number $\nabla f(\vx)^\T\vd$ is the directional derivative of $f$ at
$\vx$ in direction $\vd$ ({prf:ref}`def-01-02-directional-derivative`). Its sign has a
geometric meaning: $\nabla f(\vx)^\T\vd<0$ exactly when the angle between
$\vd$ and the gradient is obtuse, that is, when $\vd$ points into the side of
the tangent plane of the level set $\{\vy:f(\vy)=f(\vx)\}$ on which $f$ is
smaller. The name is justified by the following lemma.

:::{prf:lemma} Descent directions decrease $f$
:label: lem-06-02-descent

If $\nabla f(\vx)^\T\vd<0$, there is an $\bar\alpha>0$ such that
$f(\vx+\alpha\vd)<f(\vx)$ for all $0<\alpha\le\bar\alpha$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\varphi(\alpha)=f(\vx+\alpha\vd)$. By the chain rule
$\varphi'(0)=\nabla f(\vx)^\T\vd<0$, and by the definition of the derivative

$$
\lim_{\alpha\downarrow0}\frac{f(\vx+\alpha\vd)-f(\vx)}{\alpha}=\varphi'(0)<0 .
$$

Hence the quotient is negative for all sufficiently small $\alpha>0$, which
is the claim.
:::

The converse is false: $f$ may decrease along $\vd$ although
$\nabla f(\vx)^\T\vd=0$ (for $f(x)=-x^2$ at $x=0$, every $d\ne0$ decreases
$f$). Directions with $\nabla f(\vx)^\T\vd=0$ are undecided at first order,
and second derivatives are needed to judge them
({ref}`sec-06-02-second-order`). Among directions of unit length, the
Cauchy–Schwarz inequality $\nabla f(\vx)^\T\vd\ge-\norm{\nabla f(\vx)}$
shows that $\vd=-\nabla f(\vx)/\norm{\nabla f(\vx)}$ gives the most negative
directional derivative; this **steepest descent direction** is the basis of
the method in {doc}`06-05-unconstrained-algorithms`.

:::{prf:example} Feasible and descent directions on a triangle
:label: ex-06-02-directions

Let $\Omega=\{\vx\in\R^2:x_1+x_2\le2,\ x_1\ge0,\ x_2\ge0\}$, the triangle
with vertices $(0,0)$, $(2,0)$, $(0,2)$, and $f(\vx)=(x_1-1)^2+(x_2-2)^2$.

(a) Which directions are feasible at the point $(1,0.5)$ and at the vertex
$(2,0)$? Are $\vd=(-1,0.5)$ and $\vd=(1,1.5)$ feasible at $(2,0)$?

(b) Classify the directions $\vd_1=(-0.4,1)$, $\vd_2=(0,1)$,
$\vd_3=(0.45,1)$, $\vd_4=(-1,-0.8)$ at $(1,0.5)$ and
$\vd_5=(-1,0.5)$, $\vd_6=(0.6,1)$, $\vd_7=(1,0.13)$ at $(2,0)$ as descent or
not, and as feasible or not.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in


**(a)** The point $(1,0.5)$ satisfies all three constraints strictly. Its
distances to the three boundary lines are $1$ (to $x_1=0$), $0.5$ (to
$x_2=0$) and $\abs{1+0.5-2}/\sqrt2\approx0.354$ (to $x_1+x_2=2$), so the
ball of radius $0.354$ around it lies in $\Omega$. The point is interior,
and by {prf:ref}`prop-06-02-feasible-directions` every direction is feasible
there: whatever the direction, a step shorter than $0.354/\norm{\vd}$ stays
in the ball.

At $(2,0)$ the constraints $x_1+x_2\le2$ and $x_2\ge0$ (that is,
$-x_2\le0$) are active and $x_1\ge0$ is not. By part 3 of
{prf:ref}`prop-06-02-feasible-directions` the feasible directions are exactly

$$
D(2,0)=\{\vd:\ d_1+d_2\le0,\ -d_2\le0\}.
$$

For $\vd=(-1,0.5)$: $d_1+d_2=-0.5\le0$ and $d_2=0.5\ge0$, so it is feasible;
indeed $(2,0)+\alpha(-1,0.5)\in\Omega$ for $0\le\alpha\le2$, and $\alpha=1$
reaches $(1,0.5)$. For $\vd=(1,1.5)$: $d_1+d_2=2.5>0$, so every step
leaves $\Omega$ through the edge $x_1+x_2=2$; it is not feasible.

**(b)** The gradient is $\nabla f(\vx)=(2(x_1-1),\,2(x_2-2))^\T$, so
$\nabla f(1,0.5)=(0,-3)^\T$ and $\nabla f(2,0)=(2,-4)^\T$.

At $(1,0.5)$, $\nabla f^\T\vd=-3d_2$: a direction is a descent direction
exactly when it points upwards, $d_2>0$. Thus $\vd_1,\vd_2,\vd_3$ (each with
$\nabla f^\T\vd=-3$) are descent directions, and $\vd_4$ (with
$\nabla f^\T\vd_4=2.4$) is an ascent direction. All four are feasible since
the point is interior.

At $(2,0)$:

```{list-table} Directions at the vertex $(2,0)$.
:name: tab-06-02-directions-vertex
:header-rows: 1

* - $\vd$
  - $\nabla f(2,0)^\T\vd$
  - $d_1+d_2\le0$, $d_2\ge0$?
  - conclusion
* - $\vd_5=(-1,0.5)$
  - $-2-2=-4$
  - yes
  - feasible descent direction
* - $\vd_6=(0.6,1)$
  - $1.2-4=-2.8$
  - no ($d_1+d_2=1.6$)
  - descent, but infeasible
* - $\vd_7=(1,0.13)$
  - $2-0.52=1.48$
  - no ($d_1+d_2=1.13$)
  - ascent and infeasible
```

Only $\vd_5$ is a feasible descent direction. Since one exists, $(2,0)$ is
not a local minimum: by {prf:ref}`lem-06-02-descent` small steps along $\vd_5$
stay feasible and decrease $f$. The same is true at $(1,0.5)$. Where then is
the minimum? It is the point of $\Omega$ closest to $(1,2)$, namely
$\vx^*=(0.5,1.5)$ with $f(\vx^*)=0.5$; {ref}`sec-06-02-convex` proves this.
{numref}`fig-06-02-directions` shows both points with the contours of $f$.

```{figure} figures/06-02-directions.svg
:name: fig-06-02-directions
:width: 680px
:class: dark-light
:alt: Two panels with the triangle x1 + x2 <= 2 and circular contours of f centered at (1,2); left, directions d1 to d4 and the gradient at (1,0.5); right, directions d5 to d7 and the gradient at the vertex (2,0)

Directions of {prf:ref}`ex-06-02-directions`. Left: at the interior point
$(1,0.5)$ every direction is feasible, and $\vd_1,\vd_2,\vd_3$ make an
obtuse angle with $\nabla f=(0,-3)$ (they point into the inner side of the
contour $f=2.25$). Right: at the vertex $(2,0)$ only $\vd_5$ points into
$\Omega$ and makes an obtuse angle with $\nabla f=(2,-4)$; $\vd_6$ is a
descent direction that leaves $\Omega$.
```
:::

(sec-06-02-first-order)=
## First-order necessary conditions

A local minimum admits no feasible descent direction. Stated with the
directional derivative, this is the basic first-order condition.

:::{prf:theorem} First-order necessary condition
:label: thm-06-02-first-order

If $\vx^*$ is a local minimum of $f$ over $\Omega$, then

$$
\nabla f(\vx^*)^\T\vd\ge0\qquad\text{for every feasible direction }\vd\in D(\vx^*).
$$ (eq-06-02-first-order)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Suppose some $\vd\in D(\vx^*)$ has $\nabla f(\vx^*)^\T\vd<0$. Let
$\varepsilon>0$ be the radius in the definition of a local minimum
({prf:ref}`def-06-01-local-global`). By {eq}`eq-06-02-feasible-direction` and
{prf:ref}`lem-06-02-descent` there are $\bar\alpha_1,\bar\alpha_2>0$ with
$\vx^*+\alpha\vd\in\Omega$ for $\alpha\le\bar\alpha_1$ and
$f(\vx^*+\alpha\vd)<f(\vx^*)$ for $0<\alpha\le\bar\alpha_2$. Choose
$0<\alpha\le\min\{\bar\alpha_1,\bar\alpha_2\}$ with
$\alpha\norm{\vd}<\varepsilon$. Then $\vx^*+\alpha\vd$ is a feasible point
within distance $\varepsilon$ of $\vx^*$ with a smaller objective value,
contradicting local minimality.
:::

Geometrically, {eq}`eq-06-02-first-order` says that at a local minimum every
feasible direction makes an angle of at most $90^\circ$ with the gradient:
$-\nabla f(\vx^*)$, the direction of fastest decrease, points out of
$\Omega$ or is orthogonal to all feasible directions.

At an interior point all directions are feasible, and the condition
becomes the familiar one.

:::{prf:corollary} Stationarity at interior minima
:label: cor-06-02-stationary

If $\vx^*$ is a local minimum of $f$ over $\Omega$ and an interior point of
$\Omega$ (in particular, if $\Omega=\R^n$), then

$$
\nabla f(\vx^*)=\vzero .
$$ (eq-06-02-stationary)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

By {prf:ref}`prop-06-02-feasible-directions` both $\ve_j$ and $-\ve_j$ are
feasible directions, and {eq}`eq-06-02-first-order` gives
$\nabla f(\vx^*)^\T\ve_j=\partial f(\vx^*)/\partial x_j\ge0$ and
$\nabla f(\vx^*)^\T(-\ve_j)=-\partial f(\vx^*)/\partial x_j\ge0$. Hence
$\partial f(\vx^*)/\partial x_j=0$ for $j=1,\dots,n$.
:::

A point with $\nabla f(\vx)=\vzero$ is a **stationary point** (or critical
point) of $f$. A stationary point may be a local minimum, a local maximum or
neither; a stationary point that is neither is a **saddle point**.
{eq}`eq-06-02-stationary` is a system of $n$ equations in $n$ unknowns, and
solving it is the first step in finding unconstrained minima. In problem 1
of {prf:ref}`ex-06-01-three-feasible-sets`, $\nabla f=(2(x_1-2),2(x_2-2))$
vanishes only at $(2,2)$, the minimum found there.

When $\Omega$ is convex, every $\vy-\vx^*$ with $\vy\in\Omega$ is a
feasible direction, and {eq}`eq-06-02-first-order` implies the
**variational inequality**

$$
\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0\qquad\text{for all }\vy\in\Omega .
$$ (eq-06-02-variational-inequality)

Two special cases show how it works.

- *An interval.* For $\min\{f(x):a\le x\le b\}$ in one variable,
  {eq}`eq-06-02-variational-inequality` says: $f'(x^*)=0$ if $a<x^*<b$;
  $f'(a)\ge0$ if $x^*=a$; $f'(b)\le0$ if $x^*=b$. The function $(x^2)^x$ of
  {doc}`06-01-nlp-models` has its global minimum over $[-5,1.7]$ at the left
  end point with $f'(-5)\approx5.3\times10^{-7}>0$, in agreement with the
  condition.
- *The nonnegative orthant.* For $\Omega=\{\vx:\vx\ge\vzero\}$, the
  directions $\ve_j$ are always feasible, and $-\ve_j$ is feasible when
  $x^*_j>0$. Hence a local minimum satisfies

  $$
  \begin{aligned}
  &\frac{\partial f(\vx^*)}{\partial x_j}\ge0, && j=1,\dots,n,\\
  &x_j^*\ge0, && j=1,\dots,n,\\
  &x_j^*\,\frac{\partial f(\vx^*)}{\partial x_j}=0, && j=1,\dots,n .
  \end{aligned}
  $$ (eq-06-02-orthant)

  Each variable is either at its bound with a nonnegative partial
  derivative, or positive with a zero partial derivative. This
  *complementarity* reappears in the KKT conditions of
  {doc}`06-03-kkt-conditions`.

(sec-06-02-second-order)=
## Second-order conditions

A zero gradient does not make a minimum. The function $f(\vx)=x_1^2-x_2^2$
has $\nabla f(\vzero)=\vzero$, but it increases along the $x_1$-axis and
decreases along the $x_2$-axis: the origin is a saddle point, shaped like a
mountain pass. First derivatives cannot tell a saddle point from a minimum,
because they only describe the tangent plane, which is horizontal in both
cases. The curvature, described by the Hessian, can.

:::{prf:theorem} Second-order necessary conditions
:label: thm-06-02-second-order-necessary

Let $f$ be twice continuously differentiable and let $\vx^*$ be a local
minimum of $f$ over $\Omega$. Then for every feasible direction
$\vd\in D(\vx^*)$:

1. $\nabla f(\vx^*)^\T\vd\ge0$; and
2. if $\nabla f(\vx^*)^\T\vd=0$, then $\vd^\T\nabla^2f(\vx^*)\vd\ge0$.

In particular, if $\vx^*$ is an interior point of $\Omega$ (or
$\Omega=\R^n$), then $\nabla f(\vx^*)=\vzero$ and $\nabla^2f(\vx^*)$ is
positive semidefinite.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Part 1 is {prf:ref}`thm-06-02-first-order`. For part 2, let $\vd\in D(\vx^*)$
with $\nabla f(\vx^*)^\T\vd=0$ and suppose
$\vd^\T\nabla^2f(\vx^*)\vd<0$. Since $\nabla^2f$ is continuous, there is
$\alpha_1>0$ with $\vd^\T\nabla^2f(\vx^*+\alpha\vd)\vd<0$ for
$0\le\alpha\le\alpha_1$. Taylor's theorem with the Lagrange remainder
({prf:ref}`thm-01-02-taylor`) gives, for each such $\alpha$, a
$\vartheta\in(0,1)$ with

$$
\begin{aligned}
f(\vx^*+\alpha\vd)
&=f(\vx^*)+\alpha\nabla f(\vx^*)^\T\vd+\tfrac12\alpha^2\vd^\T\nabla^2f(\vx^*+\vartheta\alpha\vd)\vd\\
&=f(\vx^*)+\tfrac12\alpha^2\vd^\T\nabla^2f(\vx^*+\vartheta\alpha\vd)\vd ,
\end{aligned}
$$

where the first-order term vanishes because $\nabla f(\vx^*)^\T\vd=0$. Since
$\vartheta\alpha\le\alpha_1$, the last term is negative, so
$f(\vx^*+\alpha\vd)<f(\vx^*)$ for all small $\alpha>0$. These points are
feasible (for $\alpha\le\bar\alpha$) and arbitrarily close to $\vx^*$, a
contradiction.

At an interior point, $\nabla f(\vx^*)=\vzero$ by
{prf:ref}`cor-06-02-stationary`, every $\vd$ is feasible and satisfies
$\nabla f(\vx^*)^\T\vd=0$, and part 2 says $\vd^\T\nabla^2f(\vx^*)\vd\ge0$
for all $\vd\in\R^n$, which is positive semidefiniteness.
:::

Positive semidefiniteness of a symmetric matrix is equivalent to
nonnegative eigenvalues $\lambda_i(\nabla^2f(\vx^*))\ge0$
({prf:ref}`thm-01-01-eigen-test`). A stationary point at which the Hessian has a
negative eigenvalue is therefore not a local minimum; if it also has a
positive eigenvalue it is not a local maximum either, and it is a saddle
point.

Strengthening semidefiniteness to definiteness gives a sufficient
condition.

:::{prf:theorem} Second-order sufficient condition (unconstrained)
:label: thm-06-02-second-order-sufficient

Let $f$ be twice continuously differentiable on $\R^n$. If

$$
\begin{aligned}
\nabla f(\vx^*)&=\vzero,\\
\nabla^2f(\vx^*)&\succ\vzero ,
\end{aligned}
$$ (eq-06-02-sosc)

then $\vx^*$ is a strict local minimum of $f$. The same holds for
$\min_{\vx\in\Omega}f$ when $\vx^*$ is an interior point of $\Omega$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $\sigma=\lambda_{\min}(\nabla^2f(\vx^*))>0$, so that
$\vd^\T\nabla^2f(\vx^*)\vd\ge\sigma\norm{\vd}^2$ for all $\vd$. By
continuity of $\nabla^2f$ there is $\varepsilon>0$ such that the spectral
norm $\norm{\nabla^2f(\vy)-\nabla^2f(\vx^*)}\le\sigma/2$ for all
$\vy\in B(\vx^*,\varepsilon)$; for such $\vy$,

$$
\begin{aligned}
\vd^\T\nabla^2f(\vy)\vd
&\ge\vd^\T\nabla^2f(\vx^*)\vd-\norm{\nabla^2f(\vy)-\nabla^2f(\vx^*)}\,\norm{\vd}^2\\
&\ge\tfrac{\sigma}{2}\norm{\vd}^2 .
\end{aligned}
$$

Now let $\vx\in B(\vx^*,\varepsilon)$, $\vx\ne\vx^*$, and $\vd=\vx-\vx^*$.
By Taylor's theorem there is $\vartheta\in(0,1)$ with

$$
\begin{aligned}
f(\vx)&=f(\vx^*)+\nabla f(\vx^*)^\T\vd+\tfrac12\vd^\T\nabla^2f(\vx^*+\vartheta\vd)\vd\\
&\ge f(\vx^*)+\tfrac{\sigma}{4}\norm{\vx-\vx^*}^2>f(\vx^*),
\end{aligned}
$$

because $\nabla f(\vx^*)=\vzero$ and $\vx^*+\vartheta\vd\in
B(\vx^*,\varepsilon)$.
:::

There is a gap between the necessary condition (semidefinite) and the
sufficient one (definite). At a stationary point with a singular
semidefinite Hessian, second derivatives decide nothing: $x^4$, $-x^4$ and
$x^3$ all have $f'(0)=f''(0)=0$, and the origin is a minimum, a maximum and
neither, respectively. Conversely, the sufficient condition is not
necessary: $x^4$ has a strict minimum at $0$ with $f''(0)=0$.

For constrained problems, condition 2 of
{prf:ref}`thm-06-02-second-order-necessary` only involves feasible directions
along which $f$ is flat to first order; directions with
$\nabla f(\vx^*)^\T\vd>0$ already increase $f$ and need no curvature test.
When the constraints themselves are curved, feasible directions do not
capture the curvature of the boundary, and the correct second-order
conditions use the Hessian of the Lagrangian; they are given in
{doc}`06-03-kkt-conditions`.

:::{prf:example} A saddle point and its constrained version
:label: ex-06-02-saddle

Let $f(\vx)=x_1^2+x_1x_2-x_1-x_2$.

(a) Find the stationary points of $f$ on $\R^2$ and decide whether they are
local minima.

(b) Now minimize $f$ subject to $x_1\ge1$, $x_2\ge0$. Check the first- and
second-order necessary conditions at $\vx^*=(1,0)$ and decide whether it is
a minimum.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** The gradient and Hessian are

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}2x_1+x_2-1\\x_1-1\end{pmatrix},\\
\nabla^2f(\vx)&=\begin{pmatrix}2&1\\1&0\end{pmatrix}.
\end{aligned}
$$

Setting $\nabla f=\vzero$ gives $x_1=1$ and then $x_2=1-2x_1=-1$: the only
stationary point is $(1,-1)$, where $f(1,-1)=0$. The Hessian is constant
with determinant $2\cdot0-1\cdot1=-1<0$, so its eigenvalues have opposite
signs: $\lambda_{1,2}(\nabla^2f)=1\mp\sqrt2\approx-0.414,\ 2.414$. It is
indefinite, the second-order necessary condition fails, and $(1,-1)$ is not
a local minimum; by the same argument applied to $-f$ it is not a local
maximum. It is a saddle point. *(An earlier version of this example
evaluates the Hessian "at $(1,0)$"; the stationary point is $(1,-1)$, and
since the Hessian is constant the conclusion is the same.)*

The saddle can be seen directly. Factoring,

$$
\begin{aligned}
f(\vx)&=x_1(x_1+x_2)-(x_1+x_2)\\
&=(x_1-1)(x_1+x_2),
\end{aligned}
$$

so $f=0$ on the two lines $x_1=1$ and $x_1+x_2=0$, which cross at
$(1,-1)$, and $f$ changes sign across each of them. Along $\vd=(1,0)$,
$f(1+\alpha,-1)=\alpha^2>0$, while along $\vd=(1,-2)$,
$f(1+\alpha,-1-2\alpha)=\alpha\cdot(-\alpha)=-\alpha^2<0$. These agree with
the curvatures $\vd^\T\nabla^2f\vd=2$ and $-2$ (for a quadratic,
$f(\vx^*+\alpha\vd)=f(\vx^*)+\tfrac12\alpha^2\vd^\T\nabla^2f\vd$ exactly at a
stationary point).

**(b)** At $\vx^*=(1,0)$ both constraints $x_1\ge1$ and $x_2\ge0$ are
active, and since they are linear, {prf:ref}`prop-06-02-feasible-directions`
gives the feasible directions $\vd=(d_1,d_2)$ with $d_1\ge0$, $d_2\ge0$.
The gradient is $\nabla f(1,0)=(2+0-1,\,1-1)^\T=(1,0)^\T\ne\vzero$.

- First order: $\nabla f(\vx^*)^\T\vd=d_1\ge0$ for every feasible $\vd$.
  The condition holds.
- Second order: $\nabla f(\vx^*)^\T\vd=0$ means $d_1=0$, so
  $\vd=(0,d_2)$ with $d_2\ge0$, and
  $\vd^\T\nabla^2f\vd=2d_1^2+2d_1d_2=0\ge0$. The condition holds.

Both necessary conditions hold, but they do not prove optimality (the
Hessian is indefinite, and {prf:ref}`thm-06-02-second-order-sufficient` does
not apply at a boundary point anyway). The factorization settles the
question: on $\Omega$, $x_1-1\ge0$ and $x_1+x_2\ge1>0$, so
$f(\vx)\ge0=f(1,0)$. Hence $(1,0)$ is a global minimum. It is not the only
one: $f(1,x_2)=0$ for every $x_2\ge0$, so the whole ray
$\{(1,x_2):x_2\ge0\}$ consists of global minima. This explains why the
curvature along the feasible direction $(0,1)$ is exactly zero: $f$ is
constant along it.

```{figure} figures/06-02-saddle.svg
:name: fig-06-02-saddle
:width: 520px
:class: dark-light
:alt: Contour plot of f = (x1 - 1)(x1 + x2) with hyperbolic level curves around the saddle point (1,-1), the zero level lines x1 = 1 and x1 + x2 = 0, and the shaded feasible quadrant x1 >= 1, x2 >= 0 whose left edge is a ray of minima

Contours of $f(\vx)=(x_1-1)(x_1+x_2)$ in {prf:ref}`ex-06-02-saddle`. The zero
level set consists of the two dashed lines, which cross at the saddle point
$(1,-1)$; $f>0$ in the left and right sectors and $f<0$ in the upper and
lower ones. On the shaded feasible set $x_1\ge1$, $x_2\ge0$ the minimum
value $0$ is attained along the whole edge $x_1=1$ (thick line).
```
:::

(sec-06-02-procedure)=
## Using the conditions

The theorems combine into a routine for problems whose minima are interior
points, and with small changes for simple bounds.

:::{prf:algorithm} Finding and classifying unconstrained minima
:label: alg-06-02-classify

**Input** a twice continuously differentiable $f$ on an open set
$U\subseteq\R^n$ (often $U=\R^n$).

**Output** the local minima of $f$ in $U$ and, when possible, the global
minimum.

1. **Stationary points.** Solve $\nabla f(\vx)=\vzero$ in $U$
   ({prf:ref}`cor-06-02-stationary`). Every local minimum is among the
   solutions.
2. **Classify.** At each stationary point compute $\nabla^2f$ and its
   eigenvalues (or leading principal minors):
   positive definite $\Rightarrow$ strict local minimum
   ({prf:ref}`thm-06-02-second-order-sufficient`); negative definite
   $\Rightarrow$ strict local maximum; indefinite $\Rightarrow$ saddle point
   ({prf:ref}`thm-06-02-second-order-necessary` applied to $f$ and to $-f$);
   singular semidefinite $\Rightarrow$ undecided, examine $f$ directly.
3. **Global minimum.** If $f$ is convex, any stationary point is a global
   minimum ({prf:ref}`thm-06-02-convex-sufficient`). Otherwise compare the
   values at the local minima, and make sure the minimum exists at all,
   for instance by coercivity ({prf:ref}`thm-06-01-existence`); on a bounded
   domain also compare with the boundary.
:::

Step 2 rests on the second-order conditions for maxima, which follow by
applying the theorems to $-f$: at a local maximum $\nabla f=\vzero$ and
$\nabla^2f\preceq\vzero$, and $\nabla f=\vzero$ with $\nabla^2f\prec\vzero$
gives a strict local maximum. Step 3 cannot be skipped: a function can have
a single stationary point that is a strict local minimum and still be
unbounded below (for example $x_1^2+x_2^2(1-x_1)^3$, a classical
counterexample in two variables), so a local test never certifies a global
minimum without further information.

For the function $f(x)=(x^2)^x$ of {numref}`fig-06-01-x-power-2x` on the
interval $[-5,1.7]$ the routine runs as follows. Differentiating
$f'(x)=2f(x)(\ln\abs{x}+1)$ once more gives

$$
f''(x)=4f(x)\bigl(\ln\abs{x}+1\bigr)^2+\frac{2f(x)}{x},\qquad x\ne0 .
$$

At the stationary points $x=\pm1/e$ the first term vanishes, so
$f''(1/e)=2e\cdot e^{-2/e}\approx2.605>0$ (a strict local minimum) and
$f''(-1/e)=-2e\cdot e^{2/e}\approx-11.35<0$ (a strict local maximum). The
point $x=0$, where $f$ is not differentiable, and the end points $-5$ and
$1.7$, where $f'$ need not vanish, are examined separately: the end point
conditions of {eq}`eq-06-02-variational-inequality` hold at $x=-5$
($f'(-5)>0$), and comparing values, $f(-5)\approx1.0\times10^{-7}$ is
smaller than $f(1/e)\approx0.479$. The local minimum at $1/e$ is therefore
not global, as found in {doc}`06-01-nlp-models`.

For problems with linear constraints the same scheme works at boundary
points with the conditions of {prf:ref}`thm-06-02-second-order-necessary`:
check $\nabla f(\vx^*)^\T\vd\ge0$ on the cone of feasible directions, then
the curvature on its part where $\nabla f(\vx^*)^\T\vd=0$. If that
curvature is *positive* for every such nonzero $\vd$, the point is a strict
local minimum; this is the linear-constraint case of the second-order
sufficient condition proved in {doc}`06-03-kkt-conditions`
({prf:ref}`thm-06-03-second-order`). In part (b) of
{prf:ref}`ex-06-02-saddle` the curvature along the critical direction
$(0,1)$ is zero, so this test is inconclusive, which is consistent with the
ray of non-strict minima found there.

Numerical algorithms use the same conditions as stopping rules. A method
for unconstrained minimization stops when $\norm{\nabla f(\vx^{(k)})}$ falls
below a tolerance, and a method for constrained problems when the KKT
conditions of {doc}`06-03-kkt-conditions` hold approximately. Such a point is
only an approximate stationary point; whether it is a minimum is decided by
the second-order conditions or by convexity, exactly as above
({doc}`06-04-line-search`, {doc}`06-05-unconstrained-algorithms`).

(sec-06-02-convex)=
## Convex problems: the first-order condition is sufficient

For a convex program the first-order necessary condition is also
sufficient, and it certifies a *global* minimum.

:::{prf:theorem} Optimality in convex programs
:label: thm-06-02-convex-sufficient

Let $\Omega$ be a nonempty convex set and $f$ continuously differentiable
on an open set containing $\Omega$ and convex on $\Omega$. A point
$\vx^*\in\Omega$ is a global minimum of $f$ over $\Omega$ if and only if

$$
\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0\qquad\text{for all }\vy\in\Omega .
$$ (eq-06-02-convex-condition)

In particular, for $\Omega=\R^n$ (or $\vx^*$ interior), $\vx^*$ is a global
minimum if and only if $\nabla f(\vx^*)=\vzero$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

*Necessity.* A global minimum is a local minimum, and each $\vy-\vx^*$ is
a feasible direction ({prf:ref}`prop-06-02-feasible-directions`), so
{eq}`eq-06-02-convex-condition` follows from {prf:ref}`thm-06-02-first-order`.

*Sufficiency.* Suppose {eq}`eq-06-02-convex-condition` holds but some
$\vy\in\Omega$ has $f(\vy)<f(\vx^*)$. A differentiable convex function lies
above its tangent planes ({prf:ref}`thm-02-02-first-order`):

$$
f(\vy)\ge f(\vx^*)+\nabla f(\vx^*)^\T(\vy-\vx^*)\ge f(\vx^*),
$$

a contradiction. The case $\Omega=\R^n$ follows by taking
$\vy=\vx^*\pm\ve_j$, as in {prf:ref}`cor-06-02-stationary`.
:::

Together with the facts that a strictly convex function has at most one
minimum and that a twice differentiable function with positive
semidefinite Hessian on a convex set is convex there
({prf:ref}`thm-02-02-second-order`), the
theorem gives a complete recipe for convex problems: find a point that
satisfies the first-order condition; it is a global minimum, and it is the
only one if $f$ is strictly convex. For a quadratic
$f(\vx)=\vc^\T\vx+\tfrac12\vx^\T\mQ\vx$ with $\mQ\succ\vzero$, the unique
global minimum over $\R^n$ solves $\mQ\vx=-\vc$.

The objective of {prf:ref}`ex-06-02-directions` is convex and its feasible set
is a polyhedron, so the theorem applies. At $\vx^*=(0.5,1.5)$,
$\nabla f(\vx^*)=(2(0.5-1),\,2(1.5-2))^\T=(-1,-1)^\T$, and for every
$\vy\in\Omega$

$$
\begin{aligned}
\nabla f(\vx^*)^\T(\vy-\vx^*)&=-(y_1-0.5)-(y_2-1.5)\\
&=2-(y_1+y_2)\ge0,
\end{aligned}
$$

because $y_1+y_2\le2$ on $\Omega$. Hence $(0.5,1.5)$ is the global minimum,
with $f^*=0.25+0.25=0.5$, and it is unique because $f$ is strictly convex.
Geometrically, $-\nabla f(\vx^*)=(1,1)$ is perpendicular to the edge
$x_1+x_2=2$ and points out of $\Omega$; {doc}`06-03-kkt-conditions` turns this
observation into the KKT conditions.

## Summary of the conditions

```{list-table} Optimality conditions of this chapter. "Necessary" conditions hold at every local minimum; "sufficient" conditions guarantee one.
:name: tab-06-02-conditions
:header-rows: 1

* -
  - Unconstrained ($\Omega=\R^n$ or $\vx^*$ interior)
  - Constrained (feasible directions $\vd\in D(\vx^*)$)
* - First-order necessary
  - $\nabla f(\vx^*)=\vzero$
  - $\nabla f(\vx^*)^\T\vd\ge0$ for all $\vd\in D(\vx^*)$; KKT conditions
    under a constraint qualification ({doc}`06-03-kkt-conditions`)
* - Second-order necessary
  - $\nabla f(\vx^*)=\vzero$ and $\nabla^2f(\vx^*)\succeq\vzero$
  - in addition: $\nabla f(\vx^*)^\T\vd=0\Rightarrow\vd^\T\nabla^2f(\vx^*)\vd\ge0$;
    second-order KKT conditions ({doc}`06-03-kkt-conditions`)
* - Second-order sufficient
  - $\nabla f(\vx^*)=\vzero$ and $\nabla^2f(\vx^*)\succ\vzero$ (strict local
    minimum)
  - second-order KKT sufficient conditions ({doc}`06-03-kkt-conditions`)
* - Convex $f$ and $\Omega$
  - $\nabla f(\vx^*)=\vzero$ is necessary and sufficient for a global
    minimum
  - $\nabla f(\vx^*)^\T(\vy-\vx^*)\ge0$ for all $\vy\in\Omega$ is necessary
    and sufficient for a global minimum
```

## Summary and Common Mistakes

### Summary

- A feasible direction keeps a whole segment in $\Omega$
  ({prf:ref}`def-06-02-feasible-direction`); at interior points all directions
  are feasible, and at a point of a polyhedron exactly those with
  $\va_i^\T\vd\le0$ for the active constraints
  ({prf:ref}`prop-06-02-feasible-directions`).
- A descent direction has $\nabla f(\vx)^\T\vd<0$ and decreases $f$ for
  short steps ({prf:ref}`lem-06-02-descent`); $-\nabla f$ is the steepest one.
- First-order necessary condition: no feasible descent direction exists,
  $\nabla f(\vx^*)^\T\vd\ge0$ ({prf:ref}`thm-06-02-first-order`); at interior
  points $\nabla f(\vx^*)=\vzero$ ({prf:ref}`cor-06-02-stationary`); on convex
  sets the variational inequality {eq}`eq-06-02-variational-inequality`.
- Second-order necessary: nonnegative curvature along feasible directions
  with zero slope; positive semidefinite Hessian at interior minima
  ({prf:ref}`thm-06-02-second-order-necessary`).
- Second-order sufficient: zero gradient and positive definite Hessian give
  a strict local minimum ({prf:ref}`thm-06-02-second-order-sufficient`).
- {prf:ref}`alg-06-02-classify` combines the conditions: solve
  $\nabla f=\vzero$, classify by the Hessian, then settle global optimality
  by convexity, coercivity or comparison of values.
- For convex programs the first-order condition is necessary and
  sufficient for a global minimum ({prf:ref}`thm-06-02-convex-sufficient`).
- {numref}`tab-06-02-conditions` collects the conditions.

### Common mistakes

- **Stopping at $\nabla f(\vx^*)=\vzero$.** A stationary point can be a
  saddle point ({prf:ref}`ex-06-02-saddle`) or a maximum. Check the Hessian or
  convexity.
- **Requiring $\nabla f(\vx^*)=\vzero$ at a boundary point.** At a
  constrained minimum the gradient is usually not zero; only the
  directional derivatives along feasible directions must be nonnegative.
- **Concluding a minimum from a positive semidefinite Hessian.** Semidefinite
  is necessary, definite is sufficient; in between nothing can be
  concluded from second derivatives.
- **Testing definiteness of a non-symmetric matrix by its minors.** Hessians
  of twice continuously differentiable functions are symmetric; leading
  principal minors decide *definiteness* of symmetric matrices only
  ({prf:ref}`thm-01-01-sylvester`), and nonnegative leading minors do not
  imply semidefiniteness ({prf:ref}`thm-01-01-psd-minors`).
- **Treating "necessary" as "sufficient".** The point $(1,0)$ in
  {prf:ref}`ex-06-02-saddle` satisfies both necessary conditions; that alone
  does not make it a minimum. Its optimality needed a separate argument.


## References

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Section 3.4,
  minima of convex functions; Section 4.1, optimality conditions for
  unconstrained problems; Section 4.2, feasible and improving directions.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Section
  4.2.3, optimality criterion for differentiable objectives.)
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Sections 7.1–7.5, first- and second-order conditions, convex and concave
  functions, minimization of convex functions.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Section 2.1,
  recognizing a local minimum: Theorems 2.2–2.5.)
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis
  with mathematical programming methods*. Prentice-Hall. (Chapter 2, basic
  concepts in minimization problems.)
