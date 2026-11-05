# Line Search Methods

## Overview

The optimality conditions of {doc}`06-02-optimality-conditions` and
{doc}`06-03-kkt-conditions` tell us how to *recognize* a minimizer; they rarely
tell us how to *find* one, because the equations $\nabla f(\vx)=\vzero$ or the
Karush–Kuhn–Tucker (KKT) system can seldom be solved in closed form. Numerical methods replace the
single hard problem by a sequence of easy ones: from the current point they
choose a direction, move some distance along it, and repeat. This chapter
first makes the idea of an **iterative algorithm** precise (algorithmic maps,
descent functions, stopping rules and rates of convergence) and then treats
the step that every such method needs, the **line search**: minimizing a
function of one variable along a direction. Three classical methods are
developed and run on the same function, $f(x)=\sin x$ on $[3,6]$: **bisection**,
which halves an interval using the sign of $f'$; **golden-section search**,
which needs only function values and reuses one of them at every step thanks
to the ratio $\tau=(\sqrt5-1)/2$; and **Newton's method**, which fits a
quadratic and converges quadratically, but only when started close enough to
a minimizer. Quadratic interpolation and the inexact **Armijo** and **Wolfe**
conditions, which practical codes use inside the multidimensional methods of
{doc}`06-05-unconstrained-algorithms` and {doc}`06-06-constrained-algorithms`,
complete the chapter. Taylor expansions and directional derivatives from
{doc}`01-02-matrix-calculus` and convexity of functions from
{doc}`02-02-convex-functions` are used throughout.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- describe an iterative algorithm by its map, its descent function and its
  stopping rule, and classify a convergent sequence as sublinear, linear,
  superlinear or quadratic;
- reduce a line search along a direction to a one-dimensional minimization
  and state when a function of one variable is unimodal;
- run bisection and golden-section search by hand, derive the golden ratio
  $\tau$, and compute how many evaluations a given accuracy needs;
- derive Newton's method from a quadratic model, prove its local quadratic
  convergence and recognize its failure modes;
- fit a parabola through three points and use its vertex as the next trial
  point; and
- test a step length against the Armijo and Wolfe conditions and carry out a
  backtracking line search.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\vx^{(k)}$, $x_k$ | iterate $k$ of an algorithm in $\R^n$; the scalar iterate in one dimension |
| $\vd^{(k)}$, $\alpha_k$ | search direction and step length at iteration $k$ |
| $\varphi(\alpha)$ | function along a line, $\varphi(\alpha)=f(\vx+\alpha\vd)$ |
| $[a_k,b_k]$ | interval of uncertainty at iteration $k$ (one-dimensional search) |
| $p_k<q_k$ | the two trial points of golden-section search |
| $\tau$ | golden-section ratio $(\sqrt5-1)/2\approx0.6180$ |
| $e_k$ | error of iterate $k$, $e_k=\norm{\vx^{(k)}-\vx^*}$ (or $\abs{x_k-x^*}$) |
| $\varepsilon$ | tolerance of a stopping rule |
| $c_1$, $c_2$ | constants of the Armijo and curvature (Wolfe) conditions, $0<c_1<c_2<1$ |
| $\gamma$ | backtracking factor, $0<\gamma<1$ |
| $F_N$ | Fibonacci numbers, $F_0=F_1=1$, $F_N=F_{N-1}+F_{N-2}$ |

(sec-06-04-iterative-algorithms)=
## Iterative algorithms

Consider the nonlinear program

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & g_i(\vx)\le0, && i=1,\dots,m,\\
& h_j(\vx)=0, && j=1,\dots,\ell,
\end{aligned}
$$ (eq-06-04-nlp)

with feasible set $\Omega$ ({eq}`eq-06-01-feasible-set`). How can one design
a method that finds a local minimizer, or at least a KKT point, of
{eq}`eq-06-04-nlp`? Three facts shape every answer.

1. The feasible set $\Omega$ contains infinitely many points, so it cannot be
   searched exhaustively; a method must *move* from point to point.
2. At a given point the method can only compute local information: the
   values of $f$, $g_i$, $h_j$ and their first (sometimes second) derivatives.
   The direction of the next move must be built from this information.
3. The optimality conditions ({doc}`06-02-optimality-conditions`,
   {doc}`06-03-kkt-conditions`) cannot usually be solved directly, but they can be
   *checked* at a given point. They therefore serve as stopping rules.

A method built on these facts produces a sequence $\vx^{(0)},\vx^{(1)},\vx^{(2)},\dots$
in which each point is computed from the previous one by a fixed rule.

:::{prf:definition} Algorithmic map
:label: def-06-04-algorithm-map

An **iterative algorithm** is defined by a map $A$ that assigns to each point
$\vx$ of its domain a new point or a set of candidate points.

- A **point-to-point map** gives a single successor,
  $\vx^{(k+1)}=A(\vx^{(k)})$, $k=0,1,2,\dots$
- A **point-to-set map** gives a set $A(\vx^{(k)})$ of acceptable successors,
  and the algorithm chooses any $\vx^{(k+1)}\in A(\vx^{(k)})$.
:::

Point-to-set maps are the natural description of most methods, because one
iteration is usually a *composition* of simpler maps: choose a direction,
then choose a step along it, then (perhaps) round or project. When the line
search is inexact, or when the minimizer along the line is not unique, the
second map returns a set of acceptable steps rather than a single one.
Viewing the whole iteration as a point-to-set map is what allows a single
convergence theorem to cover many variants (Luenberger & Ye, 2016, Ch. 7;
Bazaraa et al., 2006, Ch. 7).

A method makes progress when some measure of quality improves from one
iterate to the next.

:::{prf:definition} Descent function
:label: def-06-04-descent-function

Let $S$ be the set of points regarded as solutions (for example, the
stationary points of $f$). A continuous function $M$ is a **descent function**
(merit function) for the algorithm $A$ and the solution set $S$ if
$M(\vy)<M(\vx)$ for every $\vx\notin S$ and every $\vy\in A(\vx)$, and
$M(\vy)\le M(\vx)$ for every $\vx\in S$ and every $\vy\in A(\vx)$. An
algorithm that has a descent function is called a **descent algorithm**.
:::

Common choices of $M$ for unconstrained minimization are

$$
\begin{gathered}
M(\vx)=f(\vx),\\
M(\vx)=\norm{\vx-\vx^*}^2,\\
M(\vx)=\norm{\nabla f(\vx)},\\
M(\vx)=-\min_{\norm{\vd}\le1}\nabla f(\vx)^\T\vd .
\end{gathered}
$$

The second is useful in proofs but not in computation, since $\vx^*$ is
unknown. The last one is the steepest rate of decrease available along a unit
direction; it is nonnegative and vanishes exactly at stationary points. The
direction must be normalized: without the bound $\norm{\vd}\le1$ the minimum
of the linear function $\nabla f(\vx)^\T\vd$ is $-\infty$ whenever
$\nabla f(\vx)\ne\vzero$. Without constraints the minimum is attained at
$\vd=-\nabla f(\vx)/\norm{\nabla f(\vx)}$ and the last measure equals
$\norm{\nabla f(\vx)}$; when $\vd$ is restricted to feasible directions of a
constrained problem it becomes a genuinely different stationarity measure.
For constrained problems the merit
function combines the objective with a measure of infeasibility, as the
penalty functions of {ref}`sec-06-06-penalty` do.

**Convergence.** An algorithm is **convergent** if, from the starting points
considered, the sequence it generates converges to a solution: a local
minimizer or a KKT point. Most theorems guarantee less, namely that every
*limit point* of a bounded sequence is a stationary point. Such a guarantee
is called **global convergence** when it holds from every starting point;
it says nothing about the speed. Theorems that need a starting point "close
enough" to $\vx^*$ are **local**.

(sec-06-04-descent-step)=
### Directions and steps

Almost every method of this chapter and the next two updates the iterate as

$$
\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)},\qquad \alpha_k>0,
$$ (eq-06-04-update)

where $\vd^{(k)}$ is a **search direction** and $\alpha_k$ a **step length**.
A direction is a **descent direction** at $\vx$
({prf:ref}`def-06-02-descent-direction`) if

$$
\nabla f(\vx)^\T\vd<0.
$$ (eq-06-04-descent-direction)

The condition is the right one because of the first-order Taylor expansion
$f(\vx+\alpha\vd)=f(\vx)+\alpha\nabla f(\vx)^\T\vd+o(\alpha)$: if
$\nabla f(\vx)^\T\vd<0$, the linear term dominates the remainder for small
$\alpha>0$, so $f(\vx+\alpha\vd)<f(\vx)$ for all sufficiently small
$\alpha>0$ (the full argument is {prf:ref}`lem-06-02-descent`). The
choice of $\vd^{(k)}$ distinguishes the methods of
{doc}`06-05-unconstrained-algorithms`; the choice of $\alpha_k$ is the **line
search**, the subject of the rest of this chapter.

**Stopping rules.** A practical algorithm stops when one of the following
holds, for user-chosen tolerances $\varepsilon$:

- the optimality condition holds approximately, $\norm{\nabla f(\vx^{(k)})}\le\varepsilon$
  (or the KKT conditions hold up to $\varepsilon$);
- the iterates or the objective values stagnate,
  $\norm{\vx^{(k+1)}-\vx^{(k)}}\le\varepsilon$ or
  $\abs{f(\vx^{(k+1)})-f(\vx^{(k)})}\le\varepsilon$, often measured relative to
  $\norm{\vx^{(k)}}$ or $\abs{f(\vx^{(k)})}$;
- a bracket that contains the solution is short enough, $b_k-a_k\le\varepsilon$
  (one-dimensional search); or
- a maximum number of iterations or function evaluations is reached (a safety
  net, not a sign of success).

Stagnation tests alone are risky: a slowly converging method can take tiny
steps far from a solution.

(sec-06-04-rates)=
### Rates of convergence

To compare methods we measure how fast the error $e_k=\norm{\vx^{(k)}-\vx^*}$
of a convergent sequence goes to zero.

:::{prf:definition} Rates of convergence
:label: def-06-04-convergence-rates

Let $\vx^{(k)}\to\vx^*$ with $e_k=\norm{\vx^{(k)}-\vx^*}>0$. The convergence is

- **linear** with rate $\beta\in(0,1)$ if $\limsup_{k\to\infty}e_{k+1}/e_k=\beta$;
- **superlinear** if $\lim_{k\to\infty}e_{k+1}/e_k=0$;
- **of order** $r>1$ if $e_{k+1}\le Ce_k^{\,r}$ for some constant $C$ and all
  large $k$; order $2$ is **quadratic** convergence;
- **sublinear** if $e_{k+1}/e_k\to1$.
:::

The sequences $e_k=2^{-k}$, $e_k=k^{-k}$, $e_k=2^{-2^k}$ and $e_k=1/k$ are
linear (rate $\tfrac12$), superlinear, quadratic and sublinear. A linearly
convergent method gains a fixed number of correct digits per iteration
(about $-\log_{10}\beta$); a quadratically convergent one roughly doubles the
number of correct digits per iteration once it is close; a sublinear one
needs about ten times more iterations for each extra digit when
$e_k\approx C/k$. Bisection (linear, $\beta=\tfrac12$), golden-section search
(linear, $\beta=\tau$) and Newton's method (quadratic) below are the first
examples; steepest descent (linear) and Newton's method in $\R^n$ (quadratic)
follow in {doc}`06-05-unconstrained-algorithms`, Frank–Wolfe (sublinear) in
{doc}`06-06-constrained-algorithms`.

(sec-06-04-one-dimensional)=
## The one-dimensional problem

Once $\vx=\vx^{(k)}$ and $\vd=\vd^{(k)}$ are fixed, the step length is found
by minimizing the function of one variable

$$
\begin{aligned}
\varphi(\alpha)&=f(\vx+\alpha\vd),\\
\varphi'(\alpha)&=\nabla f(\vx+\alpha\vd)^\T\vd,
\end{aligned}
$$ (eq-06-04-phi)

over an interval of steps, $0\le\alpha\le\bar\alpha$ (with $\bar\alpha=\infty$
in unconstrained problems). The derivative follows from the chain rule
({prf:ref}`thm-01-02-chain-rule`); in particular $\varphi'(0)=\nabla f(\vx)^\T\vd<0$
for a descent direction. Minimizing $\varphi$ exactly is an **exact line
search**; accepting any step that decreases $\varphi$ "enough" is an
**inexact line search** ({ref}`sec-06-04-inexact`).

When the move is towards a target point $\vy$, so that $\vd=\vy-\vx$, the
points of the line are convex combinations,

$$
\vx+\alpha(\vy-\vx)=(1-\alpha)\vx+\alpha\vy,
$$

and $\alpha\in[0,1]$ covers exactly the segment from $\vx$ ($\alpha=0$) to
$\vy$ ($\alpha=1$). The Frank–Wolfe method ({ref}`sec-06-06-frank-wolfe`)
works this way.

In the rest of the chapter the line search is written as the generic problem

$$
\min_{a\le x\le b} f(x)
$$ (eq-06-04-1d)

for a function $f$ of one real variable on an interval $[a,b]$, the
**interval of uncertainty**. All methods shrink this interval or move a
single point towards the minimizer $x^*$. They need some structure: an
arbitrary function can have many local minima in $[a,b]$, and no finite set of
evaluations can locate the global one. The structure used is unimodality.

:::{prf:definition} Unimodal (strictly quasiconvex) function
:label: def-06-04-unimodal

A function $f$ is **strictly quasiconvex** on $[a,b]$ if for any two points
$x_1\ne x_2$ in $[a,b]$ and any $\theta\in(0,1)$

$$
f(\theta x_1+(1-\theta)x_2)<\max\{f(x_1),f(x_2)\}.
$$

For a function that attains its minimum on $[a,b]$ (for example a
continuous one) this is the same as **unimodality**: $f$ has a unique
minimizer $x^*$ on $[a,b]$, is strictly decreasing on $[a,x^*]$ and strictly
increasing on $[x^*,b]$.
:::

Strict quasiconvexity strengthens the quasiconvexity of
{prf:ref}`def-02-02-quasiconvex`. A strictly convex function is unimodal, but a unimodal function need not be
convex or differentiable: $f(x)=\sqrt{\abs{x}}$ on $[-1,1]$ is unimodal
but neither, and $\max\{x^2,(x-1)^2\}$ is unimodal but not differentiable
at its minimizer $x=\tfrac12$. The next result is the
basis of all interval-elimination methods: two function values suffice to
discard part of the interval.

:::{prf:theorem} Interval elimination
:label: thm-06-04-interval-elimination

Let $f$ be strictly quasiconvex on $[a,b]$ and let $a\le p<q\le b$.

1. If $f(p)>f(q)$, then $f(x)>f(q)$ for every $x\in[a,p]$, so the minimizer
   lies in $(p,b]$.
2. If $f(p)<f(q)$, then $f(x)>f(p)$ for every $x\in[q,b]$, so the minimizer
   lies in $[a,q)$.
3. If $f(p)=f(q)$, the minimizer lies in $[p,q]$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

(1) Let $x\in[a,p)$. Then $p$ lies strictly between $x$ and $q$, say
$p=\theta x+(1-\theta)q$ with $\theta\in(0,1)$. Strict quasiconvexity gives
$f(p)<\max\{f(x),f(q)\}$. Since $f(p)>f(q)$, the maximum must be $f(x)$, so
$f(x)>f(p)>f(q)$; at $x=p$ the claim is the hypothesis. Hence no point of
$[a,p]$ is a minimizer. (2) is the mirror image: for $x\in(q,b]$ the point $q$
lies between $p$ and $x$, so $f(q)<\max\{f(p),f(x)\}$, and $f(q)>f(p)$ forces
$f(x)>f(q)>f(p)$. (3) If $f(p)=f(q)$, applying the argument of (1) to
$x\in[a,p)$ gives $f(p)<\max\{f(x),f(q)\}=\max\{f(x),f(p)\}$, so
$f(x)>f(p)$; similarly $f(x)>f(q)$ for $x\in(q,b]$.
:::

In words: discard the part of the interval beyond the trial point with the
**larger** value. In case 2 the proof gives the stronger bound $f(x)>f(q)$
for $x\in(q,b]$, the mirror image of case 1.

(sec-06-04-bisection)=
## Bisection search

Bisection uses the derivative. If $f$ is differentiable on $[a,b]$ and its
minimizer is interior, then $f'(x^*)=0$, and the search for $x^*$ becomes the
search for a root of $f'$. Assume

$$
f'(a)<0<f'(b),
$$ (eq-06-04-bisection-bracket)

and that $f'$ changes sign only once in $[a,b]$ (from negative to positive);
this holds, for example, when $f$ is convex and continuously differentiable.
Then the sign of $f'$ at the midpoint tells which half contains $x^*$: if
$f'(x_k)>0$ the function increases at $x_k$, so $x^*<x_k$; if $f'(x_k)<0$,
then $x^*>x_k$.

:::{prf:algorithm} Bisection search
:label: alg-06-04-bisection

**Input:** a differentiable $f$ on $[a,b]$ with $f'(a)<0<f'(b)$; a tolerance
$\varepsilon>0$.
**Output:** an estimate $\hat x$ with $\abs{\hat x-x^*}\le\varepsilon/2$.

1. **Initialize.** Set $a_0=a$, $b_0=b$, $k=0$.
2. **Stopping test.** If $b_k-a_k\le\varepsilon$, stop and return
   $\hat x=(a_k+b_k)/2$.
3. **Midpoint.** Compute $x_k=(a_k+b_k)/2$ and $f'(x_k)$.
4. **Halve.** If $f'(x_k)>0$, set $a_{k+1}=a_k$, $b_{k+1}=x_k$. If
   $f'(x_k)<0$, set $a_{k+1}=x_k$, $b_{k+1}=b_k$. If $f'(x_k)=0$, stop with
   $\hat x=x_k$.
5. **Next.** Set $k\leftarrow k+1$ and go to Step 2.
:::

Every iteration halves the interval, $b_{k+1}-a_{k+1}=\tfrac12(b_k-a_k)$, so
$b_k-a_k=(b-a)/2^k$. The stopping test holds as soon as $(b-a)/2^k\le\varepsilon$,
that is after

$$
k=\left\lceil\log_2\frac{b-a}{\varepsilon}\right\rceil
$$ (eq-06-04-bisection-count)

iterations, each with one evaluation of $f'$. Since $x^*\in[a_k,b_k]$ at all
times, the midpoint returned is within $\varepsilon/2$ of $x^*$. The error
bound halves at each step: bisection converges linearly with rate
$\tfrac12$, and it does so for every function satisfying the assumptions,
with no starting-point condition.

:::{prf:example} Bisection on $\sin x$
:label: ex-06-04-bisection-sine

Minimize $f(x)=\sin x$ on $[3,6]$ by bisection with $\varepsilon=0.05$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Here $f'(x)=\cos x$, with $f'(3)=-0.9900<0$ and $f'(6)=0.9602>0$. On
$[3,6]$ the cosine is negative to the left of $3\pi/2$ and positive to its
right, so $f'$ has the single root $x^*=3\pi/2=4.712389$ in the interval
and $f$ is unimodal there. Formula {eq}`eq-06-04-bisection-count`
gives $k=\lceil\log_2(3/0.05)\rceil=\lceil5.91\rceil=6$ halvings.
{numref}`tab-06-04-bisection` lists the iterations.

```{list-table} Bisection on $f'(x)=\cos x$, $[a_0,b_0]=[3,6]$, $\varepsilon=0.05$
:name: tab-06-04-bisection
:header-rows: 1

* - $k$
  - $a_k$
  - $b_k$
  - $b_k-a_k$
  - $x_k$
  - $f'(x_k)$
  - kept half
* - 0
  - 3
  - 6
  - 3
  - 4.5
  - $-0.2108$
  - right
* - 1
  - 4.5
  - 6
  - 1.5
  - 5.25
  - $0.5121$
  - left
* - 2
  - 4.5
  - 5.25
  - 0.75
  - 4.875
  - $0.1619$
  - left
* - 3
  - 4.5
  - 4.875
  - 0.375
  - 4.6875
  - $-0.02489$
  - right
* - 4
  - 4.6875
  - 4.875
  - 0.1875
  - 4.78125
  - $0.06881$
  - left
* - 5
  - 4.6875
  - 4.78125
  - 0.09375
  - 4.734375
  - $0.02198$
  - left
* - 6
  - 4.6875
  - 4.734375
  - 0.046875
  - stop
  -
  -
```

At $k=6$ the width $0.046875\le0.05$ triggers the stopping test, and the
estimate is the midpoint $\hat x=(4.6875+4.734375)/2=4.7109375$, with
$f'(\hat x)=-0.001451$. Its error is $\abs{\hat x-x^*}=0.001451$, well inside
the guaranteed $\varepsilon/2=0.025$; six derivative evaluations were used.
*(In some printed versions of this table the right end of the last interval
is given as $4.73475$, a typo for $4.734375$, and the last midpoint and
derivative, $4.71113$ and $-0.0013$, are computed from it; the correct values
are $4.710938$ and $-0.001451$.)*
:::

```{figure} figures/06-04-interval-reduction.svg
:name: fig-06-04-interval-reduction
:width: 640px
:class: dark-light
:alt: Two panels of nested horizontal intervals shrinking towards the minimizer 3 pi over 2 of sin x; left bisection with midpoints, right golden-section search with two trial points per row

Intervals of uncertainty for $\min\sin x$ on $[3,6]$. (a) Bisection
({prf:ref}`ex-06-04-bisection-sine`): each row halves the interval at the
midpoint (dot). (b) Golden-section search ({prf:ref}`ex-06-04-golden-sine`): each
row keeps a fraction $\tau\approx0.618$ of the interval; the filled dot is the
trial point reused from the row above, the open dot the one new evaluation.
The dashed line is $x^*=3\pi/2$.
```

(sec-06-04-golden-section)=
## Golden-section search

Bisection needs derivatives. **Golden-section search** needs only function
values and assumes only unimodality ({prf:ref}`def-06-04-unimodal`), so it also
works for non-differentiable functions. It places two trial points
$p_k<q_k$ in $[a_k,b_k]$ and discards one end by
{prf:ref}`thm-06-04-interval-elimination`. The positions are chosen with two
goals.

- **Symmetry.** The two possible new intervals, $[a_k,q_k]$ and $[p_k,b_k]$,
  should have the same length, so that the worst case is as good as the best.
  This requires $p_k-a_k=b_k-q_k$. Writing the common fraction of the interval
  kept as $\tau$,

  $$
  \begin{aligned}
  p_k&=a_k+(1-\tau)(b_k-a_k),\\
  q_k&=a_k+\tau(b_k-a_k),
  \end{aligned}
  $$ (eq-06-04-golden-points)

  and the new interval has length $b_{k+1}-a_{k+1}=\tau(b_k-a_k)$.
- **Reuse.** The surviving trial point should be one of the two trial points
  of the next interval, so that each iteration needs only **one** new
  function value.

Suppose $[a_k,q_k]$ is kept (case 2 of the theorem). Its trial points are at
fractions $1-\tau$ and $\tau$ of its length $\tau(b_k-a_k)$, measured from
$a_k$. The old point $p_k$, at distance $(1-\tau)(b_k-a_k)$ from $a_k$,
coincides with the new *right* trial point $q_{k+1}$ if

$$
\begin{aligned}
(1-\tau)(b_k-a_k)&=\tau\cdot\tau(b_k-a_k),\\
\text{that is,}\quad\tau^2&=1-\tau .
\end{aligned}
$$ (eq-06-04-golden-equation)

By symmetry the same condition makes $q_k$ the new left point $p_{k+1}$ when
$[p_k,b_k]$ is kept. The positive root of $\tau^2+\tau-1=0$ is

$$
\begin{aligned}
\tau&=\frac{\sqrt5-1}{2}\approx0.618034,\\
1-\tau&=\tau^2\approx0.381966,
\end{aligned}
$$ (eq-06-04-tau)

the reciprocal of the golden ratio $(1+\sqrt5)/2$; hence the name.
{numref}`fig-06-04-golden-section-geometry` shows the two levels.

```{figure} figures/06-04-golden-section-geometry.svg
:name: fig-06-04-golden-section-geometry
:width: 560px
:class: dark-light
:alt: Two stacked number lines; the upper interval a to b with trial points p and q, the lower interval a to q in which the old point p has become the new right trial point

Placement of the golden-section trial points. Upper line: $[a_k,b_k]$ with
$p_k$ and $q_k$ at fractions $1-\tau=\tau^2$ and $\tau$ of its length. If
$f(p_k)<f(q_k)$, the part $(q_k,b_k]$ is discarded. Lower line: in the new
interval $[a_k,q_k]$ of length $\tau(b_k-a_k)$, the old point $p_k$ sits
exactly at fraction $\tau$, so it becomes $q_{k+1}$ and only $p_{k+1}$ must
be evaluated.
```

:::{prf:algorithm} Golden-section search
:label: alg-06-04-golden-section

**Input:** a strictly quasiconvex $f$ on $[a,b]$; a tolerance
$\varepsilon>0$.
**Output:** an estimate $\hat x$ with $\abs{\hat x-x^*}\le\varepsilon/2$.

1. **Initialize.** Set $a_0=a$, $b_0=b$, $p_0=a+(1-\tau)(b-a)$,
   $q_0=a+\tau(b-a)$; evaluate $f(p_0)$ and $f(q_0)$; set $k=0$.
2. **Stopping test.** If $b_k-a_k\le\varepsilon$, stop and return
   $\hat x=(a_k+b_k)/2$.
3. **Eliminate and reuse.**
   - If $f(p_k)>f(q_k)$: set $a_{k+1}=p_k$, $b_{k+1}=b_k$,
     $p_{k+1}=q_k$ (value known), and evaluate $f$ at the new point
     $q_{k+1}=a_{k+1}+\tau(b_{k+1}-a_{k+1})$.
   - If $f(p_k)\le f(q_k)$: set $a_{k+1}=a_k$, $b_{k+1}=q_k$,
     $q_{k+1}=p_k$ (value known), and evaluate $f$ at the new point
     $p_{k+1}=a_{k+1}+(1-\tau)(b_{k+1}-a_{k+1})$.
4. **Next.** Set $k\leftarrow k+1$ and go to Step 2.
:::

A tie $f(p_k)=f(q_k)$ is placed in the second case; by case 3 of
{prf:ref}`thm-06-04-interval-elimination` the minimizer lies in
$[p_k,q_k]\subset[a_k,q_k]$, so nothing is lost.

**How many evaluations?** After $k$ reductions the interval has length
$\tau^k(b-a)$, and the first $k$ reductions cost $2+(k-1)=k+1$ function
values. With $N$ evaluations, therefore, the final interval has length
$\tau^{N-1}(b-a)$, and the accuracy $b_k-a_k\le\varepsilon$ needs

$$
\begin{aligned}
k&=\left\lceil\frac{\ln(\varepsilon/(b-a))}{\ln\tau}\right\rceil\ \text{reductions},\\
N&=k+1\ \text{evaluations.}
\end{aligned}
$$ (eq-06-04-golden-count)

The method converges linearly with rate $\tau\approx0.618$, slower per step
than bisection ($\tfrac12$): reducing an interval by a factor $10^{6}$ takes
$30$ function values against $20$ derivative values for bisection. In exchange
golden-section search uses no derivatives and tolerates kinks.

:::{prf:example} Golden-section search on $\sin x$
:label: ex-06-04-golden-sine

Minimize $f(x)=\sin x$ on $[3,6]$ by golden-section search with
$\varepsilon=0.05$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Formula {eq}`eq-06-04-golden-count` gives
$k=\lceil\ln(0.05/3)/\ln0.618\rceil=\lceil8.51\rceil=9$ reductions and
$N=10$ function values. The initial trial points are
$p_0=3+0.381966\cdot3=4.1459$ and $q_0=3+0.618034\cdot3=4.8541$. Function
values are printed with six decimals in {numref}`tab-06-04-golden` because near
the minimum they differ only in the fourth decimal.

```{list-table} Golden-section search on $\sin x$, $[a_0,b_0]=[3,6]$, $\varepsilon=0.05$ (the reused point of each row is the one shown in bold)
:name: tab-06-04-golden
:header-rows: 1

* - $k$
  - $a_k$
  - $b_k$
  - $b_k-a_k$
  - $p_k$
  - $q_k$
  - $f(p_k)$
  - $f(q_k)$
* - 0
  - 3
  - 6
  - 3
  - 4.1459
  - 4.8541
  - $-0.843789$
  - $-0.989976$
* - 1
  - 4.1459
  - 6
  - 1.8541
  - **4.8541**
  - 5.2918
  - $-0.989976$
  - $-0.836787$
* - 2
  - 4.1459
  - 5.2918
  - 1.1459
  - 4.5836
  - **4.8541**
  - $-0.991717$
  - $-0.989976$
* - 3
  - 4.1459
  - 4.8541
  - 0.7082
  - 4.4164
  - **4.5836**
  - $-0.956516$
  - $-0.991717$
* - 4
  - 4.4164
  - 4.8541
  - 0.4377
  - **4.5836**
  - 4.6869
  - $-0.991717$
  - $-0.999676$
* - 5
  - 4.5836
  - 4.8541
  - 0.2705
  - **4.6869**
  - 4.7508
  - $-0.999676$
  - $-0.999263$
* - 6
  - 4.5836
  - 4.7508
  - 0.1672
  - 4.6475
  - **4.6869**
  - $-0.997892$
  - $-0.999676$
* - 7
  - 4.6475
  - 4.7508
  - 0.1033
  - **4.6869**
  - 4.7113
  - $-0.999676$
  - $-0.999999$
* - 8
  - 4.6869
  - 4.7508
  - 0.06386
  - **4.7113**
  - 4.7264
  - $-0.999999$
  - $-0.999902$
* - 9
  - 4.6869
  - 4.7264
  - 0.03947
  - stop
  -
  -
  -
```

Row 0: $f(p_0)>f(q_0)$, so $[3,p_0]$ is discarded, the new interval is
$[4.1459,6]$, the old $q_0$ becomes $p_1$, and only $q_1=4.1459+0.618034\cdot1.8541=5.2918$
is new. Row 1: $f(p_1)<f(q_1)$, so $(q_1,6]$ goes, and $p_1$ becomes $q_2$.
Each later row follows the same pattern, and every width is $\tau$ times the
previous one ($3\tau^8=0.06386$). The search stops at $k=9$ with
$[a_9,b_9]=[4.6869,4.7264]$, which contains $x^*=4.7124$; the estimate is
$\hat x=4.7067$, error $0.0057$, after $10$ function values. Bisection reached
a smaller final interval with six derivative values; golden-section search
used no derivative at all.
:::

**Fibonacci search.** If the number $N$ of function values is fixed in
advance, the symmetric placement with reuse can be optimized over all $N$
steps instead of being the same at every step. The optimal strategy, found by
Kiefer (1953), places the trial points at ratios of consecutive Fibonacci
numbers $F_{N-1}/F_N$, $F_{N-2}/F_{N-1}$, … ($F_0=F_1=1$,
$F_N=F_{N-1}+F_{N-2}$) and reduces the interval by the factor $1/F_N$, the
best possible with $N$ evaluations (ignoring a small distinguishability
constant in the last step). Since $F_{N-1}/F_N\to\tau$, golden-section search
is the limit of Fibonacci search, and it loses little: after $N=10$
evaluations the factors are $1/F_{10}=1/89=0.01124$ and $\tau^9=0.01316$, and
for large $N$ the golden-section interval is only about $17\,\%$ longer.
Golden-section search is preferred in practice because it does not need $N$
in advance.

(sec-06-04-newton)=
## Newton's method

Bisection and golden-section search shrink an interval by a constant factor
per evaluation. If $f$ is twice differentiable, a much faster method is
obtained by **curve fitting**: replace $f$ near the current point $x_k$ by the
quadratic that agrees with it in value, slope and curvature,

$$
q(x)=f(x_k)+f'(x_k)(x-x_k)+\tfrac12f''(x_k)(x-x_k)^2,
$$ (eq-06-04-newton-model)

and take its stationary point as the next iterate. Setting
$q'(x)=f'(x_k)+f''(x_k)(x-x_k)=0$ gives

$$
x_{k+1}=x_k-\frac{f'(x_k)}{f''(x_k)}.
$$ (eq-06-04-newton-1d)

If $f''(x_k)>0$ the model is a convex parabola and $x_{k+1}$ is its minimizer.
Equation {eq}`eq-06-04-newton-1d` is also Newton's method for solving
$f'(x)=0$: $x_{k+1}$ is where the tangent to the graph of $f'$ at $x_k$ meets
the axis ({numref}`fig-06-04-newton-tangents`).

:::{prf:algorithm} Newton's method in one dimension
:label: alg-06-04-newton-1d

**Input:** a twice differentiable $f$; a starting point $x_0$; a tolerance
$\varepsilon>0$.
**Output:** an approximate stationary point of $f$.

1. **Initialize.** Set $k=0$.
2. **Stopping test.** If $\abs{f'(x_k)}\le\varepsilon$, stop and return $x_k$.
3. **Newton step.** Set $x_{k+1}=x_k-f'(x_k)/f''(x_k)$.
4. **Next.** Set $k\leftarrow k+1$ and go to Step 2.
:::

:::{prf:theorem} Local quadratic convergence of Newton's method
:label: thm-06-04-newton-quadratic

Let $f$ be three times continuously differentiable near a point $x^*$ with
$f'(x^*)=0$ and $f''(x^*)\ne0$. Then there is a $\delta>0$ such that for every
$x_0$ with $\abs{x_0-x^*}<\delta$ the iterates {eq}`eq-06-04-newton-1d` are
defined, converge to $x^*$, and satisfy

$$
\begin{aligned}
&\abs{x_{k+1}-x^*}\le C\abs{x_k-x^*}^2,\\
&\lim_{k\to\infty}\frac{x_{k+1}-x^*}{(x_k-x^*)^2}=\frac{f'''(x^*)}{2f''(x^*)}.
\end{aligned}
$$ (eq-06-04-newton-error)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Choose $\delta_0>0$ such that $\abs{f''(x)}\ge m>0$ and $\abs{f'''(x)}\le L$ on
$I=[x^*-\delta_0,x^*+\delta_0]$ (possible by continuity and $f''(x^*)\ne0$).
For $x_k\in I$, Taylor's theorem with remainder around $x_k$ gives, for some
$\xi_k$ between $x_k$ and $x^*$,

$$
\begin{aligned}
0=f'(x^*)&=f'(x_k)+f''(x_k)(x^*-x_k)\\
&\quad+\tfrac12f'''(\xi_k)(x^*-x_k)^2.
\end{aligned}
$$

Divide by $f''(x_k)$ and use {eq}`eq-06-04-newton-1d`,
$x_k-f'(x_k)/f''(x_k)=x_{k+1}$:

$$
x_{k+1}-x^*=\frac{f'''(\xi_k)}{2f''(x_k)}(x_k-x^*)^2 .
$$

Hence $\abs{x_{k+1}-x^*}\le C\abs{x_k-x^*}^2$ with $C=L/(2m)$. If in addition
$\delta=\min\{\delta_0,1/(2C)\}$ and $\abs{x_k-x^*}<\delta$, then
$\abs{x_{k+1}-x^*}\le C\delta\abs{x_k-x^*}\le\tfrac12\abs{x_k-x^*}$, so $x_{k+1}$
stays in $I$ and, by induction, $x_k\to x^*$. Since $\xi_k\to x^*$ and
$x_k\to x^*$, the ratio converges to $f'''(x^*)/(2f''(x^*))$.
:::

The theorem is local and says nothing about *which* stationary point is
found: Newton's method looks for zeros of $f'$, and a maximizer
($f''(x^*)<0$) attracts the iterates just as well as a minimizer. The
following failures are typical.

- **Negative curvature.** If $f''(x_k)<0$ the model {eq}`eq-06-04-newton-model`
  is a downward parabola, and the "Newton step" goes to its maximizer, often
  uphill.
- **Small curvature.** If $f''(x_k)\approx0$ the step $-f'(x_k)/f''(x_k)$ is
  huge and the iterate is thrown far away; if $f''(x_k)=0$ it is undefined.
- **Leaving the interval.** Equation {eq}`eq-06-04-newton-1d` ignores the bounds
  $a\le x\le b$.
- **Cycling or divergence** from poor starting points, even for convex
  functions; for example, Newton's method for the root of $f'(x)=\arctan x$
  (the derivative of a strictly convex $f$) diverges from $\abs{x_0}>1.39$.

Practical one-dimensional codes therefore **safeguard** the Newton step: they
keep a bracket $[a_k,b_k]$ with $f'(a_k)<0<f'(b_k)$ and fall back on a
bisection step whenever the Newton point is outside the bracket or the
curvature is not positive. Near the minimizer the Newton steps are then
accepted and the fast local convergence is retained.

:::{prf:example} Newton's method on $\sin x$ from two starting points
:label: ex-06-04-newton-sine

Minimize $f(x)=\sin x$ on $[3,6]$ by Newton's method with $\varepsilon=10^{-6}$,
starting (a) from $x_0=4$ and (b) from $x_0=3$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The derivatives are $f'(x)=\cos x$ and $f''(x)=-\sin x$, so
{eq}`eq-06-04-newton-1d` reads $x_{k+1}=x_k+\cos x_k/\sin x_k$.

```{list-table} Newton's method on $\sin x$ (values to 4 significant digits)
:name: tab-06-04-newton
:header-rows: 1

* - start
  - $k$
  - $x_k$
  - $f'(x_k)$
  - $f''(x_k)$
  - $x_{k+1}$
  - $\abs{x_k-3\pi/2}$
* - (a) $x_0=4$
  - 0
  - 4.000
  - $-0.6536$
  - $0.7568$
  - 4.864
  - $0.7124$
* -
  - 1
  - 4.864
  - $0.1507$
  - $0.9886$
  - 4.711
  - $0.1513$
* -
  - 2
  - 4.711
  - $-1.165\times10^{-3}$
  - $1.000$
  - 4.712
  - $1.165\times10^{-3}$
* -
  - 3
  - 4.712
  - $5.274\times10^{-10}$
  - $1.000$
  - stop
  - $5.274\times10^{-10}$
* - (b) $x_0=3$
  - 0
  - 3.000
  - $-0.9900$
  - $-0.1411$
  - $-4.015$
  - $1.712$
* -
  - 1
  - $-4.015$
  - $-0.6420$
  - $-0.7667$
  - $-4.853$
  -
* -
  - 2
  - $-4.853$
  - $0.1398$
  - $-0.9902$
  - $-4.711$
  -
* -
  - 3
  - $-4.711$
  - $-9.272\times10^{-4}$
  - $-1.000$
  - $-4.712$
  -
* -
  - 4
  - $-4.712$
  - $2.657\times10^{-10}$
  - $-1.000$
  - stop
  -
```

**(a)** From $x_0=4$ the curvature $f''(4)=0.7568$ is positive, the model is
convex, and the iterates $4$ → $4.8637$ → $4.7112$ → $4.712389$ reach
$x^*=3\pi/2$ to ten digits in three steps. The errors
$0.7124$, $0.1513$, $1.165\times10^{-3}$, $5.3\times10^{-10}$ shrink even faster
than quadratically: here $f'''(x^*)=-\cos(3\pi/2)=0$, so the limit in
{eq}`eq-06-04-newton-error` is zero, and indeed
$x_{k+1}-x^*\approx-(x_k-x^*)^3/3$ (the ratios $e_{k+1}/e_k^3$ are
$0.42,\ 0.34,\ 0.333$). Four derivative pairs gave far more accuracy than the
six or ten evaluations of the interval methods.

**(b)** From $x_0=3$ the curvature $f''(3)=-0.1411$ is negative and small. The
model is a flat downward parabola, and its stationary point
$x_1=3-(-0.9900)/(-0.1411)=-4.0153$ lies far outside $[3,6]$. From there the
iteration converges, again very fast, to $-4.712389=-3\pi/2$, where
$f''=-1<0$ and $f=\sin(-3\pi/2)=1$: a **maximizer** of $\sin x$, not a
minimizer, and not even in the interval. The stopping test
$\abs{f'(x_k)}\le\varepsilon$ is satisfied, so the algorithm reports success.
:::

Part (b) of {prf:ref}`ex-06-04-newton-sine` shows the three failure modes
at once and why the sign of $f''(x_k)$ and the bounds must be checked; a
safeguarded method started at $x_0=3$ would reject the step and bisect the
bracket $[3,6]$ instead.

```{figure} figures/06-04-newton-tangents.svg
:name: fig-06-04-newton-tangents
:width: 640px
:class: dark-light
:alt: Two panels showing the graph of cos x with tangent lines; on the left tangents from x0 equals 4 lead quickly to the root 3 pi over 2; on the right a nearly flat tangent at x0 equals 3 meets the axis near minus 4

Newton's method as tangent steps on $f'(x)=\cos x$
({prf:ref}`ex-06-04-newton-sine`). (a) From $x_0=4$ the tangents (dashed) meet
the axis at $x_1=4.864$ and $x_2=4.711$, next to the root $3\pi/2$ (a
minimizer of $\sin x$). (b) From $x_0=3$ the tangent is almost flat
($f''(3)=-0.141$) and meets the axis at $x_1=-4.015$; the iterates then
converge to $-3\pi/2$, a maximizer of $\sin x$ outside $[3,6]$.
```

(sec-06-04-curve-fitting)=
## Quadratic interpolation

Newton's method fits a quadratic using derivatives at one point. When
derivatives are unavailable or expensive, a quadratic can be fitted through
the *values* at three points $x_1<x_2<x_3$ with $f(x_2)<\min\{f(x_1),f(x_3)\}$
(a bracket of the minimum). The parabola through the three points has its
vertex at

$$
x_{\text{new}}=x_2-\frac12\,
\frac{(x_2-x_1)^2\bigl(f(x_2)-f(x_3)\bigr)-(x_2-x_3)^2\bigl(f(x_2)-f(x_1)\bigr)}
{(x_2-x_1)\bigl(f(x_2)-f(x_3)\bigr)-(x_2-x_3)\bigl(f(x_2)-f(x_1)\bigr)} ,
$$ (eq-06-04-quadratic-vertex)

which lies in $(x_1,x_3)$ because of the bracket condition. The new point
replaces one of the old ones so that the best point and its two neighbors
again form a bracket, and the process repeats. For $\sin x$ with the starting
points $4$, $4.5$, $5$ the vertices are $4.711130$, $4.712581$, $4.712393$ and
$4.712389$, with errors $1.3\times10^{-3}$, $1.9\times10^{-4}$,
$3.7\times10^{-6}$ and $1.1\times10^{-11}$: one function value per iteration
and superlinear convergence (of order about $1.3$ for smooth functions;
Luenberger & Ye, 2016, Ch. 8). Like Newton's method, quadratic
interpolation can misbehave far from the minimum, so robust codes combine it
with golden-section steps (Brent's method; Brent, 1973).

(sec-06-04-inexact)=
## Inexact line search: the Armijo and Wolfe conditions

Inside a method for $\R^n$, an exact line search is rarely worth its cost:
the direction $\vd^{(k)}$ will change at the next iteration anyway, and what
matters is a step that makes reasonable progress. An inexact line search
accepts any $\alpha$ that passes two tests. With $\varphi$ from
{eq}`eq-06-04-phi` and constants $0<c_1<c_2<1$:

$$
\varphi(\alpha)\le\varphi(0)+c_1\alpha\varphi'(0)
\qquad\text{(sufficient decrease, Armijo condition)},
$$ (eq-06-04-armijo)

$$
\varphi'(\alpha)\ge c_2\varphi'(0)
\qquad\text{(curvature condition)}.
$$ (eq-06-04-curvature)

Together they are the **Wolfe conditions** (Armijo, 1966; Wolfe, 1969). Condition {eq}`eq-06-04-armijo` says
that the actual decrease is at least a fraction $c_1$ of the decrease
predicted by the tangent line; it rules out steps that are too long. Since
$\varphi'(0)<0$, the right-hand side is a line of slope $c_1\varphi'(0)$, less
steep than the tangent, so all sufficiently small $\alpha>0$ satisfy it.
Condition {eq}`eq-06-04-curvature` rules out steps that are too short: it asks
the slope to have increased to at least $c_2\varphi'(0)$, so the step has
moved past the steep initial part. The **strong Wolfe conditions** replace
{eq}`eq-06-04-curvature` by $\abs{\varphi'(\alpha)}\le c_2\abs{\varphi'(0)}$, which
also excludes steps where $\varphi$ is already rising steeply. Typical values
are $c_1=10^{-4}$ and $c_2=0.9$ (Newton and quasi-Newton methods) or $c_2=0.1$
(conjugate gradients).

If $\varphi$ is continuously differentiable and bounded below on
$\alpha>0$, there are intervals of steps that satisfy the Wolfe and the strong
Wolfe conditions (Nocedal & Wright, 2006, Lemma 3.1). A simple way to satisfy
the Armijo condition alone is backtracking.

:::{prf:algorithm} Backtracking (Armijo) line search
:label: alg-06-04-backtracking

**Input:** $\vx$, a descent direction $\vd$ ($\nabla f(\vx)^\T\vd<0$), an initial
step $\bar\alpha>0$, constants $c_1\in(0,1)$ and $\gamma\in(0,1)$ (for
example $\gamma=\tfrac12$).
**Output:** a step $\alpha$ satisfying {eq}`eq-06-04-armijo`.

1. **Initialize.** Set $\alpha=\bar\alpha$.
2. **Test.** If $f(\vx+\alpha\vd)\le f(\vx)+c_1\alpha\nabla f(\vx)^\T\vd$, stop
   and return $\alpha$.
3. **Shrink.** Set $\alpha\leftarrow\gamma\alpha$ and go to Step 2.
:::

The loop terminates because the Armijo condition holds for all small
$\alpha$. Starting from $\bar\alpha=1$ is natural for Newton-type directions,
whose unit step is correct near the solution
({doc}`06-05-unconstrained-algorithms`). Unless it is $\bar\alpha$ itself, the
accepted step is $\gamma$ times a step that was rejected as too long, so it
cannot be needlessly short; this is why backtracking can do without the
curvature condition.

As an illustration take again $f(x)=\sin x$, now as a function on $\R$, at
$x=3$ with the descent direction $d=-f'(3)=0.9900$, so
$\varphi(\alpha)=\sin(3+0.99\alpha)$, $\varphi(0)=0.1411$ and
$\varphi'(0)=-0.9801$. With $c_1=0.25$, $\bar\alpha=4$ and $\gamma=\tfrac12$:

- $\alpha=4$: $\varphi(4)=0.6263>\varphi(0)+c_1\cdot4\cdot\varphi'(0)=-0.8390$;
  rejected (the step overshoots into the next rise of $\sin x$);
- $\alpha=2$: $\varphi(2)=-0.9644\le0.1411-0.25\cdot2\cdot0.9801=-0.3489$;
  accepted, giving $x=4.9800$.

The accepted step also satisfies the strong Wolfe conditions with $c_2=0.5$:
$\abs{\varphi'(2)}=0.2618\le0.5\cdot0.9801=0.4900$. {numref}`fig-06-04-armijo-wolfe`
shows the acceptable sets: the Armijo condition holds for
$0<\alpha\le2.749$, the curvature condition for $1.207\le\alpha\le5.426$
(the slope falls below $c_2\varphi'(0)$ again only after the next local
maximum), so the Wolfe steps form $[1.207,2.749]$ and the strong
Wolfe steps $[1.207,2.253]$; the exact minimizer $\alpha^*=1.730$ (which
reaches $x^*=3\pi/2$) lies inside both.

```{figure} figures/06-04-armijo-wolfe.svg
:name: fig-06-04-armijo-wolfe
:width: 560px
:class: dark-light
:alt: Graph of phi of alpha equals sin of 3 plus 0.99 alpha with the tangent at zero, the flatter Armijo line, shaded intervals of acceptable steps and two backtracking trials marked

The Armijo and Wolfe conditions for $\varphi(\alpha)=\sin(3+0.99\alpha)$ with
$c_1=0.25$, $c_2=0.5$. Below the Armijo line $\varphi(0)+c_1\alpha\varphi'(0)$
(dashed) the decrease is sufficient. Bars: steps that satisfy the Armijo
condition, the Wolfe conditions and the strong Wolfe conditions. Backtracking
from $\bar\alpha=4$ rejects $\alpha=4$ (open marker) and accepts $\alpha=2$
(filled marker).
```

## Summary and Common Mistakes

### Summary

- An iterative algorithm is a (point-to-point or point-to-set) map
  ({prf:ref}`def-06-04-algorithm-map`) that decreases a descent function
  ({prf:ref}`def-06-04-descent-function`); most methods take steps
  $\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$ along descent directions
  {eq}`eq-06-04-descent-direction` and stop on approximate optimality.
- Convergence speed is linear, superlinear, quadratic or sublinear
  ({prf:ref}`def-06-04-convergence-rates`).
- A line search minimizes $\varphi(\alpha)=f(\vx+\alpha\vd)$
  {eq}`eq-06-04-phi`; for unimodal functions two values suffice to discard part
  of the interval ({prf:ref}`thm-06-04-interval-elimination`).
- Bisection ({prf:ref}`alg-06-04-bisection`) halves the interval with one
  derivative per step and needs $\lceil\log_2((b-a)/\varepsilon)\rceil$ steps
  {eq}`eq-06-04-bisection-count`.
- Golden-section search ({prf:ref}`alg-06-04-golden-section`) keeps a fraction
  $\tau=(\sqrt5-1)/2$ per step, from $\tau^2=1-\tau$ {eq}`eq-06-04-golden-equation`,
  and needs one function value per step; Fibonacci search is its optimal
  finite-step version.
- Newton's method $x_{k+1}=x_k-f'(x_k)/f''(x_k)$ {eq}`eq-06-04-newton-1d`
  converges quadratically near a point with $f''(x^*)\ne0$
  ({prf:ref}`thm-06-04-newton-quadratic`) but can go to a maximizer or leave the
  interval ({numref}`fig-06-04-newton-tangents`).
- Quadratic interpolation {eq}`eq-06-04-quadratic-vertex` uses values only and
  converges superlinearly.
- Inexact line searches accept steps satisfying the Armijo
  {eq}`eq-06-04-armijo` and curvature {eq}`eq-06-04-curvature` conditions;
  backtracking ({prf:ref}`alg-06-04-backtracking`) enforces the first.

### Common mistakes

- **Writing the line as $\alpha\vx+(1-\alpha)\vy$.** With $\vd=\vy-\vx$ the
  point $\vx+\alpha\vd$ equals $(1-\alpha)\vx+\alpha\vy$; the weights are the
  other way round, and $\alpha=0$ must give the current point.
- **Discarding the wrong end.** In interval elimination the part beyond the
  trial point with the *larger* value is discarded; with equal values the
  minimizer is between the two points.
- **Recomputing both golden-section points.** With $\tau=(\sqrt5-1)/2$ one
  old point is always reused; rounding $\tau$ to $0.618$ and recomputing both
  points doubles the work and makes the trial points drift.
- **Trusting Newton's stopping test.** $\abs{f'(x_k)}\le\varepsilon$ holds at
  maximizers too; check $f''(x_k)>0$ and the bounds
  ({prf:ref}`ex-06-04-newton-sine`, start $x_0=3$).
- **Stopping on small steps only.** A slow method takes small steps far from
  the solution; combine stagnation tests with an optimality test.
- **Applying bisection without a sign change.** Bisection needs
  $f'(a)<0<f'(b)$; if both signs are equal, the minimizer is at an end point
  or the bracket is wrong.


## References

- Armijo, L. (1966). Minimization of functions having Lipschitz continuous
  first partial derivatives. *Pacific Journal of Mathematics, 16*(1), 1–3.
  <https://doi.org/10.2140/pjm.1966.16.1>
- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 7: the
  concept of an algorithm, algorithmic maps, descent functions and
  convergence; Sections 8.1–8.3: line search without and with derivatives,
  golden-section and Fibonacci search, bisection, Newton's method, quadratic
  fit and the Armijo rule.)
- Brent, R. P. (1973). *Algorithms for minimization without derivatives*.
  Prentice-Hall. (Chapter 5: one-dimensional minimization combining
  golden-section search and parabolic interpolation.)
- Kiefer, J. (1953). Sequential minimax search for a maximum. *Proceedings of
  the American Mathematical Society, 4*(3), 502–506.
  <https://doi.org/10.1090/S0002-9939-1953-0055639-3>
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming* (4th
  ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3> (Chapter 7:
  algorithms, global convergence and speed of convergence; Chapter 8: line
  search by Fibonacci and golden-section search, Newton's method, curve
  fitting and inexact line search.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Chapter 2:
  overview of line-search methods; Appendix A.2: rates of convergence;
  Sections 3.1 and 3.5: the Wolfe conditions, backtracking and step-length selection by
  interpolation.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 21: nonlinear programming algorithms, unconstrained
  one-variable search.)
- Wolfe, P. (1969). Convergence conditions for ascent methods. *SIAM Review,
  11*(2), 226–235. <https://doi.org/10.1137/1011036>
