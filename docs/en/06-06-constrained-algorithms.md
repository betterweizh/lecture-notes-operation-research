# Algorithms for Constrained Optimization

## Overview

The methods of {doc}`06-05-unconstrained-algorithms` move freely in $\R^n$. With
constraints, two things go wrong: the steepest-descent direction may leave
the feasible set, and a stationary point of $f$ need not be feasible at all.
This chapter presents the two classical ways around these difficulties. The
first keeps every iterate feasible and moves along **feasible descent
directions**; its simplest representative, the **Frank–Wolfe method**, finds
the direction by solving a *linear program* (LP) over the feasible polyhedron, so
that each iteration combines the simplex method of {doc}`03-03-simplex-method`
with a line search of {doc}`06-04-line-search`. The method also delivers a lower
bound on the optimal value at every step, but converges only sublinearly. The
second approach removes the constraints by adding a **penalty** for their
violation to the objective and solves a sequence of unconstrained problems
with a growing penalty parameter $\rho$. The minimizers $\vx(\rho)$ approach a
solution from outside the feasible set, the penalty terms estimate the Karush–Kuhn–Tucker
(KKT) multipliers of {doc}`06-03-kkt-conditions`, and the price is increasing
ill-conditioning. Barrier and augmented-Lagrangian methods, which repair
parts of this weakness, are outlined at the end. Convex sets and functions
({doc}`02-01-convex-sets`, {doc}`02-02-convex-functions`) and the feasible and
descent directions of {doc}`06-02-optimality-conditions` are used
throughout.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- explain why the negative gradient may be infeasible and what makes a
  direction feasible and descending;
- set up and solve the linear subproblem of the Frank–Wolfe method, use its
  optimal value as an optimality test and a lower bound, and carry out the
  line search on $[0,1]$;
- state and prove the $O(1/k)$ convergence bound of Frank–Wolfe and
  recognize its zig-zag between vertices;
- build the quadratic penalty function of a nonlinear program (NLP), compute the minimizers
  $\vx(\rho)$ and follow them as $\rho\to\infty$;
- estimate KKT multipliers from penalty solutions and explain the
  ill-conditioning of large penalties; and
- describe the logarithmic barrier and the augmented Lagrangian and how they
  differ from the quadratic penalty.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\vy^{(k)}$ | solution of the Frank–Wolfe linear subproblem at iteration $k$ |
| $\Gamma(\vx)$ | Frank–Wolfe gap, $\max_{\vy\in\Omega}\nabla f(\vx)^\T(\vx-\vy)\ge0$ |
| $D$, $L_f$ | diameter of $\Omega$; Lipschitz constant of $\nabla f$ |
| $\psi(\vx)$ | penalty function (measure of infeasibility), $\psi(\vx)=0$ exactly when $\vx\in\Omega$ |
| $\rho$, $\rho_k$ | penalty parameter |
| $\gamma$ | growth factor of the penalty parameter, $\rho_{k+1}=\gamma\rho_k$ |
| $\Phi_\rho(\vx)$ | penalized objective $f(\vx)+\rho\psi(\vx)$ |
| $\vx(\rho)$ | a minimizer of $\Phi_\rho$ |
| $\ell$ | number of equality constraints $h_j(\vx)=0$ |
| $\eta$ | barrier parameter |
| $L_A(\vx,\vmu;\rho)$ | augmented Lagrangian |

(sec-06-06-two-strategies)=
## Two strategies for constrained problems

The problem is

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)\\
\st\ & g_i(\vx)\le0,\quad i=1,\dots,m,\\
& h_j(\vx)=0,\quad j=1,\dots,\ell,
\end{aligned}
$$ (eq-06-06-nlp)

with feasible set $\Omega$. The goal of an algorithm is a local minimizer or,
more modestly, a KKT point ({doc}`06-03-kkt-conditions`). Two ideas are used.

1. **Feasible-direction methods.** Start from a feasible point and move only
   along directions $\vd$ that are **feasible** ($\vx+\alpha\vd\in\Omega$ for all
   small $\alpha>0$, {prf:ref}`def-06-02-feasible-direction`) and **descent**
   ($\nabla f(\vx)^\T\vd<0$). The steepest-descent direction $-\nabla f(\vx)$
   is the natural candidate, but at a point on the boundary it often points
   out of $\Omega$; the method must find a direction that is a compromise
   between descent and feasibility. Every iterate is feasible, so the method
   can be stopped at any time with a usable solution.
2. **Transformation methods.** Replace {eq}`eq-06-06-nlp` by a sequence of
   unconstrained problems whose solutions approach a solution of
   {eq}`eq-06-06-nlp`, by adding to $f$ a term that punishes infeasibility
   (penalty methods) or that blocks the way to the boundary (barrier methods),
   and solve them with the methods of {doc}`06-05-unconstrained-algorithms`.

(sec-06-06-frank-wolfe)=
## The Frank–Wolfe method

Frank and Wolfe (1956) proposed the first idea for problems with *linear*
constraints. Let

$$
\begin{aligned}
&\min_{\vx\in\Omega}f(\vx),\\
&\Omega=\{\vx\in\R^n:\mA\vx=\vb,\ \vx\ge\vzero\},
\end{aligned}
$$ (eq-06-06-fw-problem)

where $\Omega$ is nonempty and bounded (a polytope) and $f$ is convex and
continuously differentiable. Inequalities $\mA\vx\le\vb$ are allowed as well
(add slack variables, {doc}`03-02-lp-geometry`).

**Feasible directions from convexity.** Since $\Omega$ is convex, for any
$\vx^{(k)}\in\Omega$ and any $\vy\in\Omega$ the whole segment between them is
feasible:

$$
\begin{aligned}
\vx^{(k)}+\alpha(\vy-\vx^{(k)})&=(1-\alpha)\vx^{(k)}+\alpha\vy\\
&\in\Omega\qquad\text{for all }\alpha\in[0,1].
\end{aligned}
$$

So every $\vd=\vy-\vx^{(k)}$ with $\vy\in\Omega$ is a feasible direction, and
the step length is restricted to $[0,1]$.

**Descent from a linear program.** Among these directions we want the one
with the most negative directional derivative $\nabla f(\vx^{(k)})^\T(\vy-\vx^{(k)})$.
This is the **linearized subproblem**

$$
\begin{aligned}
\min_{\vy}\ & \nabla f(\vx^{(k)})^\T\vy\\
\st\ & \mA\vy=\vb,\\
& \vy\ge\vzero,
\end{aligned}
$$ (eq-06-06-fw-lp)

which has the same minimizers as $\min_{\vy\in\Omega}\nabla f(\vx^{(k)})^\T(\vy-\vx^{(k)})$
because the two objectives differ by the constant $\nabla f(\vx^{(k)})^\T\vx^{(k)}$.
Problem {eq}`eq-06-06-fw-lp` is an LP: it minimizes the first-order Taylor
approximation of $f$ over $\Omega$. Because $\Omega$ is a nonempty polytope,
it has an optimal solution at a vertex ({prf:ref}`thm-03-02-fundamental`), which the
simplex method finds ({prf:ref}`alg-03-03-simplex`). Call it $\vy^{(k)}$. The
optimal value measures how much descent is available:

:::{prf:definition} Frank–Wolfe gap
:label: def-06-06-fw-gap

For $\vx\in\Omega$, the **Frank–Wolfe gap** is

$$
\begin{aligned}
\Gamma(\vx)&=\max_{\vy\in\Omega}\nabla f(\vx)^\T(\vx-\vy)\\
&=-\min_{\vy\in\Omega}\nabla f(\vx)^\T(\vy-\vx)\ \ge0 .
\end{aligned}
$$

It is nonnegative because $\vy=\vx$ is allowed. At iteration $k$,
$\Gamma(\vx^{(k)})=-\nabla f(\vx^{(k)})^\T(\vy^{(k)}-\vx^{(k)})$.
:::

If $\Gamma(\vx^{(k)})>0$, then $\vd^{(k)}=\vy^{(k)}-\vx^{(k)}$ is a feasible
descent direction. If $\Gamma(\vx^{(k)})=0$, no feasible direction descends to
first order, and for a convex $f$ the current point is optimal. The next
theorem makes this precise and answers why a zero gap certifies optimality.

:::{prf:theorem} Optimality test and lower bound
:label: thm-06-06-fw-gap

Let $f$ be convex and differentiable on the convex set $\Omega$, and
$f^*=\min_{\vx\in\Omega}f(\vx)$. For every $\vx\in\Omega$,

$$
f(\vx)-\Gamma(\vx)\le f^*\le f(\vx).
$$ (eq-06-06-fw-bound)

In particular, $\vx$ is optimal if and only if $\Gamma(\vx)=0$, that is,
$\nabla f(\vx)^\T(\vy-\vx)\ge0$ for all $\vy\in\Omega$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

By the first-order characterization of convexity ({prf:ref}`thm-02-02-first-order`),
$f(\vy)\ge f(\vx)+\nabla f(\vx)^\T(\vy-\vx)$ for all $\vy\in\Omega$. Minimizing
both sides over $\vy\in\Omega$ gives
$f^*\ge f(\vx)+\min_{\vy\in\Omega}\nabla f(\vx)^\T(\vy-\vx)=f(\vx)-\Gamma(\vx)$.
The upper bound holds because $\vx$ is feasible. If $\Gamma(\vx)=0$, the two
bounds coincide and $f(\vx)=f^*$. Conversely, if $\Gamma(\vx)>0$, the direction
$\vy-\vx$ attaining the maximum is a feasible descent direction, so $f$
decreases along the segment and $\vx$ is not optimal.
:::

The condition $\nabla f(\vx)^\T(\vy-\vx)\ge0$ for all $\vy\in\Omega$ is the
first-order optimality condition over a convex set
({prf:ref}`thm-02-03-optimality`, {prf:ref}`thm-06-02-first-order`); for polyhedral
$\Omega$ it is equivalent to the KKT conditions ({prf:ref}`thm-06-03-kkt`), since the LP {eq}`eq-06-06-fw-lp`
then has the optimal solution $\vy=\vx$ and LP duality supplies the
multipliers. Without convexity, $\Gamma(\vx)=0$ still identifies a KKT point.
The bound {eq}`eq-06-06-fw-bound` is a by-product of the LP solve: every
iteration yields a **certified lower bound** $f(\vx^{(k)})-\Gamma(\vx^{(k)})$,
and the best of these bounds tells how far the current value can be from the
optimum.

**Line search.** The new point is the best point on the segment,

$$
\begin{aligned}
\alpha_k&\in\argmin_{0\le\alpha\le1}f\bigl(\vx^{(k)}+\alpha\vd^{(k)}\bigr),\\
\vx^{(k+1)}&=\vx^{(k)}+\alpha_k\vd^{(k)},
\end{aligned}
$$ (eq-06-06-fw-line-search)

computed by any method of {doc}`06-04-line-search`. For a quadratic $f$ with
Hessian $\mQ$, $\varphi(\alpha)=f(\vx^{(k)})-\alpha\Gamma(\vx^{(k)})+\tfrac12\alpha^2\vd^{(k)\T}\mQ\vd^{(k)}$,
and the exact step is

$$
\alpha_k=\min\Bigl\{1,\ \frac{\Gamma(\vx^{(k)})}{\vd^{(k)\T}\mQ\vd^{(k)}}\Bigr\}
$$ (eq-06-06-fw-quadratic-step)

(with $\alpha_k=1$ if $\vd^{(k)\T}\mQ\vd^{(k)}=0$). The upper bound $1$
matters: beyond $\alpha=1$ the point leaves $\Omega$ in general.

:::{prf:algorithm} Frank–Wolfe method (convex-combination method)
:label: alg-06-06-frank-wolfe

**Input:** a convex, continuously differentiable $f$; a nonempty polytope
$\Omega$; a feasible $\vx^{(0)}\in\Omega$; a tolerance $\varepsilon>0$.
**Output:** a feasible $\vx$ with $f(\vx)-f^*\le\varepsilon$.

1. **Initialize.** Set $k=0$.
2. **Linear subproblem.** Solve the LP $\min_{\vy\in\Omega}\nabla f(\vx^{(k)})^\T\vy$;
   let $\vy^{(k)}$ be an optimal vertex and $\vd^{(k)}=\vy^{(k)}-\vx^{(k)}$.
3. **Stopping test.** If
   $\Gamma(\vx^{(k)})=-\nabla f(\vx^{(k)})^\T\vd^{(k)}\le\varepsilon$, stop: by
   {eq}`eq-06-06-fw-bound`, $f(\vx^{(k)})-f^*\le\varepsilon$.
4. **Line search.** Find $\alpha_k\in[0,1]$ minimizing
   $f(\vx^{(k)}+\alpha\vd^{(k)})$.
5. **Update.** Set $\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$,
   $k\leftarrow k+1$, and go to Step 2.
:::

Each iterate is a convex combination of $\vx^{(0)}$ and the vertices
$\vy^{(0)},\dots,\vy^{(k)}$, hence the name **convex-combination method**;
it is also called the **linear-approximation** or **conditional-gradient**
method. A feasible starting point is itself the solution of an LP (for
example $\vx^{(0)}=\vy$ for the vertex returned by phase I of the simplex
method, {doc}`03-04-two-phase-big-m`).

:::{prf:example} Frank–Wolfe on a quadratic over a triangle
:label: ex-06-06-frank-wolfe

Solve

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)=-5x_1+x_1^2-8x_2+2x_2^2\\
\st\ & 3x_1+2x_2\le6,\\
& x_1,x_2\ge0,
\end{aligned}
$$

by the Frank–Wolfe method from $\vx^{(0)}=(0,0)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The gradient is $\nabla f(\vx)=(2x_1-5,\ 4x_2-8)$ and the Hessian
$\mQ=\diag(2,4)$. The feasible set is the triangle with vertices $(0,0)$,
$(2,0)$ and $(0,3)$, so each LP is solved by comparing the objective at the
three vertices. The unconstrained minimizer $(2.5,2)$ violates
$3x_1+2x_2\le6$ ($3\cdot2.5+2\cdot2=11.5$), so the constraint will be active
at the optimum.

**Iteration 0.** $\nabla f(\vx^{(0)})=(-5,-8)$. The LP $\min\,-5y_1-8y_2$
over the triangle has the values $0$, $-10$, $-24$ at the vertices, so
$\vy^{(0)}=(0,3)$, $\vd^{(0)}=(0,3)$ and $\Gamma(\vx^{(0)})=24$. Along the
segment, $\vx^{(0)}+\alpha\vd^{(0)}=(0,3\alpha)$ and

$$
\begin{aligned}
\varphi(\alpha)&=-24\alpha+18\alpha^2,\\
\varphi'(\alpha)&=-24+36\alpha=0\ \Rightarrow\ \alpha_0=\tfrac23 ,
\end{aligned}
$$

which agrees with {eq}`eq-06-06-fw-quadratic-step`: $24/(\vd^\T\mQ\vd)=24/36$. So
$\vx^{(1)}=(0,2)$ and $f(\vx^{(1)})=-8$; the lower bound is $0-24=-24$.

**Iteration 1.** $\nabla f(\vx^{(1)})=(-5,0)$. The LP $\min\,-5y_1$ gives
$\vy^{(1)}=(2,0)$, $\vd^{(1)}=(2,-2)$ and $\Gamma(\vx^{(1)})=10$. Along the
segment, $\vx^{(1)}+\alpha\vd^{(1)}=(2\alpha,\ 2-2\alpha)$ and

$$
\begin{aligned}
\varphi(\alpha)&=-8-10\alpha+12\alpha^2,\\
\varphi'(\alpha)&=-10+24\alpha=0\ \Rightarrow\ \alpha_1=\tfrac{5}{12},
\end{aligned}
$$

so $\vx^{(2)}=(0,2)+\tfrac5{12}(2,-2)=(\tfrac56,\tfrac76)$ and
$f(\vx^{(2)})=-\tfrac{121}{12}=-10.0833$; the lower bound improves to
$-8-10=-18$. *(This line function is sometimes printed as
$8+10\alpha-12\alpha^2$ with derivative $10-24\alpha$, which would make
$\alpha_1$ a maximizer; the signs are as above, and the step
$\alpha_1=\tfrac5{12}$ is a minimizer.)*

**Iteration 2.** $\nabla f(\vx^{(2)})=(-\tfrac{10}3,-\tfrac{10}3)$. The vertex
values are $0$, $-\tfrac{20}3$ and $-10$, so $\vy^{(2)}=(0,3)$,
$\vd^{(2)}=(-\tfrac56,\tfrac{11}6)$, $\Gamma(\vx^{(2)})=\tfrac{10}3$, and
{eq}`eq-06-06-fw-quadratic-step` gives $\alpha_2=\tfrac{10/3}{89/6}=\tfrac{20}{89}=0.2247$.
Hence $\vx^{(3)}=(\tfrac{115}{178},\tfrac{281}{178})=(0.6461,\ 1.5787)$ with
$f(\vx^{(3)})=-10.4579$.

**Iteration 3.** $\nabla f(\vx^{(3)})=(-3.7079,\ -1.6854)$; now $(2,0)$ wins
the LP, $\Gamma=2.3596$, $\alpha_3=0.1731$, and $\vx^{(4)}=(0.8804,\ 1.3055)$
with $f(\vx^{(4)})=-10.6620$.

{numref}`tab-06-06-frank-wolfe` continues the run. The optimum follows from the
KKT conditions with the constraint $3x_1+2x_2\le6$ active:
$2x_1-5+3\lambda=0$, $4x_2-8+2\lambda=0$, $3x_1+2x_2=6$ give $\lambda=1\ge0$,
$\vx^*=(1,\ 1.5)$ and $f^*=-11.5$.

```{list-table} Frank–Wolfe iterations for $\min\,-5x_1+x_1^2-8x_2+2x_2^2$ over the triangle (4 significant digits)
:name: tab-06-06-frank-wolfe
:header-rows: 1

* - $k$
  - $\vx^{(k)}$
  - $f(\vx^{(k)})$
  - $\vy^{(k)}$
  - $\Gamma(\vx^{(k)})$
  - $\alpha_k$
  - best lower bound
* - 0
  - $(0,\ 0)$
  - $0$
  - $(0,3)$
  - $24$
  - $0.6667$
  - $-24$
* - 1
  - $(0,\ 2)$
  - $-8$
  - $(2,0)$
  - $10$
  - $0.4167$
  - $-18$
* - 2
  - $(0.8333,\ 1.167)$
  - $-10.08$
  - $(0,3)$
  - $3.333$
  - $0.2247$
  - $-13.42$
* - 3
  - $(0.6461,\ 1.579)$
  - $-10.46$
  - $(2,0)$
  - $2.360$
  - $0.1731$
  - $-12.82$
* - 4
  - $(0.8804,\ 1.305)$
  - $-10.66$
  - $(0,3)$
  - $1.856$
  - $0.1424$
  - $-12.52$
* - 5
  - $(0.7550,\ 1.547)$
  - $-10.79$
  - $(2,0)$
  - $1.541$
  - $0.1216$
  - $-12.33$
* - 6
  - $(0.9064,\ 1.359)$
  - $-10.89$
  - $(0,3)$
  - $1.322$
  - $0.1064$
  - $-12.21$
* - 7
  - $(0.8099,\ 1.533)$
  - $-10.96$
  - $(2,0)$
  - $1.160$
  - $0.09483$
  - $-12.12$
* - 8
  - $(0.9228,\ 1.388)$
  - $-11.01$
  - $(0,3)$
  - $1.036$
  - $0.08562$
  - $-12.05$
* - 9
  - $(0.8438,\ 1.526)$
  - $-11.06$
  - $(2,0)$
  - $0.9365$
  - $0.07812$
  - $-11.99$
```

From iteration 2 on, the LP alternates between the vertices $(0,3)$ and
$(2,0)$, the two ends of the active edge, and the iterates zig-zag towards
$\vx^*=(1,1.5)$ with shrinking steps ({numref}`fig-06-06-frank-wolfe`). Progress
is slow: the gap falls below $1$ at $k=9$, below $0.1$ only at $k=105$, and
$f(\vx^{(k)})-f^*\le0.01$ needs $k=542$. Yet the bounds are always valid: at
$k=9$ the optimal value is known to lie in $[-11.99,\ -11.06]$.

```{figure} figures/06-06-frank-wolfe.svg
:name: fig-06-06-frank-wolfe
:width: 560px
:class: dark-light
:alt: Feasible triangle with elliptic contours of the objective centered outside it, and a path of Frank–Wolfe iterates that first climbs the x2 axis and then zig-zags towards the optimum on the hypotenuse

Frank–Wolfe on the triangle of {prf:ref}`ex-06-06-frank-wolfe`. Contours of $f$
(center $(2.5,2)$ outside $\Omega$), iterates $\vx^{(0)},\dots,\vx^{(9)}$ and
the optimum $\vx^*=(1,1.5)$ on the edge $3x_1+2x_2=6$. Each step moves towards
a vertex $\vy^{(k)}$ (dotted), alternately $(0,3)$ and $(2,0)$, so the path
zig-zags and the steps shrink.
```
:::

**Convergence.** The zig-zag seen in the solution of
{prf:ref}`ex-06-06-frank-wolfe` is not an accident of the example: when the
optimum lies on a face of $\Omega$ and not at a vertex, every LP solution is a
vertex, so the steps can never point straight at $\vx^*$. The rate is
sublinear.

:::{prf:theorem} Convergence rate of Frank–Wolfe
:label: thm-06-06-fw-rate

Let $\Omega$ be compact and convex with diameter
$D=\max_{\vx,\vy\in\Omega}\norm{\vx-\vy}$, and let $f$ be convex with a
Lipschitz continuous gradient, $\norm{\nabla f(\vx)-\nabla f(\vy)}\le L_f\norm{\vx-\vy}$.
With exact line search (or with the fixed steps $\alpha_k=2/(k+2)$), the
Frank–Wolfe iterates satisfy

$$
f(\vx^{(k)})-f^*\le\frac{2L_fD^2}{k+2},\qquad k\ge1 .
$$ (eq-06-06-fw-rate)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $h_k=f(\vx^{(k)})-f^*$. The Lipschitz gradient gives the quadratic upper
bound $f(\vx+\alpha\vd)\le f(\vx)+\alpha\nabla f(\vx)^\T\vd+\tfrac12L_f\alpha^2\norm{\vd}^2$.
With $\vd=\vd^{(k)}$, $\nabla f(\vx^{(k)})^\T\vd^{(k)}=-\Gamma(\vx^{(k)})\le-h_k$
by {eq}`eq-06-06-fw-bound`, and $\norm{\vd^{(k)}}\le D$:

$$
h_{k+1}\le(1-\alpha)h_k+\tfrac12L_f D^2\alpha^2\qquad\text{for every }\alpha\in[0,1],
$$

because the exact line search does at least as well as any fixed
$\alpha\in[0,1]$. For $k=0$ take $\alpha=1$: $h_1\le\tfrac12L_fD^2\le\tfrac23L_fD^2$.
If $h_k\le2L_fD^2/(k+2)$, take $\alpha=2/(k+2)$:

$$
\begin{aligned}
h_{k+1}&\le\frac{k}{k+2}\cdot\frac{2L_fD^2}{k+2}+\frac{2L_fD^2}{(k+2)^2}\\
&=\frac{2L_fD^2(k+1)}{(k+2)^2}\le\frac{2L_fD^2}{k+3},
\end{aligned}
$$

since $(k+1)(k+3)\le(k+2)^2$.
:::

In {prf:ref}`ex-06-06-frank-wolfe`, $L_f=4$ (largest eigenvalue of $\mQ$) and
$D=\sqrt{13}$ (the edge from $(2,0)$ to $(0,3)$), so the theorem guarantees
$f(\vx^{(k)})-f^*\le104/(k+2)$; the observed error is about $5/k$, much better
than the bound but of the same order $1/k$. Variants with **away steps**
(which may move away from a vertex already in the convex combination)
converge linearly on polytopes when $f$ is strongly convex.

**When is Frank–Wolfe the method of choice?** Its strength is that the
subproblem is an LP, and often a very special one. In **network flow
problems with convex arc costs**, such as traffic assignment on road networks
(Sheffi, 1985), the LP {eq}`eq-06-06-fw-lp` over the flow-conservation
constraints decomposes into **shortest-path problems** with the marginal arc
costs $\nabla f(\vx^{(k)})$ as lengths, and its solution loads all demand on
the shortest paths ("all-or-nothing" assignment). Each iteration then needs
one shortest-path computation (for example {prf:ref}`alg-04-02-dijkstra`) and no matrix at
all, which makes the method practical for networks with many thousands of
links ({prf:ref}`ex-08-09-fw-network-flow`). The sparse iterates (convex
combinations of few vertices) are also attractive in machine learning.

(sec-06-06-other-feasible-direction)=
### Other feasible-direction methods

Frank–Wolfe chooses the direction by the extreme point that is best to first
order, which wastes the information that the current point is already near
the solution. Several feasible-direction methods use the constraints more
cleverly (Bazaraa et al., 2006, Ch. 10; Luenberger & Ye, 2016, Ch. 12).

- **Zoutendijk's method** solves, at each iteration, a small LP for the
  direction $\vd$ that maximizes the rate of descent subject to staying
  feasible for the active constraints and a normalization such as
  $-1\le d_j\le1$ (Zoutendijk, 1960); the direction is then no longer forced towards a vertex.
- **Gradient projection** (Rosen, 1960) projects $-\nabla f(\vx)$ onto the subspace
  of the active linear constraints and moves along the projection; when the
  projection is zero, the signs of the multiplier estimates decide which
  constraint to release.
- **Reduced-gradient methods** (Wolfe, 1963; generalized for nonlinear constraints)
  partition the variables into basic and nonbasic ones as in the simplex
  method ({doc}`03-03-simplex-method`), eliminate the basic variables through the
  equality constraints and do steepest descent or quasi-Newton steps in the
  nonbasic ones.
- **Sequential quadratic programming** (SQP; Wilson, 1963) replaces the linear subproblem
  {eq}`eq-06-06-fw-lp` by a quadratic one that uses second-order information of
  the Lagrangian; it is Newton's method applied to the KKT conditions and
  converges superlinearly (Nocedal & Wright, 2006, Ch. 18).

All of them keep the Frank–Wolfe principle of solving a simpler constrained
subproblem at each iteration, and they all reduce to the methods of
{doc}`06-05-unconstrained-algorithms` when no constraint is active.

(sec-06-06-penalty)=
## Penalty methods

The second strategy moves the constraints into the objective. Consider
first a single inequality, $\min f(\vx)$ s.t. $g(\vx)\le0$. The idea is to
*allow* violations but to charge for them: no charge when $g(\vx)\le0$, a
charge growing with $g(\vx)$ when $g(\vx)>0$. The simplest such charge,
$\rho\max\{0,g(\vx)\}$, gives the **exact penalty** problem

$$
\min_{\vx}\ f(\vx)+\rho\max\{0,g(\vx)\}.
$$ (eq-06-06-exact-penalty)

It is called exact because, if $\vx^*$ is a strict local minimizer of the
constrained problem that satisfies the KKT conditions
({prf:ref}`thm-06-03-kkt`) with multiplier $\lambda^*$, then $\vx^*$ is also a
local minimizer of {eq}`eq-06-06-exact-penalty` for every
$\rho>\lambda^*$. Its drawback is the kink at $g(\vx)=0$, exactly where the
solution usually lies, so the smooth methods of
{doc}`06-05-unconstrained-algorithms` cannot be applied. The classical remedy
squares the violation.

:::{prf:definition} Quadratic penalty function
:label: def-06-06-quadratic-penalty

For problem {eq}`eq-06-06-nlp`, the **quadratic penalty function** is

$$
\psi(\vx)=\sum_{i=1}^{m}\bigl[\max\{0,g_i(\vx)\}\bigr]^2+\sum_{j=1}^{\ell}h_j(\vx)^2 ,
$$ (eq-06-06-quadratic-penalty)

and the **penalized objective** with parameter $\rho>0$ is
$\Phi_\rho(\vx)=f(\vx)+\rho\psi(\vx)$.
:::

The penalty $\psi$ is continuous, nonnegative, and zero exactly on $\Omega$.
If $g_i$ and $h_j$ are continuously differentiable, so is $\psi$: the function
$t\mapsto\max\{0,t\}^2$ has derivative $2\max\{0,t\}$, which is continuous at
$t=0$. More generally one can use
$\psi(\vx)=\sum_i\max\{0,g_i(\vx)\}^r+\sum_j\abs{h_j(\vx)}^r$ with an exponent
$r\ge1$: $r=1$ gives the exact (nonsmooth) penalty, and every $r>1$ gives a
continuously differentiable $\psi$. The exponent does not have to be even,
because the maximum and the absolute value already make each term
nonnegative. The quadratic penalty $r=2$ is the usual choice. It goes back to
Courant (1943); its systematic use as a sequence of unconstrained problems is
due to Fiacco and McCormick (1968).

:::{prf:example} Quadratic penalty for a one-dimensional problem
:label: ex-06-06-penalty-1d

Solve $\min x$ subject to $x\ge2$ by the quadratic penalty method.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Write the constraint as $g(x)=2-x\le0$. The penalty is
$\psi(x)=\max\{0,2-x\}^2$, equal to $0$ for $x\ge2$ and $(2-x)^2$ for $x<2$,
and

$$
\Phi_\rho(x)=x+\rho\max\{0,2-x\}^2 .
$$

On $x\ge2$, $\Phi_\rho(x)=x$ is increasing, so a minimizer must satisfy
$x\le2$. For $x<2$, $\Phi_\rho'(x)=1-2\rho(2-x)=0$ gives

$$
x(\rho)=2-\frac{1}{2\rho},
$$

and $\Phi_\rho''(x)=2\rho>0$ confirms a minimum. *(The derivative is sometimes
printed as $1-2\rho(2-x)^2$; the square does not belong there, and only the
correct derivative $1-2\rho(2-x)$ gives the formula for $x(\rho)$.)*

```{list-table} Quadratic penalty for $\min x$ s.t. $x\ge2$
:name: tab-06-06-penalty-1d
:header-rows: 1

* - $\rho$
  - $x(\rho)$
  - $\psi(x(\rho))$
  - $\rho\,\psi(x(\rho))$
  - $\Phi_\rho(x(\rho))$
  - $2\rho\max\{0,g(x(\rho))\}$
* - 1
  - 1.5
  - 0.25
  - 0.25
  - 1.75
  - 1
* - 10
  - 1.95
  - $2.5\times10^{-3}$
  - 0.025
  - 1.975
  - 1
* - 100
  - 1.995
  - $2.5\times10^{-5}$
  - $2.5\times10^{-3}$
  - 1.9975
  - 1
* - 1000
  - 1.9995
  - $2.5\times10^{-7}$
  - $2.5\times10^{-4}$
  - 1.99975
  - 1
```

As $\rho\to\infty$, $x(\rho)\to2=x^*$ **from outside** the feasible set: every
$x(\rho)$ is slightly infeasible, with violation $1/(2\rho)$. The penalized
values $\Phi_\rho(x(\rho))=2-1/(4\rho)$ increase towards $f^*=2$ and are lower
bounds on it. The last column is constant and equal to the KKT multiplier
of the problem: $1-\lambda=0$ gives $\lambda^*=1$. {numref}`fig-06-06-penalty-family`
shows the penalized functions.

```{figure} figures/06-06-penalty-family.svg
:name: fig-06-06-penalty-family
:width: 520px
:class: dark-light
:alt: Several convex curves x plus rho times the squared violation for increasing rho, all equal to the line f of x equals x on the feasible side, with their minimizers moving right towards x equals 2

Penalized objectives $\Phi_\rho(x)=x+\rho\max\{0,2-x\}^2$ of
{prf:ref}`ex-06-06-penalty-1d` for $\rho=0.5,1,2,5$. On the feasible side
$x\ge2$ (shaded) they all coincide with $f(x)=x$. The minimizers
$x(\rho)=2-1/(2\rho)=1$, $1.5$, $1.75$, $1.9$ (dots) lie outside the feasible set
and approach $x^*=2$ as $\rho$ grows, while the curves become steeper.
```
:::

For the problem of {prf:ref}`ex-06-06-penalty-1d` the exact penalty
{eq}`eq-06-06-exact-penalty`,
$x+\rho\max\{0,2-x\}$, behaves differently: for $\rho<1$ it is unbounded below
($(1-\rho)x+2\rho\to-\infty$ as $x\to-\infty$), for $\rho=1$ every $x\le2$ is a
minimizer, and for every $\rho>1=\lambda^*$ its unique minimizer is exactly
$x^*=2$. One finite $\rho$ suffices, at the price of non-smoothness.

:::{prf:algorithm} Quadratic penalty method
:label: alg-06-06-penalty

**Input:** problem {eq}`eq-06-06-nlp`; an initial parameter $\rho_0>0$; a growth
factor $\gamma>1$ (for example $\gamma=10$); a starting point
$\vx^{(-1)}$; a tolerance $\varepsilon>0$.
**Output:** an approximate solution and multiplier estimates.

1. **Initialize.** Set $k=0$.
2. **Unconstrained minimization.** Starting from $\vx^{(k-1)}$, compute a
   minimizer $\vx^{(k)}=\vx(\rho_k)$ of $\Phi_{\rho_k}(\vx)=f(\vx)+\rho_k\psi(\vx)$
   with a method of {doc}`06-05-unconstrained-algorithms`.
3. **Stopping test.** If $\rho_k\psi(\vx^{(k)})\le\varepsilon$ (or the
   constraint violation is below a tolerance), stop.
4. **Increase the penalty.** Set $\rho_{k+1}=\gamma\rho_k$,
   $k\leftarrow k+1$, and go to Step 2.
:::

Starting each minimization from the previous solution (**warm start**) is
essential: $\vx(\rho_k)$ is close to $\vx(\rho_{k+1})$, whereas a cold start
on a problem with a large $\rho$ is difficult.

**Choosing the parameters.** The initial value $\rho_0$ should make the
penalty term comparable with the objective, for example
$\rho_0\psi(\vx^{(-1)})\approx\abs{f(\vx^{(-1)})}$ at the starting point. A
large growth factor ($\gamma=10$) needs few rounds but makes each warm start
less accurate; a small one ($\gamma=2$) needs more rounds of easier
minimizations. Because the violation of an active constraint behaves like
$\lambda_i^*/(2\rho)$ ({ref}`sec-06-06-multipliers`), reaching a violation of
$10^{-6}$ with a multiplier of order $1$ needs $\rho\approx10^{6}$, which is
about where the ill-conditioning described below begins to cost digits in
double-precision arithmetic. Each subproblem need not be solved exactly
either: a gradient tolerance that decreases with $k$ is enough for
convergence. The basic convergence result needs only continuity.

:::{prf:theorem} Convergence of the penalty method
:label: thm-06-06-penalty-convergence

Let $f$, $g_i$ and $h_j$ be continuous, let {eq}`eq-06-06-nlp` have an optimal
solution with value $f^*$, let $0<\rho_k<\rho_{k+1}\to\infty$, and let each
$\vx^{(k)}$ be a global minimizer of $\Phi_{\rho_k}$. Then

1. $\Phi_{\rho_k}(\vx^{(k)})\le\Phi_{\rho_{k+1}}(\vx^{(k+1)})\le f^*$,
   $\psi(\vx^{(k+1)})\le\psi(\vx^{(k)})$ and $f(\vx^{(k)})\le f(\vx^{(k+1)})$;
2. every limit point of $\{\vx^{(k)}\}$ is an optimal solution of
   {eq}`eq-06-06-nlp`.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

(1) Write $\vx_k=\vx^{(k)}$, $\psi_k=\psi(\vx_k)$, $f_k=f(\vx_k)$. For a
feasible optimal $\vx^*$, $\Phi_{\rho_k}(\vx_k)\le\Phi_{\rho_k}(\vx^*)=f^*$. Since
$\rho_{k+1}>\rho_k$ and $\psi\ge0$,
$\Phi_{\rho_k}(\vx_k)\le\Phi_{\rho_k}(\vx_{k+1})\le\Phi_{\rho_{k+1}}(\vx_{k+1})$.
Adding the optimality inequalities $f_k+\rho_k\psi_k\le f_{k+1}+\rho_k\psi_{k+1}$
and $f_{k+1}+\rho_{k+1}\psi_{k+1}\le f_k+\rho_{k+1}\psi_k$ gives
$(\rho_{k+1}-\rho_k)(\psi_{k+1}-\psi_k)\le0$, so $\psi_{k+1}\le\psi_k$, and then
the first inequality gives $f_k\le f_{k+1}$.
(2) Let $\bar\vx$ be the limit of a subsequence $\vx_k$, $k\in\mathcal K$.
From $f_k+\rho_k\psi_k\le f^*$ and $f_k\to f(\bar\vx)$ along $\mathcal K$,
$\rho_k\psi_k$ stays bounded, so $\psi_k\to0$ and, by continuity,
$\psi(\bar\vx)=0$: $\bar\vx$ is feasible. Moreover $f(\bar\vx)=\lim f_k\le f^*$,
so $\bar\vx$ is optimal.
:::

The theorem assumes global minimizers of $\Phi_{\rho_k}$; with local
minimizers the same argument applies near an isolated local solution.
Part (1) shows that the values $\Phi_{\rho_k}(\vx^{(k)})$ are lower bounds that
increase towards $f^*$, as in {numref}`tab-06-06-penalty-1d`.

(sec-06-06-multipliers)=
### Multiplier estimates

At an unconstrained minimizer of $\Phi_\rho$ with differentiable data,

$$
\begin{aligned}
&\nabla f(\vx(\rho))+\sum_{i=1}^{m}2\rho\max\{0,g_i(\vx(\rho))\}\nabla g_i(\vx(\rho))\\
&\quad+\sum_{j=1}^{\ell}2\rho\,h_j(\vx(\rho))\nabla h_j(\vx(\rho))=\vzero .
\end{aligned}
$$ (eq-06-06-penalty-stationarity)

Compare with the KKT stationarity condition
$\nabla f+\sum_i\lambda_i\nabla g_i+\sum_j\mu_j\nabla h_j=\vzero$
({prf:ref}`thm-06-03-kkt`): the quantities

$$
\begin{aligned}
\lambda_i(\rho)&=2\rho\max\{0,g_i(\vx(\rho))\}\ge0,\\
\mu_j(\rho)&=2\rho\,h_j(\vx(\rho))
\end{aligned}
$$ (eq-06-06-multiplier-estimates)

play the role of the multipliers, and under a constraint qualification they
converge to the KKT multipliers $\lambda_i^*$, $\mu_j^*$ as $\rho\to\infty$
(Nocedal & Wright, 2006, Thm. 17.2). Inactive constraints give
$\lambda_i(\rho)=0$ once $\vx(\rho)$ is close to $\vx^*$, which is complementary
slackness. In {prf:ref}`ex-06-06-penalty-1d` the estimate is exactly $1$ for
every $\rho$. The estimates also explain *why* $\vx(\rho)$ is infeasible:
a nonzero multiplier requires a violation of size $\lambda_i^*/(2\rho)$.

(sec-06-06-ill-conditioning)=
### Ill-conditioning

The price of large $\rho$ is numerical. Near $\vx(\rho)$ the Hessian of
$\Phi_\rho$ is approximately

$$
\begin{aligned}
\nabla^2\Phi_\rho&\approx\nabla^2_{\vx\vx}L(\vx,\vlambda,\vmu)\\
&\quad+2\rho\sum_{i\in I}\nabla g_i\nabla g_i^\T+2\rho\sum_{j}\nabla h_j\nabla h_j^\T ,
\end{aligned}
$$

the Hessian of the Lagrangian plus a term of size $\rho$ in the directions of
the active constraint gradients. Some eigenvalues grow like $\rho$ and the
others stay bounded, so the condition number grows linearly in $\rho$. In
{prf:ref}`ex-06-06-penalty-1d` the curvature $2\rho$ of $\Phi_\rho$ grows without
bound; in the two-variable problem of {prf:ref}`ex-08-09-penalty-parabola` the
condition number of $\nabla^2\Phi_\rho$ grows from $3.7$ at $\rho=0.1$ to about
$1.2\times10^4$ at $\rho=10^4$. By {prf:ref}`thm-06-05-sd-rate` steepest descent
on such a function is hopelessly slow, and even Newton's method suffers from
rounding errors. This is why the penalty is increased gradually with warm
starts, and why the methods of the next two sections were developed.

(sec-06-06-barrier)=
## Barrier methods

A **barrier** (interior) method approaches the solution from *inside* the
feasible set. For inequality constraints with a nonempty interior
$\{\vx:g_i(\vx)<0\}$ it minimizes

$$
f(\vx)-\eta\sum_{i=1}^{m}\ln\bigl(-g_i(\vx)\bigr),\qquad \eta>0,
$$ (eq-06-06-log-barrier)

over the interior. The logarithmic barrier tends to $+\infty$ at the boundary,
so the iterates stay strictly feasible, and the barrier parameter $\eta$ is
decreased to zero. For $\min x$ s.t. $x\ge2$, the barrier function
$x-\eta\ln(x-2)$ has the minimizer $x(\eta)=2+\eta$, which approaches $x^*=2$
from the feasible side; the quantities $\eta/(-g_i(\vx(\eta)))$, here
$\eta/\eta=1$, estimate the multipliers. Barrier methods suffer from the same
ill-conditioning as penalty methods, but their modern form, the
**primal–dual interior-point method**, which applies Newton's method to the
perturbed KKT conditions, is among the most efficient algorithms for linear,
quadratic and convex programs (Nocedal & Wright, 2006, Ch. 19; Fiacco &
McCormick, 1968).

(sec-06-06-augmented-lagrangian)=
## Augmented Lagrangian methods

The quadratic penalty needs $\rho\to\infty$ because the multiplier
information is carried by the violation. The **augmented Lagrangian**
supplies the multipliers explicitly. For equality constraints it is

$$
L_A(\vx,\vmu;\rho)=f(\vx)+\sum_{j=1}^{\ell}\mu_jh_j(\vx)+\frac{\rho}{2}\sum_{j=1}^{\ell}h_j(\vx)^2 .
$$ (eq-06-06-augmented-lagrangian)

Each outer iteration minimizes $L_A$ over $\vx$ for fixed $\vmu$ and $\rho$, then
updates the multipliers by $\mu_j\leftarrow\mu_j+\rho\,h_j(\vx)$, the analog of
{eq}`eq-06-06-multiplier-estimates`. If $\vmu$ is the correct multiplier vector and the second-order sufficient
conditions ({prf:ref}`thm-06-03-second-order`) hold at $\vx^*$, then $\vx^*$ is a
local minimizer of $L_A$ for every sufficiently large but *finite* $\rho$
(Nocedal & Wright, 2006, Thm. 17.5), so
the ill-conditioning is avoided. For $\min x$ s.t. $h(x)=2-x=0$ with $\rho=1$
and $\mu_0=0$: minimizing $x+\mu(2-x)+\tfrac12(2-x)^2$ gives
$x=2-(1-\mu)/\rho=1$; the update gives $\mu_1=0+1\cdot(2-1)=1$, and the next
minimization gives $x=2$ exactly, with $\mu=1$ the correct multiplier.
Inequalities are handled by slack variables or by a modified update. The
method (also called the method of multipliers; Hestenes, 1969; Powell, 1969)
underlies many
large-scale NLP codes (Nocedal & Wright, 2006, Ch. 17; Bazaraa et al., 2006,
Sec. 9.3).

## Summary and Common Mistakes

### Summary

- Constrained problems are solved either by feasible-direction methods or by
  transforming them into unconstrained problems.
- Frank–Wolfe ({prf:ref}`alg-06-06-frank-wolfe`) solves the LP
  {eq}`eq-06-06-fw-lp` for a vertex $\vy^{(k)}$, moves along $\vy^{(k)}-\vx^{(k)}$
  with $\alpha\in[0,1]$, and stops when the gap $\Gamma$
  ({prf:ref}`def-06-06-fw-gap`) is small; the gap gives the lower bound
  {eq}`eq-06-06-fw-bound` ({prf:ref}`thm-06-06-fw-gap`).
- Its rate is $O(1/k)$ ({prf:ref}`thm-06-06-fw-rate`); on
  {prf:ref}`ex-06-06-frank-wolfe` the iterates zig-zag between two vertices
  towards $\vx^*=(1,1.5)$, $f^*=-11.5$ ({numref}`fig-06-06-frank-wolfe`). In
  network flow problems its LP is a shortest-path problem.
- The quadratic penalty $\psi$ ({prf:ref}`def-06-06-quadratic-penalty`) and the
  method {prf:ref}`alg-06-06-penalty` produce minimizers $\vx(\rho)$ that approach
  a solution from outside ({prf:ref}`thm-06-06-penalty-convergence`,
  {numref}`fig-06-06-penalty-family`); $2\rho\max\{0,g_i\}$ and $2\rho h_j$
  estimate the KKT multipliers {eq}`eq-06-06-multiplier-estimates`.
- Large $\rho$ makes $\Phi_\rho$ ill-conditioned; barrier methods
  {eq}`eq-06-06-log-barrier` approach from inside, and augmented Lagrangians
  {eq}`eq-06-06-augmented-lagrangian` reach the solution with a finite $\rho$.

### Common mistakes

- **Writing the segment the wrong way round.** The Frank–Wolfe update is
  $\vx^{(k)}+\alpha(\vy^{(k)}-\vx^{(k)})=(1-\alpha)\vx^{(k)}+\alpha\vy^{(k)}$; the weight
  $\alpha$ belongs to the new vertex.
- **Letting the step exceed 1.** In Frank–Wolfe, an unconstrained line
  minimizer $\alpha>1$ must be clipped to $1$; beyond $\vy^{(k)}$ the point is
  generally infeasible.
- **Using the gap of the wrong point.** The stopping test uses
  $\nabla f(\vx^{(k)})$ at the *current* point; the gradient of the previous
  iteration gives a meaningless test.
- **Expecting Frank–Wolfe to stop exactly.** Unless the optimum is a vertex,
  the gap only tends to zero; stop on $\Gamma(\vx^{(k)})\le\varepsilon$ and
  report the bound $f(\vx^{(k)})-\Gamma(\vx^{(k)})$.
- **Taking penalty solutions as feasible.** $\vx(\rho)$ violates active
  constraints by about $\lambda_i^*/(2\rho)$; a final feasibility check or
  projection is needed.
- **Starting with a huge $\rho$.** It makes the first subproblem
  ill-conditioned; increase $\rho$ gradually and warm-start.
- **Confusing slow progress with divergence.** A penalty sequence whose
  $\rho_k\psi(\vx^{(k)})$ is still above the tolerance after a few rounds has not
  failed; it has not yet met the stopping test
  ({prf:ref}`ex-08-09-penalty-parabola`).


## References

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 9: exterior
  penalty, exact penalty, augmented Lagrangian and barrier methods;
  Chapter 10: methods of feasible directions.)
- Courant, R. (1943). Variational methods for the solution of problems of
  equilibrium and vibrations. *Bulletin of the American Mathematical Society,
  49*(1), 1–23. <https://doi.org/10.1090/S0002-9904-1943-07818-4>
- Fiacco, A. V., & McCormick, G. P. (1968). *Nonlinear programming: Sequential
  unconstrained minimization techniques*. Wiley.
- Frank, M., & Wolfe, P. (1956). An algorithm for quadratic programming.
  *Naval Research Logistics Quarterly, 3*(1–2), 95–110.
  <https://doi.org/10.1002/nav.3800030109>
- Hestenes, M. R. (1969). Multiplier and gradient methods. *Journal of
  Optimization Theory and Applications, 4*(5), 303–320.
  <https://doi.org/10.1007/BF00927673>
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming* (4th
  ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3> (Chapter 12:
  primal methods, feasible directions; Chapter 13: penalty and barrier
  methods.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Chapter 17: quadratic
  penalty, exact penalty and augmented Lagrangian methods; Chapter 18:
  sequential quadratic programming; Chapter 19: interior-point methods.)
- Powell, M. J. D. (1969). A method for nonlinear constraints in minimization
  problems. In R. Fletcher (Ed.), *Optimization* (pp. 283–298). Academic
  Press.
- Rosen, J. B. (1960). The gradient projection method for nonlinear
  programming. Part I. Linear constraints. *Journal of the Society for
  Industrial and Applied Mathematics, 8*(1), 181–217.
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis with
  mathematical programming methods*. Prentice-Hall. (Chapters 4–5: the
  convex-combination (Frank–Wolfe) algorithm and its application to traffic
  assignment.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 21: constrained nonlinear algorithms.)
- Wilson, R. B. (1963). *A simplicial algorithm for concave programming*
  [Doctoral dissertation, Harvard University].
- Wolfe, P. (1963). Methods of nonlinear programming. In R. L. Graves &
  P. Wolfe (Eds.), *Recent advances in mathematical programming*
  (pp. 67–86). McGraw-Hill.
- Zoutendijk, G. (1960). *Methods of feasible directions: A study in linear
  and nonlinear programming*. Elsevier.
