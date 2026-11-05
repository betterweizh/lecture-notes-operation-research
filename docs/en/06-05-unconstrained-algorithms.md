# Algorithms for Unconstrained Optimization

## Overview

This chapter turns the one-dimensional searches of {doc}`06-04-line-search` into
methods for minimizing a smooth function of $n$ variables,
$\min_{\vx\in\R^n}f(\vx)$. All of them are **descent methods**: at the current
point they choose a direction along which $f$ decreases and then a step length
along it. The **steepest-descent method** uses the negative gradient; it is
simple and globally convergent, but its successive steps are orthogonal, its
path zig-zags across narrow valleys, and its linear rate deteriorates with the
condition number of the Hessian. **Newton's method** uses the Hessian to build
a quadratic model; it minimizes a convex quadratic in one step and converges
quadratically near a minimizer with a positive definite Hessian, but needs
second derivatives and safeguards when the Hessian is not positive definite.
Between the two lie the **quasi-Newton** (BFGS, after Broyden, Fletcher,
Goldfarb and Shanno) and **conjugate-gradient**
methods, which use only gradients and still converge superlinearly or
finitely on quadratics. All methods are run on the same ill-conditioned
quadratic, $f(\vx)=(x_1-2)^2+10(x_2-2)^2$. The chapter relies on the
optimality conditions of {doc}`06-02-optimality-conditions`, on eigenvalues
and definiteness from {doc}`01-01-linear-algebra`, and on gradients,
Hessians and Taylor expansions from {doc}`01-02-matrix-calculus`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- state the general descent algorithm and check that a direction
  $\vd=-\mD\nabla f(\vx)$ with $\mD\succ\vzero$ is a descent direction;
- run the steepest-descent method with exact line search on a quadratic,
  using the closed-form step length, and explain the zig-zag path;
- bound the rate of steepest descent by the condition number of the Hessian;
- derive Newton's method in $\R^n$, show that it minimizes a convex quadratic
  in one step, and state its local quadratic convergence;
- modify the Newton direction when the Hessian is not positive definite; and
- describe the ideas behind the BFGS and conjugate-gradient methods and
  compare the methods by cost per iteration, storage and rate.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $\vg^{(k)}$ | gradient at the iterate, $\vg^{(k)}=\nabla f(\vx^{(k)})$ |
| $\mH(\vx)$ | Hessian $\nabla^2f(\vx)$; $\mH^{(k)}=\mH(\vx^{(k)})$ |
| $\mD_k$ | positive definite matrix that scales the gradient, $\vd^{(k)}=-\mD_k\vg^{(k)}$ |
| $\mQ$, $\vb$ | Hessian and linear term of a quadratic $f(\vx)=\tfrac12\vx^\T\mQ\vx-\vb^\T\vx$, $\mQ\succ\vzero$ |
| $\kappa$ | condition number $\lambda_{\max}(\mQ)/\lambda_{\min}(\mQ)$ |
| $\vs^{(k)}$, $\vy^{(k)}$ | step and gradient change, $\vs^{(k)}=\vx^{(k+1)}-\vx^{(k)}$, $\vy^{(k)}=\vg^{(k+1)}-\vg^{(k)}$ (quasi-Newton) |
| $\theta_k$ | angle between $\vd^{(k)}$ and $-\nabla f(\vx^{(k)})$ |
| $\gamma_k$ | BFGS scalar $1/\vy^{(k)\T}\vs^{(k)}$ |
| $\beta_k$ | conjugate-gradient coefficient |
| $\Delta_k$ | trust-region radius |
| $\bar{\mH}_k$, $\nu_k$ | positive definite modification of the Hessian and the shift used to build it, $\bar{\mH}_k=\nabla^2f(\vx^{(k)})+\nu_k\mI$ |

(sec-06-05-descent-methods)=
## Descent methods

For the unconstrained problem $\min_{\vx\in\R^n}f(\vx)$ with $f$ continuously
differentiable, every method of this chapter is an instance of the following
scheme ({doc}`06-04-line-search`, {eq}`eq-06-04-update`).

:::{prf:algorithm} Descent method
:label: alg-06-05-descent

**Input:** a starting point $\vx^{(0)}$; a tolerance $\varepsilon>0$.
**Output:** an approximate stationary point.

1. **Initialize.** Set $k=0$.
2. **Stopping test.** If $\norm{\nabla f(\vx^{(k)})}\le\varepsilon$, stop.
3. **Direction.** Choose a descent direction $\vd^{(k)}$,
   $\nabla f(\vx^{(k)})^\T\vd^{(k)}<0$.
4. **Step length.** Choose $\alpha_k>0$ by a line search on
   $\varphi(\alpha)=f(\vx^{(k)}+\alpha\vd^{(k)})$ (exact, or satisfying the Wolfe
   conditions).
5. **Update.** Set $\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$,
   $k\leftarrow k+1$, and go to Step 2.
:::

A large family of directions is obtained by scaling the negative gradient
with a positive definite matrix.

:::{prf:proposition} Scaled gradient directions are descent directions
:label: prop-06-05-scaled-gradient

If $\nabla f(\vx)\ne\vzero$ and $\mD\succ\vzero$, then $\vd=-\mD\nabla f(\vx)$ is a
descent direction at $\vx$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

$\nabla f(\vx)^\T\vd=-\nabla f(\vx)^\T\mD\nabla f(\vx)<0$ because $\mD$ is
positive definite and $\nabla f(\vx)\ne\vzero$.
:::

The choice $\mD_k=\mI$ gives steepest descent, $\mD_k=[\nabla^2f(\vx^{(k)})]^{-1}$
(when positive definite) gives Newton's method, and matrices $\mD_k$ built from
gradient differences give the quasi-Newton methods.

**Global convergence.** Descent alone does not guarantee convergence: the
directions could become almost orthogonal to the gradient, or the steps too
short. Both are excluded by two mild requirements. Let $\theta_k$ be the angle
between $\vd^{(k)}$ and $-\nabla f(\vx^{(k)})$. If $f$ is bounded below, its
gradient is Lipschitz continuous, the steps satisfy the Wolfe conditions
({eq}`eq-06-04-armijo`, {eq}`eq-06-04-curvature`) and $\cos\theta_k\ge\delta>0$ for all $k$, then
$\nabla f(\vx^{(k)})\to\vzero$ (Zoutendijk's theorem; Zoutendijk, 1970; Nocedal & Wright, 2006,
Thm. 3.2). For $\vd^{(k)}=-\mD_k\vg^{(k)}$ the angle condition holds when the
condition numbers of the matrices $\mD_k$ are bounded. Consequently, if the
iterates stay bounded, every limit point is a stationary point. As for any
method based on first-order information, this is a guarantee of
stationarity, not of a minimum; in practice the descent property makes
convergence to a saddle point or maximizer unlikely.

(sec-06-05-steepest-descent)=
## The steepest-descent method

Among all directions of unit length, which one decreases $f$ fastest? The
rate of change of $f$ at $\vx$ along $\vd$ is the directional derivative
$\nabla f(\vx)^\T\vd$, and by the Cauchy–Schwarz inequality

$$
\begin{aligned}
\nabla f(\vx)^\T\vd&\ge-\norm{\nabla f(\vx)}\,\norm{\vd}\\
&=-\norm{\nabla f(\vx)}\qquad\text{for }\norm{\vd}=1,
\end{aligned}
$$

with equality for $\vd=-\nabla f(\vx)/\norm{\nabla f(\vx)}$. The negative
gradient is therefore the **steepest-descent direction**; geometrically it is
perpendicular to the level set of $f$ through $\vx$ ({prf:ref}`thm-01-02-gradient-level-set`).

:::{prf:algorithm} Steepest descent (Cauchy's method; Cauchy, 1847)
:label: alg-06-05-steepest-descent

**Input:** a starting point $\vx^{(0)}$; a tolerance $\varepsilon>0$.
**Output:** an approximate stationary point.

1. **Initialize.** Set $k=0$.
2. **Stopping test.** Compute $\vg^{(k)}=\nabla f(\vx^{(k)})$. If
   $\norm{\vg^{(k)}}\le\varepsilon$, stop.
3. **Direction.** Set $\vd^{(k)}=-\vg^{(k)}$.
4. **Line search.** Find $\alpha_k$ minimizing $f(\vx^{(k)}+\alpha\vd^{(k)})$
   over $\alpha\ge0$ (or satisfying the Wolfe conditions).
5. **Update.** Set $\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$,
   $k\leftarrow k+1$, and go to Step 2.
:::

With an exact line search the method is globally convergent in the sense of
{ref}`sec-06-05-descent-methods`, and here the argument is short enough to give
in full.

:::{prf:proposition} Limit points of steepest descent are stationary
:label: prop-06-05-sd-convergence

Let $f$ be continuously differentiable, and let steepest descent with exact
line search generate a sequence $\{\vx^{(k)}\}$ that has a limit point
$\bar\vx$. Then $\nabla f(\bar\vx)=\vzero$.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Suppose $\nabla f(\bar\vx)\ne\vzero$. Since $-\nabla f(\bar\vx)$ is a descent
direction, there is $\bar\alpha>0$ with
$f(\bar\vx-\bar\alpha\nabla f(\bar\vx))=f(\bar\vx)-2\delta$ for some
$\delta>0$. By continuity of $f$ and $\nabla f$, every $\vx$ close enough to
$\bar\vx$ satisfies $f(\vx-\bar\alpha\nabla f(\vx))\le f(\bar\vx)-\delta$. Let
$\vx^{(k)}\to\bar\vx$ along a subsequence $\mathcal K$. For $k\in\mathcal K$
large, the exact line search gives
$f(\vx^{(k+1)})\le f(\vx^{(k)}-\bar\alpha\vg^{(k)})\le f(\bar\vx)-\delta$. But
the values $f(\vx^{(k)})$ decrease monotonically and converge to $f(\bar\vx)$
along $\mathcal K$, so $f(\vx^{(j)})\ge f(\bar\vx)$ for all $j$, a
contradiction.
:::

In particular, if the level set $\{\vx:f(\vx)\le f(\vx^{(0)})\}$ is bounded,
the sequence has limit points, and all of them are stationary. The
proposition says nothing about speed; that is the subject of the rest of
this section.

**Exact step for a quadratic.** For the strictly convex quadratic

$$
\begin{aligned}
f(\vx)&=\tfrac12\vx^\T\mQ\vx-\vb^\T\vx+c, \qquad \mQ\succ\vzero,\\
\nabla f(\vx)&=\mQ\vx-\vb,
\end{aligned}
$$ (eq-06-05-quadratic)

the line function along $\vd=-\vg$ is the parabola
$\varphi(\alpha)=f(\vx)-\alpha\vg^\T\vg+\tfrac12\alpha^2\vg^\T\mQ\vg$. Setting
$\varphi'(\alpha)=-\vg^\T\vg+\alpha\vg^\T\mQ\vg=0$ gives the exact step

$$
\alpha_k=\frac{\vg^{(k)\T}\vg^{(k)}}{\vg^{(k)\T}\mQ\vg^{(k)}} .
$$ (eq-06-05-sd-step)

By the Rayleigh-quotient bounds ({prf:ref}`prop-01-01-rayleigh`),
$1/\lambda_{\max}(\mQ)\le\alpha_k\le1/\lambda_{\min}(\mQ)$.

**Orthogonal steps.** With an exact line search, $\varphi'(\alpha_k)=0$, that
is, $\nabla f(\vx^{(k+1)})^\T\vd^{(k)}=0$. For steepest descent
$\vd^{(k)}=-\vg^{(k)}$ and $\vd^{(k+1)}=-\vg^{(k+1)}$, so

$$
\vd^{(k+1)\T}\vd^{(k)}=\vg^{(k+1)\T}\vg^{(k)}=0 :
$$ (eq-06-05-orthogonal)

**successive steps are orthogonal**, whatever $f$ is. The new direction is
perpendicular to the old one, and the old direction is tangent to the level
set at the new point. In a long, narrow valley this forces the path to cross
the valley again and again in short orthogonal steps instead of following its
floor: the **zig-zag** of {numref}`fig-06-05-zigzag`.

:::{prf:example} Steepest descent on an ill-conditioned quadratic
:label: ex-06-05-steepest-descent

Minimize $f(\vx)=(x_1-2)^2+10(x_2-2)^2$ by steepest descent with exact line
search, starting from $\vx^{(0)}=(-4,-3)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The minimizer is $\vx^*=(2,2)$ with $f^*=0$. In the form
{eq}`eq-06-05-quadratic`, $\mQ=\diag(2,20)$, $\vb=(4,40)$, and
$\nabla f(\vx)=(2x_1-4,\ 20x_2-40)$. Write $u_k=4-2x_1^{(k)}$ and
$v_k=40-20x_2^{(k)}$ for the components of $\vd^{(k)}=-\vg^{(k)}$. The line
function is

$$
\varphi(\alpha)=\bigl(x_1^{(k)}+\alpha u_k-2\bigr)^2+10\bigl(x_2^{(k)}+\alpha v_k-2\bigr)^2,
$$

and $\varphi'(\alpha)=0$ gives (this is {eq}`eq-06-05-sd-step`)

$$
\alpha_k=\frac{u_k^2+v_k^2}{2u_k^2+20v_k^2}\in\Bigl[\frac1{20},\frac12\Bigr].
$$

At $\vx^{(0)}$: $\vg^{(0)}=(-12,-100)$, $\vd^{(0)}=(12,100)$,
$\alpha_0=(144+10\,000)/(288+200\,000)=0.05065$, and
$\vx^{(1)}=(-4,-3)+0.05065\,(12,100)=(-3.3922,\ 2.0647)$, with
$f(\vx^{(1)})=29.12$ against $f(\vx^{(0)})=286$. {numref}`tab-06-05-steepest-descent`
continues the iteration.

```{list-table} Steepest descent on $f(\vx)=(x_1-2)^2+10(x_2-2)^2$ from $(-4,-3)$
:name: tab-06-05-steepest-descent
:header-rows: 1

* - $k$
  - $\vx^{(k)}$
  - $\vd^{(k)}=-\nabla f(\vx^{(k)})$
  - $\alpha_k$
  - $f(\vx^{(k)})$
* - 0
  - $(-4.0000,\ -3.0000)$
  - $(12.0000,\ 100.0000)$
  - 0.05065
  - 286.0
* - 1
  - $(-3.3922,\ 2.0647)$
  - $(10.7845,\ -1.2941)$
  - 0.4434
  - 29.12
* - 2
  - $(1.3891,\ 1.4909)$
  - $(1.2217,\ 10.1811)$
  - 0.05065
  - 2.965
* - 3
  - $(1.4510,\ 2.0066)$
  - $(1.0980,\ -0.1318)$
  - 0.4434
  - 0.3018
* - 4
  - $(1.9378,\ 1.9482)$
  - $(0.1244,\ 1.0366)$
  - 0.05065
  - 0.03073
* - 5
  - $(1.9441,\ 2.0007)$
  - $(0.1118,\ -0.0134)$
  - 0.4434
  - $3.129\times10^{-3}$
* - 6
  - $(1.9937,\ 1.9947)$
  - $(0.0127,\ 0.1055)$
  - 0.05065
  - $3.185\times10^{-4}$
* - 7
  - $(1.9943,\ 2.0001)$
  - $(0.0114,\ -0.0014)$
  - 0.4434
  - $3.243\times10^{-5}$
* - 8
  - $(1.9994,\ 1.9995)$
  - $(0.0013,\ 0.0107)$
  - 0.05065
  - $3.302\times10^{-6}$
```

Three features stand out.

- **Zig-zag.** The directions alternate between "mostly $x_2$" (even $k$,
  short step $\alpha_k=0.05065$) and "mostly $x_1$" (odd $k$, long step
  $\alpha_k=0.4434$), and consecutive directions are orthogonal, for example
  $(12,100)^\T(10.7845,-1.2941)=129.4-129.4=0$ ({numref}`fig-06-05-zigzag`).
- **Linear convergence.** Every step multiplies $f(\vx^{(k)})-f^*$ by the same
  factor $0.1018$ ($29.12/286=2.965/29.12=\dots$): one decimal digit per
  iteration. The ratio is constant because in two variables the gradients
  alternate between two perpendicular directions, and the reduction factor
  $1-(\vg^\T\vg)^2/\bigl((\vg^\T\mQ\vg)(\vg^\T\mQ^{-1}\vg)\bigr)$ (derived
  in the proof of {prf:ref}`thm-06-05-sd-rate` below) takes the same value for a
  vector $\vg$ and for its perpendicular.
- **Cost.** The test $\norm{\nabla f}\le10^{-6}$ is met only at $k=17$,
  although the minimizer is at a distance of only $7.8$ from the start and
  each step is optimal along its line.
:::

```{figure} figures/06-05-zigzag.svg
:name: fig-06-05-zigzag
:width: 680px
:class: dark-light
:alt: Two contour plots of an elongated quadratic; the left shows the zig-zag steepest-descent path, the single Newton step and the two-step conjugate-gradient path to the minimizer; the right shows a slow zig-zag from the worst-case start

Contours of $f(\vx)=(x_1-2)^2+10(x_2-2)^2$, minimizer $\vx^*=(2,2)$. (a) From
$\vx^{(0)}=(-4,-3)$: steepest descent with exact line search (blue, orthogonal
zig-zag steps, {prf:ref}`ex-06-05-steepest-descent`), Newton's method (orange,
one step, {prf:ref}`ex-06-05-newton-quadratic`) and the conjugate-gradient / BFGS
path (green, two steps; its first step is the first steepest-descent step).
(b) Steepest descent from the worst-case start $(-8,1)$, drawn with the $x_2$
axis stretched: every step reduces $f$ only by the factor
$((\kappa-1)/(\kappa+1))^2=0.6694$ of {prf:ref}`thm-06-05-sd-rate`.
```

**How fast is steepest descent in general?** For quadratics the answer
depends only on the eigenvalues of $\mQ$.

:::{prf:theorem} Rate of steepest descent on a quadratic
:label: thm-06-05-sd-rate

Let $f$ be the quadratic {eq}`eq-06-05-quadratic` with $\mQ\succ\vzero$, minimizer
$\vx^*=\mQ^{-1}\vb$ and condition number
$\kappa=\lambda_{\max}(\mQ)/\lambda_{\min}(\mQ)$. Steepest descent with exact
line search satisfies

$$
f(\vx^{(k+1)})-f(\vx^*)\le\Bigl(\frac{\kappa-1}{\kappa+1}\Bigr)^2
\bigl(f(\vx^{(k)})-f(\vx^*)\bigr).
$$ (eq-06-05-sd-rate)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Let $E(\vx)=f(\vx)-f(\vx^*)=\tfrac12(\vx-\vx^*)^\T\mQ(\vx-\vx^*)$ and
$\vg=\mQ(\vx-\vx^*)$, so that $E(\vx)=\tfrac12\vg^\T\mQ^{-1}\vg$. With the step
{eq}`eq-06-05-sd-step`, $\varphi(\alpha_k)=f(\vx)-\tfrac12(\vg^\T\vg)^2/(\vg^\T\mQ\vg)$,
hence

$$
\frac{E(\vx^{(k+1)})}{E(\vx^{(k)})}
=1-\frac{(\vg^\T\vg)^2}{(\vg^\T\mQ\vg)(\vg^\T\mQ^{-1}\vg)} .
$$

The Kantorovich inequality (Luenberger & Ye, 2016, Ch. 8) states that for
$\mQ\succ\vzero$ and every $\vg\ne\vzero$

$$
\begin{aligned}
\frac{(\vg^\T\vg)^2}{(\vg^\T\mQ\vg)(\vg^\T\mQ^{-1}\vg)}
&\ge\frac{4\lambda_{\min}(\mQ)\,\lambda_{\max}(\mQ)}{\bigl(\lambda_{\min}(\mQ)+\lambda_{\max}(\mQ)\bigr)^2}\\
&=\frac{4\kappa}{(\kappa+1)^2} .
\end{aligned}
$$

Substituting,
$E(\vx^{(k+1)})/E(\vx^{(k)})\le\bigl((\kappa-1)/(\kappa+1)\bigr)^2$,
which is {eq}`eq-06-05-sd-rate`.
:::

The bound is attained. For $\mQ=\diag(2,20)$, $\kappa=10$ and the factor is
$(9/11)^2=0.6694$: from the start $(-8,1)$, whose error $(-10,-1)$ is
proportional to $(1/\lambda_1(\mQ),1/\lambda_2(\mQ))=(\tfrac12,\tfrac1{20})$,
the step length is always $\alpha_k=2/(\lambda_1(\mQ)+\lambda_2(\mQ))=1/11$ and
$f$ falls only by the factor $0.6694$ per step: $110$ → $73.64$ → $49.29$ → $33.00$ → $\dots$
({numref}`fig-06-05-zigzag`b). The start $(-4,-3)$ of
{prf:ref}`ex-06-05-steepest-descent` happens to be much more favorable
(factor $0.1018$). For $\kappa=1$ (circular contours) steepest descent reaches
$\vx^*$ in one step; for $\kappa=1000$ the worst-case factor is $0.996$, and
about $575$ iterations are needed per decimal digit of $f(\vx^{(k)})-f^*$.

**Fixed step lengths.** Many applications (machine learning among them) use
steepest descent with a constant step, $\vx^{(k+1)}=\vx^{(k)}-\alpha\nabla f(\vx^{(k)})$,
and no line search. On the quadratic {eq}`eq-06-05-quadratic` the error then
evolves linearly, $\vx^{(k+1)}-\vx^*=(\mI-\alpha\mQ)(\vx^{(k)}-\vx^*)$, and the
eigenvalues of $\mI-\alpha\mQ$ are $1-\alpha\lambda_i(\mQ)$. The iteration
converges if and only if $\abs{1-\alpha\lambda_i(\mQ)}<1$ for all $i$, that is,
$0<\alpha<2/\lambda_{\max}(\mQ)$. The best constant step,
$\alpha=2/(\lambda_{\min}(\mQ)+\lambda_{\max}(\mQ))$, balances the two extreme
eigenvalues and reduces the error norm by the factor
$(\kappa-1)/(\kappa+1)$ per step. For $\mQ=\diag(2,20)$ the steps must satisfy
$\alpha<0.1$ (with $\alpha=0.1$ the $x_2$ error just flips sign forever), and the
best one is $\alpha=\tfrac1{11}$, with factor $\tfrac9{11}=0.818$; this is the
step that the exact line search happens to choose at every iteration from the
worst-case start. A good constant step therefore requires knowledge of the
eigenvalues, which a line search supplies implicitly.

For a general smooth $f$, the same bound holds asymptotically near a
minimizer $\vx^*$ with $\nabla^2f(\vx^*)\succ\vzero$, with $\kappa$ the condition
number of $\nabla^2f(\vx^*)$ (Luenberger & Ye, 2016, Ch. 8). Steepest
descent is therefore reliable but slow on badly scaled problems.

(sec-06-05-scaling)=
### Scaling

The condition number depends on the units in which the variables are
measured, so a change of variables can help. Substitute $\vx=\mL\vy$ with an
invertible $\mL$ and apply steepest descent to $\tilde f(\vy)=f(\mL\vy)$. Since
$\nabla\tilde f(\vy)=\mL^\T\nabla f(\vx)$ ({prf:ref}`thm-01-02-chain-rule`), a step
$-\alpha\nabla\tilde f(\vy)$ in $\vy$ is the step
$-\alpha\mL\mL^\T\nabla f(\vx)$ in $\vx$: **scaled steepest descent** is the
descent method with $\vd=-\mD\nabla f$ and $\mD=\mL\mL^\T\succ\vzero$
({prf:ref}`prop-06-05-scaled-gradient`). The Hessian in the new variables is
$\mL^\T\nabla^2f\,\mL$, and the best choice of $\mL$ makes it close to a
multiple of $\mI$.

For the quadratic of {prf:ref}`ex-06-05-steepest-descent`, measuring the second
variable in units of $1/\sqrt{10}$, $y_1=x_1-2$ and $y_2=\sqrt{10}(x_2-2)$,
turns $f$ into $y_1^2+y_2^2$, whose contours are circles ($\kappa=1$). The
corresponding $\mD=\diag(1,\tfrac1{10})$ gives at $\vx^{(0)}$ the direction
$\vd=-\diag(1,\tfrac1{10})(-12,-100)=(12,10)$, and the exact step
$\alpha=\tfrac12$ along it reaches $(2,2)$ in one iteration, the same point as
the Newton step of {prf:ref}`ex-06-05-newton-quadratic`. In practice simple **diagonal scaling** (dividing each
variable by a typical magnitude) is cheap and often removes most of the
ill-conditioning caused by mixed units, for example tonnes against
kilograms or dollars against millions of dollars. Choosing $\mD$ equal to the
inverse Hessian is the ideal scaling, and leads to Newton's method.

(sec-06-05-newton)=
## Newton's method

Newton's method replaces $f$ near the current point by its second-order
Taylor model ({prf:ref}`thm-01-02-taylor`)

$$
\begin{aligned}
m_k(\vx)&=f(\vx^{(k)})+\nabla f(\vx^{(k)})^\T(\vx-\vx^{(k)})\\
&\quad+\tfrac12(\vx-\vx^{(k)})^\T\nabla^2f(\vx^{(k)})(\vx-\vx^{(k)}),
\end{aligned}
$$ (eq-06-05-newton-model)

and moves to the stationary point of the model. Setting
$\nabla m_k(\vx)=\nabla f(\vx^{(k)})+\nabla^2f(\vx^{(k)})(\vx-\vx^{(k)})=\vzero$
gives the **Newton equations** and the update

$$
\begin{aligned}
\nabla^2f(\vx^{(k)})\,\vd^{(k)}&=-\nabla f(\vx^{(k)}),\\
\vx^{(k+1)}&=\vx^{(k)}+\vd^{(k)} .
\end{aligned}
$$ (eq-06-05-newton-step)

In the framework of {prf:ref}`alg-06-05-descent` this is the direction
$\vd^{(k)}=-[\nabla^2f(\vx^{(k)})]^{-1}\nabla f(\vx^{(k)})$ with the step
$\alpha_k=1$; the generalization of the one-dimensional method
{eq}`eq-06-04-newton-1d`. When $\nabla^2f(\vx^{(k)})\succ\vzero$, the model is a
convex quadratic, $\vx^{(k+1)}$ is its minimizer, and $\vd^{(k)}$ is a descent
direction by {prf:ref}`prop-06-05-scaled-gradient`. In computation the inverse is
never formed: the linear system {eq}`eq-06-05-newton-step` is solved, by a
Cholesky factorization when the Hessian is positive definite (about $n^3/3$
operations).

:::{prf:algorithm} Newton's method
:label: alg-06-05-newton

**Input:** a twice continuously differentiable $f$; a starting point
$\vx^{(0)}$; a tolerance $\varepsilon>0$.
**Output:** an approximate stationary point.

1. **Initialize.** Set $k=0$.
2. **Stopping test.** If $\norm{\nabla f(\vx^{(k)})}\le\varepsilon$, stop.
3. **Newton direction.** Solve
   $\nabla^2f(\vx^{(k)})\vd^{(k)}=-\nabla f(\vx^{(k)})$ for $\vd^{(k)}$.
4. **Update.** Set $\vx^{(k+1)}=\vx^{(k)}+\vd^{(k)}$, $k\leftarrow k+1$, and go
   to Step 2.
:::

:::{prf:proposition} One Newton step minimizes a convex quadratic
:label: prop-06-05-newton-quadratic

For the quadratic {eq}`eq-06-05-quadratic` with $\mQ\succ\vzero$, Newton's method
reaches the minimizer $\vx^*=\mQ^{-1}\vb$ in one step from any starting point.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

The Taylor model {eq}`eq-06-05-newton-model` of a quadratic is the quadratic
itself, so its minimizer is $\vx^*$. Explicitly,
$\vx^{(1)}=\vx^{(0)}-\mQ^{-1}(\mQ\vx^{(0)}-\vb)=\mQ^{-1}\vb$.
:::

:::{prf:example} Newton's method on the same quadratic
:label: ex-06-05-newton-quadratic

Apply one Newton step to $f(\vx)=(x_1-2)^2+10(x_2-2)^2$ from
$\vx^{(0)}=(-4,-3)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The Hessian is constant, $\nabla^2f=\diag(2,20)$, and
$\nabla f(\vx^{(0)})=(-12,-100)$. The Newton equations
$\diag(2,20)\,\vd^{(0)}=(12,100)$ give $\vd^{(0)}=(6,5)$, and

$$
\begin{aligned}
\vx^{(1)}&=(-4,-3)+(6,5)=(2,2)=\vx^*,\\
f(\vx^{(1)})&=0 .
\end{aligned}
$$

The Newton direction $(6,5)$ points straight at the minimizer, whereas the
steepest-descent direction $(12,100)$ points almost along the $x_2$ axis,
where the contours are crowded ({numref}`fig-06-05-zigzag`a).
:::

As {prf:ref}`ex-06-05-newton-quadratic` shows, scaling the gradient by the
inverse Hessian, here $\mH^{-1}=\diag(\tfrac12,\tfrac1{20})$, removes
exactly the distortion that the condition number $\kappa=10$ causes.
Newton's method is unaffected by any linear change of variables (it is
**affine invariant**), so the condition number of the Hessian plays no role
in it.

For non-quadratic functions the model is only an approximation, and Newton's
method is fast only near a minimizer.

:::{prf:theorem} Local quadratic convergence of Newton's method
:label: thm-06-05-newton-quadratic

Let $f$ be twice continuously differentiable with a Lipschitz continuous
Hessian near $\vx^*$, and let $\nabla f(\vx^*)=\vzero$ and
$\nabla^2f(\vx^*)\succ\vzero$. Then there is a $\delta>0$ such that, for every
$\vx^{(0)}$ with $\norm{\vx^{(0)}-\vx^*}<\delta$, the Newton iterates
{eq}`eq-06-05-newton-step` are well defined, converge to $\vx^*$, and
$\norm{\vx^{(k+1)}-\vx^*}\le C\norm{\vx^{(k)}-\vx^*}^2$ for a constant $C$.
The gradient norms $\norm{\nabla f(\vx^{(k)})}$ converge to zero quadratically
as well.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

(Sketch; Nocedal & Wright, 2006, Thm. 3.5.) Write $\mH_k=\nabla^2f(\vx^{(k)})$,
which is invertible with $\norm{\mH_k^{-1}}\le2\norm{\nabla^2f(\vx^*)^{-1}}$
near $\vx^*$. Since $\nabla f(\vx^*)=\vzero$,

$$
\vx^{(k+1)}-\vx^*=\mH_k^{-1}\Bigl[\mH_k(\vx^{(k)}-\vx^*)-\bigl(\nabla f(\vx^{(k)})-\nabla f(\vx^*)\bigr)\Bigr].
$$

By Taylor's theorem in integral form, the bracket equals
$\int_0^1\bigl(\mH_k-\nabla^2f(\vx^*+t(\vx^{(k)}-\vx^*))\bigr)(\vx^{(k)}-\vx^*)\,dt$,
whose norm is at most $\tfrac12L\norm{\vx^{(k)}-\vx^*}^2$ when $L$ is the
Lipschitz constant of the Hessian. Hence
$\norm{\vx^{(k+1)}-\vx^*}\le L\norm{\nabla^2f(\vx^*)^{-1}}\,\norm{\vx^{(k)}-\vx^*}^2$,
and the argument of {prf:ref}`thm-06-04-newton-quadratic` completes the proof.
:::

Both assumptions matter. When $\nabla^2f(\vx^*)$ is singular, the convergence
can drop to linear: on $f(\vx)=(x_1-2)^4+(x_1-2x_2)^2$, whose Hessian at the
minimizer $(2,1)$ is singular, every Newton step reduces the error only by
the factor $\tfrac23$ ({prf:ref}`ex-08-09-newton-singular-hessian`). And far from
$\vx^*$ the Newton step can fail in all the ways of the one-dimensional method
({doc}`06-04-line-search`): the Hessian may be indefinite, the step may be huge,
or $f$ may increase.

(sec-06-05-newton-modifications)=
### Making Newton's method robust

Two modifications turn Newton's method into a globally convergent algorithm
that keeps its fast local rate.

**Line search (damped Newton).** Use $\vd^{(k)}$ as a direction and choose
$\alpha_k$ by backtracking from $\alpha=1$ ({prf:ref}`alg-06-04-backtracking`).
Near the solution the unit step satisfies the Armijo condition (with
$c_1<\tfrac12$) and is accepted, so the quadratic rate is retained; far from
it, shorter steps prevent increases of $f$.

**Hessian modification.** If $\nabla^2f(\vx^{(k)})$ is not positive definite,
the Newton direction may point uphill or towards a saddle point. It is then
replaced by $\vd^{(k)}=-\bar{\mH}_k^{-1}\nabla f(\vx^{(k)})$, where $\bar{\mH}_k\succ\vzero$ is
close to the Hessian:

- **shift:** $\bar{\mH}_k=\nabla^2f(\vx^{(k)})+\nu_k\mI$ with $\nu_k\ge0$ just large
  enough that $\lambda_{\min}(\bar{\mH}_k)\ge\delta>0$ (the Levenberg–Marquardt
  idea; Levenberg, 1944; Marquardt, 1963; large $\nu_k$ turns the step into a short steepest-descent step);
- **eigenvalue modification:** with the spectral decomposition
  $\nabla^2f(\vx^{(k)})=\mV\diag\bigl(\lambda_i(\nabla^2f(\vx^{(k)}))\bigr)\mV^\T$
  ({prf:ref}`thm-01-01-spectral`), use
  $\bar{\mH}_k=\mV\diag\bigl(\max\{\abs{\lambda_i(\nabla^2f(\vx^{(k)}))},\delta\}\bigr)\mV^\T$, which flips negative
  curvature into positive curvature;
- **modified Cholesky:** run a Cholesky factorization and enlarge diagonal
  entries whenever a pivot would be too small, which yields $\bar{\mH}_k\succ\vzero$ at
  about the cost of one factorization (the method of choice in practice).

As an illustration take $f(\vx)=\tfrac14x_1^4-\tfrac12x_1^2+x_2^2$, which has a
saddle point at $\vzero$ and minimizers $(\pm1,0)$, at $\vx=(0.2,1)$. There
$\nabla f=(x_1^3-x_1,\ 2x_2)=(-0.192,\ 2)$ and
$\nabla^2f=\diag(3x_1^2-1,\ 2)=\diag(-0.88,\ 2)$ is indefinite. The pure Newton
step $\vd=(-0.2182,\ -1)$ is still a descent direction here
($\nabla f^\T\vd=-1.958$), but in the $x_1$ coordinate it moves *towards* the
saddle point: the full step lands at $(-0.0182,0)$, next to $\vzero$. The
eigenvalue modification replaces $-0.88$ by $0.88$ and gives
$\vd=(0.2182,\ -1)$, which moves away from the saddle; with backtracking
the modified iteration then reaches the minimizer $(1,0)$, $f=-0.25$, to six
digits in five iterations, all with $\alpha_k=1$.

**Trust regions.** A different way to globalize Newton's method reverses the
order of the two decisions: instead of choosing a direction and then a step
length, it first fixes how far the quadratic model {eq}`eq-06-05-newton-model`
can be trusted, a radius $\Delta_k$, and then minimizes the model over the
ball $\norm{\vx-\vx^{(k)}}\le\Delta_k$. The ratio of the actual decrease of $f$
to the decrease predicted by the model decides whether the step is accepted
and whether $\Delta_k$ is enlarged or reduced. The subproblem has a solution
even when the Hessian is indefinite (the ball is compact), and its solution
has the form $-(\nabla^2f+\nu\mI)^{-1}\nabla f$ with $\nu\ge0$ (apart from a
degenerate "hard case"), so the shift strategy above appears automatically,
with $\nu$, the multiplier of the ball constraint, chosen by the radius.
Trust-region methods are as widely used as line-search methods (Nocedal &
Wright, 2006, Ch. 4).

(sec-06-05-quasi-newton)=
## Quasi-Newton methods

Newton's method needs second derivatives and a linear solve per iteration;
steepest descent needs neither but is slow. **Quasi-Newton methods** use
directions $\vd^{(k)}=-\mD_k\vg^{(k)}$ where $\mD_k\succ\vzero$ is an approximation
of the *inverse* Hessian that is updated from one iteration to the next using
only gradients. The update is designed so that $\mD_{k+1}$ reproduces the
curvature observed along the last step. Since
$\nabla f(\vx^{(k+1)})-\nabla f(\vx^{(k)})\approx\nabla^2f\,(\vx^{(k+1)}-\vx^{(k)})$,
it imposes the **secant condition**

$$
\begin{aligned}
\mD_{k+1}\vy^{(k)}&=\vs^{(k)},\\
\vs^{(k)}&=\vx^{(k+1)}-\vx^{(k)},\\
\vy^{(k)}&=\vg^{(k+1)}-\vg^{(k)} .
\end{aligned}
$$ (eq-06-05-secant)

Among the many matrices that satisfy {eq}`eq-06-05-secant`, the **BFGS** formula
(Broyden, 1970; Fletcher, 1970; Goldfarb, 1970; Shanno, 1970) chooses the one closest to
$\mD_k$ in a weighted norm:

$$
\begin{aligned}
\mD_{k+1}&=\Bigl(\mI-\gamma_k\vs^{(k)}\vy^{(k)\T}\Bigr)\mD_k\Bigl(\mI-\gamma_k\vy^{(k)}\vs^{(k)\T}\Bigr)\\
&\quad+\gamma_k\vs^{(k)}\vs^{(k)\T},\\
\gamma_k&=\frac1{\vy^{(k)\T}\vs^{(k)}} .
\end{aligned}
$$ (eq-06-05-bfgs)

If $\mD_k\succ\vzero$ and the **curvature condition** $\vy^{(k)\T}\vs^{(k)}>0$
holds, then $\mD_{k+1}\succ\vzero$; a line search satisfying the Wolfe conditions
{eq}`eq-06-04-curvature` guarantees $\vy^{(k)\T}\vs^{(k)}>0$, which is why
quasi-Newton codes use Wolfe line searches. The update costs $O(n^2)$
operations, no Hessian and no linear solve, and BFGS converges
**superlinearly** near a minimizer with a positive definite Hessian
(Nocedal & Wright, 2006, Sec. 6.1 and 6.4). The older DFP formula (Davidon–Fletcher–Powell; Davidon, 1959; Fletcher & Powell, 1963) updates the
same way with the roles of $\vs$ and $\vy$ exchanged in a Hessian
approximation; limited-memory BFGS (L-BFGS) stores only the last few pairs
$(\vs^{(k)},\vy^{(k)})$ and handles problems with millions of variables.

On the quadratic of {prf:ref}`ex-06-05-steepest-descent`, BFGS with $\mD_0=\mI$
and exact line searches takes the steepest-descent step first, to
$\vx^{(1)}=(-3.3922,\ 2.0647)$, and updates the matrix to

$$
\mD_1=\begin{pmatrix}1.0122&-0.006146\\-0.006146&0.05007\end{pmatrix}.
$$

The second direction $-\mD_1\vg^{(1)}=(10.9240,\ -0.1311)$ with
$\alpha_1=0.4936$ lands exactly on $\vx^*=(2,2)$, and the next update gives
$\mD_2=\diag(\tfrac12,\tfrac1{20})=\mQ^{-1}$. In general, on an $n$-variable
convex quadratic with exact line searches, BFGS terminates in at most $n$
steps and ends with $\mD_n=\mQ^{-1}$.

(sec-06-05-conjugate-gradient)=
## Conjugate-gradient methods

For very large problems even an $n\times n$ matrix is too much to store. The
**conjugate-gradient** (CG) method needs only a few vectors. Two directions
are **$\mQ$-conjugate** if $\vd^{(i)\T}\mQ\vd^{(j)}=0$ for $i\ne j$. For the
quadratic {eq}`eq-06-05-quadratic`, exact line searches along $n$ mutually
$\mQ$-conjugate directions reach the minimizer in at most $n$ steps, because
in the coordinates defined by the conjugate directions the quadratic separates
into $n$ one-dimensional problems. CG generates such directions on the fly
from the gradients:

$$
\begin{aligned}
\vd^{(0)}&=-\vg^{(0)},\\
\vd^{(k+1)}&=-\vg^{(k+1)}+\beta_k\vd^{(k)},\\
\beta_k&=\frac{\vg^{(k+1)\T}\vg^{(k+1)}}{\vg^{(k)\T}\vg^{(k)}}
\end{aligned}
$$ (eq-06-05-cg)

(the Fletcher–Reeves formula; Hestenes & Stiefel, 1952, introduced CG for
linear systems and Fletcher & Reeves, 1964, extended it to non-quadratic
functions), with $\alpha_k$ from an exact line search (for
the quadratic, $\alpha_k=-\vg^{(k)\T}\vd^{(k)}/\vd^{(k)\T}\mQ\vd^{(k)}$). The
new direction is the steepest-descent direction corrected by a multiple of
the old one; the correction is exactly what removes the zig-zag.

On the quadratic of {prf:ref}`ex-06-05-steepest-descent` the first CG step is the
steepest-descent step to $(-3.3922,\ 2.0647)$. Then
$\beta_0=\norm{\vg^{(1)}}^2/\norm{\vg^{(0)}}^2=0.01163$, the new direction
$\vd^{(1)}=(10.7845,\ -1.2941)+0.01163\,(12,100)=(10.9240,\ -0.1311)$ is
$\mQ$-conjugate to $\vd^{(0)}$, and the exact step $\alpha_1=0.4936$ reaches
$\vx^*=(2,2)$: two steps for two variables, against seventeen for steepest
descent ({numref}`fig-06-05-zigzag`a). On quadratics with exact line searches,
CG and BFGS started from $\mD_0=\mI$ generate the same iterates, as here. For
non-quadratic functions CG is restarted every $n$ iterations (or uses the
Polak–Ribière formula (Polak & Ribière, 1969) $\beta_k=\vg^{(k+1)\T}(\vg^{(k+1)}-\vg^{(k)})/\vg^{(k)\T}\vg^{(k)}$)
and converges linearly or superlinearly depending on the problem.

(sec-06-05-practice)=
## Derivatives and stopping in practice

All methods of this chapter need gradients, and Newton's method needs
Hessians. When they cannot be derived by hand, three sources are common.

- **Finite differences.** The forward difference
  $\partial f/\partial x_j\approx\bigl(f(\vx+h\ve_j)-f(\vx)\bigr)/h$ costs $n$
  extra function values per gradient. Its truncation error is of order $h$ and
  its rounding error of order $\epsilon_{\text{mach}}/h$, so
  $h\approx\sqrt{\epsilon_{\text{mach}}}\,\max\{1,\abs{x_j}\}\approx10^{-8}\max\{1,\abs{x_j}\}$
  balances the two and gives about half the digits of the function values.
  Central differences, $\bigl(f(\vx+h\ve_j)-f(\vx-h\ve_j)\bigr)/(2h)$, are more
  accurate but cost $2n$ values.
- **Automatic differentiation.** Programs that apply the chain rule to the
  code that evaluates $f$ deliver gradients exact to rounding at a small
  multiple of the cost of $f$ (the reverse mode), independent of $n$. Modeling
  languages and modern numerical libraries do this routinely ({doc}`09-06-code-nlp`).
- **Quasi-Newton updates** replace the Hessian by gradient differences
  ({ref}`sec-06-05-quasi-newton`), which is why Newton's method with exact
  Hessians is used mainly when second derivatives are cheap.

**Stopping.** The test $\norm{\nabla f(\vx^{(k)})}\le\varepsilon$ of
{prf:ref}`alg-06-05-descent` depends on the units of $f$ and $\vx$. Practical codes
use a relative version such as
$\norm{\nabla f(\vx^{(k)})}\le\varepsilon\max\{1,\abs{f(\vx^{(k)})}\}$, combined
with a test on the step length and a limit on the number of iterations, and
report which test ended the run ({doc}`06-04-line-search`,
{ref}`sec-06-04-descent-step`). A run that ends on the iteration limit, as
steepest descent easily does on an ill-conditioned problem, has not found a
solution.

(sec-06-05-comparison)=
## Comparison of the methods

```{list-table} Unconstrained descent methods for a smooth $f$ of $n$ variables
:name: tab-06-05-comparison
:header-rows: 1

* - Method
  - Direction
  - Work per iteration
  - Storage
  - Local rate
* - Steepest descent
  - $-\vg^{(k)}$
  - one gradient
  - $O(n)$
  - linear, factor up to $((\kappa-1)/(\kappa+1))^2$
* - Conjugate gradient
  - $-\vg^{(k)}+\beta_{k-1}\vd^{(k-1)}$
  - one gradient
  - $O(n)$
  - finite ($\le n$ steps) on quadratics; linear to superlinear otherwise
* - Quasi-Newton (BFGS)
  - $-\mD_k\vg^{(k)}$
  - one gradient, $O(n^2)$ update
  - $O(n^2)$
  - superlinear
* - Newton
  - $-[\nabla^2f]^{-1}\vg^{(k)}$
  - gradient, Hessian, $O(n^3)$ solve
  - $O(n^2)$
  - quadratic
```

The choice depends on what is available and what is expensive. If second
derivatives are cheap and $n$ is moderate, Newton's method with a line search
and Hessian modification is hard to beat. If only gradients are available,
BFGS is the standard choice, and L-BFGS or CG for very large $n$. Steepest
descent remains important as a building block and as the benchmark against
which the other methods are explained.

## Summary and Common Mistakes

### Summary

- Descent methods ({prf:ref}`alg-06-05-descent`) move along directions
  $\vd=-\mD\nabla f$ with $\mD\succ\vzero$, which are descent directions
  ({prf:ref}`prop-06-05-scaled-gradient`); with Wolfe steps and directions not too
  close to orthogonal to the gradient they drive $\nabla f$ to zero.
- Steepest descent ({prf:ref}`alg-06-05-steepest-descent`) uses $\vd=-\nabla f$;
  on a quadratic the exact step is {eq}`eq-06-05-sd-step`; successive steps are
  orthogonal {eq}`eq-06-05-orthogonal` and zig-zag ({numref}`fig-06-05-zigzag`).
- With exact line search every limit point is stationary
  ({prf:ref}`prop-06-05-sd-convergence`); scaling the variables
  ({ref}`sec-06-05-scaling`) changes $\kappa$ and hence the speed.
- Its convergence is linear with factor up to $((\kappa-1)/(\kappa+1))^2$
  ({prf:ref}`thm-06-05-sd-rate`); on $(x_1-2)^2+10(x_2-2)^2$ from $(-4,-3)$ the
  factor is $0.1018$ per step ({prf:ref}`ex-06-05-steepest-descent`).
- Newton's method {eq}`eq-06-05-newton-step` minimizes a convex quadratic in one
  step ({prf:ref}`prop-06-05-newton-quadratic`, {prf:ref}`ex-06-05-newton-quadratic`)
  and converges quadratically near a minimizer with a positive definite
  Hessian ({prf:ref}`thm-06-05-newton-quadratic`); line searches and Hessian
  modifications make it robust.
- Quasi-Newton (BFGS, {eq}`eq-06-05-bfgs`) and conjugate-gradient
  ({eq}`eq-06-05-cg`) methods use only gradients; both minimize an $n$-variable
  convex quadratic in at most $n$ exact-line-search steps
  ({numref}`tab-06-05-comparison`).

### Common mistakes

- **Expecting steepest descent to head for the minimizer.** The negative
  gradient is perpendicular to the contour, not directed at $\vx^*$; on
  elongated contours it is nearly orthogonal to the direction of $\vx^*$.
- **Using the unit step in steepest descent.** $\vx-\nabla f(\vx)$ has no
  reason to be a good point; the step must come from a line search (here
  $\alpha_k$ alternates between $0.05065$ and $0.4434$).
- **Inverting the Hessian.** Solve $\nabla^2f\,\vd=-\nabla f$ (Cholesky); do not
  form $[\nabla^2f]^{-1}$.
- **Running pure Newton from far away.** Without a line search and a positive
  definite (modified) Hessian, Newton's method can increase $f$, diverge, or
  converge to a saddle point or maximizer.
- **Writing the stopping test without a norm.** $\abs{\nabla f(\vx)}\le\varepsilon$
  is meaningless for a vector; use $\norm{\nabla f(\vx)}\le\varepsilon$.


## References

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Sections 8.6–8.8:
  steepest descent, Newton's method and its modifications, conjugate
  directions and quasi-Newton methods.)
- Broyden, C. G. (1970). The convergence of a class of double-rank
  minimization algorithms 1. General considerations. *Journal of the
  Institute of Mathematics and Its Applications, 6*(1), 76–90.
  <https://doi.org/10.1093/imamat/6.1.76>
- Cauchy, A. (1847). Méthode générale pour la résolution des systèmes
  d'équations simultanées. *Comptes Rendus de l'Académie des Sciences, 25*,
  536–538.
- Davidon, W. C. (1959). *Variable metric method for minimization* (Report
  ANL-5990). Argonne National Laboratory.
- Fletcher, R. (1970). A new approach to variable metric algorithms. *The
  Computer Journal, 13*(3), 317–322. <https://doi.org/10.1093/comjnl/13.3.317>
- Fletcher, R., & Powell, M. J. D. (1963). A rapidly convergent descent method
  for minimization. *The Computer Journal, 6*(2), 163–168.
  <https://doi.org/10.1093/comjnl/6.2.163>
- Fletcher, R., & Reeves, C. M. (1964). Function minimization by conjugate
  gradients. *The Computer Journal, 7*(2), 149–154.
  <https://doi.org/10.1093/comjnl/7.2.149>
- Goldfarb, D. (1970). A family of variable-metric methods derived by
  variational means. *Mathematics of Computation, 24*(109), 23–26.
- Hestenes, M. R., & Stiefel, E. (1952). Methods of conjugate gradients for
  solving linear systems. *Journal of Research of the National Bureau of
  Standards, 49*(6), 409–436. <https://doi.org/10.6028/jres.049.044>
- Levenberg, K. (1944). A method for the solution of certain nonlinear
  problems in least squares. *Quarterly of Applied Mathematics, 2*(2),
  164–168. <https://doi.org/10.1090/qam/10666>
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming* (4th
  ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3> (Chapter 8:
  steepest descent, its rate and the Kantorovich inequality, Newton's method;
  Chapter 9: conjugate-direction methods; Chapter 10: quasi-Newton methods.)
- Marquardt, D. W. (1963). An algorithm for least-squares estimation of
  nonlinear parameters. *Journal of the Society for Industrial and Applied
  Mathematics, 11*(2), 431–441. <https://doi.org/10.1137/0111030>
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Chapters 2–3:
  line-search methods, Zoutendijk's theorem, rates of steepest descent and
  Newton's method, Hessian modifications; Chapter 4: trust-region methods;
  Chapter 5: conjugate gradients; Chapter 6: quasi-Newton methods and BFGS;
  Chapter 8: finite differences and automatic differentiation.)
- Polak, E., & Ribière, G. (1969). Note sur la convergence de méthodes de
  directions conjuguées. *Revue Française d'Informatique et de Recherche
  Opérationnelle, Série Rouge, 3*(16), 35–43.
- Shanno, D. F. (1970). Conditioning of quasi-Newton methods for function
  minimization. *Mathematics of Computation, 24*(111), 647–656.
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 21: nonlinear programming algorithms, gradient method.)
- Zoutendijk, G. (1970). Nonlinear programming, computational methods. In
  J. Abadie (Ed.), *Integer and nonlinear programming* (pp. 37–86).
  North-Holland.
