# Piecewise-Linear Approximation

## Overview

Travel time on a congested road grows steeply with traffic; the cost of a
treatment plant grows more slowly than its capacity. Such nonlinear
functions can be brought into linear programs (LPs) and mixed-integer linear
programs (MILPs) by replacing
them with **piecewise-linear** functions that agree with them at chosen
**breakpoints**. This chapter shows how. A piecewise-linear function is
written as a convex combination of its breakpoints with weights $\alpha_i$
(the **convex-combination** or $\lambda$-formulation), as a sum of filled
segments (the **incremental** or $\delta$-formulation), or by choosing one
segment (the **multiple-choice** formulation). The key distinction is
convexity: when the model minimizes a convex function (or maximizes a
concave one), a plain LP is enough; otherwise the weights must satisfy an
adjacency condition, a **special ordered set of type 2** (SOS2), which needs
binary variables or solver support. We bound the approximation error for
smooth functions, extend the idea to separable programming and to functions
of two variables, and compare the formulations by size and by the strength
of their LP relaxations. Two examples, a congested two-route network and
two treatment plants with economies of scale, are solved and compared with
their exact optima. The chapter builds on the binary variables and big-$M$
constraints of {doc}`07-01-logical-constraints` and the maxima of affine
functions of {doc}`07-02-linearization-techniques`.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- write the piecewise-linear interpolant of a function from its breakpoints
  and decide whether it is convex;
- explain why minimizing a convex piecewise-linear function needs no binary
  variables, and why a nonconvex one does;
- formulate a piecewise-linear function with SOS2 weights, with binary
  adjacency constraints, with the incremental formulation and with the
  multiple-choice formulation, and compare them;
- choose breakpoints for a twice-differentiable function from the
  interpolation error bound;
- approximate separable programs and functions of two variables with
  mixed-integer linear models; and
- compare an approximate solution with the exact optimum.
:::

## Notation

| Symbol | Meaning |
|---|---|
| $f$, $f^a$ | function to be approximated and its piecewise-linear approximation |
| $a_1<a_2<\dots<a_n$ | breakpoints; $[a_k,a_{k+1}]$ is segment $k$, $k=1,\dots,n-1$ |
| $h_k=a_{k+1}-a_k$, $s_k$ | length and slope of segment $k$, $s_k=\bigl(f(a_{k+1})-f(a_k)\bigr)/h_k$ |
| $\alpha_i\ge0$ | weight of breakpoint $i$ (convex-combination formulation), $\sum_i\alpha_i=1$ |
| $\delta_k\in\{0,1\}$ | binary variable attached to segment $k$ |
| $\xi_k\in[0,1]$ | filled fraction of segment $k$ (incremental formulation) |
| $\bar x_k$ | part of $x$ in segment $k$ (segment form, multiple-choice formulation) |
| $\kappa_k$ | bound on $\abs{f''}$ on segment $k$ |
| $\beta_j\in\{0,1\}$ | binary selecting the $y$-interval $[b_j,b_{j+1}]$ (functions of two variables) |
| $a_i$, $b_j$ | grid coordinates in $x$ and $y$ for a function $f(x,y)$ |
| $\varepsilon$ | small positive tolerance |

Many textbooks write the weights as $\lambda_i$, hence the name
"$\lambda$-formulation"; this book keeps $\lambda$ for Lagrange multipliers.

(sec-07-03-pwl)=
## Piecewise-linear functions

Let $f$ be defined on $[a_1,a_n]$ and choose breakpoints
$a_1<a_2<\dots<a_n$. The **piecewise-linear interpolant** $f^a$ joins the
points $\bigl(a_i,f(a_i)\bigr)$ by straight segments. A point
$\bar x\in[a_k,a_{k+1}]$ is a convex combination of the two ends of its
segment, $\bar x=\alpha_ka_k+\alpha_{k+1}a_{k+1}$ with
$\alpha_k+\alpha_{k+1}=1$, $\alpha_k,\alpha_{k+1}\ge0$, and

$$
\begin{aligned}
f^a(\bar x)&=\alpha_kf(a_k)+\alpha_{k+1}f(a_{k+1})\\
&=f(a_k)+s_k(\bar x-a_k),\\
\alpha_{k+1}&=\frac{\bar x-a_k}{h_k} .
\end{aligned}
$$ (eq-07-03-interp)

$f^a$ is continuous and linear on each segment. It is **convex** exactly
when its slopes increase, $s_1\le s_2\le\dots\le s_{n-1}$, and then it is the
maximum of its pieces,

$$
f^a(x)=\max_{k=1,\dots,n-1}\bigl(f(a_k)+s_k(x-a_k)\bigr),
$$

which is the form handled in {ref}`sec-07-02-convex-side`. If $f$ is convex, its
interpolant is convex and lies **above** $f$ (chords of a convex function lie
above its graph); if $f$ is concave, $f^a$ is concave and lies **below** $f$.
A function with both curvatures, or a piecewise-linear cost given directly by
a tariff with a discount, is in general neither.

{numref}`fig-07-03-interpolation` shows both cases with the functions of the
examples of {ref}`sec-07-03-examples`.

```{figure} figures/07-03-interpolation.svg
:name: fig-07-03-interpolation
:width: 660px
:class: dark-light
:alt: A convex travel-time function and a concave cost function with their piecewise-linear interpolants

Left: total travel time $T_1(v)$ on route 1 of {prf:ref}`ex-07-03-congestion`
(convex) and its interpolant with breakpoints every $1000$ veh/h, which lies
above it. Right: the concave plant cost $C_1(q)=q^{0.6}$ of
{prf:ref}`ex-07-03-economies` with breakpoints every $20$ ML/d; the interpolant
lies below the curve, and the dashed chord from $0$ to $80$ is the convex
envelope that an LP without adjacency conditions would use.
```

(sec-07-03-convex-case)=
## The convex case needs no binaries

Write $x$ and $f^a(x)$ with one weight per breakpoint:

$$
\begin{aligned}
x&=\sum_{i=1}^n\alpha_ia_i,\\
f^a&=\sum_{i=1}^n\alpha_if(a_i),\\
\sum_{i=1}^n\alpha_i&=1,\qquad \alpha_i\ge0 .
\end{aligned}
$$ (eq-07-03-cc)

These conditions describe every convex combination of the breakpoints, not
only combinations of two neighbors; the pair $(x,f^a)$ can be any point of
the convex hull of the points $\bigl(a_i,f(a_i)\bigr)$. The interpolant is
recovered only when at most two weights are positive and they are
**adjacent**.

:::{prf:proposition} Adjacency is automatic in the convex case
:label: prop-07-03-convex-case

Let $f^a$ be convex and let a linear model *minimize* an objective in which
$f^a$ appears with a positive coefficient, all other constraints involving
$x$ only (not the individual weights). Then the LP with {eq}`eq-07-03-cc` has an
optimal solution with at most two positive weights, and they are adjacent;
its optimal value equals that of the model with the exact interpolant.
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Fix $x$. Among all weights with $\sum_i\alpha_ia_i=x$, $\sum_i\alpha_i=1$,
$\valpha\ge\vzero$, the value $\sum_i\alpha_if(a_i)$ is at least $f^a(x)$,
because $f^a$ is convex and agrees with $f$ at the breakpoints:
$\sum_i\alpha_if(a_i)=\sum_i\alpha_if^a(a_i)\ge f^a\bigl(\sum_i\alpha_ia_i\bigr)$
by Jensen's inequality. The two adjacent weights of {eq}`eq-07-03-interp` attain
this lower bound. Replacing any optimal weights by these does not change $x$
and does not increase the objective, so the model is optimal with adjacent
weights and its value is the value with $f^a$.
:::

The same holds for maximizing a concave $f^a$. The **segment form** gives an
equivalent LP: write $x=a_1+\sum_k\bar x_k$ with $0\le\bar x_k\le h_k$ and
$f^a=f(a_1)+\sum_ks_k\bar x_k$. With increasing slopes the minimization fills
the cheaper segments first, so a segment is never used before the previous
one is full, and no ordering constraint is needed. This is exactly the
argument behind the piecewise-linear profit of {prf:ref}`ex-03-01-piecewise-profit`
and {prf:ref}`prop-03-01-concave-split`: a contribution of \$500 per hundred
cases up to three hundred cases and \$300 beyond is concave in a
maximization, so the LP sells at \$500 first.

For a **nonconvex** case the LP misbehaves. Minimizing a concave cost, it
prefers to represent $x$ by the two extreme breakpoints, because the chord
between them lies below the curve; the LP then optimizes the **convex
envelope** of $f^a$ (the dashed line in {numref}`fig-07-03-interpolation`)
instead of $f^a$. The adjacency condition must then be imposed explicitly.

(sec-07-03-sos)=
## Special ordered sets and the convex-combination formulation

:::{prf:definition} Special ordered sets
:label: def-07-03-sos

An ordered list of nonnegative variables $(\alpha_1,\dots,\alpha_n)$ is a
**special ordered set of type 1** (SOS1) if at most one of them is nonzero,
and a **special ordered set of type 2** (SOS2) if at most two are nonzero and,
when two are, they are consecutive in the order.
:::

Special ordered sets were introduced by Beale and Tomlin (1970), together with
a branching rule that enforces them without binary variables: if the LP
solution violates the SOS2 condition, the solver branches on "$\alpha_i=0$
for all $i<r$" versus "$\alpha_i=0$ for all $i>r$" for a suitable index $r$.
A familiar SOS1 is a set of binaries with exactly one equal to $1$,
$\sum_i\delta_i=1$, $\delta_i\in\{0,1\}$: this is the multiple-choice
condition used below.

With SOS2 available, the convex-combination formulation of
$f^a$ is

$$
\begin{aligned}
x&=\sum_{i=1}^n\alpha_ia_i,\\
f^a&=\sum_{i=1}^n\alpha_if(a_i),\\
\sum_{i=1}^n\alpha_i&=1,\qquad \alpha_i\ge0,\\
&(\alpha_1,\dots,\alpha_n)\ \text{is SOS2} .
\end{aligned}
$$ (eq-07-03-sos2)

Without solver support, adjacency is enforced by a binary $\delta_k$ for each
segment ("$x$ lies in segment $k$"), with the conventions $\delta_0=\delta_n=0$:

$$
\begin{aligned}
&\alpha_i\le\delta_{i-1}+\delta_i,\quad i=1,\dots,n,\\
&\sum_{k=1}^{n-1}\delta_k=1,\\
&\delta_k\in\{0,1\},\quad k=1,\dots,n-1 .
\end{aligned}
$$ (eq-07-03-lambda)

Exactly one segment $k$ is chosen; then only $\alpha_k$ and $\alpha_{k+1}$
may be positive, since every other weight has both neighboring binaries
equal to $0$. Written out at the ends, the first and last rows are
$\alpha_1\le\delta_1$ and $\alpha_n\le\delta_{n-1}$. The objective term is
$\sum_{i=1}^n\alpha_if(a_i)$, one term per breakpoint. *(A double sum
$\sum_{i=1}^{n-1}\bigl[\alpha_if(a_i)+\alpha_{i+1}f(a_{i+1})\bigr]$ counts the
inner breakpoints twice and is incorrect. The segment indicators must be
binary; with $\delta_k\in[0,1]$ the rows {eq}`eq-07-03-lambda` are just the LP
relaxation.)*

(sec-07-03-incremental)=
## The incremental formulation

The **incremental** (or $\delta$-) formulation fills the segments from left
to right. Let $\xi_k\in[0,1]$ be the filled fraction of segment $k$:

$$
\begin{aligned}
x&=a_1+\sum_{k=1}^{n-1}h_k\,\xi_k,\\
f^a&=f(a_1)+\sum_{k=1}^{n-1}\bigl(f(a_{k+1})-f(a_k)\bigr)\xi_k,\\
\xi_{k+1}&\le\delta_k\le\xi_k,\quad k=1,\dots,n-2,\\
0&\le\xi_k\le1,\quad k=1,\dots,n-1,\\
\delta_k&\in\{0,1\},\quad k=1,\dots,n-2 .
\end{aligned}
$$ (eq-07-03-incremental)

The binary $\delta_k$ says "segment $k$ is full". If $\delta_k=1$ then
$\xi_k=1$; if $\delta_k=0$ then $\xi_{k+1}=0$. Hence the fractions have the
form $(1,\dots,1,\xi_k,0,\dots,0)$, which is the interpolant
{eq}`eq-07-03-interp` on segment $k$. The formulation needs $n-2$ binaries. In
the convex case the binaries can again be dropped (this is the segment form of
{ref}`sec-07-03-convex-case`). The name clashes with this book's use of $\delta$
for binaries; textbooks that call it the $\delta$-formulation write the
fractions as $\delta_k$ and the binaries as $y_k$.

(sec-07-03-multiple-choice)=
## The multiple-choice formulation

The third classical formulation chooses one segment and keeps a copy
$\bar x_k$ of $x$ for it:

$$
\begin{aligned}
x&=\sum_{k=1}^{n-1}\bar x_k,\\
f^a&=\sum_{k=1}^{n-1}\Bigl(s_k\bar x_k+\bigl(f(a_k)-s_ka_k\bigr)\delta_k\Bigr),\\
a_k\delta_k&\le\bar x_k\le a_{k+1}\delta_k,\quad k=1,\dots,n-1,\\
\sum_{k=1}^{n-1}\delta_k&=1,\qquad \delta_k\in\{0,1\} .
\end{aligned}
$$ (eq-07-03-multiple-choice)

If $\delta_k=1$ only $\bar x_k$ can be positive, it lies in segment $k$, and
the cost is the line of that segment. This is the disaggregated
(convex-hull) formulation of the disjunction "$x$ lies in segment 1 or 2 or
$\dots$" from {ref}`sec-07-01-either-or`.

(sec-07-03-comparison)=
## Comparing the formulations

{numref}`tab-07-03-compare` summarizes the formulations for one function with
$n$ breakpoints.

```{list-table} Formulations of a piecewise-linear function with $n$ breakpoints
:name: tab-07-03-compare
:header-rows: 1

* - Formulation
  - Continuous variables
  - Binaries
  - Remarks
* - convex combination, LP {eq}`eq-07-03-cc`
  - $n$
  - $0$
  - exact only in the convex case
* - convex combination with SOS2 {eq}`eq-07-03-sos2`
  - $n$
  - $0$
  - adjacency by SOS2 branching
* - convex combination with binaries {eq}`eq-07-03-lambda`
  - $n$
  - $n-1$
  - relaxation vertices may have fractional $\delta$
* - incremental {eq}`eq-07-03-incremental`
  - $n-1$
  - $n-2$
  - relaxation vertices have integral $\delta$
* - multiple choice {eq}`eq-07-03-multiple-choice`
  - $n-1$
  - $n-1$
  - relaxation vertices have integral $\delta$
```

Three observations, all checked numerically on the examples below.

1. **Bounds.** The LP relaxations of the three binary formulations give the
   same bound: each replaces $f^a$ by its convex envelope on $[a_1,a_n]$. In
   {prf:ref}`ex-07-03-economies` all three give $\$17.84$ million against the
   true optimum $\$19.89$ million.
2. **Integrality of the relaxation.** The incremental and multiple-choice
   formulations are *locally ideal*: every vertex of their LP relaxation has
   integral binaries (Vielma et al., 2010). The
   convex-combination formulation with binaries is not: in an experiment
   with $3000$ random objectives over a single function with five
   breakpoints, $65$ optimal vertices of its relaxation
   had fractional $\delta$, and none of those of the incremental or the
   multiple-choice formulation did.
   Branch and bound tends to need fewer nodes with locally ideal
   formulations.
3. **Size and convenience.** The convex-combination form is the easiest to
   write and to extend to several variables ({ref}`sec-07-03-two-variables`);
   with SOS2 support it adds no binaries at all.

(sec-07-03-error)=
## Approximating a nonlinear function

Approximating a smooth $f$ requires choosing the breakpoints. The error of
linear interpolation is controlled by the curvature.

:::{prf:theorem} Interpolation error
:label: thm-07-03-error

Let $f$ be twice continuously differentiable on $[a_k,a_{k+1}]$ with
$\abs{f''(x)}\le \kappa_k$ there. Then for every $x$ in the segment

$$
\abs{f(x)-f^a(x)}\le\frac{h_k^2}{8}\,\kappa_k .
$$ (eq-07-03-error)
:::

:::{dropdown} Proof
:class-container: proof
:animate: fade-in

Fix $x\in(a_k,a_{k+1})$ and let $e(t)=f(t)-f^a(t)$, which vanishes at $a_k$
and $a_{k+1}$. Define
$w(t)=e(t)-e(x)\,\dfrac{(t-a_k)(t-a_{k+1})}{(x-a_k)(x-a_{k+1})}$. Then
$w$ vanishes at $a_k$, $x$ and $a_{k+1}$, so by Rolle's theorem applied
twice there is $t^*$ in the segment with $w''(t^*)=0$. Since
$(f^a)''=0$ on the segment, $w''(t^*)=f''(t^*)-2e(x)/[(x-a_k)(x-a_{k+1})]$,
hence

$$
e(x)=\tfrac12f''(t^*)\,(x-a_k)(x-a_{k+1}) .
$$

The factor $\abs{(x-a_k)(x-a_{k+1})}$ is largest at the midpoint, where it
equals $h_k^2/4$. This gives {eq}`eq-07-03-error`.
:::

Consequences for modeling:

* To guarantee an error at most $\epsilon_f$ on a segment, take
  $h_k\le\sqrt{8\epsilon_f/\kappa_k}$. Halving the spacing divides the error by
  four.
* Place breakpoints densely where $\abs{f''}$ is large and sparsely where
  $f$ is nearly linear; equal spacing wastes breakpoints.
* Restrict the range $[a_1,a_n]$ to values the solution can actually take;
  for a steep function the error near the end of a wide range dominates the
  bound.
* For convex $f$ the interpolant overestimates ($f^a\ge f$), so a
  minimization with $f^a$ gives an upper bound on the true optimum, and its
  solution is feasible for the original problem when $f$ appears only in the
  objective; for concave $f$ the interpolant underestimates.

(sec-07-03-separable)=
## Separable programming

A function is **separable** if it is a sum of functions of one variable,
$f(\vx)=\sum_{j=1}^nf_j(x_j)$. A **separable program** has a separable
objective and separable constraints $g_i(\vx)=\sum_jg_{ij}(x_j)\le0$. Each
one-variable function is replaced by its interpolant on its own breakpoints,
with its own weights; the result is an LP with SOS2 conditions (Williams,
2013, Chapter 7; Bazaraa et al., 2006, §11.3).

* If every $f_j$ is convex (minimization) and every $g_{ij}$ is convex, the
  adjacency conditions can be dropped by {prf:ref}`prop-07-03-convex-case` and
  its analog for constraints (a convex $g_{ij}$ in a $\le$ constraint is
  over-estimated by any non-adjacent combination, which only makes the
  constraint harder to satisfy): the separable program is solved as an LP.
* Otherwise, the SOS2 conditions are enforced by SOS2 branching, by the
  binaries of {eq}`eq-07-03-lambda`, or, historically, by the **restricted basis
  entry** rule of the simplex method, which lets a weight enter the basis only
  next to a basic neighbor. The last approach finds a local optimum of the
  approximation only.

Some nonseparable terms can be made separable by a change of variables. The
classical example is a product: with $u_1=\tfrac12(x+y)$ and
$u_2=\tfrac12(x-y)$,

$$
xy=u_1^2-u_2^2,
$$

a difference of two one-variable functions; each square is then approximated
piecewise-linearly (the term $-u_2^2$ is concave and needs adjacency). For
products of bounded variables the McCormick envelope of
{ref}`sec-07-02-mccormick` is an alternative.

(sec-07-03-examples)=
## Two examples

The first example is convex and is solved as an LP; the second is concave and
needs binaries. Both are compared with the exact optimum.

:::{prf:example} System-optimal routing with congestion
:label: ex-07-03-congestion

A demand of $D=6000$ vehicles per hour travels from a suburb to the city
center on two parallel routes. Travel times follow the Bureau of Public
Roads (BPR) link performance function used in traffic assignment (Sheffi, 1985),

$$
t_r(v)=t_r^0\Bigl[1+0.15\Bigl(\frac{v}{c_r}\Bigr)^4\Bigr]\ \text{min},
$$

with free-flow times $t_1^0=10$ min, $t_2^0=15$ min and capacities
$c_1=3000$ veh/h, $c_2=4000$ veh/h. Find the split $v_1+v_2=D$ that
minimizes the total travel time $T_1(v_1)+T_2(v_2)$, where
$T_r(v)=v\,t_r(v)/60$ is in vehicle-hours per hour. Compute the exact
optimum, then solve the piecewise-linear LP with breakpoints every
$1000$ veh/h on both routes and compare the two.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Exact optimum.* $T_r$ is convex. Setting the derivatives equal (the
marginal travel times
$\mathrm d(v\,t_r)/\mathrm dv=t_r^0[1+0.75(v/c_r)^4]$) gives
$v_1^*=3070.9$ and $v_2^*=2929.1$ veh/h, with marginal times $18.235$ min on
both routes, actual times $11.647$ and $15.647$ min, and a total of
$1359.97$ veh-h/h.

*Piecewise-linear LP.* With breakpoints $0,1000,\dots,6000$ veh/h on both
routes, the interpolated values are

| $v$ (veh/h) | 0 | 1000 | 2000 | 3000 | 4000 | 5000 | 6000 |
|---|---|---|---|---|---|---|---|
| $T_1(v)$ (veh-h/h) | 0 | 166.98 | 343.21 | 575.00 | 982.72 | 1797.84 | 3400.00 |
| $T_2(v)$ (veh-h/h) | 0 | 250.15 | 504.69 | 785.60 | 1150.00 | 1707.76 | 2639.06 |

and the LP is

$$
\begin{aligned}
\min\ & \sum_{i=1}^7\bigl(\alpha_{1i}T_1(a_i)+\alpha_{2i}T_2(a_i)\bigr)\\
\st\ & \sum_{i=1}^7a_i(\alpha_{1i}+\alpha_{2i})=6000,\\
& \sum_{i=1}^7\alpha_{1i}=1,\\
& \sum_{i=1}^7\alpha_{2i}=1,\\
& \alpha_{ri}\ge0,\qquad r=1,2,\ i=1,\dots,7,
\end{aligned}
$$

with $a_i=1000(i-1)$ and no binaries, since both functions are convex. Its
solution is $\alpha_{14}=\alpha_{24}=1$, that is, $v_1=v_2=3000$ veh/h, with
value $575.00+785.60=1360.60$ veh-h/h. Check by slopes (marginal times in
minutes per vehicle): moving flow from route 2 to route 1 costs $24.46$ on
route 1's segment $[3000,4000]$ and saves $16.85$ on route 2's segment
$[2000,3000]$; moving flow the other way costs $21.86$ and saves $13.91$.
Neither move pays, so the LP is optimal.

*Comparison.* The approximate split is $71$ veh/h away from the exact one,
and its true total time, $1360.60$ veh-h/h (exact, since both flows are at
breakpoints), exceeds the optimum by only $0.63$ veh-h/h ($0.05\%$). The
interpolation itself is much less accurate: by {eq}`eq-07-03-error` with
$\kappa_k=T_1''(6000)=\tfrac{10}{60}\cdot3\cdot6000^3/3000^4=\tfrac1{750}$ the error
bound on route 1 is $1000^2/8\cdot\tfrac1{750}=166.7$ veh-h/h, and the actual
largest error is $129.2$ near $v=5523$. Near the optimum, on
$[2000,4000]$, the largest error is $33.6$. The objective is flat at its
minimum, so a coarse approximation still locates a near-optimal split; with
breakpoints every $500$ veh/h the errors fall roughly fourfold, to $36.7$
and $10.2$.
:::

:::{prf:example} Treatment plants with economies of scale
:label: ex-07-03-economies

A region needs $100$ ML/d (megaliters per day) of treatment capacity, to be
split between two plants with capacities of at most $80$ ML/d each. The
annualized costs, in \$ million per year, show economies of scale:
$C_1(q)=q^{0.6}$ for plant 1 and $C_2(q)=0.8\,q^{0.6}+0.06\,q$ for plant 2
(cheaper to build but with a pumping cost of \$0.06 million per year per
ML/d). Minimize $C_1(q_1)+C_2(q_2)$ subject to $q_1+q_2=100$,
$0\le q_1,q_2\le80$. Find the exact optimum; then approximate both costs
with breakpoints every $20$ ML/d and solve the convex-combination model
without and with adjacency conditions.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

*Exact optimum.* Both costs are concave, so the total cost along the
feasible segment $q_1\in[20,80]$ is concave and its minimum is at an end:
$C_1(80)+C_2(20)=13.863+6.027=\$19.89$ million against
$C_1(20)+C_2(80)=\$21.92$ million. The optimum is $q_1^*=80$, $q_2^*=20$.

*Breakpoints.* With $a_i=0,20,40,60,80$ ML/d,

| $q$ (ML/d) | 0 | 20 | 40 | 60 | 80 |
|---|---|---|---|---|---|
| $C_1(q)$ (\$M/yr) | 0 | 6.034 | 9.146 | 11.665 | 13.863 |
| $C_2(q)$ (\$M/yr) | 0 | 6.027 | 9.717 | 12.932 | 15.890 |

The slopes decrease (for $C_1$: $0.302$, $0.156$, $0.126$, $0.110$), so the
interpolants are concave.

*LP without adjacency.* Model {eq}`eq-07-03-cc` for both plants, with
$\sum_ia_i(\alpha_{1i}+\alpha_{2i})=100$, gives a cost of \$17.84 million
with $\alpha_{15}=1$ ($q_1=80$) and $\alpha_{21}=0.75$, $\alpha_{25}=0.25$: it
represents $q_2=20$ as a mixture of $0$ and $80$ and charges
$0.25\cdot15.890=3.97$ instead of $C_2(20)=6.03$. The weights are not
adjacent, and the cost is that of the convex envelope.

*MILP.* Adding the binaries {eq}`eq-07-03-lambda` for each plant (four per
plant) gives $\alpha_{15}=1$, $\alpha_{22}=1$, that is, $q_1=80$, $q_2=20$ at
\$19.89 million, the exact optimum (the optimum lies on breakpoints, so the
approximation error is zero there). The incremental and multiple-choice
formulations give the same solution. The LP relaxations of all three MILPs
give the same bound, \$17.84 million, the value of the LP without adjacency:
$C_1(80)+20\cdot C_2(80)/80$.
:::

{numref}`fig-07-03-sos2` illustrates the difference between adjacent and
non-adjacent weights found in {prf:ref}`ex-07-03-economies` on plant 1
alone: the point
$q=50$ is represented with SOS2 weights $\alpha_3=\alpha_4=0.5$ on the
breakpoints $40$ and $60$, giving $f^a(50)=10.41$ (true value $10.46$), or,
without adjacency, with $\alpha_1=0.375$, $\alpha_5=0.625$ on $0$ and $80$,
giving the chord value $8.66$.

```{figure} figures/07-03-sos2.svg
:name: fig-07-03-sos2
:width: 660px
:class: dark-light
:alt: Two representations of the same point as weights on breakpoints, one with adjacent weights on the curve, one with non-adjacent weights on the chord

Left: the concave cost $C_1(q)=q^{0.6}$, its interpolant, the chord from
$0$ to $80$ (dashed) and the two representations of $q=50$. Right: the weights $\alpha_i$. The SOS2 weights
(two adjacent nonzero) reproduce the interpolant; the non-adjacent weights
put the point on the chord below the curve.
```

(sec-07-03-two-variables)=
## Functions of two variables

A nonlinear function $f(x,y)$ is approximated on a grid
$a_1<\dots<a_n$ in $x$ and $b_1<\dots<b_m$ in $y$, which divides the domain
into $(n-1)(m-1)$ rectangles. Four approaches are in use (D'Ambrosio et
al., 2010).

**One-dimensional method.** Fix a representative $\tilde b_j\in[b_j,b_{j+1}]$
for each $y$-interval (usually $\tilde b_j=b_j$) and approximate
$f(x,\tilde b_j)$ in $x$ by the convex-combination formulation with weights
$\alpha_i$ (SOS2). Binaries $\beta_j$, $j=1,\dots,m-1$, select the
$y$-interval:

$$
\begin{aligned}
&\textstyle\sum_j\beta_jb_j\le y\le\sum_j\beta_jb_{j+1},\\
&\textstyle\sum_{j=1}^{m-1}\beta_j=1,\qquad \beta_j\in\{0,1\},\\
&\textstyle f^a\le\sum_i\alpha_if(a_i,\tilde b_j)+M(1-\beta_j),\quad j=1,\dots,m-1,\\
&\textstyle f^a\ge\sum_i\alpha_if(a_i,\tilde b_j)-M(1-\beta_j),\quad j=1,\dots,m-1 .
\end{aligned}
$$

For the chosen $j$ the two big-$M$ rows give
$f^a=\sum_i\alpha_if(a_i,\tilde b_j)$; the others are relaxed. The
approximation is piecewise linear in $x$ but piecewise *constant* in $y$, so
it is coarse. At a grid value $y=b_j$ two intervals qualify and the model
may pick the more favorable value; replacing the upper row by
$y\le\sum_j\beta_jb_{j+1}-\varepsilon$ removes the ambiguity, at the price of
excluding a sliver of width $\varepsilon$. *(The selection condition on the
$\beta_j$ is $\sum_j\beta_j=1$, an SOS1; a sum equal to $0$ would select no
interval.)*

**Triangle method.** Divide each rectangle by its diagonal from
$(a_i,b_j)$ to $(a_{i+1},b_{j+1})$ into an upper and a lower triangle and
interpolate linearly on each triangle, using its three vertices. With weights
$\alpha_{ij}\ge0$ on the grid points,

$$
\begin{aligned}
\sum_{i,j}\alpha_{ij}&=1,\\
x&=\sum_{i,j}\alpha_{ij}a_i,\\
y&=\sum_{i,j}\alpha_{ij}b_j,\\
f^a&=\sum_{i,j}\alpha_{ij}f(a_i,b_j),
\end{aligned}
$$

and the nonzero weights must lie on the vertices of one triangle. With
binaries $\delta^u_{ij}$ and $\delta^l_{ij}$ for the upper and lower triangle
of rectangle $(i,j)$ (zero when an index is out of range),

$$
\begin{aligned}
&\sum_{i=1}^{n-1}\sum_{j=1}^{m-1}\bigl(\delta^u_{ij}+\delta^l_{ij}\bigr)=1,\\
&\alpha_{ij}\le\delta^u_{ij}+\delta^l_{ij}+\delta^u_{i-1,j-1}+\delta^l_{i-1,j-1}\\
&\qquad+\delta^u_{i,j-1}+\delta^l_{i-1,j},\quad i=1,\dots,n,\ j=1,\dots,m .
\end{aligned}
$$

The six terms are the triangles that have $(a_i,b_j)$ as a vertex: both
triangles of the rectangle to its upper right and of the rectangle to its
lower left, the upper triangle of the rectangle to its lower right and the
lower triangle of the rectangle to its upper left ({numref}`fig-07-03-triangles`). The
approximation is continuous and exact at grid points, but it needs $nm$
weights and $2(n-1)(m-1)$ binaries. *(Some texts call the weight condition
"SOS3"; it is a triangle-adjacency condition, not a standard special ordered
set.)*

```{figure} figures/07-03-triangles.svg
:name: fig-07-03-triangles
:width: 420px
:class: dark-light
:alt: A triangulated grid with the six triangles around one grid point shaded

Triangulation of a grid by the diagonals from $(a_i,b_j)$ to
$(a_{i+1},b_{j+1})$. The weight $\alpha_{ij}$ of the marked grid point may be
positive only if the chosen triangle is one of the six shaded triangles
around it.
```

**Rectangle method.** Interpolate linearly in $x$ along the lower edge of
the chosen rectangle and add a linear correction in $y$. With
$\gamma=(y-b_j)/(b_{j+1}-b_j)$ and $\Delta(i,j)=f(a_i,b_{j+1})-f(a_i,b_j)$,

$$
\begin{aligned}
f^a(x,y)={}&\alpha f(a_i,b_j)+(1-\alpha)f(a_{i+1},b_j)\\
&+\gamma\min\bigl\{\Delta(i,j),\Delta(i+1,j)\bigr\}
\end{aligned}
$$

for $x=\alpha a_i+(1-\alpha)a_{i+1}$. If $\Delta(i,j)\le\Delta(i+1,j)$ this is
the plane through the corners $(a_i,b_j)$, $(a_{i+1},b_j)$ and
$(a_i,b_{j+1})$; otherwise the plane through $(a_i,b_j)$, $(a_{i+1},b_j)$ and
$(a_{i+1},b_{j+1})$. *(The case distinction compares the $y$-increments at the
two $x$-ends, $\Delta(i,j)$ and $\Delta(i+1,j)$.)* In the MILP, binaries
$\delta^x_i$ and $\beta_j$ select the $x$- and $y$-interval, continuous
$\gamma_j\le\beta_j$ carry the $y$-fraction through
$y=\sum_j(\beta_jb_j+\gamma_j(b_{j+1}-b_j))$, and two big-$M$ rows per
rectangle,

$$
\begin{aligned}
f^a&\le\sum_k\alpha_kf(a_k,b_j)+\gamma_j\Delta^{\min}_{ij}+M(2-\beta_j-\delta^x_i),\\
f^a&\ge\sum_k\alpha_kf(a_k,b_j)+\gamma_j\Delta^{\min}_{ij}-M(2-\beta_j-\delta^x_i),
\end{aligned}
$$

with $\Delta^{\min}_{ij}=\min\{\Delta(i,j),\Delta(i+1,j)\}$, fix $f^a$ for the chosen
rectangle. *(The relaxing term has opposite signs in the two rows: $+M$ in the
upper bound and $-M$ in the lower bound.)* Introducing
$\bar f_j=\sum_k\alpha_kf(a_k,b_j)$ as separate variables shortens these rows.
The method needs only $O(n+m)$ variables and binaries but $O(nm)$ big-$M$
rows.

**Tangent-plane method.** Wang and Lo (2010) linearize a travel-time function
of two variables by its first-order Taylor expansion at the lower-left corner
of each rectangle:

$$
\begin{aligned}
f^a(x,y)&=f(a_i,b_j)+f_x(a_i,b_j)(x-a_i)+f_y(a_i,b_j)(y-b_j)\\
&=p_{ij}x+q_{ij}y+r_{ij},
\end{aligned}
$$

with $p_{ij}=f_x(a_i,b_j)$, $q_{ij}=f_y(a_i,b_j)$ and
$r_{ij}=f(a_i,b_j)-p_{ij}a_i-q_{ij}b_j$. *(The constant uses the grid point
$(a_i,b_j)$, not the point $(x,y)$ being approximated.)* The approximation is
not continuous across rectangle boundaries. The rectangle is located by
binaries $\chi^x_i=1$ if $x<a_i$ and $\chi^y_j=1$ if $y<b_j$:

$$
\begin{aligned}
x&\ge a_i-M\chi^x_i, & i&=1,\dots,n,\\
x&\le a_i-\varepsilon+M(1-\chi^x_i), & i&=1,\dots,n,
\end{aligned}
$$

and the same rows for $y$ with $\chi^y_j$. Because the $a_i$ increase, $\chi^x_i=1$
implies $\chi^x_{i+1}=1$, so $\zeta^x_i=\chi^x_{i+1}-\chi^x_i\in\{0,1\}$ equals $1$
exactly when $a_i\le x<a_{i+1}$; likewise $\zeta^y_j=\chi^y_{j+1}-\chi^y_j$ for
$y$. *(The combination $\chi^x_i=\chi^x_{i+1}=1$ means $x<a_i$; $\chi^x_i=1$,
$\chi^x_{i+1}=0$ is impossible.)* Finally
$\abs{f^a-(p_{ij}x+q_{ij}y+r_{ij})}\le M(2-\zeta^x_i-\zeta^y_j)$, written as two
rows, fixes $f^a$ on the selected rectangle. Fixing $\chi^x_1=0$ and $\chi^x_n=1$
restricts $x$ to $[a_1,a_n)$. An equivalent model writes
$f^a=\sum_{i,j}(p_{ij}x+q_{ij}y+r_{ij})\zeta^x_i\zeta^y_j$, whose products of
binaries and continuous variables are linearized as in
{ref}`sec-07-02-binary-continuous`.

```{list-table} Size of the two-variable formulations on an $n\times m$ grid (orders of magnitude)
:name: tab-07-03-two-var
:header-rows: 1

* - Method
  - Continuous variables
  - Binaries
  - Constraints
  - Approximation
* - one-dimensional
  - $n$
  - $(n-1)+(m-1)$
  - $O(n+m)$
  - linear in $x$, constant in $y$
* - triangle
  - $nm$
  - $2(n-1)(m-1)$
  - $nm+5$
  - continuous, linear on triangles
* - rectangle
  - $n+(m-1)$
  - $(n-1)+(m-1)$
  - $O(nm)$
  - linear on half-rectangles
* - tangent plane
  - $(n-1)+(m-1)$
  - $n+m$
  - $O(nm)$
  - first-order Taylor, discontinuous
```

{numref}`tab-07-03-two-var` compares the sizes (counting the SOS2 weights of the
one-dimensional and rectangle methods as continuous variables and their
adjacency binaries among the binaries). The triangle method is the most
accurate and the largest; the one-dimensional method is the smallest and the
coarsest; the rectangle method is a compromise.

## Summary and Common Mistakes

### Summary

- The interpolant on breakpoints $a_1<\dots<a_n$ is {eq}`eq-07-03-interp`; it is
  convex when the slopes increase and then equals a maximum of affine pieces.
- Minimizing a convex (or maximizing a concave) interpolant is an LP: the
  weights are adjacent at an optimum ({prf:ref}`prop-07-03-convex-case`).
- Nonconvex cases need SOS2 weights ({prf:ref}`def-07-03-sos`,
  {eq}`eq-07-03-sos2`), binary adjacency {eq}`eq-07-03-lambda`, the incremental
  formulation {eq}`eq-07-03-incremental` or the multiple-choice formulation
  {eq}`eq-07-03-multiple-choice` ({numref}`tab-07-03-compare`,
  {numref}`fig-07-03-sos2`).
- Their LP relaxations give the convex-envelope bound; the incremental and
  multiple-choice formulations have integral relaxation vertices.
- The interpolation error is at most $h_k^2\kappa_k/8$ ({prf:ref}`thm-07-03-error`);
  place breakpoints where the curvature is large.
- Separable programs are approximated function by function; convex ones are
  LPs.
- Functions of two variables are approximated by the one-dimensional,
  triangle, rectangle or tangent-plane method ({numref}`tab-07-03-two-var`).
- In {prf:ref}`ex-07-03-congestion` a coarse LP comes within $0.05\%$ of the
  exact optimum; in {prf:ref}`ex-07-03-economies` only the MILP finds the true
  optimum.

### Common mistakes

- **Dropping adjacency in a nonconvex model.** The LP then optimizes the
  convex envelope and its "solution" mixes distant breakpoints
  ({prf:ref}`ex-07-03-economies`).
- **Adding binaries in a convex model.** They are harmless but useless:
  minimizing a convex interpolant is an LP ({prf:ref}`prop-07-03-convex-case`).
- **Summing the objective over segments** instead of over breakpoints in the
  convex-combination formulation, which double-counts inner breakpoints.
- **Declaring the segment indicators continuous.** Then {eq}`eq-07-03-lambda` is
  only a relaxation.
- **Breakpoints over too wide a range.** The error bound {eq}`eq-07-03-error` is
  governed by the largest curvature on each segment; restrict the range or
  refine where $\abs{f''}$ is large.
- **Ignoring boundary ties** in the two-variable models with big-$M$ rows; use
  $\varepsilon$ or accept the more favorable value at grid lines.

## References

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Section 11.3:
  separable programming.)
- Beale, E. M. L., & Tomlin, J. A. (1970). Special facilities in a general
  mathematical programming system for nonconvex problems using ordered sets
  of variables. In J. Lawrence (Ed.), *Proceedings of the Fifth International
  Conference on Operational Research* (pp. 447–454). Tavistock.
- D'Ambrosio, C., Lodi, A., & Martello, S. (2010). Piecewise linear
  approximation of functions of two variables in MILP models. *Operations
  Research Letters, 38*(1), 39–46.
  <https://doi.org/10.1016/j.orl.2009.09.005>
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis with
  mathematical programming methods*. Prentice-Hall. (Chapter 1: link
  performance functions; Chapter 3: the user-equilibrium and system-optimum
  formulations.)
- Vielma, J. P., Ahmed, S., & Nemhauser, G. (2010). Mixed-integer models for
  nonseparable piecewise-linear optimization: Unifying framework and
  extensions. *Operations Research, 58*(2), 303–315.
- Wang, D. Z. W., & Lo, H. K. (2010). Global optimum of the linearized
  network design problem with equilibrium flows. *Transportation Research
  Part B: Methodological, 44*(4), 482–492.
  <https://doi.org/10.1016/j.trb.2009.10.003>
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapter 7: nonlinear models and separable programming;
  Section 9.3: special ordered sets of variables.)
