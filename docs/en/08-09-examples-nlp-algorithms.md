# Nonlinear Programming Algorithms

This chapter works through the examples and problems of the last three
chapters on nonlinear programming: {doc}`06-04-line-search`,
{doc}`06-05-unconstrained-algorithms` and {doc}`06-06-constrained-algorithms`.
Iteration tables were recomputed in exact or double-precision arithmetic and
are printed to four significant digits unless a comparison needs more.

(sec-08-09-one-dimensional)=
## One-dimensional search

:::{prf:example} Three line searches on a polynomial
:label: ex-08-09-one-dim-sextic

Minimize $f(x)=2x^6+3x^4-12x$ on the interval $[0,1]$ with accuracy
$\varepsilon=0.1$ by

**(a)** bisection search (stop when $b_k-a_k\le0.1$);

**(b)** golden-section search (stop when $b_k-a_k\le0.1$);

**(c)** Newton's method from $x_0=1$ (stop when $\abs{f'(x_k)}\le0.1$).
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The derivatives are $f'(x)=12x^5+12x^3-12$ and $f''(x)=60x^4+36x^2$. Since
$f''>0$ for $x\ne0$, $f$ is strictly convex on $[0,1]$; $f'(0)=-12<0$ and
$f'(1)=12>0$, so the minimizer is interior. For reference, solving $f'(x)=0$
to machine precision gives $x^*=0.837620$ and $f^*=-7.88395$.

**(a)** Bisection ({prf:ref}`alg-06-04-bisection`) needs
$\lceil\log_2(1/0.1)\rceil=4$ halvings.

```{list-table} Bisection on $f'(x)=12x^5+12x^3-12$, $[0,1]$
:header-rows: 1

* - $k$
  - $a_k$
  - $b_k$
  - $b_k-a_k$
  - $x_k$
  - $f'(x_k)$
  - $f(x_k)$
* - 0
  - 0
  - 1
  - 1
  - 0.5
  - $-10.13$
  - $-5.781$
* - 1
  - 0.5
  - 1
  - 0.5
  - 0.75
  - $-4.090$
  - $-7.695$
* - 2
  - 0.75
  - 1
  - 0.25
  - 0.875
  - $2.194$
  - $-7.844$
* - 3
  - 0.75
  - 0.875
  - 0.125
  - 0.8125
  - $-1.314$
  - $-7.867$
* - 4
  - 0.8125
  - 0.875
  - 0.0625
  - stop
  -
  -
```

The estimate is $\hat x=0.84375$ with $f(\hat x)=-7.8829$; its error
$0.0061$ is within the guaranteed $\varepsilon/2=0.05$.

**(b)** Golden-section search ({prf:ref}`alg-06-04-golden-section`) with
$\tau=0.618034$ needs $\lceil\ln0.1/\ln\tau\rceil=5$ reductions and six
function values. Values are given to five decimals because rows 3 and 4
compare numbers that agree in the third decimal.

```{list-table} Golden-section search on $f(x)=2x^6+3x^4-12x$, $[0,1]$
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
  - 0
  - 1
  - 1
  - 0.38197
  - 0.61803
  - $-4.51352$
  - $-6.86726$
* - 1
  - 0.38197
  - 1
  - 0.61803
  - 0.61803
  - 0.76393
  - $-6.86726$
  - $-7.74793$
* - 2
  - 0.61803
  - 1
  - 0.38197
  - 0.76393
  - 0.85410
  - $-7.74793$
  - $-7.87635$
* - 3
  - 0.76393
  - 1
  - 0.23607
  - 0.85410
  - 0.90983
  - $-7.87635$
  - $-7.72778$
* - 4
  - 0.76393
  - 0.90983
  - 0.14590
  - 0.81966
  - 0.85410
  - $-7.87530$
  - $-7.87635$
* - 5
  - 0.81966
  - 0.90983
  - 0.09017
  - stop
  -
  -
  -
```

In rows 0–2 the right trial point is better, so the left part is discarded
and the old $q_k$ is reused as $p_{k+1}$; in row 3 the left point is better,
and $(q_3,1]$ is discarded; in row 4, $f(p_4)>f(q_4)$ again. The final interval
$[0.81966,\ 0.90983]$ contains $x^*$, and the estimate is its midpoint
$\hat x=0.86475$ (error $0.027$). *(Rounding $\tau$ to $0.618$ and recomputing
both trial points in every row gives slightly different points (for example
$q_4=0.854116$, although the point carried over from row 3 is $0.854134$)
and the final interval
$[0.819692,\ 0.909837]$; its midpoint is $(0.819692+0.909837)/2=0.864765$.
A printed version of this calculation has the misprint $0.81692$ for
$0.819692$ in that last sum.)*

**(c)** Newton's method ({prf:ref}`alg-06-04-newton-1d`),
$x_{k+1}=x_k-f'(x_k)/f''(x_k)$:

```{list-table} Newton's method on $f(x)=2x^6+3x^4-12x$ from $x_0=1$
:header-rows: 1

* - $k$
  - $x_k$
  - $f'(x_k)$
  - $f''(x_k)$
  - $x_{k+1}$
* - 0
  - 1
  - 12
  - 96
  - 0.875
* - 1
  - 0.875
  - 2.194
  - 62.73
  - 0.840027
* - 2
  - 0.840027
  - 0.1325
  - 55.28
  - 0.837630
* - 3
  - 0.837630
  - $5.810\times10^{-4}$
  - 54.80
  - stop
```

Newton's method stops after three steps at $x_3=0.837630$, which agrees with
$x^*=0.837620$ to five digits; the next step would give $0.83761977$. With a
derivative test the same tolerance $0.1$ is much more demanding than an
interval width of $0.1$, yet Newton's method meets it with three derivative
pairs, against four derivative values for bisection (error $0.006$) and six
function values for golden-section search (error $0.027$).
:::

:::{prf:example} Golden-section search on a nonsmooth function
:label: ex-08-09-golden-nonsmooth

Use golden-section search to find, within $\pm0.15$ of the minimizer, the
value of $x\in[0,8]$ that minimizes

$$
f(x)=\max\Bigl\{2-\frac x3,\ (x-3)^2,\ \frac x5\Bigr\}.
$$
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The function is a maximum of convex functions, hence convex and unimodal,
but it has kinks, so golden-section search (which uses no derivatives) is
appropriate. If the final interval has width at most $0.3$, its midpoint is
within $0.15$ of $x^*$; so the stopping rule is $b_k-a_k\le0.3$. This needs
$8\tau^k\le0.3$, that is $k\ge\ln(0.3/8)/\ln\tau=6.82$, so $7$ reductions.

```{list-table} Golden-section search on $\max\{2-x/3,(x-3)^2,x/5\}$, $[0,8]$
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
  - 0
  - 8
  - 8
  - 3.0557
  - 4.9443
  - 0.9814
  - 3.7802
* - 1
  - 0
  - 4.9443
  - 4.9443
  - 1.8885
  - 3.0557
  - 1.3705
  - 0.9814
* - 2
  - 1.8885
  - 4.9443
  - 3.0557
  - 3.0557
  - 3.7771
  - 0.9814
  - 0.7554
* - 3
  - 3.0557
  - 4.9443
  - 1.8885
  - 3.7771
  - 4.2229
  - 0.7554
  - 1.4955
* - 4
  - 3.0557
  - 4.2229
  - 1.1672
  - 3.5016
  - 3.7771
  - 0.8328
  - 0.7554
* - 5
  - 3.5016
  - 4.2229
  - 0.7214
  - 3.7771
  - 3.9474
  - 0.7554
  - 0.8975
* - 6
  - 3.5016
  - 3.9474
  - 0.4458
  - 3.6718
  - 3.7771
  - 0.7761
  - 0.7554
* - 7
  - 3.6718
  - 3.9474
  - 0.2755
  - stop
  -
  -
  -
```

The search stops with $[3.6718,\ 3.9474]$ and the estimate
$\hat x=3.8096$.

*Check.* Near the minimum the binding pieces are the two lines: $2-x/3$
decreases, $x/5$ increases, and they cross where $2-x/3=x/5$, at
$x^*=3.75$, with $f^*=0.75$; there $(x-3)^2=0.5625<0.75$, so the parabola is
not binding. The estimate $3.8096$ is within $0.060$ of $x^*$, inside the
required $\pm0.15$. *(A version of this solution stops as soon as
$\abs{f(p_k)-f(q_k)}\le0.15$, which happens at $k=4$ with the interval
$[3.056,\ 4.223]$ of width $1.167$, and states no final $x$. A small
difference of function values does not mean that the interval is short;
only the interval width controls the error in $x$.)*
:::

:::{prf:example} Newton's method on a cubic
:label: ex-08-09-newton-cubic

Find the minimum of $f(x)=2x^3+3x^2-12x+5$ on $[0,4]$ by Newton's method,
starting from $x_0=0$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The derivatives are $f'(x)=6x^2+6x-12=6(x+2)(x-1)$ and $f''(x)=12x+6$. The
Newton iteration is

$$
x_{k+1}=x_k-\frac{6x_k^2+6x_k-12}{12x_k+6}.
$$

```{list-table} Newton's method on $f(x)=2x^3+3x^2-12x+5$ from $x_0=0$
:header-rows: 1

* - $k$
  - $x_k$
  - $f'(x_k)$
  - $f''(x_k)$
  - $x_{k+1}$
  - $f(x_k)$
* - 0
  - 0
  - $-12$
  - 6
  - 2
  - 5
* - 1
  - 2
  - 24
  - 30
  - 1.2
  - 9
* - 2
  - 1.2
  - 3.84
  - 20.4
  - 1.011765
  - $-1.624$
* - 3
  - 1.011765
  - 0.2126
  - 18.14
  - 1.000046
  - $-1.999$
* - 4
  - 1.000046
  - $8.240\times10^{-4}$
  - 18.00
  - 1.000000
  - $-2.000$
* - 5
  - 1.000000
  - $1.257\times10^{-8}$
  - 18.00
  - stop
  - $-2.000$
```

In exact arithmetic the iterates are $0$, $2$, $\tfrac65$, $\tfrac{86}{85}$, $\tfrac{21846}{21845},\dots$
The first step overshoots to $x_1=2$ (and $f$ increases from $5$ to $9$)
because the curvature $f''(0)=6$ is small; from there the iterates converge
quadratically to the stationary point $x^*=1$, where $f''(1)=18>0$. The
errors $1$, $0.2$, $0.0118$, $4.6\times10^{-5}$, $7.0\times10^{-10}$ have ratios
$e_{k+1}/e_k^2=0.20,\ 0.29,\ 0.33$, approaching
$f'''(1)/(2f''(1))=12/36=\tfrac13$ as {prf:ref}`thm-06-04-newton-quadratic`
predicts.

*Check of the interval.* The other stationary point $x=-2$ lies outside
$[0,4]$, and $f(0)=5$, $f(1)=-2$, $f(4)=133$, so $x^*=1$ with $f^*=-2$ is the
minimum on $[0,4]$. *(A printed table for this problem has $f'(2)=12$ and
$f''(2)=13$ and continues with $f''=14,15,\dots,20$, as if the column had
been filled as a series; all rows after the first are affected. The correct
values are $f'(2)=24$, $f''(2)=30$ and those above; the final answer
$x^*=1$ is right.)*
:::

(sec-08-09-steepest-descent)=
## Steepest descent

:::{prf:example} Steepest descent in three variables
:label: ex-08-09-steepest-descent-3d

Perform four iterations of the steepest-descent method with exact line
search on

$$
\min_{\vx\in\R^3} f(\vx)=x_1^2+2x_2^2+x_3^2-2x_1-3x_2-2x_3,
$$

starting at $\vx^{(0)}=(0,0,0)$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The gradient is $\nabla f(\vx)=(2x_1-2,\ 4x_2-3,\ 2x_3-2)$, so $f$ is the
quadratic {eq}`eq-06-05-quadratic` with $\mQ=\diag(2,4,2)$ and $\vb=(2,3,2)$. With
$\vd^{(k)}=-\nabla f(\vx^{(k)})=(2-2x_1^{(k)},\ 3-4x_2^{(k)},\ 2-2x_3^{(k)})$,
setting the derivative of $\varphi(\alpha)=f(\vx^{(k)}+\alpha\vd^{(k)})$ to zero
gives the exact step {eq}`eq-06-05-sd-step`

$$
\alpha_k=\frac{\bigl(d_1^{(k)}\bigr)^2+\bigl(d_2^{(k)}\bigr)^2+\bigl(d_3^{(k)}\bigr)^2}
{2\bigl(d_1^{(k)}\bigr)^2+4\bigl(d_2^{(k)}\bigr)^2+2\bigl(d_3^{(k)}\bigr)^2}.
$$

At $\vx^{(0)}$: $\vd^{(0)}=(2,3,2)$, $\alpha_0=17/52=0.3269$ and
$\vx^{(1)}=\tfrac{17}{52}(2,3,2)=(0.6538,\ 0.9808,\ 0.6538)$.

```{list-table} Steepest descent on $x_1^2+2x_2^2+x_3^2-2x_1-3x_2-2x_3$ from $\vzero$
:header-rows: 1

* - $k$
  - $\vx^{(k)}$
  - $\vd^{(k)}$
  - $\alpha_k$
  - $f(\vx^{(k)})$
* - 0
  - $(0,\ 0,\ 0)$
  - $(2,\ 3,\ 2)$
  - 0.3269
  - 0
* - 1
  - $(0.6538,\ 0.9808,\ 0.6538)$
  - $(0.6923,\ -0.9231,\ 0.6923)$
  - 0.3400
  - $-2.779$
* - 2
  - $(0.8892,\ 0.6669,\ 0.8892)$
  - $(0.2215,\ 0.3323,\ 0.2215)$
  - 0.3269
  - $-3.087$
* - 3
  - $(0.9617,\ 0.7756,\ 0.9617)$
  - $(0.07668,\ -0.1022,\ 0.07668)$
  - 0.3400
  - $-3.121$
* - 4
  - $(0.9877,\ 0.7408,\ 0.9877)$
  -
  -
  - $-3.1245$
```

After four iterations $\vx^{(4)}=(0.9877,\ 0.7408,\ 0.9877)$ and
$f(\vx^{(4)})=-3.1245$. The exact minimizer solves $\mQ\vx=\vb$:
$\vx^*=(1,\ 0.75,\ 1)$ with $f^*=-3.125$.

*Checks.* Consecutive directions are orthogonal, as {eq}`eq-06-05-orthogonal`
requires: $(2,3,2)^\T(0.6923,-0.9231,0.6923)=1.3846-2.7692+1.3846=0$. The
steps alternate between $\tfrac{17}{52}$ and $\tfrac{17}{50}$, and every
iteration multiplies $f(\vx^{(k)})-f^*$ by $0.1108$ ($3.125$ → $0.3462$ → $0.03834$ → $0.004247$ → $0.0004705$).
With $\kappa=4/2=2$, {prf:ref}`thm-06-05-sd-rate` bounds this factor by
$((2-1)/(2+1))^2=\tfrac19=0.1111$; the start is close to the worst case.
:::

(sec-08-09-newton)=
## Newton's method

:::{prf:example} Newton's method with a singular Hessian at the solution
:label: ex-08-09-newton-singular-hessian

Apply Newton's method to

$$
\min_{\vx\in\R^2} f(\vx)=(x_1-2)^4+(x_1-2x_2)^2
$$

from $\vx^{(0)}=(0,3)$. Compute six iterations, find the minimizer, and
explain the speed of convergence.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

The gradient and Hessian are

$$
\begin{aligned}
\nabla f(\vx)&=\begin{pmatrix}4(x_1-2)^3+2(x_1-2x_2)\\-4(x_1-2x_2)\end{pmatrix},\\
\nabla^2f(\vx)&=\begin{pmatrix}12(x_1-2)^2+2&-4\\-4&8\end{pmatrix}.
\end{aligned}
$$

**Iteration 0.** At $(0,3)$: $f=16+36=52$, $\nabla f=(-44,\ 24)$ and
$\nabla^2f=\begin{pmatrix}50&-4\\-4&8\end{pmatrix}$ with determinant $384$.
Solving $\nabla^2f\,\vd=-\nabla f$,

$$
\vd^{(0)}=-\frac1{384}\begin{pmatrix}8&4\\4&50\end{pmatrix}\begin{pmatrix}-44\\24\end{pmatrix}
=\frac1{384}\begin{pmatrix}256\\-1024\end{pmatrix}
=\begin{pmatrix}0.6667\\-2.6667\end{pmatrix},
$$

so $\vx^{(1)}=(\tfrac23,\tfrac13)$ and $f(\vx^{(1)})=\tfrac{256}{81}=3.160$.

**Later iterations.** At $\vx^{(1)}$, $x_1=2x_2$, so the second gradient
component vanishes. Write $e=x_1-2$ and suppose $x_1=2x_2$. Then
$\nabla f=(4e^3,\ 0)$, and the Newton equations
$(12e^2+2)d_1-4d_2=-4e^3$, $-4d_1+8d_2=0$ give $d_2=d_1/2$ and
$12e^2d_1=-4e^3$, that is $d_1=-e/3$. The new point again satisfies
$x_1=2x_2$, and

$$
\begin{aligned}
x_1^{(k+1)}-2&=\tfrac23\bigl(x_1^{(k)}-2\bigr),\\
\vx^{(k)}&=\Bigl(2-2\bigl(\tfrac23\bigr)^k,\ 1-\bigl(\tfrac23\bigr)^k\Bigr)\quad(k\ge1),\\
f(\vx^{(k)})&=16\bigl(\tfrac23\bigr)^{4k}\quad(k\ge1).
\end{aligned}
$$

```{list-table} Newton's method on $(x_1-2)^4+(x_1-2x_2)^2$ from $(0,3)$
:header-rows: 1

* - $k$
  - $\vx^{(k)}$
  - $f(\vx^{(k)})$
  - $\nabla f(\vx^{(k)})$
  - $[\nabla^2f]_{11}$
  - $\vd^{(k)}$
* - 0
  - $(0,\ 3)$
  - 52
  - $(-44,\ 24)$
  - 50
  - $(0.6667,\ -2.6667)$
* - 1
  - $(0.6667,\ 0.3333)$
  - 3.160
  - $(-9.481,\ 0)$
  - 23.33
  - $(0.4444,\ 0.2222)$
* - 2
  - $(1.1111,\ 0.5556)$
  - 0.6243
  - $(-2.809,\ 0)$
  - 11.48
  - $(0.2963,\ 0.1481)$
* - 3
  - $(1.4074,\ 0.7037)$
  - 0.1233
  - $(-0.8324,\ 0)$
  - 6.214
  - $(0.1975,\ 0.09877)$
* - 4
  - $(1.6049,\ 0.8025)$
  - 0.02436
  - $(-0.2466,\ 0)$
  - 3.873
  - $(0.1317,\ 0.06584)$
* - 5
  - $(1.7366,\ 0.8683)$
  - $4.812\times10^{-3}$
  - $(-0.07308,\ 0)$
  - 2.832
  - $(0.08779,\ 0.04390)$
* - 6
  - $(1.8244,\ 0.9122)$
  - $9.505\times10^{-4}$
  - $(-0.02165,\ 0)$
  - 2.370
  -
```

**Minimizer and rate.** $f\ge0$ with equality only for $x_1=2$ and
$x_1=2x_2$, so $\vx^*=(2,1)$ and $f^*=0$. The iterates converge to it, but
only **linearly**: the error shrinks by the factor $\tfrac23$ per iteration,
and $\norm{\vx^{(k)}-\vx^*}<10^{-3}$ needs $k=20$. Quadratic convergence
({prf:ref}`thm-06-05-newton-quadratic`) is lost because the Hessian at the
solution, $\begin{pmatrix}2&-4\\-4&8\end{pmatrix}$, is **singular**
(eigenvalues $0$ and $10$): along the valley $x_1=2x_2$ the function behaves
like $(x_1-2)^4$, whose second derivative vanishes at the minimizer, and
Newton's method on $t^4$ gives $t_{k+1}=t_k-4t_k^3/(12t_k^2)=\tfrac23t_k$.

*(An earlier version of this table rounds every iterate to two
decimals before the next step. This creates spurious second gradient
components $\pm0.04$ and the value $f(\vx^{(1)})=3.13$ instead of $3.160$; it
also prints the step of iteration 1 as $(0.44,\ -0.23)$, whereas the step is
$(0.44,\ +0.22)$, and $f(\vx^{(4)})=0.12$, whereas $f(\vx^{(4)})=0.024$. The
exact values are those in the table.)*
:::

(sec-08-09-frank-wolfe)=
## Frank–Wolfe method

:::{prf:example} Frank–Wolfe on a quadratic over a triangle
:label: ex-08-09-fw-tutorial-qp

Consider the convex program

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)=x_1^2+x_2^2-3x_1-4x_2\\
\st\ & x_1+x_2\le1,\\
& x_1,x_2\ge0.
\end{aligned}
$$

**(a)** Starting from $\vx^{(0)}=(\tfrac14,\tfrac14)$, apply three iterations
of the Frank–Wolfe method.

**(b)** Use the Karush–Kuhn–Tucker (KKT) conditions to check whether the result of (a) is optimal.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** $\nabla f(\vx)=(2x_1-3,\ 2x_2-4)$, the Hessian is $2\mI$, and the
feasible set is the triangle with vertices $(0,0)$, $(1,0)$, $(0,1)$. For a
direction $\vd$, the line function
$\varphi(\alpha)=f(\vx+\alpha\vd)$ is a parabola with minimizer

$$
\alpha^*=\frac{-\nabla f(\vx)^\T\vd}{\vd^\T(2\mI)\vd}
=\frac{(3-2x_1)d_1+(4-2x_2)d_2}{2d_1^2+2d_2^2},
$$

and since $\varphi$ is convex the best step in $[0,1]$ is
$\alpha_k=\min\{1,\alpha^*\}$ (when $\alpha^*\ge0$), as in
{eq}`eq-06-06-fw-quadratic-step`.

*Iteration 0.* $\vx^{(0)}=(\tfrac14,\tfrac14)$, $f=-1.625$,
$\nabla f=(-\tfrac52,-\tfrac72)$. The linear program (LP) $\min\,-\tfrac52y_1-\tfrac72y_2$ has the
vertex values $0$, $-\tfrac52$, $-\tfrac72$, so $\vy^{(0)}=(0,1)$ and
$\vd^{(0)}=(-\tfrac14,\tfrac34)$. The gap is
$\Gamma=-\nabla f^\T\vd^{(0)}=-\tfrac58+\tfrac{21}8=2$. The parabola has its
minimizer at $\alpha^*=2/(2\cdot\tfrac{10}{16})=1.6>1$, so $\alpha_0=1$ and
$\vx^{(1)}=(0,1)$, $f=-3$.

*Iteration 1.* $\nabla f(\vx^{(1)})=(-3,-2)$; the LP values $0$, $-3$, $-2$
give $\vy^{(1)}=(1,0)$, $\vd^{(1)}=(1,-1)$, $\Gamma=3-2=1$, and
$\alpha_1=1/(2\cdot2)=\tfrac14$. Hence $\vx^{(2)}=(\tfrac14,\tfrac34)$ and
$f=-3.125$.

*Iteration 2.* $\nabla f(\vx^{(2)})=(-\tfrac52,-\tfrac52)$. Now the LP
objective $-\tfrac52(y_1+y_2)$ takes its minimum $-\tfrac52$ on the whole edge
$y_1+y_2=1$, which contains $\vx^{(2)}$ itself. Whichever optimal vertex is
chosen, for example $\vy^{(2)}=(1,0)$ with $\vd^{(2)}=(\tfrac34,-\tfrac34)$, the
gap is $\Gamma(\vx^{(2)})=-(-\tfrac52\cdot\tfrac34+\tfrac52\cdot\tfrac34)=0$, so
$\alpha_2=0$ and $\vx^{(3)}=\vx^{(2)}=(\tfrac14,\tfrac34)$. The zero gap is the
stopping test of {prf:ref}`alg-06-06-frank-wolfe`: by
{prf:ref}`thm-06-06-fw-gap`, $\vx^{(2)}$ is already optimal, with $f^*=-3.125$.
*(In one printed solution the direction of this iteration is labeled
$\vd^{(1)}=\vy^{(1)}-\vx^{(1)}$; it is $\vd^{(2)}=\vy^{(2)}-\vx^{(2)}$.)*

**(b)** The unconstrained minimizer $(\tfrac32,2)$ violates $x_1+x_2\le1$, so
try the constraint active with $x_1,x_2>0$ (multipliers of $x\ge0$ zero). With
multiplier $\lambda\ge0$ for $g(\vx)=x_1+x_2-1\le0$, the KKT conditions are

$$
\begin{aligned}
2x_1-3+\lambda&=0,\\
2x_2-4+\lambda&=0,\\
x_1+x_2&=1 .
\end{aligned}
$$

Subtracting the first two gives $x_2=x_1+\tfrac12$, hence $x_1=\tfrac14$,
$x_2=\tfrac34$ and $\lambda=\tfrac52\ge0$. All KKT conditions hold, and since
the problem is convex, $\vx^*=(\tfrac14,\tfrac34)$ is the global minimizer:
the Frank–Wolfe result is optimal.
:::

:::{prf:example} Frank–Wolfe with a graphically solved subproblem
:label: ex-08-09-fw-graphical

Use the Frank–Wolfe method, starting from $\vx^{(0)}=(0,-3)$, to solve

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)=4(x_1-10)^2+(x_2-4)^2\\
\st\ & x_1-x_2\le10,\\
& \tfrac15x_1-x_2\ge3,\\
& x_1\ge0.
\end{aligned}
$$

Plot the feasible region, solve each linear subproblem by inspection of the
plot, and draw the steps of the algorithm.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**Feasible region.** The lines $x_1-x_2=10$ and $\tfrac15x_1-x_2=3$ meet where
$x_1-10=\tfrac15x_1-3$, at $(\tfrac{35}4,-\tfrac54)$. With $x_1\ge0$ the
feasible set is the triangle with vertices $(0,-3)$, $(0,-10)$ and
$(\tfrac{35}4,-\tfrac54)$ ({numref}`fig-08-09-fw-graphical`a). The gradient is
$\nabla f(\vx)=(8(x_1-10),\ 2(x_2-4))$.

**Iteration 0.** $\nabla f(\vx^{(0)})=(-80,-14)$. The LP
$\min\,-80y_1-14y_2$ over the triangle is solved by sliding the line
$-80y_1-14y_2=c$ in the direction of decreasing $c$, that is towards larger
$y_1$ and $y_2$; the last vertex it touches is $(\tfrac{35}4,-\tfrac54)$ (vertex
values $42$, $140$ and $-682.5$). So $\vy^{(0)}=(\tfrac{35}4,-\tfrac54)$,
$\vd^{(0)}=(\tfrac{35}4,\tfrac74)$ and the gap is
$\Gamma=-(-80\cdot\tfrac{35}4-14\cdot\tfrac74)=724.5>0$. Along the segment,
$\vx^{(0)}+\alpha\vd^{(0)}=(\tfrac{35}4\alpha,\ -3+\tfrac74\alpha)$ and

$$
\begin{aligned}
\varphi(\alpha)&=4\bigl(\tfrac{35}4\alpha-10\bigr)^2+\bigl(\tfrac74\alpha-7\bigr)^2\\
&=309.3125\,\alpha^2-724.5\,\alpha+449 .
\end{aligned}
$$

Its unconstrained minimizer $\alpha^*=724.5/618.625=1.171$ exceeds $1$; on
$[0,1]$ the convex parabola is decreasing, so $\alpha_0=1$ and
$\vx^{(1)}=(\tfrac{35}4,-\tfrac54)$, with $f(\vx^{(1)})=\tfrac{541}{16}=33.81$
(down from $449$).

**Iteration 1.** $\nabla f(\vx^{(1)})=(-10,-\tfrac{21}2)$. The LP
$\min\,-10y_1-\tfrac{21}2y_2$ has vertex values $31.5$, $105$ and $-74.375$, so
again $\vy^{(1)}=(\tfrac{35}4,-\tfrac54)=\vx^{(1)}$, $\vd^{(1)}=\vzero$ and
$\Gamma(\vx^{(1)})=0$. The method stops: $\vx^*=(\tfrac{35}4,-\tfrac54)$ with
$f^*=33.8125$ ({numref}`fig-08-09-fw-graphical`b).

*Check by KKT.* Both constraints $g_1=x_1-x_2-10\le0$ and
$g_2=3-\tfrac15x_1+x_2\le0$ are active at $\vx^*$. Stationarity
$\nabla f+\lambda_1\nabla g_1+\lambda_2\nabla g_2=\vzero$ reads
$-10+\lambda_1-\tfrac15\lambda_2=0$ and $-\tfrac{21}2-\lambda_1+\lambda_2=0$,
giving $\lambda_2=25.625$ and $\lambda_1=15.125$, both nonnegative, so the
vertex is the global minimizer of this convex problem. The unconstrained
minimizer $(10,4)$ violates $\tfrac15x_1-x_2\ge3$, which is why the solution
lies on the boundary.

*(A printed solution writes the second linear objective as
$-80y_1-\tfrac{21}2y_2$ and evaluates the stopping test with the old gradient
$(-80,-14)$; the gradient at $\vx^{(1)}$ is $(-10,-\tfrac{21}2)$. Both give the
same vertex and a zero product with $\vd^{(1)}=\vzero$, so the conclusion is
unchanged. The line-search interval must be $0\le\alpha\le1$, and the
minimizer $1.171$ of $\varphi$ has to be clipped to $1$.)*

```{figure} figures/08-09-fw-graphical.svg
:name: fig-08-09-fw-graphical
:width: 680px
:class: dark-light
:alt: Two panels with the feasible triangle with vertices (0,-3), (0,-10) and (35/4,-5/4); the left adds level lines of the linear subproblem, the right adds contours of the objective and the single Frank–Wolfe step along the upper edge

{prf:ref}`ex-08-09-fw-graphical`. (a) Feasible region bounded by $x_1=0$,
$x_1-x_2=10$ and $\tfrac15x_1-x_2=3$, with level lines of the first linear
subproblem $-80y_1-14y_2$ (dashed); the optimal vertex is
$(\tfrac{35}4,-\tfrac54)$. (b) Contours of $f$ and the single Frank–Wolfe step
from $(0,-3)$ to the optimum $(\tfrac{35}4,-\tfrac54)$; the unconstrained
minimizer $(10,4)$ is outside the region.
```
:::

:::{prf:example} Frank–Wolfe for a minimum convex-cost flow
:label: ex-08-09-fw-network-flow

A flow of $210.5$ units per hour must be sent from node 1 to node 4 of a
network with the arcs $(1,2)$,
$(1,3)$, $(2,3)$, $(3,2)$, $(2,4)$ and $(3,4)$. The cost of a flow $x_a$ on arc
$a$ is convex and increasing in the flow (congestion):
$f_a(x_a)=A_ax_a+\tfrac12B_ax_a^2$ dollars per hour, with

| arc $a$ | $(1,2)$ | $(1,3)$ | $(2,3)$ | $(3,2)$ | $(2,4)$ | $(3,4)$ |
|---|---|---|---|---|---|---|
| $A_a$ (\$ per unit) | 4 | 6 | 1 | 1 | 7 | 3 |
| $B_a$ (\$ per unit$^2$) | 0.02 | 0.01 | 0.04 | 0.04 | 0.01 | 0.03 |

The problem is

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)=\sum_a f_a(x_a)\\
\st\ & x_{12}+x_{13}=210.5,\\
& x_{12}+x_{32}=x_{23}+x_{24},\\
& x_{13}+x_{23}=x_{32}+x_{34},\\
& x_{24}+x_{34}=210.5,\\
& \vx\ge\vzero .
\end{aligned}
$$

**(a)** Find an initial feasible flow. **(b)** Show that the linear
subproblem of the Frank–Wolfe method is a shortest-path problem, and carry
out the method. **(c)** Characterize the optimal flow.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

```{figure} figures/08-09-flow-network.svg
:name: fig-08-09-flow-network
:width: 460px
:class: dark-light
:alt: Directed network with four nodes; node 1 on the left with supply 210.5, node 4 on the right with demand 210.5, nodes 2 and 3 above and below joined by two opposite arcs; each arc labeled with its marginal cost function

Network of {prf:ref}`ex-08-09-fw-network-flow`. Each arc is labeled with its
marginal cost $f_a'(x_a)=A_a+B_ax_a$ (\$ per unit). Thick arcs carry flow at
the optimum.
```

**(a)** The network is drawn in {numref}`fig-08-09-flow-network`. Any path
from node 1 to node 4 carrying all $210.5$ units is feasible. A good choice is the path that is cheapest at zero flow, with
marginal costs $f_a'(0)=A_a$: the four paths cost
$1\text{–}2\text{–}4$: $4+7=11$, $1\text{–}3\text{–}4$: $6+3=9$,
$1\text{–}2\text{–}3\text{–}4$: $4+1+3=8$ and
$1\text{–}3\text{–}2\text{–}4$: $6+1+7=14$ (dollars per unit), so all flow is
loaded on $1\text{–}2\text{–}3\text{–}4$:
$\vx^{(0)}=(x_{12},x_{13},x_{23},x_{32},x_{24},x_{34})=(210.5,0,210.5,0,0,210.5)$,
with $f(\vx^{(0)})=\$3677.96$ per hour. *(The path must end at node 4; the
network has no node 5.)*
*(An earlier version of this example leaves the arc cost functions
unspecified; the quadratic costs in the table are an assumption made here so
that the method can be carried out with numbers.)*

**(b)** The objective is separable, so
$\nabla f(\vx^{(k)})=\bigl(f_a'(x_a^{(k)})\bigr)_a$ with
$f_a'(x_a)=A_a+B_ax_a$, the **marginal cost** of arc $a$. The linear
subproblem {eq}`eq-06-06-fw-lp` is

$$
\begin{aligned}
\min_{\vy}\ & \sum_a f_a'(x_a^{(k)})\,y_a\\
\st\ & \text{the four conservation constraints of the problem},\\
& \vy\ge\vzero ,
\end{aligned}
$$

a minimum-cost flow problem without capacities, with a single origin and
destination. Its optimal solution sends all $210.5$ units along a
**shortest path** from node 1 to node 4 with arc lengths $f_a'(x_a^{(k)})$
({prf:ref}`alg-04-02-dijkstra`): any split over several paths could not be
cheaper than the cheapest one. So each iteration (i) computes the marginal
costs, (ii) finds a shortest path, (iii) loads all flow on it to get
$\vy^{(k)}$ ("all-or-nothing" assignment), (iv) does an exact line search on
the segment (here a quadratic, so {eq}`eq-06-06-fw-quadratic-step` applies with
$\mQ=\diag(B_a)$), and (v) updates the flows.

*Iteration 0.* At $\vx^{(0)}$ the marginal costs are $8.21$, $6$, $9.42$, $1$,
$7$ and $9.315$, so the path costs are $15.21$ ($1\text{–}2\text{–}4$),
$15.315$ ($1\text{–}3\text{–}4$), $26.945$ ($1\text{–}2\text{–}3\text{–}4$) and
$14.00$ ($1\text{–}3\text{–}2\text{–}4$). The loaded path is now
$1\text{–}3\text{–}2\text{–}4$: $\vy^{(0)}=(0,210.5,0,210.5,210.5,0)$. The gap
is $\Gamma=2724.92$, the step $\alpha_0=0.41$, and
$\vx^{(1)}=(124.2,\ 86.3,\ 124.2,\ 86.3,\ 86.3,\ 124.2)$ with
$f(\vx^{(1)})=\$3119.38$. The lower bound $f(\vx^{(0)})-\Gamma=\$953.04$ is
valid already.

{numref}`tab-08-09-fw-network` lists the further iterations (best lower bound
$=\max_{i\le k}\bigl(f(\vx^{(i)})-\Gamma(\vx^{(i)})\bigr)$).

```{list-table} Frank–Wolfe on the convex-cost flow problem (costs in \$ per hour)
:name: tab-08-09-fw-network
:header-rows: 1

* - $k$
  - $f(\vx^{(k)})$
  - shortest path
  - $\Gamma(\vx^{(k)})$
  - $\alpha_k$
  - best lower bound
* - 0
  - 3677.96
  - 1–3–2–4
  - 2724.92
  - 0.4100
  - 953.04
* - 1
  - 3119.38
  - 1–3–4
  - 1176.48
  - 0.7021
  - 1942.90
* - 2
  - 2706.37
  - 1–2–4
  - 901.33
  - 0.3835
  - 1942.90
* - 3
  - 2533.53
  - 1–3–4
  - 115.59
  - 0.1593
  - 2417.94
* - 4
  - 2524.33
  - 1–2–4
  - 118.13
  - 0.1029
  - 2417.94
* - 5
  - 2518.25
  - 1–3–4
  - 74.15
  - 0.1101
  - 2444.10
* - 10
  - 2505.15
  - 1–2–4
  - 53.25
  - 0.05040
  - 2461.98
* - 20
  - 2497.30
  - 1–2–4
  - 29.22
  - 0.02853
  - 2474.25
* - 25
  - 2495.70
  - 1–2–4
  - 26.14
  - 0.02560
  - 2475.73
```

After 25 iterations the cost is known to lie between \$2475.73 and
\$2495.70 per hour. The early iterations make large gains; later the method
alternates between the paths $1\text{–}2\text{–}4$ and $1\text{–}3\text{–}4$
with shrinking steps, the typical sublinear behavior
({prf:ref}`thm-06-06-fw-rate`).

**(c)** For this convex problem a feasible flow is optimal if and only if
every path that carries flow has the **minimum marginal path cost**: otherwise
shifting a little flow from a dearer used path to a cheaper one would lower
the cost (these are the KKT conditions in path-flow form). Solving these
conditions (the quadratic program in the four path flows) gives the path
flows $89.0625$ on $1\text{–}2\text{–}4$, $115.25$ on $1\text{–}3\text{–}4$,
$6.1875$ on $1\text{–}2\text{–}3\text{–}4$ and $0$ on $1\text{–}3\text{–}2\text{–}4$,
that is the arc flows

$$
(x_{12},x_{13},x_{23},x_{32},x_{24},x_{34})=(95.25,\ 115.25,\ 6.1875,\ 0,\ 89.0625,\ 121.4375),
$$

with $f^*=\$2485.21$ per hour. *Check:* the three used paths have the same
marginal cost, $\$13.80$ per unit (for example
$1\text{–}2\text{–}4$: $4+0.02\cdot95.25+7+0.01\cdot89.0625=13.80$), and the
unused path $1\text{–}3\text{–}2\text{–}4$ costs $\$16.04$ per unit; flow is
conserved at nodes 2 ($95.25+0=6.1875+89.0625$) and 3
($115.25+6.1875=0+121.4375$). The Frank–Wolfe iterate at $k=25$ is within
\$10.49 per hour ($0.4\,\%$) of this optimum. The same scheme, with link
travel-time functions instead of costs, is the standard algorithm for
traffic assignment.
:::

(sec-08-09-penalty)=
## Penalty methods

:::{prf:example} KKT conditions and the penalty method for an equality-constrained quadratic
:label: ex-08-09-penalty-equality-qp

Solve

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)=x_1^2+x_1x_2+x_2^2-2x_2\\
\st\ & x_1+x_2=2
\end{aligned}
$$

analytically **(a)** with the KKT conditions and **(b)** with the quadratic
penalty function $\rho\,(x_1+x_2-2)^2$, letting $\rho\to\infty$.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

**(a)** With $h(\vx)=x_1+x_2-2$ and the Lagrangian
$L=f+\mu h$ ({prf:ref}`thm-06-03-kkt`), the KKT conditions are

$$
\begin{aligned}
2x_1+x_2+\mu&=0,\\
x_1+2x_2-2+\mu&=0,\\
x_1+x_2&=2 .
\end{aligned}
$$

Subtracting the first two gives $x_1-x_2+2=0$; with $x_1+x_2=2$ this yields
$x_1^*=0$, $x_2^*=2$ and $\mu^*=-2$, with $f^*=0$. The Hessian
$\begin{pmatrix}2&1\\1&2\end{pmatrix}$ is positive definite, so the problem is
convex and $(0,2)$ is the global minimizer.

**(b)** The penalized objective is
$\Phi_\rho(\vx)=x_1^2+x_1x_2+x_2^2-2x_2+\rho(x_1+x_2-2)^2$. Setting its
gradient to zero,

$$
\begin{aligned}
2x_1+x_2+2\rho(x_1+x_2-2)&=0,\\
x_1+2x_2-2+2\rho(x_1+x_2-2)&=0 .
\end{aligned}
$$

Subtracting gives $x_1(\rho)=x_2(\rho)-2$ for every $\rho$; substituting into
the first equation, $3x_2-4+4\rho x_2-8\rho=0$, so

$$
\begin{aligned}
x_2(\rho)&=\frac{4+8\rho}{3+4\rho},\\
x_1(\rho)&=-\frac{2}{3+4\rho},\\
h(\vx(\rho))&=-\frac{4}{3+4\rho}.
\end{aligned}
$$

As $\rho\to\infty$, $\vx(\rho)\to(0,2)$, the KKT solution of (a). The
multiplier estimate {eq}`eq-06-06-multiplier-estimates` is
$2\rho\,h(\vx(\rho))=-8\rho/(3+4\rho)\to-2=\mu^*$.

```{list-table} Quadratic penalty solutions for $x_1+x_2=2$
:header-rows: 1

* - $\rho$
  - $\vx(\rho)$
  - $h(\vx(\rho))$
  - $2\rho\,h(\vx(\rho))$
  - condition number of $\nabla^2\Phi_\rho$
* - 1
  - $(-0.2857,\ 1.714)$
  - $-0.5714$
  - $-1.143$
  - 7
* - 10
  - $(-0.04651,\ 1.953)$
  - $-0.09302$
  - $-1.860$
  - 43
* - 100
  - $(-0.004963,\ 1.995)$
  - $-0.009926$
  - $-1.985$
  - 403
* - 1000
  - $(-0.0004996,\ 1.9995)$
  - $-0.0009993$
  - $-1.999$
  - 4003
```

The Hessian $\nabla^2\Phi_\rho=\begin{pmatrix}2+2\rho&1+2\rho\\1+2\rho&2+2\rho\end{pmatrix}$
has eigenvalues $1$ (direction $(1,-1)$, along the constraint) and $3+4\rho$
(direction $(1,1)$, across it), so its condition number $3+4\rho$ grows
linearly with $\rho$ (see {ref}`sec-06-06-ill-conditioning`). *(A printed
solution gives $x_2(\rho)=(4+4\rho)/(3+2\rho)$; that formula belongs to the
penalty $\tfrac\rho2(x_1+x_2-2)^2$. With the penalty $\rho(x_1+x_2-2)^2$ as
stated, $x_2(\rho)=(4+8\rho)/(3+4\rho)$. Both tend to $2$.)*
:::

:::{prf:example} Quadratic penalty method with a nonlinear equality constraint
:label: ex-08-09-penalty-parabola

Apply the quadratic penalty method to

$$
\begin{aligned}
\min_{\vx}\ & f(\vx)=(x_1-2)^4+(x_1-2x_2)^2\\
\st\ & h(\vx)=x_1^2-x_2=0,
\end{aligned}
$$

with $\psi(\vx)=(x_1^2-x_2)^2$, $\rho_0=0.1$ and $\rho_{k+1}=10\rho_k$. Report
$\vx^{(k)}=\vx(\rho_k)$, $f(\vx^{(k)})$ and $\Phi_{\rho_k}(\vx^{(k)})$, and decide
whether the method converges.
:::

:::{dropdown} Solution
:class-container: solution
:animate: fade-in

Each subproblem $\min\Phi_{\rho_k}(\vx)=f(\vx)+\rho_k(x_1^2-x_2)^2$ is an
unconstrained problem in two variables, solved here by a quasi-Newton method
({ref}`sec-06-05-quasi-newton`) started from the previous solution (and,
for $k=0$, from the unconstrained minimizer $(2,1)$ of $f$).

```{list-table} Quadratic penalty method for $\min(x_1-2)^4+(x_1-2x_2)^2$ s.t. $x_1^2=x_2$
:name: tab-08-09-penalty-parabola
:header-rows: 1

* - $k$
  - $\rho_k$
  - $\vx^{(k)}$
  - $f(\vx^{(k)})$
  - $h(\vx^{(k)})$
  - $\rho_k\psi(\vx^{(k)})$
  - $\Phi_{\rho_k}(\vx^{(k)})$
  - $2\rho_kh(\vx^{(k)})$
* - 0
  - 0.1
  - $(1.4539,\ 0.7608)$
  - 0.09353
  - 1.353
  - 0.1831
  - 0.2766
  - 0.2706
* - 1
  - 1
  - $(1.1687,\ 0.7407)$
  - 0.5752
  - 0.6252
  - 0.3909
  - 0.9662
  - 1.250
* - 2
  - 10
  - $(0.9906,\ 0.8425)$
  - 1.520
  - 0.1389
  - 0.1928
  - 1.713
  - 2.777
* - 3
  - 100
  - $(0.9508,\ 0.8875)$
  - 1.891
  - 0.01648
  - 0.02717
  - 1.918
  - 3.297
* - 4
  - 1000
  - $(0.9461,\ 0.8934)$
  - 1.941
  - $1.682\times10^{-3}$
  - $2.828\times10^{-3}$
  - 1.943
  - 3.363
* - 5
  - 10000
  - $(0.9456,\ 0.8941)$
  - 1.946
  - $1.685\times10^{-4}$
  - $2.839\times10^{-4}$
  - 1.946
  - 3.370
```

**The method converges.** The iterates approach the parabola $x_2=x_1^2$
from outside ($h>0$) and settle at $(0.9456,\ 0.8941)$
({numref}`fig-08-09-penalty-path`). Solving the constrained problem directly
(KKT conditions, by a sequential quadratic programming code) gives
$\vx^*=(0.945583,\ 0.894127)$, $f^*=1.94618$ and the multiplier
$\mu^*=3.37069$ of $L=f+\mu h$. As
{prf:ref}`thm-06-06-penalty-convergence` predicts, the penalized values
$\Phi_{\rho_k}(\vx^{(k)})$ increase towards $f^*$ ($0.2766$, $0.9662$, $1.713$, $1.918$, $1.943$, $1.946$),
the violations $h(\vx^{(k)})$ decrease by about a factor $10$ per round, and
the multiplier estimates $2\rho_kh(\vx^{(k)})$ approach $\mu^*$
({ref}`sec-06-06-multipliers`). The quantity $\rho_k\psi(\vx^{(k)})$ first rises
($0.183\to0.391$) and then falls ($0.193$, $0.027$, $0.0028$, $0.00028$), so with
the stopping test $\rho_k\psi(\vx^{(k)})\le\varepsilon=0.01$ of
{prf:ref}`alg-06-06-penalty` the method stops at $k=4$; after only four rounds
($k\le3$) it has simply not yet met the test.

The price is conditioning: the condition number of $\nabla^2\Phi_{\rho_k}$ at
$\vx^{(k)}$ is $3.7$, $4.8$, $16.7$, $123$, $1183$, $1.18\times10^4$ for
$k=0,\dots,5$, growing like $\rho_k$.

*(An earlier version of this table prints the objective as
$(x_1-2)^2+(x_1-2x_2)^2$; its numbers belong to the exponent $4$ used here
(with exponent $2$ the first iterate would be $(1.5687,\ 0.8252)$). It also
gives $f(\vx^{(1)})=0.5793$, $1.5203$ and $1.8917$ where the correct values
are $0.5752$, $1.5201$ and $1.8912$; labels its last column with the previous
iterate although the values are $\Phi_{\rho_k}(\vx^{(k)})$ of the same row;
and concludes that "the method does not converge". As shown above, the
method converges; it had only not reached its stopping test.)*

```{figure} figures/08-09-penalty-path.svg
:name: fig-08-09-penalty-path
:width: 520px
:class: dark-light
:alt: Contours of the objective around the point (2,1), the parabola x2 equals x1 squared, and a sequence of penalty minimizers moving from near (1.45, 0.76) towards the constrained optimum on the parabola

{prf:ref}`ex-08-09-penalty-parabola`. Contours of $f(\vx)=(x_1-2)^4+(x_1-2x_2)^2$
(unconstrained minimizer $(2,1)$), the constraint $x_2=x_1^2$ (blue) and the
penalty minimizers $\vx(\rho)$ for $\rho=0.1,1,10,100,1000$ (dots; the last is unlabeled, next
to $\vx^*$), which approach
the constrained minimizer $\vx^*=(0.9456,\ 0.8941)$ from the infeasible side.
```
:::

## References

- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 8:
  line searches, steepest descent and Newton's method, including the example
  $(x_1-2)^4+(x_1-2x_2)^2$; Chapter 9: penalty methods.)
- Courant, R. (1943). Variational methods for the solution of problems of
  equilibrium and vibrations. *Bulletin of the American Mathematical Society,
  49*(1), 1–23. <https://doi.org/10.1090/S0002-9904-1943-07818-4>
- Frank, M., & Wolfe, P. (1956). An algorithm for quadratic programming.
  *Naval Research Logistics Quarterly, 3*(1–2), 95–110.
  <https://doi.org/10.1002/nav.3800030109>
- Kiefer, J. (1953). Sequential minimax search for a maximum. *Proceedings of
  the American Mathematical Society, 4*(3), 502–506.
  <https://doi.org/10.1090/S0002-9939-1953-0055639-3>
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming* (4th
  ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3> (Chapters 8 and
  13: descent methods; penalty methods.)
- Sheffi, Y. (1985). *Urban transportation networks: Equilibrium analysis with
  mathematical programming methods*. Prentice-Hall. (Chapters 4–5: the
  convex-combination algorithm for network flow and traffic assignment
  problems.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 21: nonlinear programming algorithms.)
