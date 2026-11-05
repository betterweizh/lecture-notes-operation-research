# Nonlinear Programming in Python

## Overview

This chapter implements the nonlinear-programming methods of
{doc}`06-01-nlp-models` to {doc}`06-06-constrained-algorithms` and uses the
solvers of SciPy and CVXPY on the book's examples. It begins with gradients
and Hessians ({doc}`01-02-matrix-calculus`), written analytically and
checked by finite differences. It
then codes the one-dimensional line searches of {doc}`06-04-line-search`
(bisection, golden-section search, Newton's method and Armijo backtracking), the
steepest-descent and Newton methods of {doc}`06-05-unconstrained-algorithms`,
and the Frank–Wolfe and quadratic-penalty methods of
{doc}`06-06-constrained-algorithms`, with `linprog` solving the linear
subproblems. A small **KKT checker** measures how far a point and a set of
multipliers are from satisfying the Karush–Kuhn–Tucker (KKT) conditions of
{doc}`06-03-kkt-conditions`: stationarity, feasibility, complementary
slackness and the sign conditions; it is
applied to the KKT example of that chapter, whose 32 active sets are also
enumerated. The general-purpose solvers of `scipy.optimize.minimize` (BFGS,
SLSQP, trust-constr) are run on the same examples and on the least-squares
estimation of {doc}`06-01-nlp-models`. The last section solves convex problems
with CVXPY, compares its dual values with the KKT multipliers, and
reproduces the Lagrange-dual examples of
{doc}`02-03-convex-optimization-duality`; it requires CVXPY (any recent
version with its default solvers). The code is not run when the book is
built: copy the blocks, in order, into a script or a notebook and compare
the printed results with the `# Expected:` comments.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- write gradients and Hessians of the book's functions and check them by
  central finite differences;
- implement bisection, golden-section search, Newton's method and
  backtracking, and the steepest-descent, Newton, Frank–Wolfe and
  quadratic-penalty methods;
- check the KKT conditions numerically and enumerate active sets;
- call `scipy.optimize.minimize` with bounds, linear and nonlinear
  constraints, read its multipliers and use multiple starts for nonconvex
  problems; and
- model a convex problem in CVXPY and interpret its dual values with the
  sign conventions of the book.
:::

(sec-09-06-setup)=
## Setup

```python
import itertools
import math

import cvxpy as cp
import numpy as np
from scipy.optimize import (
    LinearConstraint,
    linprog,
    minimize,
    minimize_scalar,
)

np.set_printoptions(precision=3, suppress=True, linewidth=120)
TAU = (math.sqrt(5) - 1) / 2
```

`TAU` is the golden-section ratio $\tau=(\sqrt5-1)/2\approx0.618$. Where
the iteration tables of {doc}`06-04-line-search` to
{doc}`06-06-constrained-algorithms` give four significant digits, the values are
rounded explicitly with `round` or `np.round`.

(sec-09-06-derivatives)=
## Gradients and Hessians

All methods below need $\nabla f$ and some need $\nabla^2f$. Analytic
derivatives are fast and exact, but easy to get wrong; a finite-difference
check catches the mistakes. The **central difference**

$$
\frac{\partial f}{\partial x_j}(\vx)\approx\frac{f(\vx+h\ve_j)-f(\vx-h\ve_j)}{2h}
$$

has an error of order $h^2$, so with $h\approx10^{-6}$ (relative to the size
of $\vx$) it agrees with a correct gradient to about eight digits. Applied to
the analytic gradient it gives the Hessian, which is then symmetrized.

```python
def fd_gradient(f, x, h: float = 1e-6) -> np.ndarray:
    """Central-difference approximation of the gradient of f at x."""
    x = np.asarray(x, dtype=float)
    if h <= 0:
        raise ValueError("step h must be positive")
    g = np.zeros_like(x)
    for j in range(x.size):
        e = np.zeros_like(x)
        e[j] = h * max(1.0, abs(x[j]))
        g[j] = (f(x + e) - f(x - e)) / (2 * e[j])
    return g


def fd_hessian(grad, x, h: float = 1e-5) -> np.ndarray:
    """Symmetrized central-difference Hessian from an analytic gradient."""
    x = np.asarray(x, dtype=float)
    H = np.column_stack(
        [fd_gradient(lambda z: grad(z)[i], x, h) for i in range(x.size)]
    ).T
    return (H + H.T) / 2
```

The quadratic of {prf:ref}`ex-06-05-steepest-descent`,
$f(\vx)=(x_1-2)^2+10(x_2-2)^2$, has $\nabla f(\vx)=(2x_1-4,\ 20x_2-40)$ and
the constant Hessian $\diag(2,20)$; the objective of
{prf:ref}`ex-06-06-frank-wolfe`, $f(\vx)=-5x_1+x_1^2-8x_2+2x_2^2$, has
$\nabla f=(2x_1-5,\ 4x_2-8)$ and $\nabla^2f=\diag(2,4)$.

```python
def f_quad(x):
    """Quadratic of ex-06-05-steepest-descent."""
    return (x[0] - 2) ** 2 + 10 * (x[1] - 2) ** 2


def grad_quad(x):
    """Gradient of f_quad."""
    return np.array([2 * x[0] - 4, 20 * x[1] - 40])


def hess_quad(x):
    """Hessian of f_quad (constant)."""
    return np.diag([2.0, 20.0])


def f_fw(x):
    """Objective of ex-06-06-frank-wolfe."""
    return -5 * x[0] + x[0] ** 2 - 8 * x[1] + 2 * x[1] ** 2


def grad_fw(x):
    """Gradient of f_fw."""
    return np.array([2 * x[0] - 5, 4 * x[1] - 8])


x_test = np.array([-4.0, -3.0])
print(grad_quad(x_test), fd_gradient(f_quad, x_test))
# Expected: [ -12. -100.] [ -12. -100.]
print(fd_hessian(grad_quad, x_test))
# Expected:
#   [[ 2.  0.]
#    [ 0. 20.]]
print(np.abs(grad_fw(x_test) - fd_gradient(f_fw, x_test)).max() < 1e-6)
# Expected: True
```

A harder case is the least-squares objective {eq}`eq-06-01-sse` of
{prf:ref}`ex-06-01-consumer-regression`, a polynomial of degree four in the six
parameters $(\vbeta,\vq)$. With the residuals $r_{ti}$ of
{eq}`eq-06-01-residual` and $m_t=E_t-\sum_jp_{tj}q_j$ (the uncommitted money of
month $t$),

$$
\begin{aligned}
\frac{\partial r_{ti}}{\partial\beta_k}&=-[i=k]\,m_t,\\
\frac{\partial r_{ti}}{\partial q_j}&=-[i=j]\,p_{ti}+\beta_ip_{tj},
\end{aligned}
$$

and $\nabla f=2\sum_{t,i}r_{ti}\nabla r_{ti}$. In NumPy the sums over $t$
and $i$ become matrix products.

```python
# Monthly budget E_t ($), prices p_ti ($/unit), quantities Q_ti (units)
E = np.array([1500.0, 2000.0, 1500.0, 1200.0])
prices = np.array(
    [
        [2.5, 4.0, 5.0],
        [2.7, 3.5, 4.9],
        [2.8, 4.0, 5.6],
        [2.6, 3.6, 4.5],
    ]
)
Q = np.array(
    [
        [198, 124, 100],
        [200, 200, 110],
        [100, 80, 60],
        [152, 115, 65],
    ]
)


def residuals(v: np.ndarray) -> np.ndarray:
    """Residuals r_ti for the parameter vector v = (beta_1..3, q_1..3)."""
    beta, q = v[:3], v[3:]
    m = E - prices @ q
    return prices * Q - (prices * q + np.outer(m, beta))


def sse(v: np.ndarray) -> float:
    """Sum of squared residuals f(beta, q)."""
    return float(np.sum(residuals(v) ** 2))


def grad_sse(v: np.ndarray) -> np.ndarray:
    """Analytic gradient of sse."""
    beta, q = v[:3], v[3:]
    r = residuals(v)
    m = E - prices @ q
    g_beta = -2 * m @ r
    g_q = -2 * np.sum(r * prices, axis=0) + 2 * prices.T @ (r @ beta)
    return np.concatenate([g_beta, g_q])


v_test = np.array([0.3, 0.4, 0.3, 50.0, 30.0, 20.0])
print(round(sse(v_test), 2))
# Expected: 148650.23
g_exact, g_fd = grad_sse(v_test), fd_gradient(sse, v_test)
print(np.abs(g_exact - g_fd).max() / np.abs(g_exact).max() < 1e-7)
# Expected: True
```

The relative difference is below $10^{-7}$, so the analytic gradient can be
trusted; it is used by SLSQP in {ref}`sec-09-06-scipy`.

(sec-09-06-line-search)=
## One-dimensional search

All four methods are tried on the problem of {doc}`06-04-line-search`:
minimize $f(x)=\sin x$ on $[3,6]$, whose minimizer is
$x^*=3\pi/2\approx4.712389$.

### Bisection

{prf:ref}`alg-06-04-bisection` halves the interval according to the sign of
$f'$ at the midpoint. The function records the rows of
{numref}`tab-06-04-bisection`.

```python
def bisection(fprime, a: float, b: float, eps: float):
    """Bisection search for a root of f' in [a, b] with f'(a) < 0 < f'(b).

    Returns the estimate (a_k + b_k)/2 and the rows (a_k, b_k, x_k, f'(x_k)).
    """
    if not (a < b and fprime(a) < 0 < fprime(b)):
        raise ValueError("need a < b and f'(a) < 0 < f'(b)")
    if eps <= 0:
        raise ValueError("eps must be positive")
    rows = []
    while b - a > eps:
        x = (a + b) / 2
        slope = fprime(x)
        rows.append((a, b, x, slope))
        if slope == 0:
            return x, rows
        if slope > 0:
            b = x
        else:
            a = x
    return (a + b) / 2, rows


x_hat, rows = bisection(math.cos, 3.0, 6.0, eps=0.05)
for a, b, x, slope in rows:
    print(f"{a:.6f} {b:.6f} {x:.6f} {slope:+.4f}")
# Expected:
#   3.000000 6.000000 4.500000 -0.2108
#   4.500000 6.000000 5.250000 +0.5121
#   4.500000 5.250000 4.875000 +0.1619
#   4.500000 4.875000 4.687500 -0.0249
#   4.687500 4.875000 4.781250 +0.0688
#   4.687500 4.781250 4.734375 +0.0220
print(round(x_hat, 7), len(rows), round(abs(x_hat - 1.5 * math.pi), 6))
# Expected: 4.7109375 6 0.001451
```

Six derivative evaluations give $\hat x=4.7109375$, with an error of
$0.001451$; the number of halvings is
$\lceil\log_2(3/0.05)\rceil=6$, as {eq}`eq-06-04-bisection-count` predicts.

### Golden-section search

{prf:ref}`alg-06-04-golden-section` keeps two interior points at the
fractions $1-\tau$ and $\tau$ of the interval and reuses one of them in
every step, so each reduction costs one new function value.

```python
def golden_section(f, a: float, b: float, eps: float):
    """Golden-section search for a minimizer of a unimodal f on [a, b].

    Returns the estimate, the final interval and the number of f values.
    """
    if not a < b or eps <= 0:
        raise ValueError("need a < b and eps > 0")
    p, q = a + (1 - TAU) * (b - a), a + TAU * (b - a)
    fp, fq = f(p), f(q)
    evaluations = 2
    while b - a > eps:
        keep_right = fp > fq
        if keep_right:
            a, p, fp = p, q, fq
            q = a + TAU * (b - a)
        else:
            b, q, fq = q, p, fp
            p = a + (1 - TAU) * (b - a)
        # the new trial point is evaluated only if the search goes on
        if b - a > eps:
            if keep_right:
                fq = f(q)
            else:
                fp = f(p)
            evaluations += 1
    return (a + b) / 2, (a, b), evaluations


x_gs, interval, n_f = golden_section(math.sin, 3.0, 6.0, eps=0.05)
print(round(x_gs, 4), np.round(interval, 4).tolist(), n_f)
# Expected: 4.7067 [4.6869, 4.7264] 10
```

The loop stops after nine reductions with the interval
$[4.6869,4.7264]$ and the estimate $4.7067$ of
{prf:ref}`ex-06-04-golden-sine`, after $N=10$ function values. The new trial
point of the last reduction is not evaluated, because the stopping test
ends the search before it would be used.

### Newton's method

{prf:ref}`alg-06-04-newton-1d` iterates $x_{k+1}=x_k-f'(x_k)/f''(x_k)$. The
function deliberately has no safeguard, so that the failure of
{prf:ref}`ex-06-04-newton-sine` (b) can be seen; it refuses only a zero
curvature.

```python
def newton_1d(
    fprime, fsecond, x0: float, eps: float = 1e-6, max_iter: int = 50
):
    """Newton's method for f'(x) = 0; returns the list of iterates."""
    xs = [x0]
    while abs(fprime(xs[-1])) > eps:
        if len(xs) > max_iter:
            raise ValueError("no convergence within max_iter iterations")
        curvature = fsecond(xs[-1])
        if curvature == 0:
            raise ValueError("zero curvature: Newton step undefined")
        xs.append(xs[-1] - fprime(xs[-1]) / curvature)
    return xs


def neg_sin(x):
    """Second derivative of sin."""
    return -math.sin(x)


print(np.round(newton_1d(math.cos, neg_sin, 4.0), 6).tolist())
# Expected: [4.0, 4.863691, 4.711224, 4.712389]
xs = newton_1d(math.cos, neg_sin, 3.0)
print(round(xs[1], 4), round(xs[-1], 6), round(math.sin(xs[-1]), 6))
# Expected: -4.0153 -4.712389 1.0
```

From $x_0=4$ three steps reach $3\pi/2$; from $x_0=3$ the first step jumps
to $-4.0153$ and the iteration converges to $-3\pi/2$, a maximizer of
$\sin x$ (value 1) outside $[3,6]$, while the stopping test reports
success.

### Backtracking (Armijo)

{prf:ref}`alg-06-04-backtracking` shrinks a trial step until the Armijo
condition {eq}`eq-06-04-armijo` holds. It works in $\R^n$; the illustration of
{ref}`sec-06-04-inexact` is the one-dimensional case $x=3$, $d=-f'(3)$, with
$c_1=0.25$, $\bar\alpha=4$ and $\gamma=\tfrac12$.

```python
def backtracking(
    f, grad, x, d, alpha0=1.0, c1=1e-4, gamma=0.5, max_halvings=60
):
    """Armijo backtracking: largest alpha0 * gamma^m with sufficient decrease.

    Returns the accepted step and the list of rejected trial steps.
    """
    x, d = np.atleast_1d(x).astype(float), np.atleast_1d(d).astype(float)
    slope = float(np.dot(grad(x), d))
    if slope >= 0:
        raise ValueError("d is not a descent direction")
    alpha, rejected = alpha0, []
    for _ in range(max_halvings):
        if f(x + alpha * d) <= f(x) + c1 * alpha * slope:
            return alpha, rejected
        rejected.append(alpha)
        alpha *= gamma
    raise ValueError("no acceptable step found")


def sin_vec(x):
    """sin of a one-element vector, for backtracking in R^n."""
    return math.sin(x[0])


def cos_vec(x):
    """Gradient of sin_vec."""
    return np.array([math.cos(x[0])])


d0 = -math.cos(3.0)
alpha, rejected = backtracking(sin_vec, cos_vec, 3.0, d0, alpha0=4.0, c1=0.25)
print(rejected, alpha, round(3.0 + alpha * d0, 4))
# Expected: [4.0] 2.0 4.98
```

The step $\alpha=4$ is rejected and $\alpha=2$ accepted, giving $x=4.98$.
For comparison, SciPy's bounded scalar minimizer (Brent's method,
parabolic interpolation safeguarded by golden-section steps; Brent, 1973) finds $3\pi/2$ to six digits.

```python
res = minimize_scalar(math.sin, bounds=(3, 6), method="bounded")
print(round(res.x, 6), res.nfev)
# Expected: 4.71239 8
```

(sec-09-06-unconstrained)=
## Steepest descent and Newton's method

{prf:ref}`alg-06-05-steepest-descent` takes the direction $-\nabla f$ and a
step from a line search. The line search is passed as a function
`step(x, d)`, so that the exact step {eq}`eq-06-05-sd-step` of a quadratic
and Armijo backtracking can be exchanged.

```python
def steepest_descent(
    grad, x0, step, tol: float = 1e-6, max_iter: int = 10_000
):
    """Steepest descent; returns the list of iterates.

    Parameters
    ----------
    grad : callable
        Gradient of f.
    x0 : array_like
        Starting point.
    step : callable
        ``step(x, d) -> alpha``, the line search.
    tol : float
        Stop when the gradient norm is at most tol.
    """
    xs = [np.asarray(x0, dtype=float)]
    while np.linalg.norm(grad(xs[-1])) > tol:
        if len(xs) > max_iter:
            raise ValueError("no convergence within max_iter iterations")
        d = -grad(xs[-1])
        xs.append(xs[-1] + step(xs[-1], d) * d)
    return xs


def exact_step(Q: np.ndarray):
    """Exact line search for a quadratic with Hessian Q: alpha = g'g / g'Qg."""
    if np.any(np.linalg.eigvalsh(Q) <= 0):
        raise ValueError("Q must be positive definite")
    return lambda x, d: float(d @ d) / float(d @ Q @ d)
```

On {prf:ref}`ex-06-05-steepest-descent` the first step reaches
$(-3.3922,\ 2.0647)$, every step multiplies $f-f^*$ by $0.1018$, and the
gradient test $\norm{\nabla f}\le10^{-6}$ is met at $k=17$.

```python
Q_quad = np.diag([2.0, 20.0])
xs = steepest_descent(grad_quad, [-4, -3], exact_step(Q_quad))
print(np.round(xs[1], 4).tolist(), round(f_quad(xs[1]), 2), len(xs) - 1)
# Expected: [-3.3922, 2.0647] 29.12 17
values = [f_quad(x) for x in xs]
print(np.round([values[k + 1] / values[k] for k in range(4)], 4).tolist())
# Expected: [0.1018, 0.1018, 0.1018, 0.1018]
print(round(float(grad_quad(xs[0]) @ grad_quad(xs[1])), 10))
# Expected: 0.0
xs_worst = steepest_descent(grad_quad, [-8, 1], exact_step(Q_quad))
print(round(f_quad(xs_worst[1]) / f_quad(xs_worst[0]), 4), len(xs_worst) - 1)
# Expected: 0.6694 86
```

Consecutive gradients are orthogonal ({eq}`eq-06-05-orthogonal`); from the
worst-case start $(-8,1)$ the ratio is $((\kappa-1)/(\kappa+1))^2=0.6694$ of
{prf:ref}`thm-06-05-sd-rate`, and many more iterations are needed. With
backtracking instead of the exact step the method still converges, at a
cost that depends on the constants.

```python
def armijo_step(f, grad, alpha0=1.0, c1=1e-4):
    """Line search for steepest_descent: Armijo backtracking from alpha0."""
    return lambda x, d: backtracking(f, grad, x, d, alpha0=alpha0, c1=c1)[0]


xs_bt = steepest_descent(grad_quad, [-4, -3], armijo_step(f_quad, grad_quad))
print(len(xs_bt) - 1, np.round(xs_bt[-1], 6).tolist())
# Expected: 64 [2.0, 2.0]
```

{prf:ref}`alg-06-05-newton` solves $\nabla^2f(\vx)\vd=-\nabla f(\vx)$ (never
inverting the Hessian) and takes the full step. On a quadratic one step
suffices ({prf:ref}`ex-06-05-newton-quadratic`).

```python
def newton(grad, hess, x0, tol: float = 1e-8, max_iter: int = 100):
    """Newton's method with unit steps; returns the list of iterates."""
    xs = [np.asarray(x0, dtype=float)]
    while np.linalg.norm(grad(xs[-1])) > tol:
        if len(xs) > max_iter:
            raise ValueError("no convergence within max_iter iterations")
        d = np.linalg.solve(hess(xs[-1]), -grad(xs[-1]))
        xs.append(xs[-1] + d)
    return xs


print(newton(grad_quad, hess_quad, [-4, -3]))
# Expected: [array([-4., -3.]), array([2., 2.])]
```

SciPy's `minimize` offers both kinds of method: `BFGS` builds a
quasi-Newton approximation from gradients (here it needs a few iterations
because its line search is inexact), and `Newton-CG` uses the Hessian.

```python
res_bfgs = minimize(f_quad, [-4, -3], jac=grad_quad, method="BFGS")
print(np.round(res_bfgs.x, 6).tolist(), res_bfgs.nit)
# Expected: [2.0, 2.0] 8
res_ncg = minimize(
    f_quad, [-4, -3], jac=grad_quad, hess=hess_quad, method="Newton-CG"
)
print(np.round(res_ncg.x, 6).tolist(), res_ncg.nit)
# Expected: [2.0, 2.0] 3
```

(sec-09-06-kkt)=
## Checking the KKT conditions

The KKT conditions of {prf:ref}`thm-06-03-kkt` for
$\min f(\vx)$ s.t. $\vg(\vx)\le\vzero$, $\vh(\vx)=\vzero$ are four groups of
conditions: stationarity {eq}`eq-06-03-kkt-stationarity`, primal feasibility
{eq}`eq-06-03-kkt-primal`, the sign condition {eq}`eq-06-03-kkt-dual` and
complementary slackness {eq}`eq-06-03-kkt-cs`. The checker returns one number
per group, the size of the worst violation; a KKT point gives four values
near zero. The constraint functions return vectors and their Jacobians have
one row per constraint ($\nabla g_i^\T$).

```python
def kkt_check(x, lam, grad_f, g, jac_g, mu=None, h=None, jac_h=None) -> dict:
    """Worst violation of each group of KKT conditions at (x, lam, mu).

    Problem: min f(x) s.t. g(x) <= 0, h(x) = 0.
    """
    x, lam = np.asarray(x, dtype=float), np.asarray(lam, dtype=float)
    gx = np.asarray(g(x), dtype=float)
    if lam.shape != gx.shape:
        raise ValueError("one multiplier per inequality constraint is needed")
    stationarity = grad_f(x) + np.asarray(jac_g(x)).T @ lam
    primal = max(0.0, gx.max())
    if h is not None:
        mu = np.asarray(mu, dtype=float)
        stationarity = stationarity + np.asarray(jac_h(x)).T @ mu
        primal = max(primal, np.abs(h(x)).max())
    return {
        "stationarity": float(np.abs(stationarity).max()),
        "primal": float(primal),
        "dual": float(max(0.0, -lam.min())),
        "complementarity": float(np.abs(lam * gx).max()),
    }
```

The example {prf:ref}`ex-06-03-kkt-example` minimizes
$x_1^2+x_2^2+x_3^2$ subject to five linear constraints, written here as
$\mA\vx-\vb\le\vzero$ with one row of $\mA$ per constraint $g_i$.

```python
A_kkt = np.array(
    [
        [2, 1, 0],
        [1, 0, 1],
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
    ],
    dtype=float,
)
b_kkt = np.array([5, 2, -1, -2, 0], dtype=float)


def f_kkt(x):
    """Objective x1^2 + x2^2 + x3^2 of ex-06-03-kkt-example."""
    return float(x @ x)


def grad_kkt(x):
    """Gradient of f_kkt."""
    return 2 * np.asarray(x, dtype=float)


def g_kkt(x):
    """Constraint values g(x) = A x - b (feasible when <= 0)."""
    return A_kkt @ x - b_kkt


def jac_kkt(x):
    """Jacobian of g_kkt, one row per constraint."""
    return A_kkt


print(
    list(
        kkt_check(
            [1, 2, 0], [0, 0, 2, 4, 0], grad_kkt, g_kkt, jac_kkt
        ).values()
    )
)
# Expected: [0.0, 0.0, 0.0, 0.0]
check = kkt_check([1.5, 2, 0], [-1.5, 0, 0, 2.5, 0], grad_kkt, g_kkt, jac_kkt)
print(list(check.values()))
# Expected: [0.0, 0.0, 1.5, 0.0]
```

The point $\vx^*=(1,2,0)$ with $\vlambda^*=(0,0,2,4,0)$ satisfies all
conditions exactly. The point $(1.5,2,0)$ of the active set $\{1,4\}$ in
{numref}`tab-06-03-active-set-cases` is stationary, feasible and
complementary, but $\lambda_1=-1.5<0$: it fails only the sign condition.

### Enumerating the active sets

{prf:ref}`alg-06-03-active-set` can be run over all $2^5=32$ subsets $S$
because the constraints are linear and $f$ is quadratic: with
$g_i=0$ for $i\in S$ and $\lambda_i=0$ otherwise, stationarity
$2\vx+\mA_S^\T\vlambda_S=\vzero$ and $\mA_S\vx=\vb_S$ form one linear
system in $(\vx,\vlambda_S)$. A least-squares solution with a nonzero
residual means that the system is **inconsistent**; otherwise the solution
is classified as infeasible (some $g_i>0$, $i\notin S$), as having a
negative multiplier, or as a KKT point.

```python
def classify_active_set(S: tuple, tol: float = 1e-9):
    """Solve the KKT system with active set S; return (outcome, x, lam)."""
    S = list(S)
    n, AS = 3, A_kkt[S]
    K = np.block([[2 * np.eye(n), AS.T], [AS, np.zeros((len(S), len(S)))]])
    rhs = np.concatenate([np.zeros(n), b_kkt[S]])
    sol = np.linalg.lstsq(K, rhs, rcond=None)[0]
    if np.abs(K @ sol - rhs).max() > tol:
        return "inconsistent", None, None
    x, lam = sol[:n], np.zeros(5)
    lam[S] = sol[n:]
    if g_kkt(x).max() > tol:
        return "infeasible", x, lam
    if lam.min() < -tol:
        return "negative multiplier", x, lam
    return "KKT point", x, lam


outcomes = {}
for size in range(6):
    for S in itertools.combinations(range(5), size):
        outcome, x, lam = classify_active_set(S)
        outcomes.setdefault(outcome, []).append(S)
kinds = ["infeasible", "negative multiplier", "inconsistent", "KKT point"]
print([len(outcomes[kind]) for kind in kinds])
# Expected: [14, 8, 8, 2]
for S in outcomes["KKT point"]:
    print([i + 1 for i in S], classify_active_set(S)[1:])
# Expected:
#   [3, 4] (array([1., 2., 0.]), array([0., 0., 2., 4., 0.]))
#   [3, 4, 5] (array([1., 2., 0.]), array([0., 0., 2., 4., 0.]))
```

The counts 14, 8, 8 and 2 and the two accepted sets $\{3,4\}$ and
$\{3,4,5\}$ are those of {prf:ref}`ex-06-03-kkt-example`; both give
$\vx^*=(1,2,0)$ with $\vlambda^*=(0,0,2,4,0)$.

### SciPy's constrained solvers

`minimize` with `method="SLSQP"` (sequential quadratic programming) takes
constraints as dictionaries. SciPy writes inequalities as
$c(\vx)\ge0$, the opposite of the book's $g(\vx)\le0$, so $c=-g$ is passed.
Recent SciPy versions (1.16 and later) also return the multipliers of
SLSQP in `res.multipliers`, in the book's sign convention.
`method="trust-constr"` accepts a `LinearConstraint`
$\vl\le\mA\vx\le\vu$ and returns its multipliers in `res.v` (one array per
constraint object); it is an interior-point method, so its solution
approaches the boundary to within the tolerance only.

```python
cons_kkt = [
    {"type": "ineq", "fun": lambda x: -g_kkt(x), "jac": lambda x: -A_kkt}
]
res_slsqp = minimize(
    f_kkt, [3, 3, 3], jac=grad_kkt, method="SLSQP", constraints=cons_kkt
)
print(
    np.round(res_slsqp.x, 6).tolist(),
    round(res_slsqp.fun, 6),
    np.round(res_slsqp.multipliers, 4).tolist(),
)
# Expected: [1.0, 2.0, 0.0] 5.0 [0.0, 0.0, 2.0, 4.0, 0.0]
linear = LinearConstraint(A_kkt, -np.inf, b_kkt)
res_tc = minimize(
    f_kkt, [3, 3, 3], jac=grad_kkt, hess=lambda x: 2 * np.eye(3),
    method="trust-constr", constraints=[linear],
    options={"gtol": 1e-10, "xtol": 1e-12},
)  # fmt: skip
print(np.round(res_tc.x, 3), np.round(res_tc.v[0], 3))
# Expected: [1. 2. 0.] [0.    0.    2.    4.    0.001]
lam_tc = np.maximum(res_tc.v[0], 0)
check = kkt_check(res_tc.x, lam_tc, grad_kkt, g_kkt, jac_kkt)
print([round(v, 3) for v in check.values()])
# Expected: [0.0, 0.0, 0.0, 0.0]
```

Both solvers find $\vx^*=(1,2,0)$, $f^*=5$ and the multipliers
$(0,0,2,4,0)$; the KKT checker quantifies how accurate the interior-point
answer is.

(sec-09-06-constrained)=
## Frank–Wolfe and the penalty method

### Frank–Wolfe with LP subproblems

{prf:ref}`alg-06-06-frank-wolfe` solves, at every iterate, the linear
program (LP)
$\min\{\nabla f(\vx^{(k)})^\T\vy:\vy\in\Omega\}$ ({eq}`eq-06-06-fw-lp`), here
with `linprog`. The optimal value of the LP gives the gap
$\Gamma(\vx^{(k)})=-\nabla f(\vx^{(k)})^\T(\vy^{(k)}-\vx^{(k)})$ and with it
the lower bound $f(\vx^{(k)})-\Gamma(\vx^{(k)})$ on $f^*$
({eq}`eq-06-06-fw-bound`). For a quadratic with Hessian $\mQ$ the exact step
along $\vd$ is $\min\{1,\Gamma/(\vd^\T\mQ\vd)\}$
({eq}`eq-06-06-fw-quadratic-step`).

```python
def frank_wolfe(f, grad, Q, A_ub, b_ub, x0, eps=1e-6, max_iter=1000):
    """Frank–Wolfe for a convex quadratic f over {y >= 0 : A_ub y <= b_ub}.

    Returns a list of (x_k, f(x_k), gap Gamma_k) for k = 0, 1, ...
    """
    x = np.asarray(x0, dtype=float)
    if np.any(A_ub @ x > b_ub + 1e-12) or np.any(x < 0):
        raise ValueError("x0 must be feasible")
    history = []
    for _ in range(max_iter):
        g = grad(x)
        lp = linprog(g, A_ub=A_ub, b_ub=b_ub, bounds=(0, None))
        if lp.status != 0:
            raise ValueError(f"direction-finding LP failed: {lp.message}")
        d = lp.x - x
        gap = float(-g @ d)
        history.append((x, f(x), gap))
        if gap <= eps:
            break
        # Exact step on the segment; f is linear along d if d'Qd = 0.
        curvature = float(d @ Q @ d)
        alpha = 1.0 if curvature <= 0 else min(1.0, gap / curvature)
        x = x + alpha * d
    return history


hist = frank_wolfe(
    f_fw,
    grad_fw,
    np.diag([2.0, 4.0]),
    np.array([[3.0, 2.0]]),
    np.array([6.0]),
    [0, 0],
    max_iter=600,
)
for k in range(5):
    x, fx, gap = hist[k]
    print(k, np.round(x, 4).tolist(), round(fx, 4), round(gap, 4))
# Expected:
#   0 [0.0, 0.0] 0.0 24.0
#   1 [0.0, 2.0] -8.0 10.0
#   2 [0.8333, 1.1667] -10.0833 3.3333
#   3 [0.6461, 1.5787] -10.4579 2.3596
#   4 [0.8804, 1.3055] -10.662 1.8559
x9, f9, gap9 = hist[9]
print(round(gap9, 4), round(f9 - gap9, 2), round(f9, 2))
# Expected: 0.9365 -11.99 -11.06
first = [
    next(k for k, (_, _, gap) in enumerate(hist) if gap < tol)
    for tol in (1, 0.1)
]
print(first, next(k for k, (_, fx, _) in enumerate(hist) if fx + 11.5 <= 0.01))
# Expected: [9, 105] 542
```

The iterates $(0,2)$, $(\tfrac56,\tfrac76)$, $(0.6461,1.5787)$,
$(0.8804,1.3055)$ and the bounds $[-11.99,-11.06]$ at $k=9$ are those of
{prf:ref}`ex-06-06-frank-wolfe`; the gap falls below 1 at $k=9$ and below 0.1
at $k=105$, and $f(\vx^{(k)})-f^*\le0.01$ needs $k=542$: the sublinear rate
of {prf:ref}`thm-06-06-fw-rate`.

### Quadratic penalty

{prf:ref}`alg-06-06-penalty` minimizes $\Phi_\rho(\vx)=f(\vx)+\rho\psi(\vx)$
with the quadratic penalty {eq}`eq-06-06-quadratic-penalty` for an increasing
sequence of $\rho$, each time starting from the previous minimizer (warm
start), and reports the multiplier estimates
$\lambda_i(\rho)=2\rho\max\{0,g_i(\vx(\rho))\}$ of
{eq}`eq-06-06-multiplier-estimates`. The unconstrained subproblems are solved
with BFGS, with the gradient of $\Phi_\rho$ supplied.

```python
def penalty_method(f, grad_f, g, jac_g, x0, rho0=1.0, factor=10.0, rounds=5):
    """Quadratic penalty method for min f s.t. g(x) <= 0.

    Returns rows (rho, x(rho), Phi(x(rho)), multiplier estimates).
    """
    if rho0 <= 0 or factor <= 1:
        raise ValueError("need rho0 > 0 and factor > 1")
    x, rho, rows = np.asarray(x0, dtype=float), rho0, []
    for _ in range(rounds):

        def phi(z, rho=rho):
            """Penalized objective Phi_rho."""
            return f(z) + rho * np.sum(np.maximum(0, g(z)) ** 2)

        def grad_phi(z, rho=rho):
            """Gradient of Phi_rho."""
            return grad_f(z) + 2 * rho * jac_g(z).T @ np.maximum(0, g(z))

        x = minimize(
            phi, x, jac=grad_phi, method="BFGS", options={"gtol": 1e-10}
        ).x
        rows.append((rho, x, phi(x), 2 * rho * np.maximum(0, g(x))))
        rho *= factor
    return rows
```

For {prf:ref}`ex-06-06-penalty-1d`, $\min x$ s.t. $2-x\le0$, the minimizers are
$x(\rho)=2-1/(2\rho)$, the penalized values $2-1/(4\rho)$ increase toward
$f^*=2$, and the multiplier estimate is 1 for every $\rho$.

```python
rows = penalty_method(
    lambda x: x[0], lambda x: np.array([1.0]),
    lambda x: np.array([2 - x[0]]), lambda x: np.array([[-1.0]]),
    x0=[0.0], rho0=0.5, factor=2, rounds=4,
)  # fmt: skip
for rho, x, phi_value, lam in rows:
    print(
        rho,
        np.round(x, 6).tolist(),
        round(phi_value, 6),
        np.round(lam, 6).tolist(),
    )
# Expected:
#   0.5 [1.0] 1.5 [1.0]
#   1.0 [1.5] 1.75 [1.0]
#   2.0 [1.75] 1.875 [1.0]
#   4.0 [1.875] 1.9375 [1.0]
```

The same code applied to the triangle of {prf:ref}`ex-06-06-frank-wolfe`
approaches $\vx^*=(1,1.5)$ from outside, and the estimate for the active
constraint $3x_1+2x_2\le6$ tends to its multiplier $\lambda=1$; the two
bounds $x_1,x_2\ge0$ are inactive and get the estimate 0.

```python
A_tri = np.array(
    [
        [3.0, 2.0],
        [-1.0, 0.0],
        [0.0, -1.0],
    ]
)
b_tri = np.array([6.0, 0.0, 0.0])
rows = penalty_method(
    f_fw, grad_fw, lambda x: A_tri @ x - b_tri, lambda x: A_tri, x0=[0.0, 0.0]
)
for rho, x, phi_value, lam in rows:
    print(
        rho,
        np.round(x, 5).tolist(),
        round(phi_value, 5),
        np.round(lam, 4).tolist(),
    )
# Expected:
#   1.0 [1.125, 1.54167] -11.72917 [0.9167, 0.0, 0.0]
#   10.0 [1.01351, 1.5045] -11.52477 [0.991, 0.0, 0.0]
#   100.0 [1.00136, 1.50045] -11.5025 [0.9991, 0.0, 0.0]
#   1000.0 [1.00014, 1.50005] -11.50025 [0.9999, 0.0, 0.0]
#   10000.0 [1.00001, 1.5] -11.50002 [1.0, 0.0, 0.0]
```

The violation $3x_1+2x_2-6$ shrinks like $1/\rho$, as
{ref}`sec-06-06-multipliers` explains.

(sec-09-06-scipy)=
## Least-squares estimation with SLSQP

The estimation problem {eq}`eq-06-01-regression-nlp` has six variables, the
equality $\beta_1+\beta_2+\beta_3=1$, bounds $0\le\beta_i\le1$ and
$0\le q_i\le\min_tQ_{ti}$, and a nonconvex objective. SLSQP finds a local
minimum from each start, so it is run from many random feasible starting
points (shares drawn uniformly from the simplex, quantities uniformly
between their bounds), and the best result is kept. The analytic gradient
of {ref}`sec-09-06-derivatives` makes each run fast and accurate.

```python
q_max = Q.min(axis=0).astype(float)
bounds = [(0, 1)] * 3 + [(0, u) for u in q_max]
shares = {
    "type": "eq",
    "fun": lambda v: v[:3].sum() - 1,
    "jac": lambda v: np.r_[1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
}
rng = np.random.default_rng(2024)
runs = []
for _ in range(200):
    v0 = np.r_[rng.dirichlet(np.ones(3)), rng.uniform(0, q_max)]
    res = minimize(
        sse, v0, jac=grad_sse, method="SLSQP", bounds=bounds,
        constraints=[shares], options={"ftol": 1e-12, "maxiter": 500},
    )  # fmt: skip
    if res.success:
        runs.append(res)
best = min(runs, key=lambda r: r.fun)
print(
    len(runs),
    f"{best.fun:.5g}",
    np.round(best.x[:3], 4).tolist(),
    np.round(best.x[3:], 2),
)
# Expected: 188 1.3547e+05 [0.2474, 0.4635, 0.2891] [89.85 16.67 37.32]
local_values = sorted({round(r.fun, -2) for r in runs})
print(
    sum(abs(r.fun - best.fun) < 1e-3 * best.fun for r in runs),
    local_values[:3],
)
# Expected: 41 [135500.0, 152300.0, 155600.0]
```

The best point is the estimate of {prf:ref}`ex-06-01-consumer-regression`:
$\vbeta^*\approx(0.2474,0.4635,0.2891)$,
$\vq^*\approx(89.85,16.67,37.32)$ units and $f^*\approx1.3547\times10^5$
$\$^2$. Only 41 of the 188 successful runs reach it; the others stop at
points with larger values (about $1.523\times10^5$ and $1.556\times10^5$
$\$^2$ and more), which shows why a nonconvex problem needs many starts.
The Hessian at the best point, computed from the analytic gradient, has a
negative eigenvalue (the objective is not convex), but projected onto the plane $\sum_i\beta_i=1$ it is positive definite, the
second-order condition of a strict local minimum.

```python
H_star = fd_hessian(grad_sse, best.x)
print(round(np.linalg.eigvalsh(H_star).min(), 2))
# Expected: -4.55
Z = np.zeros((6, 5))
Z[:3, :2] = [
    [1, 0],
    [-1, 1],
    [0, -1],
]
Z[3:, 2:] = np.eye(3)
print(np.linalg.eigvalsh(Z.T @ H_star @ Z).min() > 0)
# Expected: True
```

The columns of `Z` span the directions that keep $\sum_i\beta_i=1$
($\vd$ with $d_1+d_2+d_3=0$); every bound is inactive at the solution, so
no other direction has to be excluded.

(sec-09-06-cvxpy)=
## Convex problems with CVXPY

*This section requires CVXPY.* CVXPY lets a convex problem be written
almost as on paper: variables, an objective built from convex functions,
and a list of constraints. It checks that the problem follows its
**disciplined convex programming** (DCP) rules, transforms it into a conic
program and calls a solver. After `prob.solve()`, `prob.status`,
`prob.value`, the variables' `.value` and every constraint's `.dual_value`
are available. Results are rounded below because different solvers return
slightly different digits.

### The KKT example

For {prf:ref}`ex-06-03-kkt-example` the dual values of the vector constraint
$\mA\vx\le\vb$ are the multipliers $\lambda_i\ge0$ of the book: CVXPY's
Lagrangian for a constraint `expr <= rhs` is
$f+\vlambda^\T(\text{expr}-\text{rhs})$, the convention of
{prf:ref}`thm-06-03-kkt`.

```python
x = cp.Variable(3)
ineq = A_kkt @ x <= b_kkt
prob = cp.Problem(cp.Minimize(cp.sum_squares(x)), [ineq])
prob.solve()
print(prob.status, round(prob.value, 4), np.round(x.value, 4).tolist())
# Expected: optimal 5.0 [1.0, 2.0, 0.0]
print(np.round(ineq.dual_value, 4).tolist())
# Expected: [0.0, 0.0, 2.0, 4.0, 0.0]
```

### Economic dispatch and its multipliers

The dispatch problem with a capacity limit, {prf:ref}`ex-02-03-dispatch-capacity`,
has one inequality and one equality. For an equality `expr == rhs` CVXPY
adds $\mu(\text{expr}-\text{rhs})$ to the Lagrangian. The book writes the
demand constraint as $h(\vx)=475-x_1-x_2=0$, the negative of
`x[0] + x[1] - 475`, so the book's $\mu$ is the negative of CVXPY's dual
value.

```python
x = cp.Variable(2)
cost = 8 * x[0] + 0.01 * x[0] ** 2 + 10 * x[1] + 0.02 * x[1] ** 2
capacity = x[0] <= 300
demand = x[0] + x[1] == 475
prob = cp.Problem(cp.Minimize(cost), [capacity, demand])
prob.solve()
print(round(prob.value, 2), np.round(x.value, 2))
# Expected: 5662.5 [300. 175.]
print(
    round(float(capacity.dual_value), 3), round(-float(demand.dual_value), 3)
)
# Expected: 3.0 17.0
```

The dispatch $(300,175)$ MW costs \$5 662.50 per hour; the capacity is worth
$\lambda^*=3$ \$/MWh and the marginal cost of demand is $\mu^*=17$ \$/MWh.
Raising the capacity to 301 MW confirms the sensitivity interpretation of
{ref}`sec-02-03-sensitivity`:

```python
capacity_301 = x[0] <= 301
cp.Problem(cp.Minimize(cost), [capacity_301, demand]).solve()
print(round(5662.5 - cost.value, 2))
# Expected: 2.97
```

### Lagrange duality

Without the capacity limit ({prf:ref}`ex-02-03-dispatch-dual`) the dual
function {eq}`eq-02-03-dispatch-q` is a concave parabola in $\mu$. Maximizing it
numerically recovers $\mu^*=15$ \$/MWh and $d^*=p^*=5587.5$ \$/h, and CVXPY
returns the same price as the dual value of the demand constraint.

```python
def q_dispatch(mu: float) -> float:
    """Dual function q(mu) of the dispatch problem without capacity limit."""
    return 475 * mu - (mu - 8) ** 2 / 0.04 - (mu - 10) ** 2 / 0.08


res = minimize_scalar(lambda mu: -q_dispatch(mu))
print(
    round(res.x, 6),
    round(q_dispatch(res.x), 4),
    q_dispatch(10),
    q_dispatch(20),
)
# Expected: 15.0 5587.5 4650.0 4650.0
prob = cp.Problem(cp.Minimize(cost), [demand])
prob.solve()
print(round(prob.value, 4), round(-float(demand.dual_value), 4))
# Expected: 5587.5 15.0
```

The least-norm problem of {prf:ref}`ex-02-03-least-norm`,
$\min\{\vx^\T\vx:\mA\vx=\vb\}$, has the closed form {eq}`eq-02-03-least-norm`.
Here the book's constraint function is $\vh(\vx)=\mA\vx-\vb$, the same sign
as CVXPY's, so the dual values are $\vmu^*=(0,-2)$ directly.

```python
A_ln = np.array(
    [
        [1.0, 1.0, 1.0],
        [1.0, 2.0, 3.0],
    ]
)
b_ln = np.array([6.0, 14.0])
x = cp.Variable(3)
eq = A_ln @ x == b_ln
prob = cp.Problem(cp.Minimize(cp.sum_squares(x)), [eq])
prob.solve()
print(
    round(prob.value, 4),
    np.round(x.value, 4).tolist(),
    (np.round(eq.dual_value, 4) + 0.0).tolist(),
)
# Expected: 14.0 [1.0, 2.0, 3.0] [0.0, -2.0]
print(A_ln.T @ np.linalg.solve(A_ln @ A_ln.T, b_ln))
# Expected: [1. 2. 3.]
```

### What CVXPY refuses

The least-squares estimation problem is not convex, and CVXPY says so
before calling any solver: a product of two variables, such as
$\beta_1q_2$, violates the DCP rules. With the shares fixed, however, the
residuals are affine in $\vq$ and the problem in $\vq$ alone is a convex quadratic
program.
Fixing $\vbeta$ at the SLSQP estimate and solving for $\vq$
returns the same committed quantities, which confirms that the SLSQP point
is optimal in $\vq$ for its shares.

```python
beta, q = cp.Variable(3), cp.Variable(3)
print((beta[0] * q[1]).is_dcp())
# Expected: False
beta_star = best.x[:3]
m = E - prices @ q
fitted = cp.multiply(
    prices, np.ones((4, 1)) @ cp.reshape(q, (1, 3), order="C")
)
fitted = fitted + cp.reshape(m, (4, 1), order="C") @ beta_star.reshape(1, 3)
prob = cp.Problem(
    cp.Minimize(cp.sum_squares(prices * Q - fitted)), [q >= 0, q <= q_max]
)
prob.solve()
print(prob.is_dcp(), f"{prob.value:.5g}", np.round(q.value, 2))
# Expected: True 1.3547e+05 [89.85 16.67 37.32]
```

## Summary

- Analytic gradients and Hessians are checked by central differences
  ({ref}`sec-09-06-derivatives`), including the degree-four least-squares
  objective {eq}`eq-06-01-sse`.
- Bisection (6 derivative values, $\hat x=4.7109$), golden-section search (final
  interval $[4.6869,4.7264]$), Newton's method (fast from 4, a maximizer
  from 3) and backtracking ($\alpha=2$) reproduce {doc}`06-04-line-search`.
- Steepest descent zig-zags with the ratio 0.1018 and needs 17 iterations
  on {prf:ref}`ex-06-05-steepest-descent`; Newton's method needs one.
- The KKT checker and the enumeration of 32 active sets confirm
  $\vx^*=(1,2,0)$, $\vlambda^*=(0,0,2,4,0)$; SLSQP and trust-constr return
  the same multipliers.
- Frank–Wolfe with `linprog` subproblems reproduces
  {prf:ref}`ex-06-06-frank-wolfe` and its slow convergence; the quadratic
  penalty approaches the solution from outside with multiplier estimates
  $2\rho\max\{0,g_i\}$.
- Multi-start SLSQP finds the least-squares estimate of
  {prf:ref}`ex-06-01-consumer-regression`; CVXPY rejects the bilinear model
  but solves its convex part.
- CVXPY's dual values equal the book's multipliers for constraints written
  as $g\le0$; for equalities the sign depends on how $h$ is written
  ($\mu^*=17$ and $15$ \$/MWh in the dispatch examples, $\vmu^*=(0,-2)$ for
  the least-norm problem).

## References

- Armijo, L. (1966). Minimization of functions having Lipschitz continuous
  first partial derivatives. *Pacific Journal of Mathematics, 16*(1), 1–3.
  <https://doi.org/10.2140/pjm.1966.16.1>
- Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2006). *Nonlinear
  programming: Theory and algorithms* (3rd ed.). Wiley. (Chapter 4: KKT
  conditions; Chapters 8–10: algorithms.)
- Boyd, S., & Vandenberghe, L. (2004). *Convex optimization*. Cambridge
  University Press. <https://doi.org/10.1017/CBO9780511804441> (Chapters 4–5:
  convex problems and duality; Chapter 9: unconstrained minimization.)
- Brent, R. P. (1973). *Algorithms for minimization without derivatives*.
  Prentice-Hall. (Chapter 5: Brent's one-dimensional minimizer.)
- Cauchy, A. (1847). Méthode générale pour la résolution des systèmes
  d'équations simultanées. *Comptes Rendus de l'Académie des Sciences, 25*,
  536–538.
- Courant, R. (1943). Variational methods for the solution of problems of
  equilibrium and vibrations. *Bulletin of the American Mathematical Society,
  49*(1), 1–23. <https://doi.org/10.1090/S0002-9904-1943-07818-4>
- Fiacco, A. V., & McCormick, G. P. (1968). *Nonlinear programming:
  Sequential unconstrained minimization techniques*. Wiley.
- Frank, M., & Wolfe, P. (1956). An algorithm for quadratic programming.
  *Naval Research Logistics Quarterly, 3*(1–2), 95–110.
  <https://doi.org/10.1002/nav.3800030109>
- Kiefer, J. (1953). Sequential minimax search for a maximum. *Proceedings
  of the American Mathematical Society, 4*(3), 502–506.
  <https://doi.org/10.1090/S0002-9939-1953-0055639-3>
- Kuhn, H. W., & Tucker, A. W. (1951). Nonlinear programming. In J. Neyman
  (Ed.), *Proceedings of the Second Berkeley Symposium on Mathematical
  Statistics and Probability* (pp. 481–492). University of California
  Press.
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Chapter 8: line search and steepest descent; Chapter 13: penalty methods.)
- Nocedal, J., & Wright, S. J. (2006). *Numerical optimization* (2nd ed.).
  Springer. <https://doi.org/10.1007/978-0-387-40065-5> (Chapter 3: line
  search methods; Chapter 6: quasi-Newton methods; Chapter 8: finite
  differences; Chapter 17: penalty methods; Chapter 18: sequential
  quadratic programming.)
