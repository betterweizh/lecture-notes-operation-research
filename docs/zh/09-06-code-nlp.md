# 用 Python 求解非线性规划

## 概述

本章实现{doc}`06-01-nlp-models`至{doc}`06-06-constrained-algorithms`中的非线性规划（nonlinear programming, NLP）方法，并把 SciPy 和 CVXPY 的求解器（solver）用于本书的例题。首先给出梯度（gradient）和黑塞矩阵（Hessian matrix，见{doc}`01-02-matrix-calculus`）的解析表达式，并用有限差分加以检验。然后编写{doc}`06-04-line-search`中的一维线搜索（line search）方法（二分法（bisection search）、黄金分割法（golden-section search）、牛顿法（Newton's method）和 Armijo 回溯），{doc}`06-05-unconstrained-algorithms`中的最速下降法（steepest descent method）和牛顿法，以及{doc}`06-06-constrained-algorithms`中的 Frank–Wolfe 方法和二次罚函数法，其中线性子问题用 `linprog` 求解。一个小型的 **KKT 检验程序**（KKT checker）度量一个点和一组乘子距离满足{doc}`06-03-kkt-conditions`中的 KKT 条件（Karush–Kuhn–Tucker conditions）有多远：平稳性条件（stationarity）、可行性、互补松弛性（complementary slackness）和符号条件；它被应用于该章的 KKT 例子，该例的 32 个有效集（active set）也被全部枚举。`scipy.optimize.minimize` 的通用求解器（BFGS、SLSQP、trust-constr）被用于同样的例子以及{doc}`06-01-nlp-models`中的最小二乘估计（least-squares estimation）。最后一节用 CVXPY 求解凸问题，把它的对偶值与 KKT 乘子相比较，并重现{doc}`02-03-convex-optimization-duality`中的拉格朗日对偶（Lagrange duality）例子；这一节需要 CVXPY（任何较新版本及其默认求解器均可）。构建本书时不运行这些代码：请把各代码块按顺序复制到脚本或笔记本中，并将打印结果与 `# Expected:` 注释相比较。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 写出本书中各函数的梯度和黑塞矩阵，并用中心差分（central difference）加以检验；
- 实现二分法、黄金分割法、牛顿法和回溯法，以及最速下降法、牛顿法、Frank–Wolfe 方法和二次罚函数法；
- 数值检验 KKT 条件并枚举有效集；
- 带界约束、线性约束和非线性约束调用 `scipy.optimize.minimize`，读取其乘子，并对非凸问题使用多初始点求解；
- 在 CVXPY 中建立凸问题的模型，并按本书的符号约定解释其对偶值。
:::

(sec-09-06-setup)=
## 环境准备

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

`TAU` 是黄金分割比 $\tau=(\sqrt5-1)/2\approx0.618$。凡是{doc}`06-04-line-search`至{doc}`06-06-constrained-algorithms`中的迭代表给出四位有效数字的地方，都用 `round` 或 `np.round` 显式地舍入。

(sec-09-06-derivatives)=
## 梯度与黑塞矩阵

下面所有方法都需要 $\nabla f$，有些还需要 $\nabla^2f$。解析导数计算快且精确，但容易写错；有限差分检验可以发现这些错误。**中心差分**

$$
\frac{\partial f}{\partial x_j}(\vx)\approx\frac{f(\vx+h\ve_j)-f(\vx-h\ve_j)}{2h}
$$

的误差为 $h^2$ 阶，因此取 $h\approx10^{-6}$（相对于 $\vx$ 的大小）时，它与正确的梯度大约有八位数字相符。把它应用于解析梯度，就得到黑塞矩阵，然后再将其对称化。

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

{prf:ref}`ex-06-05-steepest-descent`中的二次函数 $f(\vx)=(x_1-2)^2+10(x_2-2)^2$ 满足 $\nabla f(\vx)=(2x_1-4,\ 20x_2-40)$，其黑塞矩阵为常数 $\diag(2,20)$；{prf:ref}`ex-06-06-frank-wolfe`的目标函数 $f(\vx)=-5x_1+x_1^2-8x_2+2x_2^2$ 满足 $\nabla f=(2x_1-5,\ 4x_2-8)$ 和 $\nabla^2f=\diag(2,4)$。

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

更难的情形是{prf:ref}`ex-06-01-consumer-regression`中的最小二乘目标函数{eq}`eq-06-01-sse`，它是六个参数 $(\vbeta,\vq)$ 的四次多项式。利用{eq}`eq-06-01-residual`中的残差 $r_{ti}$ 以及 $m_t=E_t-\sum_jp_{tj}q_j$（第 $t$ 月的可自由支配金额），

$$
\begin{aligned}
\frac{\partial r_{ti}}{\partial\beta_k}&=-[i=k]\,m_t,\\
\frac{\partial r_{ti}}{\partial q_j}&=-[i=j]\,p_{ti}+\beta_ip_{tj},
\end{aligned}
$$

且 $\nabla f=2\sum_{t,i}r_{ti}\nabla r_{ti}$。在 NumPy 中，关于 $t$ 和 $i$ 的求和变成矩阵乘积。

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

相对差小于 $10^{-7}$，因此解析梯度是可信的；{ref}`sec-09-06-scipy`中的 SLSQP 将使用它。

(sec-09-06-line-search)=
## 一维搜索

四种方法都在{doc}`06-04-line-search`的问题上试验：在 $[3,6]$ 上极小化 $f(x)=\sin x$，其极小点为 $x^*=3\pi/2\approx4.712389$。

### 二分法

{prf:ref}`alg-06-04-bisection`根据 $f'$ 在中点处的符号把区间减半。该函数记录{numref}`tab-06-04-bisection`中的各行。

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

六次导数计算给出 $\hat x=4.7109375$，误差为 $0.001451$；减半的次数为 $\lceil\log_2(3/0.05)\rceil=6$，与{eq}`eq-06-04-bisection-count`的预测一致。

### 黄金分割法

{prf:ref}`alg-06-04-golden-section`在区间的 $1-\tau$ 和 $\tau$ 分点处保留两个内点，并在每一步中重复使用其中之一，因此每次缩小区间只需计算一个新的函数值。

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

循环经过九次缩小后停止，得到{prf:ref}`ex-06-04-golden-sine`中的区间 $[4.6869,4.7264]$ 和估计值 $4.7067$，共计算了 $N=10$ 个函数值。最后一次缩小所产生的新试探点没有被计算，因为终止准则在用到它之前就结束了搜索。

### 牛顿法

{prf:ref}`alg-06-04-newton-1d`按 $x_{k+1}=x_k-f'(x_k)/f''(x_k)$ 迭代。该函数有意不加任何保护措施（safeguard），以便能看到{prf:ref}`ex-06-04-newton-sine` (b) 中的失败情形；它只拒绝零曲率。

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

从 $x_0=4$ 出发，三步就到达 $3\pi/2$；从 $x_0=3$ 出发，第一步跳到 $-4.0153$，迭代收敛到 $-3\pi/2$，这是 $\sin x$ 在 $[3,6]$ 之外的一个极大点（函数值为 1），而终止准则却报告成功。

### 回溯法（Armijo）

{prf:ref}`alg-06-04-backtracking`不断缩小试探步长，直到 Armijo 条件{eq}`eq-06-04-armijo`成立。它适用于 $\R^n$；{ref}`sec-06-04-inexact`中的示例是一维情形 $x=3$，$d=-f'(3)$，取 $c_1=0.25$，$\bar\alpha=4$，$\gamma=\tfrac12$。

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

步长 $\alpha=4$ 被拒绝，$\alpha=2$ 被接受，得到 $x=4.98$。作为比较，SciPy 的有界标量极小化程序（Brent 方法，即以黄金分割步作为保护措施的抛物线插值；Brent, 1973）求得 $3\pi/2$，精确到六位数字。

```python
res = minimize_scalar(math.sin, bounds=(3, 6), method="bounded")
print(round(res.x, 6), res.nfev)
# Expected: 4.71239 8
```

(sec-09-06-unconstrained)=
## 最速下降法与牛顿法

{prf:ref}`alg-06-05-steepest-descent`取方向 $-\nabla f$，并由线搜索确定步长。线搜索以函数 `step(x, d)` 的形式传入，这样二次函数的精确步长{eq}`eq-06-05-sd-step`与 Armijo 回溯可以互换。

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

在{prf:ref}`ex-06-05-steepest-descent`上，第一步到达 $(-3.3922,\ 2.0647)$，每一步都使 $f-f^*$ 乘以 $0.1018$，梯度判据 $\norm{\nabla f}\le10^{-6}$ 在 $k=17$ 时满足。

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

相继的梯度相互正交（{eq}`eq-06-05-orthogonal`）；从最坏情形的初始点 $(-8,1)$ 出发，比值为{prf:ref}`thm-06-05-sd-rate`中的 $((\kappa-1)/(\kappa+1))^2=0.6694$，需要多得多的迭代次数。用回溯法代替精确步长，该方法仍然收敛，所需代价取决于各个常数。

```python
def armijo_step(f, grad, alpha0=1.0, c1=1e-4):
    """Line search for steepest_descent: Armijo backtracking from alpha0."""
    return lambda x, d: backtracking(f, grad, x, d, alpha0=alpha0, c1=c1)[0]


xs_bt = steepest_descent(grad_quad, [-4, -3], armijo_step(f_quad, grad_quad))
print(len(xs_bt) - 1, np.round(xs_bt[-1], 6).tolist())
# Expected: 64 [2.0, 2.0]
```

{prf:ref}`alg-06-05-newton`求解 $\nabla^2f(\vx)\vd=-\nabla f(\vx)$（从不对黑塞矩阵求逆），并取完整的步长。对于二次函数，一步就足够了（见{prf:ref}`ex-06-05-newton-quadratic`）。

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

SciPy 的 `minimize` 提供这两类方法：`BFGS` 由梯度构造拟牛顿（quasi-Newton）近似（这里由于它的线搜索不精确，需要几次迭代），`Newton-CG` 则使用黑塞矩阵。

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
## 检验 KKT 条件

对于问题 $\min f(\vx)$，满足 $\vg(\vx)\le\vzero$，$\vh(\vx)=\vzero$，{prf:ref}`thm-06-03-kkt`中的 KKT 条件分为四组：平稳性条件{eq}`eq-06-03-kkt-stationarity`、原始可行性{eq}`eq-06-03-kkt-primal`、符号条件{eq}`eq-06-03-kkt-dual`和互补松弛性{eq}`eq-06-03-kkt-cs`。检验程序对每一组返回一个数，即最严重违反的大小；KKT 点给出四个接近零的值。约束函数返回向量，其雅可比矩阵（Jacobian matrix）每个约束一行（$\nabla g_i^\T$）。

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

{prf:ref}`ex-06-03-kkt-example`在五个线性约束下极小化 $x_1^2+x_2^2+x_3^2$，这里把约束写成 $\mA\vx-\vb\le\vzero$，$\mA$ 的每一行对应一个约束 $g_i$。

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

点 $\vx^*=(1,2,0)$ 与 $\vlambda^*=(0,0,2,4,0)$ 精确地满足全部条件。{numref}`tab-06-03-active-set-cases`中有效集 $\{1,4\}$ 对应的点 $(1.5,2,0)$ 满足平稳性、可行性和互补性，但 $\lambda_1=-1.5<0$：它只违反符号条件。

### 枚举有效集

由于约束是线性的、$f$ 是二次的，{prf:ref}`alg-06-03-active-set`可以对全部 $2^5=32$ 个子集 $S$ 运行：对 $i\in S$ 令 $g_i=0$，否则令 $\lambda_i=0$，则平稳性条件 $2\vx+\mA_S^\T\vlambda_S=\vzero$ 与 $\mA_S\vx=\vb_S$ 构成关于 $(\vx,\vlambda_S)$ 的一个线性方程组。若最小二乘解的残差不为零，则该方程组**不相容**；否则，把解归类为不可行（某个 $g_i>0$，$i\notin S$）、有负乘子，或为 KKT 点。

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

计数 14、8、8、2 以及被接受的两个集合 $\{3,4\}$ 和 $\{3,4,5\}$ 都与{prf:ref}`ex-06-03-kkt-example`相同；两者都给出 $\vx^*=(1,2,0)$ 与 $\vlambda^*=(0,0,2,4,0)$。

### SciPy 的约束求解器

带 `method="SLSQP"`（序列二次规划，sequential quadratic programming）的 `minimize` 以字典的形式接受约束。SciPy 把不等式写成 $c(\vx)\ge0$，与本书的 $g(\vx)\le0$ 相反，因此传入的是 $c=-g$。较新的 SciPy 版本（1.16 及以后）还在 `res.multipliers` 中返回 SLSQP 的乘子，其符号约定与本书相同。`method="trust-constr"` 接受 `LinearConstraint` $\vl\le\mA\vx\le\vu$，并在 `res.v` 中返回其乘子（每个约束对象一个数组）；它是一种内点法（interior-point method），因此它的解只在容差范围内逼近边界。

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

两个求解器都求得 $\vx^*=(1,2,0)$、$f^*=5$ 和乘子 $(0,0,2,4,0)$；KKT 检验程序量化了内点法答案的精度。

(sec-09-06-constrained)=
## Frank–Wolfe 方法与罚函数法

### 带 LP 子问题的 Frank–Wolfe 方法

{prf:ref}`alg-06-06-frank-wolfe`在每个迭代点求解线性规划问题（linear program, LP）$\min\{\nabla f(\vx^{(k)})^\T\vy:\vy\in\Omega\}$（{eq}`eq-06-06-fw-lp`），这里用 `linprog` 求解。LP 的最优值给出间隙 $\Gamma(\vx^{(k)})=-\nabla f(\vx^{(k)})^\T(\vy^{(k)}-\vx^{(k)})$，并由此给出 $f^*$ 的下界 $f(\vx^{(k)})-\Gamma(\vx^{(k)})$（{eq}`eq-06-06-fw-bound`）。对于黑塞矩阵为 $\mQ$ 的二次函数，沿 $\vd$ 的精确步长为 $\min\{1,\Gamma/(\vd^\T\mQ\vd)\}$（{eq}`eq-06-06-fw-quadratic-step`）。

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

迭代点 $(0,2)$、$(\tfrac56,\tfrac76)$、$(0.6461,1.5787)$、$(0.8804,1.3055)$ 以及 $k=9$ 时的界 $[-11.99,-11.06]$ 都与{prf:ref}`ex-06-06-frank-wolfe`相同；间隙在 $k=9$ 时降到 1 以下，在 $k=105$ 时降到 0.1 以下，而要使 $f(\vx^{(k)})-f^*\le0.01$ 需要 $k=542$：这就是{prf:ref}`thm-06-06-fw-rate`中的次线性收敛速度。

### 二次罚函数

{prf:ref}`alg-06-06-penalty`对递增的 $\rho$ 序列，用二次罚函数（quadratic penalty function）{eq}`eq-06-06-quadratic-penalty`极小化 $\Phi_\rho(\vx)=f(\vx)+\rho\psi(\vx)$，每次都从上一个极小点出发（热启动，warm start），并报告{eq}`eq-06-06-multiplier-estimates`中的乘子估计值 $\lambda_i(\rho)=2\rho\max\{0,g_i(\vx(\rho))\}$。无约束子问题用 BFGS 求解，并提供 $\Phi_\rho$ 的梯度。

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

对于{prf:ref}`ex-06-06-penalty-1d`，即 $\min x$，满足 $2-x\le0$，极小点为 $x(\rho)=2-1/(2\rho)$，带罚目标函数值 $2-1/(4\rho)$ 递增趋于 $f^*=2$，且对每个 $\rho$ 乘子估计值都为 1。

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

把同样的代码用于{prf:ref}`ex-06-06-frank-wolfe`中的三角形区域，迭代点从外部逼近 $\vx^*=(1,1.5)$，有效约束 $3x_1+2x_2\le6$ 的估计值趋于其乘子 $\lambda=1$；两个界约束 $x_1,x_2\ge0$ 是非有效约束，估计值为 0。

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

如{ref}`sec-06-06-multipliers`所述，违反量 $3x_1+2x_2-6$ 按 $1/\rho$ 的速度缩小。

(sec-09-06-scipy)=
## 用 SLSQP 进行最小二乘估计

估计问题{eq}`eq-06-01-regression-nlp`有六个变量、等式约束 $\beta_1+\beta_2+\beta_3=1$、界约束 $0\le\beta_i\le1$ 和 $0\le q_i\le\min_tQ_{ti}$，以及一个非凸的目标函数。SLSQP 从每个初始点出发只能找到一个局部极小点，因此从许多随机的可行初始点出发运行它（份额在单纯形上均匀抽取，数量在其上下界之间均匀抽取），并保留最好的结果。{ref}`sec-09-06-derivatives`中的解析梯度使每次运行既快又准确。

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

最好的点就是{prf:ref}`ex-06-01-consumer-regression`中的估计值：$\vbeta^*\approx(0.2474,0.4635,0.2891)$，$\vq^*\approx(89.85,16.67,37.32)$ 单位，$f^*\approx1.3547\times10^5$ $\$^2$。188 次成功运行中只有 41 次到达该点；其余的停在函数值更大的点上（约 $1.523\times10^5$ 和 $1.556\times10^5$ $\$^2$ 以及更大），这说明了为什么非凸问题需要许多初始点。在最好的点处由解析梯度计算出的黑塞矩阵有一个负特征值（目标函数不是凸的），但投影到平面 $\sum_i\beta_i=1$ 上后是正定的，这正是严格局部极小点的二阶条件。

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

`Z` 的各列张成保持 $\sum_i\beta_i=1$ 的方向（满足 $d_1+d_2+d_3=0$ 的 $\vd$）；在该解处所有的界约束都是非有效的，因此不需要排除其他方向。

(sec-09-06-cvxpy)=
## 用 CVXPY 求解凸问题

*本节需要 CVXPY。* CVXPY 让凸问题几乎可以像在纸上一样写出：变量、由凸函数构成的目标函数，以及一个约束列表。它检查问题是否遵循其**规范凸规划**（disciplined convex programming, DCP）规则，把问题转换为锥规划（conic program），然后调用求解器。执行 `prob.solve()` 之后，可以得到 `prob.status`、`prob.value`、各变量的 `.value` 以及每个约束的 `.dual_value`。由于不同求解器返回的数字略有差异，下面的结果都经过舍入。

### KKT 例子

对于{prf:ref}`ex-06-03-kkt-example`，向量约束 $\mA\vx\le\vb$ 的对偶值就是本书中的乘子 $\lambda_i\ge0$：对于约束 `expr <= rhs`，CVXPY 的拉格朗日函数为 $f+\vlambda^\T(\text{expr}-\text{rhs})$，这正是{prf:ref}`thm-06-03-kkt`的约定。

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

### 经济调度及其乘子

带容量限制的调度问题（{prf:ref}`ex-02-03-dispatch-capacity`）有一个不等式约束和一个等式约束。对于等式 `expr == rhs`，CVXPY 在拉格朗日函数中加上 $\mu(\text{expr}-\text{rhs})$。本书把需求约束写成 $h(\vx)=475-x_1-x_2=0$，它是 `x[0] + x[1] - 475` 的相反数，因此本书中的 $\mu$ 是 CVXPY 对偶值的相反数。

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

调度方案 $(300,175)$ 兆瓦的费用为每小时 5 662.50 美元；容量的价值为 $\lambda^*=3$ 美元/兆瓦时，需求的边际成本为 $\mu^*=17$ 美元/兆瓦时。把容量提高到 301 兆瓦，证实了{ref}`sec-02-03-sensitivity`中的灵敏度解释：

```python
capacity_301 = x[0] <= 301
cp.Problem(cp.Minimize(cost), [capacity_301, demand]).solve()
print(round(5662.5 - cost.value, 2))
# Expected: 2.97
```

### 拉格朗日对偶

没有容量限制时（{prf:ref}`ex-02-03-dispatch-dual`），对偶函数{eq}`eq-02-03-dispatch-q`是关于 $\mu$ 的凹抛物线。对它作数值极大化，得到 $\mu^*=15$ 美元/兆瓦时和 $d^*=p^*=5587.5$ 美元/小时，CVXPY 也把同样的价格作为需求约束的对偶值返回。

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

{prf:ref}`ex-02-03-least-norm`中的最小范数问题 $\min\{\vx^\T\vx:\mA\vx=\vb\}$ 有闭式解{eq}`eq-02-03-least-norm`。这里本书的约束函数为 $\vh(\vx)=\mA\vx-\vb$，与 CVXPY 的符号相同，因此对偶值直接就是 $\vmu^*=(0,-2)$。

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

### CVXPY 拒绝的问题

最小二乘估计问题不是凸的，CVXPY 在调用任何求解器之前就会指出这一点：两个变量之积（例如 $\beta_1q_2$）违反 DCP 规则。然而，当份额固定时，残差是 $\vq$ 的仿射函数，仅关于 $\vq$ 的问题是一个凸二次规划。把 $\vbeta$ 固定为 SLSQP 的估计值并对 $\vq$ 求解，得到相同的必需购买量，这证实了 SLSQP 所得的点在其份额下关于 $\vq$ 是最优的。

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

## 小结

- 解析梯度和黑塞矩阵用中心差分检验（见{ref}`sec-09-06-derivatives`），包括四次的最小二乘目标函数{eq}`eq-06-01-sse`。
- 二分法（6 次导数计算，$\hat x=4.7109$）、黄金分割法（最终区间 $[4.6869,4.7264]$）、牛顿法（从 4 出发收敛很快，从 3 出发得到极大点）和回溯法（$\alpha=2$）重现了{doc}`06-04-line-search`中的结果。
- 最速下降法在{prf:ref}`ex-06-05-steepest-descent`上以比值 0.1018 呈锯齿状前进，需要 17 次迭代；牛顿法只需一次。
- KKT 检验程序和 32 个有效集的枚举证实了 $\vx^*=(1,2,0)$、$\vlambda^*=(0,0,2,4,0)$；SLSQP 和 trust-constr 返回相同的乘子。
- 使用 `linprog` 子问题的 Frank–Wolfe 方法重现了{prf:ref}`ex-06-06-frank-wolfe`及其缓慢的收敛；二次罚函数法从外部逼近解，乘子估计值为 $2\rho\max\{0,g_i\}$。
- 多初始点的 SLSQP 求得{prf:ref}`ex-06-01-consumer-regression`中的最小二乘估计；CVXPY 拒绝双线性模型，但能求解其中的凸部分。
- 对于写成 $g\le0$ 的约束，CVXPY 的对偶值等于本书中的乘子；对于等式约束，符号取决于 $h$ 的写法（调度例子中 $\mu^*=17$ 和 $15$ 美元/兆瓦时，最小范数问题中 $\vmu^*=(0,-2)$）。

## 参考文献

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
