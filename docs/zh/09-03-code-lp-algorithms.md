# 用 Python 实现线性规划算法

## 概述

{doc}`09-02-code-lp-modeling`把计算交给了求解器（solver）。本章打开这个黑箱，用几百行纯 Python 代码实现{doc}`03-02-lp-geometry`至{doc}`03-07-sensitivity-analysis`中的算法：基本解（basic solution）的枚举、采用 Dantzig 规则和 Bland 规则的表格单纯形法（simplex method）、两阶段法（two-phase method）、采用逆矩阵的乘积形式（product form of the inverse）的修正单纯形法（revised simplex method）、对偶单纯形法（dual simplex method），以及灵敏度分析（sensitivity analysis）中费用与右端项（right-hand side）的变化范围分析（ranging）。代码严格遵循本书的约定：单纯形表的目标行放在最后，标记为 $-z$，其中是检验数（reduced cost）$\bar c_j$，右端项列中是 $-z$；最小化问题在所有 $\bar c_j\ge0$ 时达到最优，最大化问题在所有 $\bar c_j\le0$ 时达到最优。所有运算都用 `fractions.Fraction` **精确**进行，因此第三部分的每一张单纯形表都能逐个元素地重现：{prf:ref}`ex-03-03-simplex-run`中的产品组合计算过程、Beale 的循环（cycling）例子（{prf:ref}`ex-03-03-cycling`）、两阶段法的例子（{prf:ref}`ex-03-04-two-phase`）、修正单纯形法的计算过程（{prf:ref}`ex-03-06-revised-run`）、{doc}`03-06-revised-dual-simplex`和{doc}`03-07-sensitivity-analysis`中的对偶单纯形法例子，以及{numref}`tab-03-07-ranges`中的变化范围。每个结果都用 `linprog(method="highs")` 加以核对。

这些实现用于学习和核对手算结果。它们是稠密的、精确的、缓慢的；实用的求解程序采用修正形式、稀疏存储和浮点运算，并且在数值方面要谨慎得多（见{doc}`09-01-code-setup-solvers`）。与第九部分其余各章一样，构建本书时不运行这些代码：请把各代码块按顺序复制到脚本或笔记本中，并将输出与 `# Expected:` 注释相比较。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 用精确分数进行计算，并按本书的格式打印单纯形表；
- 枚举 $\mA\vx=\vb$ 的基本解，并找出最优的一个；
- 实现主元变换（pivot operation）、进基规则与出基规则，以及求解最小化和最大化问题的表格单纯形法，并检测循环；
- 实现第一阶段（Phase I）和第二阶段（Phase II），包括移出取值为零而仍留在基中的人工变量（artificial variable）；
- 实现采用 eta 矩阵（eta matrix）的修正单纯形法以及对偶单纯形法；
- 由最优基计算费用与右端项的变化范围，并用求解器加以验证。
:::

所有代码块都使用以下导入语句。

```python
import itertools
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

np.set_printoptions(precision=3, suppress=True, linewidth=120)
```

(sec-09-03-exact)=
## 精确运算与打印

`Fraction` 把有理数存储为一对整数，并使 `+ - * /` 保持精确：`Fraction(1, 3) + Fraction(1, 6)` 等于 `1/2`，而不是 `0.49999999999999994`。浮点单纯形程序在每一次符号判断中都需要容差（tolerance，见{ref}`sec-09-01-tolerances`）；使用分数时，判断 $\bar c_j<0$ 和 $\bar a_{is}>0$ 都是精确的，退化（degenerate）情形下的并列也是真正的并列。矩阵是行的列表，每一行是分数的列表。三个小的辅助函数分别用于转换数据、用高斯–若尔当消元法（Gauss–Jordan elimination，见{doc}`01-01-linear-algebra`）精确求解线性方程组，以及打印矩阵和单纯形表。

```python
def to_fractions(rows) -> list[list[Fraction]]:
    """Convert a matrix (list of rows of numbers or strings) to Fractions."""
    return [[Fraction(v) for v in row] for row in rows]


def gauss_jordan_inverse(M: list[list[Fraction]]) -> list[list[Fraction]]:
    """Exact inverse of a square matrix by Gauss-Jordan elimination.

    Parameters
    ----------
    M : list of list of Fraction
        Square matrix.

    Returns
    -------
    list of list of Fraction
        The inverse matrix.

    Raises
    ------
    ValueError
        If M is not square or is singular.
    """
    n = len(M)
    if any(len(row) != n for row in M):
        raise ValueError("matrix must be square")
    # Augment M with the identity and reduce the left half to I.
    work = [
        list(row) + [Fraction(int(i == j)) for j in range(n)]
        for i, row in enumerate(M)
    ]
    for col in range(n):
        pivot_row = next((r for r in range(col, n) if work[r][col] != 0), None)
        if pivot_row is None:
            raise ValueError("matrix is singular")
        work[col], work[pivot_row] = work[pivot_row], work[col]
        p = work[col][col]
        work[col] = [v / p for v in work[col]]
        for r in range(n):
            if r != col and work[r][col] != 0:
                f = work[r][col]
                work[r] = [a - f * b for a, b in zip(work[r], work[col])]
    return [row[n:] for row in work]


def mat_vec(M, v):
    """Product of a matrix (list of rows) and a vector (list)."""
    return [sum((a * x for a, x in zip(row, v)), Fraction(0)) for row in M]


def show(M) -> None:
    """Print a matrix of Fractions, one row per line, aligned."""
    cells = [[str(v) for v in row] for row in M]
    width = max(len(c) for row in cells for c in row)
    for row in cells:
        print("  ".join(c.rjust(width) for c in row))


B = to_fractions(
    [
        [1, 2, 0],
        [0, 1, 1],
        [1, 2, 1],
    ]
)
show(gauss_jordan_inverse(B))
# Expected (the inverse of the optimal basis of eq-03-07-running-lp):
# -1  -2   2
#  1   1  -1
# -1   0   1
```

**单纯形表**（simplex tableau）存储为 $m+k$ 个长度为 $n+1$ 的行组成的列表：先是 $m$ 个约束行 $[\bar a_{i1},\dots,\bar a_{in}\mid\bar b_i]$，然后是 $k\ge1$ 个目标行，$-z$ 行在最后（在第一阶段，$-\phi$ 行排在它之后，与{doc}`03-04-two-phase-big-m`相同）。基是列下标的列表 `basis`，`basis[i]` 是第 $i$ 行的基变量。代码中的下标从 0 开始，因此 $x_1$ 是第 0 列。

```python
def print_tableau(T, basis, names, labels=("-z",)) -> None:
    """Print a tableau in the book's layout (objective rows last).

    Parameters
    ----------
    T : list of list of Fraction
        Constraint rows followed by len(labels) objective rows; the last
        column is the right-hand side.
    basis : list of int
        Column index of the basic variable of each constraint row.
    names : list of str
        Column names, e.g. ["x1", "x2", "x3"].
    labels : sequence of str
        Labels of the objective rows.
    """
    m = len(basis)
    if len(T) != m + len(labels):
        raise ValueError("T must have len(basis) + len(labels) rows")
    row_names = [names[j] for j in basis] + list(labels)
    header = ["Basic"] + list(names) + ["RHS"]
    cells = [[name] + [str(v) for v in row] for name, row in zip(row_names, T)]
    width = max(len(c) for row in [header] + cells for c in row)
    for row in [header] + cells:
        print(" ".join(c.rjust(width) for c in row))
```

(sec-09-03-basic-solutions)=
## 枚举基本解

由线性规划基本定理（fundamental theorem of linear programming，{prf:ref}`thm-03-02-fundamental`），标准形（standard form）的线性规划问题（linear program, LP）若有最优解，则必有最优基本可行解（basic feasible solution, BFS），因此对小规模问题，只需尝试从 $n$ 列中选取 $m$ 列的全部 $\binom nm$ 种选法：跳过奇异的选法，求解 $\mB\vx_B=\vb$，保留满足 $\vx_B\ge\vzero$ 的解（见{prf:ref}`def-03-02-basic-solution`）。

```python
def basic_solutions(A, b):
    """All basic solutions of A x = b (A of full row rank m).

    Parameters
    ----------
    A : list of list of Fraction, shape (m, n)
    b : list of Fraction, length m

    Returns
    -------
    list of tuple
        (basis, x) for every choice of m columns, in lexicographic order;
        x is None when the columns are linearly dependent.
    """
    m, n = len(A), len(A[0])
    if len(b) != m or m > n:
        raise ValueError("need len(b) == m <= n")
    result = []
    for basis in itertools.combinations(range(n), m):
        B = [[A[i][j] for j in basis] for i in range(m)]
        try:
            x_B = mat_vec(gauss_jordan_inverse(B), b)
        except ValueError:
            result.append((basis, None))
            continue
        x = [Fraction(0)] * n
        for j, value in zip(basis, x_B):
            x[j] = value
        result.append((basis, x))
    return result


# The system of ex-03-02-basic-solutions.
A_small = to_fractions(
    [
        [2, 3, 1],
        [1, 2, -1],
    ]
)
b_small = [Fraction(5), Fraction(3)]
for basis, x in basic_solutions(A_small, b_small):
    feasible = all(v >= 0 for v in x)
    print([f"x{j + 1}" for j in basis], [str(v) for v in x], feasible)
# Expected:
# ['x1', 'x2'] ['1', '1', '0'] True
# ['x1', 'x3'] ['8/3', '0', '-1/3'] False
# ['x2', 'x3'] ['0', '8/5', '1/5'] True
```

对于标准形的产品组合问题{eq}`eq-03-02-product-mix-standard`，松弛变量为 $x_3,x_4,x_5$，由此得到{numref}`tab-03-02-product-mix-basic`中的十个候选基，其中最好的可行基就是最优基：

```python
# Product mix in standard form (slacks x3, x4, x5), maximize c^T x.
A_pm = to_fractions(
    [
        [6, 5, 1, 0, 0],
        [10, 20, 0, 1, 0],
        [1, 0, 0, 0, 1],
    ]
)
b_pm = [Fraction(60), Fraction(150), Fraction(8)]
c_pm = [Fraction(500), Fraction(450), Fraction(0), Fraction(0), Fraction(0)]

best = None
for basis, x in basic_solutions(A_pm, b_pm):
    label = ",".join(f"x{j + 1}" for j in basis)
    if x is None:
        print(f"{label:<9} singular")
        continue
    z = sum(cj * xj for cj, xj in zip(c_pm, x))
    feasible = all(v >= 0 for v in x)
    print(f"{label:<9} {str([str(v) for v in x]):<40} {feasible!s:<6} {z}")
    if feasible and (best is None or z > best[0]):
        best = (z, x)
print("best:", best[0], [str(v) for v in best[1]])
# Expected:
# x1,x2,x3  ['8', '7/2', '-11/2', '0', '0']          False  5575
# x1,x2,x4  ['8', '12/5', '0', '22', '0']            True   5080
# x1,x2,x5  ['45/7', '30/7', '0', '0', '11/7']       True   36000/7
# x1,x3,x4  ['8', '0', '12', '70', '0']              True   4000
# x1,x3,x5  ['15', '0', '-30', '0', '-7']            False  7500
# x1,x4,x5  ['10', '0', '0', '50', '-2']             False  5000
# x2,x3,x4  singular
# x2,x3,x5  ['0', '15/2', '45/2', '0', '8']          True   3375
# x2,x4,x5  ['0', '12', '0', '-90', '8']             False  5400
# x3,x4,x5  ['0', '0', '60', '150', '8']             True   0
# best: 36000/7 ['45/7', '30/7', '0', '0', '11/7']
```

有五个基可行，一个是奇异的（$x_2,x_3,x_4$：$\mB$ 的第三行为零），最好的基本可行解就是最优解 $(\tfrac{45}{7},\tfrac{30}{7})$，$z^*=\tfrac{36000}{7}$。有些不可行的基本解的 $z$ 更大，但它们不是候选解。枚举的工作量按 $\binom nm$ 增长，对于一个很小的 LP 已经达到 $\binom{20}{10}=184\,756$，这正是单纯形法改为从一个基本可行解移动到更好的相邻基本可行解的原因。

(sec-09-03-tableau)=
## 表格单纯形法

### 主元变换、进基规则与出基规则

以 $\bar a_{rs}$ 为主元（pivot element）的**主元变换**就是高斯–若尔当消元的一步{eq}`eq-03-03-pivot`：将第 $r$ 行除以主元，再从其他每一行（包括目标行）中减去它的适当倍数。对任意一个基，从原始数据 $[\mA\mid\vb]$、$[\vc^\T\mid0]$ 出发，逐行对各基变量列作主元变换，就得到该基下典式（canonical form）的单纯形表；对于松弛变量构成的基，这些主元变换不改变任何东西。

```python
def pivot(T, r: int, s: int) -> None:
    """Pivot the tableau T in place on entry (r, s), eq-03-03-pivot."""
    p = T[r][s]
    if p == 0:
        raise ValueError(f"zero pivot element at ({r}, {s})")
    T[r] = [v / p for v in T[r]]
    for i in range(len(T)):
        if i != r and T[i][s] != 0:
            f = T[i][s]
            T[i] = [a - f * b for a, b in zip(T[i], T[r])]


def make_tableau(A, b, c, basis):
    """Canonical tableau of min/max c^T x, A x = b for a given basis.

    The objective row is [c | 0]; pivoting on (i, basis[i]) for every row
    makes the basic columns unit vectors and prices out the -z row.
    """
    if len(basis) != len(A):
        raise ValueError("one basic variable per row is needed")
    T = [list(row) + [bi] for row, bi in zip(A, b)]
    T.append(list(c) + [Fraction(0)])
    for i, j in enumerate(basis):
        pivot(T, i, j)
    return T
```

主元由{numref}`tab-03-03-rules`中的规则选取。进基列根据单纯形表的**最后**一行选取，在第二阶段是 $-z$ 行，在第一阶段是 $-\phi$ 行。Dantzig 规则对最小化问题选取最负的检验数（对最大化问题选取最正的检验数），并列时取下标最小者；Bland 规则（{prf:ref}`alg-03-03-bland`）选取检验数能改进目标的最小下标。最小比值检验{eq}`eq-03-03-ratio-test`在并列时取出基变量下标最小者，与{prf:ref}`alg-03-03-simplex`相同。

```python
def entering_column(T, basis, sense="min", rule="dantzig"):
    """Entering column by Dantzig's or Bland's rule, or None if optimal."""
    if sense not in ("min", "max") or rule not in ("dantzig", "bland"):
        raise ValueError("sense must be min/max, rule dantzig/bland")
    sign = 1 if sense == "min" else -1
    # Improving columns: c_bar_j < 0 (min) or c_bar_j > 0 (max).
    cost = T[-1][:-1]
    improving = [
        j for j in range(len(cost)) if j not in basis and sign * cost[j] < 0
    ]
    if not improving:
        return None
    if rule == "bland":
        return min(improving)
    return min(improving, key=lambda j: (sign * cost[j], j))


def leaving_row(T, basis, s):
    """Minimum-ratio row for column s (ties: smallest basic index).

    Returns None if no entry of column s is positive (unbounded).
    """
    m = len(basis)
    ratios = [
        (T[i][-1] / T[i][s], basis[i], i) for i in range(m) if T[i][s] > 0
    ]
    if not ratios:
        return None
    return min(ratios)[2]
```

### 算法

现在{prf:ref}`alg-03-03-simplex`只是一个简短的循环。它把每次主元变换记录为（进基变量，出基变量，$\theta$），并以求解状态 `"optimal"`、`"unbounded"` 终止，或在基重复出现时以 `"cycling"` 终止，而基重复只可能发生在退化顶点处（见{prf:ref}`thm-03-03-finite`）。

```python
def simplex(T, basis, sense="min", rule="dantzig", max_iter=100):
    """Tableau simplex method (alg-03-03-simplex), pricing on the last row.

    Parameters
    ----------
    T : list of list of Fraction
        Tableau in canonical form for `basis` (modified in place).
    basis : list of int
        Basic column of each constraint row (modified in place).
    sense : {"min", "max"}
    rule : {"dantzig", "bland"}
    max_iter : int

    Returns
    -------
    status : str
        "optimal", "unbounded", "cycling" or "iteration limit".
    history : list of tuple
        (entering column, leaving column, theta) of every pivot.
    """
    history = []
    seen = {frozenset(basis)}
    for _ in range(max_iter):
        s = entering_column(T, basis, sense, rule)
        if s is None:
            return "optimal", history
        r = leaving_row(T, basis, s)
        if r is None:
            return "unbounded", history
        history.append((s, basis[r], T[r][-1] / T[r][s]))
        pivot(T, r, s)
        basis[r] = s
        if frozenset(basis) in seen:
            return "cycling", history
        seen.add(frozenset(basis))
    return "iteration limit", history


def solution(T, basis, n):
    """Values of the n variables in the current basic solution."""
    x = [Fraction(0)] * n
    for i, j in enumerate(basis):
        x[j] = T[i][-1]
    return x
```

### 产品组合的计算过程

从松弛变量构成的基 $\{x_3,x_4,x_5\}$ 出发，对最大化问题使用 Dantzig 规则，重现了{prf:ref}`ex-03-03-simplex-run`中的三次主元变换：$x_1$ 替换 $x_5$（$\theta=8$），$x_2$ 替换 $x_3$（$\theta=\tfrac{12}{5}$），$x_5$ 替换 $x_4$（$\theta=\tfrac{11}{7}$），路线如{numref}`fig-03-03-simplex-path`所示。

```python
names_pm = ["x1", "x2", "x3", "x4", "x5"]
basis = [2, 3, 4]
T = make_tableau(A_pm, b_pm, c_pm, basis)
status, history = simplex(T, basis, sense="max")
for s, leave, theta in history:
    print(f"{names_pm[s]} enters, {names_pm[leave]} leaves, theta = {theta}")
print(status)
print_tableau(T, basis, names_pm)
x_opt = solution(T, basis, 5)
print("x* =", [str(v) for v in x_opt], " z* =", -T[-1][-1])
# Expected:
# x1 enters, x5 leaves, theta = 8
# x2 enters, x3 leaves, theta = 12/5
# x5 enters, x4 leaves, theta = 11/7
# optimal
#    Basic       x1       x2       x3       x4       x5      RHS
#       x2        0        1     -1/7     3/35        0     30/7
#       x5        0        0     -2/7     1/14        1     11/7
#       x1        1        0      2/7    -1/14        0     45/7
#       -z        0        0   -550/7    -20/7        0 -36000/7
# x* = ['45/7', '30/7', '0', '0', '11/7']  z* = 36000/7
```

最后一张表就是{numref}`tab-03-03-pm-3`，包括第 $x_2$ 行、第 $x_3$ 列的元素 $-\tfrac17$；影子价格 $\tfrac{550}{7}$ 和 $\tfrac{20}{7}$（取负号后）出现在 $-z$ 行的松弛变量各列下方。求解器的结果与之一致：

```python
res = linprog(
    [-500, -450],
    A_ub=[[float(v) for v in row[:2]] for row in A_pm],
    b_ub=[float(v) for v in b_pm],
    method="highs",
)
print(np.isclose(-res.fun, float(-T[-1][-1])))
print(np.allclose(res.x, [float(v) for v in x_opt[:2]]))
# Expected:
# True
# True
```

### 循环与 Bland 规则

Beale 的例子（{prf:ref}`ex-03-03-cycling`）对于基 $\{x_1,x_2,x_3\}$ 已经是典式。在 Dantzig 规则下，每次主元变换都是退化的（$\theta=0$），六次主元变换之后又回到初始基。

```python
A_beale = to_fractions(
    [
        [1, 0, 0, "1/4", -8, -1, 9],
        [0, 1, 0, "1/2", -12, "-1/2", 3],
        [0, 0, 1, 0, 0, 1, 0],
    ]
)
b_beale = [Fraction(0), Fraction(0), Fraction(1)]
c_beale = [Fraction(v) for v in ["0", "0", "0", "-3/4", "20", "-1/2", "6"]]
names_beale = [f"x{j + 1}" for j in range(7)]

for rule in ["dantzig", "bland"]:
    basis = [0, 1, 2]
    T = make_tableau(A_beale, b_beale, c_beale, basis)
    status, history = simplex(T, basis, sense="min", rule=rule)
    pivots = [
        f"{names_beale[s]}/{names_beale[leave]}({theta})"
        for s, leave, theta in history
    ]
    print(rule, status, "z =", -T[-1][-1])
    print("  pivots (enters/leaves, theta):")
    print("   ", " ".join(pivots))
x_bland = solution(T, basis, 7)
print("x* =", [str(v) for v in x_bland])
# Expected:
# dantzig cycling z = 0
#   pivots (enters/leaves, theta):
#     x4/x1(0) x5/x2(0) x6/x4(0) x7/x5(0) x1/x6(0) x2/x7(0)
# bland optimal z = -5/4
#   pivots (enters/leaves, theta):
#     x4/x1(0) x5/x2(0) x6/x4(0) x1/x5(0) x2/x3(1/2) x4/x2(1)
# x* = ['3/4', '0', '0', '1', '0', '1', '0']
```

Dantzig 规则的主元序列就是{numref}`tab-03-03-cycling`中的序列；Bland 规则前三次主元变换与之相同，随后选取 $x_1$ 而不是 $x_7$ 进基，再经过三次主元变换到达最优解 $\vx^*=(\tfrac34,0,0,1,0,1,0)$，$z^*=-\tfrac54$，与{ref}`sec-03-03-degeneracy`相同。HiGHS 确认了这一最优解（它自己的防循环手段是摄动，而不是 Bland 规则）：

```python
res = linprog(
    [float(v) for v in c_beale],
    A_eq=[[float(v) for v in row] for row in A_beale],
    b_eq=[float(v) for v in b_beale],
    method="highs",
)
print(res.fun, res.x)
# Expected:
# -1.25 [0.75 0.   0.   1.   0.   1.   0.  ]
```

(sec-09-03-two-phase)=
## 两阶段法

当没有现成的松弛变量基时，{prf:ref}`alg-03-04-two-phase`给需要的行加上人工变量 $R_i$，并在第一阶段最小化它们的和 $\phi$。下面的函数（在使 $\vb\ge\vzero$ 之后）给每一行都加上人工变量，这正是本例所需要的，而且总是可行；实用程序会尽量利用已有的松弛变量列。第一阶段中 $-z$ 行一同参与计算，其下方的 $-\phi$ 行{eq}`eq-03-04-phase-one-row`决定主元变换。第一阶段结束后，仍以零值留在基中的人工变量，通过以原变量列中任意一个非零元素为主元移出基；若该行没有这样的元素，则该行冗余，予以删除（见{prf:ref}`thm-03-04-phase-one`）。

```python
def two_phase(A, b, c, sense="min", rule="dantzig", verbose=False):
    """Two-phase simplex method (alg-03-04-two-phase) with exact arithmetic.

    Parameters
    ----------
    A, b, c : lists of Fraction
        Standard-form data, A x = b, x >= 0.
    sense : {"min", "max"}
    rule : {"dantzig", "bland"}
    verbose : bool
        Print the tableau at the end of Phase I and after the clean-up.

    Returns
    -------
    status : str
        "optimal", "infeasible" or "unbounded" ("cycling" or "iteration
        limit" if a phase stops early).
    x : list of Fraction or None
        Optimal solution (None unless optimal).
    z : Fraction or None
        Optimal value.
    """
    m, n = len(A), len(A[0])
    names = [f"x{j + 1}" for j in range(n)] + [f"R{i + 1}" for i in range(m)]
    # Rows with b_i < 0 are multiplied by -1 so that b >= 0.
    rows = [
        (list(A[i]), b[i]) if b[i] >= 0 else ([-a for a in A[i]], -b[i])
        for i in range(m)
    ]
    T = [
        row + [Fraction(int(k == i)) for k in range(m)] + [rhs]
        for i, (row, rhs) in enumerate(rows)
    ]
    # Objective rows: -z (zero cost for the artificials), then -phi.
    T.append(list(c) + [Fraction(0)] * (m + 1))
    T.append([Fraction(0)] * n + [Fraction(1)] * m + [Fraction(0)])
    basis = list(range(n, n + m))
    # Price out the artificial columns of the -phi row.
    for i in range(m):
        T[-1] = [a - r for a, r in zip(T[-1], T[i])]

    # Phase I: minimize phi, pricing on the -phi row.
    status, _ = simplex(T, basis, "min", rule)
    if verbose:
        print("end of Phase I:")
        print_tableau(T, basis, names, labels=("-z", "-phi"))
    if status != "optimal":
        # Phase I stopped early: phi > 0 would not prove infeasibility.
        return status, None, None
    if -T[-1][-1] > 0:
        return "infeasible", None, None

    # Clean-up: drive basic artificial variables out of the basis.
    for i in reversed(range(len(basis))):
        if basis[i] >= n:
            s = next((j for j in range(n) if T[i][j] != 0), None)
            # No nonzero entry in an original column: the row is redundant.
            if s is None:
                del T[i], basis[i]
            else:
                pivot(T, i, s)
                basis[i] = s
    if verbose:
        print("after the clean-up:")
        print_tableau(T, basis, names, labels=("-z", "-phi"))

    # Phase II: delete the artificial columns and the -phi row.
    T = [row[:n] + [row[-1]] for row in T[:-1]]
    status, _ = simplex(T, basis, sense, rule)
    if status != "optimal":
        return status, None, None
    return status, solution(T, basis, n), -T[-1][-1]
```

{prf:ref}`ex-03-04-two-phase`在三个等式约束下最大化 $-3x_1+x_3$。第一阶段需要两次主元变换，结束时 $R_3$ 以零值留在基中；清理步骤以 $x_4$ 下方的元素 $-1$ 为主元将其移出基；第二阶段需要一次主元变换。

```python
A_tp = to_fractions(
    [
        [1, 1, 1, 1],
        [-2, 1, -1, 0],
        [0, 3, 1, 1],
    ]
)
b_tp = [Fraction(4), Fraction(1), Fraction(9)]
c_tp = [Fraction(-3), Fraction(0), Fraction(1), Fraction(0)]
status, x_tp, z_tp = two_phase(A_tp, b_tp, c_tp, sense="max", verbose=True)
print(status, [str(v) for v in x_tp], z_tp)
# Expected:
# end of Phase I:
# Basic    x1    x2    x3    x4    R1    R2    R3   RHS
#    x1     1     0   2/3   1/3   1/3  -1/3     0     1
#    x2     0     1   1/3   2/3   2/3   1/3     0     3
#    R3     0     0     0    -1    -2    -1     1     0
#    -z     0     0     3     1     1    -1     0     3
#  -phi     0     0     0     1     3     2     0     0
# after the clean-up:
# Basic    x1    x2    x3    x4    R1    R2    R3   RHS
#    x1     1     0   2/3     0  -1/3  -2/3   1/3     1
#    x2     0     1   1/3     0  -2/3  -1/3   2/3     3
#    x4     0     0     0     1     2     1    -1     0
#    -z     0     0     3     0    -1    -2     1     3
#  -phi     0     0     0     0     1     1     1     0
# optimal ['0', '5/2', '3/2', '0'] 3/2
```

打印出的两张单纯形表都与{prf:ref}`ex-03-04-two-phase`中的相同，最优解 $\vx^*=(0,\tfrac52,\tfrac32,0)$，$z^*=\tfrac32$ 与 HiGHS 一致。不可行的 LP 在第一阶段结束时 $\phi^*>0$；{doc}`03-07-sensitivity-analysis`中 $b_1=15$ 的贯穿示例（running example）就是这样一个例子（见{prf:ref}`ex-03-07-rhs`）：

```python
res = linprog(
    [3, 0, -1, 0],
    A_eq=[[float(v) for v in row] for row in A_tp],
    b_eq=[4, 1, 9],
    method="highs",
)
print(round(-res.fun, 6), res.x)

A37 = to_fractions(
    [
        [1, 2, 0, 1, 0, -6],
        [0, 1, 1, 3, -2, -1],
        [1, 2, 1, 3, -1, -5],
    ]
)
c37 = [Fraction(v) for v in [3, 2, -3, -6, 10, -5]]
b37 = [Fraction(11), Fraction(6), Fraction(13)]
print(two_phase(A37, [Fraction(15), b37[1], b37[2]], c37)[0])
status, x37, z37 = two_phase(A37, b37, c37)
print(status, [str(v) for v in x37], z37)
# Expected:
# 1.5 [0.  2.5 1.5 0. ]
# infeasible
# optimal ['3', '4', '2', '0', '0', '0'] 11
```

当 $b_1=11$ 时，同一函数从头求得{prf:ref}`ex-03-07-running-lp`的最优解 $\vx^*=(3,4,2,0,0,0)$，$z^*=11$。

(sec-09-03-revised)=
## 修正单纯形法

{prf:ref}`alg-03-06-revised-simplex`从不构造单纯形表。它保存原始数据、当前的基变量值 $\vx_B$，以及**乘积形式**{eq}`eq-03-06-product-form`的基逆矩阵，即一列 eta 矩阵 $\mE_k$，每个矩阵存储为其主元位置 $r$ 和 eta 列（eta column）$\veta$ {eq}`eq-03-06-eta`。把 $\mE$ 作用于一个列向量，会使每个元素改变第 $r$ 个元素的某个倍数（{eq}`eq-03-06-eta-action`）；把它从右边作用于一个行向量，只改变第 $r$ 个元素，该元素变为 $\vu^\T\veta$。

```python
def eta_column(a_bar, r):
    """Eta column of the pivot on entry r of the updated column a_bar."""
    if a_bar[r] == 0:
        raise ValueError("zero pivot element")
    return [
        1 / a_bar[r] if i == r else -a_bar[i] / a_bar[r]
        for i in range(len(a_bar))
    ]


def apply_eta(r, eta, v):
    """E v for the eta matrix (r, eta), by eq-03-06-eta-action."""
    return [
        eta[i] * v[r] if i == r else v[i] + eta[i] * v[r]
        for i in range(len(v))
    ]


def eta_matrix(r, eta):
    """The eta matrix as a full matrix (identity with column r = eta)."""
    m = len(eta)
    return [
        [eta[i] if j == r else Fraction(int(i == j)) for j in range(m)]
        for i in range(m)
    ]
```

算法从松弛变量构成的基出发，其基矩阵为单位矩阵，因此 eta 列表开始时为空。定价时把各 eta 矩阵从最后一个到第一个依次作用于 $\vc_B^\T$；进基列则从第一个到最后一个依次更新。进基规则采用 Dantzig 规则，但可以用列下标的列表 `first_entering` 规定最初的几次选择，以便重演手算过程。

```python
def revised_simplex(A, b, c, basis, first_entering=(), max_iter=50):
    """Revised simplex method in product form (alg-03-06-revised-simplex).

    Minimization; the starting basis must be a slack basis (B = I).

    Returns
    -------
    status : str
        "optimal" or "unbounded".
    basis : list of int
    x_B : list of Fraction
        Values of the basic variables.
    pi : list of Fraction
        Simplex multipliers at the end.
    etas : list of tuple
        The eta factors (r, eta) in the order they were created.
    """
    m, n = len(A), len(A[0])
    if [[A[i][j] for j in basis] for i in range(m)] != [
        [Fraction(int(i == k)) for k in range(m)] for i in range(m)
    ]:
        raise ValueError("the starting basis must be the identity matrix")
    basis, x_B, etas = list(basis), list(b), []
    for it in range(max_iter):
        # Step 1, price: pi^T = c_B^T E_k ... E_1 (etas applied last first).
        pi = [c[j] for j in basis]
        for r, eta in reversed(etas):
            pi[r] = sum((u * e for u, e in zip(pi, eta)), Fraction(0))
        reduced = {
            j: c[j] - sum((pi[i] * A[i][j] for i in range(m)), Fraction(0))
            for j in range(n)
            if j not in basis
        }
        # Step 2, entering variable.
        if it < len(first_entering):
            s = first_entering[it]
        else:
            negative = [j for j, cj in reduced.items() if cj < 0]
            if not negative:
                return "optimal", basis, x_B, pi, etas
            s = min(negative, key=lambda j: (reduced[j], j))
        a_bar = [A[i][s] for i in range(m)]
        for r, eta in etas:
            a_bar = apply_eta(r, eta, a_bar)
        # Step 3, ratio test.
        ratios = [
            (x_B[i] / a_bar[i], basis[i], i) for i in range(m) if a_bar[i] > 0
        ]
        if not ratios:
            return "unbounded", basis, x_B, pi, etas
        r = min(ratios)[2]
        # Step 4, update: new eta factor, x_B <- E x_B, basis change.
        eta = eta_column(a_bar, r)
        etas.append((r, eta))
        x_B = apply_eta(r, eta, x_B)
        basis[r] = s
    raise RuntimeError("iteration limit reached")
```

{prf:ref}`ex-03-06-revised-run`先让 $x_2$ 进基，随后遵循 Dantzig 规则。三个 eta 矩阵和最终的逆矩阵 $\mB_4^{-1}=\mE_3\mE_2\mE_1$ 都与该例相同：

```python
A_rev = to_fractions(
    [
        [2, 1, 1, 1, 0, 0],
        [1, 2, 3, 0, 1, 0],
        [2, 2, 1, 0, 0, 1],
    ]
)
b_rev = [Fraction(2), Fraction(5), Fraction(6)]
c_rev = [Fraction(v) for v in [-3, -1, -3, 0, 0, 0]]

status, basis, x_B, pi, etas = revised_simplex(
    A_rev, b_rev, c_rev, [3, 4, 5], first_entering=[1]
)
for k, (r, eta) in enumerate(etas, start=1):
    print(f"E{k}:")
    show(eta_matrix(r, eta))
# B^{-1} = E3 E2 E1: apply the etas to the columns of the identity.
columns = []
for j in range(3):
    e_j = [Fraction(int(i == j)) for i in range(3)]
    for r, eta in etas:
        e_j = apply_eta(r, eta, e_j)
    columns.append(e_j)
print("B4^-1:")
show([list(row) for row in zip(*columns)])
print(status, [f"x{j + 1}" for j in basis], [str(v) for v in x_B])
print("pi =", [str(v) for v in pi])
print("z* =", sum(c_rev[j] * v for j, v in zip(basis, x_B)))
# Expected:
# E1:
#  1   0   0
# -2   1   0
# -2   0   1
# E2:
#  1  -1   0
#  0   1   0
#  0   1   1
# E3:
# 1/5    0    0
# 3/5    1    0
#   1    0    1
# B4^-1:
#  3/5  -1/5     0
# -1/5   2/5     0
#   -1     0     1
# optimal ['x1', 'x3', 'x6'] ['1/5', '8/5', '4']
# pi = ['-6/5', '-3/5', '0']
# z* = -27/5
```

最优解 $\vx^*=(\tfrac15,0,\tfrac85,0,0,4)$，$z^*=-\tfrac{27}{5}$ 以及单纯形乘子 $\vpi=(-\tfrac65,-\tfrac35,0)$ 与该例一致，也满足 $z^*=\vpi^\T\vb$。若从一开始就用 Dantzig 规则，则两次主元变换就能到达同一个基，正如该例之后所指出的：

```python
status, basis, x_B, pi, etas = revised_simplex(A_rev, b_rev, c_rev, [3, 4, 5])
print(status, len(etas), [f"x{j + 1}" for j in basis], [str(v) for v in x_B])
res = linprog(
    [float(v) for v in c_rev],
    A_eq=[[float(v) for v in row] for row in A_rev],
    b_eq=[2, 5, 6],
    method="highs",
)
print(round(res.fun, 6), res.x)
# Expected:
# optimal 2 ['x1', 'x3', 'x6'] ['1/5', '8/5', '4']
# -5.4 [0.2 0.  1.6 0.  0.  4. ]
```

(sec-09-03-dual-simplex)=
## 对偶单纯形法

{prf:ref}`alg-03-06-dual-simplex`作用于这样的单纯形表：检验数全部非负（对偶可行，最小化），但右端项有负元素。它先选取出基行（最负的 $\bar b_r$）；若该行没有负元素，则以不可行终止；否则让对偶比值 $\bar c_j/(-\bar a_{rj})$（$\bar a_{rj}<0$）最小的列进基（{eq}`eq-03-06-dual-ratio`）。

```python
def dual_simplex(T, basis, max_iter=50):
    """Dual simplex method (alg-03-06-dual-simplex) on a minimization tableau.

    T must be dual feasible: every reduced cost in the last row is >= 0.
    T and basis are modified in place.

    Returns
    -------
    status : str
        "optimal" or "infeasible".
    history : list of tuple
        (leaving column, entering column, dual ratio) of every pivot.
    """
    m, n = len(basis), len(T[0]) - 1
    if any(v < 0 for v in T[-1][:n]):
        raise ValueError("the tableau is not dual feasible")
    history = []
    for _ in range(max_iter):
        r = min(range(m), key=lambda i: (T[i][-1], basis[i]))
        if T[r][-1] >= 0:
            return "optimal", history
        candidates = [j for j in range(n) if T[r][j] < 0]
        if not candidates:
            return "infeasible", history
        s = min(candidates, key=lambda j: (T[-1][j] / -T[r][j], j))
        history.append((basis[r], s, T[-1][s] / -T[r][s]))
        pivot(T, r, s)
        basis[r] = s
    raise RuntimeError("iteration limit reached")
```

{prf:ref}`ex-03-06-dual-simplex`直接给出了单纯形表。一次主元变换（$x_2$ 出基、$x_4$ 进基，对偶比值为 $\tfrac12$）就使其达到最优，$z^*=1$：

```python
names6 = [f"x{j + 1}" for j in range(6)]
T = to_fractions(
    [
        [1, 0, 0, 4, -5, 7, 8],
        [0, 1, 0, -2, 4, -2, -2],
        [0, 0, 1, 1, -3, 2, 2],
        [0, 0, 0, 1, 3, 2, 0],
    ]
)
basis = [0, 1, 2]
status, history = dual_simplex(T, basis)
print(status, [(names6[l], names6[e], str(t)) for l, e, t in history])
print_tableau(T, basis, names6)
# Expected:
# optimal [('x2', 'x4', '1/2')]
# Basic    x1    x2    x3    x4    x5    x6   RHS
#    x1     1     2     0     0     3     3     4
#    x4     0  -1/2     0     1    -2     1     1
#    x3     0   1/2     1     0    -1     1     1
#    -z     0   1/2     0     0     5     1    -1
```

对偶单纯形法是 $\vb$ 改变或新增约束之后重新优化（re-optimization）的工具（见{ref}`sec-03-07-reoptimisation`）。对于{prf:ref}`ex-03-07-rhs`，$b_1=15$：原最优基 $(x_1,x_2,x_3)$ 仍然对偶可行，有两个基变量的值变为负数，该方法经过一次主元变换就证明了问题不可行。对于{prf:ref}`ex-03-07-new-constraint`，添加带新松弛变量 $x_7$ 的行 $x_1+x_2+x_7=6$，减去各基变量所在行使之化为典式，一次主元变换就恢复了可行性。

```python
# (b) b1 = 15: optimal basis of the running example, new right-hand side.
T = make_tableau(A37, [Fraction(15), b37[1], b37[2]], c37, [0, 1, 2])
basis = [0, 1, 2]
status, history = dual_simplex(T, basis)
print(status, [(names6[l], names6[e], str(t)) for l, e, t in history])
print("row of x1:", [str(v) for v in T[basis.index(0)]])

# New constraint x1 + x2 <= 6 with slack x7.
T = make_tableau(A37, b37, c37, [0, 1, 2])
# Insert a zero column for x7 before the right-hand side.
T = [row[:-1] + [Fraction(0)] + row[-1:] for row in T]
new_row = [Fraction(v) for v in [1, 1, 0, 0, 0, 0, 1, 6]]
# Canonical form: subtract multiples of the rows of x1, x2 (and x3).
for i, j in enumerate([0, 1, 2]):
    new_row = [a - new_row[j] * t for a, t in zip(new_row, T[i])]
T.insert(3, new_row)
basis = [0, 1, 2, 6]
print("new row:", [str(v) for v in new_row])
status, history = dual_simplex(T, basis)
names7 = names6 + ["x7"]
print(status, [(names7[l], names7[e], str(t)) for l, e, t in history])
x_new = solution(T, basis, 7)
print("x* =", [str(v) for v in x_new[:6]], " z* =", -T[-1][-1])
# Expected:
# infeasible [('x3', 'x5', '3')]
# row of x1: ['1', '0', '2', '3', '0', '0', '-5']
# new row: ['0', '0', '0', '0', '-1', '4', '1', '-1']
# optimal [('x7', 'x5', '3')]
# x* = ['1', '5', '3', '0', '1', '0']  z* = 14
```

两次运行都与{doc}`03-07-sensitivity-analysis`一致：当 $b_1=15$ 时，行 $x_1+2x_3+3x_4=-5$ 没有负元素，证明 LP 不可行；加入新约束后，最优解移到 $\vx^*=(1,5,3,0,1,0)$，$z^*=14$。HiGHS 的结果一致：

```python
A37f = np.array([[float(v) for v in row] for row in A37])
c37f = np.array([float(v) for v in c37])
res_b = linprog(c37f, A_eq=A37f, b_eq=[15, 6, 13], method="highs")
res_c = linprog(
    c37f,
    A_ub=[[1, 1, 0, 0, 0, 0]],
    b_ub=[6],
    A_eq=A37f,
    b_eq=[11, 6, 13],
    method="highs",
)
print(res_b.status, round(res_c.fun, 6), res_c.x)
# Expected:
# 2 14.0 [1. 5. 3. 0. 1. 0.]
```

(sec-09-03-sensitivity)=
## 灵敏度变化范围分析

给定最优基，{doc}`03-07-sensitivity-analysis`中的变化范围只需要 $\mB^{-1}$、最优单纯形表和检验数：

- **非基变量的费用** $c_r$ 可以无限增大，最多可减小到盈亏平衡成本（break-even cost）$\vpi^\T\mA_{\cdot r}=c_r-\bar c_r$（{eq}`eq-03-07-nonbasic-range`）；
- 位于第 $i$ 行的**基变量的费用** $c_r$ 的改变量 $\delta$ 可以在区间{eq}`eq-03-07-basic-range`内变化，该区间由各非基变量列上的比值 $\bar c_j/\bar a_{ij}$ 确定；
- **右端项** $b_r$ 的改变量 $\delta$ 只要满足 $\vx_B+\delta\,\mB^{-1}\ve_r\ge\vzero$ 即可（{eq}`eq-03-07-rhs-range`）。

下面的函数对最小化问题返回全部三类变化范围，无穷端点用 `None` 表示。

```python
def sensitivity_ranges(A, b, c, basis):
    """Cost and right-hand-side ranges of an optimal basis (minimization).

    Parameters
    ----------
    A, b, c : lists of Fraction
        Standard-form data.
    basis : list of int
        An optimal basis (basis[i] is basic in row i).

    Returns
    -------
    cost_ranges : list of tuple
        (low, high) for every c_j; None for an infinite end.
    rhs_ranges : list of tuple
        (low, high) for every b_i.
    """
    m, n = len(A), len(A[0])
    B_inv = gauss_jordan_inverse([[A[i][j] for j in basis] for i in range(m)])
    x_B = mat_vec(B_inv, b)
    pi = [
        sum((c[basis[k]] * B_inv[k][i] for k in range(m)), Fraction(0))
        for i in range(m)
    ]
    A_bar = [mat_vec(B_inv, [A[i][j] for i in range(m)]) for j in range(n)]
    c_bar = [
        c[j] - sum((pi[i] * A[i][j] for i in range(m)), Fraction(0))
        for j in range(n)
    ]
    if min(x_B) < 0 or min(c_bar) < 0:
        raise ValueError("the basis is not optimal")
    nonbasic = [j for j in range(n) if j not in basis]

    cost_ranges = []
    for j in range(n):
        # Nonbasic cost: eq-03-07-nonbasic-range.
        if j not in basis:
            cost_ranges.append((c[j] - c_bar[j], None))
            continue
        # Basic cost: eq-03-07-basic-range with row i of the tableau.
        i = basis.index(j)
        up = [c_bar[k] / A_bar[k][i] for k in nonbasic if A_bar[k][i] > 0]
        down = [c_bar[k] / A_bar[k][i] for k in nonbasic if A_bar[k][i] < 0]
        cost_ranges.append(
            (
                c[j] + max(down) if down else None,
                c[j] + min(up) if up else None,
            )
        )

    rhs_ranges = []
    # Right-hand sides: eq-03-07-rhs-range with column r of B^{-1}.
    for r in range(m):
        col = [B_inv[i][r] for i in range(m)]
        down = [-x_B[i] / col[i] for i in range(m) if col[i] > 0]
        up = [-x_B[i] / col[i] for i in range(m) if col[i] < 0]
        rhs_ranges.append(
            (
                b[r] + max(down) if down else None,
                b[r] + min(up) if up else None,
            )
        )
    return cost_ranges, rhs_ranges


def fmt_range(low, high):
    """Interval as text, with infinite ends."""
    lo = "-inf" if low is None else str(low)
    hi = "inf" if high is None else str(high)
    return f"[{lo}, {hi}]"


cost_ranges, rhs_ranges = sensitivity_ranges(A37, b37, c37, [0, 1, 2])
for j, (low, high) in enumerate(cost_ranges):
    print(f"c{j + 1} = {c37[j]}: {fmt_range(low, high)}")
for i, (low, high) in enumerate(rhs_ranges):
    print(f"b{i + 1} = {b37[i]}: {fmt_range(low, high)}")
# Expected:
# c1 = 3: [2, 9/2]
# c2 = 2: [-1, 3]
# c3 = -3: [-6, -5/2]
# c4 = -6: [-7, inf]
# c5 = 10: [7, inf]
# c6 = -5: [-13, inf]
# b1 = 11: [7, 13]
# b2 = 6: [2, 15/2]
# b3 = 13: [23/2, 17]
```

这些就是{numref}`tab-03-07-ranges`中的变化范围。用求解器检验可以使其含义更加具体：费用刚好在范围之内时，最优解保持为 $(3,4,2,0,0,0)$，刚好在范围之外时最优解发生移动；右端项刚好在范围之内时，最优基（取正值的变量集合）保持为 $\{x_1,x_2,x_3\}$，刚好在范围之外时最优基改变或 LP 变为不可行。

```python
def optimum(c_vec, b_vec):
    """HiGHS optimum of min c^T x, A37 x = b, x >= 0 (None if not optimal)."""
    res = linprog(c_vec, A_eq=A37f, b_eq=b_vec, method="highs")
    return res.x if res.status == 0 else None


x_star = np.array([3.0, 4.0, 2.0, 0.0, 0.0, 0.0])
b37f = np.array([float(v) for v in b37])
eps = 0.01
for j, (low, high) in enumerate(cost_ranges):
    checks = []
    for end, inward in [(low, eps), (high, -eps)]:
        if end is None:
            continue
        for shift, expect_same in [(inward, True), (-inward, False)]:
            c_test = c37f.copy()
            c_test[j] = float(end) + shift
            x_test = optimum(c_test, b37f)
            same = x_test is not None and np.allclose(x_test, x_star)
            checks.append(same == expect_same)
    print(f"c{j + 1}: range confirmed {all(checks)}")
for i, (low, high) in enumerate(rhs_ranges):
    checks = []
    for end, inward in [(low, eps), (high, -eps)]:
        for shift, expect_same in [(inward, True), (-inward, False)]:
            b_test = b37f.copy()
            b_test[i] = float(end) + shift
            x_test = optimum(c37f, b_test)
            support = (
                None if x_test is None else set(np.flatnonzero(x_test > 1e-9))
            )
            checks.append((support == {0, 1, 2}) == expect_same)
    print(f"b{i + 1}: range confirmed {all(checks)}")
# Expected:
# c1: range confirmed True
# c2: range confirmed True
# c3: range confirmed True
# c4: range confirmed True
# c5: range confirmed True
# c6: range confirmed True
# b1: range confirmed True
# b2: range confirmed True
# b3: range confirmed True
```

在右端项的变化范围内，最优值按 $\pi_i$ 的速率变化；再加上{eq}`eq-03-07-T-running`中的 $\mT$ 矩阵，即把任一数据列变为最优单纯形表中相应列的加边逆矩阵（{eq}`eq-03-06-T-column`），这就是一份灵敏度报告的全部内容。$\mT$ 矩阵只需一次精确求逆即可得到：

```python
bordered = [[A37[i][j] for j in [0, 1, 2]] + [Fraction(0)] for i in range(3)]
bordered.append([c37[j] for j in [0, 1, 2]] + [Fraction(1)])
T_matrix = gauss_jordan_inverse(bordered)
show(T_matrix)
# Column of the new activity of ex-03-07-new-activity (b) with its cost.
print([str(v) for v in mat_vec(T_matrix, [3, -1, 1, 4])])
# Expected:
# -1  -2   2   0
#  1   1  -1   0
# -1   0   1   0
# -2   4  -1   1
# ['1', '1', '-2', '-7']
```

$\mT$ 的最后一行是 $(-\vpi^\T,1)=(-2,4,-1,1)$，{prf:ref}`ex-03-07-new-activity` (b) 中的新增活动的检验数为 $-7$，更新后的列为 $(1,1,-2)^\T$，与该例中的计算一致。

## 本章回顾

### 小结

- 精确的 `Fraction` 运算可以逐位重现本书中的单纯形表（见{ref}`sec-09-03-exact`）。
- 枚举 $\binom nm$ 个基可以求得小规模 LP 的最优 BFS；产品组合问题的十个基中有五个可行（见{ref}`sec-09-03-basic-solutions`）。
- 表格单纯形法就是主元变换{eq}`eq-03-03-pivot`加上进基规则与出基规则；Dantzig 规则三次主元变换解出产品组合问题，在 Beale 的例子上发生循环，而 Bland 规则到达 $z^*=-\tfrac54$（见{ref}`sec-09-03-tableau`）。
- 两阶段法添加人工变量，最小化它们的和，把退化的人工变量移出基，然后用一同更新的 $-z$ 行继续计算（见{ref}`sec-09-03-two-phase`）。
- 修正单纯形法只保存 $\vx_B$ 和 $\mB^{-1}$ 的 eta 因子（见{ref}`sec-09-03-revised`）；对偶单纯形法保持对偶可行性，在 $\vb$ 改变或新增约束后重新优化（见{ref}`sec-09-03-dual-simplex`）。
- 费用与右端项的变化范围以及 $\mT$ 矩阵都由最优基得到，并与用 HiGHS 重新求解的结果一致（见{ref}`sec-09-03-sensitivity`）。

### 常见错误

- **在错误的行上定价。** 第一阶段的进基变量由 $-\phi$ 行决定，而不是 $-z$ 行；每次主元变换都要更新这两行。
- **忽略优化方向。** 最大化问题的单纯形表在所有 $\bar c_j\le0$ 时最优；使用最小化的判别准则会使方法反向运行。
- **浮点数下的并列。** 使用浮点数时，像 $3/3$ 与 $6/6$ 这样的退化并列可能被舍入噪声打破；需要使用精确分数，或在每次比较中使用容差。
- **忘记第一阶段后以零值留在基中的人工变量**，或删除了并不冗余的行。
- **以错误的顺序应用 eta 矩阵。** 更新列向量时先用 $\mE_1$，更新行向量（$\vpi^\T$）时先用最后一个 eta 矩阵。
- **从非对偶可行的单纯形表出发使用对偶单纯形法**；本章的函数会检查这一点并拒绝执行。
- **从最优单纯形表中基变量所在的列而不是所在的行读取其费用变化范围。**

## 参考文献

- Bertsimas, D., & Tsitsiklis, J. N. (1997). *Introduction to linear
  optimization*. Athena Scientific. (Chapter 3: the simplex method, the full
  tableau and revised implementations, anticycling, Phase I; Chapter 4: the
  dual simplex method; Chapter 5: sensitivity analysis.)
- Bland, R. G. (1977). New finite pivoting rules for the simplex method.
  *Mathematics of Operations Research, 2*(2), 103–107.
  <https://doi.org/10.1287/moor.2.2.103>
- Dantzig, G. B. (1963). *Linear programming and extensions*. Princeton
  University Press. (Chapter 5: the simplex method; Chapters 8–9: inverses in
  product form and the simplex method using multipliers, that is, the
  revised simplex method; Chapter 11: variants of the simplex algorithm,
  including the dual simplex method.)
- Lemke, C. E. (1954). The dual method of solving the linear programming
  problem. *Naval Research Logistics Quarterly, 1*(1), 36–47.
  <https://doi.org/10.1002/nav.3800010107>
- Luenberger, D. G., & Ye, Y. (2016). *Linear and nonlinear programming*
  (4th ed.). Springer. <https://doi.org/10.1007/978-3-319-18842-3>
  (Chapters 3–4: the simplex method, the revised simplex method, duality and
  the dual simplex method.)
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 3: the simplex method and two-phase method; Chapter 4: dual
  simplex and post-optimal analysis; Chapter 7: the revised simplex method.)
