# Linear Programming Algorithms in Python

## Overview

{doc}`09-02-code-lp-modeling` let a solver do the work. This chapter opens
the box and implements the algorithms of {doc}`03-02-lp-geometry` to
{doc}`03-07-sensitivity-analysis` in a few hundred lines of plain Python:
the enumeration of basic solutions, the tableau simplex method with
Dantzig's and Bland's rules, the two-phase method, the revised simplex
method with the product form of the inverse, the dual simplex method, and
the cost and right-hand-side ranging of sensitivity analysis. The code
follows the book's conventions exactly: tableaux have the objective row
last, labeled $-z$, holding the reduced costs $\bar c_j$ and, in the
right-hand-side column, $-z$; minimization is optimal when all
$\bar c_j\ge0$, maximization when all $\bar c_j\le0$. All arithmetic is
**exact**, with `fractions.Fraction`, so every tableau of Part III is
reproduced entry by entry: the product-mix run of
{prf:ref}`ex-03-03-simplex-run`, Beale's cycling example
{prf:ref}`ex-03-03-cycling`, the two-phase example {prf:ref}`ex-03-04-two-phase`,
the revised run {prf:ref}`ex-03-06-revised-run`, the dual simplex examples of
{doc}`03-06-revised-dual-simplex` and {doc}`03-07-sensitivity-analysis`, and the
ranges of {numref}`tab-03-07-ranges`. Every result is checked against
`linprog(method="highs")`.

These implementations are for learning and for checking hand computations.
They are dense, exact and slow; production codes are revised, sparse,
floating-point and far more careful about numerics
({doc}`09-01-code-setup-solvers`). As in the rest of Part IX, the code is not
run when the book is built: copy the blocks in order into a script or
notebook and compare the output with the `# Expected:` comments.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- compute with exact fractions and print tableaux in the book's layout;
- enumerate the basic solutions of $\mA\vx=\vb$ and find an optimal one;
- implement the pivot, the entering and leaving rules and the tableau
  simplex method for minimization and maximization, and detect cycling;
- implement Phase I and Phase II, including the removal of artificial
  variables that stay basic at zero;
- implement the revised simplex method with eta matrices and the dual
  simplex method; and
- compute cost and right-hand-side ranges from an optimal basis and verify
  them with a solver.
:::

All blocks use the following imports.

```python
import itertools
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

np.set_printoptions(precision=3, suppress=True, linewidth=120)
```

(sec-09-03-exact)=
## Exact arithmetic and printing

A `Fraction` stores a rational number as a pair of integers and keeps
`+ - * /` exact: `Fraction(1, 3) + Fraction(1, 6)` is `1/2`, not
`0.49999999999999994`. Floating-point simplex codes need tolerances in
every sign test ({ref}`sec-09-01-tolerances`); with fractions the tests
$\bar c_j<0$ and $\bar a_{is}>0$ are exact, and degenerate ties are real
ties. A matrix is a list of rows, each a list of fractions. Three small
helpers convert data, solve linear systems exactly by Gauss–Jordan
elimination ({doc}`01-01-linear-algebra`), and print matrices and tableaux.

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

A **tableau** is stored as a list of $m+k$ rows of length $n+1$: the $m$
constraint rows $[\bar a_{i1},\dots,\bar a_{in}\mid\bar b_i]$, followed by
$k\ge1$ objective rows, the $-z$ row last (in Phase I the $-\phi$ row comes
after it, as in {doc}`03-04-two-phase-big-m`). The basis is a list `basis` of
column indices, `basis[i]` being the basic variable of row $i$. Indices are
0-based in the code, so $x_1$ is column 0.

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
## Enumerating basic solutions

By the fundamental theorem of linear programming
({prf:ref}`thm-03-02-fundamental`) a linear program (LP) in standard form
with an optimal solution has an optimal basic feasible solution (BFS), so
for a small problem it suffices to try all $\binom nm$ choices of $m$
columns: skip the singular ones, solve $\mB\vx_B=\vb$, keep the
solutions with $\vx_B\ge\vzero$ ({prf:ref}`def-03-02-basic-solution`).

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

For the product mix in standard form {eq}`eq-03-02-product-mix-standard`,
with slacks $x_3,x_4,x_5$, the ten candidate bases of
{numref}`tab-03-02-product-mix-basic` follow, and the best feasible one is
optimal:

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

Five bases are feasible, one is singular ($x_2,x_3,x_4$: the third row of
$\mB$ is zero), and the best feasible basic solution is the optimum
$(\tfrac{45}{7},\tfrac{30}{7})$ with $z^*=\tfrac{36000}{7}$. Some
infeasible basic solutions have larger $z$; they are not candidates.
Enumeration grows like $\binom nm$, already $\binom{20}{10}=184\,756$ for a
tiny LP, which is why the simplex method moves from one basic feasible
solution to a better neighbor instead.

(sec-09-03-tableau)=
## The tableau simplex method

### Pivot, entering and leaving rules

The **pivot** on $\bar a_{rs}$ is the Gauss–Jordan step {eq}`eq-03-03-pivot`:
divide row $r$ by the pivot element and subtract multiples of it from every
other row, including the objective rows. A tableau in canonical form for
any basis is obtained from the plain data $[\mA\mid\vb]$, $[\vc^\T\mid0]$
by pivoting on the basic columns one row at a time; for a slack basis
these pivots change nothing.

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

The rules of {numref}`tab-03-03-rules` choose the pivot. The entering column
is chosen from the **last** row of the tableau, which is the $-z$ row in
Phase II and the $-\phi$ row in Phase I. Dantzig's rule takes the most
negative reduced cost for a minimization (most positive for a
maximization), ties going to the smallest index; Bland's rule
({prf:ref}`alg-03-03-bland`) takes the smallest index with an improving
reduced cost. The ratio test {eq}`eq-03-03-ratio-test` breaks ties by the
smallest index of the leaving basic variable, as in
{prf:ref}`alg-03-03-simplex`.

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

### The algorithm

{prf:ref}`alg-03-03-simplex` is now a short loop. It records every pivot as
(entering, leaving, $\theta$) and stops with status `"optimal"`,
`"unbounded"`, or `"cycling"` when a basis repeats, which can only happen
at a degenerate vertex ({prf:ref}`thm-03-03-finite`).

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

### The product-mix run

Starting from the slack basis $\{x_3,x_4,x_5\}$, Dantzig's rule for the
maximization reproduces the three pivots of {prf:ref}`ex-03-03-simplex-run`:
$x_1$ replaces $x_5$ ($\theta=8$), $x_2$ replaces $x_3$
($\theta=\tfrac{12}{5}$) and $x_5$ replaces $x_4$
($\theta=\tfrac{11}{7}$), along the path of {numref}`fig-03-03-simplex-path`.

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

The last tableau is {numref}`tab-03-03-pm-3`, including the entry
$-\tfrac17$ in row $x_2$, column $x_3$, and the shadow prices
$\tfrac{550}{7}$ and $\tfrac{20}{7}$ appear (negated) under the slack
columns of the $-z$ row. The solver agrees:

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

### Cycling and Bland's rule

Beale's example {prf:ref}`ex-03-03-cycling` is already in canonical form for
the basis $\{x_1,x_2,x_3\}$. Under Dantzig's rule every pivot is degenerate
($\theta=0$), and after six pivots the starting basis returns.

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

The Dantzig sequence is that of {numref}`tab-03-03-cycling`; Bland's rule
makes the same first three pivots, then takes $x_1$ instead of $x_7$ and
reaches the optimum $\vx^*=(\tfrac34,0,0,1,0,1,0)$, $z^*=-\tfrac54$, in
three more pivots, as in {ref}`sec-03-03-degeneracy`. HiGHS confirms the
optimum (its own anti-cycling devices are perturbations, not Bland's rule):

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
## The two-phase method

When no slack basis is available, {prf:ref}`alg-03-04-two-phase` adds an
artificial variable $R_i$ to the rows that need one and minimizes their sum
$\phi$ in Phase I. The function below adds an artificial variable to every
row (after making $\vb\ge\vzero$), which is what the example needs and
always works; a production code would reuse slack columns where they exist.
The $-z$ row is carried along during Phase I, and the $-\phi$ row
{eq}`eq-03-04-phase-one-row` below it drives the pivots. After Phase I, an
artificial variable still basic at level zero is pivoted out on any
nonzero entry of an original column; if its row has none, the row is
redundant and is deleted ({prf:ref}`thm-03-04-phase-one`).

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

{prf:ref}`ex-03-04-two-phase` maximizes $-3x_1+x_3$ over three equations.
Phase I needs two pivots and ends with $R_3$ basic at zero; the clean-up
pivots it out on the entry $-1$ under $x_4$; Phase II needs one pivot.

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

Both printed tableaux are those of {prf:ref}`ex-03-04-two-phase`, and the
optimum $\vx^*=(0,\tfrac52,\tfrac32,0)$, $z^*=\tfrac32$ agrees with HiGHS.
An infeasible LP ends Phase I with $\phi^*>0$; the running example of
{doc}`03-07-sensitivity-analysis` with $b_1=15$ is one ({prf:ref}`ex-03-07-rhs`):

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

With $b_1=11$ the same function finds the optimum $\vx^*=(3,4,2,0,0,0)$,
$z^*=11$ of {prf:ref}`ex-03-07-running-lp` from scratch.

(sec-09-03-revised)=
## The revised simplex method

{prf:ref}`alg-03-06-revised-simplex` never forms the tableau. It keeps the
original data, the current basic values $\vx_B$ and the inverse of the
basis in **product form** {eq}`eq-03-06-product-form`, a list of eta matrices
$\mE_k$, each stored as its pivot position $r$ and eta column
$\veta$ {eq}`eq-03-06-eta`. Applying $\mE$ to a column changes
every entry by a multiple of entry $r$ ({eq}`eq-03-06-eta-action`); applying
it to a row vector from the right changes only entry $r$, which becomes
$\vu^\T\veta$.

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

The algorithm starts from a slack basis, whose matrix is the identity, so
the list of etas starts empty. Pricing applies the etas to $\vc_B^\T$ from
the last to the first; the entering column is updated from the first to
the last. The entering rule is Dantzig's, but a list `first_entering` of
column indices can prescribe the first choices, to replay a hand
computation.

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

{prf:ref}`ex-03-06-revised-run` lets $x_2$ enter first and then follows
Dantzig's rule. The three eta matrices and the final inverse
$\mB_4^{-1}=\mE_3\mE_2\mE_1$ are those of the example:

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

The optimum $\vx^*=(\tfrac15,0,\tfrac85,0,0,4)$, $z^*=-\tfrac{27}{5}$ and
the multipliers $\vpi=(-\tfrac65,-\tfrac35,0)$ agree with the example and
with $z^*=\vpi^\T\vb$. With Dantzig's rule from the start, the same basis
is reached in two pivots, as remarked after the example:

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
## The dual simplex method

{prf:ref}`alg-03-06-dual-simplex` works on a tableau whose reduced costs are
all nonnegative (dual feasible, minimization) but whose right-hand side
has negative entries. It chooses the leaving row first (the most negative
$\bar b_r$), stops with infeasibility if that row has no negative entry,
and otherwise lets the column with the smallest dual ratio
$\bar c_j/(-\bar a_{rj})$, $\bar a_{rj}<0$, enter
({eq}`eq-03-06-dual-ratio`).

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

{prf:ref}`ex-03-06-dual-simplex` gives the tableau directly. One pivot, $x_2$
leaving and $x_4$ entering with dual ratio $\tfrac12$, makes it optimal
with $z^*=1$:

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

The dual simplex method is the tool for re-optimizing after a change of
$\vb$ or a new constraint ({ref}`sec-03-07-reoptimisation`). For
{prf:ref}`ex-03-07-rhs`, $b_1=15$: the old optimal basis $(x_1,x_2,x_3)$
stays dual feasible, two basic values become negative, and the method
proves infeasibility after one pivot. For {prf:ref}`ex-03-07-new-constraint`,
the row $x_1+x_2+x_7=6$ with a new slack $x_7$ is appended, brought into
canonical form by subtracting the rows of the basic variables, and one
pivot restores feasibility.

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

Both runs match {doc}`03-07-sensitivity-analysis`: with $b_1=15$ the row
$x_1+2x_3+3x_4=-5$ has no negative entry and proves the LP infeasible; with
the new constraint the optimum moves to $\vx^*=(1,5,3,0,1,0)$, $z^*=14$.
HiGHS agrees:

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
## Sensitivity ranging

Given an optimal basis, the ranges of {doc}`03-07-sensitivity-analysis` need
only $\mB^{-1}$, the optimal tableau and the reduced costs:

- a **nonbasic cost** $c_r$ may grow without limit and fall to the
  break-even cost $\vpi^\T\mA_{\cdot r}=c_r-\bar c_r$
  ({eq}`eq-03-07-nonbasic-range`);
- a **basic cost** $c_r$, in row $i$, may change by $\delta$ in the interval
  {eq}`eq-03-07-basic-range`, from the ratios $\bar c_j/\bar a_{ij}$ over the
  nonbasic columns;
- a **right-hand side** $b_r$ may change by $\delta$ while
  $\vx_B+\delta\,\mB^{-1}\ve_r\ge\vzero$ ({eq}`eq-03-07-rhs-range`).

The function below returns all three kinds of ranges for a minimization,
with `None` for an infinite end.

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

These are the ranges of {numref}`tab-03-07-ranges`. A solver check makes the
meaning concrete: just inside a cost range the optimal solution stays
$(3,4,2,0,0,0)$, just outside it moves; just inside a right-hand-side range
the optimal basis (the set of positive variables) stays
$\{x_1,x_2,x_3\}$, just outside it changes or the LP becomes infeasible.

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

Inside the right-hand-side ranges the optimal value changes at the rate
$\pi_i$; together with the $\mT$-matrix of {eq}`eq-03-07-T-running`, the
bordered inverse that turns any column of data into a column of the
optimal tableau ({eq}`eq-03-06-T-column`), this is all a sensitivity report
contains. The $\mT$-matrix is one exact inverse away:

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

The last row of $\mT$ is $(-\vpi^\T,1)=(-2,4,-1,1)$, and the new activity
of {prf:ref}`ex-03-07-new-activity` (b) has reduced cost $-7$ with updated
column $(1,1,-2)^\T$, as computed there.

## Summary and Common Mistakes

### Summary

- Exact `Fraction` arithmetic reproduces the book's tableaux digit for
  digit ({ref}`sec-09-03-exact`).
- Enumerating the $\binom nm$ bases finds the optimal BFS of small LPs;
  the product mix has five feasible bases out of ten
  ({ref}`sec-09-03-basic-solutions`).
- The tableau simplex method is a pivot {eq}`eq-03-03-pivot` plus entering
  and leaving rules; Dantzig's rule solves the product mix in three pivots
  and cycles on Beale's example, Bland's rule reaches $z^*=-\tfrac54$
  ({ref}`sec-09-03-tableau`).
- The two-phase method adds artificial variables, minimizes their sum,
  pivots out degenerate artificials and continues with the carried $-z$
  row ({ref}`sec-09-03-two-phase`).
- The revised simplex method keeps only $\vx_B$ and eta factors of
  $\mB^{-1}$ ({ref}`sec-09-03-revised`); the dual simplex method keeps dual
  feasibility and re-optimizes after changes of $\vb$ or new constraints
  ({ref}`sec-09-03-dual-simplex`).
- Cost and right-hand-side ranges and the $\mT$-matrix follow from the
  optimal basis and agree with re-solving by HiGHS
  ({ref}`sec-09-03-sensitivity`).

### Common mistakes

- **Pricing on the wrong row.** In Phase I the entering variable comes from
  the $-\phi$ row, not from the $-z$ row; both rows are updated by every
  pivot.
- **Ignoring the sense.** A maximization tableau is optimal when all
  $\bar c_j\le0$; using the minimization test runs the method backwards.
- **Floating-point ties.** With floats, a degenerate tie such as $3/3$
  against $6/6$ may be broken by rounding noise; exact fractions or a
  tolerance in every comparison are needed.
- **Forgetting the artificial variables that stay basic** at zero after
  Phase I, or deleting a row that is not redundant.
- **Applying eta matrices in the wrong order.** Columns are updated with
  $\mE_1$ first, row vectors ($\vpi^\T$) with the last eta first.
- **Starting the dual simplex method from a tableau that is not dual
  feasible**; the function checks this and refuses.
- **Reading a cost range for a basic variable from its column instead of
  its row** of the optimal tableau.

## References

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
