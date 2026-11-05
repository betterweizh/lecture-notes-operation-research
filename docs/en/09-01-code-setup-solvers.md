# Software Setup and Solvers

## Overview

Part IX turns the methods of Parts II–VII into Python code. This chapter
prepares the ground. It sets up a Python environment with the few packages
that the code chapters use, surveys the **solvers** that do the numerical
work (open-source and commercial solvers for linear programs (LPs) and
mixed-integer linear programs (MILPs), conic solvers for convex problems) and the **modeling languages** in which models are written for
them, and explains how to read what a solver returns: a status, an optimal
value, primal values, dual values with their sign conventions, slacks, and
the tolerances that make all of these approximate. It ends with a short
installation check that solves the product-mix LP of
{prf:ref}`ex-03-01-product-mix` with SciPy, PuLP and CVXPY and compares the
three answers with the exact optimum $\vx^*=(\tfrac{45}{7},\tfrac{30}{7})$,
$z^*=\tfrac{36000}{7}\approx\$5142.86$ per week. The chapter assumes the
product-mix model ({doc}`03-01-lp-model-building`), standard form and slacks
({doc}`03-02-lp-geometry`), shadow prices ({doc}`03-05-lp-duality`), and
basic Python with NumPy arrays.

The code in Part IX is **not run when the book is built**. Copy the blocks
of a chapter, in order, into a script or a notebook (each chapter is one
namespace, like the cells of a notebook) and compare the printed results
with the `# Expected:` comments, which repeat the numbers of the worked
examples. Small differences in the last printed digit, or `-0.` instead of
`0.`, come from floating-point arithmetic and are harmless.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- create a virtual environment, install the packages of Part IX and check
  that they work;
- name the main open-source and commercial LP, MILP and conic solvers and
  say which problem classes each one handles;
- compare modeling languages and frameworks (AMPL, GAMS, Pyomo, JuMP,
  OR-Tools, PuLP, CVXPY, YALMIP) and read a small AMPL model;
- read a solver's status, optimal value, primal values, duals and slacks,
  and convert dual values between the sign conventions of SciPy, PuLP and
  CVXPY and the shadow prices of {doc}`03-05-lp-duality`; and
- explain feasibility, optimality and integrality tolerances and compare
  floating-point results with exact values accordingly.
:::

(sec-09-01-environment)=
## Setting up Python

The code needs Python 3.11 or later (the oldest version that recent SciPy
releases support). It is good practice to install the
packages of a project into a **virtual environment**, a folder that holds
its own copy of the interpreter's package directory, so that the versions
used here do not interfere with other projects. From the root of the book's
repository:

```bash
python -m venv .venv                 # create the environment once
source .venv/bin/activate            # Linux and macOS
.venv\Scripts\activate               # Windows (PowerShell or cmd)
pip install -r requirements.txt      # numpy, scipy, networkx, pulp, cvxpy, ...
```

The file `requirements.txt` lists the packages with minimum versions; the
same file also installs the tools that build the book. Without the
repository, the code chapters need only

```bash
pip install numpy scipy networkx matplotlib pulp cvxpy
```

Jupyter (`pip install jupyterlab`) is convenient for running the blocks
interactively but not required. Every later block assumes the imports of the
first block of its chapter. For this chapter they are:

```python
import importlib.metadata
import math
import sys

import cvxpy as cp
import numpy as np
import pulp
from scipy.optimize import linprog

np.set_printoptions(precision=3, suppress=True, linewidth=120)
```

The second block checks the interpreter and the installed versions against
the minimum versions of `requirements.txt`. Exact version numbers differ
from machine to machine; what matters is that every check prints `True`.

```python
MINIMUM_VERSIONS = {
    "numpy": (1, 25),
    "scipy": (1, 16),
    "networkx": (3, 0),
    "matplotlib": (3, 7),
    "pulp": (2, 7),
    "cvxpy": (1, 4),
}


def major_minor(version: str) -> tuple[int, int]:
    """Return the (major, minor) numbers of a version string like '1.26.4'.

    Parameters
    ----------
    version : str
        Version string as reported by ``importlib.metadata.version``.

    Returns
    -------
    tuple of int
        The first two numeric components.
    """
    parts = version.split(".")
    if len(parts) < 2:
        raise ValueError(f"unexpected version string {version!r}")
    # Leading digits only: "0rc1" gives 0, not 1.
    minor = ""
    for ch in parts[1]:
        if not ch.isdigit():
            break
        minor += ch
    return int(parts[0]), int(minor or "0")


print("Python >= 3.11:", sys.version_info >= (3, 11))
for name, minimum in MINIMUM_VERSIONS.items():
    try:
        version = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        print(f"{name:<11} not installed")
        continue
    ok = major_minor(version) >= minimum
    print(f"{name:<11} {version:<10} ok: {ok}")
# Expected (version numbers differ):
# Python >= 3.11: True
# numpy       2.x.x      ok: True
# scipy       1.x.x      ok: True
# ... one line per package, each ending in "ok: True"
```

(sec-09-01-packages)=
## The packages of Part IX

Part IX uses six packages besides the standard library. {numref}`tab-09-01-packages`
lists them with their role. The standard library's `fractions.Fraction`
deserves a mention of its own: the hand-sized tableaux of
{doc}`09-03-code-lp-algorithms` are computed with exact rational numbers, so
that they reproduce the printed tableaux of Part III digit for digit.

```{list-table} Packages used in Part IX
:name: tab-09-01-packages
:header-rows: 1

* - Package
  - What it is used for
  - Chapters
* - NumPy
  - arrays, matrix–vector products, linear systems (`numpy.linalg`)
  - all
* - SciPy
  - `scipy.optimize.linprog` and `milp` (both call HiGHS) for LPs and MILPs;
    `minimize` for nonlinear programs (NLPs)
  - {doc}`09-02-code-lp-modeling`, {doc}`09-03-code-lp-algorithms`,
    {doc}`09-04-code-networks`, {doc}`09-06-code-nlp`
* - PuLP
  - algebraic modeling of LPs and MILPs with named variables and
    constraints; ships with the CBC solver
  - {doc}`09-02-code-lp-modeling`
* - CVXPY
  - modeling of convex problems (LPs, quadratic programs (QPs) and
    second-order cone programs (SOCPs)) in vector form, with
    automatic convexity checking; calls conic solvers
  - {doc}`09-02-code-lp-modeling`, {doc}`09-06-code-nlp`
* - NetworkX
  - graphs and reference implementations of network algorithms
  - {doc}`09-04-code-networks`
* - Matplotlib
  - plots, only where a picture is essential
  - {doc}`09-06-code-nlp`
```

Two remarks help to choose among them. **SciPy** expects a model in
*matrix form*, $\min\{\vc^\T\vx:\mA_{ub}\vx\le\vb_{ub},\
\mA_{eq}\vx=\vb_{eq},\ \vl\le\vx\le\vu\}$: compact and fast, but the
modeler must number the variables and build the matrices. **PuLP** and
**CVXPY** let the modeler write constraints as expressions in named
variables and build the matrices themselves. PuLP is limited to linear and
mixed-integer linear models; CVXPY handles any convex model that follows its
rules (disciplined convex programming), but integer variables only with a
solver that supports them.

(sec-09-01-solvers)=
## The solver landscape

A **solver** is the program that runs an algorithm on a model given as
data: matrices, bounds, variable types. A **modeling layer** (PuLP, CVXPY,
AMPL, …) translates the model as the modeler writes it into that data and
the solver's answer back into named quantities. Most modeling layers can
call several solvers, so the two choices are largely independent.

### Algorithms inside LP and MILP solvers

LP solvers implement three families of methods. The **primal simplex
method** of {doc}`03-03-simplex-method` and the **dual simplex method** of
{doc}`03-06-revised-dual-simplex` move between bases, in the revised form with
an LU factorization of $\mB$ instead of a tableau; the dual simplex method
is usually the default because it re-optimizes quickly after a change.
**Interior-point** (barrier) methods follow a path through the interior of
the feasible region and need far fewer, but more expensive, iterations; a
final **crossover** step recovers an optimal basis when one is wanted (for
example for sensitivity analysis). MILP solvers run **branch and bound** on
the LP relaxation, strengthened by **cutting planes**, presolve reductions
and primal heuristics. All of them work in floating-point arithmetic with
tolerances ({ref}`sec-09-01-tolerances`).

### Open-source LP and MILP solvers

- **HiGHS** (Huangfu & Hall, 2018) is a high-performance solver for LP,
  MILP and convex QP developed at the University of Edinburgh and
  distributed under the MIT license. It has dual and primal revised simplex
  codes, an interior-point code and a branch-and-cut MILP solver. SciPy's
  `linprog` uses HiGHS by default (`method="highs"`), and
  `scipy.optimize.milp` calls the HiGHS MILP solver; HiGHS also has its own
  Python interface `highspy`.
- **CBC** (COIN-OR Branch and Cut) is the MILP solver of the COIN-OR
  project; it solves LPs with the COIN-OR simplex code CLP. A CBC binary is
  bundled with PuLP, which is why PuLP works without any further
  installation.
- **GLPK** (GNU Linear Programming Kit) is an older, slower but very
  portable LP/MILP package with its own modeling language (GNU MathProg, a
  subset of AMPL); PuLP calls it through `GLPK_CMD` when `glpsol` is on the
  path.
- **SCIP** is a constraint-integer-programming framework that solves MILPs
  and mixed-integer nonlinear programs (MINLPs); it is distributed under the Apache
  2.0 license (since version 8.0.3) and has the Python interface PySCIPOpt.

### Conic solvers behind CVXPY

CVXPY reformulates a convex model as a **conic program**: a linear
objective, linear equations and the requirement that some affine expressions
lie in cones (the nonnegative orthant for LP, second-order cones for SOCP,
the cone of positive semidefinite matrices for semidefinite programs, or
SDPs; see
{doc}`02-03-convex-optimization-duality`). It then calls one of these solvers:

- **Clarabel**, an interior-point solver written in Rust for LP, QP, SOCP,
  SDP and exponential and power cones; recent CVXPY versions install it by
  default and prefer it for most problems;
- **SCS** (O'Donoghue et al., 2016), a first-order splitting method (the
  alternating direction method of multipliers, ADMM)
  for large conic problems, accurate to about $10^{-4}$ to $10^{-6}$ rather
  than to machine precision;
- **OSQP** (Stellato et al., 2020), an ADMM solver for convex QPs, fast and
  warm-startable, popular in control;
- **ECOS** (Domahidi et al., 2013), an interior-point SOCP solver for
  embedded systems, the former default of CVXPY, now an optional install.

CVXPY can also call HiGHS, GLPK, CBC and SCIP for LPs and MILPs, and the
commercial solvers below. `cp.installed_solvers()` lists what is available.

### Commercial solvers

The commercial solvers are faster on large and hard MILPs, have better
presolve and more robust numerics, and report sensitivity ranges directly.
All four below give **free academic licenses** to students and staff of
degree-granting institutions; registration with a university e-mail address
gives a license key that is activated on the user's machine (for Gurobi
with the command `grbgetkey <key>` from a university network, for CPLEX
through IBM's academic program). The Python packages of Gurobi, CPLEX and
Xpress install with `pip` together with a size-limited license for small
models; MOSEK's package also installs with `pip` but needs a trial or
academic license.

- **Gurobi** (LP, QP, SOCP, MILP, mixed-integer QP, and nonconvex quadratic
  problems);
  Python API `gurobipy`.
- **IBM ILOG CPLEX** (LP, QP, SOCP, MILP), part of CPLEX Optimization
  Studio with the OPL modeling language; Python API `cplex` and the
  modeling layer `docplex`.
- **FICO Xpress** (LP, QP, MILP, nonlinear); Python API `xpress`.
- **MOSEK** (LP, QP, conic problems including SDP, mixed-integer conic);
  particularly strong on conic problems; Python API `mosek`.

Switching solver does not change the model. In PuLP the solver is an
argument of `solve`, for example `prob.solve(pulp.GUROBI_CMD())` or
`pulp.CPLEX_CMD()`, `pulp.XPRESS()`, `pulp.MOSEK()`, `pulp.GLPK_CMD()`;
`pulp.listSolvers(onlyAvailable=True)` shows which ones PuLP can find. In
CVXPY it is a keyword, `problem.solve(solver=cp.GUROBI)` (or `cp.CPLEX`,
`cp.XPRESS`, `cp.MOSEK`, `cp.HIGHS`, `cp.SCS`, …). {numref}`tab-09-01-solvers`
summarizes the landscape.

```{list-table} Solvers and the problem classes they handle (✓ = supported)
:name: tab-09-01-solvers
:header-rows: 1

* - Solver
  - License
  - LP
  - MILP
  - QP / SOCP
  - SDP
  - Reached from
* - HiGHS
  - open (MIT)
  - ✓
  - ✓
  - convex QP
  -
  - SciPy, PuLP, CVXPY, Pyomo, JuMP
* - CBC / CLP
  - open (EPL)
  - ✓
  - ✓
  -
  -
  - PuLP (bundled), CVXPY, Pyomo, JuMP
* - GLPK
  - open (GPL)
  - ✓
  - ✓
  -
  -
  - PuLP, CVXPY, Pyomo
* - SCIP
  - open (Apache)
  - ✓
  - ✓
  - ✓ (and MINLP)
  -
  - PySCIPOpt, CVXPY, Pyomo
* - Clarabel, ECOS
  - open
  - ✓
  -
  - ✓
  - Clarabel only
  - CVXPY
* - SCS
  - open (MIT)
  - ✓
  -
  - ✓
  - ✓
  - CVXPY
* - OSQP
  - open
  -
  -
  - QP only
  -
  - CVXPY
* - Gurobi
  - commercial, free academic
  - ✓
  - ✓
  - ✓
  -
  - gurobipy, PuLP, CVXPY, Pyomo, AMPL, GAMS
* - CPLEX
  - commercial, free academic
  - ✓
  - ✓
  - ✓
  -
  - cplex/docplex, PuLP, CVXPY, Pyomo, AMPL, GAMS
* - Xpress
  - commercial, free academic
  - ✓
  - ✓
  - ✓
  -
  - xpress, PuLP, CVXPY, Pyomo, AMPL, GAMS
* - MOSEK
  - commercial, free academic
  - ✓
  - ✓
  - ✓
  - ✓
  - mosek, PuLP, CVXPY, Pyomo, AMPL, GAMS
```

For the models of this book every solver in the table is more than fast
enough; the differences matter for MILPs with thousands of binaries and for
LPs with millions of variables.

(sec-09-01-modelling-languages)=
## Modeling languages and frameworks

An **algebraic modeling language** lets the modeler write
$\sum_{j\in J}a_{ij}x_j\le b_i$ for all $i\in I$ almost as on paper, with
index sets, parameters, variables, objective and constraints as separate
declarations, and with the data kept apart from the model. The best known
is **AMPL** (Fourer et al., 2003). An AMPL project usually has
three files: the model (`.mod`), the data (`.dat`) and a script of
commands (`.run`). The product-mix model of {prf:ref}`ex-03-01-product-mix`
reads, as an illustration,

```text
# product_mix.mod -- the model, without numbers
set PRODUCTS;
set RESOURCES;
param profit {PRODUCTS};                  # $ per hundred cases
param use {RESOURCES, PRODUCTS};          # resource use per hundred cases
param avail {RESOURCES};                  # amount available
param demand {PRODUCTS} default Infinity; # demand limit, hundred cases
var Make {p in PRODUCTS} >= 0, <= demand[p];
maximize Contribution: sum {p in PRODUCTS} profit[p] * Make[p];
subject to Capacity {r in RESOURCES}:
    sum {p in PRODUCTS} use[r, p] * Make[p] <= avail[r];
```

```text
# product_mix.dat -- the numbers
set PRODUCTS := juice cocktail;
set RESOURCES := machine storage;
param profit := juice 500  cocktail 450;
param avail  := machine 60  storage 150;
param demand := juice 8;
param use:    juice  cocktail :=
    machine     6       5
    storage    10      20 ;
```

```text
# product_mix.run -- the commands
reset;
model product_mix.mod;
data product_mix.dat;
option solver highs;     # or cbc, gurobi, cplex, ...
solve;
display Make, Capacity.dual, Capacity.slack;
```

The declarations follow a few conventions worth copying in any language:
index sets in capitals, parameters in lower case, variables and constraint
names capitalized; a parameter's admissible range (`param avail {RESOURCES}
>= 0;`) is a check on the data, not a constraint of the model; every
constraint and the objective carry a name, through which their dual values
(`Capacity.dual`) and slacks are displayed. Changing the data file solves a
different instance of the same model.

The other main options are:

- **GAMS**, an algebraic modeling system as old as AMPL, widespread in
  economics and energy modeling, with its own language and a large set of
  solver links.
- **Pyomo** (Bynum et al., 2021), an open-source algebraic modeling
  library for Python with sets, parameters and indexed variables in the AMPL
  style; it supports LP, MILP and nonlinear models and almost every solver.
- **JuMP** (Lubin et al., 2023), the algebraic modeling package of the
  Julia language, fast to build large models and solver-independent.
- **Google OR-Tools**, a C++ library with Python bindings; it includes the
  LP solver GLOP, the constraint-programming and SAT solver CP-SAT (very
  strong on scheduling and pure-integer problems), vehicle-routing code, and
  a wrapper that calls HiGHS, SCIP, Gurobi and others.
- **PuLP** (Mitchell et al., 2011), a light Python library for LP and MILP
  models, used in {doc}`09-02-code-lp-modeling`; and **CVXPY** (Diamond &
  Boyd, 2016) for convex models.
- **MATLAB** with its Optimization Toolbox (`linprog`, `intlinprog`,
  `quadprog`, `fmincon`), and **YALMIP** (Löfberg, 2004), a free MATLAB
  toolbox that plays the role of CVXPY or PuLP: symbolic decision variables
  (`sdpvar`, `binvar`), constraints written as inequalities, and
  `optimize(constraints, objective)` with a choice of solver.

```{list-table} Modeling languages and frameworks
:name: tab-09-01-modelling
:header-rows: 1

* - Framework
  - Language
  - License
  - Problem classes
  - Typical use
* - AMPL
  - own language (`.mod`/`.dat`/`.run`), Python API `amplpy`
  - commercial (free community and academic editions)
  - LP, MILP, NLP, MINLP
  - teaching, industrial models, nonlinear models
* - GAMS
  - own language
  - commercial
  - LP, MILP, NLP, MINLP, equilibrium models
  - economics, energy systems
* - Pyomo
  - Python
  - open
  - LP, MILP, NLP, MINLP
  - research models in Python
* - JuMP
  - Julia
  - open
  - LP, MILP, conic, NLP
  - large models, research on algorithms
* - OR-Tools
  - C++, Python, Java, C#
  - open
  - LP, MILP, constraint programming, routing
  - scheduling and routing applications
* - PuLP
  - Python
  - open
  - LP, MILP
  - small and medium models, teaching
* - CVXPY
  - Python
  - open
  - convex (LP, QP, SOCP, SDP), mixed-integer convex with suitable solvers
  - convex optimization, statistics, control
* - YALMIP
  - MATLAB
  - free
  - LP, MILP, conic, NLP
  - control engineering in MATLAB
```

The code chapters use SciPy, PuLP and CVXPY because they are free, install
with `pip` and are enough for every model of the book. The concepts carry
over: a model built with named sets and parameters in PuLP is translated
into AMPL or Pyomo line by line.

(sec-09-01-reading-output)=
## Reading solver output

### Status

A solver's first answer is a **status**: did it find an optimal solution,
prove the problem infeasible or unbounded, or stop for another reason (an
iteration or time limit, numerical trouble)? Nothing else in the result is
meaningful until the status has been checked. {numref}`tab-09-01-status`
compares the three Python interfaces.

```{list-table} Status reports of the three interfaces
:name: tab-09-01-status
:header-rows: 1

* - Outcome
  - `linprog` / `milp` (`res.status`)
  - PuLP (`pulp.LpStatus[prob.status]`)
  - CVXPY (`problem.status`)
* - optimal
  - `0`
  - `"Optimal"`
  - `"optimal"`
* - iteration or time limit reached
  - `1`
  - `"Not Solved"`
  - `"user_limit"` or an error
* - infeasible
  - `2`
  - `"Infeasible"`
  - `"infeasible"`
* - unbounded
  - `3`
  - `"Unbounded"`
  - `"unbounded"`
* - numerical difficulties
  - `4`
  - `"Undefined"`
  - `"optimal_inaccurate"` and similar
```

Infeasible and unbounded results are the two outcomes of {doc}`03-02-lp-geometry`
that have no optimal solution, and in a model under construction they almost
always mean a modeling error: an infeasible model usually has a wrong sign
or unit, an unbounded one a missing constraint ({ref}`sec-03-01-modelling-steps`).
Some solvers cannot tell the two apart after presolve and report
"infeasible or unbounded"; solving again with the objective set to zero
decides which.

### Optimal value and primal values

With status optimal, the solver reports the optimal objective value and a
primal solution. Two details matter. First, `linprog` and `milp` always
*minimize*: a maximization $\max\vc^\T\vx$ is passed as $\min(-\vc)^\T\vx$,
and the optimal value is $-$`res.fun` ({eq}`eq-03-01-max-min`). Second, an LP
solver returns a basic optimal solution when it uses the simplex method or
crossover, but when the optimum is not unique ({ref}`sec-03-03-special-cases`)
which optimal vertex it returns depends on the solver, its version and its
settings. Two correct solvers may report different optimal solutions with
the same optimal value.

### Dual values and their sign conventions

For an optimal basis the solver also knows the simplex multipliers
$\vpi^\T=\vc_B^\T\mB^{-1}$ of {doc}`03-03-simplex-method`, which are optimal
dual values ({doc}`03-05-lp-duality`) and, at a nondegenerate optimum, the
shadow prices $\partial z^*/\partial b_i$ ({prf:ref}`thm-03-05-shadow-price`).
Solvers call them duals, dual values, marginals or shadow prices, and the
interfaces do not agree on their signs, because a sign depends on whether
the problem is a minimization or a maximization, on whether the row is
$\le$ or $\ge$, and on whether the interface reports the multiplier of a
Lagrangian or a derivative of the optimal value. The safe rule is to know
what each interface reports and to convert it to the book's shadow price,
the rate of change of the optimal value of the problem *as stated* per unit
increase of $b_i$:

```{list-table} Dual values and slacks in the three interfaces
:name: tab-09-01-duals
:header-rows: 1

* - Quantity
  - `linprog`
  - PuLP with CBC
  - CVXPY
* - dual value of a row
  - `res.ineqlin.marginals`, `res.eqlin.marginals`: derivative of the
    *minimized* objective `res.fun` with respect to the right-hand side
  - `prob.constraints[name].pi`: derivative of the objective of the
    problem as stated (max or min)
  - `constraint.dual_value`: Lagrange multiplier, nonnegative for
    inequality constraints in a minimization and in a maximization
* - shadow price of a $\le$ row in $\max\vc^\T\vx$
  - `-res.ineqlin.marginals` ($\ge0$)
  - `pi` ($\ge0$)
  - `dual_value` ($\ge0$)
* - shadow price of a $\le$ row in $\min\vc^\T\vx$
  - `res.ineqlin.marginals` ($\le0$)
  - `pi` ($\le0$)
  - `-dual_value` ($\le0$)
* - reduced costs $\bar c_j$
  - `res.lower.marginals` (for variables at their lower bound)
  - `variable.dj`
  - dual value of the constraint `x >= 0`
* - slack of a row
  - `res.slack` $=\vb_{ub}-\mA_{ub}\vx$; `res.con` for equations
  - `constraint.slack` $=$ right-hand side $-$ left-hand side
  - `constraint.expr.value` or the residual computed by hand
```

`linprog` has only $\le$ rows. A $\ge$ row $\va^\T\vx\ge\beta$ is passed as
$-\va^\T\vx\le-\beta$, which flips the sign of its marginal once more: the
derivative with respect to $\beta$ is $-1$ times the reported marginal.
Whatever the interface, a dual value can always be checked by its meaning:
solve again with $b_i$ increased by a small amount $\Delta$ and compare the
change of the optimal value with $\pi_i\Delta$, as
{ref}`sec-09-02-duals` does. Exact ranges over which the shadow prices hold
({doc}`03-07-sensitivity-analysis`) are reported by the commercial solvers
(and by AMPL with them); with the open-source interfaces used here they are
computed from the optimal basis, as in {ref}`sec-09-03-sensitivity`.

### Slacks

The **slack** of a $\le$ row is $b_i-\va_i^\T\vx^*\ge0$, the unused amount
of a resource; for a $\ge$ row the **surplus** $\va_i^\T\vx^*-b_i\ge0$ is
the amount by which a requirement is exceeded ({doc}`03-02-lp-geometry`). A row
with zero slack is **binding**. By complementary slackness
({doc}`03-05-lp-duality`) a row with positive slack has dual value zero, which
is a quick consistency check of any solver output. PuLP reports the slack of
every row as right-hand side minus left-hand side, so the surplus of a
$\ge$ row appears with a minus sign.

(sec-09-01-tolerances)=
## Numerical tolerances

Solvers work in double-precision floating-point arithmetic, with about 16
significant decimal digits, and every test inside them uses a
**tolerance**:

- the **primal feasibility tolerance** (typically $10^{-7}$) is the largest
  violation of a constraint or bound that still counts as feasible;
- the **dual feasibility** or **optimality tolerance** (typically
  $10^{-7}$) is the largest wrong-signed reduced cost that still counts as
  optimal ({prf:ref}`thm-03-03-optimality`);
- the **integrality tolerance** (typically $10^{-6}$ to $10^{-5}$) is how far
  an integer variable may be from an integer, so a binary may be reported as
  $0.9999999$ or $10^{-9}$;
- the **relative MIP gap** (typically $10^{-4}$, that is $0.01\%$) is how close
  the best integer solution must be to the best bound before a MILP solver
  stops; an "optimal" MILP result is optimal only up to this gap.

Two consequences follow for the code of this book. Results are compared
with expected values by `math.isclose` or `numpy.allclose` with a tolerance,
never with `==`, and printed with a fixed number of decimals. And badly
scaled models, with coefficients spanning many orders of magnitude (a big-$M$
of $10^9$ next to coefficients of order 1, {doc}`07-01-logical-constraints`),
make tolerances act in unexpected ways: a constraint violated by
$10^{-7}\cdot10^9=100$ units may still count as satisfied. Choosing units
that keep coefficients between about $10^{-3}$ and $10^{4}$, and big-$M$
constants as small as the bounds allow, avoids most numerical trouble.
SciPy passes tolerances to HiGHS through the `options` dictionary of
`linprog` and `milp`, for example the key `primal_feasibility_tolerance`
(value `1e-9`) for `linprog` and the key `mip_rel_gap` (value `1e-6`) for
`milp`.

Exact rational arithmetic avoids tolerances altogether, at the price of
speed. It is the right tool for checking small examples, and
{doc}`09-03-code-lp-algorithms` uses it for every tableau.

(sec-09-01-check)=
## Checking the installation: the product mix three ways

The product-mix LP {eq}`eq-03-01-product-mix` is small enough to know by
heart: maximize $500x_1+450x_2$ subject to machine time $6x_1+5x_2\le60$ h,
storage $10x_1+20x_2\le150$ m³ and demand $x_1\le8$ hundred cases, with
$\vx\ge\vzero$. Its optimum is $x_1^*=\tfrac{45}{7}\approx6.429$,
$x_2^*=\tfrac{30}{7}\approx4.286$ hundred cases per week, the optimal value
is $z^*=\tfrac{36000}{7}\approx\$5142.86$ per week, the shadow prices are
$\pi_1^*=\tfrac{550}{7}\approx\$78.57$ per machine hour,
$\pi_2^*=\tfrac{20}{7}\approx\$2.86$ per m³ and $\pi_3^*=0$, and only the
demand row has slack, $\tfrac{11}{7}\approx1.571$ hundred cases
({prf:ref}`ex-03-03-simplex-run`, {prf:ref}`ex-03-05-product-mix-dual`). The
data are stored once:

```python
# Product-mix data: contributions ($ per hundred cases), resource use per
# hundred cases (rows: machine h, storage m3, demand) and availabilities.
c = np.array([500.0, 450.0])
A = np.array(
    [
        [6.0, 5.0],
        [10.0, 20.0],
        [1.0, 0.0],
    ]
)
b = np.array([60.0, 150.0, 8.0])
z_exact = 36000 / 7
pi_exact = np.array([550 / 7, 20 / 7, 0.0])
```

**SciPy.** `linprog` minimizes, so the contributions enter with a minus
sign, and the shadow prices are the negated marginals
({numref}`tab-09-01-duals`).

```python
res = linprog(-c, A_ub=A, b_ub=b, bounds=(0, None), method="highs")
print("status:", res.status, res.success)
print("x* =", res.x)
print("z* =", round(-res.fun, 2))
print("slacks =", res.slack)
print("shadow prices =", -res.ineqlin.marginals)
# Expected:
# status: 0 True
# x* = [6.429 4.286]
# z* = 5142.86
# slacks = [0.    0.    1.571]
# shadow prices = [78.571  2.857  0.   ]
```

**PuLP.** The model is written with named variables and constraints and
solved by the bundled CBC; `msg=False` silences the solver log.

```python
prob = pulp.LpProblem("product_mix", pulp.LpMaximize)
x1 = pulp.LpVariable("juice", lowBound=0)
x2 = pulp.LpVariable("cocktail", lowBound=0)
prob += 500 * x1 + 450 * x2, "contribution"
prob += 6 * x1 + 5 * x2 <= 60, "machine"
prob += 10 * x1 + 20 * x2 <= 150, "storage"
prob += x1 <= 8, "demand"
prob.solve(pulp.PULP_CBC_CMD(msg=False))

print("status:", pulp.LpStatus[prob.status])
print("x* =", round(x1.value(), 3), round(x2.value(), 3))
print("z* =", round(pulp.value(prob.objective), 2))
for name, con in prob.constraints.items():
    print(f"{name:<8} pi = {con.pi:8.3f}   slack = {con.slack:6.3f}")
# Expected:
# status: Optimal
# x* = 6.429 4.286
# z* = 5142.86
# machine  pi =   78.571   slack = -0.000
# storage  pi =    2.857   slack = -0.000
# demand   pi =   -0.000   slack =  1.571
```

CBC writes its solution with about eight significant digits, which is why
`-0.000` can appear and why the results are rounded before printing. The
values agree with the exact ones within the tolerances of
{ref}`sec-09-01-tolerances`:

```python
pulp_x = np.array([x1.value(), x2.value()])
print(np.allclose(pulp_x, res.x, atol=1e-6))
print(math.isclose(pulp.value(prob.objective), z_exact, abs_tol=1e-4))
print(math.isclose(-res.fun, z_exact, rel_tol=1e-9))
# Expected:
# True
# True
# True
```

**CVXPY** (requires CVXPY). The model is written in vector form with the
same arrays; the constraint object is kept so that its dual values can be
read after the solve. CVXPY's dual values of inequality constraints are
nonnegative multipliers, which for a maximization are the shadow prices
directly.

```python
x = cp.Variable(2, nonneg=True)
resources = A @ x <= b
problem = cp.Problem(cp.Maximize(c @ x), [resources])
problem.solve()
print("status:", problem.status)
print("x* =", np.round(x.value, 3))
print("z* =", round(problem.value, 2))
print("shadow prices =", np.round(resources.dual_value, 3))
print("matches exact:", np.allclose(resources.dual_value, pi_exact, atol=1e-5))
print("solver used:", problem.solver_stats.solver_name)
# Expected:
# status: optimal
# x* = [6.429 4.286]
# z* = 5142.86
# shadow prices = [78.571  2.857  0.   ]
# matches exact: True
# solver used: CLARABEL   (or another installed solver, for example ECOS)
```

Interior-point solvers such as Clarabel stop close to, but not exactly at,
the optimal vertex, so their values agree with the simplex results to
about six digits and a non-binding dual may print as `0.` or as a tiny
number.

**Two failures on purpose.** The status must be read before anything else.
Requiring at least 9 hundred cases of juice glasses, which contradicts the
demand limit $x_1\le8$, makes the model infeasible; forgetting the machine
and storage rows makes it unbounded, because nothing then limits cocktail
glasses.

```python
STATUS = {
    0: "optimal",
    1: "iteration limit",
    2: "infeasible",
    3: "unbounded",
    4: "numerical difficulties",
}

# Infeasible: x1 >= 9 is passed to linprog as -x1 <= -9.
A_inf = np.vstack([A, [-1.0, 0.0]])
b_inf = np.append(b, -9.0)
res_inf = linprog(-c, A_ub=A_inf, b_ub=b_inf, method="highs")
print("x1 >= 9 added:", STATUS[res_inf.status])

# Unbounded: only the demand row x1 <= 8 is kept.
res_unb = linprog(-c, A_ub=A[2:], b_ub=b[2:], method="highs")
print("machine and storage rows dropped:", STATUS[res_unb.status])
# Expected:
# x1 >= 9 added: infeasible
# machine and storage rows dropped: unbounded
```

In both cases `res.x` and `res.fun` must not be used. Finally,
`pulp.listSolvers(onlyAvailable=True)` and (with CVXPY)
`cp.installed_solvers()` show which solvers each interface can call.

```python
print(pulp.listSolvers(onlyAvailable=True))
# Expected: ['PULP_CBC_CMD']   (plus any other solver installed)
```

```python
print(cp.installed_solvers())
# Expected (requires CVXPY; the list depends on the installation):
# ['CLARABEL', 'ECOS', 'HIGHS', 'OSQP', 'SCIPY', 'SCS', ...]
```

If all blocks of this section print their expected values, the environment
is ready for the rest of Part IX.

## Summary and Common Mistakes

### Summary

- A virtual environment with `pip install -r requirements.txt` provides
  NumPy, SciPy, NetworkX, Matplotlib, PuLP (with CBC) and CVXPY, the
  packages of Part IX ({numref}`tab-09-01-packages`).
- Open-source LP/MILP solvers are HiGHS (behind SciPy's `linprog` and
  `milp`), CBC (bundled with PuLP), GLPK and SCIP; CVXPY calls conic solvers
  such as Clarabel, SCS, OSQP and ECOS; Gurobi, CPLEX, Xpress and MOSEK are
  commercial with free academic licenses and plug into PuLP and CVXPY as a
  solver argument ({numref}`tab-09-01-solvers`).
- Modeling languages separate model and data; AMPL uses `.mod`, `.dat` and
  `.run` files, and Pyomo, JuMP, OR-Tools, PuLP, CVXPY and YALMIP play the
  same role in Python, Julia and MATLAB ({numref}`tab-09-01-modelling`).
- Read the status first ({numref}`tab-09-01-status`); then the optimal value,
  primal values, duals and slacks, converting duals to shadow prices with
  {numref}`tab-09-01-duals`.
- Every solver result is exact only up to feasibility, optimality,
  integrality and gap tolerances ({ref}`sec-09-01-tolerances`).
- SciPy, PuLP and CVXPY all reproduce the product-mix optimum
  $\vx^*=(6.429,4.286)$, $z^*=\$5142.86$ and the shadow prices
  $(78.57,\,2.86,\,0)$ ({ref}`sec-09-01-check`).

### Common mistakes

- **Using `res.x` without checking the status.** For an infeasible or
  unbounded model the arrays hold meaningless numbers, or nothing.
- **Forgetting that `linprog` minimizes.** A maximization must be passed as
  $\min(-\vc)^\T\vx$; the optimal value is then $-$`res.fun`, and the
  shadow prices are the negated marginals.
- **Mixing dual-sign conventions.** SciPy marginals, PuLP's `pi` and
  CVXPY's `dual_value` differ in sign for some combinations of sense and
  row type ({numref}`tab-09-01-duals`). Check a dual value by re-solving with a
  slightly larger right-hand side.
- **Comparing floats with `==`.** Use tolerances; a binary reported as
  $0.9999999$ is $1$.
- **Expecting the same optimal vertex from two solvers.** With alternative
  optima the solution reported depends on the solver; only the optimal
  value is unique.
- **Badly scaled models and huge big-$M$ constants**, which turn
  tolerances into real constraint violations.

## References

- Bynum, M. L., Hackebeil, G. A., Hart, W. E., Laird, C. D., Nicholson,
  B. L., Siirola, J. D., Watson, J.-P., & Woodruff, D. L. (2021).
  *Pyomo — Optimization modeling in Python* (3rd ed.). Springer. (Chapters
  1–3.)
- Diamond, S., & Boyd, S. (2016). CVXPY: A Python-embedded modeling
  language for convex optimization. *Journal of Machine Learning Research,
  17*(83), 1–5.
- Domahidi, A., Chu, E., & Boyd, S. (2013). ECOS: An SOCP solver for
  embedded systems. In *Proceedings of the European Control Conference*
  (pp. 3071–3076).
- Fourer, R., Gay, D. M., & Kernighan, B. W. (2003). *AMPL: A modeling
  language for mathematical programming* (2nd ed.). Duxbury Press. (Chapter 1:
  production models; Chapter 12: display commands, including duals and
  slacks; Chapter 13: command scripts.)
- Huangfu, Q., & Hall, J. A. J. (2018). Parallelizing the dual revised
  simplex method. *Mathematical Programming Computation, 10*(1), 119–142.
- Löfberg, J. (2004). YALMIP: A toolbox for modeling and optimization in
  MATLAB. In *Proceedings of the IEEE International Symposium on Computer
  Aided Control System Design* (pp. 284–289).
- Lubin, M., Dowson, O., Dias Garcia, J., Huchette, J., Legat, B., &
  Vielma, J. P. (2023). JuMP 1.0: Recent improvements to a modeling language
  for mathematical optimization. *Mathematical Programming Computation,
  15*, 581–589.
- Mitchell, S., O'Sullivan, M., & Dunning, I. (2011). *PuLP: A linear
  programming toolkit for Python*. Department of Engineering Science, The
  University of Auckland.
- O'Donoghue, B., Chu, E., Parikh, N., & Boyd, S. (2016). Conic optimization
  via operator splitting and homogeneous self-dual embedding. *Journal of
  Optimization Theory and Applications, 169*(3), 1042–1068.
- Stellato, B., Banjac, G., Goulart, P., Bemporad, A., & Boyd, S. (2020).
  OSQP: An operator splitting solver for quadratic programs. *Mathematical
  Programming Computation, 12*(4), 637–672.
- Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T.,
  Cournapeau, D., Burovski, E., Peterson, P., Weckesser, W., Bright, J.,
  van der Walt, S. J., Brett, M., Wilson, J., Millman, K. J., Mayorov, N.,
  Nelson, A. R. J., Jones, E., Kern, R., Larson, E., … SciPy 1.0
  Contributors. (2020). SciPy 1.0: Fundamental algorithms for scientific
  computing in Python. *Nature Methods, 17*(3), 261–272.
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapter 2: solving mathematical programming models;
  Section 3.5: the use of modeling languages.)
