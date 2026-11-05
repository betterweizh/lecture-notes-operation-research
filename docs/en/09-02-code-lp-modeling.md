# Modeling LPs and MILPs in Python

## Overview

This chapter writes the linear programs (LPs) and mixed-integer linear
programs (MILPs) of Parts III, IV, VII and VIII
in Python and lets a solver do
the arithmetic. It starts with the product-mix LP of
{prf:ref}`ex-03-01-product-mix` in the three styles of {doc}`09-01-code-setup-solvers`:
**matrix form** with `scipy.optimize.linprog`, **algebraic form** with named
variables and constraints in PuLP, and **vector form** in CVXPY. It then
reads the solver's dual values and slacks and checks them against the
shadow prices of {doc}`03-05-lp-duality` and the optimal tableaux of
{doc}`03-03-simplex-method` and {doc}`03-07-sensitivity-analysis`. Larger models
are built from data held in dictionaries and indexed by sets, as in an
algebraic modeling language: a blending model, a multi-period
production–inventory model and the transportation model. The second half
turns to **mixed-integer** models with `scipy.optimize.milp` and PuLP: the
facility-opening problem with fixed charges, an either–or constraint, and a
concave piecewise-linear cost with the $\lambda$-formulation, each compared
with its LP relaxation. The chapter ends with ways to check a model: reading
infeasible and unbounded statuses, measuring constraint violations, and
enumerating small instances by brute force. The models come from
{doc}`03-01-lp-model-building`, {doc}`08-01-examples-lp-modeling`,
{doc}`04-04-minimum-cost-flow`, {doc}`07-01-logical-constraints` and
{doc}`07-03-piecewise-linear-approximation`.

As everywhere in Part IX, the code is not run when the book is built. Copy
the blocks in order into a script or notebook and compare the output with
the `# Expected:` comments, which repeat the numbers of the worked examples.

:::{admonition} Learning objectives
:class: objectives
After completing this chapter, readers should be able to:

- write an LP in matrix form for `linprog`, in algebraic form with PuLP and
  in vector form with CVXPY, and obtain the same optimum from all three;
- read dual values, reduced costs and slacks from each interface, convert
  them to shadow prices and verify them by re-solving;
- build models from data stored in dictionaries, with variables and
  constraints indexed by sets;
- formulate and solve MILPs with binary variables in `milp` and PuLP, and
  compare formulations by the bounds of their LP relaxations; and
- check a model by its status, by the residuals of its constraints and by
  enumeration on a small instance.
:::

All blocks use the following imports.

```python
import itertools

import cvxpy as cp
import numpy as np
import pulp
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

np.set_printoptions(precision=3, suppress=True, linewidth=120)

# The CBC solver bundled with PuLP, without its log.
CBC = pulp.PULP_CBC_CMD(msg=False)
```

(sec-09-02-product-mix)=
## The product-mix LP three ways

### Matrix form with SciPy

`linprog` solves $\min\{\vc^\T\vx:\mA_{ub}\vx\le\vb_{ub},\
\mA_{eq}\vx=\vb_{eq},\ \vl\le\vx\le\vu\}$. The product-mix LP
{eq}`eq-03-01-product-mix` is a maximization with three $\le$ rows, so it is
passed as $\min(-\vc)^\T\vx$. A small wrapper hides the two sign changes,
in the objective and in the dual values ({numref}`tab-09-01-duals`), and stops
with an error when the LP has no optimal solution.

```python
def solve_max_lp(
    c: np.ndarray, A_ub: np.ndarray, b_ub: np.ndarray
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    """Solve max c^T x s.t. A_ub x <= b_ub, x >= 0 with HiGHS.

    Parameters
    ----------
    c : ndarray, shape (n,)
        Objective coefficients (contributions).
    A_ub : ndarray, shape (m, n)
        Constraint matrix of the <= rows.
    b_ub : ndarray, shape (m,)
        Right-hand sides (amounts available).

    Returns
    -------
    x : ndarray, shape (n,)
        Optimal solution.
    z : float
        Optimal value c^T x.
    pi : ndarray, shape (m,)
        Shadow prices dz*/db_i (nonnegative for a maximization).
    slack : ndarray, shape (m,)
        Unused amounts b_ub - A_ub x.
    """
    c, A_ub, b_ub = np.asarray(c), np.asarray(A_ub), np.asarray(b_ub)
    if A_ub.shape != (len(b_ub), len(c)):
        raise ValueError("A_ub must have shape (len(b_ub), len(c))")
    res = linprog(-c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
    if res.status != 0:
        raise ValueError(f"no optimal solution: {res.message}")
    return res.x, -res.fun, -res.ineqlin.marginals, res.slack


# Product mix: x1 juice, x2 cocktail glasses (hundreds of cases per week).
c = np.array([500.0, 450.0])
A = np.array(
    [
        [6.0, 5.0],
        [10.0, 20.0],
        [1.0, 0.0],
    ]
)
b = np.array([60.0, 150.0, 8.0])

x, z, pi, slack = solve_max_lp(c, A, b)
print("x* =", x, " z* =", round(z, 2))
print("pi =", pi, " slack =", slack)
# Expected:
# x* = [6.429 4.286]  z* = 5142.86
# pi = [78.571  2.857  0.   ]  slack = [0.    0.    1.571]
```

The matrix form is compact, but the meaning of each column and row lives
only in comments. That is fine for a model with two variables and
unmanageable for one with two thousand.

### Algebraic form with PuLP

An algebraic model separates the **data** (sets and parameters) from the
**model** (variables, objective, constraints), exactly as the `.mod` and
`.dat` files of AMPL do ({ref}`sec-09-01-modelling-languages`). In Python the
data are dictionaries keyed by the names of products and resources, and the
model is a function of the data.

```python
# Data of the product mix, keyed by names instead of positions:
# contribution ($ per hundred cases), resource use per hundred cases
# (machine h, storage m3), availability (h per week, m3) and the demand
# limit (hundred cases per week; none for cocktail glasses).
profit = {"juice": 500, "cocktail": 450}
use = {
    "machine": {"juice": 6, "cocktail": 5},
    "storage": {"juice": 10, "cocktail": 20},
}
avail = {"machine": 60, "storage": 150}
demand = {"juice": 8}


def build_product_mix(
    profit: dict, use: dict, avail: dict, demand: dict
) -> tuple[pulp.LpProblem, dict]:
    """Build the resource-allocation LP (eq-03-01-resource-allocation) in PuLP.

    Parameters
    ----------
    profit : dict
        Contribution per unit of each product.
    use : dict of dict
        use[r][p] = amount of resource r used per unit of product p.
    avail : dict
        Amount available of each resource.
    demand : dict
        Upper limit on the output of some products.

    Returns
    -------
    prob : pulp.LpProblem
        The model (not yet solved).
    make : dict
        The decision variables, one per product.
    """
    if set(use) != set(avail):
        raise ValueError("use and avail must have the same resources")
    prob = pulp.LpProblem("product_mix", pulp.LpMaximize)
    make = {
        p: pulp.LpVariable(f"make_{p}", lowBound=0, upBound=demand.get(p))
        for p in profit
    }
    prob += pulp.lpSum(profit[p] * make[p] for p in profit), "contribution"
    for r in avail:
        prob += (
            pulp.lpSum(use[r][p] * make[p] for p in profit) <= avail[r],
            r,
        )
    return prob, make


prob, make = build_product_mix(profit, use, avail, demand)
prob.solve(CBC)
print(pulp.LpStatus[prob.status], round(pulp.value(prob.objective), 2))
for p, var in make.items():
    print(f"{p:<9} {var.value():.3f}  reduced cost {var.dj:.3f}")
for name, con in prob.constraints.items():
    print(f"{name:<8} pi = {con.pi:.3f}  slack = {con.slack:.3f}")
# Expected:
# Optimal 5142.86
# juice     6.429  reduced cost -0.000
# cocktail  4.286  reduced cost -0.000
# machine  pi = 78.571  slack = -0.000
# storage  pi = 2.857  slack = -0.000
```

Here the demand limit is written as an upper **bound** on the variable
rather than as a constraint row, which is how solvers treat simple bounds
anyway; as a consequence it has no row dual. Its slack is
$8-6.429=1.571$, and its dual value would be zero, because it is not
binding.

PuLP can also print the model in the LP file format, which is the
quickest way to see whether the data were read as intended:

```python
print(prob)
# Expected (spacing may differ):
# product_mix:
# MAXIMIZE
# 450*make_cocktail + 500*make_juice + 0.0
# SUBJECT TO
# machine: 5 make_cocktail + 6 make_juice <= 60
#
# storage: 20 make_cocktail + 10 make_juice <= 150
#
# VARIABLES
# make_cocktail Continuous
# make_juice <= 8 Continuous
```

### Vector form with CVXPY

CVXPY (requires CVXPY) uses the same arrays as SciPy but keeps the
maximization and the $\le$ signs as they are. The constraint object returned
by `A @ x <= b` carries the dual values after the solve.

```python
x_cp = cp.Variable(2, nonneg=True)
rows = A @ x_cp <= b
mix = cp.Problem(cp.Maximize(c @ x_cp), [rows])
mix.solve()
print(mix.status, round(mix.value, 2))
print("x* =", np.round(x_cp.value, 3))
print("pi =", np.round(rows.dual_value, 3))
# Expected:
# optimal 5142.86
# x* = [6.429 4.286]
# pi = [78.571  2.857  0.   ]
```

All three interfaces give $\vx^*\approx(6.429,4.286)$ hundred cases per
week and $z^*\approx\$5142.86$, the optimum
$(\tfrac{45}{7},\tfrac{30}{7})$, $\tfrac{36000}{7}$ of
{prf:ref}`ex-03-03-simplex-run`.

(sec-09-02-duals)=
## Duals, shadow prices and slacks

### Matching the solver with the theory

The dual values $\vpi^*=(78.571,\ 2.857,\ 0)$ are the shadow prices
$\tfrac{550}{7}$ \$/h, $\tfrac{20}{7}$ \$/m³ and $0$ of
{prf:ref}`ex-03-05-product-mix-dual`, and also the negated entries of the
$-z$ row under the slack columns of the optimal tableau
{numref}`tab-03-03-pm-3`. Complementary slackness
({prf:ref}`thm-03-05-complementary-slackness`) can be read off the output:
the two binding rows have positive prices, the demand row with slack
$1.571$ has price zero, and both products, being produced, have reduced
cost zero. Strong duality says that the resources are worth exactly the
optimal contribution, $\vb^\T\vpi^*=z^*$:

```python
print(round(b @ pi, 2), round(c @ x, 2))
print(np.allclose(pi, [550 / 7, 20 / 7, 0.0]))
# Expected:
# 5142.86 5142.86
# True
```

### Checking a shadow price by re-solving

A shadow price predicts the change of the optimal value per extra unit of a
resource, but only within its range of validity
({prf:ref}`thm-03-05-shadow-price`). Re-solving the LP for several amounts of
machine time shows both the prediction and its limits. By
{ref}`sec-03-05-shadow-prices`, $\pi_1^*=\tfrac{550}{7}$ holds for machine
time between 37.5 and 65.5 h per week.

```python
def contribution_for_hours(hours: float) -> float:
    """Optimal contribution of the product mix with `hours` machine hours."""
    if hours < 0:
        raise ValueError("hours must be nonnegative")
    b_h = np.array([hours, 150.0, 8.0])
    return solve_max_lp(c, A, b_h)[1]


for hours in [30.0, 37.5, 50.0, 60.0, 61.0, 65.5, 70.0]:
    z_h = contribution_for_hours(hours)
    slope = contribution_for_hours(hours + 0.5) - z_h
    print(f"{hours:5.1f} h  z* = {z_h:8.2f}  next half hour: {slope:6.2f}")
# Expected:
#  30.0 h  z* =  2700.00  next half hour:  45.00
#  37.5 h  z* =  3375.00  next half hour:  39.29
#  50.0 h  z* =  4357.14  next half hour:  39.29
#  60.0 h  z* =  5142.86  next half hour:  39.29
#  61.0 h  z* =  5221.43  next half hour:  39.29
#  65.5 h  z* =  5575.00  next half hour:   0.00
#  70.0 h  z* =  5575.00  next half hour:   0.00
```

Between 37.5 and 65.5 h each half hour is worth $\tfrac12\cdot78.57=\$39.29$,
and one extra hour raises the contribution from \$5142.86 to \$5221.43, as
{ref}`sec-03-05-shadow-prices` predicts. Above 65.5 h the demand limit binds
and more machine time is worth nothing; below 37.5 h only cocktail glasses
are made, each hour yields $\tfrac{450}{5}=\$90$ ($\$45$ per half hour),
and the price changes. The optimal value is a concave, piecewise-linear
function of the resource, as {prf:ref}`thm-03-07-value-convex` states for the
minimization form.

### Duals of equality constraints and reduced costs

For a minimization in standard form, $\min\{\vc^\T\vx:\mA\vx=\vb,\
\vx\ge\vzero\}$, `linprog` returns the multipliers of the equations in
`res.eqlin.marginals` (free in sign) and the reduced costs
$\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$ of the variables at their lower bound
zero in `res.lower.marginals`. For the running example
{eq}`eq-03-07-running-lp` of {doc}`03-07-sensitivity-analysis` they should be
$\vpi=(2,-4,1)$ and $\bar\vc=(0,0,0,1,3,8)$, the $-z$ row of
{numref}`tab-03-07-optimal-tableau`.

```python
# Running example of chapter 3.7: min c^T x, A x = b, x >= 0.
c37 = np.array([3.0, 2.0, -3.0, -6.0, 10.0, -5.0])
A37 = np.array(
    [
        [1.0, 2.0, 0.0, 1.0, 0.0, -6.0],
        [0.0, 1.0, 1.0, 3.0, -2.0, -1.0],
        [1.0, 2.0, 1.0, 3.0, -1.0, -5.0],
    ]
)
b37 = np.array([11.0, 6.0, 13.0])

res37 = linprog(c37, A_eq=A37, b_eq=b37, bounds=(0, None), method="highs")
print("x* =", res37.x, " z* =", round(res37.fun, 6))
print("pi =", res37.eqlin.marginals)
print("reduced costs =", res37.lower.marginals)
# Expected:
# x* = [3. 4. 2. 0. 0. 0.]  z* = 11.0
# pi = [ 2. -4.  1.]
# reduced costs = [0. 0. 0. 1. 3. 8.]
```

The same model in PuLP returns the same numbers in `pi` and `dj`; for a
minimization they are again derivatives of the optimal value:

```python
lp37 = pulp.LpProblem("running_example", pulp.LpMinimize)
x37 = [pulp.LpVariable(f"x{j + 1}", lowBound=0) for j in range(6)]
lp37 += pulp.lpSum(c37[j] * x37[j] for j in range(6))
for i in range(3):
    lp37 += pulp.lpSum(A37[i, j] * x37[j] for j in range(6)) == b37[i]
lp37.solve(CBC)
print([round(con.pi, 3) for con in lp37.constraints.values()])
print([round(var.dj, 3) for var in x37])
# Expected:
# [2.0, -4.0, 1.0]
# [0.0, 0.0, 0.0, 1.0, 3.0, 8.0]
```

The sign of $\pi_2=-4$ is meaningful: raising $b_2$ *lowers* the minimum
cost. Simultaneous changes can be predicted the same way as long as the
basis stays optimal. Raising $b_1$ to 12 and $b_3$ to 14 uses 75 % of the
allowed increases ({prf:ref}`thm-03-07-hundred-percent`), so the predicted
optimum is $11+2\cdot1+1\cdot1=14$:

```python
res_new = linprog(
    c37, A_eq=A37, b_eq=[12.0, 6.0, 14.0], bounds=(0, None), method="highs"
)
print(round(res_new.fun, 6), round(11 + res37.eqlin.marginals @ [1, 0, 1], 6))
print("x* =", res_new.x)
# Expected:
# 14.0 14.0
# x* = [4. 4. 2. 0. 0. 0.]
```

Neither SciPy nor PuLP with CBC reports the **ranges** over which the duals
and the optimal basis remain valid (the cost and right-hand-side ranges of
{numref}`tab-03-07-ranges`). {ref}`sec-09-03-sensitivity` computes them from the
optimal basis with the formulas of {doc}`03-07-sensitivity-analysis`.

(sec-09-02-data-driven)=
## Data-driven models

### A blending model

The foundry of {prf:ref}`ex-08-01-foundry` must deliver 1000 kg of castings
with at least 0.45 % manganese and between 3.5 % and 6 % silicon. It can
melt three pig irons, whose content and cost (purchase plus \$9 melting per
1000 kg) are given per thousand kilograms, and add pure manganese at \$14
per kg. With $x_k$ the thousands of kilograms of pig iron $k$ and $x_{Mn}$
the kilograms of pure manganese, the model is the blending model of
{ref}`sec-03-01-blending` with a fixed total. Written from data, every
constraint is one line, and a new pig iron is one more dictionary entry.

```python
# Pig irons: kg of manganese and silicon, and cost in $, per 1000 kg.
pig_irons = {
    "A": {"manganese": 4.5, "silicon": 40.0, "cost": 40.0 + 9.0},
    "B": {"manganese": 5.5, "silicon": 13.0, "cost": 45.0 + 9.0},
    "C": {"manganese": 4.2, "silicon": 5.0, "cost": 25.0 + 9.0},
}
# Pure manganese ($ per kg), the order (kg) and the content limits (kg).
pure_mn_cost = 14.0
order_kg = 1000.0
mn_min = 0.0045 * order_kg
si_min, si_max = 0.035 * order_kg, 0.06 * order_kg

foundry = pulp.LpProblem("foundry", pulp.LpMinimize)
melt = {k: pulp.LpVariable(f"pig_{k}", lowBound=0) for k in pig_irons}
pure_mn = pulp.LpVariable("pure_manganese", lowBound=0)

foundry += (
    pulp.lpSum(d["cost"] * melt[k] for k, d in pig_irons.items())
    + pure_mn_cost * pure_mn
)
foundry += 1000 * pulp.lpSum(melt.values()) + pure_mn == order_kg, "weight"
manganese = pulp.lpSum(d["manganese"] * melt[k] for k, d in pig_irons.items())
silicon = pulp.lpSum(d["silicon"] * melt[k] for k, d in pig_irons.items())
foundry += manganese + pure_mn >= mn_min, "manganese"
foundry += silicon >= si_min, "silicon_min"
foundry += silicon <= si_max, "silicon_max"
foundry.solve(CBC)

cost = pulp.value(foundry.objective)
print(pulp.LpStatus[foundry.status])
for k, var in melt.items():
    print(f"pig iron {k}: {1000 * var.value():6.1f} kg")
print(f"pure manganese: {pure_mn.value():.1f} kg")
print(f"cost = ${cost:.2f}, profit = ${0.85 * order_kg - cost:.2f}")
for name, con in foundry.constraints.items():
    print(f"{name:<12} pi = {con.pi:7.3f}  slack = {con.slack:6.2f}")
# Expected:
# Optimal
# pig iron A:  849.2 kg
# pig iron B:   34.8 kg
# pig iron C:  116.0 kg
# pure manganese: 0.0 kg
# cost = $47.43, profit = $802.57
# weight       pi =  -0.024  slack =  -0.00
# manganese    pi =  13.457  slack =  -0.00
# silicon_min  pi =   0.313  slack =  -0.00
# silicon_max  pi =   0.000  slack =  25.00
```

The solution is $\vx^*=(\tfrac{366}{431},\tfrac{15}{431},\tfrac{50}{431},0)$
thousand kg at a cost of $\tfrac{20444}{431}\approx\$47.43$, as in
{prf:ref}`ex-08-01-foundry`. The duals explain it. The manganese requirement
costs \$13.46 per extra kilogram when it is met through the pig-iron mix,
less than the \$14 of pure manganese, so no pure manganese is bought. The
upper silicon limit has 25 kg of slack and price zero. For a $\ge$ row,
PuLP's slack (right-hand side minus left-hand side) is $\le0$, so a surplus
would appear with a minus sign.

### A multi-period model

Multi-period models are where indexing pays off. The production–inventory
problem of {prf:ref}`ex-08-01-production-inventory` has two products, three
months and the balance equations {eq}`eq-03-01-inventory-balance`, one per
product and month, with no initial stock. The variables are dictionaries
keyed by (product, month).

```python
PRODUCTS = ["A", "B"]
MONTHS = ["Jun", "Jul", "Aug"]
demand_pt = {
    ("A", "Jun"): 500,
    ("A", "Jul"): 5000,
    ("A", "Aug"): 750,
    ("B", "Jun"): 1000,
    ("B", "Jul"): 1200,
    ("B", "Aug"): 1200,
}
# Capacity (h per month), labor (h per unit), production cost ($ per
# unit) and holding cost ($ per unit per month).
capacity_h = {"Jun": 3000, "Jul": 3500, "Aug": 3000}
hours_per_unit = {"A": 0.75, "B": 1.0}
production_cost = {"A": 30.0, "B": 28.0}
holding_cost = {"A": 0.90, "B": 0.75}

plan = pulp.LpProblem("production_inventory", pulp.LpMinimize)
produce = pulp.LpVariable.dicts("p", (PRODUCTS, MONTHS), lowBound=0)
stock = pulp.LpVariable.dicts("I", (PRODUCTS, MONTHS), lowBound=0)

plan += pulp.lpSum(
    production_cost[k] * produce[k][t] + holding_cost[k] * stock[k][t]
    for k in PRODUCTS
    for t in MONTHS
)
for t in MONTHS:
    plan += (
        pulp.lpSum(hours_per_unit[k] * produce[k][t] for k in PRODUCTS)
        <= capacity_h[t],
        f"capacity_{t}",
    )
for k in PRODUCTS:
    for s, t in enumerate(MONTHS):
        # Stock at the end of the previous month (zero before June).
        previous = stock[k][MONTHS[s - 1]] if s > 0 else 0
        plan += (
            stock[k][t] == previous + produce[k][t] - demand_pt[k, t],
            f"balance_{k}_{t}",
        )
plan.solve(CBC)

print(pulp.LpStatus[plan.status], round(pulp.value(plan.objective), 2))
for k in PRODUCTS:
    made = [round(produce[k][t].value(), 2) for t in MONTHS]
    held = [round(stock[k][t].value(), 2) for t in MONTHS]
    print(f"{k}: produce {made}  stock {held}")
for t in MONTHS:
    con = plan.constraints[f"capacity_{t}"]
    print(f"capacity {t}: pi = {con.pi:.2f} $/h, slack = {con.slack:.1f} h")
# Expected:
# Optimal 283900.0
# A: produce [833.33, 4666.67, 750.0]  stock [333.33, 0.0, 0.0]
# B: produce [2200.0, 0.0, 1200.0]  stock [1200.0, 0.0, 0.0]
# capacity Jun: pi = 0.00 $/h, slack = 175.0 h
# capacity Jul: pi = -1.20 $/h, slack = -0.0 h
# capacity Aug: pi = 0.00 $/h, slack = 1237.5 h
```

The plan and the cost $z^*=\$283\,900$ are those of
{prf:ref}`ex-08-01-production-inventory`. July is the bottleneck. Its shadow
price, $-\$1.20$ per hour, is the saving from one more July hour: it would
let $\tfrac43$ units of A be made in July instead of June and save
$\tfrac43\cdot\$0.90=\$1.20$ of holding cost. The extra hour would *not* be
used for B: moving one unit of B (one hour) from June to July saves only
\$0.75 of holding cost. The two other months have slack and price zero.

(sec-09-02-transportation)=
## The transportation model

The transportation problem of {prf:ref}`ex-04-04-transportation` ships 40 and
50 t per week from plants A and B to three markets that need 20, 25 and
45 t. Its supply and demand rows are the balance equations
{eq}`eq-04-04-transportation`; the arcs are the keys of the cost dictionary.

```python
# Supplies and demands (t per week) and shipping costs ($ per t).
supply = {"A": 40, "B": 50}
market_demand = {1: 20, 2: 25, 3: 45}
ship_cost = {
    ("A", 1): 4,
    ("A", 2): 6,
    ("A", 3): 9,
    ("B", 1): 2,
    ("B", 2): 7,
    ("B", 3): 5,
}

transport = pulp.LpProblem("transportation", pulp.LpMinimize)
ship = {
    (i, j): pulp.LpVariable(f"x_{i}{j}", lowBound=0) for (i, j) in ship_cost
}
transport += pulp.lpSum(ship_cost[a] * ship[a] for a in ship_cost)
for i in supply:
    transport += (
        pulp.lpSum(ship[i, j] for j in market_demand) == supply[i],
        f"supply_{i}",
    )
for j in market_demand:
    transport += (
        pulp.lpSum(ship[i, j] for i in supply) == market_demand[j],
        f"demand_{j}",
    )
transport.solve(CBC)

print(pulp.LpStatus[transport.status], pulp.value(transport.objective))
print({a: var.value() for a, var in ship.items() if var.value() > 0})
# Expected:
# Optimal 445.0
# {('A', 1): 15.0, ('A', 2): 25.0, ('B', 1): 5.0, ('B', 3): 45.0}
```

The plan costs \$445 per week, as found by successive shortest paths in
{ref}`sec-04-04-ssp`, and it is integral although nothing required it: the
constraint matrix is totally unimodular ({prf:ref}`thm-04-04-integrality`).

**Duals as node potentials.** With $u_i$ the dual of the supply row of
plant $i$ and $v_j$ that of the demand row of market $j$, the reduced cost
of arc $(i,j)$ is $c_{ij}-u_i-v_j$, zero on every arc used. One of the five
equations is redundant (total supply equals total demand), so the duals are
determined only up to adding a constant to all $u_i$ and subtracting it from
all $v_j$, and different solvers report different, equally correct, values.
The node potentials of {ref}`sec-04-04-ssp` are $\pi_i=u_i+k$ for the plants
and $\pi_j=k-v_j$ for the markets, with $k=v_3$ so that $\pi_3=0$:

```python
u = {i: transport.constraints[f"supply_{i}"].pi for i in supply}
v = {j: transport.constraints[f"demand_{j}"].pi for j in market_demand}
k = v[3]
potential = {i: u[i] + k for i in supply}
potential.update({j: k - v[j] for j in market_demand})
print(potential)
reduced = {
    (i, j): ship_cost[i, j] - potential[i] + potential[j]
    for (i, j) in ship_cost
}
for arc, value in reduced.items():
    print(arc, value)
# Expected:
# {'A': 7.0, 'B': 5.0, 1: 3.0, 2: 1.0, 3: 0.0}
# ('A', 1) 0.0
# ('A', 2) 0.0
# ('A', 3) 2.0
# ('B', 1) 0.0
# ('B', 2) 3.0
# ('B', 3) 0.0
```

These are the potentials $\pi_A=7$, $\pi_B=5$, $\pi_1=3$, $\pi_2=1$,
$\pi_3=0$ of {ref}`sec-04-04-ssp`, and the two unused arcs have the positive
reduced costs $\bar c_{A3}=2$ and $\bar c_{B2}=3$ that prove the plan
optimal and unique ({prf:ref}`thm-04-04-reduced-cost`).

For `linprog` the same model needs the constraint matrix explicitly: one
column per arc, a 1 in the row of its plant and in the row of its market.
Building it from the arc list is the incidence-matrix construction of
{doc}`04-04-minimum-cost-flow`, and the result is the same:

```python
arcs = list(ship_cost)
plants, markets = list(supply), list(market_demand)
A_tr = np.zeros((len(plants) + len(markets), len(arcs)))
for col, (i, j) in enumerate(arcs):
    A_tr[plants.index(i), col] = 1.0
    A_tr[len(plants) + markets.index(j), col] = 1.0
b_tr = np.array(
    [supply[i] for i in plants] + [market_demand[j] for j in markets]
)
res_tr = linprog(
    [ship_cost[a] for a in arcs], A_eq=A_tr, b_eq=b_tr, method="highs"
)
print(A_tr)
print(res_tr.fun, res_tr.x)
# Expected:
# [[1. 1. 1. 0. 0. 0.]
#  [0. 0. 0. 1. 1. 1.]
#  [1. 0. 0. 1. 0. 0.]
#  [0. 1. 0. 0. 1. 0.]
#  [0. 0. 1. 0. 0. 1.]]
# 445.0 [15. 25.  0.  5.  0. 45.]
```

(sec-09-02-milp)=
## Mixed-integer models

### Binaries in `scipy.optimize.milp`

`milp` solves $\min\{\vc^\T\vx:\vb_l\le\mA\vx\le\vb_u,\ \vl\le\vx\le\vu,\
x_j\in\mathbb Z\ \text{for}\ j\in J\}$ with the HiGHS branch-and-cut code.
Constraints are passed as `LinearConstraint(A, lb, ub)` (an equation has
`lb = ub`, a $\le$ row `lb = -inf`), bounds as `Bounds(l, u)`, and the
integer variables are marked by an `integrality` vector with 1 for integer
and 0 for continuous; a binary is an integer variable with bounds $[0,1]$.

The facility-opening problem of {prf:ref}`ex-07-01-facility` has four
candidate sites with fixed costs $F_i$ (\$ thousand per year) and
capacities $u_i$ (thousand tonnes per year), and five districts with waste
$d_j$; $c_{ij}$ is in dollars per tonne, hence in \$ thousand per thousand
tonnes. The variables are the 20 flows $x_{ij}$, stored row by row
($x_{11},\dots,x_{15},x_{21},\dots$), followed by the four binaries
$\delta_i$. The function below builds the three formulations compared in
{prf:ref}`ex-07-01-facility` and solves either the MILP or its LP relaxation;
it can also fix the binaries, which is used for enumeration below.

```python
# Facility data (tab-07-01-facility).
fac_cost = np.array(
    [
        [4.0, 6.0, 9.0, 10.0, 8.0],
        [7.0, 3.0, 5.0, 9.0, 10.0],
        [10.0, 8.0, 4.0, 3.0, 6.0],
        [9.0, 10.0, 7.0, 5.0, 3.0],
    ]
)
# Fixed costs ($ thousand per year), capacities and waste (thousand t per
# year).
fixed = np.array([300.0, 250.0, 280.0, 220.0])
cap = np.array([90.0, 70.0, 80.0, 60.0])
waste = np.array([30.0, 25.0, 35.0, 20.0, 30.0])


def solve_facility(
    link: str = "capacity", relax: bool = False, open_sites=None
):
    """Solve the facility-opening MILP of ex-07-01-facility with milp.

    Parameters
    ----------
    link : {"capacity", "loose", "disaggregated"}
        How flows are linked to the binaries: sum_j x_ij <= u_i delta_i
        ("capacity"); sum_j x_ij <= 140 delta_i plus sum_j x_ij <= u_i
        ("loose"); or the capacity link plus x_ij <= min(d_j, u_i) delta_i
        ("disaggregated").
    relax : bool
        Solve the LP relaxation (delta_i continuous in [0, 1]).
    open_sites : sequence of int, optional
        If given, fix delta_i = 1 for these sites (1-based) and 0 otherwise.

    Returns
    -------
    scipy.optimize.OptimizeResult
        Result of milp; x[:20] are the flows, x[20:] the binaries.
    """
    if link not in ("capacity", "loose", "disaggregated"):
        raise ValueError(f"unknown link {link!r}")
    n_sites, n_districts = fac_cost.shape
    n_flow = n_sites * n_districts
    n_var = n_flow + n_sites
    objective = np.concatenate([fac_cost.ravel(), fixed])
    rows, lower, upper = [], [], []

    def flow(i, j):
        """Column of x_ij."""
        return i * n_districts + j

    def add_row(row, lo, hi):
        """Append the constraint lo <= row @ x <= hi."""
        rows.append(row)
        lower.append(lo)
        upper.append(hi)

    # Every district's waste is collected: sum_i x_ij = d_j.
    for j in range(n_districts):
        row = np.zeros(n_var)
        row[[flow(i, j) for i in range(n_sites)]] = 1.0
        add_row(row, waste[j], waste[j])
    # Link of the flows into site i with its binary.
    big_m = 140.0 if link == "loose" else None
    for i in range(n_sites):
        row = np.zeros(n_var)
        row[flow(i, 0) : flow(i, 0) + n_districts] = 1.0
        row[n_flow + i] = -(big_m if big_m else cap[i])
        add_row(row, -np.inf, 0.0)
        if link == "loose":
            row = np.zeros(n_var)
            row[flow(i, 0) : flow(i, 0) + n_districts] = 1.0
            add_row(row, -np.inf, cap[i])
    if link == "disaggregated":
        for i in range(n_sites):
            for j in range(n_districts):
                row = np.zeros(n_var)
                row[flow(i, j)] = 1.0
                row[n_flow + i] = -min(waste[j], cap[i])
                add_row(row, -np.inf, 0.0)

    lb, ub = np.zeros(n_var), np.full(n_var, np.inf)
    ub[n_flow:] = 1.0
    if open_sites is not None:
        if any(not 1 <= s <= n_sites for s in open_sites):
            raise ValueError(f"sites are numbered 1, ..., {n_sites}")
        chosen = np.zeros(n_sites)
        chosen[[s - 1 for s in open_sites]] = 1.0
        lb[n_flow:] = ub[n_flow:] = chosen
    integrality = np.zeros(n_var)
    if not relax:
        integrality[n_flow:] = 1
    return milp(
        objective,
        constraints=LinearConstraint(np.array(rows), lower, upper),
        integrality=integrality,
        bounds=Bounds(lb, ub),
    )


best = solve_facility()
print(best.status, round(best.fun, 2))
print("open sites:", [i + 1 for i in range(4) if best.x[20 + i] > 0.5])
# Adding 0.0 turns -0.0 into 0.0 for printing.
print(np.round(best.x[:20].reshape(4, 5), 3) + 0.0)
# Expected:
# 0 1200.0
# open sites: [2, 3]
# [[ 0.  0.  0.  0.  0.]
#  [30. 25.  5.  0.  0.]
#  [ 0.  0. 30. 20. 30.]
#  [ 0.  0.  0.  0.  0.]]
```

The optimum opens sites 2 and 3 at \$1 200 thousand per year: \$530
thousand of fixed cost and \$670 thousand of transport. (The flow of
district 3 may be split differently between the two open sites by another
solver version, at the same cost.) The LP relaxations of the three
formulations give the bounds of {prf:ref}`ex-07-01-facility`; the tighter the
link, the closer the bound to the integer optimum:

```python
for link in ["loose", "capacity", "disaggregated"]:
    bound = solve_facility(link, relax=True).fun
    gap = (best.fun - bound) / best.fun
    print(f"{link:<14} LP bound {bound:8.2f}   gap {100 * gap:5.1f} %")
# Expected:
# loose          LP bound   751.07   gap  37.4 %
# capacity       LP bound   976.79   gap  18.6 %
# disaggregated  LP bound  1168.46   gap   2.6 %
```

### The same MILP in PuLP

In PuLP a binary is `cat="Binary"`, and the model reads like
{prf:ref}`ex-07-01-facility`; the solver is CBC instead of HiGHS.

```python
SITES = range(1, 5)
DISTRICTS = range(1, 6)
facility = pulp.LpProblem("facility", pulp.LpMinimize)
send = pulp.LpVariable.dicts("x", (SITES, DISTRICTS), lowBound=0)
is_open = pulp.LpVariable.dicts("delta", SITES, cat="Binary")
facility += pulp.lpSum(fixed[i - 1] * is_open[i] for i in SITES) + pulp.lpSum(
    fac_cost[i - 1, j - 1] * send[i][j] for i in SITES for j in DISTRICTS
)
for j in DISTRICTS:
    facility += pulp.lpSum(send[i][j] for i in SITES) == waste[j - 1]
for i in SITES:
    facility += (
        pulp.lpSum(send[i][j] for j in DISTRICTS) <= cap[i - 1] * is_open[i]
    )
facility.solve(CBC)
print(pulp.LpStatus[facility.status], pulp.value(facility.objective))
print("open sites:", [i for i in SITES if is_open[i].value() > 0.5])
# Expected:
# Optimal 1200.0
# open sites: [2, 3]
```

### Checking by enumeration

With four sites there are only $2^4=16$ choices of open sites, so the MILP
can be checked by fixing the binaries and solving the remaining transport
LP for each choice, which reproduces the table of {prf:ref}`ex-07-01-facility`.

```python
totals = {}
for k in range(1, 5):
    for sites in itertools.combinations(range(1, 5), k):
        res_s = solve_facility(open_sites=sites)
        totals[sites] = res_s.fun if res_s.status == 0 else None
for sites, total in totals.items():
    if total is None:
        print(sites, "infeasible")
    else:
        print(sites, round(total, 2))
print("best:", min((t, s) for s, t in totals.items() if t is not None))
# Expected:
# (1,) infeasible
# (2,) infeasible
# (3,) infeasible
# (4,) infeasible
# (1, 2) 1350.0
# (1, 3) 1240.0
# (1, 4) 1275.0
# (2, 3) 1200.0
# (2, 4) infeasible
# (3, 4) 1260.0
# (1, 2, 3) 1410.0
# (1, 2, 4) 1330.0
# (1, 3, 4) 1360.0
# (2, 3, 4) 1325.0
# (1, 2, 3, 4) 1535.0
# best: (1200.0, (2, 3))
```

Single sites and the pair $\{2,4\}$ (130 thousand tonnes of capacity for
140 of waste) are infeasible, and $\{2,3\}$ is the cheapest choice.
Enumeration grows as $2^n$ and is useless for real instances, but on a
small one it is the most convincing test of a formulation.

### An either–or constraint

{prf:ref}`ex-07-01-either-or` maximizes $-x_1+2x_2$ over the box $[0,4]^2$
when $\vx$ must lie in the triangle $x_1+x_2\le2$ or in the square
$x_1\ge3$, $x_2\ge3$. With $\delta=1$ selecting the triangle, the big-$M$
rows are $x_1+x_2-2\le M_T(1-\delta)$, $3-x_1\le M_1\delta$,
$3-x_2\le M_2\delta$. The integer optimum does not depend on the constants
(as long as they are valid), but the LP relaxation does:

```python
def solve_either_or(m_t: float, m_1: float, m_2: float, relax: bool):
    """Either-or model of ex-07-01-either-or; variables (x1, x2, delta).

    Returns the optimal value of max -x1 + 2 x2 and the solution.
    """
    if min(m_t, m_1, m_2) < 0:
        raise ValueError("big-M constants must be nonnegative")
    # Rows: x1 + x2 + M_T delta <= 2 + M_T; -x1 - M_1 delta <= -3;
    # -x2 - M_2 delta <= -3.
    rows = np.array(
        [
            [1.0, 1.0, m_t],
            [-1.0, 0.0, -m_1],
            [0.0, -1.0, -m_2],
        ]
    )
    upper = np.array([2.0 + m_t, -3.0, -3.0])
    res = milp(
        [1.0, -2.0, 0.0],
        constraints=LinearConstraint(rows, -np.inf, upper),
        integrality=[0, 0, 0 if relax else 1],
        bounds=Bounds([0, 0, 0], [4, 4, 1]),
    )
    if res.status != 0:
        # For example infeasible when the constants are too small.
        raise ValueError(f"no optimal solution: {res.message}")
    return -res.fun, res.x


for m_values in [(6.0, 3.0, 3.0), (10.0, 10.0, 10.0)]:
    z_int, x_int = solve_either_or(*m_values, relax=False)
    z_lp, x_lp = solve_either_or(*m_values, relax=True)
    print(m_values, "MILP", z_int, x_int, " LP", round(z_lp, 3))
# Expected:
# (6.0, 3.0, 3.0) MILP 5.0 [3. 4. 0.]  LP 6.0
# (10.0, 10.0, 10.0) MILP 5.0 [3. 4. 0.]  LP 8.0
```

Both MILPs return $\vx^*=(3,4)$, $\delta^*=0$ with value 5. The tight
constants $M_T=6$, $M_1=M_2=3$ from {eq}`eq-07-01-m-from-bounds` give the
bound 6, the loose $M=10$ only 8, as {numref}`fig-07-01-either-or` shows.

### A concave piecewise-linear cost

The two treatment plants of {prf:ref}`ex-07-03-economies` must provide
100 ML/d of capacity, at most 80 ML/d each, with concave costs
$C_1(q)=q^{0.6}$ and $C_2(q)=0.8q^{0.6}+0.06q$ in \$ million per year.
With breakpoints $a_i=0,20,40,60,80$ ML/d, each cost is replaced by its
interpolant in the convex-combination form {eq}`eq-07-03-sos2`, and because the
costs are concave, adjacency of the weights must be enforced by the segment
binaries of {eq}`eq-07-03-lambda`: $\alpha_i\le\delta_{i-1}+\delta_i$ and
$\sum_k\delta_k=1$. The function below builds the model with or without the
binaries.

```python
# Breakpoints and required capacity (ML/d); plant costs ($ million per
# year) as functions of the capacity q.
breakpoints = [0.0, 20.0, 40.0, 60.0, 80.0]
plant_cost = {
    1: lambda q: q**0.6,
    2: lambda q: 0.8 * q**0.6 + 0.06 * q,
}
required = 100.0


def solve_plants(adjacency: bool) -> tuple[float, dict]:
    """Lambda-formulation of ex-07-03-economies in PuLP.

    Parameters
    ----------
    adjacency : bool
        Add the segment binaries (eq-07-03-lambda); without them the model is
        an LP that may combine non-adjacent breakpoints.

    Returns
    -------
    cost : float
        Optimal value of the model ($ million per year).
    q : dict
        Capacity of each plant (ML/d).
    """
    n = len(breakpoints)
    prob = pulp.LpProblem("plants", pulp.LpMinimize)
    alpha = pulp.LpVariable.dicts("alpha", (plant_cost, range(n)), lowBound=0)
    prob += pulp.lpSum(
        plant_cost[p](breakpoints[i]) * alpha[p][i]
        for p in plant_cost
        for i in range(n)
    )
    q = {
        p: pulp.lpSum(breakpoints[i] * alpha[p][i] for i in range(n))
        for p in plant_cost
    }
    prob += pulp.lpSum(q.values()) == required
    for p in plant_cost:
        prob += pulp.lpSum(alpha[p][i] for i in range(n)) == 1
        if adjacency:
            # delta[p][k] = 1 if plant p's capacity lies in segment k.
            delta = pulp.LpVariable.dicts(
                f"delta{p}", range(n - 1), cat="Binary"
            )
            prob += pulp.lpSum(delta.values()) == 1
            for i in range(n):
                neighbors = [delta[k] for k in (i - 1, i) if 0 <= k < n - 1]
                prob += alpha[p][i] <= pulp.lpSum(neighbors)
    prob.solve(CBC)
    if prob.status != pulp.LpStatusOptimal:
        raise ValueError(pulp.LpStatus[prob.status])
    return pulp.value(prob.objective), {p: q[p].value() for p in q}


for adjacency in [False, True]:
    total, q = solve_plants(adjacency)
    print(f"adjacency={adjacency}: cost {total:.2f}, q = {q}")
exact = plant_cost[1](80.0) + plant_cost[2](20.0)
print(f"true cost of q1 = 80, q2 = 20: {exact:.2f}")
# Expected:
# adjacency=False: cost 17.84, q = {1: 80.0, 2: 20.0}
# adjacency=True: cost 19.89, q = {1: 80.0, 2: 20.0}
# true cost of q1 = 80, q2 = 20: 19.89
```

Both models report $q_1=80$, $q_2=20$ ML/d, but only the MILP prices it
correctly. The LP without adjacency writes $q_2=20$ as $0.75$ of the
breakpoint 0 plus $0.25$ of the breakpoint 80 and charges the chord value
$0.25\cdot15.89$ instead of $C_2(20)=6.03$; its \$17.84 million is the cost
of the convex envelope, a lower bound and not a cost. With the binaries
the model returns \$19.89 million, the exact optimum
({prf:ref}`ex-07-03-economies`). Some solvers accept the adjacency condition
directly as an SOS2 constraint (Gurobi, CPLEX, Xpress through their own
interfaces); the CBC command-line interface of PuLP does not pass SOS
information to CBC, so with PuLP the binaries of {eq}`eq-07-03-lambda` are the
safe choice.

(sec-09-02-checking)=
## Checking a model

### Status first

{doc}`09-01-code-setup-solvers` lists the status codes. Two models of
{doc}`03-07-sensitivity-analysis` show the two failures. With $b_1=15$ the
running example is infeasible ({prf:ref}`ex-03-07-rhs`), and with the new
activity $\mA_{\cdot7}=(3,-1,1)^\T$, $c_7=4$ it is unbounded
({prf:ref}`ex-03-07-new-activity`):

```python
res_inf = linprog(
    c37, A_eq=A37, b_eq=[15.0, 6.0, 13.0], bounds=(0, None), method="highs"
)
A_new = np.column_stack([A37, [3.0, -1.0, 1.0]])
c_new = np.append(c37, 4.0)
res_unb = linprog(
    c_new, A_eq=A_new, b_eq=b37, bounds=(0, None), method="highs"
)
print(res_inf.status, res_inf.message.split(".")[0])
print(res_unb.status, res_unb.message.split(".")[0])
# Expected:
# 2 The problem is infeasible
# 3 The problem is unbounded
```

An infeasible model is debugged by relaxing constraints one family at a
time (or by adding nonnegative slack variables with a large cost to every
row and seeing which slacks become positive); an unbounded one by bounding
every variable with a large box and seeing which variables run to the box.

### Residuals

A reported optimum can be checked independently of the solver by
evaluating the constraints. The function below returns the largest
violation of $\mA_{ub}\vx\le\vb_{ub}$, $\mA_{eq}\vx=\vb_{eq}$ and
$\vx\ge\vzero$; anything above the feasibility tolerance of about $10^{-7}$
signals an error in the model or the data handling.

```python
def max_violation(x, A_ub=None, b_ub=None, A_eq=None, b_eq=None) -> float:
    """Largest violation of A_ub x <= b_ub, A_eq x = b_eq and x >= 0.

    Parameters
    ----------
    x : array_like, shape (n,)
        Candidate solution.
    A_ub, b_ub, A_eq, b_eq : array_like, optional
        Constraint data as for linprog.

    Returns
    -------
    float
        The maximum violation (0.0 if x is feasible).
    """
    x = np.asarray(x, dtype=float)
    worst = max(0.0, -x.min())
    if A_ub is not None:
        if b_ub is None:
            raise ValueError("b_ub is required with A_ub")
        worst = max(worst, (np.asarray(A_ub) @ x - b_ub).max())
    if A_eq is not None:
        if b_eq is None:
            raise ValueError("b_eq is required with A_eq")
        worst = max(worst, np.abs(np.asarray(A_eq) @ x - b_eq).max())
    return float(worst)


print(max_violation(x, A_ub=A, b_ub=b) < 1e-9)
print(max_violation(res37.x, A_eq=A37, b_eq=b37) < 1e-9)
print(max_violation(res_tr.x, A_eq=A_tr, b_eq=b_tr) < 1e-9)
print(max_violation([8.0, 4.0], A_ub=A, b_ub=b))
# Expected:
# True
# True
# True
# 10.0
```

The last line evaluates the infeasible plan $(8,4)$ of the product mix,
which needs 68 machine hours and 160 m³: the storage row is violated by
10 m³ and the machine row by 8 h, and the larger violation is reported.

### Enumeration of small instances

Rounding an LP optimum is not the same as solving the integer program
({ref}`sec-03-01-when-lp-fails`). For the product mix with whole hundreds of
cases, $x_1\in\{0,\dots,8\}$ and $x_2\in\{0,\dots,7\}$ (by the storage row),
so all 72 integer points can be tried; 54 of them are feasible:

```python
feasible_points = [
    (x1, x2)
    for x1, x2 in itertools.product(range(9), range(8))
    if np.all(A @ np.array([x1, x2]) <= b)
]
best_point = max(feasible_points, key=lambda p: c @ np.array(p))
print(len(feasible_points), best_point, c @ np.array(best_point))
print("rounded down:", (6, 4), c @ np.array([6, 4]))

res_int = milp(
    -c,
    constraints=LinearConstraint(A, -np.inf, b),
    integrality=np.ones(2),
    bounds=Bounds(0, np.inf),
)
print("milp:", res_int.x, -res_int.fun)
# Expected:
# 54 (8, 2) 4900.0
# rounded down: (6, 4) 4800.0
# milp: [8. 2.] 4900.0
```

The best integer plan $(8,2)$ earns \$4900 per week, \$100 more than the
rounded LP optimum $(6,4)$, as noted in {ref}`sec-03-01-when-lp-fails`; `milp`
agrees with the enumeration. A model that passes such a test on a small
instance, and whose duals make economic sense, is ready for real data.

## Summary and Common Mistakes

### Summary

- `linprog` takes the matrix form; PuLP builds the same LP from named
  variables and data dictionaries; CVXPY writes it in vector form. All three
  give the product-mix optimum $\vx^*\approx(6.429,4.286)$,
  $z^*\approx\$5142.86$ ({ref}`sec-09-02-product-mix`).
- Dual values are shadow prices after the sign conventions of
  {numref}`tab-09-01-duals` are applied: $(78.57,2.86,0)$ for the product mix
  and $\vpi=(2,-4,1)$ with reduced costs $(0,0,0,1,3,8)$ for the running
  example of {doc}`03-07-sensitivity-analysis`; re-solving confirms them within
  their ranges ({ref}`sec-09-02-duals`).
- Indexed models (blending, production–inventory, transportation) are
  written once and filled from data; their duals explain the solutions
  (\$13.46 per kg of manganese, \$1.20 per July hour, node potentials
  $(7,5,3,1,0)$) ({ref}`sec-09-02-data-driven`, {ref}`sec-09-02-transportation`).
- `milp` and PuLP solve MILPs with binaries; the facility problem costs
  \$1 200 thousand per year with sites 2 and 3, and its LP bounds 751.07,
  976.79 and 1168.46 measure the strength of three formulations; the
  $\lambda$-formulation with binaries prices concave costs correctly
  (\$19.89 million against \$17.84 million without adjacency)
  ({ref}`sec-09-02-milp`).
- Check every model by its status, its residuals and, on a small instance,
  by enumeration ({ref}`sec-09-02-checking`).

### Common mistakes

- **Forgetting the sign changes of `linprog`** for a maximization and for
  $\ge$ rows; wrap them in one function, as `solve_max_lp` does.
- **Reading duals of a degenerate or non-unique optimum as unique.** The
  transportation duals are determined only up to a constant; normalize them
  before comparing.
- **Using a shadow price outside its range.** Re-solve to check, as in
  {ref}`sec-09-02-duals`.
- **Hard-coding indices.** Build constraints from the data dictionaries;
  a new product, month or site should not require editing the model.
- **Off-by-one errors in balance equations.** State the initial stock
  explicitly (here zero) and test the model on a period with known data.
- **Treating an LP relaxation value as a cost.** For nonconvex
  piecewise-linear costs the relaxation charges the convex envelope; only
  the MILP with adjacency gives true costs.
- **Trusting `prob.sos2` in PuLP with the bundled CBC**: the sets are not
  passed to the solver, and the result is the LP without adjacency.

## References

- Diamond, S., & Boyd, S. (2016). CVXPY: A Python-embedded modeling
  language for convex optimization. *Journal of Machine Learning Research,
  17*(83), 1–5.
- Hillier, F. S., & Lieberman, G. J. (2021). *Introduction to operations
  research* (11th ed.). McGraw-Hill. (Chapter 3: LP models and their
  formulation in spreadsheets and modeling languages; Chapter 7: sensitivity
  output; Chapter 12: integer programming formulations.)
- Mitchell, S., O'Sullivan, M., & Dunning, I. (2011). *PuLP: A linear
  programming toolkit for Python*. Department of Engineering Science, The
  University of Auckland.
- Taha, H. A. (2017). *Operations research: An introduction* (10th ed.).
  Pearson. (Chapter 2: LP applications; Chapter 5: transportation model;
  Chapter 9: integer linear programming.)
- Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T.,
  Cournapeau, D., Burovski, E., Peterson, P., Weckesser, W., Bright, J.,
  van der Walt, S. J., Brett, M., Wilson, J., Millman, K. J., Mayorov, N.,
  Nelson, A. R. J., Jones, E., Kern, R., Larson, E., … SciPy 1.0
  Contributors. (2020). SciPy 1.0: Fundamental algorithms for scientific
  computing in Python. *Nature Methods, 17*(3), 261–272.
- Williams, H. P. (2013). *Model building in mathematical programming*
  (5th ed.). Wiley. (Chapters 3–5: building LP models; Chapter 6:
  interpreting solutions and dual values; Chapters 9–10: integer programming
  formulations.)
