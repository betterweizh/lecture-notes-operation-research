# 用 Python 建立 LP 与 MILP 模型

## 概述

本章用 Python 写出第三、四、七、八部分中的线性规划问题（linear program, LP）和混合整数线性规划（mixed-integer linear program, MILP），把算术运算交给求解器（solver）。首先以{doc}`09-01-code-setup-solvers`中的三种风格写出{prf:ref}`ex-03-01-product-mix`的产品组合（product mix）LP：用 `scipy.optimize.linprog` 的**矩阵形式**（matrix form），在 PuLP 中使用具名变量与具名约束的**代数形式**（algebraic form），以及 CVXPY 中的**向量形式**（vector form）。随后读取求解器给出的对偶值与松弛量，并与{doc}`03-05-lp-duality`中的影子价格（shadow price）以及{doc}`03-03-simplex-method`和{doc}`03-07-sensitivity-analysis`中的最优单纯形表相对照。较大的模型则像代数建模语言（algebraic modeling language）那样，由存放在字典中、按集合索引的数据构建：配料模型、多周期生产–库存模型（production–inventory model）和运输模型。本章后半部分转向**混合整数**模型，使用 `scipy.optimize.milp` 和 PuLP：带固定费用（fixed charge）的设施选址问题、一个择一约束（either–or constraint），以及用 $\lambda$ 形式表示的凹分段线性费用，每个模型都与其线性规划松弛（LP relaxation）作比较。本章最后介绍检查模型的方法：解读不可行与无界的求解状态（status）、度量约束违反量，以及对小规模实例作穷举。这些模型取自{doc}`03-01-lp-model-building`、{doc}`08-01-examples-lp-modeling`、{doc}`04-04-minimum-cost-flow`、{doc}`07-01-logical-constraints`和{doc}`07-03-piecewise-linear-approximation`。

与第九部分其他各章一样，构建本书时不运行这些代码。请把各代码块按顺序复制到脚本或笔记本中，并将输出与 `# Expected:` 注释相比较；这些注释重复了例题中的数值。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 把一个 LP 分别写成 `linprog` 所用的矩阵形式、PuLP 中的代数形式和 CVXPY 中的向量形式，并从三者得到相同的最优解；
- 从每种接口读取对偶值、检验数（reduced cost）和松弛量，把它们换算为影子价格，并通过重新求解加以验证；
- 由存放在字典中的数据构建模型，变量与约束按集合索引；
- 在 `milp` 和 PuLP 中建立并求解含 0–1 变量（binary variable）的 MILP，并用线性规划松弛的界比较不同的建模形式；
- 通过求解状态、约束残差以及小规模实例上的枚举（enumeration）来检查模型。
:::

所有代码块都使用以下导入语句。

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
## 产品组合 LP 的三种写法

### 用 SciPy 写成矩阵形式

`linprog` 求解 $\min\{\vc^\T\vx:\mA_{ub}\vx\le\vb_{ub},\ \mA_{eq}\vx=\vb_{eq},\ \vl\le\vx\le\vu\}$。产品组合 LP {eq}`eq-03-01-product-mix` 是含三个 $\le$ 行的最大化问题，因此以 $\min(-\vc)^\T\vx$ 的形式传入。一个简短的包装函数隐藏了目标函数和对偶值中的两处符号变换（见{numref}`tab-09-01-duals`），并在 LP 没有最优解时报错终止。

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

矩阵形式很紧凑，但每一列、每一行的含义只能写在注释里。对于两个变量的模型这样做没有问题，对于两千个变量的模型则无法管理。

### 用 PuLP 写成代数形式

代数模型把**数据**（集合与参数）与**模型**（变量、目标函数、约束）分开，正如 AMPL 的 `.mod` 文件和 `.dat` 文件那样（见{ref}`sec-09-01-modelling-languages`）。在 Python 中，数据是以产品名和资源名为键的字典，模型则是数据的函数。

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

这里需求上限写成变量的上**界**，而不是一个约束行，求解器本来也是这样处理简单界的；因此它没有行对偶值。它的松弛量为 $8-6.429=1.571$，由于它不是紧约束（binding constraint），其对偶值为零。

PuLP 还可以把模型以 LP 文件格式打印出来，这是检查数据是否按预期读入的最快方法：

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

### 用 CVXPY 写成向量形式

CVXPY（需要安装 CVXPY）使用与 SciPy 相同的数组，但保留最大化和 $\le$ 号的原样。由 `A @ x <= b` 返回的约束对象在求解后带有对偶值。

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

三种接口都给出 $\vx^*\approx(6.429,4.286)$（每周百箱）和 $z^*\approx\$5142.86$，即{prf:ref}`ex-03-03-simplex-run`中的最优解 $(\tfrac{45}{7},\tfrac{30}{7})$ 与最优值 $\tfrac{36000}{7}$。

(sec-09-02-duals)=
## 对偶值、影子价格与松弛量

### 求解器输出与理论的对照

对偶值 $\vpi^*=(78.571,\ 2.857,\ 0)$ 就是{prf:ref}`ex-03-05-product-mix-dual`中的影子价格 $\tfrac{550}{7}$ 美元/小时、$\tfrac{20}{7}$ 美元/立方米和 $0$，也是最优单纯形表（{numref}`tab-03-03-pm-3`）中松弛变量各列下方 $-z$ 行元素的相反数。互补松弛性（complementary slackness，见{prf:ref}`thm-03-05-complementary-slackness`）可以直接从输出中读出：两个紧约束行的价格为正，松弛量为 $1.571$ 的需求行价格为零，两种产品都在生产，检验数均为零。强对偶性（strong duality）表明，资源的价值恰好等于最优贡献，即 $\vb^\T\vpi^*=z^*$：

```python
print(round(b @ pi, 2), round(c @ x, 2))
print(np.allclose(pi, [550 / 7, 20 / 7, 0.0]))
# Expected:
# 5142.86 5142.86
# True
```

### 通过重新求解检验影子价格

影子价格预测的是某种资源每增加一个单位时最优值的变化，但只在其有效范围（validity range）内成立（见{prf:ref}`thm-03-05-shadow-price`）。对若干不同的机时数重新求解 LP，既能显示这一预测，也能显示它的局限。由{ref}`sec-03-05-shadow-prices`，当每周机时在 37.5 小时到 65.5 小时之间时，$\pi_1^*=\tfrac{550}{7}$ 成立。

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

在 37.5 小时到 65.5 小时之间，每半小时价值 $\tfrac12\cdot78.57=\$39.29$，增加一小时使贡献由 5142.86 美元增至 5221.43 美元，与{ref}`sec-03-05-shadow-prices`的预测一致。超过 65.5 小时后需求上限成为紧约束，更多的机时不再有价值；低于 37.5 小时时只生产鸡尾酒杯，每小时产出 $\tfrac{450}{5}=\$90$（每半小时 $\$45$），价格随之改变。最优值是资源量的凹分段线性函数，这正是{prf:ref}`thm-03-07-value-convex`对最小化形式所述的结论。

### 等式约束的对偶值与检验数

对于标准形（standard form）的最小化问题 $\min\{\vc^\T\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$，`linprog` 在 `res.eqlin.marginals` 中返回各等式的乘子（符号不限），在 `res.lower.marginals` 中返回处于下界零的变量的检验数 $\bar c_j=c_j-\vpi^\T\mA_{\cdot j}$。对于{doc}`03-07-sensitivity-analysis`中的贯穿示例（running example）{eq}`eq-03-07-running-lp`，它们应当是 $\vpi=(2,-4,1)$ 和 $\bar\vc=(0,0,0,1,3,8)$，即{numref}`tab-03-07-optimal-tableau`中的 $-z$ 行。

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

同一模型在 PuLP 中通过 `pi` 和 `dj` 返回相同的数值；对于最小化问题，它们同样是最优值的导数：

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

$\pi_2=-4$ 的符号是有意义的：增大 $b_2$ 会*降低*最小费用。只要基保持最优，同时改变多个数据的效果也可以用同样的方法预测。把 $b_1$ 增至 12、$b_3$ 增至 14，用掉了允许增量的 75%（见{prf:ref}`thm-03-07-hundred-percent`），因此预测的最优值为 $11+2\cdot1+1\cdot1=14$：

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

无论 SciPy 还是使用 CBC 的 PuLP，都不报告对偶值和最优基保持有效的**范围**（即{numref}`tab-03-07-ranges`中的费用与右端项变化范围）。{ref}`sec-09-03-sensitivity`用{doc}`03-07-sensitivity-analysis`中的公式，由最优基计算这些范围。

(sec-09-02-data-driven)=
## 数据驱动的模型

### 配料模型

{prf:ref}`ex-08-01-foundry`中的铸造厂须交付 1000 千克铸件，含锰量至少 0.45%，含硅量在 3.5% 到 6% 之间。它可以熔炼三种生铁，其成分和费用（购价加每 1000 千克 9 美元的熔炼费）按每千千克给出，还可以按每千克 14 美元的价格加入纯锰。以 $x_k$ 表示生铁 $k$ 的千千克数，$x_{Mn}$ 表示纯锰的千克数，该模型就是{ref}`sec-03-01-blending`中总量固定的配料模型。由数据写出时，每个约束只占一行，增加一种生铁只需在字典中多加一项。

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

最优解为 $\vx^*=(\tfrac{366}{431},\tfrac{15}{431},\tfrac{50}{431},0)$ 千千克，费用为 $\tfrac{20444}{431}\approx\$47.43$，与{prf:ref}`ex-08-01-foundry`一致。对偶值解释了这个解。通过生铁配比满足锰含量要求时，每多一千克锰的费用为 13.46 美元，低于纯锰的 14 美元，所以不购买纯锰。硅含量上限有 25 千克的松弛量，价格为零。对于 $\ge$ 行，PuLP 的松弛量（右端减左端）$\le0$，因此剩余量会以负号出现。

### 多周期模型

多周期模型最能体现索引的好处。{prf:ref}`ex-08-01-production-inventory`中的生产–库存问题有两种产品、三个月，以及每种产品每月一个的平衡方程（balance equation）{eq}`eq-03-01-inventory-balance`，期初无库存。变量是以（产品，月份）为键的字典。

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

生产计划和费用 $z^*=\$283\,900$ 与{prf:ref}`ex-08-01-production-inventory`相同。七月是瓶颈。其影子价格为每小时 $-\$1.20$，即七月多一小时所节省的费用：它可以使 $\tfrac43$ 单位的 A 改在七月而不是六月生产，节省 $\tfrac43\cdot\$0.90=\$1.20$ 的库存持有费用。这多出的一小时*不会*用于 B：把一单位 B（一小时）从六月移到七月只节省 0.75 美元的持有费用。另外两个月有松弛，价格为零。

(sec-09-02-transportation)=
## 运输模型

{prf:ref}`ex-04-04-transportation`中的运输问题（transportation problem）每周从工厂 A 和 B 分别运出 40 吨和 50 吨，送往需求量为 20 吨、25 吨和 45 吨的三个市场。其供应行与需求行就是平衡方程{eq}`eq-04-04-transportation`；弧就是费用字典的键。

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

该运输方案每周费用 445 美元，与{ref}`sec-04-04-ssp`中用连续最短路算法（successive shortest paths）求得的结果相同；而且它是整数解，尽管模型中并没有这样要求：约束矩阵是全单模的（见{prf:ref}`thm-04-04-integrality`）。

**对偶值即节点势。** 以 $u_i$ 表示工厂 $i$ 供应行的对偶值，$v_j$ 表示市场 $j$ 需求行的对偶值，则弧 $(i,j)$ 的检验数为 $c_{ij}-u_i-v_j$，在每条被使用的弧上为零。五个等式中有一个是冗余的（总供应等于总需求），因此对偶值只确定到相差一个常数：可以给所有 $u_i$ 加上一个常数，同时从所有 $v_j$ 中减去它；不同的求解器会报告不同但同样正确的值。{ref}`sec-04-04-ssp`中的节点势（node potential）对工厂为 $\pi_i=u_i+k$，对市场为 $\pi_j=k-v_j$，其中取 $k=v_3$ 使 $\pi_3=0$：

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

这正是{ref}`sec-04-04-ssp`中的节点势 $\pi_A=7$、$\pi_B=5$、$\pi_1=3$、$\pi_2=1$、$\pi_3=0$，两条未使用的弧具有正的检验数 $\bar c_{A3}=2$ 和 $\bar c_{B2}=3$，这证明该方案最优且唯一（见{prf:ref}`thm-04-04-reduced-cost`）。

对于 `linprog`，同一模型需要显式给出约束矩阵：每条弧一列，在其所属工厂的行和所属市场的行上各有一个 1。由弧表构建该矩阵就是{doc}`04-04-minimum-cost-flow`中关联矩阵的构造方法，结果相同：

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
## 混合整数模型

### `scipy.optimize.milp` 中的 0–1 变量

`milp` 用 HiGHS 的分支切割程序求解 $\min\{\vc^\T\vx:\vb_l\le\mA\vx\le\vb_u,\ \vl\le\vx\le\vu,\ x_j\in\mathbb Z\ \text{对}\ j\in J\}$。约束以 `LinearConstraint(A, lb, ub)` 的形式传入（等式取 `lb = ub`，$\le$ 行取 `lb = -inf`），界以 `Bounds(l, u)` 传入，整数变量由 `integrality` 向量标记，整数取 1、连续取 0；0–1 变量就是界为 $[0,1]$ 的整数变量。

{prf:ref}`ex-07-01-facility`中的设施选址问题有四个候选地点，固定费用为 $F_i$（千美元/年），容量为 $u_i$（千吨/年），还有五个废物量为 $d_j$ 的区；$c_{ij}$ 的单位为美元/吨，也就是千美元/千吨。变量是 20 个流量 $x_{ij}$，按行存放（$x_{11},\dots,x_{15},x_{21},\dots$），其后是四个 0–1 变量 $\delta_i$。下面的函数构建{prf:ref}`ex-07-01-facility`中比较的三种建模形式，并求解 MILP 或其线性规划松弛；它还可以固定 0–1 变量的取值，下文的枚举就要用到这一功能。

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

最优解开设地点 2 和 3，费用为每年 120 万美元：固定费用 53 万美元，运输费用 67 万美元。（换一个求解器版本，区 3 的流量在两个开设地点之间的分配可能不同，但费用相同。）三种建模形式的线性规划松弛给出{prf:ref}`ex-07-01-facility`中的界；联系约束越紧，界越接近整数最优值：

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

### 用 PuLP 写同一个 MILP

在 PuLP 中，0–1 变量写作 `cat="Binary"`，模型读起来与{prf:ref}`ex-07-01-facility`一样；求解器是 CBC 而不是 HiGHS。

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

### 用枚举法检验

只有四个地点时，开设地点的选择只有 $2^4=16$ 种，因此可以固定 0–1 变量，对每种选择求解剩下的运输 LP，以此检验 MILP，这样就重现了{prf:ref}`ex-07-01-facility`中的表格。

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

只开设单个地点以及地点对 $\{2,4\}$（容量 130 千吨，废物 140 千吨）都不可行，$\{2,3\}$ 是费用最低的选择。枚举的工作量按 $2^n$ 增长，对实际规模的实例毫无用处，但在小实例上，它是检验建模形式最有说服力的方法。

### 择一约束

{prf:ref}`ex-07-01-either-or`在方形区域 $[0,4]^2$ 上最大化 $-x_1+2x_2$，要求 $\vx$ 位于三角形 $x_1+x_2\le2$ 中，或位于正方形 $x_1\ge3$, $x_2\ge3$ 中。令 $\delta=1$ 表示选择三角形，则大 $M$ 约束（big-M constraint）为 $x_1+x_2-2\le M_T(1-\delta)$、$3-x_1\le M_1\delta$、$3-x_2\le M_2\delta$。整数最优解与这些常数无关（只要它们有效），但线性规划松弛与之有关：

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

两个 MILP 都返回 $\vx^*=(3,4)$、$\delta^*=0$，目标函数值为 5。由{eq}`eq-07-01-m-from-bounds`得到的紧常数 $M_T=6$, $M_1=M_2=3$ 给出界 6，而宽松的 $M=10$ 只给出界 8，如{numref}`fig-07-01-either-or`所示。

### 凹分段线性费用

{prf:ref}`ex-07-03-economies`中的两座处理厂须提供 100 百万升/日的处理能力，每座至多 80 百万升/日，费用为凹函数 $C_1(q)=q^{0.6}$ 和 $C_2(q)=0.8q^{0.6}+0.06q$（百万美元/年）。取分段点（breakpoint）$a_i=0,20,40,60,80$ 百万升/日，把每个费用函数替换为凸组合形式{eq}`eq-07-03-sos2`中的分段线性插值函数；由于费用是凹的，必须用{eq}`eq-07-03-lambda`中的区段 0–1 变量强制权重相邻：$\alpha_i\le\delta_{i-1}+\delta_i$ 且 $\sum_k\delta_k=1$。下面的函数构建带或不带这些 0–1 变量的模型。

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

两个模型都报告 $q_1=80$、$q_2=20$ 百万升/日，但只有 MILP 给出了正确的费用。不带相邻条件的 LP 把 $q_2=20$ 写成分段点 0 的 $0.75$ 倍加分段点 80 的 $0.25$ 倍，收取的是弦上的值 $0.25\cdot15.89$，而不是 $C_2(20)=6.03$；它给出的 1 784 万美元是凸包络（convex envelope）的费用，只是一个下界，而不是真实费用。加上 0–1 变量后，模型返回 1 989 万美元，即精确的最优值（见{prf:ref}`ex-07-03-economies`）。有些求解器可以直接以 SOS2 约束的形式接受相邻条件（Gurobi、CPLEX、Xpress，通过它们各自的接口）；PuLP 的 CBC 命令行接口不会把 SOS 信息传给 CBC，因此使用 PuLP 时，{eq}`eq-07-03-lambda`中的 0–1 变量才是稳妥的选择。

(sec-09-02-checking)=
## 检查模型

### 先看求解状态

{doc}`09-01-code-setup-solvers`列出了各种状态码。{doc}`03-07-sensitivity-analysis`中的两个模型展示了两种失败情形。当 $b_1=15$ 时，贯穿示例不可行（见{prf:ref}`ex-03-07-rhs`）；加入新增活动 $\mA_{\cdot7}=(3,-1,1)^\T$、$c_7=4$ 后，它无界（见{prf:ref}`ex-03-07-new-activity`）：

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

调试不可行模型的方法是一次放松一族约束（或者给每一行加上费用很大的非负松弛变量，看哪些松弛变量变为正值）；调试无界模型的方法是用一个很大的箱形区域约束每个变量，看哪些变量跑到箱子的边界上。

### 残差

求解器报告的最优解可以通过计算约束值来独立检验。下面的函数返回 $\mA_{ub}\vx\le\vb_{ub}$、$\mA_{eq}\vx=\vb_{eq}$ 和 $\vx\ge\vzero$ 的最大违反量；任何超过约 $10^{-7}$ 的可行性容差的值，都表明模型或数据处理中有错误。

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

最后一行计算产品组合的不可行方案 $(8,4)$，它需要 68 机时和 160 立方米：仓储行违反 10 立方米，机时行违反 8 小时，报告的是较大的违反量。

### 小规模实例的枚举

把 LP 最优解取整与求解整数规划并不是一回事（见{ref}`sec-03-01-when-lp-fails`）。对于以整百箱计的产品组合，$x_1\in\{0,\dots,8\}$，$x_2\in\{0,\dots,7\}$（由仓储行），因此可以尝试全部 72 个整数点，其中 54 个可行：

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

最好的整数方案 $(8,2)$ 每周赚 4 900 美元，比取整后的 LP 最优解 $(6,4)$ 多 100 美元，与{ref}`sec-03-01-when-lp-fails`中指出的一致；`milp` 的结果与枚举一致。一个在小实例上通过了这类检验、且对偶值在经济上讲得通的模型，就可以用于实际数据了。

## 本章回顾

### 小结

- `linprog` 接受矩阵形式；PuLP 由具名变量和数据字典构建同一个 LP；CVXPY 将其写成向量形式。三者都给出产品组合的最优解 $\vx^*\approx(6.429,4.286)$、$z^*\approx\$5142.86$（见{ref}`sec-09-02-product-mix`）。
- 按{numref}`tab-09-01-duals`的符号约定换算后，对偶值就是影子价格：产品组合为 $(78.57,2.86,0)$；{doc}`03-07-sensitivity-analysis`的贯穿示例为 $\vpi=(2,-4,1)$，检验数为 $(0,0,0,1,3,8)$；在各自的有效范围内重新求解可以证实它们（见{ref}`sec-09-02-duals`）。
- 带索引的模型（配料、生产–库存、运输）只需写一次，再由数据填充；它们的对偶值解释了所得的解（每千克锰 13.46 美元、七月每小时 1.20 美元、节点势 $(7,5,3,1,0)$）（见{ref}`sec-09-02-data-driven`、{ref}`sec-09-02-transportation`）。
- `milp` 和 PuLP 求解含 0–1 变量的 MILP；设施选址问题开设地点 2 和 3，费用为每年 120 万美元，其 LP 界 751.07、976.79 和 1168.46 衡量了三种建模形式的强弱；带 0–1 变量的 $\lambda$ 形式能正确计算凹费用（1 989 万美元，而不带相邻条件时为 1 784 万美元）（见{ref}`sec-09-02-milp`）。
- 每个模型都要通过求解状态、残差以及小实例上的枚举加以检查（见{ref}`sec-09-02-checking`）。

### 常见错误

- **忘记 `linprog` 的符号变换**：最大化问题和 $\ge$ 行都需要变号；应像 `solve_max_lp` 那样把它们封装在一个函数里。
- **把退化或非唯一最优解的对偶值当作唯一的。** 运输问题的对偶值只确定到相差一个常数；比较之前要先归一化。
- **在有效范围之外使用影子价格。** 应像{ref}`sec-09-02-duals`那样重新求解加以检验。
- **把索引写死在代码里。** 应由数据字典构建约束；增加一种产品、一个月份或一个地点，不应需要修改模型。
- **平衡方程中的差一错误。** 应显式写出期初库存（此处为零），并在数据已知的周期上测试模型。
- **把线性规划松弛的值当作费用。** 对于非凸的分段线性费用，松弛收取的是凸包络的费用；只有带相邻条件的 MILP 才给出真实费用。
- **在使用自带 CBC 的 PuLP 中依赖 `prob.sos2`**：这些集合不会传给求解器，得到的是不带相邻条件的 LP。

## 参考文献

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
