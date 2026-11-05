# 软件环境与求解器

## 概述

第九部分把第二至第七部分中的方法转化为 Python 代码。本章为此做好准备：搭建一个包含代码各章所需少量软件包的 Python 环境；概览完成数值计算的**求解器**（solver）——包括线性规划（linear program, LP）与混合整数线性规划（mixed-integer linear program, MILP）的开源求解器（open-source solver）和商业求解器（commercial solver），以及面向凸问题的锥规划求解器——和用来为它们编写模型的**建模语言**（modeling language）；并说明如何解读求解器的返回结果：求解状态、最优值、原始变量值、对偶值及其符号约定、松弛量，以及使所有这些量都只是近似值的容差。本章最后给出一个简短的安装检验：分别用 SciPy、PuLP 和 CVXPY 求解{prf:ref}`ex-03-01-product-mix`的产品组合线性规划，并把三个结果与精确最优解 $\vx^*=(\tfrac{45}{7},\tfrac{30}{7})$、每周 $z^*=\tfrac{36000}{7}\approx\$5142.86$ 进行比较。本章假定读者已经熟悉产品组合模型（{doc}`03-01-lp-model-building`）、标准形与松弛变量（{doc}`03-02-lp-geometry`）、影子价格（{doc}`03-05-lp-duality`），以及基本的 Python 和 NumPy 数组。

第九部分的代码**在构建本书时并不运行**。请把一章中的代码块按顺序复制到一个脚本或笔记本中（每章是一个命名空间，就像笔记本中的各个单元格），并把打印结果与 `# Expected:` 注释进行比较，这些注释重复了例题中的数值。最后一位打印数字的微小差别，或者出现 `-0.` 而不是 `0.`，都源于浮点运算，并无妨碍。

:::{admonition} 学习目标
:class: objectives
学完本章后，读者应能够：

- 创建虚拟环境（virtual environment），安装第九部分所用的软件包并检验它们能否正常工作；
- 说出主要的开源与商业 LP、MILP 及锥规划求解器，并指出每种求解器能处理哪些问题类型；
- 比较各种建模语言与建模框架（AMPL、GAMS、Pyomo、JuMP、OR-Tools、PuLP、CVXPY、YALMIP），并读懂一个小型 AMPL 模型；
- 读取求解器给出的求解状态、最优值、原始变量值、对偶值和松弛量，并在 SciPy、PuLP、CVXPY 的符号约定与{doc}`03-05-lp-duality`中的影子价格之间换算对偶值；
- 解释可行性容差、最优性容差和整数性容差，并据此把浮点结果与精确值进行比较。
:::

(sec-09-01-environment)=
## Python 环境准备

代码需要 Python 3.11 或更高版本（这是近期 SciPy 版本所支持的最低版本）。良好的做法是把一个项目的软件包安装到一个**虚拟环境**中。虚拟环境是一个文件夹，其中保存了解释器软件包目录的独立副本，这样这里使用的版本就不会干扰其他项目。在本书代码仓库的根目录下执行：

```bash
python -m venv .venv                 # create the environment once
source .venv/bin/activate            # Linux and macOS
.venv\Scripts\activate               # Windows (PowerShell or cmd)
pip install -r requirements.txt      # numpy, scipy, networkx, pulp, cvxpy, ...
```

文件 `requirements.txt` 列出了各软件包及其最低版本；同一文件还会安装构建本书所需的工具。如果没有代码仓库，代码各章只需要

```bash
pip install numpy scipy networkx matplotlib pulp cvxpy
```

Jupyter（`pip install jupyterlab`）便于交互式地运行代码块，但并非必需。后面的每个代码块都假定已经执行了所在章第一个代码块中的导入语句。本章的导入语句为：

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

第二个代码块检查解释器和已安装的版本是否满足 `requirements.txt` 中的最低版本。具体的版本号因机器而异；重要的是每项检查都打印 `True`。

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
## 第九部分使用的软件包

除标准库外，第九部分使用六个软件包。{numref}`tab-09-01-packages`列出了它们及其用途。标准库中的 `fractions.Fraction` 值得单独一提：{doc}`09-03-code-lp-algorithms`中可以手算规模的单纯形表都用精确有理数计算，因而能逐位重现第三部分中印出的单纯形表。

```{list-table} 第九部分使用的软件包
:name: tab-09-01-packages
:header-rows: 1

* - 软件包
  - 用途
  - 章节
* - NumPy
  - 数组、矩阵–向量乘积、线性方程组（`numpy.linalg`）
  - 全部
* - SciPy
  - `scipy.optimize.linprog` 和 `milp`（两者都调用 HiGHS），用于 LP 和 MILP；`minimize` 用于非线性规划（nonlinear program, NLP）
  - {doc}`09-02-code-lp-modeling`、{doc}`09-03-code-lp-algorithms`、{doc}`09-04-code-networks`、{doc}`09-06-code-nlp`
* - PuLP
  - 用具名变量和具名约束对 LP 和 MILP 进行代数建模；自带 CBC 求解器
  - {doc}`09-02-code-lp-modeling`
* - CVXPY
  - 以向量形式对凸问题（LP、二次规划（quadratic program, QP）和二阶锥规划（second-order cone program, SOCP））建模，并自动检查凸性；调用锥规划求解器
  - {doc}`09-02-code-lp-modeling`、{doc}`09-06-code-nlp`
* - NetworkX
  - 图以及网络算法的参考实现
  - {doc}`09-04-code-networks`
* - Matplotlib
  - 绘图，仅在图形必不可少时使用
  - {doc}`09-06-code-nlp`
```

在它们之间作选择时，有两点值得注意。**SciPy** 要求模型以*矩阵形式*给出，即 $\min\{\vc^\T\vx:\mA_{ub}\vx\le\vb_{ub},\ \mA_{eq}\vx=\vb_{eq},\ \vl\le\vx\le\vu\}$：紧凑而快速，但建模者必须自己为变量编号并构造矩阵。**PuLP** 和 **CVXPY** 则允许建模者用具名变量的表达式写出约束，由它们自己构造矩阵。PuLP 只限于线性模型和混合整数线性模型；CVXPY 能处理任何遵循其规则（规范凸规划，disciplined convex programming）的凸模型，但只有在所用求解器支持时才能处理整数变量。

(sec-09-01-solvers)=
## 求解器概览

**求解器**是对以数据形式给出的模型（矩阵、界、变量类型）运行某种算法的程序。**建模层**（modeling layer）（PuLP、CVXPY、AMPL 等）把建模者所写的模型翻译成这些数据，再把求解器的答案翻译回具名的量。大多数建模层都能调用多个求解器，因此这两项选择在很大程度上是相互独立的。

### LP 与 MILP 求解器内部的算法

LP 求解器实现了三类方法。{doc}`03-03-simplex-method`中的**原始单纯形法**（primal simplex method）和{doc}`03-06-revised-dual-simplex`中的**对偶单纯形法**（dual simplex method）在基之间移动，采用修正形式，用 $\mB$ 的 LU 分解代替单纯形表；对偶单纯形法通常是默认方法，因为在问题改变后它能迅速重新优化。**内点法**（interior-point method，即障碍函数法）沿可行域内部的一条路径前进，所需迭代次数少得多，但每次迭代的代价更高；需要最优基时（例如为了做灵敏度分析），最后用一个**交叉（crossover）**步骤恢复出最优基。MILP 求解器对线性规划松弛（LP relaxation）运行**分支定界法**（branch and bound），并辅以**割平面**（cutting plane）、预处理化简和原始启发式算法加以强化。所有这些求解器都在带容差的浮点运算中工作（{ref}`sec-09-01-tolerances`）。

### 开源 LP 与 MILP 求解器

- **HiGHS**（Huangfu & Hall, 2018）是爱丁堡大学开发的面向 LP、MILP 和凸 QP 的高性能求解器，以 MIT 许可证发布。它包含对偶与原始修正单纯形法程序、内点法程序以及分支切割 MILP 求解器。SciPy 的 `linprog` 默认使用 HiGHS（`method="highs"`），`scipy.optimize.milp` 调用 HiGHS 的 MILP 求解器；HiGHS 也有自己的 Python 接口 `highspy`。
- **CBC**（COIN-OR Branch and Cut）是 COIN-OR 项目的 MILP 求解器；它用 COIN-OR 的单纯形法程序 CLP 求解 LP。PuLP 自带一个 CBC 可执行文件，因此 PuLP 无需任何额外安装即可使用。
- **GLPK**（GNU Linear Programming Kit）是一个较老、较慢但可移植性极好的 LP/MILP 软件包，有自己的建模语言（GNU MathProg，AMPL 的一个子集）；当 `glpsol` 在搜索路径中时，PuLP 通过 `GLPK_CMD` 调用它。
- **SCIP** 是一个约束整数规划框架，可求解 MILP 和混合整数非线性规划（mixed-integer nonlinear program, MINLP）；它以 Apache 2.0 许可证发布（自 8.0.3 版起），并有 Python 接口 PySCIPOpt。

### CVXPY 背后的锥规划求解器

CVXPY 把凸模型改写为**锥规划**（conic program）：线性目标函数、线性方程，以及要求某些仿射表达式属于某些锥（LP 对应非负卦限，SOCP 对应二阶锥，半定规划（semidefinite program, SDP）对应半正定矩阵锥；见{doc}`02-03-convex-optimization-duality`）。然后它调用下列求解器之一：

- **Clarabel**，用 Rust 编写的内点法求解器，适用于 LP、QP、SOCP、SDP 以及指数锥和幂锥；近期的 CVXPY 版本默认安装它，并对大多数问题优先使用它；
- **SCS**（O'Donoghue et al., 2016），面向大规模锥问题的一阶分裂方法（交替方向乘子法，alternating direction method of multipliers, ADMM），精度约为 $10^{-4}$ 到 $10^{-6}$，而达不到机器精度；
- **OSQP**（Stellato et al., 2020），面向凸 QP 的 ADMM 求解器，速度快且支持热启动，在控制领域很流行；
- **ECOS**（Domahidi et al., 2013），面向嵌入式系统的内点法 SOCP 求解器，曾是 CVXPY 的默认求解器，现在需要另行安装。

CVXPY 还可以调用 HiGHS、GLPK、CBC 和 SCIP 求解 LP 和 MILP，也可以调用下面介绍的商业求解器。`cp.installed_solvers()` 会列出当前可用的求解器。

### 商业求解器

商业求解器在大规模、高难度的 MILP 上速度更快，预处理更好，数值上更稳健，并能直接报告灵敏度范围。下面四种求解器都向授予学位的院校的学生和教职工提供**免费学术许可**（academic license）；用大学电子邮箱注册即可获得许可证密钥，并在用户自己的机器上激活（Gurobi 需在校园网内执行命令 `grbgetkey <key>`，CPLEX 则通过 IBM 的学术计划激活）。Gurobi、CPLEX 和 Xpress 的 Python 软件包可用 `pip` 安装，并附带一个只能求解小规模模型的受限许可证；MOSEK 的软件包也可用 `pip` 安装，但需要试用许可证或学术许可证。

- **Gurobi**（LP、QP、SOCP、MILP、混合整数 QP 以及非凸二次问题）；Python API 为 `gurobipy`。
- **IBM ILOG CPLEX**（LP、QP、SOCP、MILP），是 CPLEX Optimization Studio 的一部分，该套件还包含 OPL 建模语言；Python API 为 `cplex`，建模层为 `docplex`。
- **FICO Xpress**（LP、QP、MILP、非线性问题）；Python API 为 `xpress`。
- **MOSEK**（LP、QP、包括 SDP 在内的锥问题、混合整数锥问题）；尤其擅长锥问题；Python API 为 `mosek`。

更换求解器不会改变模型。在 PuLP 中，求解器是 `solve` 的一个参数，例如 `prob.solve(pulp.GUROBI_CMD())`，或 `pulp.CPLEX_CMD()`、`pulp.XPRESS()`、`pulp.MOSEK()`、`pulp.GLPK_CMD()`；`pulp.listSolvers(onlyAvailable=True)` 显示 PuLP 能找到哪些求解器。在 CVXPY 中，求解器是一个关键字参数，`problem.solve(solver=cp.GUROBI)`（或 `cp.CPLEX`、`cp.XPRESS`、`cp.MOSEK`、`cp.HIGHS`、`cp.SCS` 等）。{numref}`tab-09-01-solvers`对上述内容作了总结。

```{list-table} 求解器及其能处理的问题类型（✓ = 支持）
:name: tab-09-01-solvers
:header-rows: 1

* - 求解器
  - 许可证
  - LP
  - MILP
  - QP / SOCP
  - SDP
  - 调用途径
* - HiGHS
  - 开源（MIT）
  - ✓
  - ✓
  - 凸 QP
  -
  - SciPy, PuLP, CVXPY, Pyomo, JuMP
* - CBC / CLP
  - 开源（EPL）
  - ✓
  - ✓
  -
  -
  - PuLP（自带）, CVXPY, Pyomo, JuMP
* - GLPK
  - 开源（GPL）
  - ✓
  - ✓
  -
  -
  - PuLP, CVXPY, Pyomo
* - SCIP
  - 开源（Apache）
  - ✓
  - ✓
  - ✓（以及 MINLP）
  -
  - PySCIPOpt, CVXPY, Pyomo
* - Clarabel, ECOS
  - 开源
  - ✓
  -
  - ✓
  - 仅 Clarabel
  - CVXPY
* - SCS
  - 开源（MIT）
  - ✓
  -
  - ✓
  - ✓
  - CVXPY
* - OSQP
  - 开源
  -
  -
  - 仅 QP
  -
  - CVXPY
* - Gurobi
  - 商业，学术免费
  - ✓
  - ✓
  - ✓
  -
  - gurobipy, PuLP, CVXPY, Pyomo, AMPL, GAMS
* - CPLEX
  - 商业，学术免费
  - ✓
  - ✓
  - ✓
  -
  - cplex/docplex, PuLP, CVXPY, Pyomo, AMPL, GAMS
* - Xpress
  - 商业，学术免费
  - ✓
  - ✓
  - ✓
  -
  - xpress, PuLP, CVXPY, Pyomo, AMPL, GAMS
* - MOSEK
  - 商业，学术免费
  - ✓
  - ✓
  - ✓
  - ✓
  - mosek, PuLP, CVXPY, Pyomo, AMPL, GAMS
```

对本书中的模型而言，表中任何一个求解器都绰绰有余；只有在含数千个 0–1 变量的 MILP 和含数百万个变量的 LP 上，它们之间的差别才显得重要。

(sec-09-01-modelling-languages)=
## 建模语言与建模框架

**代数建模语言**（algebraic modeling language）允许建模者几乎像在纸上那样写出“对所有 $i\in I$，$\sum_{j\in J}a_{ij}x_j\le b_i$”：下标集、参数、变量、目标函数和约束分别声明，数据与模型分开存放。最著名的是 **AMPL**（Fourer et al., 2003）。一个 AMPL 项目通常有三个文件：模型（`.mod`）、数据（`.dat`）和命令脚本（`.run`）。作为示例，{prf:ref}`ex-03-01-product-mix`的产品组合模型写成

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

这些声明遵循若干值得在任何语言中借鉴的约定：下标集用大写，参数用小写，变量名和约束名首字母大写；参数的允许范围（`param avail {RESOURCES} >= 0;`）是对数据的检查，而不是模型的约束；每个约束和目标函数都有名字，通过名字可以显示它们的对偶值（`Capacity.dual`）和松弛量。更换数据文件，就能求解同一模型的另一个实例。

其他主要选择有：

- **GAMS**，与 AMPL 一样历史悠久的代数建模系统，在经济学和能源建模中应用广泛，有自己的语言和大量求解器接口。
- **Pyomo**（Bynum et al., 2021），开源的 Python 代数建模库，采用 AMPL 风格的集合、参数和带下标的变量；支持 LP、MILP 和非线性模型以及几乎所有求解器。
- **JuMP**（Lubin et al., 2023），Julia 语言的代数建模包，构建大规模模型速度快，且与求解器无关。
- **Google OR-Tools**，带有 Python 绑定的 C++ 库；它包含 LP 求解器 GLOP、约束规划与 SAT 求解器 CP-SAT（在排程问题和纯整数问题上非常强）、车辆路径程序，以及可调用 HiGHS、SCIP、Gurobi 等求解器的封装接口。
- **PuLP**（Mitchell et al., 2011），面向 LP 和 MILP 模型的轻量级 Python 库，用于{doc}`09-02-code-lp-modeling`；以及面向凸模型的 **CVXPY**（Diamond & Boyd, 2016）。
- **MATLAB** 及其 Optimization Toolbox（`linprog`、`intlinprog`、`quadprog`、`fmincon`），以及 **YALMIP**（Löfberg, 2004）——一个免费的 MATLAB 工具箱，其作用相当于 CVXPY 或 PuLP：符号决策变量（`sdpvar`、`binvar`）、以不等式写出的约束，以及可选择求解器的 `optimize(constraints, objective)`。

```{list-table} 建模语言与建模框架
:name: tab-09-01-modelling
:header-rows: 1

* - 框架
  - 语言
  - 许可证
  - 问题类型
  - 典型用途
* - AMPL
  - 自有语言（`.mod`/`.dat`/`.run`），Python API `amplpy`
  - 商业（有免费的社区版和学术版）
  - LP、MILP、NLP、MINLP
  - 教学、工业模型、非线性模型
* - GAMS
  - 自有语言
  - 商业
  - LP、MILP、NLP、MINLP、均衡模型
  - 经济学、能源系统
* - Pyomo
  - Python
  - 开源
  - LP、MILP、NLP、MINLP
  - Python 中的研究模型
* - JuMP
  - Julia
  - 开源
  - LP、MILP、锥规划、NLP
  - 大规模模型、算法研究
* - OR-Tools
  - C++, Python, Java, C#
  - 开源
  - LP、MILP、约束规划、路径问题
  - 排程与路径规划应用
* - PuLP
  - Python
  - 开源
  - LP、MILP
  - 中小规模模型、教学
* - CVXPY
  - Python
  - 开源
  - 凸问题（LP、QP、SOCP、SDP），配合合适的求解器可处理混合整数凸问题
  - 凸优化、统计、控制
* - YALMIP
  - MATLAB
  - 免费
  - LP、MILP、锥规划、NLP
  - MATLAB 中的控制工程
```

代码各章使用 SciPy、PuLP 和 CVXPY，因为它们免费、可用 `pip` 安装，并且足以处理本书的所有模型。其中的概念是通用的：在 PuLP 中用具名集合和参数构建的模型，可以逐行翻译成 AMPL 或 Pyomo。

(sec-09-01-reading-output)=
## 解读求解器输出

### 求解状态

求解器的第一个回答是**求解状态**（status）：它是找到了最优解，证明了问题不可行或无界，还是因其他原因（迭代次数或时间限制、数值困难）而停止？在检查求解状态之前，结果中的其他内容都没有意义。{numref}`tab-09-01-status`比较了三种 Python 接口。

```{list-table} 三种接口报告的求解状态
:name: tab-09-01-status
:header-rows: 1

* - 结果
  - `linprog` / `milp`（`res.status`）
  - PuLP（`pulp.LpStatus[prob.status]`）
  - CVXPY（`problem.status`）
* - 最优
  - `0`
  - `"Optimal"`
  - `"optimal"`
* - 达到迭代次数或时间限制
  - `1`
  - `"Not Solved"`
  - `"user_limit"` 或报错
* - 不可行
  - `2`
  - `"Infeasible"`
  - `"infeasible"`
* - 无界
  - `3`
  - `"Unbounded"`
  - `"unbounded"`
* - 数值困难
  - `4`
  - `"Undefined"`
  - `"optimal_inaccurate"` 等
```

不可行和无界是{doc}`03-02-lp-geometry`中没有最优解的两种结局；对于正在构建的模型，它们几乎总是意味着建模错误：不可行的模型通常有符号或单位错误，无界的模型通常缺少某个约束（{ref}`sec-03-01-modelling-steps`）。有些求解器在预处理之后无法区分这两种情形，只报告“不可行或无界”；把目标函数设为零再求解一次，即可判定是哪一种。

### 最优值与原始变量值

当求解状态为最优时，求解器报告最优目标函数值和一个原始解。有两个细节需要注意。第一，`linprog` 和 `milp` 总是*最小化*：最大化问题 $\max\vc^\T\vx$ 要以 $\min(-\vc)^\T\vx$ 的形式传入，最优值为 $-$`res.fun`（{eq}`eq-03-01-max-min`）。第二，LP 求解器在使用单纯形法或交叉步骤时返回一个最优基本解；但当最优解不唯一时（{ref}`sec-03-03-special-cases`），返回哪一个最优顶点取决于求解器、其版本及其设置。两个都正确的求解器可能报告不同的最优解，而最优值相同。

### 对偶值及其符号约定

对于最优基，求解器还知道{doc}`03-03-simplex-method`中的单纯形乘子 $\vpi^\T=\vc_B^\T\mB^{-1}$，它们是最优对偶值（{doc}`03-05-lp-duality`），在非退化最优解处也是影子价格 $\partial z^*/\partial b_i$（{prf:ref}`thm-03-05-shadow-price`）。求解器把它们称为 duals、dual values、marginals 或 shadow prices，而各接口对其符号的规定并不一致，因为符号取决于问题是最小化还是最大化、该行是 $\le$ 还是 $\ge$，以及接口报告的是拉格朗日函数的乘子还是最优值的导数。稳妥的做法是弄清每个接口报告的是什么，并把它换算成本书的影子价格，即*按原样表述*的问题的最优值在 $b_i$ 每增加一个单位时的变化率：

```{list-table} 三种接口中的对偶值与松弛量
:name: tab-09-01-duals
:header-rows: 1

* - 量
  - `linprog`
  - PuLP（CBC）
  - CVXPY
* - 一行的对偶值
  - `res.ineqlin.marginals`、`res.eqlin.marginals`：*最小化*目标函数 `res.fun` 对右端项的导数
  - `prob.constraints[name].pi`：按原样表述（max 或 min）的问题的目标函数的导数
  - `constraint.dual_value`：拉格朗日乘子，对最小化和最大化问题中的不等式约束都非负
* - $\max\vc^\T\vx$ 中 $\le$ 行的影子价格
  - `-res.ineqlin.marginals`（$\ge0$）
  - `pi`（$\ge0$）
  - `dual_value`（$\ge0$）
* - $\min\vc^\T\vx$ 中 $\le$ 行的影子价格
  - `res.ineqlin.marginals`（$\le0$）
  - `pi`（$\le0$）
  - `-dual_value`（$\le0$）
* - 检验数 $\bar c_j$
  - `res.lower.marginals`（对处于下界的变量）
  - `variable.dj`
  - 约束 `x >= 0` 的对偶值
* - 一行的松弛量
  - `res.slack` $=\vb_{ub}-\mA_{ub}\vx$；等式行用 `res.con`
  - `constraint.slack` $=$ 右端项 $-$ 左端项
  - `constraint.expr.value` 或手工计算的残差
```

`linprog` 只接受 $\le$ 行。$\ge$ 行 $\va^\T\vx\ge\beta$ 要以 $-\va^\T\vx\le-\beta$ 的形式传入，这会使其边际值（marginal）的符号再翻转一次：关于 $\beta$ 的导数等于所报告边际值的 $-1$ 倍。无论使用哪种接口，都可以按对偶值的含义来检验它：把 $b_i$ 增加一个小量 $\Delta$ 后重新求解，并把最优值的变化与 $\pi_i\Delta$ 比较，正如{ref}`sec-09-02-duals`中所做的那样。影子价格保持有效的精确范围（{doc}`03-07-sensitivity-analysis`）由商业求解器（以及配合商业求解器使用的 AMPL）报告；对于这里使用的开源接口，这些范围要从最优基计算出来，如{ref}`sec-09-03-sensitivity`所示。

### 松弛量

$\le$ 行的**松弛量**为 $b_i-\va_i^\T\vx^*\ge0$，即资源的未用量；对 $\ge$ 行，**剩余量** $\va_i^\T\vx^*-b_i\ge0$ 是需求被超出的量（{doc}`03-02-lp-geometry`）。松弛量为零的行是**紧约束**（binding constraint）。由互补松弛性（{doc}`03-05-lp-duality`），松弛量为正的行对偶值为零，这是检验任何求解器输出是否自洽的快捷方法。PuLP 把每一行的松弛量报告为右端项减左端项，因此 $\ge$ 行的剩余量会带负号出现。

(sec-09-01-tolerances)=
## 数值容差

求解器在双精度浮点运算中工作，约有 16 位有效十进制数字，其内部的每项检验都使用一个**容差**（tolerance）：

- **原始可行性容差**（primal feasibility tolerance）（通常为 $10^{-7}$）是仍被视为可行的约束或界的最大违反量；
- **对偶可行性容差**或**最优性容差**（optimality tolerance）（通常为 $10^{-7}$）是仍被视为最优的、符号错误的检验数的最大值（{prf:ref}`thm-03-03-optimality`）；
- **整数性容差**（integrality tolerance）（通常为 $10^{-6}$ 到 $10^{-5}$）是整数变量与整数之间允许的最大距离，因此 0–1 变量可能被报告为 $0.9999999$ 或 $10^{-9}$；
- **相对 MIP 间隙**（relative MIP gap）（通常为 $10^{-4}$，即 $0.01\%$）规定了在 MILP 求解器停止之前，最好的整数解与最好的界必须接近到什么程度；“最优”的 MILP 结果只是在这个间隙之内最优。

由此对本书的代码有两点推论。结果与期望值的比较使用带容差的 `math.isclose` 或 `numpy.allclose`，绝不用 `==`，并以固定的小数位数打印。此外，尺度不良的模型——系数跨越许多数量级（例如 $10^9$ 的大 $M$ 与量级为 1 的系数并存，{doc}`07-01-logical-constraints`）——会使容差以意想不到的方式起作用：一个被违反了 $10^{-7}\cdot10^9=100$ 个单位的约束仍可能被视为满足。选择适当的单位，使系数大致保持在 $10^{-3}$ 到 $10^{4}$ 之间，并在界允许的范围内把大 $M$ 常数取得尽可能小，就能避免大多数数值问题。SciPy 通过 `linprog` 和 `milp` 的 `options` 字典把容差传给 HiGHS，例如对 `linprog` 使用键 `primal_feasibility_tolerance`（取值 `1e-9`），对 `milp` 使用键 `mip_rel_gap`（取值 `1e-6`）。

精确的有理数运算完全避免了容差，代价是速度。它是检验小型例子的合适工具，{doc}`09-03-code-lp-algorithms`中的每个单纯形表都用它计算。

(sec-09-01-check)=
## 检验安装：用三种方式求解产品组合问题

产品组合线性规划 {eq}`eq-03-01-product-mix` 小到可以熟记于心：极大化 $500x_1+450x_2$，约束为机器工时 $6x_1+5x_2\le60$ 小时、仓储 $10x_1+20x_2\le150$ 立方米、需求 $x_1\le8$ 百箱，以及 $\vx\ge\vzero$。其最优解为每周 $x_1^*=\tfrac{45}{7}\approx6.429$、$x_2^*=\tfrac{30}{7}\approx4.286$ 百箱，每周的最优值为 $z^*=\tfrac{36000}{7}\approx\$5142.86$；影子价格为每机器工时 $\pi_1^*=\tfrac{550}{7}\approx\$78.57$、每立方米 $\pi_2^*=\tfrac{20}{7}\approx\$2.86$，以及 $\pi_3^*=0$；只有需求行有松弛量，为 $\tfrac{11}{7}\approx1.571$ 百箱（{prf:ref}`ex-03-03-simplex-run`、{prf:ref}`ex-03-05-product-mix-dual`）。数据只存储一次：

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

**SciPy。** `linprog` 做最小化，因此各项贡献带负号输入，影子价格是 `marginals` 的相反数（{numref}`tab-09-01-duals`）。

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

**PuLP。** 模型用具名变量和具名约束写出，并由自带的 CBC 求解；`msg=False` 关闭求解器日志。

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

CBC 输出的解约有八位有效数字，这就是可能出现 `-0.000` 的原因，也是打印前先对结果取整的原因。这些值与精确值在{ref}`sec-09-01-tolerances`的容差范围内一致：

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

**CVXPY**（需要安装 CVXPY）。模型用同样的数组以向量形式写出；保留约束对象，以便在求解后读取其对偶值。CVXPY 给出的不等式约束对偶值是非负乘子，对最大化问题而言它们就直接是影子价格。

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

Clarabel 这类内点法求解器停在最优顶点附近而非恰好在最优顶点上，因此它们的结果与单纯形法结果约有六位数字相同，非紧约束的对偶值可能打印为 `0.` 或一个极小的数。

**两个有意制造的失败。** 必须先读取求解状态，再看其他内容。要求果汁杯至少生产 9 百箱，这与需求限制 $x_1\le8$ 矛盾，会使模型不可行；漏掉机器工时行和仓储行则会使模型无界，因为此时没有任何约束限制鸡尾酒杯的产量。

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

在这两种情形下都不能使用 `res.x` 和 `res.fun`。最后，`pulp.listSolvers(onlyAvailable=True)` 和（安装了 CVXPY 时）`cp.installed_solvers()` 显示每个接口能调用哪些求解器。

```python
print(pulp.listSolvers(onlyAvailable=True))
# Expected: ['PULP_CBC_CMD']   (plus any other solver installed)
```

```python
print(cp.installed_solvers())
# Expected (requires CVXPY; the list depends on the installation):
# ['CLARABEL', 'ECOS', 'HIGHS', 'OSQP', 'SCIPY', 'SCS', ...]
```

如果本节所有代码块都打印出预期的值，环境就已为第九部分的其余内容准备就绪。

## 本章回顾

### 小结

- 用 `pip install -r requirements.txt` 建立的虚拟环境提供了第九部分所用的软件包：NumPy、SciPy、NetworkX、Matplotlib、PuLP（含 CBC）和 CVXPY（{numref}`tab-09-01-packages`）。
- 开源 LP/MILP 求解器有 HiGHS（SciPy 的 `linprog` 和 `milp` 背后的求解器）、CBC（随 PuLP 提供）、GLPK 和 SCIP；CVXPY 调用 Clarabel、SCS、OSQP 和 ECOS 等锥规划求解器；Gurobi、CPLEX、Xpress 和 MOSEK 是提供免费学术许可的商业求解器，可作为求解器参数接入 PuLP 和 CVXPY（{numref}`tab-09-01-solvers`）。
- 建模语言把模型与数据分开；AMPL 使用 `.mod`、`.dat` 和 `.run` 文件，Pyomo、JuMP、OR-Tools、PuLP、CVXPY 和 YALMIP 在 Python、Julia 和 MATLAB 中扮演同样的角色（{numref}`tab-09-01-modelling`）。
- 先读取求解状态（{numref}`tab-09-01-status`）；再读取最优值、原始变量值、对偶值和松弛量，并用{numref}`tab-09-01-duals`把对偶值换算为影子价格。
- 任何求解器结果都只在可行性、最优性、整数性和间隙容差的意义下是精确的（{ref}`sec-09-01-tolerances`）。
- SciPy、PuLP 和 CVXPY 都重现了产品组合问题的最优解 $\vx^*=(6.429,4.286)$、$z^*=\$5142.86$ 以及影子价格 $(78.57,\,2.86,\,0)$（{ref}`sec-09-01-check`）。

### 常见错误

- **不检查求解状态就使用 `res.x`。** 对于不可行或无界的模型，这些数组中是毫无意义的数，或者什么也没有。
- **忘记 `linprog` 做的是最小化。** 最大化问题必须以 $\min(-\vc)^\T\vx$ 的形式传入；此时最优值为 $-$`res.fun`，影子价格是 `marginals` 的相反数。
- **混淆对偶值的符号约定。** 对于某些优化方向与行类型的组合，SciPy 的 `marginals`、PuLP 的 `pi` 和 CVXPY 的 `dual_value` 符号不同（{numref}`tab-09-01-duals`）。可以把右端项稍微增大后重新求解，以此检验对偶值。
- **用 `==` 比较浮点数。** 应使用容差；报告为 $0.9999999$ 的 0–1 变量就是 $1$。
- **期望两个求解器给出同一个最优顶点。** 存在多重最优解时，报告的解取决于求解器；只有最优值是唯一的。
- **尺度不良的模型和过大的大 $M$ 常数**，它们会使容差变成真正的约束违反。

## 参考文献

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
