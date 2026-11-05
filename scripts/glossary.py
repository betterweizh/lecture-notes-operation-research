#!/usr/bin/env python3
"""English–Chinese terminology of the book; writes both glossary pages.

This module is the single source of the technical terms of *Lecture Notes on
Operations Research* and of their Simplified Chinese equivalents. The terms
were collected from the English chapters: bold terms at their first
definition, the titles of definitions, theorems and algorithms, section
headings, and the symbols of ``docs/en/notation.md``. The Chinese terms follow
the standard Chinese textbooks of operations research and optimization
(《运筹学》清华大学出版社; 袁亚湘《最优化理论与方法》; 刁在筠《运筹学》) and
the national terminology standards (全国科学技术名词审定委员会: 数学名词,
运筹学名词) where they exist. Translators of the Chinese edition use exactly
these terms.

Data layout
-----------
``SECTIONS`` is a list of sections ``(English title, 中文标题, terms)``; each
term is a tuple ``(English, 中文, symbol)``. The symbol is LaTeX written with
the book's MathJax macros (``\\vx``, ``\\bar c_j``, ...) and without dollar
signs, or ``""`` when the term has no symbol. Every English term appears only
once in the whole list.

Order
-----
Within each section the terms are sorted **alphabetically by the English
term** (case-insensitive, ignoring accents, dollar signs, backslashes and a
leading "the"). The data may therefore be written in any order. The Chinese
page uses the same order, so that both pages list the same term in the same
place and the English column of the Chinese page serves as the index.

Usage::

    python3 scripts/glossary.py     # writes docs/en/glossary.md and docs/zh/glossary.md

Both pages are unnumbered appendix pages: they hold tables only and no
equation labels.
"""

from __future__ import annotations

import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"

Term = tuple[str, str, str]  # (English, 中文, LaTeX symbol or "")
Section = tuple[str, str, list[Term]]  # (English title, 中文标题, terms)


# ---------------------------------------------------------------------------
# General words: page furniture, kinds of statements, chapter structure
# ---------------------------------------------------------------------------
GENERAL: list[Term] = [
    ("algorithm", "算法", ""),
    ("appendix", "附录", ""),
    ("assumption", "假设", ""),
    ("chapter", "章", ""),
    ("chapter review", "本章回顾", ""),
    ("common mistakes", "常见错误", ""),
    ("corollary", "推论", ""),
    ("definition", "定义", ""),
    ("enumeration", "枚举", ""),
    ("example", "例", ""),
    ("examples and problems", "例题与习题", ""),
    ("existence and uniqueness", "存在性与唯一性", ""),
    ("glossary", "术语表", ""),
    ("implementation", "实现", ""),
    ("input", "输入", ""),
    ("iteration", "迭代", ""),
    ("learning objectives", "学习目标", ""),
    ("lemma", "引理", ""),
    ("necessary condition", "必要条件", ""),
    ("notation", "符号说明", ""),
    ("operation count", "运算量", ""),
    ("output", "输出", ""),
    ("overview", "概述", ""),
    ("part (of the book)", "部分", ""),
    ("problem (exercise)", "习题", ""),
    ("proof", "证明", ""),
    ("proposition", "命题", ""),
    ("references", "参考文献", ""),
    ("remark", "注", ""),
    ("running example", "贯穿示例", ""),
    ("section", "节", ""),
    ("solution (of an example or problem)", "解答", ""),
    ("sufficient condition", "充分条件", ""),
    ("summary", "小结", ""),
    ("summary and common mistakes", "本章回顾", ""),
    ("theorem", "定理", ""),
    ("worked example", "例题", ""),
]

# ---------------------------------------------------------------------------
# Linear algebra (chapter 1.1)
# ---------------------------------------------------------------------------
LINEAR_ALGEBRA: list[Term] = [
    ("all-ones vector", "全一向量", r"\vone"),
    ("angle", "夹角", ""),
    ("augmented matrix", "增广矩阵", ""),
    ("back substitution", "回代", ""),
    ("basis (of a subspace)", "基", ""),
    ("block inverse", "分块矩阵的逆", ""),
    ("Cauchy–Schwarz inequality", "柯西–施瓦茨不等式", ""),
    ("characteristic polynomial", "特征多项式", ""),
    ("Cholesky factorization", "Cholesky 分解", ""),
    ("cofactor", "代数余子式", ""),
    ("cofactor expansion", "代数余子式展开", ""),
    ("column space", "列空间", r"\col(\mA)"),
    ("component (of a vector)", "分量", r"x_j"),
    ("condition number", "条件数", r"\kappa(\mA)"),
    ("definiteness", "有定性", ""),
    ("determinant", "行列式", r"\det\mA"),
    ("diagonal matrix", "对角矩阵", r"\diag(\cdot)"),
    ("diagonalization", "对角化", ""),
    ("dimension", "维数", ""),
    ("distance", "距离", ""),
    ("eigenvalue", "特征值", r"\lambda_i(\mA)"),
    ("eigenvector", "特征向量", ""),
    ("elementary matrix", "初等矩阵", r"\mE"),
    ("elementary row operations", "初等行变换", ""),
    ("entry (of a matrix)", "元素", r"a_{ij}"),
    ("eta column", "eta 列", r"\veta"),
    ("eta matrix", "eta 矩阵", r"\mE_k"),
    ("Euclidean length", "欧氏长度", r"\norm{\vx}_2"),
    ("Farkas' lemma", "Farkas 引理", ""),
    ("forward elimination", "前向消元", ""),
    ("four fundamental subspaces", "四个基本子空间", ""),
    ("Frobenius norm", "Frobenius 范数", r"\norm{\mA}_F"),
    ("full column rank", "列满秩", ""),
    ("full row rank", "行满秩", ""),
    ("Gauss–Jordan elimination", "高斯–若尔当消元法", ""),
    ("Gaussian elimination", "高斯消元法", ""),
    ("Gram–Schmidt process", "格拉姆–施密特正交化", ""),
    ("identity matrix", "单位矩阵", r"\mI"),
    ("ill-conditioned", "病态的", ""),
    ("indefinite", "不定", ""),
    ("induced matrix norm", "诱导矩阵范数", ""),
    ("inner product", "内积", r"\vx^\T\vy"),
    ("inverse matrix", "逆矩阵", r"\mA^{-1}"),
    ("invertible matrix", "可逆矩阵", ""),
    ("leading principal minor", "顺序主子式", r"D_k"),
    ("least-norm solution", "最小范数解", ""),
    ("left null space", "左零空间", r"\nullsp(\mA^\T)"),
    ("linear combination", "线性组合", ""),
    ("linear system", "线性方程组", r"\mA\vx=\vb"),
    ("linearly dependent", "线性相关", ""),
    ("linearly independent", "线性无关", ""),
    ("lower triangular matrix", "下三角矩阵", r"\mL"),
    ("LU factorization", "LU 分解", r"\mA=\mL\mU"),
    ("matrix", "矩阵", r"\mA"),
    ("minor", "子式", ""),
    ("negative definite", "负定", ""),
    ("negative semidefinite", "半负定", ""),
    ("nonsingular matrix", "非奇异矩阵", ""),
    ("norm", "范数", r"\norm{\vx}"),
    ("null space", "零空间", r"\nullsp(\mA)"),
    ("orthogonal", "正交", ""),
    ("orthogonal complement", "正交补", ""),
    ("orthonormal", "标准正交", ""),
    ("partial pivoting", "部分选主元", ""),
    ("partitioned (block) matrix", "分块矩阵", ""),
    ("pivot column", "主元列", r"s"),
    ("pivot element", "主元", r"\bar a_{rs}"),
    ("pivot matrix", "主元变换矩阵", r"\mE_k"),
    ("pivot operation", "主元变换", ""),
    ("pivot row", "主元行", r"r"),
    ("positive definite", "正定", r"\mA\succ\vzero"),
    ("positive semidefinite", "半正定", r"\mA\succeq\vzero"),
    ("positive semidefinite cone", "半正定锥", r"\Sym^n_+"),
    ("principal minor", "主子式", ""),
    ("principal submatrix", "主子矩阵", ""),
    ("product form of the inverse", "逆矩阵的乘积形式", ""),
    ("projection", "投影", ""),
    ("projection matrix", "投影矩阵", r"\mP"),
    ("quadratic form", "二次型", r"\vx^\T\mA\vx"),
    ("rank", "秩", r"\rank\mA"),
    ("rank–nullity theorem", "秩–零化度定理", ""),
    ("Rayleigh quotient", "瑞利商", ""),
    ("reduced row echelon form", "简化行阶梯形", ""),
    ("row space", "行空间", r"\row(\mA)"),
    ("scalar", "标量", ""),
    ("singular matrix", "奇异矩阵", ""),
    ("span", "张成", ""),
    ("spectral decomposition", "谱分解", r"\mA=\mV\mD\mV^\T"),
    ("spectral norm", "谱范数", r"\norm{\mA}_2"),
    ("spectral theorem", "谱定理", ""),
    ("square matrix", "方阵", ""),
    ("subspace", "子空间", ""),
    ("Sylvester's criterion", "西尔维斯特判据", ""),
    ("symmetric matrix", "对称矩阵", r"\Sym^n"),
    ("trace", "迹", r"\tr\mA"),
    ("transpose", "转置", r"\vx^\T,\ \mA^\T"),
    ("triangle inequality", "三角不等式", ""),
    ("unit vector", "单位向量", r"\ve_i"),
    ("upper triangular matrix", "上三角矩阵", r"\mU"),
    ("vector", "向量", r"\vx"),
    ("zero vector", "零向量", r"\vzero"),
    ("adjugate", "伴随矩阵", ""),
    ("Cramer's rule", "克拉默法则", ""),
    ("orthogonal matrix", "正交矩阵", ""),
    ("permutation matrix", "置换矩阵", ""),
    ("rectangle rule (pivoting)", "矩形法则", ""),
    ("Schur complement", "Schur 补", ""),
]

# ---------------------------------------------------------------------------
# Calculus (chapter 1.2)
# ---------------------------------------------------------------------------
CALCULUS: list[Term] = [
    ("affine map", "仿射映射", ""),
    ("central difference", "中心差分", ""),
    ("chain rule", "链式法则", ""),
    ("continuously differentiable", "连续可微", r"C^1"),
    ("contour (level line)", "等值线", ""),
    ("derivative with respect to a matrix", "对矩阵的导数", r"\partial f/\partial\mX"),
    ("differentiable", "可微", ""),
    ("directional derivative", "方向导数", r"f'(\vx;\vd)"),
    ("gradient", "梯度", r"\nabla f(\vx)"),
    ("Hessian matrix", "黑塞矩阵", r"\nabla^2 f(\vx)"),
    ("Jacobian matrix", "雅可比矩阵", r"D\vg(\vx)"),
    ("Landau symbols", "朗道符号", r"o(\cdot),\ O(\cdot)"),
    ("layout convention", "布局约定", ""),
    ("least squares", "最小二乘法", ""),
    ("level set", "水平集", ""),
    ("Lipschitz continuous", "Lipschitz 连续", ""),
    ("mean value theorem", "中值定理", ""),
    ("normal equations", "正规方程", r"\mA^\T\mA\vx=\mA^\T\vb"),
    ("partial derivative", "偏导数", r"\partial f/\partial x_j"),
    ("quadratic function", "二次函数", r"\tfrac12\vx^\T\mQ\vx+\vc^\T\vx"),
    ("restriction to a line", "在直线上的限制", r"\varphi(\alpha)=f(\vx+\alpha\vd)"),
    ("scalar function", "标量函数", ""),
    ("Taylor expansion", "泰勒展开", ""),
    ("Taylor's theorem", "泰勒定理", ""),
    ("twice continuously differentiable", "二阶连续可微", r"C^2"),
    ("vector function", "向量值函数", r"\vg(\vx)"),
    ("automatic differentiation", "自动微分", ""),
    ("finite difference", "有限差分", ""),
    ("Lagrange remainder", "拉格朗日余项", ""),
    ("numerator layout", "分子布局", ""),
    ("Rolle's theorem", "罗尔定理", ""),
]

# ---------------------------------------------------------------------------
# Optimization models (chapters 1.3 and 6.1)
# ---------------------------------------------------------------------------
OPTIMIZATION: list[Term] = [
    ("coercive function", "强制函数", ""),
    ("constraint", "约束", ""),
    ("decision variable", "决策变量", r"\vx"),
    ("deterministic", "确定性的", ""),
    ("equality constraint", "等式约束", r"h_j(\vx)=0"),
    ("feasible region", "可行域", ""),
    ("feasible set", "可行集", r"\Omega"),
    ("feasible solution", "可行解", ""),
    ("global minimum", "全局极小点", ""),
    ("inequality constraint", "不等式约束", r"g_i(\vx)\le0"),
    ("infeasible", "不可行", ""),
    ("infimum", "下确界", r"\inf"),
    ("iterate", "迭代点", r"\vx^{(k)}"),
    ("local maximum", "局部极大点", ""),
    ("local minimum", "局部极小点", ""),
    ("mathematical programming", "数学规划", ""),
    ("maximization", "最大化", ""),
    ("minimization", "最小化", ""),
    ("minimizer", "极小点", r"\argmin"),
    ("modeling process", "建模过程", ""),
    ("multiple objectives", "多目标", ""),
    ("objective function", "目标函数", r"f"),
    ("operations research (OR)", "运筹学", ""),
    ("optimal solution", "最优解", r"\vx^*"),
    ("optimal value", "最优值", r"p^*"),
    ("optimization", "最优化", ""),
    ("optimization model", "最优化模型", ""),
    ("parameter", "参数", ""),
    ("stochastic", "随机的", ""),
    ("stopping rule", "终止准则", r"\varepsilon"),
    ("strict local minimum", "严格局部极小点", ""),
    ("subject to", "满足（s.t.）", r"\st"),
    ("supremum", "上确界", r"\sup"),
    ("unbounded", "无界", ""),
    ("uncertainty", "不确定性", ""),
    ("all-or-nothing assignment", "全有全无分配", ""),
    ("decision science", "决策科学", ""),
    ("digital twin", "数字孪生", ""),
    ("free-flow travel time", "自由流行程时间", ""),
    ("linear expenditure system", "线性支出系统", ""),
    ("link performance function", "路段阻抗函数", ""),
    ("management science", "管理科学", ""),
    ("marginal propensity to spend", "边际支出倾向", ""),
    ("queueing theory", "排队论", ""),
    ("system-optimal assignment", "系统最优交通分配", ""),
    ("system-optimal routing", "系统最优路径分配", ""),
    ("traffic assignment", "交通分配", ""),
    ("unit commitment", "机组组合", ""),
]

# ---------------------------------------------------------------------------
# Convex analysis and Lagrange duality (part II)
# ---------------------------------------------------------------------------
CONVEX: list[Term] = [
    ("affine combination", "仿射组合", ""),
    ("affine hull", "仿射包", r"\aff S"),
    ("affine set", "仿射集", ""),
    ("ball", "球", ""),
    ("Bolzano–Weierstrass theorem", "波尔查诺–魏尔斯特拉斯定理", ""),
    ("boundary", "边界", ""),
    ("boundary point", "边界点", ""),
    ("bounded set", "有界集", ""),
    ("Carathéodory's theorem", "Carathéodory 定理", ""),
    ("chord", "弦", ""),
    ("closed ball", "闭球", ""),
    ("closed set", "闭集", ""),
    ("closure", "闭包", ""),
    ("compact set", "紧集", ""),
    ("concave function", "凹函数", ""),
    ("cone", "锥", r"K"),
    ("conic combination", "锥组合", ""),
    ("conic hull", "锥包", r"\cone S"),
    ("convergence (of a sequence)", "收敛", ""),
    ("convex combination", "凸组合", r"\theta_i"),
    ("convex cone", "凸锥", ""),
    ("convex function", "凸函数", ""),
    ("convex hull", "凸包", r"\conv S"),
    ("convex optimization problem", "凸优化问题", ""),
    ("convex program", "凸规划", ""),
    ("convex set", "凸集", r"C"),
    ("domain", "定义域", r"\dom f"),
    ("dual optimal value", "对偶最优值", r"d^*"),
    ("duality gap", "对偶间隙", r"p^*-d^*"),
    ("ellipsoid", "椭球", ""),
    ("epigraph", "上图", r"\operatorname{epi} f"),
    ("extended-value function", "扩展值函数", ""),
    ("extreme point", "极点", ""),
    ("half-space", "半空间", r"\{\vx:\va^\T\vx\le b\}"),
    ("hyperplane", "超平面", r"\{\vx:\va^\T\vx=b\}"),
    ("indicator function", "示性函数", ""),
    ("interior", "内部", r"\intr S"),
    ("interior point", "内点", ""),
    ("Jensen's inequality", "Jensen 不等式", ""),
    ("Lagrange dual function", "拉格朗日对偶函数", r"q(\vlambda,\vmu)"),
    ("Lagrange dual problem", "拉格朗日对偶问题", ""),
    ("Lagrange duality", "拉格朗日对偶", ""),
    ("Lagrange multiplier", "拉格朗日乘子", r"\lambda_i,\ \mu_j"),
    ("Lagrangian", "拉格朗日函数", r"L(\vx,\vlambda,\vmu)"),
    ("limit point", "极限点", ""),
    ("line", "直线", ""),
    ("line segment", "线段", ""),
    ("Minkowski–Weyl theorem", "Minkowski–Weyl 定理", ""),
    ("nonnegative orthant", "非负卦限", r"\R^n_+"),
    ("norm cone", "范数锥", ""),
    ("normal vector", "法向量", r"\va"),
    ("open ball", "开球", r"B(\vx,r)"),
    ("open set", "开集", ""),
    ("operations that preserve convexity", "保凸运算", ""),
    ("polyhedron", "多面体", r"P"),
    ("polytope", "多胞形", ""),
    ("projection onto a convex set", "到凸集上的投影", ""),
    ("quasiconcave function", "拟凹函数", ""),
    ("quasiconvex function", "拟凸函数", ""),
    ("relative interior", "相对内部", r"\relint S"),
    ("second-order cone", "二阶锥", ""),
    ("separating hyperplane", "分离超平面", ""),
    ("sequence", "序列", ""),
    ("Slater's condition", "Slater 条件", ""),
    ("strictly concave", "严格凹", ""),
    ("strictly convex", "严格凸", ""),
    ("strictly feasible", "严格可行", ""),
    ("strictly quasiconvex", "严格拟凸", ""),
    ("strong duality", "强对偶性", ""),
    ("sublevel set", "下水平集", r"\{\vx:f(\vx)\le\alpha\}"),
    ("subsequence", "子序列", ""),
    ("supporting hyperplane", "支撑超平面", ""),
    ("unimodal function", "单峰函数", ""),
    ("unit simplex", "单位单纯形", ""),
    ("weak duality", "弱对偶性", ""),
    ("Weierstrass theorem", "魏尔斯特拉斯定理", ""),
    ("arithmetic–geometric mean inequality", "算术–几何平均不等式", ""),
    ("hypograph", "亚图", ""),
    ("ice-cream cone", "冰淇淋锥", ""),
    ("Krein–Milman theorem", "Krein–Milman 定理", ""),
    ("linear-fractional map", "线性分式映射", ""),
    ("perspective map", "透视映射", ""),
    ("quasilinear", "拟线性的", ""),
    ("second-order cone program (SOCP)", "二阶锥规划", ""),
    ("semidefinite program (SDP)", "半定规划", ""),
]

# ---------------------------------------------------------------------------
# Linear programming and the simplex method (chapters 3.1–3.4, 3.6)
# ---------------------------------------------------------------------------
LINEAR_PROGRAMMING: list[Term] = [
    ("additivity", "可加性", ""),
    ("adjacent basic feasible solutions", "相邻基本可行解", ""),
    ("alternative optima", "多重最优解", ""),
    ("artificial variable", "人工变量", r"R_i"),
    ("assignment problem", "指派问题", ""),
    ("assumptions of linear programming", "线性规划的基本假设", ""),
    ("balance equation", "平衡方程", ""),
    ("basic feasible solution (BFS)", "基本可行解", ""),
    ("basic solution", "基本解", ""),
    ("basic variable", "基变量", r"\vx_B"),
    ("basis", "基", r"\mathcal B"),
    ("basis matrix", "基矩阵", r"\mB"),
    ("big-M method", "大 M 法", r"M"),
    ("binding constraint", "紧约束", ""),
    ("Bland's rule", "Bland 规则", ""),
    ("blending problem", "配料问题", ""),
    ("bordered basis matrix", "加边基矩阵", ""),
    ("canonical form", "典式", ""),
    ("certainty", "确定性", ""),
    ("constraint matrix", "约束矩阵", r"\mA"),
    ("cost vector", "费用向量", r"\vc"),
    ("cycling", "循环", ""),
    ("Dantzig's rule", "Dantzig 规则", ""),
    ("degenerate", "退化", ""),
    ("degenerate pivot", "退化的主元变换", ""),
    ("demand constraint", "需求约束", ""),
    ("diet problem", "营养配餐问题", ""),
    ("divisibility", "可分割性", ""),
    ("end-of-horizon effect", "期末效应", ""),
    ("entering variable", "进基变量", r"x_s"),
    ("finite termination", "有限终止", ""),
    ("free variable", "自由变量", r"x_j=x_j^+-x_j^-"),
    ("fundamental theorem of linear programming", "线性规划基本定理", ""),
    ("graphical method", "图解法", ""),
    ("leaving variable", "出基变量", r"x_{B(r)}"),
    ("lexicographic rule", "字典序规则", ""),
    ("linear program", "线性规划问题", ""),
    ("linear programming (LP)", "线性规划", ""),
    ("minimum-ratio test", "最小比值检验", r"\theta"),
    ("multi-period model", "多周期模型", ""),
    ("non-binding constraint", "非紧约束", ""),
    ("nonbasic matrix", "非基矩阵", r"\mN"),
    ("nonbasic variable", "非基变量", r"\vx_N"),
    ("nondegenerate", "非退化", ""),
    ("nonnegativity constraint", "非负约束", r"\vx\ge\vzero"),
    ("objective row", "目标行（$-z$ 行）", r"-z"),
    ("objective value", "目标函数值", r"z"),
    ("optimal basic feasible solution", "最优基本可行解", ""),
    ("optimality test", "最优性检验", ""),
    ("perturbation", "摄动", ""),
    ("Phase I", "第一阶段", r"\phi"),
    ("Phase II", "第二阶段", ""),
    ("pricing out", "检验数消元", ""),
    ("product mix", "产品组合", ""),
    ("production–inventory model", "生产–库存模型", ""),
    ("proportionality", "比例性", ""),
    ("redundant constraint", "冗余约束", ""),
    ("reduced cost", "检验数", r"\bar c_j"),
    ("reduced form", "约化形式", ""),
    ("reinversion", "重新求逆", ""),
    ("revised simplex method", "修正单纯形法", ""),
    ("right-hand side", "右端项", r"\vb"),
    ("simplex method", "单纯形法", ""),
    ("simplex multipliers", "单纯形乘子", r"\vpi"),
    ("simplex tableau", "单纯形表", r"\bar\mA"),
    ("slack variable", "松弛变量", r"s_i"),
    ("sparse", "稀疏的", ""),
    ("specification constraint", "规格约束", ""),
    ("standard form", "标准形", r"\mA\vx=\vb,\ \vx\ge\vzero"),
    ("supply constraint", "供应约束", ""),
    ("surplus variable", "剩余变量", r"s_i"),
    ("T-matrix (bordered inverse)", "T 矩阵（加边逆矩阵）", r"\mT"),
    ("tie (in a ratio test)", "并列", ""),
    ("transportation problem", "运输问题", ""),
    ("two-phase method", "两阶段法", ""),
    ("unbounded objective", "目标函数无界", ""),
    ("contribution (per unit)", "贡献", ""),
    ("detached coefficient tableau", "分离系数表", ""),
    ("edge direction", "棱方向", ""),
    ("partial pricing", "部分定价", ""),
    ("percentage specification", "百分比规格", ""),
    ("pricing (revised simplex)", "定价", ""),
    ("slack (of a row)", "松弛量", ""),
    ("split variables", "拆分变量", ""),
    ("support (of a vector)", "支撑集", "J(\\vx)"),
    ("surplus (of a row)", "剩余量", ""),
    ("temporary degeneracy", "暂时退化", ""),
    ("$z$-row", "$z$ 行", ""),
]

# ---------------------------------------------------------------------------
# Duality and sensitivity analysis (chapters 3.5–3.7)
# ---------------------------------------------------------------------------
DUALITY: list[Term] = [
    ("100 percent rule", "百分之百法则", r"\kappa_r"),
    ("activity", "活动", ""),
    ("artificial constraint", "人工约束", ""),
    ("break-even cost", "盈亏平衡成本", ""),
    ("canonical primal–dual pair", "对称形式的原始–对偶问题对", ""),
    ("change of a cost or right-hand side", "费用或右端项的改变量", r"\delta"),
    ("complementary slackness", "互补松弛性", ""),
    ("dual constraint", "对偶约束", ""),
    ("dual feasible", "对偶可行", ""),
    ("dual feasible basis", "对偶可行基", ""),
    ("dual objective value", "对偶目标函数值", r"w"),
    ("dual problem", "对偶问题", r"(\mathrm D)"),
    ("dual ratio test", "对偶比值检验", r"\theta_D"),
    ("dual simplex method", "对偶单纯形法", ""),
    ("dual slack", "对偶松弛变量", ""),
    ("dual variable", "对偶变量", r"\vpi"),
    ("global sensitivity inequality", "全局灵敏度不等式", ""),
    ("imputed cost", "隐含成本", r"\vpi^\T\mA_{\cdot j}"),
    ("marginal cost", "边际成本", ""),
    ("natural type (of a constraint)", "自然类型", ""),
    ("new activity", "新增活动", ""),
    ("new constraint", "新增约束", ""),
    ("optimality certificate", "最优性证书", ""),
    ("parametric right-hand side", "参数化右端项", ""),
    ("primal feasible", "原始可行", ""),
    ("primal feasible basis", "原始可行基", ""),
    ("primal problem", "原问题", r"(\mathrm P)"),
    ("primal simplex method", "原始单纯形法", ""),
    ("primal–dual pair", "原始–对偶问题对", ""),
    ("ranging", "变化范围分析", ""),
    ("re-optimization", "重新优化", ""),
    ("sensitivity analysis", "灵敏度分析", ""),
    ("shadow price", "影子价格", r"\pi_i"),
    ("theorem of the alternative", "择一定理", ""),
    ("validity range", "有效范围", ""),
    ("competitive prices", "有竞争力的价格", ""),
    ("dual value", "对偶值", ""),
    ("parametric analysis", "参数分析", ""),
]

# ---------------------------------------------------------------------------
# Network optimization (part IV)
# ---------------------------------------------------------------------------
NETWORKS: list[Term] = [
    ("$s$–$t$ cut", "$s$–$t$ 割", r"[S,\bar S]"),
    ("acyclic network", "无圈网络", ""),
    ("adjacency list", "邻接表", r"\mathcal A(i)"),
    ("adjacent nodes", "相邻节点", ""),
    ("all-pairs shortest paths", "所有点对最短路", ""),
    ("arc", "弧", r"(i,j)"),
    ("arc list", "弧表", ""),
    ("augmenting path", "增广路", r"\mathcal P"),
    ("backward arc", "后向弧", ""),
    ("Bellman's equations", "Bellman 方程", ""),
    ("Bellman–Ford algorithm", "Bellman–Ford 算法", r"V^{(k)}(j)"),
    ("bipartite matching", "二部图匹配", ""),
    ("Borůvka's algorithm", "Borůvka 算法", ""),
    ("breadth-first search", "广度优先搜索", ""),
    ("capacitated network", "容量网络", ""),
    ("capacity", "容量", r"u_{ij}"),
    ("chain", "链", ""),
    ("clustering", "聚类", ""),
    ("component (of a network)", "连通分支", ""),
    ("connected", "连通的", ""),
    ("cut", "割", r"[S,\bar S]"),
    ("cut capacity", "割容量", r"u[S,\bar S]"),
    ("cut property", "割性质", ""),
    ("cycle", "圈", r"\mathcal W"),
    ("cycle canceling", "消圈算法", ""),
    ("cycle property", "圈性质", ""),
    ("degree", "度", r"\deg(i)"),
    ("demand node", "需求节点", ""),
    ("destination", "终点", r"t"),
    ("Dijkstra's algorithm", "Dijkstra 算法", ""),
    ("directed cycle", "有向圈", ""),
    ("directed graph", "有向图", ""),
    ("directed network", "有向网络", ""),
    ("directed path", "有向路", ""),
    ("distance label", "距离标号", r"V(j)"),
    ("distance matrix", "距离矩阵", r"\mD_k"),
    ("dummy node", "虚拟节点", ""),
    ("edge", "边", ""),
    ("Edmonds–Karp rule", "Edmonds–Karp 规则", ""),
    ("feasible flow", "可行流", ""),
    ("Floyd–Warshall algorithm", "Floyd–Warshall 算法", ""),
    ("flow", "流", r"x_{ij}"),
    ("flow conservation", "流量守恒", ""),
    ("flow value", "流值", r"v"),
    ("flow-balance constraint", "流量平衡约束", ""),
    ("Ford–Fulkerson algorithm", "Ford–Fulkerson 算法", ""),
    ("forest", "森林", ""),
    ("forward arc", "前向弧", ""),
    ("forward star", "前向星", ""),
    ("graph", "图", ""),
    ("head (of an arc)", "弧头", r"j"),
    ("in-degree", "入度", ""),
    ("incident", "关联", ""),
    ("integrality of flows", "流的整数性", ""),
    ("Kirchhoff's matrix-tree theorem", "基尔霍夫矩阵树定理", ""),
    ("Kruskal's algorithm", "Kruskal 算法", ""),
    ("label", "标号", ""),
    ("label-correcting algorithm", "标号修正算法", ""),
    ("label-setting algorithm", "标号设定算法", ""),
    ("labeling procedure", "标号法", ""),
    ("Laplacian matrix", "拉普拉斯矩阵", r"\mL"),
    ("leaf", "叶节点", ""),
    ("length matrix", "长度矩阵", ""),
    ("longest path", "最长路", ""),
    ("lower bound (of an arc flow)", "流量下界", r"l_{ij}"),
    ("max-flow min-cut theorem", "最大流最小割定理", ""),
    ("maximum flow value", "最大流值", r"v^*"),
    ("maximum matching", "最大匹配", ""),
    ("maximum-flow problem", "最大流问题", ""),
    ("minimum cut", "最小割", ""),
    ("minimum spanning tree (MST)", "最小生成树", ""),
    ("minimum-cost flow problem", "最小费用流问题", ""),
    ("negative cycle", "负圈", ""),
    ("negative-cycle optimality conditions", "负圈最优性条件", ""),
    ("net supply", "净供应量", r"b_i"),
    ("network", "网络", r"G=(\mathcal N,\mathcal A)"),
    ("network simplex method", "网络单纯形法", ""),
    ("node", "节点", r"i\in\mathcal N"),
    ("node potential", "节点势", r"\pi_i"),
    ("node–arc incidence matrix", "节点–弧关联矩阵", r"\mA"),
    ("origin", "起点", r"s"),
    ("out-degree", "出度", ""),
    ("path", "路", ""),
    ("permanent label", "永久标号", r"S"),
    ("predecessor", "前驱节点", r"P(j)"),
    ("predecessor matrix", "前驱矩阵", r"\mP_k"),
    ("Prim's algorithm", "Prim 算法", r"[\gamma(j),P(j)]"),
    ("reduced cost of an arc", "弧的检验数", r"\bar c_{ij}"),
    ("reduced length", "约化长度", r"c_{ij}+V(i)-V(j)"),
    ("reduced-cost optimality conditions", "检验数最优性条件", ""),
    ("residual capacity", "残余容量", r"r_{ij}"),
    ("residual network", "残量网络", r"G(\vx)"),
    ("sequence matrix", "序列矩阵", ""),
    ("shortest-path problem", "最短路问题", ""),
    ("shortest-path tree", "最短路树", ""),
    ("single-source shortest paths", "单源最短路", ""),
    ("sink", "汇点", r"t"),
    ("source", "源点", r"s"),
    ("spanning tree", "生成树", ""),
    ("strongly connected", "强连通的", ""),
    ("sub-network", "子网络", ""),
    ("successive shortest paths", "连续最短路算法", ""),
    ("super-sink", "超级汇点", ""),
    ("super-source", "超级源点", ""),
    ("supply node", "供应节点", ""),
    ("tail (of an arc)", "弧尾", r"i"),
    ("temporary label", "临时标号", r"T"),
    ("topological order", "拓扑序", ""),
    ("totally unimodular matrix", "全单模矩阵", ""),
    ("transshipment node", "转运节点", ""),
    ("tree", "树", r"\mathcal T"),
    ("triple operation", "三元运算", ""),
    ("undirected network", "无向网络", ""),
    ("union–find", "并查集", ""),
    ("vertex", "顶点", ""),
    ("walk", "途径", ""),
    ("adjacency matrix", "邻接矩阵", ""),
    ("binary heap", "二叉堆", ""),
    ("breakthrough path", "突破路", ""),
    ("capacity scaling", "容量尺度化", ""),
    ("Cayley's formula", "Cayley 公式", ""),
    ("circulation", "循环流", ""),
    ("cross-dock", "越库配送中心", ""),
    ("depth-first search", "深度优先搜索", ""),
    ("dissimilarity", "相异度", ""),
    ("Hungarian method", "匈牙利算法", ""),
    ("König's theorem", "König 定理", ""),
    ("lazy deletion", "惰性删除", ""),
    ("matrix-tree theorem", "矩阵树定理", ""),
    ("minimum-cost arborescence", "最小费用树形图", ""),
    ("path halving", "路径减半", ""),
    ("preflow–push", "预流推进", ""),
    ("single-linkage clustering", "单链接聚类", ""),
    ("spacing (of a clustering)", "间隔", ""),
    ("strongly feasible tree", "强可行树", ""),
    ("time-expanded network", "时间扩展网络", ""),
]

# ---------------------------------------------------------------------------
# Dynamic programming (part V)
# ---------------------------------------------------------------------------
DYNAMIC_PROGRAMMING: list[Term] = [
    ("0–1 knapsack problem", "0–1 背包问题", ""),
    ("age (of a machine)", "役龄", r"i"),
    ("artificial stage variable", "人工阶段变量", ""),
    ("backward recursion", "逆序递推", ""),
    ("boundary condition", "边界条件", ""),
    ("buy (decision)", "购买", r"\mathrm B"),
    ("capital budgeting", "资本预算", ""),
    ("curse of dimensionality", "维数灾难", ""),
    ("Dantzig's recursion", "Dantzig 递推", r"\tilde D^{(k)}(i,j)"),
    ("decision", "决策", r"d"),
    ("decision set", "决策集合", r"\mathcal D_k(y)"),
    ("deterministic dynamic programming", "确定性动态规划", ""),
    ("dynamic programming (DP)", "动态规划", ""),
    ("equipment replacement", "设备更新", ""),
    ("feedback rule", "反馈规则", ""),
    ("Floyd's recursion", "Floyd 递推", r"D^{(k)}(i,j)"),
    ("forward recursion", "顺序递推", ""),
    ("forward value", "前向值", r"F_k(y)"),
    ("integer knapsack problem", "整数背包问题", ""),
    ("keep (decision)", "保留", r"\mathrm K"),
    ("knapsack problem", "背包问题", ""),
    ("lattice", "格点网络", r"(k,y)"),
    ("lease (decision)", "租赁", r"\mathrm L"),
    ("memoized recursion", "记忆化递归", ""),
    ("multistage decision problem", "多阶段决策问题", ""),
    ("number of stages", "阶段数", r"N"),
    ("one-pass recursion", "单遍递推", ""),
    ("operating cost", "运行费用", r"c(i)"),
    ("optimal policy function", "最优策略函数", r"P_k(y)"),
    ("optimal value function", "最优值函数", r"V_k(y)"),
    ("overhaul (decision)", "大修", r"\mathrm O"),
    ("policy", "策略", ""),
    ("principle of optimality", "最优性原理", ""),
    ("recurrence relation", "递推关系", ""),
    ("resource allocation", "资源分配", r"x,\ X"),
    ("return function", "收益函数", r"r_k(x)"),
    ("salvage value", "残值", r"s(i)"),
    ("separability", "可分性", ""),
    ("stage", "阶段", r"k"),
    ("stage cost", "阶段费用", r"c_k(y,d)"),
    ("staged network", "阶段网络", ""),
    ("state", "状态", r"y"),
    ("state augmentation", "状态扩充", ""),
    ("state set", "状态集合", r"\mathcal Y_k"),
    ("stochastic dynamic programming", "随机动态规划", ""),
    ("successive approximation", "逐次逼近法", ""),
    ("terminal cost", "终端费用", r"g(y)"),
    ("trade-in value", "旧设备折价", r"t(i)"),
    ("transition function", "状态转移函数", r"t_k(y,d)"),
    ("turn penalty", "转弯惩罚", ""),
    ("value function", "值函数", ""),
    ("answer cell", "所求值", ""),
    ("decrease (successive approximation)", "下降", ""),
    ("goodwill", "商誉", ""),
    ("Markov (no-memory) property", "Markov 性（无记忆性）", ""),
    ("reversal (successive approximation)", "转向", ""),
    ("run (monotone piece of a walk)", "单调段", ""),
    ("sunk cost", "沉没成本", ""),
    ("Tabourier's modification", "Tabourier 修正", ""),
]

# ---------------------------------------------------------------------------
# Nonlinear programming and optimality conditions (chapters 6.1–6.3)
# ---------------------------------------------------------------------------
NONLINEAR: list[Term] = [
    ("Abadie's constraint qualification", "Abadie 约束规范", ""),
    ("active constraint", "有效约束", ""),
    ("active set", "有效集", r"I(\vx)"),
    ("active-set case analysis", "有效集分情形讨论", ""),
    ("bound-constrained problem", "界约束问题", ""),
    ("cone of feasible directions", "可行方向锥", r"D(\vx)"),
    ("cone of linearized feasible directions", "线性化可行方向锥", r"L(\vx^*)"),
    ("constrained problem", "约束问题", ""),
    ("constraint qualification (CQ)", "约束规范", ""),
    ("critical cone", "临界锥", r"C(\vx^*,\vlambda^*)"),
    ("descent direction", "下降方向", r"\vd"),
    ("dual feasibility", "对偶可行性", r"\vlambda\ge\vzero"),
    ("economic dispatch", "经济调度", ""),
    ("feasible descent direction", "可行下降方向", ""),
    ("feasible direction", "可行方向", ""),
    ("first-order necessary condition", "一阶必要条件", ""),
    ("Fritz John conditions", "Fritz John 条件", r"\lambda_0"),
    ("gradient of the Lagrangian", "拉格朗日函数的梯度", r"\nabla_{\vx}L"),
    ("Hessian of the Lagrangian", "拉格朗日函数的黑塞矩阵", r"\nabla^2_{\vx\vx}L"),
    ("inactive constraint", "非有效约束", ""),
    ("Karush–Kuhn–Tucker (KKT) conditions", "KKT 条件", ""),
    ("KKT point", "KKT 点", ""),
    ("least-squares estimation", "最小二乘估计", ""),
    ("linear independence constraint qualification (LICQ)", "线性无关约束规范", ""),
    ("linearly constrained problem", "线性约束问题", ""),
    ("nonlinear program", "非线性规划问题", ""),
    ("nonlinear programming (NLP)", "非线性规划", ""),
    ("primal feasibility", "原始可行性", ""),
    ("quadratic programming (QP)", "二次规划", ""),
    ("regular point", "正则点", ""),
    ("saddle point", "鞍点", ""),
    ("second-order necessary condition", "二阶必要条件", ""),
    ("second-order sufficient condition", "二阶充分条件", ""),
    ("stationarity", "平稳性条件", r"\nabla_{\vx}L=\vzero"),
    ("stationary point", "驻点", ""),
    ("strict complementarity", "严格互补性", ""),
    ("tangent cone", "切锥", r"T(\vx^*)"),
    ("tangent direction", "切方向", ""),
    ("unconstrained problem", "无约束问题", ""),
    ("variational inequality", "变分不等式", ""),
    ("degenerate constraint", "退化约束", ""),
    ("egg-crate surface", "蛋托曲面", ""),
    ("Mangasarian–Fromovitz condition", "Mangasarian–Fromovitz 条件", ""),
    ("omega constant", "欧米加常数", ""),
]

# ---------------------------------------------------------------------------
# Line search and algorithms (chapters 6.4–6.6)
# ---------------------------------------------------------------------------
ALGORITHMS: list[Term] = [
    ("affine invariance", "仿射不变性", ""),
    ("algorithmic map", "算法映射", ""),
    ("Armijo condition", "Armijo 条件", r"c_1"),
    ("augmented Lagrangian", "增广拉格朗日函数", r"L_A(\vx,\vmu;\rho)"),
    ("augmented Lagrangian method", "增广拉格朗日法", ""),
    ("away step", "远离步", ""),
    ("backtracking line search", "回溯线搜索", r"\gamma"),
    ("barrier method", "障碍函数法", ""),
    ("barrier parameter", "障碍参数", r"\eta"),
    ("BFGS update", "BFGS 修正公式", ""),
    ("bisection search", "二分法", ""),
    ("Cauchy's method", "柯西法", ""),
    ("conditional-gradient method", "条件梯度法", ""),
    ("conjugate directions", "共轭方向", ""),
    ("conjugate-gradient method", "共轭梯度法", ""),
    ("convex-combination method", "凸组合法", ""),
    ("curvature condition", "曲率条件", r"c_2"),
    ("curve fitting", "曲线拟合法", ""),
    ("descent algorithm", "下降算法", ""),
    ("descent function", "下降函数", ""),
    ("descent method", "下降法", ""),
    ("DFP update", "DFP 修正公式", ""),
    ("diagonal scaling", "对角尺度化", ""),
    ("exact line search", "精确线搜索", ""),
    ("exact penalty function", "精确罚函数", ""),
    ("feasible-direction method", "可行方向法", ""),
    ("Fibonacci search", "斐波那契法", r"F_N"),
    ("Frank–Wolfe gap", "Frank–Wolfe 间隙", r"\Gamma(\vx)"),
    ("Frank–Wolfe method", "Frank–Wolfe 方法", ""),
    ("global convergence", "全局收敛性", ""),
    ("golden-section search", "黄金分割法", r"\tau"),
    ("gradient projection method", "梯度投影法", ""),
    ("Hessian modification", "黑塞矩阵修正", r"\bar\mH_k"),
    ("inexact line search", "非精确线搜索", ""),
    ("interior-point method", "内点法", ""),
    ("interval elimination", "区间消去", ""),
    ("interval of uncertainty", "不确定区间", r"[a_k,b_k]"),
    ("iterative algorithm", "迭代算法", ""),
    ("line search", "线搜索", ""),
    ("linear convergence", "线性收敛", ""),
    ("linearized subproblem", "线性化子问题", ""),
    ("local convergence", "局部收敛性", ""),
    ("Newton equations", "牛顿方程", r"\nabla^2 f(\vx^{(k)})\,\vd=-\vg^{(k)}"),
    ("Newton step", "牛顿步", ""),
    ("Newton's method", "牛顿法", ""),
    ("one-dimensional search", "一维搜索", ""),
    ("order of convergence", "收敛阶", ""),
    ("penalized objective", "带罚目标函数", r"\Phi_\rho(\vx)"),
    ("penalty function", "罚函数", r"\psi(\vx)"),
    ("penalty method", "罚函数法", ""),
    ("penalty parameter", "罚参数", r"\rho"),
    ("point-to-point map", "点到点映射", ""),
    ("point-to-set map", "点到集映射", ""),
    ("primal–dual interior-point method", "原始–对偶内点法", ""),
    ("quadratic convergence", "二次收敛", ""),
    ("quadratic interpolation", "二次插值法", ""),
    ("quadratic penalty function", "二次罚函数", ""),
    ("quasi-Newton method", "拟牛顿法", ""),
    ("quasi-Newton update", "拟牛顿修正", r"\vs^{(k)},\ \vy^{(k)}"),
    ("rate of convergence", "收敛速度", ""),
    ("reduced-gradient method", "既约梯度法", ""),
    ("safeguard", "保护措施", ""),
    ("scaled steepest descent", "尺度化最速下降法", ""),
    ("scaling matrix", "尺度矩阵", r"\mD_k"),
    ("search direction", "搜索方向", r"\vd^{(k)}"),
    ("secant condition", "割线条件", r"\mD_{k+1}\vy^{(k)}=\vs^{(k)}"),
    ("sequential quadratic programming (SQP)", "序列二次规划", ""),
    ("steepest descent method", "最速下降法", ""),
    ("steepest-descent direction", "最速下降方向", r"-\nabla f(\vx)"),
    ("step length", "步长", r"\alpha_k"),
    ("strong Wolfe conditions", "强 Wolfe 条件", ""),
    ("sublinear convergence", "次线性收敛", ""),
    ("superlinear convergence", "超线性收敛", ""),
    ("trust region", "信赖域", ""),
    ("unimodality", "单峰性", ""),
    ("warm start", "热启动", ""),
    ("Wolfe conditions", "Wolfe 条件", ""),
    ("zigzagging", "锯齿现象", ""),
    ("Zoutendijk's method", "Zoutendijk 方法", ""),
    ("bracket (of a minimum)", "括区间", ""),
    ("Brent's method", "Brent 方法", ""),
    ("certified lower bound", "可证明的下界", ""),
    ("damped Newton method", "阻尼牛顿法", ""),
    ("Kantorovich inequality", "Kantorovich 不等式", ""),
    ("limited-memory BFGS (L-BFGS)", "有限内存 BFGS", ""),
    ("linear-approximation method", "线性近似法", ""),
    ("merit function", "价值函数", ""),
    ("method of multipliers", "乘子法", ""),
    ("modified Cholesky factorization", "修正 Cholesky 分解", ""),
    ("multi-start", "多初始点", ""),
    ("multiplier estimate", "乘子估计", ""),
    ("$\\mQ$-conjugate", "$\\mQ$-共轭", ""),
]

# ---------------------------------------------------------------------------
# Linearization and mixed-integer programming (part VII)
# ---------------------------------------------------------------------------
LINEARIZATION: list[Term] = [
    ("absolute value", "绝对值", r"\abs{x}"),
    ("aggregated formulation", "聚合形式", ""),
    ("at least k of m constraints", "m 个约束中至少 k 个成立", ""),
    ("auxiliary variable", "辅助变量", r"\eta"),
    ("big-M constant", "大 M 常数", r"M_i"),
    ("big-M constraint", "大 M 约束", ""),
    ("bilinear term", "双线性项", r"xy"),
    ("binary variable", "0–1 变量", r"\delta_i\in\{0,1\}"),
    ("breakpoint", "分段点", r"a_i"),
    ("Charnes–Cooper transformation", "Charnes–Cooper 变换", ""),
    ("conjunctive normal form", "合取范式", ""),
    ("convex-combination formulation", "凸组合形式", ""),
    ("disaggregated formulation", "分解形式", ""),
    ("economies of scale", "规模经济", ""),
    ("either–or constraint", "择一约束", ""),
    ("epigraph form", "上图形式", ""),
    ("fixed charge", "固定费用", r"F_i"),
    ("if–then constraint", "若–则约束", ""),
    ("incremental formulation", "增量形式", r"\xi_k"),
    ("indicator constraint", "指示约束", ""),
    ("indicator variable", "指示变量", r"\delta_i"),
    ("integer programming", "整数规划", ""),
    ("interpolation error", "插值误差", ""),
    ("interpolation weights", "插值权重", r"\alpha_i"),
    ("least absolute deviations", "最小绝对偏差", ""),
    ("linear-fractional program", "线性分式规划", ""),
    ("linearization", "线性化", ""),
    ("logical constraint", "逻辑约束", ""),
    ("LP relaxation", "线性规划松弛", ""),
    ("maximin objective", "极大化极小目标", ""),
    ("McCormick envelope", "McCormick 包络", ""),
    ("McCormick inequalities", "McCormick 不等式", ""),
    ("minimax objective", "极小化极大目标", ""),
    ("mixed-integer linear program (MILP)", "混合整数线性规划", ""),
    ("multiple-choice formulation", "多重选择形式", ""),
    ("over-estimator", "上估计函数", ""),
    ("piecewise-linear approximation", "分段线性逼近", r"f^a"),
    ("piecewise-linear function", "分段线性函数", ""),
    ("piecewise-linear interpolant", "分段线性插值函数", ""),
    ("positive and negative parts", "正部与负部", r"\vx^+,\ \vx^-"),
    ("precedence constraint", "先后次序约束", ""),
    ("product of binary variables", "0–1 变量之积", r"\omega=\delta_A\delta_B"),
    ("quadratic 0–1 program", "0–1 二次规划", ""),
    ("segment (of a piecewise-linear function)", "区段", ""),
    ("semi-continuous variable", "半连续变量", ""),
    ("separable function", "可分函数", ""),
    ("separable programming", "可分规划", ""),
    ("special ordered set (SOS)", "特殊有序集", ""),
    ("special ordered set of type 1 (SOS1)", "第一类特殊有序集", ""),
    ("special ordered set of type 2 (SOS2)", "第二类特殊有序集", ""),
    ("split form", "拆分形式", ""),
    ("under-estimator", "下估计函数", ""),
    ("valid inequality", "有效不等式", ""),
    ("weighted median", "加权中位数", ""),
    ("clause (logic)", "子句", ""),
    ("concave envelope", "凹包络", ""),
    ("convex envelope", "凸包络", ""),
    ("disjunctive programming", "析取规划", ""),
    ("exclusive or", "异或", ""),
    ("facility location problem", "设施选址问题", ""),
    ("literal (logic)", "文字", ""),
    ("locally ideal (formulation)", "局部理想的", ""),
    ("one-dimensional method (bivariate interpolation)", "一维方法", ""),
    ("partitioning (McCormick envelope)", "分割", ""),
    ("rectangle method", "矩形方法", ""),
    ("restricted basis entry", "限制进基", ""),
    ("tangent-plane method", "切平面方法", ""),
    ("triangle method", "三角形方法", ""),
]

# ---------------------------------------------------------------------------
# Software and solvers (part IX). Names of programs (Python, SciPy, HiGHS,
# PuLP, CVXPY, NetworkX, CBC, Gurobi, ...) are not translated.
# ---------------------------------------------------------------------------
SOFTWARE: list[Term] = [
    ("academic license", "学术许可", ""),
    ("algebraic form", "代数形式", ""),
    ("algebraic modeling language", "代数建模语言", ""),
    ("array", "数组", ""),
    ("branch and bound", "分支定界法", ""),
    ("commercial solver", "商业求解器", ""),
    ("conic program", "锥规划", ""),
    ("crossover", "交叉", ""),
    ("cutting plane", "割平面", ""),
    ("disciplined convex programming (DCP)", "规范凸规划", ""),
    ("floating-point arithmetic", "浮点运算", ""),
    ("heap", "堆", ""),
    ("integrality tolerance", "整数性容差", ""),
    ("KKT checker", "KKT 检验程序", ""),
    ("matrix form", "矩阵形式", ""),
    ("modeling language", "建模语言", ""),
    ("modeling layer", "建模层", ""),
    ("open-source solver", "开源求解器", ""),
    ("optimality tolerance", "最优性容差", ""),
    ("primal feasibility tolerance", "原始可行性容差", ""),
    ("relative MIP gap", "相对 MIP 间隙", ""),
    ("running time", "运行时间", ""),
    ("solver", "求解器", ""),
    ("status (of a solve)", "求解状态", ""),
    ("tolerance", "容差", ""),
    ("vector form", "向量形式", ""),
    ("virtual environment", "虚拟环境", ""),
    ("alternating direction method of multipliers (ADMM)", "交替方向乘子法", ""),
    ("branch and cut", "分支切割", ""),
    ("marginal (dual value in SciPy)", "边际值", ""),
    ("mixed-integer nonlinear program (MINLP)", "混合整数非线性规划", ""),
    ("modeling framework", "建模框架", ""),
    ("presolve", "预处理", ""),
    ("primal heuristic", "原始启发式算法", ""),
    ("residual (of a constraint)", "残差", ""),
]

# The sections in page order: (English title, 中文标题, terms).
SECTIONS: list[Section] = [
    ("General terms", "通用术语", GENERAL),
    ("Linear algebra", "线性代数", LINEAR_ALGEBRA),
    ("Calculus", "微积分", CALCULUS),
    ("Optimization models", "最优化模型", OPTIMIZATION),
    ("Convex analysis and Lagrange duality", "凸分析与拉格朗日对偶", CONVEX),
    ("Linear programming and the simplex method", "线性规划与单纯形法", LINEAR_PROGRAMMING),
    ("Duality and sensitivity analysis", "对偶理论与灵敏度分析", DUALITY),
    ("Network optimization", "网络优化", NETWORKS),
    ("Dynamic programming", "动态规划", DYNAMIC_PROGRAMMING),
    ("Nonlinear programming and optimality conditions", "非线性规划与最优性条件", NONLINEAR),
    ("Line search and algorithms", "线搜索与算法", ALGORITHMS),
    ("Linearization and mixed-integer programming", "线性化与混合整数规划", LINEARIZATION),
    ("Software and solvers", "软件与求解器", SOFTWARE),
]

INTRO_EN = """\
This glossary lists the technical terms of the book with the Simplified
Chinese terms used in the Chinese edition and, where the book has one, the
symbol. The terms are grouped by topic and sorted alphabetically within each
group. Symbols that are used across chapters are explained on the
{doc}`notation` page."""

INTRO_ZH = """\
本术语表按主题列出本书的专业术语、对应的英文术语以及（如有）所用符号。各组内按英文术语的字母顺序排列，与英文版术语表的顺序一致。跨章节使用的符号见{doc}`notation`页。"""


def sort_key(term: Term) -> str:
    """Return the alphabetical sort key of a term.

    Parameters
    ----------
    term : Term
        A tuple ``(English, 中文, symbol)``.

    Returns
    -------
    str
        The English term in lower case, without accents, dollar signs and
        backslashes, and without a leading "the ".
    """
    text = unicodedata.normalize("NFKD", term[0])
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    text = text.replace("$", "").replace("\\", "").lower()
    return text.removeprefix("the ")


def check(sections: list[Section]) -> None:
    """Check the data before writing anything.

    Parameters
    ----------
    sections : list of Section
        The glossary sections.

    Raises
    ------
    ValueError
        If an English term occurs twice, a field is empty or carries
        surrounding blanks, or a cell contains ``|`` (which would break the
        Markdown table).
    """
    seen: dict[str, str] = {}
    for title, _, terms in sections:
        for english, chinese, symbol in terms:
            key = english.lower()
            if key in seen:
                raise ValueError(f"term {english!r} is in {seen[key]!r} and in {title!r}")
            seen[key] = title
            if not english or not chinese:
                raise ValueError(f"empty field in {(english, chinese, symbol)!r}")
            for cell in (english, chinese, symbol):
                if cell != cell.strip():
                    raise ValueError(f"surrounding blanks in {cell!r}")
                if "|" in cell:
                    raise ValueError(f"'|' in {cell!r}; use \\abs or \\norm")


def math(symbol: str) -> str:
    """Return the symbol as inline math, or an empty cell when there is none."""
    return f"${symbol}$" if symbol else ""


def render(sections: list[Section], lang: str) -> str:
    """Render the glossary page of one edition.

    Parameters
    ----------
    sections : list of Section
        The glossary sections.
    lang : {"en", "zh"}
        The edition. English pages use the columns ``English | 中文 |
        Symbol``; Chinese pages use ``中文 | English | 符号``.

    Returns
    -------
    str
        The MyST Markdown source of the page.
    """
    if lang == "en":
        lines = ["# Glossary", "", INTRO_EN, ""]
        header = ["| English | 中文 | Symbol |", "|---|---|---|"]
    else:
        lines = ["# 术语表", "", INTRO_ZH, ""]
        header = ["| 中文 | English | 符号 |", "|---|---|---|"]

    for title_en, title_zh, terms in sections:
        lines += [f"## {title_en if lang == 'en' else title_zh}", ""]
        lines += header
        for english, chinese, symbol in sorted(terms, key=sort_key):
            if lang == "en":
                lines.append(f"| {english} | {chinese} | {math(symbol)} |")
            else:
                lines.append(f"| {chinese} | {english} | {math(symbol)} |")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    """Check the data and write ``docs/en/glossary.md`` and ``docs/zh/glossary.md``."""
    check(SECTIONS)
    for lang in ("en", "zh"):
        outdir = DOCS / lang
        outdir.mkdir(parents=True, exist_ok=True)  # creates docs/zh/ if missing
        path = outdir / "glossary.md"
        # newline="\n": the same bytes on Windows (no CRLF), so running the
        # script never changes a checked-out page.
        path.write_text(render(SECTIONS, lang), encoding="utf-8", newline="\n")
        print(f"wrote {path.relative_to(ROOT)}")
    counts = ", ".join(f"{title}: {len(terms)}" for title, _, terms in SECTIONS)
    print(f"{sum(len(t) for _, _, t in SECTIONS)} terms ({counts})")


if __name__ == "__main__":
    main()
