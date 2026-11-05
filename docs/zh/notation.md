# 符号表

全书使用统一的符号。向量用粗体小写字母表示（$\vx$, $\vc$, $\vpi$），矩阵用粗体大写字母表示（$\mA$, $\mB$）；标量以及向量和矩阵的分量用斜体表示（$x_j$, $a_{ij}$）。所有向量都是列向量：$\vx^\T$ 表示转置，因此 $\vc^\T\vx$ 是一个数量积，梯度 $\nabla f(\vx)$ 也是列向量。向量之间的不等式，如 $\vx\ge\vzero$ 或 $\mA\vx\le\vb$，按分量成立；符号 $\succeq$ 与 $\succ$ 专用于表示对称矩阵的（半）正定性。特征值总是连同其矩阵一起写出，即 $\lambda_i(\mA)$，以便单独的 $\lambda_i$ 留给拉格朗日乘子使用。下列各表列出在多个章节中使用、或读者可能需要查阅的符号；只属于某一章或某一例题的符号在引入之处定义（理论章节也会在各自的“符号说明”一节中列出）。最后一列给出符号首次引入的章节。

## 一般符号与线性代数

| 符号 | 含义 | 首次出现 |
|---|---|---|
| $\R^n$, $\R^{m\times n}$ | 实 $n$ 维向量，实 $m\times n$ 矩阵 | {doc}`01-01-linear-algebra` |
| $\vx$, $x_j$ | 向量及其第 $j$ 个分量；$\vx_1,\dots,\vx_k$（带下标的粗体）表示一组向量 | {doc}`01-01-linear-algebra` |
| $\vx^\T$, $\mA^\T$, $\mA^{-\T}$ | 转置；转置的逆 | {doc}`01-01-linear-algebra` |
| $\vzero$, $\vone$, $\ve_i$ | 零向量，全一向量，第 $i$ 个单位向量 | {doc}`01-01-linear-algebra` |
| $\mI$, $\mI_m$ | （$m$ 阶）单位矩阵 | {doc}`01-01-linear-algebra` |
| $a_{ij}$, $\mA_{\cdot j}$, $\mA_{i\cdot}$ | $\mA$ 的元素、第 $j$ 列和第 $i$ 行 | {doc}`01-01-linear-algebra` |
| $\col(\mA)$, $\row(\mA)$, $\nullsp(\mA)$ | 列空间，行空间，零空间 | {doc}`01-01-linear-algebra` |
| $\rank\mA$, $\det\mA$, $\tr\mA$ | 秩，行列式，迹 | {doc}`01-01-linear-algebra` |
| $\diag(\cdot)$ | 以给定元素为对角线的对角矩阵 | {doc}`01-01-linear-algebra` |
| $\mE$, $\mE_k$; $\veta$ | 初等矩阵（主元变换矩阵、eta 矩阵）；其 eta 列 | {doc}`01-01-linear-algebra` |
| $\bar\mA$, $\bar a_{rs}$ | 当前单纯形表；位于第 $r$ 行、第 $s$ 列的主元 | {doc}`01-01-linear-algebra` |
| $\mL$, $\mU$ | 下三角因子与上三角因子，$\mA=\mL\mU$ | {doc}`01-01-linear-algebra` |
| $\lambda_i(\mA)$, $\lambda_{\min}(\mA)$, $\lambda_{\max}(\mA)$ | $\mA$ 的特征值 | {doc}`01-01-linear-algebra` |
| $\mV$, $\mD$ | 特征向量矩阵与特征值对角矩阵，$\mA=\mV\mD\mV^\T$ | {doc}`01-01-linear-algebra` |
| $\Sym^n$, $\Sym^n_+$, $\Sym^n_{++}$ | $n\times n$ 对称矩阵、半正定矩阵、正定矩阵的集合 | {doc}`01-01-linear-algebra` |
| $\mA\succeq\vzero$, $\mA\succ\vzero$ | $\mA$ 半正定，正定 | {doc}`01-01-linear-algebra` |
| $D_k$ | 第 $k$ 个顺序主子式 | {doc}`01-01-linear-algebra` |
| $\norm{\cdot}_1$, $\norm{\cdot}_2$, $\norm{\cdot}_\infty$, $\norm{\cdot}_F$ | 向量范数及其诱导矩阵范数，Frobenius 范数；$\norm{\cdot}=\norm{\cdot}_2$ | {doc}`01-01-linear-algebra` |
| $\kappa(\mA)$ | 条件数 $\norm{\mA}\,\norm{\mA^{-1}}$ | {doc}`01-01-linear-algebra` |
| $\mP$ | 投影矩阵 | {doc}`01-01-linear-algebra` |
| $\abs{S}$ | 有限集 $S$ 的元素个数 | {doc}`03-03-simplex-method` |
| $\vx^{(k)}$ | 序列的第 $k$ 项；算法的第 $k$ 个迭代点 | {doc}`02-01-convex-sets` |
| $\varepsilon$ | 小正数；终止准则的容差 | {doc}`01-02-matrix-calculus` |

## 微积分

| 符号 | 含义 | 首次出现 |
|---|---|---|
| $\nabla f(\vx)$, $\nabla^2 f(\vx)$ | $f$ 的梯度（列向量）与黑塞矩阵 | {doc}`01-02-matrix-calculus` |
| $D\vg(\vx)$ | $\vg:\R^n\to\R^m$ 的雅可比矩阵，即以 $\nabla g_i(\vx)^\T$ 为行的 $m\times n$ 矩阵 | {doc}`01-02-matrix-calculus` |
| $f'(\vx;\vd)$ | $f$ 在 $\vx$ 处沿 $\vd$ 的方向导数 | {doc}`01-02-matrix-calculus` |
| $C^1$, $C^2$ | 一阶、二阶连续可微 | {doc}`01-02-matrix-calculus` |
| $\bar\vx$ | 泰勒公式的展开点；中间点为 $\bar\vx+t\vd$，$0<t<1$（自{doc}`06-02-optimality-conditions`起记作 $\bar\vx+\vartheta\vd$） | {doc}`01-02-matrix-calculus` |
| $o(\cdot)$, $O(\cdot)$ | 朗道符号 | {doc}`01-02-matrix-calculus` |
| $\mX$, $\partial f/\partial\mX$ | 矩阵变量及对它的导数 | {doc}`01-02-matrix-calculus` |
| $\mQ$ | 二次函数 $\tfrac12\vx^\T\mQ\vx+\vc^\T\vx$ 的对称矩阵 | {doc}`01-02-matrix-calculus` |
| $\varphi(\alpha)$ | $f$ 在直线上的限制，$\varphi(\alpha)=f(\vx+\alpha\vd)$ | {doc}`01-02-matrix-calculus` |

## 凸分析与对偶

| 符号 | 含义 | 首次出现 |
|---|---|---|
| $f$, $g_i$, $h_j$ | 目标函数，不等式约束 $g_i(\vx)\le0$（$i=1,\dots,m$），等式约束 $h_j(\vx)=0$（$j=1,\dots,\ell$） | {doc}`01-02-matrix-calculus` |
| $\vg(\vx)$, $\vh(\vx)$ | 约束函数构成的向量 | {doc}`01-02-matrix-calculus` |
| $m$, $\ell$ | 不等式约束的个数与等式约束的个数 | {doc}`01-02-matrix-calculus` |
| $\Omega$ | 最优化问题的可行集 | {doc}`01-03-introduction` |
| $\vx^*$ | 最优解（极小点） | {doc}`01-02-matrix-calculus` |
| $p^*$ | 最小化问题的最优值 | {doc}`01-03-introduction` |
| $C$, $K$ | 一个（凸）集；一个锥 | {doc}`02-01-convex-sets` |
| $\theta$, $\theta_i$ | 仿射组合、凸组合或锥组合的权重 | {doc}`02-01-convex-sets` |
| $\aff S$, $\conv S$, $\cone S$ | $S$ 的仿射包、凸包与锥包 | {doc}`02-01-convex-sets` |
| $\R^n_+$ | 非负卦限 | {doc}`02-01-convex-sets` |
| $B(\vx,r)$ | 以 $\vx$ 为中心、半径为 $r$ 的开球 | {doc}`02-01-convex-sets` |
| $\intr S$ | $S$ 的内部 | {doc}`02-01-convex-sets` |
| $\dom f$ | $f$ 的定义域 | {doc}`02-02-convex-functions` |
| $\operatorname{epi} f$ | $f$ 的上图 | {doc}`02-02-convex-functions` |
| $L(\vx,\vlambda,\vmu)$ | 拉格朗日函数 $f(\vx)+\vlambda^\T\vg(\vx)+\vmu^\T\vh(\vx)$ | {doc}`01-02-matrix-calculus` |
| $\lambda_i$, $\vlambda$ | 不等式约束 $g_i\le0$ 的乘子（$\lambda_i\ge0$） | {doc}`01-02-matrix-calculus` |
| $\mu_j$, $\vmu$ | 等式约束 $h_j=0$ 的乘子（符号不限） | {doc}`01-02-matrix-calculus` |
| $q(\vlambda,\vmu)$ | 拉格朗日对偶函数 | {doc}`02-03-convex-optimization-duality` |
| $d^*$ | 对偶问题的最优值 | {doc}`02-03-convex-optimization-duality` |

## 线性规划

| 符号 | 含义 | 首次出现 |
|---|---|---|
| $\mA\in\R^{m\times n}$, $\vb\in\R^m$, $\vc\in\R^n$ | 约束矩阵，右端项，费用（或利润）向量；$m$ 个约束，$n$ 个变量 | {doc}`03-01-lp-model-building` |
| $z$, $z^*$ | （原）线性规划的目标函数值 $\vc^\T\vx$；其最优值 | {doc}`03-01-lp-model-building` |
| $s_i$ | 约束 $i$ 的松弛变量或剩余变量 | {doc}`03-02-lp-geometry` |
| $x_j^+$, $x_j^-$ | 自由变量的正部与负部 | {doc}`03-02-lp-geometry` |
| $P$ | 可行多面体 $\{\vx:\mA\vx=\vb,\ \vx\ge\vzero\}$ | {doc}`03-02-lp-geometry` |
| $\mB$, $\mN$ | 基矩阵与非基矩阵 | {doc}`01-01-linear-algebra` |
| $\vx_B$, $\vx_N$, $\vc_B$, $\vc_N$ | $\vx$ 与 $\vc$ 的基部分与非基部分 | {doc}`01-01-linear-algebra` |
| $\mathcal B$, $B(i)$ | 基变量的下标集；第 $i$ 行的基变量 | {doc}`03-03-simplex-method` |
| $\bar\vb=\mB^{-1}\vb$, $\bar a_{ij}$ | 更新后的右端项与单纯形表元素 | {doc}`03-03-simplex-method` |
| $\vpi$, $\pi_i$ | 单纯形乘子、对偶变量、影子价格：$\vpi^\T=\vc_B^\T\mB^{-1}$ | {doc}`01-01-linear-algebra` |
| $\bar c_j$ | 检验数（相对费用）$c_j-\vpi^\T\mA_{\cdot j}$ | {doc}`03-02-lp-geometry` |
| $-z$ 行 | 单纯形表的目标行：各 $\bar c_j$，以及右端的 $-z$ | {doc}`03-03-simplex-method` |
| $s$, $r$ | 单纯形迭代的进基列与主元行 | {doc}`03-03-simplex-method` |
| $\theta$ | 主元变换的步长，即最小比值 $\min\{\bar b_i/\bar a_{is}:\bar a_{is}>0\}$ | {doc}`03-03-simplex-method` |
| $R_i$ | 第 $i$ 行的人工变量 | {doc}`03-04-two-phase-big-m` |
| $\phi$ | 第一阶段的目标函数，即人工变量之和 | {doc}`03-04-two-phase-big-m` |
| $M$ | 大 M 罚系数（一个很大的常数） | {doc}`03-04-two-phase-big-m` |
| (P), (D) | 一对原始–对偶问题中的原问题与对偶问题 | {doc}`03-05-lp-duality` |
| $w$, $w^*$ | 对偶目标函数值 $\vb^\T\vpi$；其最优值 | {doc}`03-05-lp-duality` |
| $I(\vx)$ | 在 $\vx$ 处有效（紧）的约束的下标集 | {doc}`03-05-lp-duality` |
| $\mT$ | 加边逆矩阵 $\begin{pmatrix}\mB&\vzero\\\vc_B^\T&1\end{pmatrix}^{-1}$（$T$ 矩阵） | {doc}`01-01-linear-algebra` |
| $\mB_k$, $\mE_k$ | 修正单纯形法第 $k$ 次迭代的基矩阵与 eta 矩阵 | {doc}`03-06-revised-dual-simplex` |
| $\theta_D$ | 对偶单纯形法的对偶比值 | {doc}`03-06-revised-dual-simplex` |
| $\vl$, $\vu$ | 变量的下界向量与上界向量 | {doc}`04-04-minimum-cost-flow` |

## 网络

| 符号 | 含义 | 首次出现 |
|---|---|---|
| $G=(\mathcal N,\mathcal A)$ | 节点集为 $\mathcal N$、弧（边）集为 $\mathcal A$ 的网络 | {doc}`04-01-networks-mst` |
| $(i,j)$ | 由弧尾 $i$ 指向弧头 $j$ 的弧；在无向网络中为边 $i$–$j$ | {doc}`04-01-networks-mst` |
| $\abs{\mathcal N}$, $\abs{\mathcal A}$ | 节点数，弧数 | {doc}`04-01-networks-mst` |
| $\mathcal A(i)$ | 从节点 $i$ 出发的弧（邻接表） | {doc}`04-01-networks-mst` |
| $\deg(i)$ | 节点 $i$ 的度 | {doc}`04-01-networks-mst` |
| $c_{ij}$ | 弧 $(i,j)$ 的长度或单位费用 | {doc}`04-01-networks-mst` |
| $\mathcal T$ | 一棵（生成）树，视为边的集合 | {doc}`04-01-networks-mst` |
| $\mathcal W$, $c(\mathcal W)$ | 一个圈及其长度或费用 | {doc}`04-01-networks-mst` |
| $[S,\bar S]$ | $S$ 与 $\bar S=\mathcal N\setminus S$ 之间的割 | {doc}`04-01-networks-mst` |
| $[\gamma(j),P(j)]$ | Prim 算法的标号：连向树的最便宜边的费用，及该边在树中的端点 | {doc}`04-01-networks-mst` |
| $\mA$ | 节点–弧关联矩阵；弧 $(i,j)$ 对应的列为 $\ve_i-\ve_j$ | {doc}`04-01-networks-mst` |
| $s$, $t$ | 起点与终点（最短路）；源点与汇点（网络流） | {doc}`04-02-shortest-paths` |
| $V(j)$, $P(j)$ | 最短路标号：从 $s$ 到 $j$ 的最短路长度，以及 $j$ 在该路上的前驱节点 | {doc}`04-02-shortest-paths` |
| $S$, $T$ | 具有永久标号和临时标号的节点集合（Dijkstra 算法） | {doc}`04-02-shortest-paths` |
| $\mD_k=(d^{(k)}_{ij})$, $\mP_k=(p^{(k)}_{ij})$ | 第 $k$ 次迭代后的 Floyd–Warshall 距离矩阵与前驱矩阵 | {doc}`04-02-shortest-paths` |
| $\bar c_{ij}$ | 弧 $(i,j)$ 的约化长度或检验数：用距离标号时为 $c_{ij}+V(i)-V(j)$，用节点势时为 $c_{ij}-\pi_i+\pi_j$ | {doc}`04-02-shortest-paths` |
| $\pi_i$ | 节点势，即节点 $i$ 的平衡方程所对应的对偶变量 | {doc}`04-02-shortest-paths` |
| $x_{ij}$ | 弧 $(i,j)$ 上的流量 | {doc}`04-03-maximum-flow` |
| $u_{ij}$ | 弧 $(i,j)$ 的容量（上界） | {doc}`04-03-maximum-flow` |
| $l_{ij}$ | 弧 $(i,j)$ 上流量的下界 | {doc}`04-04-minimum-cost-flow` |
| $v$, $v^*$ | 流值；最大流值 | {doc}`04-03-maximum-flow` |
| $r_{ij}$, $G(\vx)$ | 残余容量；流 $\vx$ 的残量网络 | {doc}`04-03-maximum-flow` |
| $\mathcal P$, $\Delta$ | 增广路；沿该路（或沿某个圈）输送的流量 | {doc}`04-03-maximum-flow` |
| $u[S,\bar S]$ | $s$–$t$ 割的容量 | {doc}`04-03-maximum-flow` |
| $b_i$ | 节点 $i$ 的净供应量（$b_i>0$ 为供应，$b_i<0$ 为需求） | {doc}`04-04-minimum-cost-flow` |

## 动态规划

| 符号 | 含义 | 首次出现 |
|---|---|---|
| $k$, $N$ | 阶段下标；阶段数 | {doc}`05-01-dp-fundamentals` |
| $y$, $\mathcal Y_k$ | 状态；第 $k$ 阶段的状态集合 | {doc}`05-01-dp-fundamentals` |
| $d$, $\mathcal D_k(y)$ | 决策；第 $k$ 阶段处于状态 $y$ 时允许的决策 | {doc}`05-01-dp-fundamentals` |
| $t_k(y,d)$ | 状态转移函数：下一个状态 | {doc}`05-01-dp-fundamentals` |
| $c_k(y,d)$, $g(y)$ | 阶段费用；终端费用 | {doc}`05-01-dp-fundamentals` |
| $V_k(y)$ | 最优值函数：从第 $k$ 阶段的状态 $y$ 到终点的最优值 | {doc}`05-01-dp-fundamentals` |
| $P_k(y)$ | 最优策略函数：第 $k$ 阶段处于状态 $y$ 时的最优决策 | {doc}`05-01-dp-fundamentals` |
| $F_k(y)$ | 前向值：从起点到第 $k$ 阶段状态 $y$ 的最优值 | {doc}`05-01-dp-fundamentals` |
| $(k,y)$, $c^{\uparrow}(k,y)$, $c^{\downarrow}(k,y)$ | 格点网络的顶点（列，层）；其向上弧与向下弧的费用 | {doc}`05-01-dp-fundamentals` |
| $\mathrm U$, $\mathrm D$ | 格点网络上的决策“向上”与“向下” | {doc}`05-01-dp-fundamentals` |
| $\mathrm B$, $\mathrm K$, $\mathrm O$, $\mathrm L$ | 决策：购买、保留、大修、租赁（设备更新） | {doc}`05-02-dp-equipment-replacement` |
| $i$; $c(i)$, $p$, $t(i)$, $s(i)$ | 机器役龄；役龄为 $i$ 时的运行费用，新机器的价格，旧设备折价，残值 | {doc}`05-02-dp-equipment-replacement` |
| $U_k$ | 单参数设备更新模型的值函数 | {doc}`05-02-dp-equipment-replacement` |
| $x$, $X$ | 分配给某项活动的资源单位数；资源总单位数 | {doc}`05-03-dp-resource-allocation` |
| $r_k(x)$ | 活动 $k$ 获得 $x$ 单位资源时的收益 | {doc}`05-03-dp-resource-allocation` |
| $V^{(k)}(j)$, $P^{(k)}(j)$ | 至多使用 $k$ 条弧时的 Bellman–Ford 值与前驱节点（过程 1） | {doc}`05-04-dp-shortest-paths` |
| $\tilde V^{(k)}(j)$, $\hat V^{(k)}(j)$ | 过程 2 与过程 3 的值：至多含 $k-1$ 次下降、至多含 $k-1$ 次反向的最短途径 | {doc}`05-04-dp-shortest-paths` |
| $D^{(k)}(i,j)$, $P^{(k)}(i,j)$ | Floyd 值：中间节点取自 $\{1,\dots,k\}$ 时 $i$–$j$ 的最短距离；$j$ 的前驱节点 | {doc}`05-04-dp-shortest-paths` |
| $\tilde D^{(k)}(i,j)$ | Dantzig 值：起点与终点均取自 $\{1,\dots,k\}$ 时的同一量 | {doc}`05-04-dp-shortest-paths` |
| $\mathcal S_k$, $i_k$, $T^{(k)}(j)$ | 动态规划形式的 Dijkstra 方法：距节点 1 最近的 $k$ 个节点，第 $k$ 阶段加入的节点，$j$ 的临时标号 | {doc}`05-04-dp-shortest-paths` |
| $\bar V(i)$ | 从 $i$ 到最后一个节点的最短路长度（逆序递推） | {doc}`05-04-dp-shortest-paths` |

## 非线性规划

| 符号 | 含义 | 首次出现 |
|---|---|---|
| $f^*$ | $f$ 在局部或全局极小点 $\vx^*$ 处的值 | {doc}`06-01-nlp-models` |
| $\vd$, $\alpha$ | 方向；沿该方向的步长 | {doc}`06-02-optimality-conditions` |
| $D(\vx)$ | $\Omega$ 在 $\vx$ 处的可行方向锥 | {doc}`06-02-optimality-conditions` |
| $L(\vx^*)$, $T(\vx^*)$ | 线性化可行方向锥；$\Omega$ 在 $\vx^*$ 处的切锥 | {doc}`06-03-kkt-conditions` |
| $\nabla_{\vx}L$, $\nabla^2_{\vx\vx}L$ | 拉格朗日函数关于 $\vx$ 的梯度与黑塞矩阵 | {doc}`06-03-kkt-conditions` |
| $C(\vx^*,\vlambda^*)$ | 临界锥 | {doc}`06-03-kkt-conditions` |
| $\lambda_0$ | Fritz John 条件中 $f$ 的乘子 | {doc}`06-03-kkt-conditions` |
| $\vd^{(k)}$, $\alpha_k$ | 第 $k$ 次迭代的搜索方向与步长：$\vx^{(k+1)}=\vx^{(k)}+\alpha_k\vd^{(k)}$ | {doc}`06-04-line-search` |
| $[a_k,b_k]$ | 一维搜索的不确定区间 | {doc}`06-04-line-search` |
| $\tau$ | 黄金分割比 $(\sqrt5-1)/2\approx0.618$ | {doc}`06-04-line-search` |
| $F_N$ | 斐波那契数 | {doc}`06-04-line-search` |
| $c_1$, $c_2$ | Armijo 常数与曲率（Wolfe）常数，$0<c_1<c_2<1$ | {doc}`06-04-line-search` |
| $\gamma$ | 回溯的缩减因子 | {doc}`06-04-line-search` |
| $\vg^{(k)}$ | 第 $k$ 个迭代点处的梯度 $\nabla f(\vx^{(k)})$ | {doc}`06-05-unconstrained-algorithms` |
| $\mH(\vx)$, $\bar\mH_k$ | 黑塞矩阵 $\nabla^2 f(\vx)$；它的一个正定修正 | {doc}`06-05-unconstrained-algorithms` |
| $\mD_k$ | 正定尺度矩阵，$\vd^{(k)}=-\mD_k\vg^{(k)}$（拟牛顿法中为黑塞矩阵之逆的近似） | {doc}`06-05-unconstrained-algorithms` |
| $\vs^{(k)}$, $\vy^{(k)}$ | 拟牛顿修正中的步长向量与梯度变化量 | {doc}`06-05-unconstrained-algorithms` |
| $\Gamma(\vx)$ | Frank–Wolfe 间隙 | {doc}`06-06-constrained-algorithms` |
| $\psi(\vx)$, $\rho$ | 罚函数与罚参数 | {doc}`06-06-constrained-algorithms` |
| $\Phi_\rho(\vx)$, $\vx(\rho)$ | 带罚目标函数 $f+\rho\psi$ 及其极小点 | {doc}`06-06-constrained-algorithms` |
| $L_A(\vx,\vmu;\rho)$ | 增广拉格朗日函数 | {doc}`06-06-constrained-algorithms` |

## 线性化

| 符号 | 含义 | 首次出现 |
|---|---|---|
| $\delta_i\in\{0,1\}$, $\vdelta$ | 0–1 变量（指示变量）及其向量 | {doc}`07-01-logical-constraints` |
| $M$, $M_i$ | 大 M 常数：$g$、$g_i$ 在 $\Omega_0$ 上的上界 | {doc}`07-01-logical-constraints` |
| $\underline g$, $\underline g_i$ | $g$、$g_i$ 在 $\Omega_0$ 上的下界 | {doc}`07-01-logical-constraints` |
| $\Omega_0$ | 模型中所有其他约束所允许的集合 | {doc}`07-01-logical-constraints` |
| $F$, $F_i$ | 固定费用（准备费用、开设费用） | {doc}`07-01-logical-constraints` |
| $\omega$, $\omega_{ij}$ | 替代乘积（$\delta_A\delta_B$, $\delta x$, $xy$）的变量 | {doc}`07-01-logical-constraints` |
| $\underline x$, $\overline x$ | 标量变量的下界与上界 | {doc}`07-02-linearization-techniques` |
| $\vx^+$, $\vx^-$ | $\vx$ 的正部与负部 | {doc}`07-02-linearization-techniques` |
| $\eta$ | 用于界定绝对值、最大值或最小值的辅助变量 | {doc}`07-02-linearization-techniques` |
| $f^a$ | $f$ 的分段线性逼近 | {doc}`07-03-piecewise-linear-approximation` |
| $a_1<\dots<a_n$ | 分段线性逼近的分段点 | {doc}`07-03-piecewise-linear-approximation` |
| $\alpha_i$ | 分段点的插值权重（SOS2 权重） | {doc}`07-03-piecewise-linear-approximation` |
| $\xi_k$ | 区段 $k$ 被填充的比例（增量形式） | {doc}`07-03-piecewise-linear-approximation` |

## 一符多义的符号

少数字母被有意地用于两种不同的角色，而这两种角色从不在同一论证中同时出现。根据上下文以及符号的自变量，总能将它们区分开来。

| 符号 | 含义 |
|---|---|
| $L$ | 拉格朗日函数 $L(\vx,\vlambda,\vmu)$（{doc}`01-02-matrix-calculus`、{doc}`02-03-convex-optimization-duality`）与只有一个自变量的线性化可行方向锥 $L(\vx^*)$（{doc}`06-03-kkt-conditions`）；粗体 $\mL$ 是三角因子（{doc}`01-01-linear-algebra`）或图的拉普拉斯矩阵（{doc}`09-04-code-networks`） |
| $\theta$ | 单纯形主元变换的最小比值（{doc}`03-03-simplex-method`）与凸组合的权重（{doc}`02-01-convex-sets`、{doc}`02-02-convex-functions`） |
| $s$, $t$ | 源点与汇点（{doc}`04-02-shortest-paths`–{doc}`04-04-minimum-cost-flow`）；$s$ 也表示主元变换的进基列（{doc}`03-03-simplex-method`），$s_i$ 表示松弛变量（{doc}`03-02-lp-geometry`）；$t_k(y,d)$ 是动态规划的状态转移（{doc}`05-01-dp-fundamentals`），$t(i)$ 是旧设备折价（{doc}`05-02-dp-equipment-replacement`） |
| $S$, $T$ | Dijkstra 算法的永久标号集与临时标号集（{doc}`04-02-shortest-paths`）；$S$ 也表示割的一侧（{doc}`04-01-networks-mst`、{doc}`04-03-maximum-flow`），$T^{(k)}(j)$ 表示动态规划形式的 Dijkstra 方法中的临时标号（{doc}`05-04-dp-shortest-paths`），$T(\vx^*)$ 表示切锥（{doc}`06-03-kkt-conditions`），粗体 $\mT$ 表示加边逆矩阵（{doc}`03-06-revised-dual-simplex`） |
| $P$ | 前驱节点 $P(j)$（{doc}`04-02-shortest-paths`）与最优策略 $P_k(y)$（{doc}`05-01-dp-fundamentals`），是同一思想在两种情形中的体现；{doc}`03-02-lp-geometry`中的多面体 $P$；粗体 $\mP$ 表示投影矩阵（{doc}`01-01-linear-algebra`）或 Floyd–Warshall 前驱矩阵 $\mP_k$（{doc}`04-02-shortest-paths`） |
| $\mD_k$ | Floyd–Warshall 距离矩阵（{doc}`04-02-shortest-paths`）与下降法的尺度矩阵（{doc}`06-05-unconstrained-algorithms`） |
| $D$ | 雅可比算子 $D\vg$（{doc}`01-02-matrix-calculus`）、顺序主子式 $D_k$（{doc}`01-01-linear-algebra`）与可行方向锥 $D(\vx)$（{doc}`06-02-optimality-conditions`） |
| $\vy^{(k)}$ | 拟牛顿修正中的梯度变化量（{doc}`06-05-unconstrained-algorithms`）与 Frank–Wolfe 方法线性子问题返回的顶点（{doc}`06-06-constrained-algorithms`） |
| $\alpha$ | 步长（{doc}`06-02-optimality-conditions`、{doc}`06-04-line-search`），以及带下标时的 $\alpha_i$，即{doc}`07-03-piecewise-linear-approximation`中的 SOS2 权重 |
| $\eta$ | 极大化极小模型的保证收益（{doc}`03-01-lp-model-building`）、障碍参数（{doc}`06-06-constrained-algorithms`）与起界定作用的辅助变量（{doc}`07-02-linearization-techniques`）；粗体 $\veta$ 表示 eta 列（{doc}`01-01-linear-algebra`） |
| $\omega$ | 离开格点网络顶点的方向（{doc}`05-01-dp-fundamentals`）与乘积变量（{doc}`07-01-logical-constraints`、{doc}`07-02-linearization-techniques`） |
| $\gamma$ | Prim 标号 $\gamma(j)$（{doc}`04-01-networks-mst`）、回溯因子（{doc}`06-04-line-search`）与罚参数的增长因子（{doc}`06-06-constrained-algorithms`） |
| $\delta$ | 0–1 变量 $\delta_i$（{doc}`07-01-logical-constraints`），以及不带下标时，灵敏度分析中费用或右端项的改变量（{doc}`03-07-sensitivity-analysis`） |
| $\kappa$ | 条件数 $\kappa(\mA)$（{doc}`01-01-linear-algebra`、{doc}`06-05-unconstrained-algorithms`）、百分之百法则中的比例 $\kappa_r$（{doc}`03-07-sensitivity-analysis`）与曲率界 $\kappa_k$（{doc}`07-03-piecewise-linear-approximation`） |
| $c_1$, $c_2$ | Armijo 常数与曲率常数（{doc}`06-04-line-search`），以及两种资源的单价（{doc}`05-03-dp-resource-allocation`） |
