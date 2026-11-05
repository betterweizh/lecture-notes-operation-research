"""Figures of chapter 5.3, Resource Allocation and Capital Budgeting.

* ``05-03-return-functions``: the return tables r_k(x), k = 1..4, of the
  one-resource allocation example, with the optimal allocation (2, 6, 0, 0)
  marked.

Labels are symbols and numbers (language neutral).
"""

import matplotlib.pyplot as plt

from .common import INK, PALETTE, save_mpl

#: returns r_k(x) for x = 0..8 (Example 5.3.1)
RETURNS = [
    [0, 3, 7, 10, 12, 13, 14, 14, 14],
    [0, 1, 2, 4, 8, 13, 17, 19, 20],
    [0, 2, 4, 6, 8, 10, 12, 14, 16],
    [0, 1, 3, 6, 9, 12, 14, 16, 17],
]
#: optimal allocation found by the DP
X_OPT = [2, 6, 0, 0]
MARKERS = ["o", "s", "^", "D"]


def return_functions(outdir, lang="en"):
    """Return functions of the four activities (fig-05-03-returns).

    Each table r_k(x) is drawn as points joined by thin lines (only integer
    allocations exist; the lines guide the eye). The optimal allocation
    x = (2, 6, 0, 0) of Example 5.3.1 is marked by rings on
    r_1(2) = 7 and r_2(6) = 17 (the zero allocations sit at the origin).
    """
    fig, ax = plt.subplots(figsize=(5.5, 4))
    xs = range(9)
    for k, (r, mk) in enumerate(zip(RETURNS, MARKERS)):
        # one color and marker per activity
        ax.plot(
            xs,
            r,
            marker=mk,
            markersize=5,
            linewidth=1.3,
            color=PALETTE[k],
            label=f"$r_{k + 1}(x)$",
        )
    # rings on the chosen points with a positive allocation
    for k, x in enumerate(X_OPT):
        if x > 0:
            ax.plot(
                x,
                RETURNS[k][x],
                marker="o",
                markersize=15,
                markerfacecolor="none",
                markeredgecolor=INK,
                markeredgewidth=1.4,
                zorder=5,
            )
            ax.annotate(
                f"$r_{k + 1}({x}) = {RETURNS[k][x]}$",
                xy=(x, RETURNS[k][x]),
                xytext=(-14, 8),
                textcoords="offset points",
                ha="right",
                fontsize=10,
            )
    ax.set_xlabel("$x$")
    ax.set_ylabel("$r_k(x)$")
    ax.set_xticks(list(xs))
    ax.set_yticks(range(0, 21, 5))
    ax.set_xlim(-0.3, 8.4)
    ax.set_ylim(-0.5, 21.5)
    ax.grid(True)
    ax.legend(loc="upper left", frameon=False)
    return save_mpl(fig, "05-03-return-functions", outdir)


FIGURES = {
    "05-03-return-functions": return_functions,
}
