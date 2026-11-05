"""Figures of chapter 6.4, Line Search Methods.

All figures are computed from the running example f(x) = sin x on [3, 6]
(minimizer x* = 3 pi / 2):

* ``06-04-interval-reduction``: intervals of bisection (eps = 0.05) and of
  golden-section search (eps = 0.05), examples ex-06-04-bisection-sine and
  ex-06-04-golden-sine;
* ``06-04-golden-section-geometry``: placement of the golden-section trial
  points and the reuse of one point (section sec-06-04-golden-section);
* ``06-04-newton-tangents``: Newton's method as tangent steps on f' = cos x
  from x0 = 4 and x0 = 3 (example ex-06-04-newton-sine);
* ``06-04-armijo-wolfe``: phi(alpha) = sin(3 + 0.99 alpha), the Armijo line,
  the Wolfe intervals and two backtracking trials (section sec-06-04-inexact).
"""

from math import cos, pi, sin, sqrt

import matplotlib.pyplot as plt
import numpy as np

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl

TAU = (sqrt(5) - 1) / 2  # golden-section ratio tau
XSTAR = 3 * pi / 2  # minimizer of sin x on [3, 6]


def _t(lang, en, zh):
    """Label in the language of the edition: ``zh`` for "zh", else ``en``."""
    return zh if lang == "zh" else en


def _bisection_rows(a=3.0, b=6.0, eps=0.05):
    """Intervals [a_k, b_k] and midpoints x_k of bisection on f' = cos x."""
    rows = []
    while b - a > eps:
        x = (a + b) / 2
        rows.append((a, b, x))
        if cos(x) > 0:  # f increasing at x: minimizer to the left
            b = x
        else:
            a = x
    rows.append((a, b, None))  # stopping row (no new midpoint)
    return rows


def _golden_rows(a=3.0, b=6.0, eps=0.05):
    """Intervals and trial points of golden-section search on f = sin x.

    Each row is (a, b, p, q, reused) where ``reused`` is "p" or "q" (the
    point carried over from the previous row) or None in row 0.
    """
    p, q = a + (1 - TAU) * (b - a), a + TAU * (b - a)
    rows, reused = [], None
    while b - a > eps:
        rows.append((a, b, p, q, reused))
        if sin(p) > sin(q):  # discard [a, p]; old q becomes new p
            a, p = p, q
            q = a + TAU * (b - a)
            reused = "p"
        else:  # discard (q, b]; old p becomes new q
            b, q = q, p
            p = a + (1 - TAU) * (b - a)
            reused = "q"
    rows.append((a, b, None, None, None))  # stopping row
    return rows


def interval_reduction(outdir, lang="en"):
    """Nested intervals of bisection (left) and golden-section search (right).

    Used in chapter 6.4 (fig-06-04-interval-reduction). Row k shows the
    interval [a_k, b_k] as a horizontal bar at height -k; dots mark the evaluated
    points; the dashed vertical line is x* = 3 pi / 2.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.4), sharey=False)
    panels = [
        (axes[0], _bisection_rows(), _t(lang, "(a) bisection", "(a) 二分法")),
        (
            axes[1],
            _golden_rows(),
            _t(lang, "(b) golden-section search", "(b) 黄金分割法"),
        ),
    ]
    for ax, rows, title in panels:
        for k, row in enumerate(rows):
            a, b = row[0], row[1]
            # the interval of uncertainty as a thick bar with end ticks
            ax.plot([a, b], [-k, -k], color=BLUE, lw=3, solid_capstyle="butt")
            ax.plot([a, a], [-k - 0.18, -k + 0.18], color=BLUE, lw=1.2)
            ax.plot([b, b], [-k - 0.18, -k + 0.18], color=BLUE, lw=1.2)
            if len(row) == 3 and row[2] is not None:
                # bisection: the midpoint where f' is evaluated
                ax.plot(row[2], -k, "o", color=ORANGE, ms=6, zorder=5)
            if len(row) == 5 and row[2] is not None:
                _, _, p, q, reused = row
                for name, x in (("p", p), ("q", q)):
                    filled = reused == name  # reused point: filled dot
                    ax.plot(
                        x, -k, "o", ms=6, zorder=5,
                        color=ORANGE,
                        markerfacecolor=ORANGE if filled else "none",
                        markeredgewidth=1.5,
                    )  # fmt: skip
        n = len(rows)
        # x* = 3 pi / 2
        ax.axvline(XSTAR, color=MUTED, ls="--", lw=1)
        ax.text(XSTAR, 0.9, r"$x^*=3\pi/2$", ha="center", va="bottom",
                color=MUTED, fontsize=10)  # fmt: skip
        ax.set_xlim(2.9, 6.1)
        ax.set_ylim(-n + 0.4, 1.6)
        ax.set_yticks([-k for k in range(n)])
        ax.set_yticklabels([str(k) for k in range(n)])
        ax.set_ylabel("$k$")
        ax.set_xlabel("$x$")
        ax.set_title(title, fontsize=11, loc="left")
        ax.grid(axis="y", visible=False)
    fig.tight_layout()
    return save_mpl(fig, "06-04-interval-reduction", outdir)


def golden_section_geometry(outdir, lang="en"):
    """Schematic of the golden-section trial points on two levels.

    Used in chapter 6.4 (fig-06-04-golden-section-geometry). The upper line is
    [a_k, b_k] = [0, 1] with p_k = 1 - tau and q_k = tau; the lower line is
    the kept interval [a_k, q_k] in which the old p_k becomes q_{k+1}.
    """
    fig, ax = plt.subplots(figsize=(6.4, 2.9))
    ax.set_axis_off()
    t2 = 1 - TAU  # = tau^2
    y1, y2 = 1.0, 0.0  # heights of the two levels

    def point(x, y, label, filled=True, color=ORANGE, above=True):
        # a trial point with its label above or below the line
        ax.plot(x, y, "o", ms=8, color=color,
                markerfacecolor=color if filled else "white",
                markeredgewidth=1.6, zorder=5)  # fmt: skip
        ax.text(x, y + (0.13 if above else -0.15), label, ha="center",
                va="bottom" if above else "top", fontsize=12)  # fmt: skip

    def dim(x0, x1, y, label):
        # dimension arrow between x0 and x1 at height y, label above it
        ax.annotate("", xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1,
                                    shrinkA=0, shrinkB=0))  # fmt: skip
        ax.text((x0 + x1) / 2, y + 0.03, label, ha="center", va="bottom",
                color=MUTED, fontsize=10.5)  # fmt: skip

    # upper level: [a_k, b_k] with p_k and q_k
    ax.plot([0, 1], [y1, y1], color=INK, lw=1.6)
    for x in (0, 1):
        ax.plot([x, x], [y1 - 0.06, y1 + 0.06], color=INK, lw=1.6)
    ax.text(0, y1 - 0.12, "$a_k$", ha="center", va="top", fontsize=12)
    ax.text(1, y1 - 0.12, "$b_k$", ha="center", va="top", fontsize=12)
    point(t2, y1, "$p_k$")
    point(TAU, y1, "$q_k$")
    # discarded part (q_k, b_k] drawn muted and hatched-looking
    ax.plot([TAU, 1], [y1, y1], color=MUTED, lw=5, alpha=0.35,
            solid_capstyle="butt")  # fmt: skip
    dim(0, t2, y1 + 0.42, r"$(1-\tau)(b_k-a_k)$")
    dim(TAU, 1, y1 + 0.42, r"$(1-\tau)(b_k-a_k)$")
    dim(0, TAU, y1 + 0.72, r"$\tau(b_k-a_k)$")

    # lower level: kept interval [a_{k+1}, b_{k+1}] = [a_k, q_k]
    ax.plot([0, TAU], [y2, y2], color=INK, lw=1.6)
    for x in (0, TAU):
        ax.plot([x, x], [y2 - 0.06, y2 + 0.06], color=INK, lw=1.6)
    ax.text(-0.01, y2 - 0.12, "$a_{k+1}$", ha="right", va="top", fontsize=12)
    ax.text(
        TAU, y2 - 0.12, "$b_{k+1}=q_k$", ha="center", va="top", fontsize=12
    )
    new = _t(lang, "$p_{k+1}$ (new)", "$p_{k+1}$（新）")
    point(TAU * t2, y2, new, filled=False, above=False)
    point(t2, y2, r"$q_{k+1}=p_k$")
    # the reused point: dashed connector from p_k down to q_{k+1}
    ax.plot([t2, t2], [y1 - 0.1, y2 + 0.35], color=MUTED, ls=":", lw=1.2)
    keep = _t(lang, "keep $[a_k,q_k]$", "保留 $[a_k,q_k]$")
    ax.text(1.0, y2 + 0.05, r"$f(p_k)<f(q_k)$" + "\n" + keep,
            ha="right", va="center", fontsize=10.5, color=MUTED)  # fmt: skip
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.55, 2.05)
    return save_mpl(fig, "06-04-golden-section-geometry", outdir)


def newton_tangents(outdir, lang="en"):
    """Newton's method as tangent steps on f'(x) = cos x, two starts.

    Used in chapter 6.4 (fig-06-04-newton-tangents). Panel (a): x0 = 4,
    iterates 4 -> 4.864 -> 4.711 (converges to 3 pi / 2). Panel (b): x0 = 3, the
    almost flat tangent sends x1 to -4.015, then the iterates converge to
    -3 pi / 2, a maximizer of sin x.
    """
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))

    def newton_path(x0, n):
        # iterates x_{k+1} = x_k - cos(x_k) / (-sin(x_k))
        xs = [x0]
        for _ in range(n):
            x = xs[-1]
            xs.append(x - cos(x) / (-sin(x)))
        return xs

    for ax, x0, xlim, title, nsteps in (
        (axes[0], 4.0, (3.6, 5.6), r"(a) $x_0=4$", 2),
        (axes[1], 3.0, (-5.6, 3.6), r"(b) $x_0=3$", 1),
    ):
        xx = np.linspace(*xlim, 400)
        ax.plot(
            xx, np.cos(xx), color=BLUE, lw=2, label="$f$\u2032$(x)=\\cos x$"
        )
        ax.axhline(0, color=INK, lw=0.8)
        xs = newton_path(x0, nsteps)
        for k in range(nsteps):
            xk, xn = xs[k], xs[k + 1]
            # tangent to f' at x_k reaches zero at x_{k+1}
            ax.plot([xk, xn], [cos(xk), 0.0], color=ORANGE, ls="--", lw=1.4)
            ax.plot([xk, xk], [0, cos(xk)], color=MUTED, ls=":", lw=1)
            ax.plot(xk, cos(xk), "o", color=ORANGE, ms=5, zorder=5)
        for xk in xs:
            ax.plot(xk, 0, "o", color=INK, ms=4, zorder=5)
        # labels of the iterates
        if x0 == 4.0:
            ax.annotate("$x_0$", (4.0, 0), xytext=(0, 8), textcoords="offset points",
                        ha="center")  # fmt: skip
            ax.annotate("$x_1$", (xs[1], 0), xytext=(4, -16), textcoords="offset points",
                        ha="left")  # fmt: skip
            ax.annotate("$x_2$", (xs[2], 0), xytext=(-4, -16), textcoords="offset points",
                        ha="right")  # fmt: skip
            ax.axvline(XSTAR, color=GREEN, lw=1, ls="-.")
            ax.text(
                XSTAR, 0.85, r"$3\pi/2$", color=GREEN, ha="right", fontsize=10
            )
        else:
            ax.annotate("$x_0=3$", (3.0, 0), xytext=(0, 8), textcoords="offset points",
                        ha="center")  # fmt: skip
            ax.annotate("$x_1=-4.015$", (xs[1], 0), xytext=(4, 8),
                        textcoords="offset points", ha="left")  # fmt: skip
            ax.axvline(-XSTAR, color=GREEN, lw=1, ls="-.")
            # label at the top, where cos x is far from the line
            ax.text(
                -XSTAR + 0.15,
                0.85,
                r"$-3\pi/2$",
                color=GREEN,
                ha="left",
                fontsize=10,
            )
            # the interval [3, 6] that the iterates should stay in
            ax.axvspan(3.0, 3.6, color=BLUE, alpha=0.08)
        ax.set_xlim(*xlim)
        ax.set_ylim(-1.15, 1.15)
        ax.set_xlabel("$x$")
        ax.set_title(title, fontsize=11, loc="left")
    axes[0].legend(loc="lower right", fontsize=9, frameon=False)
    fig.tight_layout()
    return save_mpl(fig, "06-04-newton-tangents", outdir)


def armijo_wolfe(outdir, lang="en"):
    """The Armijo and Wolfe conditions for phi(alpha) = sin(3 + 0.99 alpha).

    Used in chapter 6.4 (fig-06-04-armijo-wolfe). Data: x = 3, d = -cos 3,
    c1 = 0.25, c2 = 0.5. Intervals of acceptable steps are found on a fine
    grid (Armijo: (0, 2.749]; Wolfe: [1.207, 2.749]; strong Wolfe:
    [1.207, 2.253]); backtracking trials alpha = 4 (rejected), 2 (accepted).
    """
    x0 = 3.0
    d = -cos(x0)
    c1, c2 = 0.25, 0.5
    A = np.linspace(0, 4.5, 4501)
    phi = np.sin(x0 + A * d)
    dphi = np.cos(x0 + A * d) * d
    phi0, dphi0 = sin(x0), cos(x0) * d

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(A, phi, color=BLUE, lw=2, label=r"$\varphi(\alpha)$")
    # tangent at 0 and the flatter Armijo line
    ax.plot(A[A <= 1.25], phi0 + dphi0 * A[A <= 1.25], color=MUTED, lw=1,
            label=_t(lang, "tangent", "切线")
            + " $\\varphi(0)+\\alpha\\varphi$\u2032$(0)$")  # fmt: skip
    ax.plot(A, phi0 + c1 * dphi0 * A, color=ORANGE, ls="--", lw=1.5,
            label=_t(lang, r"Armijo line, $c_1=0.25$",
                     r"Armijo 直线，$c_1=0.25$"))  # fmt: skip

    def spans(mask):
        # maximal runs of True in mask, as (start, end) alpha values
        out, start = [], None
        for i, m in enumerate(mask):
            if m and start is None:
                start = A[i]
            if (not m or i == len(mask) - 1) and start is not None:
                out.append((start, A[i - 1] if not m else A[i]))
                start = None
        return out

    armijo = phi <= phi0 + c1 * A * dphi0
    curv = dphi >= c2 * dphi0
    strong = np.abs(dphi) <= c2 * abs(dphi0)
    bars = [
        (armijo & (A > 0), "Armijo", -1.25, ORANGE),
        (armijo & curv, "Wolfe", -1.4, GREEN),
        (armijo & strong, _t(lang, "strong Wolfe", "强 Wolfe"), -1.55, BLUE),
    ]
    for mask, name, y, col in bars:
        for s, e in spans(mask):
            ax.plot([s, e], [y, y], color=col, lw=5, solid_capstyle="butt")
        ax.text(4.55, y, name, va="center", ha="left", fontsize=9.5, color=col)
    # backtracking trials
    for al, ok in ((4.0, False), (2.0, True)):
        v = sin(x0 + al * d)
        ax.plot(al, v, "o", ms=8, color=INK,
                markerfacecolor=INK if ok else "white", zorder=6)  # fmt: skip
        ax.annotate(
            r"$\alpha=%g$ %s"
            % (
                al,
                (
                    _t(lang, "accepted", "接受")
                    if ok
                    else _t(lang, "rejected", "拒绝")
                ),
            ),
            (al, v),
            # accepted: below right of the point, clear of the curve
            xytext=(8, -14 if ok else -2),
            textcoords="offset points",
            fontsize=9.5,
        )
    # exact minimizer alpha* = (3 pi / 2 - 3) / d
    ast = (XSTAR - x0) / d
    ax.plot([ast, ast], [-1.0, -1.06], color=INK, lw=1)
    ax.text(ast, -1.065, r"$\alpha^*$", ha="center", va="top", fontsize=10)
    ax.set_xlim(0, 4.5)
    ax.set_ylim(-1.65, 1.25)
    ax.set_xlabel(r"$\alpha$")
    ax.legend(loc="upper left", fontsize=9, frameon=False)
    fig.tight_layout()
    return save_mpl(fig, "06-04-armijo-wolfe", outdir)


FIGURES = {
    "06-04-interval-reduction": interval_reduction,
    "06-04-golden-section-geometry": golden_section_geometry,
    "06-04-newton-tangents": newton_tangents,
    "06-04-armijo-wolfe": armijo_wolfe,
}
