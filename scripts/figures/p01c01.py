"""Figures of chapter 1.1, Linear Algebra Essentials.

Five matplotlib figures, all drawn from closed-form data (no randomness):

* ``01-01-span-combination``: a linear combination of two vectors in R^2 and the
  span of a dependent pair (section on linear combinations and span);
* ``01-01-projection-line``: orthogonal projection of b onto the line through a
  (example ex-01-01-projection-line);
* ``01-01-circle-ellipse``: the symmetric matrix [[2, 1], [1, 2]] maps the unit
  circle to an ellipse whose axes are its eigenvectors (ex-01-01-spectral-2x2);
* ``01-01-quadratic-forms``: contours of x^T A x for a positive definite, a
  singular positive semidefinite and an indefinite matrix;
* ``01-01-unit-balls``: unit balls of the 1-, 2- and infinity-norms in R^2.
"""

import matplotlib.pyplot as plt
import numpy as np

from .common import BLUE, GREEN, INK, MUTED, ORANGE, save_mpl


def _arrow(ax, start, end, color=INK, lw=1.8, ls="-", zorder=4):
    """Draw a vector as an arrow from ``start`` to ``end`` (data coords)."""
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=lw,
            linestyle=ls,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=14,
        ),
        zorder=zorder,
    )


def _plane_axes(ax, xlim, ylim):
    """Equal-aspect axes through the origin with x1, x2 titles."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axhline(0, color=MUTED, lw=0.9, zorder=1)
    ax.axvline(0, color=MUTED, lw=0.9, zorder=1)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$x_2$")


def span_combination(outdir, lang="en"):
    """Linear combination and span in R^2 (two panels).

    (a) u = (2, 0.5) and v = (0.5, 1.5) are independent; w = 1.5u + v is
    built by the parallelogram rule; the thin skewed grid shows the lines
    alpha*u + t*v and beta*v + t*u for integer alpha, beta (coordinates in
    the basis {u, v}), which cover the whole plane: span{u, v} = R^2.
    (b) u = (2, 1) and v = -0.5u are dependent: every combination lies on
    the line through the origin with direction u, so span{u, v} is a line.
    Used in section sec-01-01-span of 01-01-linear-algebra.md.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.2))

    # --- panel (a): independent pair -----------------------------------
    u = np.array([2.0, 0.5])
    v = np.array([0.5, 1.5])
    w = 1.5 * u + v  # = (3.5, 2.25)
    _plane_axes(ax1, (-1.6, 4.6), (-1.6, 3.6))
    ax1.grid(False)  # the skewed grid below replaces the square one
    # skewed coordinate grid: lines alpha*u + t*v and beta*v + t*u for
    # integer alpha, beta; every point of the plane lies on this grid's
    # coordinate system, so span{u, v} = R^2
    tt = np.array([-6.0, 6.0])
    for k in range(-4, 5):
        ax1.plot(
            k * u[0] + tt * v[0],
            k * u[1] + tt * v[1],
            color=MUTED,
            lw=0.5,
            alpha=0.6,
            zorder=1,
        )
        ax1.plot(
            k * v[0] + tt * u[0],
            k * v[1] + tt * u[1],
            color=MUTED,
            lw=0.5,
            alpha=0.6,
            zorder=1,
        )
    # parallelogram: 1.5u then v, dashed helper edges
    ax1.plot(
        [1.5 * u[0], w[0]], [1.5 * u[1], w[1]], ls="--", color=MUTED, lw=1.1
    )
    ax1.plot([v[0], w[0]], [v[1], w[1]], ls="--", color=MUTED, lw=1.1)
    _arrow(ax1, (0, 0), 1.5 * u, color=BLUE, lw=1.2, ls="--")
    _arrow(ax1, (0, 0), u, color=BLUE)
    _arrow(ax1, (0, 0), v, color=GREEN)
    _arrow(ax1, (0, 0), w, color=ORANGE)
    ax1.text(*(u + [-0.15, -0.42]), r"$\mathbf{u}$", color=BLUE, fontsize=13)
    ax1.text(
        *(1.5 * u + [0.05, -0.45]),
        r"$1.5\,\mathbf{u}$",
        color=BLUE,
        fontsize=12
    )
    ax1.text(*(v + [-0.45, 0.05]), r"$\mathbf{v}$", color=GREEN, fontsize=13)
    ax1.text(
        *(w + [-1.55, 0.18]),
        r"$1.5\,\mathbf{u}+\mathbf{v}$",
        color=ORANGE,
        fontsize=12
    )
    ax1.set_title("(a)", loc="left")

    # --- panel (b): dependent pair ---------------------------------------
    u2 = np.array([2.0, 1.0])
    v2 = -0.5 * u2  # dependent: v = -0.5 u
    _plane_axes(ax2, (-2.6, 3.6), (-1.6, 3.6))
    t = np.array([-2.0, 2.0])
    ax2.plot(
        t * u2[0],
        t * u2[1],
        color=BLUE,
        lw=5,
        alpha=0.18,
        solid_capstyle="butt",
    )
    _arrow(ax2, (0, 0), u2, color=BLUE)
    _arrow(ax2, (0, 0), v2, color=GREEN)
    ax2.text(*(u2 + [0.0, -0.5]), r"$\mathbf{u}$", color=BLUE, fontsize=13)
    ax2.text(
        *(v2 + [0.05, -0.62]),
        r"$\mathbf{v}=-0.5\,\mathbf{u}$",
        color=GREEN,
        fontsize=12
    )
    # label of the line, written along it (slope 1/2 = 26.57 degrees)
    ax2.text(
        1.95,
        1.3,
        r"span$\{\mathbf{u},\mathbf{v}\}$",
        color=BLUE,
        fontsize=12,
        rotation=np.degrees(np.arctan(0.5)),
        rotation_mode="anchor",
    )
    ax2.set_title("(b)", loc="left")
    fig.tight_layout()
    return save_mpl(fig, "01-01-span-combination", outdir)


def projection_line(outdir, lang="en"):
    """Orthogonal projection of b = (1, 2) onto the line through a = (2, 1).

    p = (a^T b / a^T a) a = 0.8 a = (1.6, 0.8); the residual r = b - p =
    (-0.6, 1.2) is orthogonal to a (right-angle mark at p). Used in
    example ex-01-01-projection-line of 01-01-linear-algebra.md.
    """
    fig, ax = plt.subplots(figsize=(5.4, 3.9))
    a = np.array([2.0, 1.0])
    b = np.array([1.0, 2.0])
    p = (a @ b) / (a @ a) * a  # (1.6, 0.8)
    e = b - p  # residual r = (-0.6, 1.2)
    _plane_axes(ax, (-0.6, 3.0), (-0.4, 2.4))
    # the line span{a}
    t = np.array([-0.3, 1.45])
    ax.plot(t * a[0], t * a[1], color=MUTED, lw=1.0)
    ax.text(
        2.62,
        1.6,
        r"span$\{\mathbf{a}\}$",
        color=MUTED,
        fontsize=11,
        ha="center",
    )
    _arrow(ax, (0, 0), a, color=BLUE)
    _arrow(ax, (0, 0), b, color=INK)
    _arrow(ax, (0, 0), p, color=ORANGE, lw=2.4)
    _arrow(ax, p, b, color=GREEN, ls="--")
    # right-angle mark at p: small square spanned by unit a and unit e
    ua = a / np.linalg.norm(a) * 0.12
    ue = e / np.linalg.norm(e) * 0.12
    sq = np.array([p + ua, p + ua + ue, p + ue])
    ax.plot(sq[:, 0], sq[:, 1], color=INK, lw=0.9)
    ax.text(*(a + [0.06, -0.2]), r"$\mathbf{a}$", color=BLUE, fontsize=13)
    ax.text(*(b + [-0.12, 0.1]), r"$\mathbf{b}$", color=INK, fontsize=13)
    ax.text(*(p + [0.02, -0.33]), r"$\mathbf{p}$", color=ORANGE, fontsize=13)
    ax.text(
        *((p + b) / 2 + [0.12, -0.17]),
        r"$\mathbf{r}=\mathbf{b}-\mathbf{p}$",
        color=GREEN,
        fontsize=12
    )
    return save_mpl(fig, "01-01-projection-line", outdir)


def circle_ellipse(outdir, lang="en"):
    """Image of the unit circle under A = [[2, 1], [1, 2]].

    The eigenvectors v1 = (1, 1)/sqrt2 (eigenvalue 3) and v2 = (1, -1)/sqrt2
    (eigenvalue 1) are only stretched: A v1 = 3 v1, A v2 = v2, so the image
    is an ellipse with semi-axes 3 and 1 along v1 and v2. A generic vector
    such as x = (1, 0) is also turned: A x = (2, 1). Used in example
    ex-01-01-spectral-2x2 of 01-01-linear-algebra.md.
    """
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    A = np.array([[2.0, 1.0], [1.0, 2.0]])
    s = np.linspace(0, 2 * np.pi, 400)
    circle = np.vstack([np.cos(s), np.sin(s)])
    ellipse = A @ circle
    _plane_axes(ax, (-3.0, 3.0), (-2.6, 2.9))
    ax.plot(*circle, color=MUTED, lw=1.2)
    ax.plot(*ellipse, color=BLUE, lw=1.8)
    v1 = np.array([1.0, 1.0]) / np.sqrt(2)
    v2 = np.array([1.0, -1.0]) / np.sqrt(2)
    # eigenvectors (on the circle) and their images (on the ellipse)
    _arrow(ax, (0, 0), A @ v1, color=BLUE, lw=1.4, ls="--")
    _arrow(ax, (0, 0), v1, color=BLUE)
    _arrow(ax, (0, 0), v2, color=GREEN, lw=2.6)
    # a generic vector x = e1 and its image A x = (2, 1)
    x = np.array([1.0, 0.0])
    _arrow(ax, (0, 0), x, color=ORANGE)
    _arrow(ax, (0, 0), A @ x, color=ORANGE, lw=1.4, ls="--")
    ax.text(*(v1 + [-0.62, 0.02]), r"$\mathbf{v}_1$", color=BLUE, fontsize=13)
    ax.text(
        *(A @ v1 + [-1.5, 0.42]),
        r"$\mathbf{A}\mathbf{v}_1=3\mathbf{v}_1$",
        color=BLUE,
        fontsize=12
    )
    ax.text(
        *(v2 + [0.08, -0.32]),
        r"$\mathbf{v}_2=\mathbf{A}\mathbf{v}_2$",
        color=GREEN,
        fontsize=12
    )
    ax.text(*(x + [0.06, 0.04]), r"$\mathbf{x}$", color=ORANGE, fontsize=13)
    ax.text(
        *(A @ x + [0.2, -0.36]),
        r"$\mathbf{A}\mathbf{x}$",
        color=ORANGE,
        fontsize=12
    )
    return save_mpl(fig, "01-01-circle-ellipse", outdir)


def quadratic_forms(outdir, lang="en"):
    """Contours of q(x) = x^T A x for three symmetric 2x2 matrices.

    (a) A = [[2, 1], [1, 2]], eigenvalues 1, 3: positive definite, closed
        ellipses around the unique minimiser 0.
    (b) A = [[1, -1], [-1, 1]], eigenvalues 0, 2: positive semidefinite and
        singular; q = (x1 - x2)^2 vanishes on the whole line x1 = x2 (green).
    (c) A = [[1, 2], [2, 1]], eigenvalues -1, 3: indefinite; q = 0 on the two
        lines x2 = (-2 +/- sqrt3) x1 (green), positive levels blue, negative
        levels orange (a saddle).
    Used in section sec-01-01-quadratic-forms of 01-01-linear-algebra.md.
    """
    mats = [
        (
            np.array([[2.0, 1.0], [1.0, 2.0]]),
            r"(a)  $\lambda_i(\mathbf{A})=1,\ 3$",
        ),
        (
            np.array([[1.0, -1.0], [-1.0, 1.0]]),
            r"(b)  $\lambda_i(\mathbf{A})=0,\ 2$",
        ),
        (
            np.array([[1.0, 2.0], [2.0, 1.0]]),
            r"(c)  $\lambda_i(\mathbf{A})=-1,\ 3$",
        ),
    ]
    g = np.linspace(-2, 2, 241)
    X1, X2 = np.meshgrid(g, g)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.7))
    for ax, (A, title) in zip(axes, mats):
        Q = A[0, 0] * X1**2 + 2 * A[0, 1] * X1 * X2 + A[1, 1] * X2**2
        pos = [0.5, 1.5, 3, 5, 8]
        neg = [-5, -3, -1.5, -0.5]
        ax.contour(X1, X2, Q, levels=pos, colors=BLUE, linewidths=1.2)
        if Q.min() < -0.1:
            ax.contour(X1, X2, Q, levels=neg, colors=ORANGE, linewidths=1.2)
        # zero level set: the origin (a), a line (b), two lines (c)
        if np.isclose(np.linalg.det(A), 0):
            ax.plot(g, g, color=GREEN, lw=2.2)
        elif np.linalg.det(A) < 0:
            for m in (-2 + np.sqrt(3), -2 - np.sqrt(3)):
                ax.plot(g, m * g, color=GREEN, lw=2.2)
        else:
            ax.plot(0, 0, "o", color=GREEN, ms=6)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ax.set_aspect("equal")
        ax.set_xticks([-2, -1, 0, 1, 2])
        ax.set_yticks([-2, -1, 0, 1, 2])
        ax.set_xlabel("$x_1$")
        ax.set_title(title, loc="left", fontsize=11.5)
    axes[0].set_ylabel("$x_2$")
    fig.tight_layout()
    return save_mpl(fig, "01-01-quadratic-forms", outdir)


def unit_balls(outdir, lang="en"):
    """Boundaries of the unit balls of the 1-, 2- and infinity-norms in R^2.

    |x1| + |x2| = 1 (diamond), x1^2 + x2^2 = 1 (circle) and
    max(|x1|, |x2|) = 1 (square); the balls are nested,
    B_1 inside B_2 inside B_inf. Used in section sec-01-01-norms.
    """
    fig, ax = plt.subplots(figsize=(4.6, 4.3))
    _plane_axes(ax, (-1.6, 1.6), (-1.6, 1.6))
    s = np.linspace(0, 2 * np.pi, 400)
    ax.plot([1, 0, -1, 0, 1], [0, 1, 0, -1, 0], color=BLUE)  # 1-norm
    ax.plot(np.cos(s), np.sin(s), color=ORANGE)  # 2-norm
    ax.plot([1, -1, -1, 1, 1], [1, 1, -1, -1, 1], color=GREEN)  # inf-norm
    # labels outside the square, with thin leader lines to each curve
    leader = dict(arrowstyle="-", lw=0.8, shrinkA=2, shrinkB=0)
    ax.annotate(
        r"$\|\mathbf{x}\|_1=1$",
        xy=(-0.5, 0.5),
        xytext=(-1.55, 1.3),
        color=BLUE,
        fontsize=11,
        arrowprops=dict(color=BLUE, **leader),
    )
    ax.annotate(
        r"$\|\mathbf{x}\|_2=1$",
        xy=(np.cos(0.9), np.sin(0.9)),
        xytext=(0.35, 1.3),
        color=ORANGE,
        fontsize=11,
        arrowprops=dict(color=ORANGE, **leader),
    )
    ax.annotate(
        r"$\|\mathbf{x}\|_\infty=1$",
        xy=(1.0, -0.75),
        xytext=(0.62, -1.42),
        color=GREEN,
        fontsize=11,
        arrowprops=dict(color=GREEN, **leader),
    )
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    return save_mpl(fig, "01-01-unit-balls", outdir)


FIGURES = {
    "01-01-span-combination": span_combination,
    "01-01-projection-line": projection_line,
    "01-01-circle-ellipse": circle_ellipse,
    "01-01-quadratic-forms": quadratic_forms,
    "01-01-unit-balls": unit_balls,
}
