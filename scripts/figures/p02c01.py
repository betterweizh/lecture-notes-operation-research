"""Figures of chapter 2.1, Convex Sets.

Six matplotlib figures, all from closed-form data (no randomness):

* ``02-01-convex-nonconvex``: a convex hexagon and three nonconvex sets, each
  with a segment (section sec-02-01-convex-sets);
* ``02-01-hulls``: convex hull of nine points and conic hull of two vectors;
* ``02-01-hyperplane``: a hyperplane, its normal vector and its half-spaces;
* ``02-01-balls``: unit balls of three norms and an ellipsoid;
* ``02-01-example-sets``: the sets Omega_1 ... Omega_5 of ex-02-01-seven-sets;
* ``02-01-separation``: separating hyperplane of a disc and a point, and a
  supporting hyperplane of a polygon.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.path import Path as MplPath
from matplotlib.patches import Circle, PathPatch, Polygon
from scipy.spatial import ConvexHull

from .common import BLUE, INK, MUTED, ORANGE, save_mpl

FILL = BLUE  # color of sets
FILL_ALPHA = 0.18
# Face color with the transparency built in, so that edges stay opaque.
FACE = to_rgba(FILL, FILL_ALPHA)


def _bare(ax, xlim, ylim):
    """Equal-aspect panel without axes, ticks or grid (pure geometry)."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.axis("off")


def _plane(ax, xlim, ylim, labels=("$x_1$", "$x_2$")):
    """Equal-aspect panel with coordinate axes, light grid and axis titles."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])


def _arrow(ax, start, end, color=INK, lw=1.8, zorder=5):
    """A vector drawn as an arrow from ``start`` to ``end`` (data coords)."""
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=lw,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=14,
        ),
        zorder=zorder,
    )


def _segment_in_out(ax, p, q, inside, color=ORANGE):
    """Draw the segment [p, q]: solid where ``inside`` holds, dashed elsewhere.

    ``inside`` maps an (N, 2) array of points to N booleans. The segment is
    sampled at 801 points and drawn as runs of equal membership.
    """
    t = np.linspace(0.0, 1.0, 801)
    pts = np.outer(1 - t, p) + np.outer(t, q)
    flags = inside(pts)
    start = 0
    for k in range(1, len(t) + 1):
        if k == len(t) or flags[k] != flags[start]:
            run = pts[start : k + 1] if k < len(t) else pts[start:k]
            ax.plot(run[:, 0], run[:, 1], color=color, lw=2.0,
                    ls="-" if flags[start] else (0, (3, 2.5)), zorder=4)  # fmt: skip
            start = k
    ax.plot(*np.array([p, q]).T, "o", color=color, ms=5, zorder=5)


def _hexagon(center, radius, phase=0.0):
    """Vertices of a regular hexagon (phase 0: a vertex points right)."""
    ang = phase + np.arange(6) * np.pi / 3
    return np.c_[
        center[0] + radius * np.cos(ang), center[1] + radius * np.sin(ang)
    ]


def convex_nonconvex(outdir, lang="en"):
    """One convex and three nonconvex sets, each with a segment.

    (a) regular hexagon, segment inside; (b) kidney-shaped set given in
    polar form r(t) = 1 - 0.5 max(0, sin t)^3 (a dent at the top), segment
    across the dent; (c) ring 0.55 <= ||x|| <= 1, segment across the hole;
    (d) two hexagons touching at the vertex (1, 0), segment across the gap.
    Parts of a segment outside the set are dashed. Used in sec-02-01-convex-sets.
    """
    fig, axes = plt.subplots(1, 4, figsize=(9.6, 2.7))
    lim = (-1.15, 1.15)

    # (a) convex hexagon: the segment stays inside
    ax = axes[0]
    hexa = _hexagon((0, 0), 1.0)
    ax.add_patch(Polygon(hexa, fc=FACE, ec="none"))
    ax.add_patch(Polygon(hexa, fill=False, ec=INK, lw=1.4))
    _segment_in_out(
        ax, (-0.55, -0.45), (0.6, 0.55), MplPath(hexa).contains_points
    )

    # (b) kidney: polar curve with a dent at t = pi/2
    ax = axes[1]
    t = np.linspace(0, 2 * np.pi, 400)
    r = 1 - 0.5 * np.maximum(0, np.sin(t)) ** 3
    kidney = np.c_[r * np.cos(t), r * np.sin(t)]
    ax.add_patch(Polygon(kidney, fc=FACE, ec="none"))
    ax.add_patch(Polygon(kidney, fill=False, ec=INK, lw=1.4))
    # both end points lie in the set (|x| = 0.83 < r = 0.85 there); the
    # middle of the segment passes above the dent (r = 0.5 at the top)
    _segment_in_out(
        ax, (-0.62, 0.55), (0.62, 0.55), MplPath(kidney).contains_points
    )

    # (c) ring between radii 0.55 and 1
    ax = axes[2]
    # the ring as one path: outer circle counter-clockwise, inner clockwise
    tt = np.linspace(0, 2 * np.pi, 200)
    outer = np.c_[np.cos(tt), np.sin(tt)]
    inner = 0.55 * np.c_[np.cos(tt), -np.sin(tt)]
    verts = np.vstack([outer, inner])
    codes = np.r_[
        [MplPath.MOVETO],
        [MplPath.LINETO] * 199,
        [MplPath.MOVETO],
        [MplPath.LINETO] * 199,
    ]
    ax.add_patch(PathPatch(MplPath(verts, codes), fc=FACE, ec="none"))
    ax.add_patch(Circle((0, 0), 1.0, fill=False, ec=INK, lw=1.4))
    ax.add_patch(Circle((0, 0), 0.55, fill=False, ec=INK, lw=1.4))

    def in_ring(pts):
        n = np.hypot(pts[:, 0], pts[:, 1])
        return (n >= 0.55) & (n <= 1.0)

    _segment_in_out(ax, (-0.78, -0.1), (0.78, 0.1), in_ring)

    # (d) two hexagons with a common vertex (1, 0); scaled to fit the panel
    ax = axes[3]
    h1 = _hexagon((-0.5, 0), 0.5)
    h2 = _hexagon((0.5, 0), 0.5)
    for h in (h1, h2):
        ax.add_patch(Polygon(h, fc=FACE, ec="none"))
        ax.add_patch(Polygon(h, fill=False, ec=INK, lw=1.4))

    def in_two(pts):
        return MplPath(h1).contains_points(pts) | MplPath(h2).contains_points(
            pts
        )

    _segment_in_out(ax, (-0.45, 0.3), (0.45, 0.3), in_two)

    for ax, tag in zip(axes, "abcd"):
        _bare(ax, lim, lim)
        ax.set_title(f"({tag})")
    return save_mpl(fig, "02-01-convex-nonconvex", outdir)


def hulls(outdir, lang="en"):
    """Convex hull of a finite set and conic hull of two vectors.

    (a) nine points; the hull (scipy ConvexHull) is shaded, its vertices
    (the extreme points) are filled, the other points are open circles.
    (b) x_1 = (2, 0.6) and x_2 = (0.7, 2); the wedge of their nonnegative
    combinations is shaded up to the window edge. Used in sec-02-01-convex-sets
    and sec-02-01-cones.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.6))

    # (a) nine points: six on the hull, three inside / on an edge
    pts = np.array(
        [[0.2, 0.6], [1.4, 0.1], [2.8, 0.5], [3.3, 1.8], [2.1, 2.9],
         [0.6, 2.4], [1.3, 1.2], [2.2, 1.5], [1.6, 2.1]]
    )  # fmt: skip
    hull = ConvexHull(pts)
    poly = pts[hull.vertices]
    ax1.add_patch(Polygon(poly, fc=FACE, ec=INK, lw=1.4))
    on_hull = np.zeros(len(pts), bool)
    on_hull[hull.vertices] = True
    ax1.plot(*pts[on_hull].T, "o", color=INK, ms=6, zorder=5)
    ax1.plot(*pts[~on_hull].T, "o", mfc="white", mec=INK, ms=6, zorder=5)
    ax1.text(1.9, 1.0, r"$\operatorname{conv} S$", ha="center", color=INK)
    _bare(ax1, (-0.1, 3.6), (-0.2, 3.2))
    ax1.set_title("(a)")

    # (b) conic hull of x_1, x_2: wedge from the origin, clipped to the window
    x1 = np.array([2.0, 0.6])
    x2 = np.array([0.7, 2.0])
    big = 10.0
    wedge = np.array([[0, 0], big * x1, big * x2])
    ax2.add_patch(Polygon(wedge, fc=FACE, ec="none"))
    for v in (x1, x2):  # boundary rays, thin and muted
        ax2.plot([0, big * v[0]], [0, big * v[1]], color=MUTED, lw=1.0)
    _arrow(ax2, (0, 0), x1)
    _arrow(ax2, (0, 0), x2)
    ax2.text(x1[0] + 0.05, x1[1] - 0.32, r"$\mathbf{x}_1$")
    ax2.text(x2[0] - 0.42, x2[1] + 0.02, r"$\mathbf{x}_2$")
    ax2.text(
        2.35,
        2.25,
        r"$\operatorname{cone}\{\mathbf{x}_1,\mathbf{x}_2\}$",
        ha="center",
    )
    ax2.plot(0, 0, "o", color=INK, ms=4)
    ax2.text(-0.12, -0.25, r"$\mathbf{0}$")
    _bare(ax2, (-0.4, 3.6), (-0.4, 3.2))
    ax2.set_title("(b)")
    return save_mpl(fig, "02-01-hulls", outdir)


def hyperplane(outdir, lang="en"):
    """Hyperplane a^T x = b with normal a = (1, 2), b = 4, through x0 = (2, 1).

    The half-space a^T x <= b (away from a) is shaded. A point x = (0.6, 0.9)
    of that half-space is joined to x0 by a dashed vector x - x0, which makes
    an obtuse angle with a. Used in sec-02-01-hyperplanes.
    """
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    a = np.array([1.0, 2.0])
    b = 4.0
    x0 = np.array([2.0, 1.0])  # a^T x0 = 4
    xs = np.array([-0.5, 4.5])
    line = (b - a[0] * xs) / a[1]  # x2 on the hyperplane
    # shade a^T x <= b: below the line inside the window
    ax.fill_between(xs, -0.5, line, color=FILL, alpha=FILL_ALPHA, lw=0)
    ax.plot(xs, line, color=INK, lw=1.8)
    # normal vector at x0 (scaled by 0.55 for the picture)
    _arrow(ax, x0, x0 + 0.55 * a, color=ORANGE)
    ax.text(
        *(x0 + 0.55 * a + np.array([0.08, 0.0])), r"$\mathbf{a}$", color=ORANGE
    )
    # a point x of the shaded half-space and the vector x - x0
    x = np.array([0.6, 0.9])
    _arrow(ax, x0, x, color=MUTED, lw=1.4)
    ax.plot(*x0, "o", color=INK, ms=5, zorder=6)
    ax.plot(*x, "o", color=INK, ms=5, zorder=6)
    ax.text(x0[0] + 0.1, x0[1] - 0.3, r"$\mathbf{x}_0$")
    ax.text(x[0] - 0.15, x[1] - 0.32, r"$\mathbf{x}$")
    # labels of the hyperplane and the two half-spaces
    ax.text(3.45, 0.62, r"$\mathbf{a}^\mathsf{T}\mathbf{x}=b$", ha="left")
    ax.text(0.7, 0.05, r"$\mathbf{a}^\mathsf{T}\mathbf{x}\leq b$", ha="center")
    ax.text(2.8, 2.6, r"$\mathbf{a}^\mathsf{T}\mathbf{x}\geq b$", ha="center")
    _plane(ax, (-0.5, 4.5), (-0.5, 3.2))
    return save_mpl(fig, "02-01-hyperplane", outdir)


def balls(outdir, lang="en"):
    """Unit balls of three norms and an ellipsoid.

    (a) {||x||_1 <= 1} (diamond), {||x||_2 <= 1} (disc), {||x||_inf <= 1}
    (square). (b) the ellipsoid (x - xc)^T P^{-1} (x - xc) <= 1 with
    xc = (2, 1), P = V diag(4, 1) V^T, V the rotation by 30 degrees; the
    semi-axes 2 v_1 and 1 v_2 are drawn from the center. Used in
    sec-02-01-balls.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.4, 3.7))

    # (a) unit balls, drawn from the largest to the smallest
    sq = np.array([[1, 1], [-1, 1], [-1, -1], [1, -1]])
    dia = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]])
    ax1.add_patch(Polygon(sq, fc="none", ec=ORANGE, lw=1.8))
    ax1.add_patch(Circle((0, 0), 1, fc="none", ec=BLUE, lw=1.8))
    ax1.add_patch(Polygon(dia, fc=FACE, ec=INK, lw=1.8))
    ax1.text(
        1.1, 0.95, r"$\|\mathbf{x}\|_\infty\leq 1$", color=ORANGE, va="center"
    )
    # the 2-norm label sits right of the square, with a pointer to the circle
    ax1.annotate(
        r"$\|\mathbf{x}\|_2\leq 1$", xy=(np.cos(0.35), np.sin(0.35)), xytext=(1.1, 0.25),
        color=BLUE, va="center", arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.9),
    )  # fmt: skip
    ax1.text(0.0, 0.0, r"$\|\mathbf{x}\|_1\leq 1$", ha="center", va="center")
    _plane(ax1, (-1.5, 2.4), (-1.5, 1.5))
    ax1.set_title("(a)")

    # (b) ellipsoid: P = V diag(4, 1) V^T, semi-axes sqrt(4) = 2 and 1
    phi = np.radians(30)
    V = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    lam = np.array([4.0, 1.0])  # eigenvalues of P
    xc = np.array([2.0, 1.0])
    t = np.linspace(0, 2 * np.pi, 300)
    # E = {xc + P^{1/2} u : ||u|| <= 1}, P^{1/2} = V diag(sqrt(lam)) V^T
    Phalf = V @ np.diag(np.sqrt(lam)) @ V.T
    pts = xc[:, None] + Phalf @ np.vstack([np.cos(t), np.sin(t)])
    ax2.add_patch(Polygon(pts.T, fc=FACE, ec=INK, lw=1.8))
    ax2.plot(*xc, "o", color=INK, ms=5, zorder=6)
    ax2.text(xc[0] - 0.05, xc[1] - 0.3, r"$\mathbf{x}_c$", ha="right")
    for k, col in zip(range(2), (ORANGE, BLUE)):
        tip = xc + np.sqrt(lam[k]) * V[:, k]
        _arrow(ax2, xc, tip, color=col)
    # semi-axes 2 v_1 and 1 v_2 (v_k = unit eigenvectors, columns of V)
    ax2.text(3.8, 2.05, r"$2\,\mathbf{v}_1$", color=ORANGE, ha="left")
    ax2.text(1.3, 2.05, r"$\mathbf{v}_2$", color=BLUE, ha="center")
    _plane(ax2, (-0.2, 5.2), (-0.6, 2.6))
    ax2.set_title("(b)")
    return save_mpl(fig, "02-01-balls", outdir)


def example_sets(outdir, lang="en"):
    """The planar sets Omega_1 ... Omega_5 of ex-02-01-seven-sets.

    Omega_1: polygon (0,0), (1,0), (3/2,1/2), (0,2). Omega_2: closed disc,
    center (1,1), radius 1. Omega_3: nonconvex polygon (0,0), (1,1), (3,0),
    (0,3) with the segment (0,0)-(3,0) and its midpoint (1.5, 0) outside.
    Omega_4: Omega_1 without its two slanted edges (dashed) and the three
    corners on them (open circles). Omega_5: half-plane x1 + x2 >= 2.
    """
    fig, axes = plt.subplots(1, 5, figsize=(11.5, 2.9))
    win = (-0.3, 3.2)

    def region(ax, verts):
        ax.add_patch(Polygon(verts, fc=FACE, ec="none"))

    # Omega_1: closed polygon, all edges solid
    o1 = np.array([[0, 0], [1, 0], [1.5, 0.5], [0, 2]])
    region(axes[0], o1)
    axes[0].add_patch(Polygon(o1, fill=False, ec=INK, lw=1.6))

    # Omega_2: closed disc
    axes[1].add_patch(Circle((1, 1), 1, fc=FACE, ec=INK, lw=1.6))

    # Omega_3: nonconvex quadrilateral and the segment that leaves it
    o3 = np.array([[0, 0], [1, 1], [3, 0], [0, 3]])
    region(axes[2], o3)
    axes[2].add_patch(Polygon(o3, fill=False, ec=INK, lw=1.6))
    axes[2].plot(
        [0, 3], [0, 0], color=ORANGE, lw=2.0, ls=(0, (3, 2.5)), zorder=4
    )
    axes[2].plot([0, 1.5, 3], [0, 0, 0], "o", color=ORANGE, ms=5, zorder=5)

    # Omega_4: solid edges on the axes, dashed open edges, open corners
    region(axes[3], o1)
    axes[3].plot([0, 0, 1], [2, 0, 0], color=INK, lw=1.6)
    axes[3].plot([1, 1.5, 0], [0, 0.5, 2], color=INK, lw=1.6, ls=(0, (4, 3)))
    axes[3].plot(
        [1, 1.5, 0], [0, 0.5, 2], "o", mfc="white", mec=INK, ms=6, zorder=5
    )

    # Omega_5: half-plane above the line x1 + x2 = 2, clipped to the window
    hp = np.array(
        [[2 - win[1], win[1]], [win[1], 2 - win[1]], [win[1], win[1]]]
    )
    region(axes[4], hp)
    xs = np.array(win)
    axes[4].plot(xs, 2 - xs, color=INK, lw=1.6)

    for k, ax in enumerate(axes, start=1):
        ax.set_xlim(*win)
        ax.set_ylim(*win)
        ax.set_aspect("equal")
        ax.set_xticks([0, 1, 2, 3])
        ax.set_yticks([0, 1, 2, 3])
        ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
        ax.axvline(0, color=MUTED, lw=0.8, zorder=1)
        ax.set_title(rf"$\Omega_{k}$")
    return save_mpl(fig, "02-01-example-sets", outdir)


def separation(outdir, lang="en"):
    """Separating and supporting hyperplanes.

    (a) C = disc with center (1, 1) and radius 1, y = (3, 2). The projection
    x* = c + (y - c)/||y - c|| = (1.894, 1.447); the normal a = y - x*; the
    hyperplane a^T x = b passes through the midpoint of x* and y. (b) A
    convex polygon and the supporting line at its vertex x0 = (2.5, 1.5)
    with outward normal a = (1, 1): the polygon lies in a^T x <= a^T x0.
    Used in sec-02-01-separation.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.8))

    # (a) disc and point
    c = np.array([1.0, 1.0])
    y = np.array([3.0, 2.0])
    xs = c + (y - c) / np.linalg.norm(y - c)  # projection onto the disc
    a = y - xs
    mid = (xs + y) / 2
    ax1.add_patch(Circle(c, 1, fc=FACE, ec=INK, lw=1.6))
    ax1.text(0.85, 0.9, "$C$", ha="center", fontsize=13)
    ax1.plot(*np.array([xs, y]).T, color=MUTED, lw=1.2, ls="--")
    # hyperplane through mid, direction perpendicular to a
    d = np.array([-a[1], a[0]]) / np.linalg.norm(a)
    seg = np.array([mid - 1.9 * d, mid + 1.5 * d])
    ax1.plot(*seg.T, color=INK, lw=1.8)
    # normal vector a drawn at a point of the hyperplane above the midpoint
    foot = mid + 0.9 * d
    tip = foot + 0.8 * a / np.linalg.norm(a)
    _arrow(ax1, foot, tip, color=ORANGE)
    ax1.text(tip[0] + 0.05, tip[1] + 0.05, r"$\mathbf{a}$", color=ORANGE)
    ax1.plot(*xs, "o", color=INK, ms=5, zorder=6)
    ax1.plot(*y, "o", color=INK, ms=5, zorder=6)
    ax1.text(xs[0] - 0.1, xs[1] + 0.12, r"$\mathbf{x}^*$", ha="right")
    ax1.text(y[0] + 0.12, y[1] - 0.22, r"$\mathbf{y}$")
    ax1.text(
        seg[1][0] + 0.12,
        seg[1][1] - 0.08,
        r"$\mathbf{a}^\mathsf{T}\mathbf{x}=b$",
    )
    _plane(ax1, (-0.2, 3.8), (-0.2, 3.2))
    ax1.set_title("(a)")

    # (b) polygon and a supporting line at a vertex
    poly = np.array(
        [[0.4, 0.3], [2.2, 0.4], [2.5, 1.5], [1.2, 2.3], [0.3, 1.6]]
    )
    x0 = poly[2]
    n = np.array([1.0, 1.0]) / np.sqrt(2)  # outward normal at x0
    ax2.add_patch(Polygon(poly, fc=FACE, ec=INK, lw=1.6))
    ax2.text(1.3, 1.2, "$C$", ha="center", fontsize=13)
    t = np.array([-1.6, 1.4])
    dline = np.array([-n[1], n[0]])
    ax2.plot(x0[0] + t * dline[0], x0[1] + t * dline[1], color=INK, lw=1.8)
    _arrow(ax2, x0, x0 + 0.8 * n, color=ORANGE)
    ax2.text(
        *(x0 + 0.8 * n + np.array([0.05, 0.02])), r"$\mathbf{a}$", color=ORANGE
    )
    ax2.plot(*x0, "o", color=INK, ms=5, zorder=6)
    ax2.text(x0[0] + 0.12, x0[1] - 0.2, r"$\mathbf{x}_0$")
    _plane(ax2, (-0.2, 3.8), (-0.2, 3.2))
    ax2.set_title("(b)")
    return save_mpl(fig, "02-01-separation", outdir)


FIGURES = {
    "02-01-convex-nonconvex": convex_nonconvex,
    "02-01-hulls": hulls,
    "02-01-hyperplane": hyperplane,
    "02-01-balls": balls,
    "02-01-example-sets": example_sets,
    "02-01-separation": separation,
}
