"""Render a deterministic square style specimen."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from catplot import apply_style, mm_to_inches, save_panel


def main() -> None:
    theme = apply_style("paper", "default")
    rng = np.random.default_rng(24082026)
    edge_in = mm_to_inches(170.0)
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(edge_in, edge_in),
        layout="constrained",
    )

    x = np.linspace(0, 10, 250)
    y = np.exp(-0.12 * x) * np.sin(1.45 * x)
    axes[0, 0].plot(x, y, label="Calculation")
    axes[0, 0].fill_between(x, y - 0.08, y + 0.08, alpha=0.20, linewidth=0)
    axes[0, 0].scatter(x[::25], y[::25] + rng.normal(0, 0.045, 10),
                       color=theme["semantic"]["experiment"], label="Experiment",
                       zorder=3)
    axes[0, 0].axhline(0, color=theme["semantic"]["reference"], linewidth=0.7)
    axes[0, 0].set(xlabel="Time / ps", ylabel="Correlation / -", title="Line and uncertainty")
    axes[0, 0].legend(ncols=2, loc="upper right")

    sx = rng.normal(size=90)
    sy = 0.65 * sx + rng.normal(scale=0.75, size=sx.size)
    sc = axes[0, 1].scatter(sx, sy, c=np.hypot(sx, sy), cmap="viridis",
                            edgecolors="white", linewidths=0.35)
    fig.colorbar(sc, ax=axes[0, 1], label=r"Magnitude / $\AA$")
    axes[0, 1].set(xlabel="Descriptor A", ylabel="Descriptor B", title="Continuous colour")

    categories = ["P-S", "P-Se", "Sb-S", "Sb-Se"]
    values = np.array([0.40, 0.32, 0.34, 0.24])
    errors = np.array([0.025, 0.018, 0.020, 0.015])
    positions = np.arange(len(categories))
    axes[1, 0].errorbar(positions, values, yerr=errors, fmt="o", linestyle="none")
    axes[1, 0].set_xticks(positions, categories)
    axes[1, 0].set(xlabel="Composition", ylabel=r"$E_a$ / eV", title="Categorical estimate")

    gx = np.linspace(-3, 3, 160)
    gy = np.linspace(-3, 3, 160)
    xx, yy = np.meshgrid(gx, gy)
    zz = np.sin(xx**2 + yy**2) * np.exp(-0.10 * (xx**2 + yy**2))
    mesh = axes[1, 1].pcolormesh(gx, gy, zz, shading="auto", cmap="coolwarm", rasterized=True)
    axes[1, 1].contour(gx, gy, zz, levels=7, colors="0.15", linewidths=0.45, alpha=0.65)
    fig.colorbar(mesh, ax=axes[1, 1], label="Response / a.u.")
    axes[1, 1].set_box_aspect(1)
    axes[1, 1].set(xlabel=r"$q_x$ / $\AA^{-1}$", ylabel=r"$q_y$ / $\AA^{-1}$", title="Raster map, vector axes")

    for label, ax in zip("abcd", axes.flat):
        ax.text(-0.14, 1.04, label, transform=ax.transAxes, fontsize=9,
                fontweight="bold", va="bottom", ha="left")

    root = Path(__file__).resolve().parents[1]
    metadata = {"Title": "CATPlot style specimen", "Creator": "catplot 0.2.0.dev1"}
    save_panel(fig, root / "output/pdf/style_specimen.pdf", metadata=metadata)
    save_panel(fig, root / "output/preview/style_specimen.png", metadata=metadata)
    plt.close(fig)


if __name__ == "__main__":
    main()
