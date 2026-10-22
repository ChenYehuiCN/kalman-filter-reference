#!/usr/bin/env python3
"""
Visualize the output of example/eskf_localization/main.c.

Usage:
    python3 plot.py [CSV] [--outdir DIR] [--show]

Figures written next to the CSV by default:

    eskf_trajectory.png  -- path, landmarks, dead reckoning and 2-sigma ellipses
    eskf_errors.png      -- ESKF error and 2-sigma envelopes
    eskf_vs_dr.png       -- ESKF versus dead reckoning (log error scale)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

if "--show" not in sys.argv[1:]:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Patch

# ---------------------------------------------------------------------------
# Style and helpers
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent

C_TRUE = "#262626"
C_EST = "#1f77b4"
C_DR = "#e07b39"
C_LM = "#d62728"
C_GRID = "#e7e7e7"

# Must match the `model` initializer in main.c.
LANDMARKS = np.array(
    [
        [14.0, 9.0],
        [-11.0, 16.0],
        [-15.0, -18.0],
    ]
)

plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "#fcfcfc",
        "axes.edgecolor": "#c9c9c9",
        "axes.linewidth": 0.8,
        "axes.grid": True,
        "grid.color": C_GRID,
        "grid.linewidth": 0.8,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
        "legend.framealpha": 0.9,
        "savefig.dpi": 160,
        "figure.dpi": 110,
    }
)


def wrap_angle(angle):
    return (angle + np.pi) % (2.0 * np.pi) - np.pi


def rmse(values):
    return float(np.sqrt(np.mean(np.square(values))))


def add_covariance_ellipse(ax, x, y, pxx, pxy, pyy, nstd=2.0, **kwargs):
    cov = np.array([[pxx, pxy], [pxy, pyy]], dtype=float)
    if not np.all(np.isfinite(cov)):
        return
    values, vectors = np.linalg.eigh(cov)
    if values[0] <= 0.0:
        return
    order = values.argsort()[::-1]
    values = values[order]
    vectors = vectors[:, order]
    angle = np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0]))
    width, height = 2.0 * nstd * np.sqrt(values)
    ax.add_patch(
        Ellipse((x, y), width=width, height=height, angle=angle, **kwargs)
    )


def load_csv(path: Path) -> np.ndarray:
    if not path.is_file():
        raise SystemExit(
            f"CSV not found: {path}\n"
            "Run `make run` in example/eskf_localization first."
        )
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.shape == ():
        raise SystemExit(f"CSV contains no samples: {path}")
    return data


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def plot_trajectory(data, outdir: Path):
    fig, ax = plt.subplots(figsize=(8.9, 8.0), constrained_layout=True)

    ax.plot(
        data["true_x"], data["true_y"], color=C_TRUE, lw=2.0, label="truth",
        zorder=4,
    )
    ax.plot(
        data["est_x"], data["est_y"], color=C_EST, lw=1.5,
        label="ESKF estimate", zorder=5,
    )
    ax.plot(
        data["dr_x"], data["dr_y"], color=C_DR, lw=1.2, ls="--",
        label="dead reckoning", zorder=3,
    )

    # Covariance: dense during the initial convergence, sparse afterwards.
    indices = list(range(0, min(201, len(data)), 10))
    indices += list(range(300, len(data), 500))
    for i in indices:
        add_covariance_ellipse(
            ax,
            data["est_x"][i],
            data["est_y"][i],
            data["Pxx"][i],
            data["Pxy"][i],
            data["Pyy"][i],
            nstd=2.0,
            facecolor=C_EST,
            edgecolor=C_EST,
            alpha=0.10,
            linewidth=0.6,
            zorder=1,
        )

    for i, (lx, ly) in enumerate(LANDMARKS):
        ax.scatter(
            lx,
            ly,
            marker="*",
            s=170,
            color=C_LM,
            edgecolor="white",
            linewidth=0.8,
            zorder=7,
            label="landmarks" if i == 0 else None,
        )
        ax.annotate(
            f"L{i + 1}",
            (lx, ly),
            textcoords="offset points",
            xytext=(9, 6),
            color=C_LM,
            fontsize=10,
            weight="bold",
        )

    ax.scatter(
        data["true_x"][0],
        data["true_y"][0],
        s=75,
        color="#2a9d4b",
        edgecolor="white",
        linewidth=0.8,
        zorder=8,
        label="start (truth)",
    )
    ax.scatter(
        data["est_x"][0],
        data["est_y"][0],
        marker="s",
        s=65,
        facecolor="none",
        edgecolor=C_EST,
        linewidth=1.2,
        zorder=8,
        label="initial estimate",
    )

    handles, labels = ax.get_legend_handles_labels()
    handles.append(
        Patch(facecolor=C_EST, edgecolor=C_EST, alpha=0.18, label=r"2$\sigma$ position")
    )
    labels.append(r"2$\sigma$ position")

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(
        "ESKF localization: unicycle motion, body-frame landmark observations"
    )
    ax.set_xlim(-20.0, 22.0)
    ax.set_ylim(-24.0, 19.0)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(handles, labels, loc="upper right", ncol=2)
    fig.savefig(outdir / "eskf_trajectory.png", bbox_inches="tight")
    plt.close(fig)


def plot_errors(data, outdir: Path):
    t = data["t"]
    ex = data["est_x"] - data["true_x"]
    ey = data["est_y"] - data["true_y"]
    eyaw = np.degrees(wrap_angle(data["est_yaw"] - data["true_yaw"]))
    sx = np.sqrt(np.maximum(data["Pxx"], 0.0))
    sy = np.sqrt(np.maximum(data["Pyy"], 0.0))
    syaw = np.degrees(np.sqrt(np.maximum(data["Pyaw"], 0.0)))

    rows = [
        ("Position x", ex, sx, "m"),
        ("Position y", ey, sy, "m"),
        ("Yaw", eyaw, syaw, "deg"),
    ]

    fig, axes = plt.subplots(
        3, 2, figsize=(12.6, 8.4), constrained_layout=True
    )
    zoom_start = 1.0

    for row, (title, error, sigma, unit) in enumerate(rows):
        for col, mask in enumerate(
            [np.ones_like(t, dtype=bool), t >= zoom_start]
        ):
            ax = axes[row, col]
            tt = t[mask]
            ee = error[mask]
            ss = sigma[mask]

            ax.fill_between(
                tt,
                -2.0 * ss,
                2.0 * ss,
                color=C_EST,
                alpha=0.18,
                linewidth=0,
                label=r"$\pm 2\sigma$",
            )
            ax.plot(tt, ee, color=C_EST, lw=1.1, label="ESKF error")
            ax.axhline(0.0, color="#9a9a9a", lw=0.8, zorder=0)
            if col == 0:
                ax.axvline(
                    zoom_start, color="#777777", lw=0.9, ls=":", zorder=0
                )
                ax.set_ylabel(f"error [{unit}]")
            ax.set_title(
                f"{title} — {'full run' if col == 0 else f't ≥ {zoom_start:g} s'}",
                fontsize=10,
            )
            ax.set_xlim(tt[0], tt[-1])
            if row < 2:
                ax.tick_params(labelbottom=False)

            inside = np.mean(np.abs(ee) <= 2.0 * np.maximum(ss, 1e-12))
            ax.text(
                0.02,
                0.94,
                f"RMSE {rmse(ee):.4f} {unit}\nwithin 2σ: {100 * inside:4.1f}%",
                transform=ax.transAxes,
                va="top",
                fontsize=8.2,
                bbox=dict(
                    boxstyle="round,pad=0.25", fc="white", ec="#dddddd", alpha=0.9
                ),
            )

    axes[0, 0].legend(loc="upper right")
    axes[-1, 0].set_xlabel("time [s]")
    axes[-1, 1].set_xlabel("time [s]")
    fig.suptitle(
        "ESKF estimation error (blue band = 2σ; right column zooms the steady state)",
        fontsize=11.5,
    )
    fig.savefig(outdir / "eskf_errors.png", bbox_inches="tight")
    plt.close(fig)


def plot_vs_dead_reckoning(data, outdir: Path):
    t = data["t"]
    eskf_pos = np.hypot(
        data["est_x"] - data["true_x"], data["est_y"] - data["true_y"]
    )
    dr_pos = np.hypot(data["dr_x"] - data["true_x"], data["dr_y"] - data["true_y"])
    eskf_yaw = np.abs(
        np.degrees(wrap_angle(data["est_yaw"] - data["true_yaw"]))
    )
    dr_yaw = np.abs(np.degrees(wrap_angle(data["dr_yaw"] - data["true_yaw"])))

    fig, axes = plt.subplots(
        2, 1, figsize=(11.0, 7.0), sharex=True, constrained_layout=True
    )

    axes[0].semilogy(
        t,
        np.maximum(eskf_pos, 1.0e-4),
        color=C_EST,
        lw=1.4,
        label=f"ESKF (RMSE {rmse(eskf_pos):.3f} m)",
    )
    axes[0].semilogy(
        t,
        np.maximum(dr_pos, 1.0e-4),
        color=C_DR,
        lw=1.2,
        ls="--",
        label=f"dead reckoning (RMSE {rmse(dr_pos):.3f} m)",
    )
    axes[0].set_ylabel("position error [m]")
    axes[0].set_title("Position error magnitude")
    axes[0].legend(loc="upper left")

    axes[1].semilogy(
        t,
        np.maximum(eskf_yaw, 1.0e-4),
        color=C_EST,
        lw=1.4,
        label=f"ESKF (RMSE {rmse(eskf_yaw):.3f}°)",
    )
    axes[1].semilogy(
        t,
        np.maximum(dr_yaw, 1.0e-4),
        color=C_DR,
        lw=1.2,
        ls="--",
        label=f"dead reckoning (RMSE {rmse(dr_yaw):.3f}°)",
    )
    axes[1].set_ylabel("yaw error [deg]")
    axes[1].set_xlabel("time [s]")
    axes[1].set_title("Yaw error magnitude")
    axes[1].legend(loc="upper left")

    fig.suptitle(
        "ESKF vs dead reckoning (log scale; both start from the same prior)",
        fontsize=11.5,
    )
    fig.savefig(outdir / "eskf_vs_dr.png", bbox_inches="tight")
    plt.close(fig)


def print_report(data):
    ex = data["est_x"] - data["true_x"]
    ey = data["est_y"] - data["true_y"]
    eyaw = np.degrees(wrap_angle(data["est_yaw"] - data["true_yaw"]))
    drx = data["dr_x"] - data["true_x"]
    dry = data["dr_y"] - data["true_y"]
    dryaw = np.degrees(wrap_angle(data["dr_yaw"] - data["true_yaw"]))

    eskf_pos = np.hypot(ex, ey)
    dr_pos = np.hypot(drx, dry)
    print(f"samples        : {len(data)}")
    print(f"duration       : {data['t'][-1]:.2f} s")
    print(f"ESKF updates   : {int(data['updated'].sum())}")
    print(f"ESKF position  : RMSE {rmse(eskf_pos):.4f} m, final {eskf_pos[-1]:.4f} m")
    print(f"ESKF yaw       : RMSE {rmse(eyaw):.4f}°, final {abs(eyaw[-1]):.4f}°")
    print(f"DR position    : RMSE {rmse(dr_pos):.4f} m, final {dr_pos[-1]:.4f} m")
    print(f"DR yaw         : RMSE {rmse(dryaw):.4f}°, final {abs(dryaw[-1]):.4f}°")
    if rmse(eskf_pos) > 0:
        print(f"position improvement over DR: {rmse(dr_pos) / rmse(eskf_pos):.1f}x")

    for name, error, var, unit in [
        ("position x", ex, data["Pxx"], "m"),
        ("position y", ey, data["Pyy"], "m"),
        ("yaw", eyaw, data["Pyaw"], "deg"),
    ]:
        sigma = np.sqrt(np.maximum(var, 0.0))
        if unit == "deg":
            sigma = np.degrees(sigma)
        inside = 100.0 * np.mean(np.abs(error) <= 2.0 * np.maximum(sigma, 1e-12))
        print(
            f"  {name:12s}: RMSE {rmse(error):8.4f} {unit:4s}   "
            f"within 2σ {inside:5.1f}%"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv",
        nargs="?",
        default=SCRIPT_DIR / "output" / "eskf_localization.csv",
        help="CSV produced by eskf_demo (default: output/eskf_localization.csv)",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="directory for PNG files (default: same directory as the CSV)",
    )
    parser.add_argument("--show", action="store_true", help="open figures interactively")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    data = load_csv(csv_path)
    outdir = Path(args.outdir) if args.outdir else csv_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    plot_trajectory(data, outdir)
    plot_errors(data, outdir)
    plot_vs_dead_reckoning(data, outdir)

    print_report(data)
    print(f"figures written to {outdir}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
