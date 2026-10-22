#!/usr/bin/env python3
"""
Visualize the output of example/kf_tracking/main.c.

Usage:
    python3 plot.py [CSV] [--outdir DIR] [--show]

The CSV is produced by the C program (make run).  Figures are written next to
the CSV by default:

    kf_trajectory.png  -- 2D path, measurements and 2-sigma ellipses
    kf_states.png      -- truth, estimate and 2-sigma envelopes
    kf_errors.png      -- estimation error and 2-sigma envelopes
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

# Keep the script usable on a headless machine unless the user asks to show.
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
C_MEAS = "#8c8c8c"
C_BAD = "#e07b39"
C_GRID = "#e7e7e7"
C_SHADE = "#f2c14e"

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

STATE_COLUMNS = [
    # (title, truth, estimate, covariance, unit)
    ("Position x", "true_px", "est_px", "Pxx", "m"),
    ("Position y", "true_py", "est_py", "Pyy", "m"),
    ("Velocity x", "true_vx", "est_vx", "Pvx", "m/s"),
    ("Velocity y", "true_vy", "est_vy", "Pvy", "m/s"),
    ("Acceleration x", "true_ax", "est_ax", "Pax", "m/s²"),
    ("Acceleration y", "true_ay", "est_ay", "Pay", "m/s²"),
]

POSITION_COLUMNS = {"est_px", "est_py"}


def add_covariance_ellipse(ax, x, y, pxx, pxy, pyy, nstd=2.0, **kwargs):
    """Draw a covariance ellipse.  Silently skip degenerate covariances."""
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


def degraded_mask(data):
    """Samples whose measurement noise is clearly above the normal level."""
    sigma_min = float(np.min(data["sigma"]))
    return data["sigma"] > 1.5 * sigma_min


def degradation_label(data):
    """Describe the time window in which the GPS noise is visibly larger."""
    mask = degraded_mask(data)
    if not np.any(mask):
        return "constant GPS noise"
    times = data["t"][mask]
    step = data["t"][1] - data["t"][0] if len(data) > 1 else 0.0
    return (
        f"orange window: t = {times.min():g}–{times.max() + step:g} s, "
        f"GPS σ = {data['sigma'].max():g} m"
    )


def shade_bad_window(ax, data):
    mask = degraded_mask(data)
    if not np.any(mask):
        return
    times = data["t"][mask]
    ax.axvspan(
        times.min(), times.max(), color=C_SHADE, alpha=0.12, linewidth=0, zorder=0
    )


def rmse(values):
    return float(np.sqrt(np.mean(np.square(values))))


def load_csv(path: Path) -> np.ndarray:
    if not path.is_file():
        raise SystemExit(
            f"CSV not found: {path}\nRun `make run` in example/kf_tracking first."
        )
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.shape == ():
        raise SystemExit(f"CSV contains no samples: {path}")
    return data


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def plot_trajectory(data, outdir: Path):
    fig, ax = plt.subplots(figsize=(10.0, 6.8), constrained_layout=True)

    bad = degraded_mask(data)
    good = ~bad
    good_sigma = data["sigma"][good].min() if np.any(good) else 0.0
    ax.scatter(
        data["meas_px"][good],
        data["meas_py"][good],
        s=5,
        color=C_MEAS,
        alpha=0.35,
        linewidths=0,
        label=rf"GPS, $\sigma={good_sigma:g}$ m",
        rasterized=True,
        zorder=2,
    )
    if np.any(bad):
        bad_sigma = data["sigma"][bad].max()
        ax.scatter(
            data["meas_px"][bad],
            data["meas_py"][bad],
            s=9,
            color=C_BAD,
            alpha=0.55,
            linewidths=0,
            label=rf"GPS, $\sigma={bad_sigma:g}$ m",
            rasterized=True,
            zorder=3,
        )

    ax.plot(
        data["true_px"],
        data["true_py"],
        color=C_TRUE,
        lw=2.0,
        label="truth",
        zorder=4,
    )
    ax.plot(
        data["est_px"],
        data["est_py"],
        color=C_EST,
        lw=1.4,
        label="KF estimate",
        zorder=5,
    )

    # 2-sigma position covariance every 2 seconds.
    stride = 100
    for i in range(0, len(data), stride):
        add_covariance_ellipse(
            ax,
            data["est_px"][i],
            data["est_py"][i],
            data["Pxx"][i],
            data["Pxy"][i],
            data["Pyy"][i],
            nstd=2.0,
            facecolor=C_EST,
            edgecolor=C_EST,
            alpha=0.13,
            linewidth=0.7,
            zorder=1,
        )

    ax.scatter(
        data["true_px"][0],
        data["true_py"][0],
        marker="o",
        s=70,
        color="#2a9d4b",
        edgecolor="white",
        linewidth=0.8,
        zorder=6,
        label="start (truth)",
    )
    ax.scatter(
        data["est_px"][0],
        data["est_py"][0],
        marker="s",
        s=60,
        facecolor="none",
        edgecolor=C_EST,
        linewidth=1.2,
        zorder=6,
        label="start (initial estimate)",
    )

    handles, labels = ax.get_legend_handles_labels()
    handles.append(
        Patch(facecolor=C_EST, edgecolor=C_EST, alpha=0.2, label=r"2$\sigma$ position")
    )
    labels.append(r"2$\sigma$ position")

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("KF tracking: 6-state CA model, known jerk + 2D position")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(handles, labels, loc="upper center", ncol=3)
    fig.savefig(outdir / "kf_trajectory.png", bbox_inches="tight")
    plt.close(fig)


def plot_states(data, outdir: Path):
    t = data["t"]
    fig, axes = plt.subplots(
        3, 2, figsize=(12.5, 8.2), sharex=True, constrained_layout=True
    )

    for idx, (title, truth_key, est_key, var_key, unit) in enumerate(STATE_COLUMNS):
        ax = axes.flat[idx]
        shade_bad_window(ax, data)

        if est_key in POSITION_COLUMNS:
            bad = degraded_mask(data)
            good = ~bad
            ax.scatter(
                t[good],
                data["meas_px" if est_key == "est_px" else "meas_py"][good],
                s=4,
                color=C_MEAS,
                alpha=0.3,
                linewidths=0,
                rasterized=True,
                zorder=1,
            )
            if np.any(bad):
                ax.scatter(
                    t[bad],
                    data["meas_px" if est_key == "est_px" else "meas_py"][bad],
                    s=7,
                    color=C_BAD,
                    alpha=0.5,
                    linewidths=0,
                    rasterized=True,
                    zorder=2,
                )

        mean = data[est_key]
        sigma = np.sqrt(np.maximum(data[var_key], 0.0))
        ax.fill_between(
            t,
            mean - 2.0 * sigma,
            mean + 2.0 * sigma,
            color=C_EST,
            alpha=0.17,
            linewidth=0,
            zorder=2,
        )
        ax.plot(t, data[truth_key], color=C_TRUE, lw=1.1, ls="--", label="truth")
        ax.plot(t, mean, color=C_EST, lw=1.5, label="KF estimate")
        ax.set_title(title)
        ax.set_ylabel(unit)
        if idx >= 4:
            ax.set_xlabel("time [s]")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles.append(Patch(facecolor=C_EST, alpha=0.25, label=r"estimate $\pm 2\sigma$"))
    labels.append(r"estimate $\pm 2\sigma$")
    axes[0, 0].legend(handles, labels, loc="upper left")

    fig.suptitle(
        f"KF state estimates ({degradation_label(data)})",
        fontsize=12,
    )
    fig.savefig(outdir / "kf_states.png", bbox_inches="tight")
    plt.close(fig)


def plot_errors(data, outdir: Path):
    # The initial convergence is already visible in kf_states.png; zoom the
    # error figure on the steady state so the 2-sigma bands stay readable.
    zoom_start = 2.0
    mask = data["t"] >= zoom_start
    t = data["t"][mask]
    fig, axes = plt.subplots(
        3, 2, figsize=(12.5, 8.2), sharex=True, constrained_layout=True
    )

    for idx, (title, truth_key, est_key, var_key, unit) in enumerate(STATE_COLUMNS):
        ax = axes.flat[idx]
        shade_bad_window(ax, data)

        error = (data[est_key] - data[truth_key])[mask]
        sigma = np.sqrt(np.maximum(data[var_key], 0.0))[mask]
        ax.fill_between(
            t,
            -2.0 * sigma,
            2.0 * sigma,
            color=C_EST,
            alpha=0.18,
            linewidth=0,
            label=r"$\pm 2\sigma$",
        )
        ax.plot(t, error, color=C_EST, lw=1.2, label="error")
        ax.axhline(0.0, color="#9a9a9a", lw=0.8, zorder=0)
        ax.set_title(title)
        ax.set_ylabel(unit)
        if idx >= 4:
            ax.set_xlabel("time [s]")

        inside = np.mean(np.abs(error) <= 2.0 * np.maximum(sigma, 1e-12))
        ax.text(
            0.02,
            0.96,
            f"RMSE {rmse(error):.3f} {unit}\nwithin 2σ: {100 * inside:4.1f}%",
            transform=ax.transAxes,
            va="top",
            fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#dddddd", alpha=0.9),
        )

    axes[0, 0].legend(loc="lower right")
    fig.suptitle(
        f"KF steady-state errors (t ≥ {zoom_start:g} s; blue band = 2σ; "
        "orange window = degraded GPS)",
        fontsize=12,
    )
    fig.savefig(outdir / "kf_errors.png", bbox_inches="tight")
    plt.close(fig)


def print_report(data):
    pos_error = np.hypot(
        data["est_px"] - data["true_px"], data["est_py"] - data["true_py"]
    )
    vel_error = np.hypot(
        data["est_vx"] - data["true_vx"], data["est_vy"] - data["true_vy"]
    )
    acc_error = np.hypot(
        data["est_ax"] - data["true_ax"], data["est_ay"] - data["true_ay"]
    )
    print(f"samples           : {len(data)}")
    print(f"duration          : {data['t'][-1]:.2f} s")
    print(f"position RMSE     : {rmse(pos_error):.3f} m (2D)")
    print(f"velocity RMSE     : {rmse(vel_error):.3f} m/s (2D)")
    print(f"acceleration RMSE : {rmse(acc_error):.3f} m/s² (2D)")
    for key, truth, var, name, unit in [
        ("x", "true_px", "Pxx", "position x", "m"),
        ("y", "true_py", "Pyy", "position y", "m"),
        ("vx", "true_vx", "Pvx", "velocity x", "m/s"),
        ("vy", "true_vy", "Pvy", "velocity y", "m/s"),
        ("ax", "true_ax", "Pax", "accel x", "m/s²"),
        ("ay", "true_ay", "Pay", "accel y", "m/s²"),
    ]:
        err = data[f"est_p{key}" if key in ("x", "y") else f"est_{key}"] - data[truth]
        sigma = np.sqrt(np.maximum(data[var], 0.0))
        inside = 100.0 * np.mean(np.abs(err) <= 2.0 * np.maximum(sigma, 1e-12))
        print(f"  {name:12s}: RMSE {rmse(err):8.3f} {unit:5s}   within 2σ {inside:5.1f}%")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv",
        nargs="?",
        default=SCRIPT_DIR / "output" / "kf_tracking.csv",
        help="CSV produced by kf_demo (default: output/kf_tracking.csv)",
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
    plot_states(data, outdir)
    plot_errors(data, outdir)

    print_report(data)
    print(f"figures written to {outdir}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
