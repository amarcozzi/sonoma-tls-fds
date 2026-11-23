#!/usr/bin/env python
"""Figures for sonoma-tls-fds, rendered from the saved postprocess outputs.

Plotting is kept apart from the analysis: every figure is drawn from the small
files postprocess.py writes (metrics.json, the CSVs, arrival_map.npz), so the
figures can be restyled and regenerated locally without touching the raw FDS
output on the HPC.

Per simulation, in <out>/<sim_id>/:
    summary.{pdf,png}   (a) arrival map, then against downwind distance:
                        (b) rate of spread, (c) flame length/height,
                        (d) flame tilt, (e) flame depth
Across simulations, in <out>/:
    comparison_ros.{pdf,png}
    comparison_distributions.{pdf,png}
    comparison_flame_length.{pdf,png}

Captions for all of these are in POSTPROCESSING.md.

Usage:
    python scripts/figures.py postprocess              # everything
    python scripts/figures.py postprocess --only c4_p25
"""

import argparse
import glob
import json
import logging
import os
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotstyle as S

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

UNBURNED_S = 1.0e5     # arrival >= this means never burned (FDS writes 1e6)


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------

def read_csv(path):
    """CSV with a header row -> {column: float array}, or None."""
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        header = fh.readline().strip().split(",")
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    if data.size == 0:
        return None
    data = np.atleast_2d(data)
    return {name: data[:, i] for i, name in enumerate(header)}


def load_record(sim_dir):
    path = os.path.join(sim_dir, "metrics.json")
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def front_columns(front):
    """Name of the front-time column (t_median_s or t_min_s)."""
    return next(c for c in front if c.startswith("t_")
                and not c.endswith(("_raw_s", "leading_s", "trailing_s")))


def finite(v):
    return v is not None and np.isfinite(v)


def roi_times(rec):
    sc = rec.get("scalars", {})
    t0, t1 = sc.get("t_enter_roi_s"), sc.get("t_exit_roi_s")
    return (t0 if finite(t0) else None), (t1 if finite(t1) else None)


def short_status(rec):
    return {"no_roi_arrival": "did not reach ROI",
            "stalled_in_roi": "stalled in ROI",
            "no_spread": "no spread",
            "error": "postprocessing failed"}.get(rec.get("status"), "")


# ---------------------------------------------------------------------------
# shared drawing helpers
# ---------------------------------------------------------------------------

def shade_roi_time(ax, rec):
    t0, t1 = roi_times(rec)
    if t0 is not None:
        ax.axvspan(t0, t1 if t1 is not None else ax.get_xlim()[1],
                   color=S.REGION, lw=0, zorder=0)


def headroom(ax, frac=0.42):
    """Zero-based y axis with free space on top for an in-axes legend."""
    ax.set_ylim(0, ax.get_ylim()[1] * (1 + frac))


def outlined(color="white", stroke=S.INK, width=1.6):
    return [pe.Stroke(linewidth=width, foreground=stroke), pe.Normal()]


# ---------------------------------------------------------------------------
# per-simulation summary figure
# ---------------------------------------------------------------------------

def front_x_at(front, t):
    """Median front position at times t (NaN outside the tracked record)."""
    tc = front_columns(front)
    ok = np.isfinite(front[tc])
    if ok.sum() < 2:
        return np.full(np.shape(t), np.nan)
    return np.interp(t, front[tc][ok], front["x_m"][ok],
                     left=np.nan, right=np.nan)


def x_limits(rec):
    """Shared downwind axis: ignition line to the x=+30 m kill devices."""
    x_ign = rec.get("scalars", {}).get("x_ignition_line_m")
    return ((x_ign if finite(x_ign) else -80.0) - 2.0, 30.0)


def per_sim(sim_dir):
    """Five-panel summary figure for one simulation. Returns the stem written.

    Everything except the map is plotted against downwind distance x, so the
    panels line up with the map and the ROI is the same band in all of them.
    """
    rec = load_record(sim_dir)
    if rec is None or rec.get("status") == "error":
        return None
    color = S.group_color(rec.get("group"))
    meta, sc = rec.get("meta", {}), rec.get("scalars", {})
    roi = meta.get("roi_x_m", [-10.0, 10.0])
    y_band = meta.get("y_band_m", 10.0)

    front = read_csv(os.path.join(sim_dir, "front_trajectory.csv"))
    ros_t = read_csv(os.path.join(sim_dir, "ros_time.csv"))
    ros_s = read_csv(os.path.join(sim_dir, "ros_segments.csv"))
    flame = read_csv(os.path.join(sim_dir, "flame_geometry.csv"))
    amap_path = os.path.join(sim_dir, "arrival_map.npz")

    fig = plt.figure(figsize=(S.DOUBLE, 5.0))
    outer = fig.add_gridspec(2, 1, height_ratios=(1.3, 1), hspace=0.38)
    top = outer[0].subgridspec(1, 2, wspace=0.3, width_ratios=(1.05, 1))
    bot = outer[1].subgridspec(1, 3, wspace=0.42)
    ax_map = fig.add_subplot(top[0, 0])
    ax_ros = fig.add_subplot(top[0, 1])
    ax_len = fig.add_subplot(bot[0, 0], sharex=ax_ros)
    ax_tilt = fig.add_subplot(bot[0, 1], sharex=ax_ros)
    ax_dep = fig.add_subplot(bot[0, 2], sharex=ax_ros)
    xlim = x_limits(rec)

    t_ign = sc.get("t_ignition_s", 0.0)
    t_end = sc.get("t_track_end_s")

    # (a) arrival map ------------------------------------------------------
    ax = ax_map
    if os.path.isfile(amap_path):
        d = np.load(amap_path)
        arr = np.ma.masked_where(d["arrival_s"] >= UNBURNED_S, d["arrival_s"])
        cmap = plt.get_cmap("magma_r").copy()
        cmap.set_bad(S.UNBURNED)
        vmax = t_end if finite(t_end) else float(arr.max())
        im = ax.pcolormesh(d["x_m"], d["y_m"], arr.T, cmap=cmap, vmin=t_ign,
                           vmax=vmax, shading="auto", rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlim(d["x_m"].min(), d["x_m"].max())
        ax.set_ylim(d["y_m"].min(), d["y_m"].max())
        ax.add_patch(Rectangle((roi[0], -y_band), roi[1] - roi[0], 2 * y_band,
                               fill=False, ec="white", lw=0.9, ls=(0, (3, 2)),
                               path_effects=outlined(), zorder=3))
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, aspect=22)
        cb.set_label("Arrival time (s)")
        cb.outline.set_linewidth(0.6)
        cb.ax.tick_params(width=0.6, length=2.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(True)
    else:
        ax.text(0.5, 0.5, "arrival_map.npz missing", ha="center",
                va="center", transform=ax.transAxes, color=S.MUTED)
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$y$ (m)")
    S.panel_label(ax, "a", x=-0.16)

    # (b) rate of spread ---------------------------------------------------
    ax = ax_ros
    ax.axvspan(roi[0], roi[1], color=S.REGION, lw=0, zorder=0)
    if ros_t is not None:
        ax.plot(ros_t["x_front_m"], ros_t["ros_m_per_s"], color=color,
                lw=1.0, label=f"Smoothed ({meta.get('ros_smooth_s', 5):g} s)")
    seg = meta.get("ros_segment_m", 2.5)
    if ros_s is not None:
        ax.hlines(ros_s["ros_m_per_s"], ros_s["x_mid_m"] - seg / 2,
                  ros_s["x_mid_m"] + seg / 2, color=S.INK, lw=1.4,
                  label=f"{seg:g} m segments")
    bulk = sc.get("ros_bulk_roi_m_per_s")
    if finite(bulk):
        ax.hlines(bulk, roi[0], roi[1], color=S.INK_2, lw=0.9,
                  ls=(0, (3, 2)), zorder=4,
                  label=f"ROI mean, {bulk:.2f} m s$^{{-1}}$")
    status = short_status(rec)
    if status:
        ax.text(0.98, 0.04, status, transform=ax.transAxes, ha="right",
                va="bottom", fontsize=7, color=S.INK_2, style="italic")
    ax.set_xlim(*xlim)
    headroom(ax, 0.35)
    ax.text(np.mean(roi), 0.99, "ROI", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=7, color=S.INK_2)
    ax.set_xlabel("Distance downwind, $x$ (m)")
    ax.set_ylabel("ROS (m s$^{-1}$)")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.93))
    S.panel_label(ax, "b", x=-0.12)

    # (c)-(e) flame geometry, mapped onto the front position at each time --
    for ax, letter in ((ax_len, "c"), (ax_tilt, "d"), (ax_dep, "e")):
        ax.axvspan(roi[0], roi[1], color=S.REGION, lw=0, zorder=0)
        ax.set_xlabel("Distance downwind, $x$ (m)")
        S.panel_label(ax, letter, x=-0.2)
    if flame is not None and front:
        xf = front_x_at(front, flame["t_s"])
        slices = [c for c in flame if c.startswith("flame_length_y")]
        if slices:
            stack = np.vstack([flame[c] for c in slices])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                lo, hi = np.nanmin(stack, 0), np.nanmax(stack, 0)
            ax_len.fill_between(xf, lo, hi, color=color, alpha=0.18, lw=0,
                                label="Range across slices")
        ax_len.plot(xf, flame["flame_length_mean"], color=color, lw=1.0,
                    label="Length, $L$")
        ax_len.plot(xf, flame["flame_height_mean"], color=S.INK, lw=0.9,
                    ls=(0, (3, 1.5)), label="Height, $h_F$")
        ax_tilt.plot(xf, flame["flame_tilt_mean"], color=color, lw=1.0)
        ax_dep.plot(xf, flame["flame_depth_mean"], color=color, lw=1.0)
    else:
        for ax in (ax_len, ax_tilt, ax_dep):
            ax.text(0.5, 0.5, "no flame-geometry data", ha="center",
                    va="center", transform=ax.transAxes, color=S.MUTED)
    headroom(ax_len, 0.55)
    ax_len.set_ylabel("Flame length, height (m)")
    ax_len.legend(loc="upper left", fontsize=6.5, handlelength=1.4)
    ax_tilt.set_ylim(0, 90)
    ax_tilt.set_yticks([0, 30, 60, 90])
    ax_tilt.set_ylabel("Flame tilt from vertical (°)")
    ax_dep.set_ylim(bottom=0)
    ax_dep.set_ylabel("Flame depth (m)")
    ax_ros.set_xlim(*xlim)
    ax_ros.xaxis.set_major_locator(MaxNLocator(5))

    stem = os.path.join(sim_dir, "summary")
    S.save(fig, stem)
    plt.close(fig)
    return stem


# ---------------------------------------------------------------------------
# cross-simulation figures
# ---------------------------------------------------------------------------

def _sorted(records):
    return sorted(records, key=lambda r: r.get("sim_id", ""))


def _spread_labels(ys, gap):
    """Nudge sorted label positions apart by at least `gap`."""
    order = np.argsort(ys)
    placed = np.array(ys, dtype=float)
    for k in range(1, len(order)):
        i, j = order[k], order[k - 1]
        placed[i] = max(placed[i], placed[j] + gap)
    return placed


PLOT_STYLES = ["-", (0, (4, 1.6)), (0, (1, 1.3)), (0, (5, 1.5, 1, 1.5))]


def rolling_mean(t, v, width):
    """Centred moving average of width `width` (s), ignoring NaNs."""
    out = np.full_like(v, np.nan, dtype=float)
    for i, ti in enumerate(t):
        w = np.abs(t - ti) <= width / 2
        if np.isfinite(v[w]).any() and np.isfinite(v[i]):
            out[i] = np.nanmean(v[w])
    return out


def _series_ros(out_root, rec):
    r = read_csv(os.path.join(out_root, rec["sim_id"], "ros_time.csv"))
    if not r:
        return None
    return r["x_front_m"], r["ros_m_per_s"]


def _series_flame_length(out_root, rec):
    """Slice-mean flame length, 5 s moving average, placed at the front."""
    d = os.path.join(out_root, rec["sim_id"])
    f = read_csv(os.path.join(d, "flame_geometry.csv"))
    front = read_csv(os.path.join(d, "front_trajectory.csv"))
    if not f or not front:
        return None
    width = rec.get("meta", {}).get("ros_smooth_s", 5.0)
    return (front_x_at(front, f["t_s"]),
            rolling_mean(f["t_s"], f["flame_length_mean"], width))


def comparison_lines(records, out_root, series, ylabel, stem):
    """Every run on one axis against downwind distance.

    Colour is the canopy group and line style the plot, so identity never
    rests on colour alone; the legend sits outside the data.
    """
    roi = next((r["meta"]["roi_x_m"] for r in records if r.get("meta")),
               [-10.0, 10.0])
    xlim = x_limits(next((r for r in records if r.get("scalars")), {}))
    plots = sorted({r.get("plot") for r in records if r.get("plot")})
    style = {p: PLOT_STYLES[i % len(PLOT_STYLES)] for i, p in enumerate(plots)}

    fig, ax = plt.subplots(figsize=(S.DOUBLE * 0.6, 2.7))
    ax.axvspan(roi[0], roi[1], color=S.REGION, lw=0, zorder=0)
    for rec in _sorted(records):
        xy = series(out_root, rec)
        if xy is None:
            continue
        x, v = xy
        ok = np.isfinite(x) & np.isfinite(v)
        if not ok.any():
            continue
        x, v = x[ok], v[ok]
        color = S.group_color(rec.get("group"))
        stopped = rec.get("status") != "propagated"
        ax.plot(x, v, color=color, lw=1.1, ls=style.get(rec.get("plot"), "-"),
                label=rec["sim_id"].replace("_", " ")
                + (" (stopped)" if stopped else ""))
        if stopped:
            ax.plot(x[-1], v[-1], marker="x", ms=4, mew=1.0, color=color)
    ax.set_xlim(*xlim)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), handlelength=2.6)
    ax.text(np.mean(roi), 0.99, "ROI", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=7, color=S.INK_2)
    ax.set_xlabel("Distance downwind, $x$ (m)")
    ax.set_ylabel(ylabel)
    S.save(fig, os.path.join(out_root, stem))
    plt.close(fig)


def comparison_ros(records, out_root):
    comparison_lines(records, out_root, _series_ros, "ROS (m s$^{-1}$)",
                     "comparison_ros")


def comparison_flame_length(records, out_root):
    comparison_lines(records, out_root, _series_flame_length,
                     "Flame length, $L$ (m)", "comparison_flame_length")


def comparison_distributions(records, out_root, window):
    """(a) segment ROS samples and (b) flame-length distribution per run."""
    flag = {"roi": "in_roi", "steady": "in_steady"}.get(window)
    recs = _sorted(records)
    n = len(recs)
    fig, (ax_r, ax_l) = plt.subplots(
        1, 2, figsize=(S.DOUBLE, 0.36 * n + 0.75), sharey=True,
        gridspec_kw={"wspace": 0.08})
    ypos = np.arange(n)[::-1]          # first sim at the top
    rng = np.random.default_rng(0)

    def windowed(rec, fname, col):
        s = read_csv(os.path.join(out_root, rec["sim_id"], fname))
        if not s or col not in s:
            return np.array([])
        v = s[col]
        if flag and flag in s:
            v = v[s[flag] > 0.5]
        return v[np.isfinite(v)]

    for y, rec in zip(ypos, recs):
        color = S.group_color(rec.get("group"))
        # (a) segment ROS: few samples (8 across the ROI), so show every one
        v = windowed(rec, "ros_segments.csv", "ros_m_per_s")
        if v.size:
            jit = rng.uniform(-0.13, 0.13, v.size)
            ax_r.plot(v, y + jit, "o", ms=3.2, mfc=color, mec="white",
                      mew=0.4, zorder=3)
            ax_r.plot([np.median(v)] * 2, [y - 0.28, y + 0.28], color=S.INK,
                      lw=1.3, zorder=4)
        else:
            ax_r.text(0.02, y, short_status(rec) or "no data",
                      transform=ax_r.get_yaxis_transform(), va="center",
                      fontsize=7, color=S.INK_2, style="italic")
        # (b) flame length: many (correlated) samples, so a box
        v = windowed(rec, "flame_geometry.csv", "flame_length_mean")
        if v.size >= 5:
            p5, q1, med, q3, p95 = np.percentile(v, (5, 25, 50, 75, 95))
            ax_l.plot([p5, q1], [y, y], color=color, lw=0.8)
            ax_l.plot([q3, p95], [y, y], color=color, lw=0.8)
            ax_l.add_patch(Rectangle((q1, y - 0.26), q3 - q1, 0.52,
                                     fc="white", ec=color, lw=0.9, zorder=2))
            ax_l.plot([med, med], [y - 0.26, y + 0.26], color=color, lw=1.6,
                      zorder=3)
        else:
            ax_l.text(0.02, y, short_status(rec) or "no data",
                      transform=ax_l.get_yaxis_transform(), va="center",
                      fontsize=7, color=S.INK_2, style="italic")

    ax_r.set_yticks(ypos)
    ax_r.set_yticklabels([r["sim_id"].replace("_", " ") for r in recs])
    ax_r.set_ylim(-0.6, n - 0.4)
    ax_r.tick_params(axis="y", length=0)
    ax_l.tick_params(axis="y", length=0)
    for lab, rec in zip(ax_r.get_yticklabels(), recs):
        lab.set_color(S.group_color(rec.get("group")))
    ax_r.set_xlim(left=0)
    ax_l.set_xlim(left=0)
    ax_r.set_xlabel("ROS (m s$^{-1}$)")
    ax_l.set_xlabel("Flame length, $L$ (m)")
    S.panel_label(ax_r, "a", x=-0.13)
    S.panel_label(ax_l, "b", x=-0.03)
    S.save(fig, os.path.join(out_root, "comparison_distributions"))
    plt.close(fig)


def comparison(records, out_root, window="roi"):
    records = [r for r in records if r.get("sim_id")]
    comparison_ros(records, out_root)
    comparison_flame_length(records, out_root)
    comparison_distributions(records, out_root, window)


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", nargs="?", default="postprocess",
                    help="postprocess output directory (default ./postprocess)")
    ap.add_argument("--only", nargs="+", help="per-simulation figures for "
                    "these sim_ids only; skips the comparison figures")
    ap.add_argument("--window", default="roi", choices=("roi", "steady", "full"))
    args = ap.parse_args()

    S.apply()
    dirs = sorted(d for d in glob.glob(os.path.join(args.out, "*"))
                  if os.path.isfile(os.path.join(d, "metrics.json")))
    if args.only:
        dirs = [d for d in dirs if os.path.basename(d) in args.only]
    for d in dirs:
        stem = per_sim(d)
        print(f"  {stem or os.path.basename(d) + ' (skipped)'}")
    if not args.only:
        records = [r for r in (load_record(d) for d in dirs) if r]
        comparison(records, args.out, args.window)
        print(f"  comparison figures in {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
