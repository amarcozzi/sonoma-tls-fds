#!/usr/bin/env python
"""Merge per-simulation metrics.json into cross-simulation tables and figures.

Run after the postprocess.py array job. Reads `<out>/*/metrics.json` plus the
per-simulation time-series CSVs (all small), and writes:

  metrics_summary.csv     one row per simulation, wide. The headline table:
                          sim_id, group, plot, status, then every statistic of
                          every metric for ONE analysis window (--window,
                          default 'roi'). This is the file to join against the
                          lacunarity curves.
  metrics_long.csv        tidy: sim_id, group, plot, window, metric, stat,
                          value, units. Every window, for analysis.
  metrics_dictionary.csv  column -> description + units for the wide table.
  comparison_*.png        cross-simulation figures.

Usage:
    python scripts/aggregate_metrics.py --out postprocess
    python scripts/aggregate_metrics.py --out postprocess --window steady
"""

import argparse
import csv
import glob
import json
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotstyle
from postprocess import METRICS, STAT_NAMES, WINDOWS

SCALARS = [
    ("t_ignition_s", "s", "burner ignition time detected from the HRR record"),
    ("t_enter_roi_s", "s", "time the front reaches x=-10 m"),
    ("t_exit_roi_s", "s", "time the front reaches x=+10 m"),
    ("t_enter_roi_leading_edge_s", "s", "leading-edge (not median) ROI entry time"),
    ("t_track_end_s", "s", "last tracked front arrival time"),
    ("x_front_max_m", "m", "furthest downwind x the front reached"),
    ("ros_bulk_roi_m_per_s", "m/s", "20 m / (t_exit_roi - t_enter_roi); single "
                                    "chord across the ROI, no smoothing"),
    ("x_ignition_line_m", "m", "detected ignition line; the backing fire "
                               "upwind of it is excluded from the front"),
    ("front_flattened_fraction", "", "QC: fraction of the front the monotonic "
                                     "pass had to lift. A healthy front is a "
                                     "few percent; above ~0.2 tau(x) is "
                                     "largely flat and every ROS in the row "
                                     "is untrustworthy"),
]

STAT_DESC = {
    "n": "number of finite samples",
    "mean": "mean", "std": "sample standard deviation (ddof=1)",
    "min": "minimum", "max": "maximum",
    "p02": "2nd percentile", "p05": "5th percentile", "p25": "25th percentile",
    "p50": "median", "p75": "75th percentile", "p95": "95th percentile",
    "p98": "98th percentile",
}


def load_metrics(out_root):
    records = []
    for path in sorted(glob.glob(os.path.join(out_root, "*", "metrics.json"))):
        with open(path) as fh:
            try:
                records.append(json.load(fh))
            except json.JSONDecodeError as exc:
                print(f"  skipping unreadable {path}: {exc}", file=sys.stderr)
    return records


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6g}"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def write_summary(records, path, window):
    id_cols = ["sim_id", "group", "plot", "status", "fds_completed", "window"]
    scalar_cols = [s[0] for s in SCALARS]
    stat_cols = [f"{m}_{s}" for m in METRICS for s in STAT_NAMES]
    header = id_cols + scalar_cols + stat_cols

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for rec in records:
            row = [rec.get("sim_id", ""), rec.get("group", ""),
                   rec.get("plot", ""), rec.get("status", ""),
                   fmt(rec.get("fds_completed")), window]
            scalars = rec.get("scalars", {})
            row += [fmt(scalars.get(c)) for c in scalar_cols]
            stats = rec.get("stats", {})
            for m in METRICS:
                block = stats.get(m, {}).get(window, {})
                row += [fmt(block.get(s)) for s in STAT_NAMES]
            w.writerow(row)
    return header


def write_long(records, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sim_id", "group", "plot", "status", "window", "metric",
                    "stat", "value", "units"])
        for rec in records:
            base = [rec.get("sim_id", ""), rec.get("group", ""),
                    rec.get("plot", ""), rec.get("status", "")]
            for metric, (units, _) in METRICS.items():
                for window in WINDOWS:
                    block = rec.get("stats", {}).get(metric, {}).get(window, {})
                    for stat in STAT_NAMES:
                        if stat in block:
                            w.writerow(base + [window, metric, stat,
                                               fmt(block[stat]), units])
            for name, units, _ in SCALARS:
                val = rec.get("scalars", {}).get(name)
                if val is not None:
                    w.writerow(base + ["", "scalar", name, fmt(val), units])


def write_dictionary(path, window, records):
    meta = next((r.get("meta", {}) for r in records if r.get("meta")), {})
    window_desc = {
        "roi": f"front inside the ROI, x in {meta.get('roi_x_m', [-10, 10])} m",
        "steady": f"front inside x in {meta.get('steady_x_m', [-30, 20])} m "
                  "(past the stand-edge entrance jet, short of the kill devices)",
        "full": "the whole tracked propagation record, from ignition to the "
                "last arrival",
    }
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["column", "units", "description"])
        w.writerow(["sim_id", "", "simulation identifier, <group>_<plot>"])
        w.writerow(["group", "", "TLS canopy field group (c4 / c8)"])
        w.writerow(["plot", "", "TLS plot identifier"])
        w.writerow(["status", "", "propagated | stalled_in_roi | no_roi_arrival "
                                  "| no_spread | error"])
        w.writerow(["fds_completed", "", "true if FDS wrote its .end file "
                                         "(otherwise the run was still going "
                                         "or was killed)"])
        w.writerow(["window", "", f"analysis window for every *_<stat> column: "
                                  f"{window} = {window_desc.get(window, '')}"])
        for name, units, desc in SCALARS:
            w.writerow([name, units, desc])
        for metric, (units, mdesc) in METRICS.items():
            pooled = (" Pooled over all vertical HRRPUV slices and all "
                      "in-window frames, so n is ~5x the frame count."
                      if metric.startswith("flame") else "")
            for stat in STAT_NAMES:
                w.writerow([f"{metric}_{stat}",
                            "" if stat == "n" else units,
                            f"{STAT_DESC[stat]} of {metric} ({mdesc})."
                            f"{pooled}"])
        w.writerow([])
        w.writerow(["# provenance", "", json.dumps(meta, sort_keys=True)])


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def read_series(out_root, sim_id, filename):
    path = os.path.join(out_root, sim_id, filename)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            header = fh.readline().strip().split(",")
        data = np.genfromtxt(path, delimiter=",", skip_header=1)
    except Exception:
        return None
    if data.size == 0:
        return None
    data = np.atleast_2d(data)
    return {name: data[:, i] for i, name in enumerate(header)}


def declutter_labels(ax, points, color, min_gap_frac=0.052, fontsize=7):
    """Direct-label line ends, pushed apart vertically so they never overprint.

    points: iterable of (x, y, text) at the end of each line. Labels that would
    collide are nudged along y and given a hairline leader back to the anchor.
    """
    points = [p for p in points if np.isfinite(p[0]) and np.isfinite(p[1])]
    if not points:
        return
    ax.relim()
    ax.autoscale_view()
    lo, hi = ax.get_ylim()
    gap = abs(hi - lo) * min_gap_frac
    ordered = sorted(points, key=lambda p: p[1])
    placed = [ordered[0][1]]
    for _, y, _ in ordered[1:]:
        placed.append(max(y, placed[-1] + gap))
    shift = max(0.0, placed[-1] - hi + gap)      # keep the stack inside the axes
    for (x, y, text), y_lab in zip(ordered, placed):
        y_lab -= shift
        if abs(y_lab - y) > 0.25 * gap:
            ax.plot([x, x], [y, y_lab], color=color, lw=0.5, alpha=0.45,
                    zorder=1)
        ax.annotate(text, (x, y_lab), fontsize=fontsize, color=color,
                    textcoords="offset points", xytext=(4, 0), va="center")


def comparison_figures(records, out_root, window):
    S = plotstyle
    groups = sorted({r.get("group", "") for r in records})
    win_flag = {"roi": "in_roi", "steady": "in_steady"}.get(window)

    # --- 1. front trajectories, small multiples by group ------------------
    fig, axes = plt.subplots(1, max(1, len(groups)),
                             figsize=(4.2 * max(1, len(groups)), 4.0),
                             sharex=True, sharey=True, squeeze=False)
    for ax, group in zip(axes[0], groups):
        labels = []
        for rec in records:
            if rec.get("group") != group:
                continue
            s = read_series(out_root, rec["sim_id"], "front_trajectory.csv")
            if not s:
                continue
            tcol = next((c for c in s if c.startswith("t_") and
                         not c.endswith(("_raw_s", "leading_s", "trailing_s"))),
                        None)
            if tcol is None:
                continue
            ax.plot(s[tcol], s["x_m"], lw=1.2, alpha=0.85,
                    color=S.group_color(group))
            good = np.isfinite(s[tcol])
            if good.any():
                i = np.argmax(np.where(good, s[tcol], -np.inf))
                labels.append((s[tcol][i], s["x_m"][i], rec["plot"]))
        declutter_labels(ax, labels, color=S.INK_2)
        S.roi_marker(ax, plotstyle_roi(records), axis="y")
        ax.set_title(f"{group}", color=S.group_color(group))
        ax.set_xlabel("time (s)")
    axes[0][0].set_ylabel("front position x (m)")
    fig.suptitle("Fire front trajectories (dashed: ROI bounds)", y=1.02)
    fig.savefig(os.path.join(out_root, "comparison_front_position.png"))
    plt.close(fig)

    # --- 2. ROS and flame-length distributions, one row per simulation ----
    for metric, filename, column, xlabel in (
        ("ros_seg", "ros_segments.csv", "ros_m_per_s",
         f"chord ROS over segments (m/s) — {window} window"),
        ("flame_length", "flame_geometry.csv", "flame_length_mean",
         f"flame length $L$ (m) — {window} window"),
    ):
        samples, labels, colors = [], [], []
        for rec in sorted(records, key=lambda r: r.get("sim_id", "")):
            s = read_series(out_root, rec["sim_id"], filename)
            if not s or column not in s:
                continue
            v = s[column]
            if win_flag and win_flag in s:
                v = v[s[win_flag] > 0.5]
            v = v[np.isfinite(v)]
            if v.size < 2:
                continue
            samples.append(v)
            labels.append(rec["sim_id"])
            colors.append(S.group_color(rec.get("group", "")))
        if not samples:
            continue

        order = np.argsort([np.median(v) for v in samples])
        samples = [samples[i] for i in order]
        labels = [labels[i] for i in order]
        colors = [colors[i] for i in order]

        fig, ax = plt.subplots(figsize=(6.5, 0.32 * len(samples) + 1.6))
        bp = ax.boxplot(samples, orientation="horizontal", widths=0.6,
                        showfliers=False,
                        patch_artist=True, medianprops=dict(color=S.INK, lw=1.6),
                        whiskerprops=dict(color=S.MUTED, lw=1.0),
                        capprops=dict(color=S.MUTED, lw=1.0),
                        whis=(5, 95))
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.35)
            patch.set_edgecolor(color)
            patch.set_linewidth(1.4)
        for i, v in enumerate(samples, start=1):
            ax.plot([np.percentile(v, 98)], [i], "|", ms=9, color=S.INK_2)
        ax.set_yticks(range(1, len(labels) + 1))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel(xlabel)
        ax.set_title(f"{metric} distribution by simulation "
                     f"(box: IQR, whiskers: p5–p95, tick: p98)")
        handles = [plt.Line2D([], [], color=S.group_color(g), lw=6, alpha=0.5,
                              label=g) for g in groups if g]
        if len(handles) > 1:
            ax.legend(handles=handles, loc="lower right", title=None)
        ax.grid(axis="y", visible=False)
        fig.savefig(os.path.join(out_root, f"comparison_{metric}_distribution.png"))
        plt.close(fig)

    # --- 3. ROS vs flame length ------------------------------------------
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    plotted = False
    for rec in records:
        stats = rec.get("stats", {})
        xm = stats.get("ros_seg", {}).get(window, {})
        ym = stats.get("flame_length", {}).get(window, {})
        if not xm or not ym:
            continue
        xv, yv = xm.get("mean"), ym.get("mean")
        if xv is None or yv is None or not (np.isfinite(xv) and np.isfinite(yv)):
            continue
        color = S.group_color(rec.get("group", ""))
        ax.errorbar(xv, yv,
                    xerr=[[max(0, xv - xm.get("p05", xv))],
                          [max(0, xm.get("p95", xv) - xv)]],
                    yerr=[[max(0, yv - ym.get("p05", yv))],
                          [max(0, ym.get("p95", yv) - yv)]],
                    fmt="o", ms=7, color=color, mec=S.SURFACE, mew=1.4,
                    ecolor=color, elinewidth=1.0, alpha=0.9, capsize=0)
        ax.annotate(rec["sim_id"], (xv, yv), fontsize=7, color=S.INK_2,
                    textcoords="offset points", xytext=(8, 4))
        plotted = True
    if plotted:
        ax.set_xlabel("mean chord ROS (m/s)")
        ax.set_ylabel("mean flame length $L$ (m)")
        ax.set_title(f"ROS vs flame length ({window} window; bars: p5–p95)")
        handles = [plt.Line2D([], [], color=S.group_color(g), marker="o", ls="",
                              ms=7, label=g) for g in groups if g]
        if len(handles) > 1:
            ax.legend(handles=handles, loc="best")
        fig.savefig(os.path.join(out_root, "comparison_ros_vs_flame_length.png"))
    plt.close(fig)


def plotstyle_roi(records):
    meta = next((r.get("meta", {}) for r in records if r.get("meta")), {})
    return meta.get("roi_x_m", [-10.0, 10.0])


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default="postprocess",
                    help="directory holding the per-simulation subdirectories")
    ap.add_argument("--window", choices=WINDOWS, default="roi",
                    help="analysis window used for the wide summary table")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    records = load_metrics(args.out)
    if not records:
        print(f"no metrics.json found under {args.out}", file=sys.stderr)
        return 1
    print(f"aggregating {len(records)} simulation(s) "
          f"(window: {args.window})")

    summary = os.path.join(args.out, "metrics_summary.csv")
    write_summary(records, summary, args.window)
    write_long(records, os.path.join(args.out, "metrics_long.csv"))
    write_dictionary(os.path.join(args.out, "metrics_dictionary.csv"),
                     args.window, records)

    if not args.no_figures:
        plotstyle.apply()
        try:
            comparison_figures(records, args.out, args.window)
        except Exception as exc:
            print(f"  WARNING comparison figures failed: {exc}", file=sys.stderr)

    by_status = {}
    for rec in records:
        by_status[rec.get("status", "?")] = by_status.get(rec.get("status", "?"), 0) + 1
    print("  status: " + ", ".join(f"{k}={v}" for k, v in sorted(by_status.items())))
    print(f"  wrote {summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
