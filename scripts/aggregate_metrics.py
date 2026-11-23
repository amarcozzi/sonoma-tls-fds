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
  comparison_*.{pdf,png}  cross-simulation figures (drawn by figures.py).

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


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figures
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
                                    "segment across the ROI, no smoothing"),
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
            figures.comparison(records, args.out, args.window)
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
