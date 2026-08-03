# Postprocessing on Ceres

The raw FDS output is far too large to move, so the analysis runs next to the
data and only tables and figures come home.

```
run_postprocess.sh          driver: uv sync on the login node, then submits both jobs
  submit_postprocess.sh     Slurm array, one task per simulation  -> per-sim CSV/PNG
  submit_aggregate.sh       chained after the array               -> cross-sim tables/figures
scripts/postprocess.py      the per-simulation analysis
scripts/aggregate_metrics.py  merges metrics.json into the tables
scripts/plotstyle.py        shared figure styling
```

## Quickstart

```sh
cd ~/.../sonoma-tls-fds
./run_postprocess.sh                    # everything in identifiers.txt
./run_postprocess.sh c4_p10 c8_p25      # just these two
WINDOW=steady ./run_postprocess.sh      # summary table over the steady window
```

The driver resolves the environment once on the login node (`module load uv;
uv sync`) so 21 array tasks don't race to build `.venv` on a shared
filesystem; the array tasks then run `uv run --no-sync`. Overrides:
`SIM_ROOT` (default `simulations`), `OUT_DIR` (default `postprocess`),
`WINDOW` (default `roi`), `CONCURRENT` (default 10).

Runs that are still going are fine to process — the analysis uses whatever FDS
has flushed, and `fds_completed` records whether the `.end` file exists. A
simulation that fails outright still gets a row with `status=error`, so nothing
silently disappears from the table.

Pulling results back:

```sh
rsync -av --include='*/' --include='*.csv' --include='*.png' --exclude='*' \
      ceres:~/.../sonoma-tls-fds/postprocess/ ./postprocess/
```

## Outputs

Per simulation, in `postprocess/<sim_id>/`:

| file | contents |
|---|---|
| `metrics.json` | everything below, machine-readable; the aggregation input |
| `front_trajectory.csv` | `tau(x)` — median, leading- and trailing-edge arrival time |
| `ros_time.csv` | `t`, `x_front`, instantaneous ROS, window flags |
| `ros_segments.csv` | chord ROS per x-segment, window flags |
| `flame_geometry.csv` | flame length / height / tilt / depth, slice mean and per slice |
| `front_position.png`, `arrival_map.png`, `ros.png`, `ros_segments.png`, `flame_geometry.png` | |

Across simulations, in `postprocess/`:

| file | contents |
|---|---|
| **`metrics_summary.csv`** | **one row per simulation, wide — the table to join against the lacunarity curves** |
| `metrics_summary_{roi,steady,full}.csv` | the same table for each of the other windows |
| `metrics_long.csv` | tidy `sim_id, group, plot, window, metric, stat, value, units` |
| `metrics_dictionary.csv` | every column explained, with units and the run parameters |
| `comparison_front_position.png` | front trajectories, small multiples by group |
| `comparison_ros_seg_distribution.png` | ROS distribution per simulation (box + p98) |
| `comparison_flame_length_distribution.png` | flame-length distribution per simulation |
| `comparison_ros_vs_flame_length.png` | mean ROS vs mean flame length, p5–p95 bars |

`metrics_summary.csv` columns are `sim_id, group, plot, status, fds_completed,
window`, then the scalars, then `<metric>_<stat>` for every combination of

* metric: `ros`, `ros_seg`, `flame_length`, `flame_height`, `flame_tilt`, `flame_depth`
* stat: `n, mean, std, min, p02, p05, p25, p50, p75, p95, p98, max`

`sim_id` splits as `<group>_<plot>` (e.g. `c4_p10` -> `c4`, `p10`), matching the
`Auxiliary_Files/*.bdf` TLS bulk-density fields the runs were built from.

## What the metrics mean

**Rate of spread.** The front comes from the cumulative `FIRE ARRIVAL TIME`
boundary file on the ground (FDS UG Sec. 22.10.34), reduced across the
`|y| <= 10 m` band to `tau(x)` by the median (the mean fireline, robust to
fingering) and forced monotonic in `+x`. Two estimators:

* `ros` — d/dt of the smoothed `X_front(t)`, sampled every 0.5 s with a 5 s
  smoothing window. Dense, but successive samples are correlated by that window.
* `ros_seg` — `2.5 m / (tau(x+2.5) - tau(x))`. Sparser, but the samples are
  close to independent, so this is the defensible one to quote *spread* from.
* `ros_bulk_roi_m_per_s` — a single chord, `20 m / (t_exit - t_enter)`.

**Flame geometry**, from the five vertical HRRPUV slices at `y = -2..2`, per
Fig. 3 of Alexander (1982):

* `flame_length` **L** — straight-line distance from the *middle of the flame
  base* to the flame tip. This is the wind-tilted length that feeds
  Byram-style intensity relations, and is the "flame length" you asked for.
* `flame_height` **h_F** — vertical fuel-bed-to-tip distance. This is what the
  previous version of `postprocess.py` reported; under a tilted flame
  `h_F < L`, so the two are not interchangeable.
* `flame_tilt` — angle of the base→tip vector from vertical (deg).
* `flame_depth` — along-wind extent of the flaming zone at the base.

A cell is flaming when `HRRPUV >= min(200, 20/dx)` kW/m³ (UG Eq. 22.18) — the
same threshold FDS uses internally for `FIRE ARRIVAL TIME`, so ROS and flame
geometry share one definition of "flame". Each slice is assembled across meshes
by hand because `Slice.to_global()` fails on the stretched TRNZ grid. The
flaming field is then restricted to `|x - X_front(t)| <= 15 m` and the flame is
taken to be the **ground-connected component nearest the front**, so residual
burning back at the igniter, detached flamelets, and canopy flare-ups outside
the head-fire band cannot inflate the length.

**Analysis windows.**

| window | definition |
|---|---|
| `roi` | front inside x ∈ (−10, 10) — the region of interest |
| `steady` | front inside x ∈ (−30, 20) — past the stand-edge entrance jet, short of the kill devices at x=+30 |
| `full` | the whole tracked record, ignition to last arrival |

Flame statistics pool all five slices at all in-window frames (so `n` ≈ 5× the
frame count); ROS statistics use the single front track.

## Read the `n` columns before quoting a percentile

This matters for the ROS numbers specifically:

* **`flame_*`** are well sampled even in the ROI — 5 slices × ~90 frames ≈ 450
  samples. `p95`/`p98` are meaningful.
* **`ros_seg`** gives only **8** samples across the 20 m ROI. At that count
  `p95` and `p98` are just interpolations near the maximum and carry no
  information beyond `max`. Quote ROI ROS as `mean` / `p50` / `ros_bulk_roi`,
  and take the upper percentiles from `metrics_summary_steady.csv` (20 samples)
  or `metrics_summary_full.csv` instead.
* **`ros`** has ~80 samples in the ROI but a 5 s smoothing window, so its `std`
  understates the true variability. Use `ros_seg` for spread, `ros` for shape.

Raising `--segment` trades sample count for independence; the per-simulation
CSVs are written so any of this can be recomputed without touching Ceres again.

## Useful knobs

```sh
python scripts/postprocess.py --root simulations --out postprocess \
    --front-stat median \      # or 'min' for the leading edge
    --y-band 10 \              # half-width of the front averaging band
    --segment 2.5 \            # chord-ROS segment length (m)
    --smooth 5 \               # ROS derivative smoothing (s)
    --steady-x -30 20 \        # bounds of the 'steady' window
    --flame-window 15 \        # half-width of the head-fire search band (m)
    --no-figures               # tables only
```

## Caveat

The numerics (front tracking, both ROS estimators, the flame-geometry
extraction, the statistics and windowing) are validated against synthetic
fields with known analytic answers. The fdsreader I/O layer —
`assemble_arrival_map` and `assemble_slice` — could **not** be exercised
locally, because no FDS output was on this machine when the scripts were
written. Check the first array task's log before trusting a full sweep.
