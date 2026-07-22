# Analysis: why c4_p25 / c8_p25 do not sustain spread (2026-07-22)

Diagnostics of the three example runs (`c4_p10`, `c4_p25`, `c8_p25`) at
`SPEED=3.75` wind. The three input decks are **byte-identical except for the
tree `BULK_DENSITY_FILE`** — "p10/p25" are different TLS plots (canopy
structures), not moisture levels; litter is 10 % MC in all three.

## Key findings

1. **The canopy never burns.** < 2 % of canopy mass consumed, < 0.03 % ever
   exceeds 250 °C in any run. These are pure surface fires under a canopy
   that only removes energy and momentum.
2. **The p25 stands stall because of canopy structure.** They carry
   ~1.5–2x more canopy mass below 15 m over the fetch (c4_p25 mean 87 kg/m²
   vs c4_p10 58; c4_p10 is tall-crowned and patchy, with an open corridor at
   x ≈ −30…−10). Two mechanisms follow:
   - *Drag*: mid-fetch sub-canopy wind ~0.5–1.5 m/s under the p25 canopies
     vs ~2–2.5 m/s under c4_p10 (the empirical sustainment threshold on this
     litter bed).
   - *Moisture sink*: 125 %-MC foliage; c4_p25 evaporates 163 kg/s of water
     on a 1.6 GW fire vs c4_p10's 145 kg/s on 2.8 GW (~2x the drying burden
     per unit HRR).
3. **All runs flash to 1.6–2.8 GW at t ≈ 40 s on the ignition pulse**, then
   c4_p10 stabilizes at ~1 GW while both p25 runs decay monotonically to
   extinction at ~100–115 s. There is no plateau: they were never in a
   self-sustained regime, so a longer ignition pulse would not change the
   outcome.
4. **Wind-field checks** (`wind_checks/`):
   - The inflow/stand-edge adjustment jet (6–11 m/s, full height) dies by
     x ≈ −50; vertical profiles at x = −40/−20/0 collapse → the ROI sits in
     canopy-equilibrated flow. No domain change needed.
   - The mean field is developed by t ≈ 15–20 s → the 30 s spin-up before
     ignition is adequate. (Statistics, not development, limit 14-s
     averages: canopy gust cycles are 10–30 s, so short-window means carry
     ±0.3–0.5 m/s.)
   - The near-zero ROI sub-canopy wind is *real equilibrium*, verified
     against c4_p25's fire-free tail (pre-fire 0.61 m/s vs fully developed
     0.69 m/s over t = 145–320 s). The ROI is the weakest wind environment
     on the fetch: 60–80 m into the stand, with local canopy clusters at
     x ≈ 0–10 in all three plots.

## Figures

`extinction_analysis/` (from `scripts/extinction_analysis.py`):
- `canopy_structure.png` — initial canopy load along x and vertical profile
- `subcanopy_wind.png` — pre-fire mean U below 3 m vs x (note: doubled
  lines/diagonal segments are a stitching artifact of the y=0 slice lying on
  a mesh interface; the signal is real)
- `canopy_heating.png` — % of canopy mass drying (90–110 °C) / pyrolyzing (>250 °C)
- `canopy_consumption.png` — canopy mass over time and consumption along x
- `hrr_water.png` — total HRR and water-vapor MLR (drying heat sink)
- `summary.txt` — per-run numbers

`wind_checks/` (from `scripts/wind_checks.py` and `scripts/roi_wind_check.py`):
- `u_cross_sections.png` — pre-fire mean u(x,z); entrance jet and canopy
  shear layer
- `u_profiles_vs_x.png` — vertical profile self-similarity along the fetch
- `stationarity.png` — pre-ignition equilibration of the sub-canopy flow
- `roi_wind_equilibrium.png` — ROI wind over the full record; pre-fire vs
  fire-free developed state (also shows fire indraft reversal in c4_p10)

## Regenerating

Scripts read simulation output from the external drive (paths hardcoded at
the top of each script: `/Volumes/T7 Shield/sonoma-tls-fds/simulations`) and
use the repo venv (fdsreader 1.11.9):

```sh
.venv/bin/python scripts/extinction_analysis.py   # reads ~15 GB of prt5, few minutes
.venv/bin/python scripts/wind_checks.py
.venv/bin/python scripts/roi_wind_check.py
```

`extinction_analysis.py` contains a minimal seek-based PRT5 reader (frame
subsampling, centre-row meshes only) because fdsreader would load all 41 GB
of c4_p25's particle data; it is also the workaround for two fdsreader
1.11.9 limitations (Slice.to_global() fails on the stretched TRNZ grid;
Boundary.data keeps only one patch per orientation).

General ROS/flame-height postprocessing (arrival-time front tracking, flame
height from HRRPUV slices) lives in `scripts/postprocess.py`.
