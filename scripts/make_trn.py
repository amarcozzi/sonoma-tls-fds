"""Generate a vertical mesh-stretching file (TRNZ) for the Sonoma FDS template.

The vertical mesh is built with the same design used throughout this project:

    * a near-surface layer of `n_uniform` cells of constant height `dz0`
      (to resolve the litter bed, canopy particles, and near-ground fire), then
    * a geometrically growing region whose ratio is solved so the column
      reaches `ztop` in exactly `nz` cells.

FDS applies the map via `TRNZ` control points: `CC` is the cell-face position on
the *uniform* (computational) grid, `PC` is the physical height it maps to. The
CC list therefore spans [0, ztop] with spacing ztop/nz -- i.e. the file is tied
to BOTH the domain height and the vertical cell count. Change either and you
must regenerate this file, then keep `IJK`'s 3rd value and the mesh `XB` z-range
in `template.fds` consistent with `--nz` and `--ztop`.

Examples
--------
Reproduce the current mesh (0-40 m, 50 cells, 0.125 m ground cell):
    python scripts/make_trn.py --ztop 40 --nz 50 --dz0 0.125 --n-uniform 4 \
        --output Auxiliary_Files/TRN.fds

Taller 60 m domain, keeping ground resolution, adding cells aloft:
    python scripts/make_trn.py --ztop 60 --nz 56 --output Auxiliary_Files/TRN_60m.fds
"""
import argparse
import sys
from pathlib import Path


def solve_ratio(nz, ztop, dz0, n_uniform):
    """Bisection-solve the geometric growth ratio so the column reaches `ztop`.

    The first `n_uniform` cells are `dz0` tall; the remaining `n_geo` cells grow
    geometrically (first geometric cell = dz0 * ratio). Returns the ratio.
    """
    n_geo = nz - n_uniform
    if n_geo < 1:
        raise ValueError("nz must be greater than n_uniform")
    remaining = ztop - n_uniform * dz0
    if remaining <= n_geo * dz0:
        raise ValueError(
            f"Domain too short: {n_uniform} uniform cells of {dz0} m already "
            f"fill {n_uniform * dz0} m of a {ztop} m domain with no room to grow."
        )
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(500):
        r = 0.5 * (lo + hi)
        total = dz0 * r * (r ** n_geo - 1) / (r - 1)
        if total > remaining:
            hi = r
        else:
            lo = r
    return 0.5 * (lo + hi)


def build_cell_heights(nz, ztop, dz0, n_uniform):
    """Return the list of `nz` physical cell heights summing exactly to `ztop`."""
    ratio = solve_ratio(nz, ztop, dz0, n_uniform)
    dz = [dz0] * n_uniform + [dz0 * ratio ** k for k in range(1, nz - n_uniform + 1)]
    # Correct the tiny bisection residual so the column lands exactly on ztop.
    scale = ztop / sum(dz)
    dz = [d * scale for d in dz]
    return dz, ratio * scale


def trnz_lines(nz, ztop, dz0, n_uniform, trn_id):
    """Build the TRNZ control-point lines for the stretched column."""
    dz, ratio = build_cell_heights(nz, ztop, dz0, n_uniform)
    faces = [0.0]
    for d in dz:
        faces.append(faces[-1] + d)
    dz_comp = ztop / nz  # uniform computational face spacing
    # Interior faces only (indices 1..nz-1); endpoints 0 and ztop are implicit.
    lines = [
        f"&TRNZ ID='{trn_id}', CC={i * dz_comp:19.15f}, PC={faces[i]:19.15f} /"
        for i in range(1, nz)
    ]
    return lines, dz, ratio


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Generate a TRNZ vertical mesh-stretching file for FDS.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--ztop", type=float, default=40.0, help="Domain top height (m).")
    p.add_argument("--nz", type=int, default=50, help="Number of vertical cells.")
    p.add_argument("--dz0", type=float, default=0.125, help="Ground cell height (m).")
    p.add_argument("--n-uniform", type=int, default=4,
                   help="Number of uniform near-surface cells of height dz0.")
    p.add_argument("--id", dest="trn_id", default="TRNZ", help="TRNZ ID string.")
    p.add_argument("--output", type=Path, default=None,
                   help="Output path (default: Auxiliary_Files/TRN.fds relative to repo root).")
    args = p.parse_args(argv)

    if args.output is None:
        repo_root = Path(__file__).resolve().parent.parent
        args.output = repo_root / "Auxiliary_Files" / "TRN.fds"

    try:
        lines, dz, ratio = trnz_lines(
            args.nz, args.ztop, args.dz0, args.n_uniform, args.trn_id
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")

    print(f"Wrote {args.output}  ({len(lines)} TRNZ lines)")
    print(f"  domain      : 0 - {args.ztop:g} m")
    print(f"  cells (nz)  : {args.nz}   ({args.n_uniform} uniform + "
          f"{args.nz - args.n_uniform} geometric)")
    print(f"  ground cell : {dz[0]:.4f} m")
    print(f"  top cell    : {dz[-1]:.4f} m")
    print(f"  growth ratio: {ratio:.4f}")
    print()
    print("  Remember to set template.fds MESH consistently:")
    print(f"    IJK = 40,40,{args.nz}   and   XB z-range = 0,{args.ztop:g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
