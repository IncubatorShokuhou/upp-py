"""Command line: ``upp-py wrfout -o out.nc --config fields.yaml``."""

from __future__ import annotations

import argparse
import sys

from upp_py import __version__
from upp_py.config import dump_example, load_config
from upp_py.grib2 import write_grib2
from upp_py.post import post_wrfout, write_netcdf


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="upp-py",
        description="Post a WRF-ARW wrfout NetCDF file with UPP v4.1 algorithms.",
    )
    parser.add_argument("wrfout", nargs="?", help="path to wrfout_d0*")
    parser.add_argument("-o", "--output", help="output NetCDF path")
    parser.add_argument("--grib2", help="also write GRIB2 to this path (needs the eccodes extra)")
    parser.add_argument("--config", help="YAML field list (default: standard core set)")
    parser.add_argument("--time", type=int, default=0, help="Time index in the wrfout file")
    parser.add_argument(
        "--times",
        help="Comma-separated time indices. Overrides --time. APCP at each time is the run total.",
    )
    parser.add_argument(
        "--accum-hours",
        type=float,
        default=None,
        help="Also write apcp_accum over this many hours, ending at each selected time.",
    )
    parser.add_argument("--write-example-config", metavar="PATH", help="write an example YAML and exit")
    parser.add_argument("--version", action="version", version=f"upp-py {__version__}")
    args = parser.parse_args(argv)
    if args.write_example_config:
        dump_example(args.write_example_config)
        return 0
    if not args.wrfout or not args.output:
        parser.error("wrfout and -o/--output are required")
    cfg = load_config(args.config)
    indices = None
    if args.times:
        indices = tuple(int(part) for part in args.times.split(",") if part.strip() != "")
    ds = post_wrfout(
        args.wrfout,
        cfg,
        time_index=args.time,
        time_indices=indices,
        accum_hours=args.accum_hours,
    )
    write_netcdf(ds, args.output)
    print(f"wrote {args.output} ({len(ds.data_vars)} variables)", file=sys.stderr)
    if args.grib2:
        skipped = write_grib2(ds, args.grib2)
        print(f"wrote {args.grib2}", file=sys.stderr)
        for line in skipped:
            print(f"grib2 skip: {line}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
