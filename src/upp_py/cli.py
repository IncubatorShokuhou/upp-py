"""Command line: ``upp-py wrfout -o out.nc --config fields.yaml``."""

from __future__ import annotations

import argparse
import sys

from upp_py import __version__
from upp_py.config import dump_example, load_config
from upp_py.post import post_wrfout, write_netcdf


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="upp-py",
        description="Post a WRF-ARW wrfout NetCDF file with UPP v4.1 algorithms.",
    )
    parser.add_argument("wrfout", nargs="?", help="path to wrfout_d0*")
    parser.add_argument("-o", "--output", help="output NetCDF path")
    parser.add_argument("--config", help="YAML field list (default: standard core set)")
    parser.add_argument("--time", type=int, default=0, help="Time index in the wrfout file")
    parser.add_argument("--write-example-config", metavar="PATH", help="write an example YAML and exit")
    parser.add_argument("--version", action="version", version=f"upp-py {__version__}")
    args = parser.parse_args(argv)
    if args.write_example_config:
        dump_example(args.write_example_config)
        return 0
    if not args.wrfout or not args.output:
        parser.error("wrfout and -o/--output are required")
    cfg = load_config(args.config)
    ds = post_wrfout(args.wrfout, cfg, time_index=args.time)
    write_netcdf(ds, args.output)
    print(f"wrote {args.output} ({len(ds.data_vars)} variables)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
