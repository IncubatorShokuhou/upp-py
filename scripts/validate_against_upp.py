#!/usr/bin/env python3
"""Compare upp-py to a unipost.exe GRIB2 file.

Reads the reference with cfgrib (optional extra ``grib``). Writes
``docs/validation.md`` and ``docs/validation_generated.rst``.

Example::

    python scripts/validate_against_upp.py \\
        --wrfout /path/wrfout_d01_2018-07-23_21_00_00 \\
        --grib /path/WRFPRS_d01.00
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.stats import pearsonr

from upp_py.post import post_wrfout


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3 or np.std(a[mask]) == 0 or np.std(b[mask]) == 0:
        return float("nan")
    return float(pearsonr(a[mask], b[mask])[0])


def _metrics(model: np.ndarray, ref: np.ndarray) -> dict:
    a = np.asarray(model, dtype=np.float64)
    b = np.asarray(ref, dtype=np.float64)
    d = a - b
    finite = np.isfinite(d)
    dd = d[finite]
    return {
        "bias": float(np.mean(dd)),
        "rmse": float(np.sqrt(np.mean(dd**2))),
        "maxabs": float(np.max(np.abs(dd))),
        "corr": _corr(a, b),
    }


def _grib(path: str, **keys):
    ds = xr.open_dataset(path, engine="cfgrib", backend_kwargs={"filter_by_keys": keys})
    name = [n for n in ds.data_vars if n not in ("crs",)][0]
    return np.asarray(ds[name].values)


def collect(wrfout: str, grib: str) -> list[dict]:
    ds = post_wrfout(wrfout)
    rows = []

    def add(name, model, ref, unit):
        m = _metrics(model, ref)
        m.update(field=name, unit=unit)
        rows.append(m)
        print(f"{name:22s} bias {m['bias']:10.4g} rmse {m['rmse']:10.4g} max {m['maxabs']:10.4g} corr {m['corr']:.5f}")

    add("PRMSL", ds["prmsl"].values[0], _grib(grib, shortName="prmsl"), "Pa")
    add("MSLET", ds["mslet"].values[0], _grib(grib, shortName="mslet"), "Pa")
    add("TMP 2 m", ds["t2m"].values[0], _grib(grib, shortName="2t"), "K")
    add("DPT 2 m", ds["dpt2m"].values[0], _grib(grib, shortName="2d"), "K")
    add("RH 2 m", ds["rh2m"].values[0], _grib(grib, shortName="2r"), "%")
    add("UGRD 10 m", ds["u10"].values[0], _grib(grib, shortName="10u"), "m/s")
    add("VGRD 10 m", ds["v10"].values[0], _grib(grib, shortName="10v"), "m/s")
    for p in (850, 500, 250):
        add(
            f"HGT {p} hPa",
            ds["gh"].sel(isobaric=p * 100).values[0],
            _grib(grib, shortName="gh", typeOfLevel="isobaricInhPa", level=p),
            "m",
        )
        add(
            f"TMP {p} hPa",
            ds["t"].sel(isobaric=p * 100).values[0],
            _grib(grib, shortName="t", typeOfLevel="isobaricInhPa", level=p),
            "K",
        )
        add(
            f"RH {p} hPa",
            ds["r"].sel(isobaric=p * 100).values[0],
            _grib(grib, shortName="r", typeOfLevel="isobaricInhPa", level=p),
            "%",
        )
    add("REFC", ds["refc"].values[0], _grib(grib, shortName="refc"), "dBZ")
    add(
        "REFD 1 km",
        ds["refd"].sel(height_agl=1000).values[0],
        _grib(grib, shortName="refd", typeOfLevel="heightAboveGround", level=1000),
        "dBZ",
    )
    add("CAPE surface", ds["cape"].sel(parcel="surface").values[0], _grib(grib, shortName="cape", typeOfLevel="surface"), "J/kg")
    add("CIN surface", ds["cin"].sel(parcel="surface").values[0], _grib(grib, shortName="cin", typeOfLevel="surface"), "J/kg")
    for label, top in (("mixed_layer", 90), ("best", 180), ("most_unstable", 255)):
        add(
            f"CAPE {label}",
            ds["cape"].sel(parcel=label).values[0],
            _grib(grib, shortName="cape", typeOfLevel="pressureFromGroundLayer", topLevel=top * 100, bottomLevel=0),
            "J/kg",
        )
        add(
            f"CIN {label}",
            ds["cin"].sel(parcel=label).values[0],
            _grib(grib, shortName="cin", typeOfLevel="pressureFromGroundLayer", topLevel=top * 100, bottomLevel=0),
            "J/kg",
        )
    hl = xr.open_dataset(grib, engine="cfgrib", backend_kwargs={"filter_by_keys": {"shortName": "hlcy"}})
    # cfgrib orders the layer coordinate ascending: 1000 m then 3000 m.
    add("HLCY 0-1 km", ds["hlcy"].sel(hel_depth=1000).values[0], hl["hlcy"].values[0], "m2/s2")
    add("HLCY 0-3 km", ds["hlcy"].sel(hel_depth=3000).values[0], hl["hlcy"].values[1], "m2/s2")
    add("PWAT", ds["pwat"].values[0], _grib(grib, shortName="pwat"), "kg/m2")
    add("APCP", ds["apcp"].values[0], _grib(grib, shortName="tp"), "kg/m2")
    return rows


def write_reports(rows: list[dict], md_path: Path, rst_path: Path) -> None:
    header = "| 场 | 单位 | bias | RMSE | max abs | 相关系数 |"
    sep = "| --- | --- | ---: | ---: | ---: | ---: |"
    lines = [header, sep]
    rst = [
        ".. list-table::",
        "   :header-rows: 1",
        "   :widths: 22 10 14 14 14 12",
        "",
        "   * - 场",
        "     - 单位",
        "     - bias",
        "     - RMSE",
        "     - max abs",
        "     - 相关系数",
    ]
    for r in rows:
        corr = "nan" if not np.isfinite(r["corr"]) else f"{r['corr']:.5f}"
        lines.append(
            f"| {r['field']} | {r['unit']} | {r['bias']:.4g} | {r['rmse']:.4g} | {r['maxabs']:.4g} | {corr} |"
        )
        rst.append(f"   * - {r['field']}")
        rst.append(f"     - {r['unit']}")
        rst.append(f"     - {r['bias']:.4g}")
        rst.append(f"     - {r['rmse']:.4g}")
        rst.append(f"     - {r['maxabs']:.4g}")
        rst.append(f"     - {corr}")
    md_path.write_text(
        "\n".join(
            [
                "# upp-py 与 unipost.exe 的对比",
                "",
                "Carr Fire `wrfout_d01_2018-07-23_21_00_00`（WRF-ARW 4.7.1，199×199×79，3 km Lambert）",
                "对 UPP v4.1 `unipost.exe` 的 `WRFPRS_d01.00`（549 条 GRIB2，起报 2018072321，fhr 00）。",
                "差值为 upp-py 减参考场。GRIB 由 cfgrib 读出，打包误差会留在小数位上。",
                "",
                *lines,
                "",
                "## 读表时需要知道的事",
                "",
                "- **PRMSL** 在这份 `postxconfig-NT-WRF.txt` 里是可用场 105，对应 `NGMSLP.f` 的 Shuell 订正，不是膜法。",
                "  膜法是可用场 23，短名 **MSLET**（`SLP_new.f` 的 `MEMSLP`）。两套都做了，并分别和对应的 GRIB 记录比。",
                "- 2 m 温度、露点、相对湿度来自 `TH2`/`Q2` 和 `SURFCE.f` 的遮蔽层气压；10 m 风直接用 `U10`/`V10`，不旋转到地球坐标。",
                "  这些场和 250/500 hPa 的高度、温度几乎只剩 GRIB 打包误差。",
                "- 850 hPa 有少量点落在地下外推里，温度最大差约 3 K，高度最大差约 11 m，相关系数仍高于 0.999。",
                "- **CAPE/CIN**：地面、混合层（0–90 hPa）和最不稳定（搜索 300 hPa，GRIB 层标成 255 hPa）与 `CALCAPE.f` 很接近。",
                "  最佳气块（六个 30 hPa 层里 Bolton θe 最大的一层）和 unipost.exe 的剩余差别是预期的，也是有意保留的。upp-py 直接解湿绝热和 Bolton θe，UPP 的 `CALCAPE` 用查表。θe 极接近时两边可能选到相邻一层。约 0.2% 的格点选层不同，RMSE 约 100 J/kg，99 分位绝对差约 20 J/kg。",
                "  地面、混合层和最不稳定 CIN 的平均偏差约 1–2 J/kg，最大绝对差约 680–750 J/kg。离群点上直接求解没有正浮力层（约 350–400 hPa 只比环境冷 0.01–0.04 K），按 `CALCAPE.f` 把 CIN 留在 0。参考场有约 1–9 J/kg 的 CAPE 和数百 J/kg 的负 CIN，是查表把该层报成略暖后，从抬升凝结高度积到这个假平衡高度。浮力符号以直接求解为准，求解器不改为对齐这几百 J/kg。例子：(173, 51) 42.688N 124.422W 地面 0 对 -748；(151, 55) 42.076N 124.253W 0 对 -746；(92, 53) 40.432N 124.282W 混合层 0 对 -685。",
                "- **HLCY** 0–1 km 和 0–3 km 与 `CALHEL.f`（Bunkers）一致，差在打包误差。",
                "- **PWAT** 按 `CALPW.f`（`IDECID<=1`，只积分水汽）计算。",
                "- **APCP** 在这个时次是模式起始时刻，`RAINC+RAINNC` 为 0，和 GRIB 的 0–0 累积一致。",
                "- **REFC**：NCAR + Thompson（`mp_physics=8`）写出的是模式变量 `REFL_10CM` 的柱最大（`MDLFLD.f`）。这个 wrfout 里该变量全是 0，所以两边都是 0，相关系数没有定义。Stoelinga 诊断仍会算（unipost 日志里的 `dbze_max` 约 1 dBZ），但不会写进这条 REFC。",
                "- **REFD 1 km**：Stoelinga 公式里对流项不低于 0.1 mm6/m3，晴空裸算约 -10 dBZ。写出时把晴空和低于 0 dBZ 的值记为 0，这是有意的。参考 GRIB 的 1 km 场整场为 0，两边因此一致。",
                "- 整场为常数时（REFC、APCP）相关系数没有定义，表里写成 nan。",
                "",
            ]
        ),
        encoding="utf-8",
    )
    rst_path.write_text("\n".join(rst) + "\n", encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--wrfout", required=True)
    p.add_argument("--grib", required=True)
    p.add_argument("--md", default="docs/validation.md")
    p.add_argument("--rst", default="docs/validation_generated.rst")
    args = p.parse_args()
    rows = collect(args.wrfout, args.grib)
    write_reports(rows, Path(args.md), Path(args.rst))
    print(f"wrote {args.md} and {args.rst}")


if __name__ == "__main__":
    main()
