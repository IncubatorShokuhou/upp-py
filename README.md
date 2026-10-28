# upp-py

WRF-ARW 后处理的纯 Python 实现，算法对齐 NOAA/DTC UPP v4.1（`EMC_post` 标签 `dtc_post_v4.1.0`，`sorc/ncep_post.fd`）。数值用 NumPy，读写 NetCDF 用 xarray。不调用 `unipost.exe`，也没有 Fortran/C 扩展。

不追求和 GRIB 逐比特相同。公式和常数按 Fortran 来，对数气压插值与地下外推也是。Shuell / 膜法海平面气压同样照 Fortran，相对湿度和露点、CAPE 气块以及 Bunkers 螺旋度也是。

**English.** upp-py is a pure-Python WRF-ARW post processor. Numerics follow UPP v4.1 (not bit-for-bit). See `docs/` for the Sphinx site (`cd docs && make html`) and `docs/validation.md` for the Carr Fire comparison against `unipost.exe`.

## 安装

```bash
pip install -e ".[dev]"
# 读参考 GRIB，以及写出 GRIB2。eccodes 的 wheel 自带 C 库。
pip install -e ".[grib]"
```

## 用法

```bash
upp-py wrfout_d01_2018-07-23_21_00_00 -o post.nc
upp-py wrfout_d01_2018-07-23_21_00_00 -o post.nc --config configs/example.yaml
upp-py wrfout_d01 -o post.nc --grib2 post.grib2 --times 0,1,2 --accum-hours 1
```

```python
from upp_py import load_config, post_wrfout, write_grib2, write_netcdf

ds = post_wrfout("wrfout_d01_2018-07-23_21_00_00", load_config("configs/example.yaml"))
write_netcdf(ds, "post.nc")
write_grib2(ds, "post.grib2")
```

## 做了什么

| 产品 | Fortran |
| --- | --- |
| 读 wrfout 与去交错，含 `P_HYD`、位温转温度和地面气压的水汽订正 | `INITPOST.F`, `getVariable.f` |
| 气压层 `ln p` 插值与地下外推 | `MDL2P.f` |
| `PRMSL`（Shuell）和 `MSLET`（膜法） | `NGMSLP.f`, `SLP_new.f` |
| 相对湿度、露点与可降水量 | `CALRH.f`, `CALDWP.f`, `DEWPOINT.f`, `CALPW.f` |
| 2 m 温湿、10 m 风 | `SURFCE.f` |
| 累积降水、`PREC_ACC_*` 桶、`bucket_mm` 计数器 | `INITPOST.F`, `SURFCE.f` |
| 能见度 `VIS`、结冰类别 `TIPD` | `CALVIS.f`, `AVIATION.f` `CALICING` |
| 气溶胶数浓度与烟、550 nm 光学厚度（文件里有才写） | `INITPOST.F`, `MDLFLD.f` |
| Thompson 反射率；NCAR 的合成反射率取 `REFL_10CM` 柱最大 | `MDLFLD.f`, `PARAMR.f` |
| 地面 / 混合层 / 最佳 / 最不稳定 CAPE 与 CIN，LCL、LFC | `CALCAPE.f`, `BNDLYR.f`, `CLDRAD.f`, `MISCLN.f` |
| 0–1 km、0–3 km 风暴相对螺旋度 | `CALHEL.f` |
| 低中高云量 | `INITPOST.F`（`COMMCODE`） |

## 文档

```bash
pip install -e ".[docs]"
cd docs && make html
```

用户指南是中文。API 文档字符串是英文。图集默认不读 wrfout；要画 Carr Fire，先 `export UPP_PY_WRFOUT=/path/to/wrfout`。

## 验证

Carr Fire（2018-07-23 21:00，199×199×79，3 km）对同一文件上的 `unipost.exe` 输出。完整表格在 [`docs/validation.md`](docs/validation.md)。

| 场 | RMSE | 相关系数 |
| --- | ---: | ---: |
| PRMSL（Shuell） | 0.013 Pa | 1.000 |
| MSLET（膜法） | 4.0 Pa | 1.000 |
| HGT/TMP 850 hPa | 0.002 m / 0.029 K | 1.000 |
| TMP 2 m | 0.005 K | 1.000 |
| HGT/TMP 500 hPa | 0.002 m / 0.003 K | 1.000 |
| CAPE 地面 | 11 J/kg | 1.000 |
| HLCY 0–3 km | 0.03 m²/s² | 1.000 |
| PWAT | 0.03 kg/m² | 1.000 |

参考 GRIB 的合成反射率是 0，因为这个 wrfout 的 `REFL_10CM` 全是 0，而 NCAR+Thompson 的写出路径用的就是它。1 km 反射率的公式下限裸算约 −10 dBZ，写出时晴空记为 0，和参考场一致。最佳 CAPE 有约 0.2% 的格点因查表和直接求解选层不同，RMSE 约 100 J/kg。

亮温 `BRTMP` 在独立包 `upp-crtm` 里算。没装时请求该场只记 `skipped`，本仓库不带 CRTM 源码或系数。`CALVIS_GSD` 和 NCAR GFIP 不算。能见度用 `CALVIS`，结冰类别用 `CALICING`。细节在 `docs/`。

```bash
python scripts/validate_against_upp.py --wrfout wrfout_d01 --grib WRFPRS_d01.00
```

wrfout 和参考 GRIB 都不要提交进仓库。
