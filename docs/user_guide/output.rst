输出
====

默认写 NetCDF，变量带 CF 风格的长名和单位，另外有 ``grib2_shortName``。水平格点挂在 Lambert（或极射、墨卡托）的 ``crs`` 上，同时给出二维 ``latitude`` 和 ``longitude``，所以不依赖投影库也能画图。

.. list-table::
   :header-rows: 1

   * - 短名
     - 数据集里的名字
     - 单位
   * - PRMSL
     - ``prmsl``
     - Pa
   * - MSLET
     - ``mslet``
     - Pa
   * - HGT / TMP / DPT / RH / SPFH / UGRD / VGRD / VVEL / ABSV
     - ``gh`` ``t`` ``dpt`` ``r`` ``q`` ``u`` ``v`` ``w`` ``absv``
     - m, K, %, kg/kg, m/s, Pa/s, 1/s
   * - TMP, DPT, RH, UGRD, VGRD 的 2 m 或 10 m
     - ``t2m`` ``dpt2m`` ``rh2m`` ``u10`` ``v10``
     - K, %, m/s
   * - REFC, REFD
     - ``refc`` ``refd``
     - dBZ
   * - CAPE, CIN
     - ``cape`` ``cin``
     - J/kg
   * - HLCY
     - ``hlcy``
     - m\ :sup:`2` s\ :sup:`−2`
   * - PWAT, APCP, PRATE
     - ``pwat`` ``apcp`` ``prate``
     - kg/m2, kg/m2, kg/m2/s

云量 ``LCDC`` / ``MCDC`` / ``HCDC`` / ``TCDC`` 来自 ``INITPOST.F`` 的社区代码分支：``qc+qi`` 超过 ``10^{-7}`` kg/kg 算有云，再取 30 km 方框里的云点比例。低云是 642 hPa 以下的最大比例，中云到 350 hPa，高云到 150 hPa。总云量是这三段的最大。单位是百分数。

气压坐标 ``isobaric`` 的单位是 Pa，向下为正。时间是 wrfout 的有效时间，预报参考时间写在数据集属性 ``forecast_reference_time``。

GRIB2 不在默认依赖里。要对参考场做对比时安装 ``cfgrib``（``pip install -e ".[grib]"``），用 ``scripts/validate_against_upp.py``。upp-py 自己不写 GRIB。

.. code-block:: python

   import xarray as xr
   ds = xr.open_dataset("post.nc")
   ds["t"].sel(isobaric=50000).plot()

和 UPP 的差别
-------------

UPP 的 GRIB 用 ``params_grib2_tbl_new`` 里的十进制比例尺打包，温度常常只保留约 0.01 K。NetCDF 用 float32 和 zlib，保留诊断值本身。和 GRIB 相减时，0.01 量级的差来自打包，不来自公式。
