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
   * - APCP 的 PREC_ACC 桶、窗口累积
     - ``apcp_bucket`` ``apcp_accum``
     - kg/m2
   * - VIS, TIPD
     - ``vis`` ``tipd``
     - m, 0 或 1

云量 ``LCDC`` / ``MCDC`` / ``HCDC`` / ``TCDC`` 来自 ``INITPOST.F`` 的社区代码分支：``qc+qi`` 超过 ``10^{-7}`` kg/kg 算有云，再取 30 km 方框里的云点比例。低云是 642 hPa 以下的最大比例，中云到 350 hPa，高云到 150 hPa。总云量是这三段的最大。单位是百分数。

气压坐标 ``isobaric`` 的单位是 Pa，向下为正。时间是 wrfout 的有效时间，预报参考时间写在数据集属性 ``forecast_reference_time``。

GRIB2
-----

写出用 ECMWF 的 ``eccodes`` Python 绑定，和读参考场是同一个可选依赖：

.. code-block:: bash

   pip install -e ".[grib]"
   upp-py wrfout_d01 -o post.nc --grib2 post.grib2 --config configs/example.yaml

.. code-block:: python

   from upp_py import post_wrfout, write_grib2

   ds = post_wrfout("wrfout_d01")
   write_grib2(ds, "post.grib2")

中心设为 NCEP（``kwbc``）。短名走 eccodes 的 kwbc 表，例如温度是 ``t`` / ``2t``，累积降水是 ``tp``，能见度是 ``vis``，结冰类别是 ``tipd``。网格是规则经纬或 Lambert，经度按 0–360 写。打包是 grid_simple、24 bit，不是 ``params_grib2_tbl_new`` 的十进制比例，也不和 ``unipost.exe`` 逐比特相同。

模式层上的气溶胶数浓度、以及 kwbc 表里没有的短名，留在 NetCDF 里。``write_grib2`` 的返回值是没写进 GRIB 的原因。``-o`` 仍然是 NetCDF。

读参考 GRIB 仍是 ``pip install -e ".[grib]"`` 之后的 ``scripts/validate_against_upp.py``。

.. code-block:: python

   import xarray as xr
   ds = xr.open_dataset("post.nc")
   ds["t"].sel(isobaric=50000).plot()

和 UPP 的差别
-------------

UPP 的 GRIB 用 ``params_grib2_tbl_new`` 里的十进制比例尺打包，温度常常只保留约 0.01 K。NetCDF 用 float32 和 zlib，保留诊断值本身。和 GRIB 相减时，0.01 量级的差来自打包，不来自公式。
