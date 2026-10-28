upp-py
======

upp-py 用纯 Python 复现 NOAA/DTC UPP v4.1（``dtc_post_v4.1.0``）里 WRF-ARW 这一支的后处理。数值用 NumPy，读 ``wrfout``、写 NetCDF 用 xarray。它不调用 ``unipost.exe``，也不包一层 Fortran。

算法一致，依据同一套公式与常数，也依据插值和下地下外推。输出不追求和 GRIB 逐比特相同。GRIB 打包、查表量化会留下可以量出来的小数差。

.. toctree::
   :maxdepth: 2
   :caption: 开始

   getting_started/install
   getting_started/quickstart

.. toctree::
   :maxdepth: 2
   :caption: 用户指南

   user_guide/reading
   user_guide/interpolation
   user_guide/mslp
   user_guide/moisture
   user_guide/surface
   user_guide/precip
   user_guide/visibility
   user_guide/icing
   user_guide/aerosol
   user_guide/reflectivity
   user_guide/cape
   user_guide/helicity
   user_guide/config
   user_guide/output

.. toctree::
   :maxdepth: 2
   :caption: API

   api

.. toctree::
   :maxdepth: 1
   :caption: 图集

   auto_examples/index

.. toctree::
   :maxdepth: 2
   :caption: 验证

   validation

源码对照的是 `EMC_post` 标签 ``dtc_post_v4.1.0`` 的 ``sorc/ncep_post.fd``。每个用户指南章节都写明对应的 Fortran 子程序。
