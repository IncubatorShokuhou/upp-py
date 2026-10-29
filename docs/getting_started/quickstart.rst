快速开始
========

对一个 WRF-ARW 的 ``wrfout`` 做后处理，写出 NetCDF：

.. code-block:: bash

   upp-py wrfout_d01_2018-07-23_21_00_00 -o post.nc

不给 ``--config`` 时，标准气压层上会输出高度和温度，以及露点、相对湿度和比湿。风、垂直速度和绝对涡度也在这些层上。海平面气压、2 m 温湿和 10 m 风一并写出。反射率和几类 CAPE/CIN 也在输出里，螺旋度、可降水量和累积降水同样如此。场很多，第一次可以先缩成自己要的几项。

``configs/example.yaml`` 长这样：

.. code-block:: yaml

   fields:
     - PRMSL
     - MSLET
     - HGT
     - TMP
     - RH
     - TMP2M
     - UGRD10M
     - VGRD10M
     - REFC
     - CAPE
     - CIN
     - PWAT
   pressure_hpa: [250, 500, 850]

.. code-block:: bash

   upp-py wrfout_d01_2018-07-23_21_00_00 -o post.nc --config configs/example.yaml

在 Python 里等价的调用是：

.. code-block:: python

   from upp_py import load_config, post_wrfout, write_netcdf

   ds = post_wrfout("wrfout_d01_2018-07-23_21_00_00", load_config("configs/example.yaml"))
   write_netcdf(ds, "post.nc")

``ds`` 是一个 xarray ``Dataset``。气压层在坐标 ``isobaric``（Pa）上，水平坐标是 ``latitude`` / ``longitude``，投影写在 ``crs``。每个变量的 ``grib2_shortName`` 属性是 wgrib2 里的短名。

多个时次用 ``--times 0,1,2``。每个时次的 ``apcp`` 是从起报累到该时次的连续总量；``--accum-hours 1`` 再给出以该时次为终点的 1 小时累积。只处理一个时次时，``--time`` 仍然可用，默认第 0 个。

.. code-block:: bash

   upp-py wrfout_d01 -o post.nc --grib2 post.grib2
