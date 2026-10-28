配置与命令行
============

UPP 的 ``postxconfig-NT-WRF.txt`` 是一份给代码生成器用的场清单。这里用一份短 YAML 代替它：短名列表加上要的气压层。不解析 Fortran 控制文件的位串。

.. code-block:: yaml

   fields:
     - PRMSL
     - MSLET
     - HGT
     - TMP
     - RH
     - UGRD
     - VGRD
     - TMP2M
     - DPT2M
     - RH2M
     - UGRD10M
     - VGRD10M
     - REFC
     - REFD
     - CAPE
     - CIN
     - HLCY
     - PWAT
     - APCP
   pressure_hpa: [250, 500, 700, 850, 925, 1000]
   reflectivity_agl_m: [1000, 4000]

``CAPE`` 和 ``CIN`` 会把四种气块都算出来，用坐标 ``parcel`` 区分。``HLCY`` 带 1000 m 和 3000 m 两个深度。``REFD`` 目前只做 1000 m 和 4000 m，和 UPP 的两档一致。

命令行：

.. code-block:: bash

   upp-py wrfout_d01 -o out.nc --config fields.yaml --time 0
   upp-py wrfout_d01 -o out.nc --times 0,1,2 --accum-hours 1 --grib2 out.grib2

``--times`` 覆盖 ``--time``。``--accum-hours`` 覆盖 YAML 里的 ``accum_hours``，多写一个变量 ``apcp_accum``。``--grib2`` 在 NetCDF 之外再写一份 GRIB2，需要 ``eccodes``。

``--write-example-config`` 会把一份示例写到指定路径。未知短名会直接报错。算不了的场（没有 ``PREC_ACC_*``、没有 CRTM 系数、微物理选项对不上）不报错，写进数据集属性 ``skipped``。

短名和工作里的对应关系见 :doc:`output`。完整列表在 ``upp_py.config.KNOWN_FIELDS``。
