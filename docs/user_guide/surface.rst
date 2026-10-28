地面与近地面场
================

对应 ``SURFCE.f`` 里遮蔽层（2 m）和风速计高度（10 m）那一段，读数来自 ``INITPOST.F``。

2 m 温度、露点与相对湿度
--------------------------

WRF 的 ``TH2`` 是 2 m 位温，不是温度。遮蔽层气压在 ``SURFCE.f`` 里用最低质量层温度估：

.. math::

   p_{2m} = p_{sfc}\,\exp(-0.068283 / T_\mathrm{low})

然后

.. math::

   T_{2m} = \theta_{2m}\,(p_{2m}\times 10^{-5})^\kappa

``Q2`` 是混合比，先换成比湿。露点和相对湿度用上一章的 ``CALDWP`` / ``CALRH``，气压用 ``p_{2m}`` 而不是 ``p_{sfc}``。

10 m 风
-------

``U10``、``V10`` 直接读出，不再从模式层诊断。它们是格点相对的，UPP 对 Lambert 没有再乘 ``COSALPHA`` 转到地球坐标。NetCDF 里的长名写了 grid-relative，短名仍是 ``UGRD`` / ``VGRD``，和 GRIB 一致。

地形和地面气压
--------------

地形是 ``HGT``。写出的地面气压是经过水汽柱订正的那个，不一定等于 wrfout 里的 ``PSFC`` 诊断量。平均而言两者很接近，单点可以差到 1 hPa 量级，因为平滑水汽的做法和 WRF 自己的诊断不完全一样。

.. code-block:: python

   state.t2          # K
   state.u10, state.v10
   state.p_shelter   # 2 m 气压

和 UPP 的差别
-------------

RAP/HRRR 会改用最低模式层的比湿，或者另做时间平均风。``MODELNAME=='NCAR'`` 不走那些分支。阵风、摩擦速度和感热通量这一类地面通量没有做。
